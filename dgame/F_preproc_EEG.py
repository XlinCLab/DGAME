import argparse
import json
import os
from collections import Counter
from dataclasses import dataclass
from logging import Logger

import matplotlib.pyplot as plt
import mne
import numpy as np
import pandas as pd
from meegkit.asr import ASR
from mne_icalabel import label_components
from pyprep import NoisyChannels
from scipy.stats import kurtosis, trim_mean
from scipy.stats.mstats import trimmed_std

from dgame.amica_utils import run_amica
from dgame.constants import BLOCK_IDS, STEP_F_KEY
from dgame.matlab_scripts.dependencies import EEGLAB_PLUGIN_PATH
from experiment.load_experiment import Experiment
from utils.utils import _safe_float
from utils.xdf_utils import extract_eeg_stream_samples, get_xdf_stream_by_type

EEG_REMOVE_LABELS = {"ACC128", "ACC129", "ACC130", "Packet Counter", "TRIGGER"}

_UV_LABELS = {"microvolt", "microvolts", "µv", "uv", "μv"}
_V_LABELS = {"volt", "volts", "v"}


@dataclass
class EEGPreprocParams:
    flatline_seconds: float
    neighbor_corr_threshold: float
    line_noise_z_threshold: float
    kurtosis_z_threshold: float
    high_pass_filter_min_hz: float
    low_pass_filter_max_hz: float
    notch_filter_hz: float
    asr_cutoff: float
    ica_downsample_hz: float
    use_amica: bool
    ica_max_iter: int
    ica_max_threads: int


def load_eeg_preproc_params(experiment: Experiment) -> EEGPreprocParams:
    """Load EEG preprocessing parameters from Experiment config."""
    return EEGPreprocParams(
        flatline_seconds = experiment.get_dgame_step_parameter(STEP_F_KEY, "channel_rejection", "flatline_seconds"),
        neighbor_corr_threshold = experiment.get_dgame_step_parameter(STEP_F_KEY, "channel_rejection", "neighbor_corr_threshold"),
        line_noise_z_threshold = experiment.get_dgame_step_parameter(STEP_F_KEY, "channel_rejection", "line_noise_z_threshold"),
        kurtosis_z_threshold = experiment.get_dgame_step_parameter(STEP_F_KEY, "channel_rejection", "kurtosis_z_threshold"),
        high_pass_filter_min_hz = experiment.get_dgame_step_parameter(STEP_F_KEY, "cleaning", "high_pass_filter_min_hz"),
        low_pass_filter_max_hz = experiment.get_dgame_step_parameter(STEP_F_KEY, "cleaning", "low_pass_filter_max_hz"),
        notch_filter_hz = experiment.get_dgame_step_parameter(STEP_F_KEY, "cleaning", "notch_filter_hz"),
        asr_cutoff = experiment.get_dgame_step_parameter(STEP_F_KEY, "cleaning", "asr_cutoff"),
        ica_downsample_hz = experiment.get_dgame_step_parameter(STEP_F_KEY, "cleaning", "ica_downsample_hz"),
        use_amica = experiment.get_dgame_step_parameter(STEP_F_KEY, "cleaning", "use_amica"),
        ica_max_iter = experiment.get_dgame_step_parameter(STEP_F_KEY, "cleaning", "ica_max_iter"),
        ica_max_threads = experiment.get_dgame_step_parameter(STEP_F_KEY, "cleaning", "ica_max_threads"),
    )


def make_events_from_words(words_df: pd.DataFrame) -> pd.DataFrame:
    df = words_df.copy()
    if "trial_time" in df.columns:
        df["trial_time"] = df["trial_time"].apply(_safe_float)
    df["type"] = df.get("pos", "")
    df["duration"] = np.nan
    if "tmax" in df.columns and "time" in df.columns:
        df["duration"] = df["tmax"] - df["time"]
    df["saccAmpl"] = np.nan
    df["fix_at"] = None
    return df


def make_events_from_fixations(fix_df: pd.DataFrame) -> pd.DataFrame:
    df = fix_df.copy()
    df["type"] = "fixation"
    df["duration"] = pd.to_numeric(df.get("duration", np.nan), errors="coerce")
    # DGAME fixation CSVs have duration in milliseconds; MNE annotations expect seconds
    df["duration"] = df["duration"] / 1000.0
    if "saccAmpl" in df.columns:
        df["saccAmpl"] = df["saccAmpl"].where(df["saccAmpl"] > 0, np.nan)
    return df


def _extract_xdf_unit(eeg_stream: dict) -> str | None:
    """Return the channel unit string from XDF stream metadata, or None if absent/inconsistent."""
    try:
        desc = eeg_stream.get("info", {}).get("desc", [])
        if isinstance(desc, list):
            desc = desc[0] if desc else {}
        channels = desc.get("channels", {}) if isinstance(desc, dict) else {}
        if isinstance(channels, list):
            channels = channels[0] if channels else {}
        channel_list = channels.get("channel", []) if isinstance(channels, dict) else []
        units = set()
        for ch in channel_list:
            if not isinstance(ch, dict):
                continue
            u = ch.get("unit")
            if isinstance(u, list):
                u = u[0] if u else None
            if u:
                units.add(str(u).strip())
        if len(units) == 1:
            return units.pop()
    except Exception:
        pass
    return None


def build_raw_from_xdf(xdf_file: str, logger) -> tuple["mne.io.Raw", list[str]]:
    eeg_stream = get_xdf_stream_by_type(stream_type="EEG", xdf_file=xdf_file)
    data, srate, labels = extract_eeg_stream_samples(eeg_stream)
    if data.ndim != 2:
        raise RuntimeError(f"Unexpected EEG data shape in {xdf_file}: {data.shape}")
    if srate <= 0:
        raise RuntimeError(f"Invalid sampling rate in {xdf_file}: {srate}")
    if len(labels) != data.shape[0]:
        labels = [f"EEG{idx+1:03d}" for idx in range(data.shape[0])]

    keep_mask = [label not in EEG_REMOVE_LABELS for label in labels]
    data = data[keep_mask, :]
    labels = [label for label, keep in zip(labels, keep_mask) if keep]

    # Determine the data unit from XDF stream metadata and convert to V for MNE
    unit = _extract_xdf_unit(eeg_stream)
    if unit is None:
        logger.warning(
            f"XDF stream in {xdf_file} has no channel unit metadata — "
            "assuming µV and converting to V"
        )
        data = data * 1e-6
    elif unit.lower() in _UV_LABELS:
        logger.info(f"XDF stream unit is '{unit}' — converting µV → V")
        data = data * 1e-6
    elif unit.lower() in _V_LABELS:
        logger.info(f"XDF stream unit is '{unit}' — no unit conversion needed")
    else:
        logger.warning(
            f"XDF stream unit '{unit}' in {xdf_file} is unrecognized — "
            "assuming µV and converting to V"
        )
        data = data * 1e-6

    info = mne.create_info(ch_names=labels, sfreq=srate, ch_types="eeg")
    raw = mne.io.RawArray(data, info, verbose="ERROR")
    return raw, labels



def apply_kurtosis_rejection(raw: mne.io.Raw, z_threshold: float = 2.0) -> list[str]:
    data = raw.get_data()
    k = kurtosis(data, axis=1, fisher=True, bias=False)
    # Trimmed normalization matching MATLAB's rejkurt with normval=2 (called via
    # pop_rejchan with 'norm','on'): sort channels by kurtosis, remove the bottom
    # and top 10%, compute mean and std of the remaining 80%, then z-score all channels
    # against those trimmed statistics. This prevents extreme outlier channels from
    # inflating the std and masking other bad channels.
    k_trimmed_mean = trim_mean(k, proportiontocut=0.1)
    k_trimmed_std = trimmed_std(k, limits=(0.1, 0.1))
    z = (k - k_trimmed_mean) / k_trimmed_std
    # Two-tailed rejection matching MATLAB's abs(kurto) > threshold
    bad_idx = np.where(np.abs(z) > z_threshold)[0].tolist()
    bad_channels = sorted([raw.ch_names[i] for i in bad_idx])
    return bad_channels


def clean_rawdata_channel_rejection(
    raw: mne.io.Raw,
    flatline_seconds: float = 5.0,
    neighbor_corr_threshold: float = 0.8,
    line_noise_z_threshold: float = 4.0,
    logger: Logger | None = None, 
) -> list[str]:
    """
    Approximate the *channel rejection* part of EEGLAB clean_rawdata (pop_clean_rawdata)
    used in MATLAB Step F *before* kurtosis rejection.
    The three criteria match the EEGLAB source (clean_flatlines.m + clean_channels.m)
    as faithfully as possible:

    FlatlineCriterion (clean_flatlines.m):
      Flag channels with any consecutive run of near-zero first differences longer than
      flatline_seconds. Uses EEGLAB's jitter threshold of 20*eps.

    LineNoiseCriterion (clean_channels.m lines 88-100):
      pyprep.NoisyChannels.find_bad_by_hfnoise() — a direct port of EEGLAB's criterion.
      Same formula: MAD(residual)/MAD(lowpassed), robust z-score (median + MAD * 1.4826).
      Flatline channels are pre-marked so they are excluded from the noise computation.

    ChannelCriterion (clean_channels.m lines 121-141):
      pyprep.NoisyChannels.find_bad_by_ransac() — a direct port of EEGLAB's RANSAC criterion.
      Parameters match EEGLAB's clean_rawdata defaults: 50 random subsets of 25% of channels,
      5-second windows, MaxBrokenTime=0.4. A fresh NoisyChannels object is created after
      LineNoiseCriterion so that flatline + HFnoise channels are both excluded from RANSAC
      predictions (pyprep's RANSAC only auto-excludes bad_by_deviation/correlation/dropout,
      not bad_by_hf_noise, so pre-marking in raw.info['bads'] is required).
    """
    sfreq = float(raw.info["sfreq"])
    x = raw.get_data()  # (n_channels, n_samples)
    n_channels, n_samples = x.shape
    bads: set[str] = set()

    # FlatlineCriterion: consecutive run of near-zero first differences (cf. EEGLAB clean_flatlines.m)
    max_jitter = 20.0 * np.finfo(float).eps
    flatline_bads = set()
    for ch_idx, ch_name in enumerate(raw.ch_names):
        is_flat = np.abs(np.diff(x[ch_idx])) < max_jitter
        if not np.any(is_flat):
            continue
        padded = np.concatenate([[False], is_flat, [False]])
        changes = np.diff(padded.astype(int))
        run_starts = np.where(changes == 1)[0]
        run_ends = np.where(changes == -1)[0]
        if len(run_starts) > 0 and np.max(run_ends - run_starts) > flatline_seconds * sfreq:
            flatline_bads.add(ch_name)
    if logger:
        logger.info(f"{len(flatline_bads)} bad channel(s) identified by FlatlineCriterion: {', '.join(sorted(list(flatline_bads)))}")
    bads.update(flatline_bads)

    # LineNoiseCriterion + ChannelCriterion via pyprep.
    # raw.info["bads"] is temporarily updated between steps so that each stage excludes
    # channels already found bad. It is restored to its original value in the finally block.
    original_bads = list(raw.info["bads"])

    # LineNoiseCriterion: find_bad_by_hfnoise uses the same formula as EEGLAB's
    # clean_channels.m: MAD(residual)/MAD(lowpassed), robust z-score with 1.4826 scaling.
    # Flatline channels are pre-marked so they are excluded from the noise computation.
    raw.info["bads"] = sorted(bads)
    nc_noise = NoisyChannels(raw, do_detrend=True, random_state=0, matlab_strict=True)
    nc_noise.find_bad_by_hfnoise(HF_zscore_threshold=line_noise_z_threshold)
    if logger:
        logger.info(f"{len(nc_noise.bad_by_hf_noise)} bad channel(s) identified by LineNoiseCriterion: {', '.join(sorted(list(nc_noise.bad_by_hf_noise)))}")
    bads.update(nc_noise.bad_by_hf_noise)

    # ChannelCriterion: RANSAC on a fresh NoisyChannels object with all prior bads
    # (flatline + HFnoise) pre-marked so they are excluded from RANSAC predictions.
    # A second object is needed because pyprep's RANSAC only auto-excludes
    # bad_by_deviation/correlation/dropout — not bad_by_hf_noise.
    raw.info["bads"] = sorted(bads)
    nc_ransac = NoisyChannels(raw, do_detrend=True, random_state=0, matlab_strict=True)
    nc_ransac.find_bad_by_ransac(
        n_samples=50,                         # EEGLAB: num_samples = 50
        sample_prop=0.25,                     # EEGLAB: subset_size = 0.25
        corr_thresh=neighbor_corr_threshold,  # EEGLAB: ChannelCriterion (default 0.8)
        frac_bad=0.4,                         # EEGLAB: MaxBrokenTime = 0.4
        corr_window_secs=5.0,                 # EEGLAB: window_len = 5 s
    )
    if logger:
        logger.info(f"{len(nc_ransac.bad_by_ransac)} bad channel(s) identified by ChannelCriterion: {', '.join(sorted(list(nc_ransac.bad_by_ransac)))}")
    bads.update(nc_ransac.bad_by_ransac)
    raw.info["bads"] = original_bads

    return sorted(bads)


def restore_missing_channels(raw: mne.io.Raw, missing_chs: list[str], montage: mne.channels.DigMontage) -> mne.io.Raw:
    if not missing_chs:
        return raw
    info = mne.create_info(ch_names=missing_chs, sfreq=raw.info["sfreq"], ch_types="eeg")
    missing = mne.io.RawArray(np.zeros((len(missing_chs), raw.n_times)), info, verbose="ERROR")
    missing.set_montage(montage, match_case=False, on_missing="ignore")
    raw.add_channels([missing], force_update_info=True)
    raw.info["bads"].extend([ch for ch in missing_chs if ch not in raw.info["bads"]])
    raw.interpolate_bads(reset_bads=True)
    return raw


def apply_asr(raw: mne.io.Raw, cutoff: float = 10.0) -> mne.io.Raw:
    # Artifact Subspace Reconstruction (ASR) using meegkit.
    # Per meegkit docs: X is shaped (n_channels, n_samples) (or 3D with trials).
    # We use the 2D path here and write the cleaned samples back onto the preloaded MNE Raw.
    sfreq = float(raw.info["sfreq"])
    x = raw.get_data()  # (n_channels, n_samples)

    asr = ASR(sfreq=sfreq, cutoff=cutoff)
    asr.fit(x)
    cleaned = asr.transform(x)

    if cleaned.shape != x.shape:
        raise RuntimeError(f"ASR output shape mismatch: expected {x.shape}, got {cleaned.shape}")

    raw._data = np.asarray(cleaned, dtype=float)
    return raw


def plot_montage(montage: mne.channels.DigMontage, outplot: str, kind: str = "topomap"):
    """Plot montage to visualize electrode positions on head."""
    fig = montage.plot(kind=kind, show=False)
    if kind == "3d":
        ax = fig.gca()
        ax.view_init(azim=70, elev=15)
    fig.savefig(outplot, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main(experiment: str | dict | Experiment) -> Experiment:
    # Initialize DGAME experiment from config
    if not isinstance(experiment, Experiment):
        from dgame.dgame import DGAME
        experiment = DGAME.from_input(experiment)
    logger = experiment.logger

    # Retrieve list of electrodes which were removed to fit eyetracking glasses
    removed_electrodes = experiment.get_dgame_step_parameter(STEP_F_KEY, "removed_electrodes")
    if removed_electrodes:
        logger.info(f"Electrodes removed to fit eyetracking glasses: {', '.join(sorted(removed_electrodes))}")
    else:
        logger.warning("No electrodes specified as removed to fit eyetracking glasses!")
    # Get any override to-remove channels
    # (e.g. due to broken electrodes or exceptional noise in specific channels for specific participants)
    channels_to_remove = experiment.get_dgame_step_parameter(STEP_F_KEY, "channels_to_remove")

    # Load parameters from config for EEG bad channel rejection and cleaning
    eeg_preproc_params = load_eeg_preproc_params(experiment)

    # Load montage from preprocessed version of standard-10-5-cap385.elp (omit first line only)
    # Matches previous handling in MATLAB that references this file from standard_BESA
    montage = mne.channels.read_custom_montage(experiment.montage_file)

    # Plot montage to visualize electrode positions on head
    montage_plot_2d = os.path.join(experiment.eeg_outdir, "montage_2d.png")
    plot_montage(montage, montage_plot_2d, kind="topomap")
    logger.info(f"Plotted 2D montage topomap to {montage_plot_2d}")
    montage_plot_3d = os.path.join(experiment.eeg_outdir, "montage_3d.png")
    plot_montage(montage, montage_plot_3d, kind="3d")
    logger.info(f"Plotted 3D montage to {montage_plot_3d}")

    for subject_id in experiment.subject_ids:
        subject_xdf_dir = os.path.join(experiment.xdf_indir, subject_id)
        outpath = os.path.join(experiment.outdir, "eeg", subject_id)
        subj_ica_outdir = os.path.join(experiment.eeg_ica_outdir, subject_id)
        os.makedirs(outpath, exist_ok=True)
        os.makedirs(subj_ica_outdir, exist_ok=True)

        raws = []
        all_events = []
        total_offset = 0.0
        for block in BLOCK_IDS:
            logger.info(f"Building EEG events for subject <{subject_id}> in block <{block}>...")
            xdf_file = os.path.join(
                subject_xdf_dir,
                "Director",
                f"dgame{experiment.dgame_version}_{subject_id}_Director_{block}.xdf",
            )
            raw_block, _ = build_raw_from_xdf(xdf_file, logger=logger)
            raw_block.set_montage(montage, match_case=False, on_missing="ignore")

            # Load events
            trialtime_filename = f"{subject_id}_words2erp_{block}_trialtime.csv"
            event_file = os.path.join(experiment.outdir, "audio", subject_id, trialtime_filename)
            words_df = pd.read_csv(event_file)
            words_events = make_events_from_words(words_df)

            fix_filename = f"fixations_times_{block}_trials.csv"
            fix_file = os.path.join(
                experiment.input_dir, "preproc", "eyetracking", "fixations", subject_id, fix_filename
            )
            fix_df = pd.read_csv(fix_file)
            fix_events = make_events_from_fixations(fix_df)

            block_events = pd.concat([words_events, fix_events], axis=0, ignore_index=True)
            block_events["onset"] = block_events["time"].astype(float) + total_offset
            block_events["duration"] = block_events.get("duration", np.nan).astype(float)
            block_events["block"] = block
            all_events.append(block_events)

            # Add annotations for basic timing
            ann = mne.Annotations(
                onset=block_events["time"].astype(float).to_numpy(),
                duration=block_events["duration"].fillna(0).to_numpy(),
                description=block_events["type"].astype(str).to_numpy(),
            )
            raw_block.set_annotations(ann)

            raw_block.resample(250, npad="auto")
            raws.append(raw_block)
            # MNE uses seconds as its time base; `total_offset` is in seconds across concatenated blocks
            total_offset += raw_block.n_times / raw_block.info["sfreq"]

        raw = mne.concatenate_raws(raws)
        events_df = pd.concat(all_events, axis=0, ignore_index=True)

        # Save raw before filtering (optional)
        raw_before_filtering = os.path.join(outpath, f"{subject_id}_raw_before_filtering_raw.fif")
        raw.save(raw_before_filtering, overwrite=True)

        # Make copy of raw EEG data before cleaning/filtering
        raw_backup = raw.copy()

        # Pre-cleaning high-pass filter
        logger.info(f"Pre-cleaning EEG data with high-pass filter at {eeg_preproc_params.high_pass_filter_min_hz} Hz...")
        raw.filter(l_freq=eeg_preproc_params.high_pass_filter_min_hz, h_freq=None, verbose="ERROR")

        # Remove channels (temporarily) and track for interpolation later
        subject_channels_to_remove = channels_to_remove.get(subject_id, [])
        if subject_channels_to_remove:
            logger.info(f"Removing {len(subject_channels_to_remove)} custom channel(s) for subject <{subject_id}>: {', '.join(sorted(subject_channels_to_remove))}")
        channels_to_drop = removed_electrodes + subject_channels_to_remove
        channels_to_drop = sorted([ch for ch in channels_to_drop if ch in raw.ch_names])
        if channels_to_drop:
            logger.info(
                f"Subject <{subject_id}>: "
                f"Dropping {len(channels_to_drop)} total channel(s) "
                f"({len(removed_electrodes)} experiment-defined, "
                f"{len(subject_channels_to_remove)} subject-specific)"
            )
            raw.drop_channels(channels_to_drop)

        # Channel rejection (approximation of clean_rawdata's channel criteria, prior to kurtosis)
        # This corresponds to the first pop_clean_rawdata call in the MATLAB pipeline where
        # BurstCriterion is off (i.e., no ASR burst cleaning yet)
        clean_raw_data_params = {
            "flatline_seconds": eeg_preproc_params.flatline_seconds,
            "neighbor_corr_threshold": eeg_preproc_params.neighbor_corr_threshold,
            "line_noise_z_threshold": eeg_preproc_params.line_noise_z_threshold,
        }
        logger.info(f"Identifying bad channels using parameters:\n{json.dumps(clean_raw_data_params, indent=4)}")
        clean_rawdata_bads = clean_rawdata_channel_rejection(raw, logger=logger, **clean_raw_data_params)
        if clean_rawdata_bads:
            logger.info(
                f"Subject <{subject_id}>: {len(clean_rawdata_bads)} bad channel(s) rejected via clean_rawdata-like criteria: "
                f"{', '.join(clean_rawdata_bads)}"
            )
            raw.drop_channels([b for b in clean_rawdata_bads if b in raw.ch_names])

        # Kurtosis-based rejection
        logger.info(f"Applying kurtosis rejection with z_threshold={eeg_preproc_params.kurtosis_z_threshold}...")
        bads = apply_kurtosis_rejection(raw, z_threshold=eeg_preproc_params.kurtosis_z_threshold)
        logger.info(f"Subject <{subject_id}>: Bad channels rejected via kurtosis criterion: {', '.join(bads)}")
        raw.drop_channels([b for b in bads if b in raw.ch_names])

        # Low-pass filter at 100 Hz and notch at 50 Hz
        logger.info(f"Cleaning EEG data with low-pass filter at {eeg_preproc_params.low_pass_filter_max_hz} Hz...")
        raw.filter(l_freq=None, h_freq=eeg_preproc_params.low_pass_filter_max_hz, verbose="ERROR")
        # NB: notch_filter is best Python/MNE equivalent to CleanLine in MATLAB EEGLAB 
        logger.info(f"Applying notch filter at {eeg_preproc_params.notch_filter_hz} Hz...")
        raw.notch_filter(freqs=[eeg_preproc_params.notch_filter_hz], verbose="ERROR")

        # ASR
        logger.info(f"Applying Artifact Subspace Reconstruction (ASR) with cutoff={eeg_preproc_params.asr_cutoff} ...")
        raw = apply_asr(raw, cutoff=eeg_preproc_params.asr_cutoff)

        # Interpolate all removed/rejected channels back in before average reference,
        # so the reference is computed over a full, symmetric head coverage, then drop them again
        all_channels = set(raw_backup.ch_names)
        kept_channels = set(raw.ch_names)
        bad_channels = sorted(all_channels - kept_channels)
        raw_for_ref = raw.copy()
        raw_for_ref = restore_missing_channels(raw_for_ref, bad_channels, montage)
        raw_for_ref = mne.add_reference_channels(raw_for_ref, ref_channels=["initialReference"], copy=False)
        raw_for_ref.set_eeg_reference("average", projection=False)
        raw_for_ref.drop_channels(["initialReference"])
        if bad_channels:
            raw_for_ref.drop_channels([ch for ch in bad_channels if ch in raw_for_ref.ch_names])
        raw = raw_for_ref

        pre_ica_file = os.path.join(outpath, f"{subject_id}_director_preICA_raw.fif")
        raw.save(pre_ica_file, overwrite=True)

        # ICA on downsampled copy
        logger.info(f"Running ICA (downsampled to {eeg_preproc_params.ica_downsample_hz} Hz)...")
        ica_raw = raw.copy().resample(eeg_preproc_params.ica_downsample_hz, npad="auto")
        # Compute rank as n_channels - 1 to account for the average reference applied earlier.
        # Average reference makes one channel a linear combination of the others, reducing
        # the true mathematical rank to n-1. MNE's rank="info" and rank="full" both ignore
        # this (returning 109 for 109 channels) because the average reference was applied as
        # projection=False (direct subtraction, no SSP projector recorded in info).
        # MATLAB achieves the same result via: sum(eig(cov(data')) > 1e-6), where the
        # near-zero eigenvalue from the average reference falls below the threshold.
        rank = ica_raw.info["nchan"] - 1
        logger.info(
            f"Subject {subject_id}: ICA rank={rank} "
            f"({ica_raw.info['nchan']} channels − 1 for average reference)"
        )

        if eeg_preproc_params.use_amica:
            amica_plugin_dir = os.path.join(
                experiment.matlab_root, EEGLAB_PLUGIN_PATH, "amica"
            )
            ica = run_amica(
                raw=ica_raw,
                n_pcs=rank,
                amica_plugin_dir=amica_plugin_dir,
                work_dir=outpath,
                max_iter=eeg_preproc_params.ica_max_iter,
                max_threads=eeg_preproc_params.ica_max_threads,
                log_prefix=f"Subject {subject_id}: ",
            )
        else:
            ica = mne.preprocessing.ICA(
                n_components=rank,
                method="infomax",
                fit_params={"extended": True},
                random_state=97,
                max_iter="auto",
            )
            ica.fit(ica_raw, verbose="ERROR")

        ica.save(os.path.join(subj_ica_outdir, f"{subject_id}_ica.fif"), overwrite=True)

        post_ica_file = os.path.join(outpath, f"{subject_id}_director_postIC_raw.fif")
        ica_raw.save(post_ica_file, overwrite=True)

        # Rebuild from the unfiltered backup with the final 0.3–20 Hz filter
        final_raw = raw_backup.copy()
        if bad_channels:
            final_raw.drop_channels([ch for ch in bad_channels if ch in final_raw.ch_names])
        final_raw.filter(l_freq=0.3, h_freq=20.0, verbose="ERROR")

        # ASR on the 0.3–20 Hz data
        logger.info(f"Applying Artifact Subspace Reconstruction (ASR) with cutoff={eeg_preproc_params.asr_cutoff} ...")
        final_raw = apply_asr(final_raw, cutoff=eeg_preproc_params.asr_cutoff)

        # Interpolate bad channels for avg ref, compute reference, then drop bad channels
        # (channel set must match the ICA decomposition)
        final_raw_for_ref = final_raw.copy()
        final_raw_for_ref = restore_missing_channels(final_raw_for_ref, bad_channels, montage)
        final_raw_for_ref = mne.add_reference_channels(final_raw_for_ref, ref_channels=["initialReference"], copy=False)
        final_raw_for_ref.set_eeg_reference("average", projection=False)
        final_raw_for_ref.drop_channels(["initialReference"])
        if bad_channels:
            final_raw_for_ref.drop_channels([ch for ch in bad_channels if ch in final_raw_for_ref.ch_names])

        # Run ICLabel on the 0.3–20 Hz final data
        labeled_ica = label_components(final_raw_for_ref, ica, method="iclabel")
        ica_labels = labeled_ica["labels"]
        label_counts = Counter(ica_labels)
        label_summary = ", ".join(f"{lbl}: {cnt}" for lbl, cnt in sorted(label_counts.items()))
        logger.info(f"Subject {subject_id}: ICLabel classification — {label_summary}")
        ic_exclude = [i for i, label in enumerate(ica_labels) if label != "brain"]
        ica_excluded_percent = round((len(ic_exclude) / len(ica_labels) * 100), 2)
        logger.info(f"Subject {subject_id}: Excluding {len(ic_exclude)} ({ica_excluded_percent}%) non-brain ICs")

        # Save all-ICs file before removal
        # Use a copy so that final_raw_for_ref is not mutated before ica.apply()
        all_ics_raw = restore_missing_channels(final_raw_for_ref.copy(), bad_channels, montage)
        all_ics_file = os.path.join(outpath, f"{subject_id}_director_allICs_raw.fif")
        all_ics_raw.save(all_ics_file, overwrite=True)

        # Remove non-brain ICs and restore bad channels for final output
        ica.apply(final_raw_for_ref, exclude=ic_exclude)
        # Final interpolation: restore the missing channels in the saved outputs
        final_raw = restore_missing_channels(final_raw_for_ref, bad_channels, montage)

        cleaned_file = os.path.join(outpath, f"{subject_id}_director_cleaned_raw.fif")
        final_raw.save(cleaned_file, overwrite=True)

        # Save events
        events_file = os.path.join(outpath, f"{subject_id}_director_events.csv")
        # TODO write CSV files in MATLAB version and compare with version created by Python implementation
        events_df.to_csv(events_file, index=False)
        logger.info(f"Saved EEG events to {events_file}")

    return experiment


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Preprocess EEG data.")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(os.path.abspath(args.config))
