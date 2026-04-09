import argparse
import json
import os

import matplotlib.pyplot as plt
import mne
import numpy as np
import pandas as pd
from meegkit.asr import ASR
from mne_icalabel import label_components
from scipy.stats import kurtosis

from dgame.constants import BLOCK_IDS, STEP_F_KEY
from experiment.load_experiment import Experiment
from utils.utils import _safe_float
from utils.xdf_utils import extract_eeg_stream_samples, get_xdf_stream_by_type

SCALE_FACTOR = 104.1178  # NB: seemed to be done in MoBILAB, empirically determined as 1.041177792474590e+02
EEG_REMOVE_LABELS = {"ACC128", "ACC129", "ACC130", "Packet Counter", "TRIGGER"}


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
    df["duration"] = df.get("duration", np.nan)
    if "saccAmpl" in df.columns:
        df["saccAmpl"] = df["saccAmpl"].where(df["saccAmpl"] > 0, np.nan)
    return df


def build_raw_from_xdf(xdf_file: str) -> tuple["mne.io.Raw", list[str]]:
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

    data = data / SCALE_FACTOR
    info = mne.create_info(ch_names=labels, sfreq=srate, ch_types="eeg")
    raw = mne.io.RawArray(data, info, verbose="ERROR")
    return raw, labels


def apply_kurtosis_rejection(raw: mne.io.Raw, z_threshold: float = 2.0) -> list[str]:
    data = raw.get_data()
    k = kurtosis(data, axis=1, fisher=True, bias=False)
    z = (k - np.nanmean(k)) / np.nanstd(k)
    bad_idx = np.where(z > z_threshold)[0].tolist()
    bads = sorted([raw.ch_names[i] for i in bad_idx])
    return bads


def clean_rawdata_channel_rejection(
    raw: mne.io.Raw,
    flatline_seconds: float = 5.0,
    neighbor_corr_threshold: float = 0.8,
    line_noise_z_threshold: float = 4.0,
) -> list[str]:
    """
    Approximate the *channel rejection* part of EEGLAB clean_rawdata (pop_clean_rawdata)
    used in MATLAB Step F *before* kurtosis rejection:
      - FlatlineCriterion=5 (seconds)
      - ChannelCriterion=0.8 (correlation with nearby channels)
      - LineNoiseCriterion=4 (high-frequency noise outliers)

    This does NOT attempt to do burst (ASR) cleaning; that is handled separately later.
    """
    sfreq = float(raw.info["sfreq"])
    x = raw.get_data()  # (n_channels, n_samples)
    n_channels, n_samples = x.shape

    bads: set[str] = set()

    # FlatlineCriterion: mark channels whose peak-to-peak stays ~0 for an entire 5s window.
    win = int(round(flatline_seconds * sfreq))
    if win > 1 and win <= n_samples:
        step = max(int(round(1.0 * sfreq)), 1)  # 1s stride
        eps = np.finfo(float).eps * 10
        for ch_idx, ch_name in enumerate(raw.ch_names):
            max_ptp = 0.0
            for start in range(0, n_samples - win + 1, step):
                seg = x[ch_idx, start:start + win]
                ptp = float(np.ptp(seg))
                if ptp > max_ptp:
                    max_ptp = ptp
                if max_ptp > eps:
                    break
            if max_ptp <= eps:
                bads.add(ch_name)

    # LineNoiseCriterion: identify channels with unusually high high-frequency RMS (z > thresh).
    # This is a pragmatic proxy for clean_rawdata's line-noise criterion.
    x_hf = mne.filter.filter_data(
        x,
        sfreq=sfreq,
        l_freq=45.0,
        h_freq=min(100.0, sfreq / 2 - 1.0),
        verbose="ERROR",
    )
    hf_rms = np.sqrt(np.mean(x_hf ** 2, axis=1))
    if np.nanstd(hf_rms) > 0:
        z = (hf_rms - np.nanmean(hf_rms)) / np.nanstd(hf_rms)
        for ch_name, zval in zip(raw.ch_names, z):
            if float(zval) > line_noise_z_threshold:
                bads.add(ch_name)

    # ChannelCriterion: correlation with nearby channels < threshold.
    adjacency, adj_ch_names = mne.channels.find_ch_adjacency(raw.info, ch_type="eeg")
    adjacency = adjacency.tocsr()
    name_to_idx = {name: i for i, name in enumerate(adj_ch_names)}
    for ch_name in raw.ch_names:
        if ch_name in bads or ch_name not in name_to_idx:
            continue
        i = name_to_idx[ch_name]
        neigh = adjacency.indices[adjacency.indptr[i] : adjacency.indptr[i + 1]]
        if len(neigh) == 0:
            continue
        xi = x[raw.ch_names.index(ch_name)]
        best = -1.0
        for j in neigh:
            nn = adj_ch_names[j]
            if nn not in raw.ch_names:
                continue
            xj = x[raw.ch_names.index(nn)]
            if np.std(xi) == 0 or np.std(xj) == 0:
                continue
            r = float(np.corrcoef(xi, xj)[0, 1])
            if r > best:
                best = r
        if best >= 0 and best < neighbor_corr_threshold:
            bads.add(ch_name)

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

    # Load parameters from config for channel rejection
    flatline_seconds = experiment.get_dgame_step_parameter(STEP_F_KEY, "channel_rejection", "flatline_seconds")
    neighbor_corr_threshold = experiment.get_dgame_step_parameter(STEP_F_KEY, "channel_rejection", "neighbor_corr_threshold")
    line_noise_z_threshold = experiment.get_dgame_step_parameter(STEP_F_KEY, "channel_rejection", "line_noise_z_threshold")
    kurtosis_z_threshold = experiment.get_dgame_step_parameter(STEP_F_KEY, "channel_rejection", "kurtosis_z_threshold")
    asr_cutoff = experiment.get_dgame_step_parameter(STEP_F_KEY, "asr_cutoff")

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
            xdf_file = os.path.join(
                subject_xdf_dir,
                "Director",
                f"dgame{experiment.dgame_version}_{subject_id}_Director_{block}.xdf",
            )
            raw_block, _ = build_raw_from_xdf(xdf_file)
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
            # TODO why this warning about so many annotations outside time range? e.g.
            # Limited 171 annotation(s) that were expanding outside the data range.

            raw_block.resample(250, npad="auto")
            raws.append(raw_block)
            # NB: division by sfreq (=500) if EEG times are in samples, otherwise not necessary if EEG times are in seconds
            total_offset += raw_block.n_times / raw_block.info["sfreq"]

        raw = mne.concatenate_raws(raws)
        events_df = pd.concat(all_events, axis=0, ignore_index=True)

        # Save raw before filtering (optional)
        raw_before_filtering = os.path.join(outpath, f"{subject_id}_raw_before_filtering_raw.fif")
        raw.save(raw_before_filtering, overwrite=True)

        raw_raw = raw.copy()

        # Pre-cleaning
        raw.filter(l_freq=2.0, h_freq=None, verbose="ERROR")

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
            "flatline_seconds": flatline_seconds,
            "neighbor_corr_threshold": neighbor_corr_threshold,
            "line_noise_z_threshold": line_noise_z_threshold,
        }
        logger.info(f"Identifying bad channels using parameters:\n{json.dumps(clean_raw_data_params, indent=4)}")
        clean_rawdata_bads = clean_rawdata_channel_rejection(raw, **clean_raw_data_params)
        if clean_rawdata_bads:
            logger.info(
                f"Subject <{subject_id}>: {len(clean_rawdata_bads)} bad channel(s) rejected via clean_rawdata-like criteria: "
                f"{', '.join(clean_rawdata_bads)}"
            )
            raw.drop_channels([b for b in clean_rawdata_bads if b in raw.ch_names])

        # Kurtosis-based rejection
        logger.info(f"Applying kurtosis rejection with z_thresh={kurtosis_z_threshold}...")
        bads = apply_kurtosis_rejection(raw, z_threshold=kurtosis_z_threshold)
        logger.info(f"Subject <{subject_id}>: Bad channels rejected via kurtosis criterion: {', '.join(bads)}")
        raw.drop_channels([b for b in bads if b in raw.ch_names])
        missing_chs = list(dict.fromkeys(channels_to_drop + clean_rawdata_bads + bads))

        # Low-pass filter at 100 Hz and notch at 50 Hz
        raw.filter(l_freq=None, h_freq=100.0, verbose="ERROR")
        raw.notch_filter(freqs=[50.0], verbose="ERROR")

        # ASR
        logger.info(f"Applying Artifact Subspace Reconstruction (ASR) with cutoff={asr_cutoff} ...")
        raw = apply_asr(raw, cutoff=asr_cutoff)

        # Interpolate missing channels
        # TODO not correctly interpolated here
        raw.set_eeg_reference("average", projection=False)

        pre_ica_file = os.path.join(outpath, f"{subject_id}_director_preICA_raw.fif")
        raw.save(pre_ica_file, overwrite=True)

        # ICA on downsampled copy
        logger.info("Running ICA...")
        ica_raw = raw.copy().resample(100, npad="auto")
        rank = mne.compute_rank(ica_raw, rank="info").get("eeg", None)
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

        # Classify ICs
        ic_exclude = []
        labels = label_components(ica_raw, ica, method="iclabel")
        label_df = labels["labels"]
        probs = labels["y_pred_proba"]
        brain_idx = label_df == "brain"
        ic_exclude = [i for i, is_brain in enumerate(brain_idx) if not is_brain]
        logger.info(f"Subject {subject_id}: Excluding {len(ic_exclude)} non-brain ICs")

        # Apply ICA to raw data (re-filtered)
        final_raw = raw_raw.copy()
        if missing_chs:
            final_raw.drop_channels([ch for ch in missing_chs if ch in final_raw.ch_names])
        final_raw.filter(l_freq=0.3, h_freq=20.0, verbose="ERROR")
        final_raw.notch_filter(freqs=[50.0], verbose="ERROR")

        final_raw = apply_asr(final_raw, cutoff=10.0, logger=logger)

        final_raw.set_eeg_reference("average", projection=False)
        ica.apply(final_raw, exclude=ic_exclude)
        final_raw = restore_missing_channels(final_raw, missing_chs, montage)

        all_ics_file = os.path.join(outpath, f"{subject_id}_director_allICs_raw.fif")
        final_raw.save(all_ics_file, overwrite=True)

        cleaned_file = os.path.join(outpath, f"{subject_id}_director_cleaned_raw.fif")
        final_raw.save(cleaned_file, overwrite=True)

        # Save events
        events_file = os.path.join(outpath, f"{subject_id}_director_events.csv")
        # TODO write CSV files in MATLAB version and compare with version created by Python implementation
        events_df.to_csv(events_file, index=False)

    return experiment


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Preprocess EEG data.")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(os.path.abspath(args.config))
