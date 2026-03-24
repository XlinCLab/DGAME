import argparse
import os

import asrpy
import mne
import numpy as np
import pandas as pd
from mne_icalabel import label_components
from scipy.stats import kurtosis

from dgame.constants import BLOCK_IDS, STEP_F_KEY
from experiment.load_experiment import Experiment
from utils.utils import _ensure_list, _safe_float
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


def apply_kurtosis_rejection(raw: "mne.io.Raw", z_thresh: float = 2.0) -> list[str]:
    data = raw.get_data()
    k = kurtosis(data, axis=1, fisher=True, bias=False)
    z = (k - np.nanmean(k)) / np.nanstd(k)
    bad_idx = np.where(z > z_thresh)[0].tolist()
    bads = sorted([raw.ch_names[i] for i in bad_idx])
    return bads


def restore_missing_channels(raw: "mne.io.Raw", missing_chs: list[str], montage: "mne.channels.DigMontage") -> "mne.io.Raw":
    if not missing_chs:
        return raw
    info = mne.create_info(ch_names=missing_chs, sfreq=raw.info["sfreq"], ch_types="eeg")
    missing = mne.io.RawArray(np.zeros((len(missing_chs), raw.n_times)), info, verbose="ERROR")
    missing.set_montage(montage, match_case=False, on_missing="ignore")
    raw.add_channels([missing], force_update_info=True)
    raw.info["bads"].extend([ch for ch in missing_chs if ch not in raw.info["bads"]])
    raw.interpolate_bads(reset_bads=True)
    return raw


def apply_asr(raw: "mne.io.Raw", cutoff: float = 10.0, logger=None) -> "mne.io.Raw":
    # TODO try with meegkit to see if better/faster
    asr = asrpy.ASR(sfreq=raw.info["sfreq"], cutoff=cutoff)
    asr.fit(raw)
    return asr.transform(raw)


def main(experiment: str | dict | Experiment) -> Experiment:
    # Initialize DGAME experiment from config
    if not isinstance(experiment, Experiment):
        from dgame.dgame import DGAME
        experiment = DGAME.from_input(experiment)
    logger = experiment.logger

    # Retrieve list of electrodes which were removed to fit eyetracking glasses
    removed_electrodes = experiment.get_dgame_step_parameter(STEP_F_KEY, "removed_electrodes")
    # Get any override to-remove channels (e.g. due to broken electrodes or exceptional noise in specific channels for specific participants)
    channels_to_remove = experiment.get_dgame_step_parameter(STEP_F_KEY, "channels_to_remove")

    montage = mne.channels.make_standard_montage("standard_1005")
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
        to_drop = _ensure_list(removed_electrodes) + _ensure_list(subject_channels_to_remove)
        to_drop = sorted([ch for ch in to_drop if ch in raw.ch_names])
        if to_drop:
            logger.info(f"Subject <{subject_id}>: Removing channels: {', '.join(to_drop)}")
            raw.drop_channels(to_drop)

        # Kurtosis-based rejection  # TODO ASR before kurtosis
        bads = apply_kurtosis_rejection(raw, z_thresh=2.0)
        logger.info(f"Subject <{subject_id}>: Bad channels rejected via kurtosis criterion: {', '.join(bads)}")
        raw.drop_channels([b for b in bads if b in raw.ch_names])
        missing_chs = list(dict.fromkeys(to_drop + bads))

        # Low-pass filter at 100 Hz and notch at 50 Hz
        raw.filter(l_freq=None, h_freq=100.0, verbose="ERROR")
        raw.notch_filter(freqs=[50.0], verbose="ERROR")

        # ASR
        logger.info("Applying ASR...")
        raw = apply_asr(raw, cutoff=10.0, logger=logger)

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
        ica.save(os.path.join(subj_ica_outdir, f"{subject_id}_ica.fif"))

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
