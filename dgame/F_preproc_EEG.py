import argparse
import json
import os
from collections import Counter
from dataclasses import dataclass
from typing import Iterable

import matplotlib.pyplot as plt
import mne
import numpy as np
import pandas as pd
from meegkit.asr import ASR
from mne.channels import DigMontage
from mne_icalabel import label_components
from pyprep import NoisyChannels
from scipy.stats import kurtosis, trim_mean
from scipy.stats.mstats import trimmed_std

from dgame.amica_utils import run_amica
from dgame.constants import BLOCK_IDS
from dgame.pipeline import STEP_F_KEY
from experiment.input_validation import InputValidationError
from experiment.load_experiment import Experiment
from utils.utils import _safe_float
from utils.xdf_utils import (extract_eeg_stream_samples, fill_stream_gaps,
                             get_xdf_stream_by_type)

EEG_REMOVE_LABELS = {"ACC128", "ACC129", "ACC130", "Packet Counter", "TRIGGER"}

_UV_LABELS = {"microvolt", "microvolts", "µv", "uv", "μv"}
_V_LABELS = {"volt", "volts", "v"}


@dataclass
class EEGPreprocParams:
    block_resample: float
    flatline_seconds: float
    neighbor_corr_threshold: float
    line_noise_z_threshold: float
    kurtosis_z_threshold: float
    high_pass_filter_min_hz: float
    low_pass_filter_max_hz: float
    notch_filter_hz: float
    asr_cutoff: float
    ica_downsample_hz: float
    ica_method: str
    ica_max_iter: int
    ica_max_threads: int
    iclabel_high_pass_filter_min_hz: float
    iclabel_low_pass_filter_max_hz: float


class ExperimentEEGHandler:
    def __init__(self, experiment):
        # Initialize experiment, logger, and EEG output directory
        self.experiment = experiment
        self.logger = experiment.logger
        self.eeg_outdir = experiment.eeg_outdir

    def debug(self, msg: str):
        self.logger.debug(msg)

    def info(self, msg: str):
        self.logger.info(msg)

    def warning(self, msg: str):
        self.logger.warning(msg)

    def error(self, msg: str):
        self.logger.error(msg)


class EEGPipeline(ExperimentEEGHandler):
    def __init__(self, experiment: Experiment):
        super().__init__(experiment)

        # Validate that all expected per-subject/per-block input files exist
        self.validate_inputs()

        # Load experiment montage file
        self.montage = self.load_montage()

        # Parse pre-designated electrodes to exclude
        self.removed_electrodes, self.channels_to_remove = self.parse_excluded_channels()

        # Load parameters from config for EEG bad channel rejection and cleaning
        self.params = self.load_eeg_preproc_params()

    def get_xdf_file(self, subject_id: str, block: int) -> str:
        return os.path.join(
            self.experiment.xdf_indir,
            subject_id,
            "Director",
            f"dgame{self.experiment.dgame_version}_{subject_id}_Director_{block}.xdf",
        )

    def get_trialtime_file(self, subject_id: str, block: int) -> str:
        return os.path.join(
            self.experiment.audio_outdir,
            subject_id,
            f"{subject_id}_words2erp_{block}_trialtime.csv",
        )

    def get_fixation_file(self, subject_id: str, block: int) -> str:
        return os.path.join(
            self.experiment.fixations_outdir,
            subject_id,
            f"fixations_times_{block}_trials.csv",
        )

    def validate_inputs(self) -> None:
        """Validate that all expected EEG pipeline input files exist for every subject and block.
        Collects every missing file before raising error, rather than failing on the first one found."""
        missing_files = []
        for subject_id in self.experiment.subject_ids:
            for block in BLOCK_IDS:
                for filepath in (
                    self.get_xdf_file(subject_id, block),
                    self.get_trialtime_file(subject_id, block),
                    self.get_fixation_file(subject_id, block),
                ):
                    if not os.path.exists(filepath):
                        missing_files.append(filepath)
        if missing_files:
            missing_list = "\n".join(f"  - {filepath}" for filepath in missing_files)
            raise InputValidationError(
                f"Missing {len(missing_files)} expected input file(s) for EEG preprocessing:\n{missing_list}"
            )

    def parse_excluded_channels(self) -> tuple[list, dict[str, list]]:
        """Retrieve list of electrodes/channels which were removed to fit eyetracking glasses
        as well as subject-specific electrodes/channels to exclude."""
        removed_electrodes = self.experiment.get_dgame_step_parameter(STEP_F_KEY, "removed_electrodes")
        if removed_electrodes:
            self.info(f"Electrodes removed to fit eyetracking glasses: {', '.join(sorted(removed_electrodes))}")
        else:
            self.warning("No electrodes specified as removed to fit eyetracking glasses!")
        # Get any override to-remove channels
        # (e.g. due to broken electrodes or exceptional noise in specific channels for specific participants)
        channels_to_remove = self.experiment.get_dgame_step_parameter(STEP_F_KEY, "channels_to_remove")
        return removed_electrodes, channels_to_remove

    def load_montage(self) -> DigMontage:
        # Load montage from preprocessed version of standard-10-5-cap385.elp (omit first line only)
        # Matches previous handling in MATLAB that references this file from standard_BESA
        return mne.channels.read_custom_montage(self.experiment.montage_file)
    
    def plot_montage(self, outdir: str = None) -> None:
        """Plot montage to visualize electrode positions on head."""
        if outdir is None:
            outdir = self.eeg_outdir
        montage_plot_2d = os.path.join(outdir, "montage_2d.png")
        plot_montage(self.montage, montage_plot_2d, kind="topomap")
        self.info(f"Plotted 2D montage topomap to {montage_plot_2d}")
        montage_plot_3d = os.path.join(outdir, "montage_3d.png")
        plot_montage(self.montage, montage_plot_3d, kind="3d")
        self.info(f"Plotted 3D montage to {montage_plot_3d}")

    def load_eeg_preproc_params(self) -> EEGPreprocParams:
        """Load EEG preprocessing parameters from Experiment config."""
        experiment = self.experiment
        return EEGPreprocParams(
            block_resample = experiment.get_dgame_step_parameter(STEP_F_KEY, "block_resample"),
            flatline_seconds = experiment.get_dgame_step_parameter(STEP_F_KEY, "channel_rejection", "flatline_seconds"),
            neighbor_corr_threshold = experiment.get_dgame_step_parameter(STEP_F_KEY, "channel_rejection", "neighbor_corr_threshold"),
            line_noise_z_threshold = experiment.get_dgame_step_parameter(STEP_F_KEY, "channel_rejection", "line_noise_z_threshold"),
            kurtosis_z_threshold = experiment.get_dgame_step_parameter(STEP_F_KEY, "channel_rejection", "kurtosis_z_threshold"),
            high_pass_filter_min_hz = experiment.get_dgame_step_parameter(STEP_F_KEY, "cleaning", "high_pass_filter_min_hz"),
            low_pass_filter_max_hz = experiment.get_dgame_step_parameter(STEP_F_KEY, "cleaning", "low_pass_filter_max_hz"),
            notch_filter_hz = experiment.get_dgame_step_parameter(STEP_F_KEY, "cleaning", "notch_filter_hz"),
            asr_cutoff = experiment.get_dgame_step_parameter(STEP_F_KEY, "cleaning", "asr_cutoff"),
            ica_downsample_hz = experiment.get_dgame_step_parameter(STEP_F_KEY, "ica", "ica_downsample_hz"),
            ica_method = experiment.get_dgame_step_parameter(STEP_F_KEY, "ica", "method").lower(),
            ica_max_iter = experiment.get_dgame_step_parameter(STEP_F_KEY, "ica", "amica", "ica_max_iter"),
            ica_max_threads = experiment.get_dgame_step_parameter(STEP_F_KEY, "ica", "amica", "ica_max_threads"),
            iclabel_high_pass_filter_min_hz = experiment.get_dgame_step_parameter(STEP_F_KEY, "ica", "iclabel", "high_pass_filter_min_hz"),
            iclabel_low_pass_filter_max_hz = experiment.get_dgame_step_parameter(STEP_F_KEY, "ica", "iclabel", "low_pass_filter_max_hz"),
        )

    def run(self):
        """Run EEG preprocessing pipeline for all available subjects in experiment."""
        for subject_id in self.experiment.subject_ids:
            self.info(f"Preprocessing EEG data for subject <{subject_id}>...")
            subject_eeg_processor = SubjectEEGPreprocessor(
                experiment=self.experiment,
                subject_id=subject_id,
            )
            subject_eeg_processor.run_pipeline()


class SubjectEEGPreprocessor(EEGPipeline):
    def __init__(self, experiment: Experiment, subject_id: str):
        super().__init__(experiment)
        self.subject_id = subject_id
        # Create subject output directories
        self.outdir = os.path.join(self.eeg_outdir, subject_id)
        os.makedirs(self.outdir, exist_ok=True)
        # Initialize lookup dictionary of EEG versions at given steps
        self.original_eeg_key = "START"
        self.eeg_versions =  {}
        # Sets of channels
        self.channel_names = set()
        self.bad_channels = set()

    def save_eeg_version(self, raw_eeg: mne.io.Raw, label: str) -> None:
        """Save a given EEG data version under a specified label."""
        if label in self.eeg_versions:
            raise ValueError(f"Label <{label}> has already used for an EEG version!")
        self.eeg_versions[label] = raw_eeg

    def set_original_eeg(self, raw_eeg: mne.io.Raw, label="START") -> None:
        """Save original version of EEG data before any processing and its set of included channels."""
        self.save_eeg_version(raw_eeg, label=label)
        self.original_eeg_key = label
        self.channel_names = set(raw_eeg.ch_names)
    
    def get_original_eeg(self):
        """Retrieve the original EEG data before processing."""
        return self.eeg_versions.get(self.original_eeg_key)

    def get_missing_channels(self, raw_eeg: mne.io.Raw) -> list:
        """Returns a list of channels included in the saved original EEG data but absent from the specified EEG object."""
        kept_channels = set(raw_eeg.ch_names)
        missing_channels = sorted(self.channel_names - kept_channels)
        return missing_channels
    
    def set_missing_channels_as_bad(self, raw_eeg: mne.io.Raw) -> set:
        """Identifies missing channels with respect to the original EEG data and marks these as bad channels."""
        self.bad_channels.update(self.get_missing_channels(raw_eeg))
        return self.bad_channels
    
    def drop_bad_channels(self, raw_eeg: mne.io.Raw, bad_channels: Iterable = None) -> mne.io.Raw:
        """Drop designated bad channels from a raw EEG dataset."""
        if bad_channels is None:
            bad_channels = self.bad_channels
        if not bad_channels:
            self.warning("No bad channels found to drop!")
            return raw_eeg
        raw_eeg.drop_channels(
            [ch for ch in bad_channels if ch in raw_eeg.ch_names]
        )
        return raw_eeg

    def restore_missing_channels(self, raw_eeg: mne.io.Raw, missing_chs: list[str]= None) -> mne.io.Raw:
        if missing_chs is None:
            missing_chs = self.get_missing_channels(raw_eeg)
        if not missing_chs:
            return raw_eeg
        info = mne.create_info(ch_names=list(missing_chs), sfreq=raw_eeg.info["sfreq"], ch_types="eeg")
        missing = mne.io.RawArray(np.zeros((len(missing_chs), raw_eeg.n_times)), info, verbose="ERROR")
        missing.set_montage(self.montage, match_case=False, on_missing="ignore")
        raw_eeg.add_channels([missing], force_update_info=True)
        raw_eeg.info["bads"].extend([ch for ch in missing_chs if ch not in raw_eeg.info["bads"]])
        raw_eeg.interpolate_bads(reset_bads=True)
        return raw_eeg

    def run_pipeline(self):
        """Run EEG preprocessing pipeline for single subject."""

        # Build raw subject EEG blocks
        raw, events_df = self.build_subject_eeg_blocks()
        # Save events CSV
        # TODO write CSV files in MATLAB version and compare with version created by Python implementation
        self.write_events(events_df)
        # Save raw EEG before filtering
        raw.save(
            os.path.join(self.outdir, f"{self.subject_id}_raw_before_filtering.fif"),
            overwrite=True,
        )

        # Create backup of raw EEG data before cleaning/filtering
        self.set_original_eeg(raw.copy())

        # Pre-cleaning high-pass filter
        self.info(f"Pre-cleaning EEG data with high-pass filter at {self.params.high_pass_filter_min_hz} Hz...")
        raw.filter(l_freq=self.params.high_pass_filter_min_hz, h_freq=None, verbose="ERROR")

        # Remove channels (temporarily) and track for interpolation later
        raw, _ = self.drop_subject_channels(raw_eeg=raw)

        # Channel rejection (approximation of clean_rawdata's channel criteria, prior to kurtosis)
        raw = self.clean_rawdata_channel_rejection(raw_eeg=raw)

        # Kurtosis-based rejection
        raw = self.kurtosis_channel_rejection(raw_eeg=raw)

        # Identify bad channels up to this point
        self.set_missing_channels_as_bad(raw)

        # Low-pass filter (default = 100 Hz)
        self.info(f"Cleaning EEG data with low-pass filter at {self.params.low_pass_filter_max_hz} Hz...")
        raw.filter(l_freq=None, h_freq=self.params.low_pass_filter_max_hz, verbose="ERROR")

        # Notch filter (default = 50 Hz)
        # NB: notch_filter is best Python/MNE equivalent to CleanLine in MATLAB EEGLAB 
        self.info(f"Applying notch filter at {self.params.notch_filter_hz} Hz...")
        raw.notch_filter(freqs=[self.params.notch_filter_hz], verbose="ERROR")

        # ASR
        self.info(f"Applying Artifact Subspace Reconstruction (ASR) with cutoff={self.params.asr_cutoff} ...")
        raw = apply_asr(raw, cutoff=self.params.asr_cutoff)

        # Interpolate all removed/rejected channels back in before average reference,
        # so the reference is computed over a full, symmetric head coverage, then drop them again
        raw = self.interpolate_bads_and_set_average_reference(raw_eeg=raw)

        # Run ICA
        pre_ica_file = os.path.join(self.outdir, f"{self.subject_id}_director_preICA_raw.fif")
        raw.save(pre_ica_file, overwrite=True)
        ica = self.run_ica(raw_eeg=raw)

        # Rebuild from the unfiltered backup with the final filter (default: 0.3–20 Hz)
        self.info(f"Filtering original raw EEG data to {self.params.iclabel_high_pass_filter_min_hz}-{self.params.iclabel_low_pass_filter_max_hz} Hz...")
        final_raw = self.get_original_eeg().copy()
        final_raw = self.drop_bad_channels(raw_eeg=final_raw)
        final_raw.filter(
            l_freq=self.params.iclabel_high_pass_filter_min_hz,
            h_freq=self.params.iclabel_low_pass_filter_max_hz,
            verbose="ERROR"
        )

        # ASR on the final filtered (default = 0.3-20 Hz) data
        self.info(f"Applying Artifact Subspace Reconstruction (ASR) with cutoff={self.params.asr_cutoff} ...")
        final_raw = apply_asr(final_raw, cutoff=self.params.asr_cutoff)

        # Interpolate bad channels for avg ref, compute reference, then drop bad channels
        # (channel set must match the ICA decomposition)
        final_raw_for_ref = self.interpolate_bads_and_set_average_reference(raw_eeg=final_raw)

        # Save all-ICs file before removal of non-brain ICs
        # Use a copy so that raw_eeg is not mutated before ica.apply()
        all_ics_raw = self.restore_missing_channels(
            raw_eeg=final_raw_for_ref.copy(),
            missing_chs=self.bad_channels,
        )
        all_ics_file = os.path.join(self.outdir, f"{self.subject_id}_director_allICs_raw.fif")
        all_ics_raw.save(all_ics_file, overwrite=True)

        # Run ICLabel on the final filtered (default = 0.3-20 Hz) data
        final_raw_for_ref = self.label_ica(raw_eeg=final_raw_for_ref, ica=ica)

        # Final interpolation: restore the missing channels in the saved outputs
        final_raw = self.restore_missing_channels(
            raw_eeg=final_raw_for_ref,
            missing_chs=self.bad_channels,
        )

        # Save final cleaned EEG file
        cleaned_file = os.path.join(self.outdir, f"{self.subject_id}_director_cleaned_raw.fif")
        final_raw.save(cleaned_file, overwrite=True)

    def build_subject_eeg_blocks(self) -> tuple[mne.io.Raw, pd.DataFrame]:
        raws = []
        all_events = []
        total_offset = 0.0
        for block in BLOCK_IDS:
            self.info(f"Building EEG events for subject <{self.subject_id}> in block <{block}>...")
            xdf_file = self.get_xdf_file(self.subject_id, block)
            raw_block, _, gap_events = self.build_raw_from_xdf(xdf_file)
            raw_block.set_montage(self.montage, match_case=False, on_missing="ignore")
            if len(gap_events) > 0:
                self.warning(
                    f"Subject <{self.subject_id}> block <{block}>: found {len(gap_events)} "
                    f"gap(s) totalling {gap_events['duration'].sum():.2f}s of gap-filled samples"
                )

            # Load events
            event_file = self.get_trialtime_file(self.subject_id, block)
            words_df = pd.read_csv(event_file)
            words_events = make_events_from_words(words_df)

            fix_file = self.get_fixation_file(self.subject_id, block)
            fix_df = pd.read_csv(fix_file)
            fix_events = make_events_from_fixations(fix_df)

            block_events = pd.concat([words_events, fix_events, gap_events], axis=0, ignore_index=True)
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
            self.info(f"Resampling block {block} at {self.params.block_resample} Hz...")
            raw_block.resample(self.params.block_resample, npad="auto")
            raws.append(raw_block)
            # MNE uses seconds as its time base; `total_offset` is in seconds across concatenated blocks
            total_offset += raw_block.n_times / raw_block.info["sfreq"]

        raw = mne.concatenate_raws(raws)
        events_df = pd.concat(all_events, axis=0, ignore_index=True)
        return raw, events_df

    def build_raw_from_xdf(self, xdf_file: str) -> tuple["mne.io.Raw", list[str], pd.DataFrame]:
        eeg_stream = get_xdf_stream_by_type(stream_type="EEG", xdf_file=xdf_file)
        # Insert zero-valued samples at dropped-sample gaps so sample count/timing stays
        # correct (otherwise a gap silently shifts everything after it earlier than its real
        # time, desyncing this stream from the separately-timed word/fixation events below).
        # The filled-in samples are tracked via "is_gap" and returned as "BAD_gap" events so
        # they end up as MNE annotations and get excluded from epoching/ICA downstream.
        eeg_stream = fill_stream_gaps(eeg_stream)
        data, srate, labels = extract_eeg_stream_samples(eeg_stream)
        gap_events = make_events_from_gaps(
            eeg_stream["is_gap"],
            np.asarray(eeg_stream["time_stamps"], dtype=np.float64),
            srate=srate,
        )
        if data.ndim != 2:
            raise RuntimeError(f"Unexpected EEG data shape in {xdf_file}: {data.shape}")
        if srate <= 0:
            raise RuntimeError(f"Invalid sampling rate in {xdf_file}: {srate}")

        # Some recordings run at a sustained rate that deviates from the nominal_srate in the
        # XDF metadata (e.g. Bluetooth-based amplifiers under packet loss), without tripping
        # fill_stream_gaps's discrete-gap detection (no single interval is anomalously large -
        # the whole stream is just uniformly slower). Building the Raw object at the nominal
        # rate would then silently drift out of sync with the separately-timed word/fixation
        # events over the course of a several-minute block. Using the actual measured rate
        # keeps whatever samples we do have correctly time-aligned; it cannot recover the
        # missing samples themselves (their exact position within the stream can't be
        # determined once pyxdf's jitter-removal has smoothed the real, uneven arrival times
        # into a uniform series) - that data loss is real and is not fixed by this.
        timestamps = np.asarray(eeg_stream["time_stamps"], dtype=np.float64)
        effective_srate = (len(timestamps) - 1) / (timestamps[-1] - timestamps[0])
        if abs(effective_srate - srate) / srate > 0.01:
            self.warning(
                f"XDF stream in {xdf_file} has effective rate {effective_srate:.2f}Hz, "
                f"deviating >1% from nominal {srate}Hz - using effective rate for MNE timing"
            )
        srate = effective_srate
        if len(labels) != data.shape[0]:
            labels = [f"EEG{idx+1:03d}" for idx in range(data.shape[0])]

        keep_mask = [label not in EEG_REMOVE_LABELS for label in labels]
        data = data[keep_mask, :]
        labels = [label for label, keep in zip(labels, keep_mask) if keep]

        # Determine the data unit from XDF stream metadata and convert to V for MNE
        unit = _extract_xdf_unit(eeg_stream)
        if unit is None:
            self.warning(
                f"XDF stream in {xdf_file} has no channel unit metadata — "
                "assuming µV and converting to V"
            )
            data = data * 1e-6
        elif unit.lower() in _UV_LABELS:
            self.info(f"XDF stream unit is '{unit}' — converting µV → V")
            data = data * 1e-6
        elif unit.lower() in _V_LABELS:
            self.info(f"XDF stream unit is '{unit}' — no unit conversion needed")
        else:
            self.warning(
                f"XDF stream unit '{unit}' in {xdf_file} is unrecognized — "
                "assuming µV and converting to V"
            )
            data = data * 1e-6

        info = mne.create_info(ch_names=labels, sfreq=srate, ch_types="eeg")
        raw = mne.io.RawArray(data, info, verbose="ERROR")
        return raw, labels, gap_events

    def write_events(self, events_df: pd.DataFrame) -> None:
        """Write a Pandas DataFrame containing annotated events from XDF file to CSV."""
        events_file = os.path.join(self.outdir, f"{self.subject_id}_director_events.csv")
        events_df.to_csv(events_file, index=False)
        self.info(f"Wrote EEG events to {events_file}")

    def drop_subject_channels(self, raw_eeg: mne.io.Raw) -> tuple[mne.io.Raw, list]:
        subject_channels_to_remove = self.channels_to_remove.get(self.subject_id, [])
        if subject_channels_to_remove:
            self.info(f"Removing {len(subject_channels_to_remove)} custom channel(s) for subject <{self.subject_id}>: {', '.join(sorted(subject_channels_to_remove))}")
        channels_to_drop = set(self.removed_electrodes + subject_channels_to_remove)
        channels_to_drop = sorted([ch for ch in channels_to_drop if ch in raw_eeg.ch_names])
        if channels_to_drop:
            self.info(
                f"Subject <{self.subject_id}>: "
                f"Dropping {len(channels_to_drop)} total channel(s) "
                f"({len(self.removed_electrodes)} experiment-defined, "
                f"{len(subject_channels_to_remove)} subject-specific)"
            )
            raw_eeg.drop_channels(channels_to_drop)
        return raw_eeg, channels_to_drop

    def clean_rawdata_channel_rejection(self, raw_eeg: mne.io.Raw) -> mne.io.Raw:
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
        clean_rawdata_params = {
            "flatline_seconds": self.params.flatline_seconds,
            "neighbor_corr_threshold": self.params.neighbor_corr_threshold,
            "line_noise_z_threshold": self.params.line_noise_z_threshold,
        }
        self.info(f"Identifying bad channels using parameters:\n{json.dumps(clean_rawdata_params, indent=4)}")
        sfreq = float(raw_eeg.info["sfreq"])
        x = raw_eeg.get_data()  # (n_channels, n_samples)
        bads: set[str] = set()

        # FlatlineCriterion: consecutive run of near-zero first differences (cf. EEGLAB clean_flatlines.m)
        max_jitter = 20.0 * np.finfo(float).eps
        flatline_bads = set()
        for ch_idx, ch_name in enumerate(raw_eeg.ch_names):
            is_flat = np.abs(np.diff(x[ch_idx])) < max_jitter
            if not np.any(is_flat):
                continue
            padded = np.concatenate([[False], is_flat, [False]])
            changes = np.diff(padded.astype(int))
            run_starts = np.where(changes == 1)[0]
            run_ends = np.where(changes == -1)[0]
            if len(run_starts) > 0 and np.max(run_ends - run_starts) > self.params.flatline_seconds * sfreq:
                flatline_bads.add(ch_name)
        self.info(f"{len(flatline_bads)} bad channel(s) identified by FlatlineCriterion: {', '.join(sorted(list(flatline_bads)))}")
        bads.update(flatline_bads)

        # LineNoiseCriterion + ChannelCriterion via pyprep.
        # raw.info["bads"] is temporarily updated between steps so that each stage excludes
        # channels already found bad. It is restored to its original value in the finally block.
        original_bads = list(raw_eeg.info["bads"])

        # LineNoiseCriterion: find_bad_by_hfnoise uses the same formula as EEGLAB's
        # clean_channels.m: MAD(residual)/MAD(lowpassed), robust z-score with 1.4826 scaling.
        # Flatline channels are pre-marked so they are excluded from the noise computation.
        raw_eeg.info["bads"] = sorted(bads)
        nc_noise = NoisyChannels(raw_eeg, do_detrend=True, random_state=0, matlab_strict=True)
        nc_noise.find_bad_by_hfnoise(HF_zscore_threshold=self.params.line_noise_z_threshold)
        self.info(f"{len(nc_noise.bad_by_hf_noise)} bad channel(s) identified by LineNoiseCriterion: {', '.join(sorted(list(nc_noise.bad_by_hf_noise)))}")
        bads.update(nc_noise.bad_by_hf_noise)

        # ChannelCriterion: RANSAC on a fresh NoisyChannels object with all prior bads
        # (flatline + HFnoise) pre-marked so they are excluded from RANSAC predictions.
        # A second object is needed because pyprep's RANSAC only auto-excludes
        # bad_by_deviation/correlation/dropout — not bad_by_hf_noise.
        raw_eeg.info["bads"] = sorted(bads)
        nc_ransac = NoisyChannels(raw_eeg, do_detrend=True, random_state=0, matlab_strict=True)
        nc_ransac.find_bad_by_ransac(
            n_samples=50,                         # EEGLAB: num_samples = 50
            sample_prop=0.25,                     # EEGLAB: subset_size = 0.25
            corr_thresh=self.params.neighbor_corr_threshold,  # EEGLAB: ChannelCriterion (default 0.8)
            frac_bad=0.4,                         # EEGLAB: MaxBrokenTime = 0.4
            corr_window_secs=5.0,                 # EEGLAB: window_len = 5 s
        )
        self.info(f"{len(nc_ransac.bad_by_ransac)} bad channel(s) identified by ChannelCriterion: {', '.join(sorted(list(nc_ransac.bad_by_ransac)))}")
        bads.update(nc_ransac.bad_by_ransac)
        raw_eeg.info["bads"] = original_bads

        # Drop identified bad channels
        if bads:
            self.info(
                f"{len(bads)} bad channel(s) rejected via clean_rawdata-like criteria: "
                f"{', '.join(sorted(bads))}"
            )
            raw_eeg = self.drop_bad_channels(raw_eeg, bad_channels=bads)

        return raw_eeg

    def kurtosis_channel_rejection(self, raw_eeg: mne.io.Raw) -> mne.io.Raw:
        self.info(f"Applying kurtosis rejection with z_threshold={self.params.kurtosis_z_threshold}...")
        bad_channels = apply_kurtosis_rejection(raw_eeg, z_threshold=self.params.kurtosis_z_threshold)
        self.info(f"Bad channels rejected via kurtosis criterion: {', '.join(bad_channels)}")
        raw_eeg = self.drop_bad_channels(raw_eeg, bad_channels=bad_channels)
        return raw_eeg
    
    def interpolate_bads_and_set_average_reference(self, raw_eeg: mne.io.Raw) -> mne.io.Raw:
        """Interpolate all removed/rejected channels back in before average reference so that
        the reference is computed over a full, symmetric head coverage, then drop them again."""
        bad_channels = self.get_missing_channels(raw_eeg)
        raw_for_ref = raw_eeg.copy()
        raw_for_ref = self.restore_missing_channels(raw_eeg=raw_for_ref, missing_chs=bad_channels)
        raw_for_ref = mne.add_reference_channels(raw_for_ref, ref_channels=["initialReference"], copy=False)
        raw_for_ref.set_eeg_reference("average", projection=False)
        raw_for_ref.drop_channels(["initialReference"])
        raw_for_ref = self.drop_bad_channels(raw_for_ref, bad_channels=bad_channels)
        return raw_for_ref

    def run_ica(self, raw_eeg: mne.io.Raw) -> mne.preprocessing.ICA:
        # Create output directory for ICA results
        subj_ica_outdir = os.path.join(self.experiment.eeg_ica_outdir, self.subject_id)
        os.makedirs(subj_ica_outdir, exist_ok=True)

        # Run ICA on downsampled copy
        self.info(f"Downsampling EEG to {self.params.ica_downsample_hz} Hz for ICA...")
        ica_raw = raw_eeg.copy().resample(self.params.ica_downsample_hz, npad="auto")
        # Compute rank as n_channels - 1 to account for the average reference applied earlier.
        # Average reference makes one channel a linear combination of the others, reducing
        # the true mathematical rank to n-1. MNE's rank="info" and rank="full" both ignore
        # this (returning 109 for 109 channels) because the average reference was applied as
        # projection=False (direct subtraction, no SSP projector recorded in info).
        # MATLAB achieves the same result via: sum(eig(cov(data')) > 1e-6), where the
        # near-zero eigenvalue from the average reference falls below the threshold.
        rank = ica_raw.info["nchan"] - 1
        self.info(f"ICA rank={rank} ({ica_raw.info['nchan']} channels − 1 for average reference)")

        if self.params.ica_method == "amica":
            self.info("Running AMICA...")
            amica_plugin_dir = os.path.abspath(
                self.experiment.get_analysis_parameter("dependencies", "amica", "dir")
            )
            if not os.path.isdir(amica_plugin_dir):
                raise FileNotFoundError(
                    f"AMICA directory not found at '{amica_plugin_dir}'. "
                    f"Run ./install_amica.sh to install it, or set "
                    f"analysis.dependencies.amica.dir in your config."
                )
            amica_log_path = os.path.join(self.experiment.logdir, "steps", f"{STEP_F_KEY}_amica", f"{self.subject_id}_amica.log")
            self.info(f"AMICA output log: {amica_log_path}")
            ica = run_amica(
                raw=ica_raw,
                n_pcs=rank,
                amica_plugin_dir=amica_plugin_dir,
                work_dir=self.outdir,
                max_iter=self.params.ica_max_iter,
                max_threads=self.params.ica_max_threads,
                log_prefix=f"Subject {self.subject_id}: ",
                amica_log_path=amica_log_path,
            )
        else:
            self.info("Running ICA with Extended Infomax...")
            ica = mne.preprocessing.ICA(
                n_components=rank,
                method="infomax",
                fit_params={"extended": True},
                random_state=97,
                max_iter="auto",
            )
            ica.fit(ica_raw, verbose="ERROR")

        # Save ICA artifacts
        ica.save(os.path.join(subj_ica_outdir, f"{self.subject_id}_ica.fif"), overwrite=True)
        post_ica_file = os.path.join(self.outdir, f"{self.subject_id}_director_postIC_raw.fif")
        ica_raw.save(post_ica_file, overwrite=True)

        return ica

    def label_ica(self,raw_eeg: mne.io.Raw, ica: mne.preprocessing.ICA) -> mne.io.Raw:
        # Run ICLabel to yield IC labels within EEG signal (brain, heartbeat, muscle, etc.)
        self.info("Running ICLabel to separate brain from non-brain components...")
        labeled_ica = label_components(raw_eeg, ica, method="iclabel")
        ica_labels = labeled_ica["labels"]
        label_counts = Counter(ica_labels)
        label_summary = ", ".join(f"{label}: {count}" for label, count in sorted(label_counts.items()))
        self.info(f"ICLabel classification — {label_summary}")

        # Get non-brain ICs to exclude
        ic_exclude = [i for i, label in enumerate(ica_labels) if label != "brain"]
        ica_excluded_percent = round((len(ic_exclude) / len(ica_labels) * 100), 2)
        self.info(f"Excluding {len(ic_exclude)} ({ica_excluded_percent}%) non-brain ICs")

        # Remove non-brain ICs
        ica.apply(raw_eeg, exclude=ic_exclude)
        return raw_eeg


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


def make_events_from_gaps(is_gap: np.ndarray, timestamps: np.ndarray, srate: float) -> pd.DataFrame:
    """Build one "BAD_gap" event per contiguous run of gap-filled (zero-inserted) samples,
    in the same {"time", "duration", "type"} schema as make_events_from_words/_fixations, so
    they flow into the same MNE annotations and get excluded from epoching/ICA by MNE's
    standard "BAD_"-prefix convention (e.g. via reject_by_annotation=True)."""
    is_gap = np.asarray(is_gap, dtype=bool)
    if not is_gap.any():
        return pd.DataFrame(columns=["time", "duration", "type"])

    padded = np.concatenate(([False], is_gap, [False]))
    run_starts = np.flatnonzero(np.diff(padded.astype(int)) == 1)
    run_ends = np.flatnonzero(np.diff(padded.astype(int)) == -1) - 1

    period = 1.0 / srate if srate > 0 else 0.0
    onset = timestamps[run_starts] - timestamps[0]
    # +period so the interval fully covers the last gap-filled sample, not just up to its onset
    duration = timestamps[run_ends] - timestamps[run_starts] + period

    return pd.DataFrame({"time": onset, "duration": duration, "type": "BAD_gap"})


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

    # Initialize EEGPipeline instance from loaded experiment and its config
    eeg_pipeline = EEGPipeline(experiment)

    # Plot montage to visualize electrode positions on head
    eeg_pipeline.plot_montage()

    # Run full EEG preprocessing pipeline for all subjects
    eeg_pipeline.run()

    return experiment


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Preprocess EEG data.")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(os.path.abspath(args.config))
