import argparse
import os
from typing import Any

import mne
import numpy as np
import pandas as pd

from dgame.constants import (CONFLICT_LABEL, NO_CONFLICT_LABEL, STEP_F_KEY,
                             TRIAL_TIME_OFFSET)
from experiment.input_validation import InputValidationError, ensure_columns_exist
from experiment.load_experiment import Experiment


def update_fixation_events_df(events: pd.DataFrame) -> pd.DataFrame:
    """
    Setup the fixations for analysis: loop through all events and find fixations to add
    condition information retrieved from the noun (type == 'N') of the corresponding trial.

    The logic is based on:
    - `trial_time`: time relative to noun onset (negative = before noun, positive = after noun)
    - if within +/- TRIAL_TIME_OFFSET seconds of noun onset, search for the corresponding noun
    - search window is +/- 200 events (arbitrary but intended to cover all fixations per trial)
    """
    events = events.copy()
    events["type"] = events["type"].astype(str)

    if "fix_at" not in events.columns:
        events["fix_at"] = "NA"
    events["fix_at"] = events["fix_at"].astype(str)

    if "trial_time" not in events.columns:
        events["trial_time"] = np.nan
    events["trial_time"] = pd.to_numeric(events["trial_time"], errors="coerce")

    for idx in range(len(events)):
        if events.at[idx, "type"] != "fixation":
            continue

        # Summarize non-target fixations.
        # Any fixation not explicitly labeled target/other is treated as "elsewhere".
        fix_at = events.at[idx, "fix_at"]
        if fix_at not in ("target", "other"):
            events.at[idx, "fix_at"] = "elsewhere"

        tt = events.at[idx, "trial_time"]
        if pd.isna(tt):
            continue

        if tt < 0 and tt >= -TRIAL_TIME_OFFSET:
            # If it has trial_time < 0, i.e. occurs before noun onset:
            # go forward a few events and see if you find a noun; stop at the first one.
            # (200 is arbitrary so we are sure to get all nouns in all trials)
            for find_idx in range(1, 201):
                cand_idx = idx + find_idx
                if cand_idx >= len(events):
                    break
                if events.at[cand_idx, "type"] == "N":
                    for col in ("condition", "set", "pattern"):
                        if col in events.columns:
                            events.at[idx, col] = events.at[cand_idx, col]
                    break
        elif tt > 0 and tt <= TRIAL_TIME_OFFSET:
            # Do the same for fixations after noun onset:
            # go backward a few events and see if you find a noun; stop at the first one.
            for find_idx in range(1, 201):
                cand_idx = idx - find_idx
                if cand_idx < 0:
                    break
                if events.at[cand_idx, "type"] == "N":
                    for col in ("condition", "set", "pattern"):
                        if col in events.columns:
                            events.at[idx, col] = events.at[cand_idx, col]
                    break

        raw_condition = events.at[idx, "condition"] if "condition" in events.columns else None
        condition = str(raw_condition) if pd.notna(raw_condition) else ""
        if condition not in (CONFLICT_LABEL, NO_CONFLICT_LABEL):
            # Set all fixations without a condition to a different type.
            events.at[idx, "type"] = "other_fixation"
            if "condition_old" not in events.columns:
                events["condition_old"] = "NA"
            events.at[idx, "condition_old"] = condition if condition else "NA"
            events.at[idx, "condition"] = "NA"

    # Ensure fixation fix_at labeling is consistent
    fix_mask = events["type"] == "fixation"
    events.loc[fix_mask & ~events["fix_at"].isin(["target", "other"]), "fix_at"] = "elsewhere"

    return events


def main(experiment: str | dict | Experiment) -> Experiment:
    # Initialize DGAME experiment from config
    if not isinstance(experiment, Experiment):
        from dgame.dgame import DGAME
        experiment = DGAME.from_input(experiment)
    logger = experiment.logger

    # Get experiment subject IDs and their corresponding EEG directory paths
    subject_eeg_dirs_dict = experiment.get_subject_dirs_dict(experiment.eeg_outdir)

    # Pre-unfold processing from MNE .fif + events CSV
    per_subject_inputs: dict[str, dict[str, Any]] = {}
    for subject_id, subject_dirs in subject_eeg_dirs_dict.items():
        subject_eeg_dir = subject_dirs[0]
        cleaned_fif = f"{subject_id}_director_cleaned_raw.fif"
        cleaned_fif_path = os.path.join(subject_eeg_dir, cleaned_fif)
        if not os.path.exists(cleaned_fif_path):
            raise InputValidationError(f"Missing cleaned EEG file (from step {STEP_F_KEY}): {cleaned_fif_path}")

        events_csv = os.path.join(subject_eeg_dir, f"{subject_id}_director_events.csv")
        if not os.path.exists(events_csv):
            raise InputValidationError(f"Missing events CSV (from step {STEP_F_KEY}): {events_csv}")

        logger.info(f"Loading cleaned EEG for subject {subject_id}: {cleaned_fif_path}")
        raw = mne.io.read_raw_fif(cleaned_fif_path, preload=True, verbose="ERROR")

        events = pd.read_csv(events_csv)
        # Timestamp normalization:
        # Step F (Python/MNE) is expected to provide `onset` in seconds on a single absolute
        # timeline across concatenated blocks (i.e., per-block offsets already applied).
        ensure_columns_exist(events, ["type", "onset"], events_csv)
        events["onset"] = pd.to_numeric(events["onset"], errors="coerce")
        events = events.dropna(subset=["onset"]).sort_values("onset").reset_index(drop=True)

        # Ensure required columns exist for Julia step
        for col in ("condition", "trial", "trial_time", "saccAmpl", "fix_at", "set", "pattern"):
            if col not in events.columns:
                logger.warning(f"Events CSV file {events_csv} is missing column <{col}>; adding column with all NA values.")
                events[col] = np.nan

        # For fixation events close to noun onset (trial_time within +/- TRIAL_TIME_OFFSET),
        # search forward/backward up to 200 events to find the corresponding noun (N) and copy
        # its condition metadata onto the fixation.
        events = update_fixation_events_df(events)

        # Unfold.jl expects `latency` in 1-indexed samples (EEGLAB convention).
        # MNE onset times are in seconds on a 0-indexed timeline (first sample = t=0),
        # so onset * sfreq gives a 0-indexed sample number. Adding 1 converts to the
        # 1-indexed convention that Julia arrays and Unfold.jl use internally.
        events["latency"] = np.round(events["onset"] * raw.info["sfreq"]).astype(float) + 1

        # Inject boundary events from MNE annotations into the pre-unfold events CSV.
        # mne.concatenate_raws inserts a 'BAD boundary' annotation at each block junction;
        # these are preserved in the .fif file and signal that the EEG is discontinuous
        # at that point. The Julia artifact detection function respects these by refusing
        # to create a rejection window that straddles a boundary — matching EEGLAB's
        # uf_continuousArtifactDetect behaviour. Without them the detector would compare
        # signal segments from different recording blocks, producing spurious rejections
        # or missing real artifacts at the joins. 'EDGE boundary' falls at the same sample
        # as 'BAD boundary', so only one entry per junction is needed.
        boundary_onsets = sorted({
            ann["onset"] for ann in raw.annotations
            if ann["description"] == "BAD boundary"
        })
        if boundary_onsets:
            boundary_rows = pd.DataFrame({
                "type": "boundary",
                "onset": boundary_onsets,
                "latency": np.round(np.array(boundary_onsets) * raw.info["sfreq"]).astype(float) + 1,
            })
            events = pd.concat([events, boundary_rows], ignore_index=True).sort_values("onset").reset_index(drop=True)
        logger.info(f"Adding {len(boundary_onsets)} boundary event(s) to pre-unfold events for subject {subject_id}")

        # Create output directory for unfold results
        unfold_outpath = os.path.join(subject_eeg_dir, "unfold_out")
        os.makedirs(unfold_outpath, exist_ok=True)
        # Save events CSV prior to unfold
        pre_unfold_events_csv = os.path.join(unfold_outpath, f"{subject_id}_events_pre_unfold.csv")
        events.to_csv(pre_unfold_events_csv, index=False)
        logger.info(f"Saved pre-unfold EEG event CSV for subject {subject_id}: {pre_unfold_events_csv}")

        per_subject_inputs[subject_id] = {
            # MNE stores EEG in Volts; Unfold.jl and EEGLAB both expect microvolts.
            # The artifact threshold (150 µV) and all output betas depend on this scale.
            "data": raw.get_data() * 1e6,
            "srate": float(raw.info["sfreq"]),
            "events_csv": pre_unfold_events_csv,
            "chan_names": raw.ch_names,
            "outpath": unfold_outpath,
        }

    # Run Unfold analysis in Julia
    logger.info("Running unfold analysis in Julia...")
    julia_path = experiment.julia_dir
    unfold_julia_script = os.path.join(julia_path, "unfold_step_g.jl")
    end_log_cmd = experiment.set_julia_logfile(unfold_julia_script)
    jl = experiment.julia_interface
    # NB: workaround to silence non-critical error about missing extension CategoricalArraysExt
    jl.seval('import Logging; Logging.disable_logging(Logging.Error)')
    jl.seval('Base.retry_load_extensions()')
    jl.seval(f'include("{unfold_julia_script}")')
    # Reset logging defaults after workaround
    jl.seval('Logging.disable_logging(Logging.BelowMinLevel)')
    for subject_id, inputs in per_subject_inputs.items():
        logger.info(f"Running Unfold.jl for subject {subject_id}...")
        chan_names = jl.Vector[jl.String](inputs["chan_names"])
        jl.run_unfold_step_g_from_arrays(
            inputs["data"],
            inputs["srate"],
            inputs["events_csv"],
            chan_names,
            inputs["outpath"],
            subject_id,
        )
    # Close the Julia logfile and reset logging
    jl.seval(end_log_cmd)

    return experiment


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Compute rERPs.")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(os.path.abspath(args.config))
