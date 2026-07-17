import glob
import logging
import os
import re

import numpy as np
import pandas as pd

from dgame.constants import (COLUMN_DATA_TYPES, FIXATION_ID_FIELD,
                             FIXATION_TIMES_TRIALS_SUFFIX,
                             FIXATIONS_FILE_SUFFIX, GAZE_TIMESTAMP_FIELD,
                             ROUND_N, WORD_END_FIELD, WORD_ONSET_FIELD)
from dgame.eyetracking.saccades import compute_saccade_amplitude
from experiment.load_experiment import Experiment
from experiment.test_subjects import list_subject_files


def load_filtered_gaze_data(experiment: Experiment) -> pd.DataFrame:
    """Load the aggregate gaze positions file (output of step Ca) and filter to valid trial
    data points with non-negative word durations.

    Shared by Da_compute_trialtime and Da_gaze_language_timecourse_stats, which both derive
    their input from this same filtered dataframe but otherwise have no dependency on one another.
    """
    # TODO confirm that this should use the aggregated file, not individual per subject
    gaze_infile = os.path.join(experiment.gaze_outdir, "gaze_positions_all_4analysis.csv")

    # Ensure that the subj column is read as a string; e.g. '02' will be read in as an integer unless otherwise specified in dtype arg
    gaze_data = pd.read_csv(gaze_infile, dtype=COLUMN_DATA_TYPES)
    gaze_data = gaze_data.loc[
        gaze_data["condition"].notna() &
        gaze_data["trial_time"].notna() &
        (~gaze_data["trackloss"]) &
        gaze_data["subj"].isin(experiment.subject_ids)
    ].drop_duplicates()

    # Add column "duration": tmax - time (rounded to ROUND_N digits)
    gaze_data["duration"] = (gaze_data[WORD_END_FIELD] - gaze_data[WORD_ONSET_FIELD]).round(ROUND_N)

    # Filter valid durations and then drop duration column
    gaze_data = gaze_data[gaze_data["duration"] >= 0].drop(columns="duration")

    return gaze_data


def load_and_combine_surface_files(surface_file_list: list) -> pd.DataFrame:
    """Load a series of surface fixation CSV files and combine into single dataframe."""
    surface_pos_data = None
    for surface_file in surface_file_list:
        surface_id = re.search(r"_(\d+)\.csv", os.path.basename(surface_file)).group(1)
        tmp = pd.read_csv(surface_file)
        # Drop all columns except GAZE_TIMESTAMP_FIELD ("gaze_timestamp") and "on_surf"
        tmp = tmp.drop(columns=[
            "confidence",
            "world_index",
            "x_norm",
            "y_norm",
            "x_scaled",
            "y_scaled",
            "world_timestamp",
        ])
        # Rename "on_surf" column to the surface ID
        tmp = tmp.rename(columns={"on_surf": surface_id})
        # Merge each successive file into combined dataframe
        if surface_pos_data is not None:
            surface_pos_data = surface_pos_data.merge(tmp, on=GAZE_TIMESTAMP_FIELD, how='left')
        else:
            surface_pos_data = tmp
    return surface_pos_data


def load_fixation_files(surface_dir):
    """Retrieve and load surface fixation files."""
    fixation_file_regex = re.compile(FIXATIONS_FILE_SUFFIX)
    # NB: fixation files sequence needs to be (deterministically) sorted
    fixation_files = sorted([
        filepath for filepath in glob.glob(os.path.join(surface_dir, "*"))
        if fixation_file_regex.match(os.path.basename(filepath))
    ])

    # Load all files into single dataframe
    fixation_positions = None
    for fix_file in fixation_files:
        tmp_fixation_positions = pd.read_csv(fix_file)
        # Group by fixation_id column and take first row of each
        tmp_fixation_positions = tmp_fixation_positions.groupby(FIXATION_ID_FIELD, as_index=False).first()
        # Rename "on_surf" column to surface ID
        surface_id = fixation_file_regex.search(fix_file).group(1)
        tmp_fixation_positions = tmp_fixation_positions.rename(columns={"on_surf": surface_id})
        # Merge each successive file into combined dataframe, merging on fixation_id column
        if fixation_positions is not None:
            # Keep only surface and "fixation_id" column
            tmp_fixation_positions = tmp_fixation_positions[[FIXATION_ID_FIELD, surface_id]]
            fixation_positions = fixation_positions.merge(tmp_fixation_positions, on=FIXATION_ID_FIELD, how='left')
        else:
            fixation_positions = tmp_fixation_positions

    # Drop "world_timestamp" column
    fixation_positions = fixation_positions.drop(columns=["world_timestamp"])

    # Rename 'start_timestamp' to 'gaze_timestamp'
    fixation_positions = fixation_positions.rename(columns={'start_timestamp': GAZE_TIMESTAMP_FIELD})

    # Add "saccAmpl" column
    fixation_positions = fixation_positions.reset_index(drop=True)
    fixation_positions["saccAmpl"] = 0.0
    for idx, row in fixation_positions.iterrows():
        if idx == 0:  # start from second row in order to look one row backward
            continue
        previous_row = fixation_positions.loc[idx - 1]
        start_coords = np.array([previous_row["norm_pos_x"], previous_row["norm_pos_y"]], dtype=float)
        end_coords = np.array([row["norm_pos_x"], row["norm_pos_y"]], dtype=float)
        fixation_positions.loc[idx, "saccAmpl"] = compute_saccade_amplitude(start_coords, end_coords)

    return fixation_positions


def load_fixation_times_trials_files(
        subj_fixation_dirs_dict: dict,
        logger: logging.Logger,
        ) -> pd.DataFrame:
    """Loads fixation times trials files from selected subjects into a single dataframe."""
    fixation_times_trials_df = pd.DataFrame()
    for subject_id, subj_fixation_dirs in subj_fixation_dirs_dict.items():
        if len(subj_fixation_dirs) > 1:
            logger.warning(f">1 matching directory found for subject ID '{subject_id}'")
        subj_fixation_dir = subj_fixation_dirs[0]
        fixation_times_trials_files = list_subject_files(
            dir=subj_fixation_dir,
            subject_regex=r"^",
            suffix=FIXATION_TIMES_TRIALS_SUFFIX
        )
        for fixation_time_trial_file in fixation_times_trials_files:
            # Ensure that the subj column is read as a string, e.g. '02' will be read in as an integer
            data = pd.read_csv(fixation_time_trial_file, dtype={"subj": object})
            # Ensure that the subj column has a single entry and matches subject_id
            assert len(data["subj"].unique()) == 1
            assert data["subj"].unique()[0] == subject_id
            fixation_times_trials_df = pd.concat([fixation_times_trials_df, data], axis=0, ignore_index=True)
    return fixation_times_trials_df
