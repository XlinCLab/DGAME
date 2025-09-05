import argparse
import glob
import logging
import os
import re
from typing import Iterable

import numpy as np
import pandas as pd

from dgame.constants import (AOI_COLUMNS, FIXATION_ID_FIELD,
                             FIXATIONS_FILE_SUFFIX, GAZE_TIMESTAMP_FIELD,
                             SURFACE_LIST)
from experiment.load_experiment import Experiment

logger = logging.getLogger(__name__)


def calculate_saccade_amplitude(start_coords: Iterable, end_coords: Iterable):
    """Compute the amplitude of a saccade from staring and ending coordinates."""
    # Convert to NumPy arrays if necessary
    if not isinstance(start_coords, np.ndarray):
        start_coords = np.array(start_coords, dtype=float)
    if not isinstance(end_coords, np.ndarray):
        end_coords = np.array(end_coords, dtype=float)

    numerator = sum(start_coords * end_coords)
    denominator = np.sqrt(sum(start_coords * start_coords)) * np.sqrt(sum(end_coords * end_coords))
    saccade_amplitude = np.acos(numerator / denominator)
    return saccade_amplitude


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
        fixation_positions.loc[idx, "saccAmpl"] = calculate_saccade_amplitude(start_coords, end_coords)

    return fixation_positions


def main(experiment: str | dict | Experiment) -> Experiment:
    # Initialize DGAME experiment from config
    if not isinstance(experiment, Experiment):
        from dgame.dgame import DGAME
        experiment = DGAME.from_input(experiment)

    # Iterate over subject directories
    subject_gaze_dirs_dict = experiment.get_subject_dirs_dict(experiment.gaze_outdir)
    subject_ids = sorted(list(subject_gaze_dirs_dict.keys()))
    logger.info(f"Processing {len(subject_ids)} subject ID(s): {', '.join(subject_ids)}")
    for subject_id, subject_gaze_dirs in subject_gaze_dirs_dict.items():
        logger.info(f"Processing subject '{subject_id}'...")

        # Load per-subject fixation position files
        fixation_data = load_fixation_files(os.path.join(experiment.surface_indir, subject_id))

        # Get per-subject gaze and fixation directories/files
        if len(subject_gaze_dirs) > 1:
            logger.warning(f">1 matching directory found for subject ID '{subject_id}'")
        subject_gaze_dir = subject_gaze_dirs[0]
        gaze_file = os.path.join(subject_gaze_dir, "gaze_positions_4analysis.csv")
        fix_subj_out = os.path.join(experiment.fixations_outdir, subject_id, "fixations_4analysis.csv")
        # Create per-subject fixation outdir in case not already done
        os.makedirs(os.path.dirname(fix_subj_out), exist_ok=True)

        # Load per-subject gaze file
        # Ensure that subject ID is read in as a string
        gaze_data = pd.read_csv(gaze_file, dtype={"subj": object})
        # Drop all surface-label columns (11, 12, 13, ...) and some other unneeded columns
        # TODO 'director' column should also be dropped, but not currently in the input file
        gaze_data = gaze_data.drop(columns=SURFACE_LIST + ["norm_pos_x", "norm_pos_y", "base_data"])
        # Combine with fixation data, merging on gaze timestamp field
        # NB: rounding to 7 digits here yields different results than in the R script; no rounding matches R script
        gaze_and_fixation_data = gaze_data.merge(fixation_data, how="left", on=GAZE_TIMESTAMP_FIELD)
        # Add "end_time" column, computed from gaze_timestamp and duration
        gaze_and_fixation_data["end_time"] = gaze_and_fixation_data[GAZE_TIMESTAMP_FIELD] + (gaze_and_fixation_data["duration"] / 1000)

        # Initialize AOI (area of interest) columns
        for aoi_column, lookup_column in AOI_COLUMNS.items():
            # Skip "aoi_other", "aoi_goal", "aoi_empty" altogether
            if lookup_column is None:
                continue

            # Initialize AOI column value as False
            gaze_and_fixation_data[aoi_column] = False

            # Retrieve surface code ("11", "23", etc.)
            codes = gaze_and_fixation_data[lookup_column].apply(lambda x: str(int(x)) if pd.notna(x) else pd.NA)

            # Filter valid surface codes
            valid_codes = gaze_and_fixation_data[FIXATION_ID_FIELD].notna() & codes.isin(SURFACE_LIST)
            valid_indices = gaze_and_fixation_data.index[valid_codes]

            # For each valid row, update gaze_and_fixation_data[aoi_name] by looking up gaze_and_fixation_data[code][row]
            def get_aoi_lookup_column_value(row_idx):
                code = codes.at[row_idx]  # noqa: B023
                # If the code column exists, get its value in this row; else NaN
                return gaze_and_fixation_data.at[row_idx, code] if code in gaze_and_fixation_data.columns else pd.NA  # noqa: B023

            gaze_and_fixation_data.loc[valid_indices, aoi_column] = [get_aoi_lookup_column_value(i) for i in valid_indices]

            # (Should already be the case but) ensure aoi_column values are boolean
            gaze_and_fixation_data[aoi_column] = gaze_and_fixation_data[aoi_column].fillna(False).astype(bool)

        # Compute "aoi_empty": True if all AOI columns are False in a row
        # but only in rows which have a fixation
        has_fixation_mask = gaze_and_fixation_data[FIXATION_ID_FIELD].notna()
        aoi_cols_present = [col for col, lookup_col in AOI_COLUMNS.items() if lookup_col is not None]
        # Mask for rows with no AOI match (all AOI cols are False or NaN)
        no_aoi_match_mask = ~gaze_and_fixation_data[aoi_cols_present].any(axis=1)
        # Set aoi_empty to True for rows with a fixation but without any AOI match
        gaze_and_fixation_data["aoi_empty"] = False
        target_rows = has_fixation_mask & no_aoi_match_mask
        gaze_and_fixation_data.loc[target_rows, "aoi_empty"] = True

        # Set aoi_fillerA and aoi_fillerB to False where aoi_target is True
        gaze_and_fixation_data.loc[gaze_and_fixation_data["aoi_target"], "aoi_fillerA"] = False
        gaze_and_fixation_data.loc[gaze_and_fixation_data["aoi_target"], "aoi_fillerB"] = False

        # Save to CSV
        gaze_and_fixation_data.to_csv(fix_subj_out, index=False)
        logger.info(f"Wrote subject {subject_id} gaze fixation CSV to {fix_subj_out}")

    return experiment


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Process fixation files from pupil player and prepare them for further processing.")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(os.path.abspath(args.config))
