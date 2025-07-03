import argparse
import logging
import os
import re
import time
from collections import defaultdict

import numpy as np
import pandas as pd
from tqdm import tqdm

from constants import (AOI_COLUMNS, AUDIO_ERP_FILE_SUFFIX, CONDITIONS,
                       DEFAULT_CONFIDENCE, ERROR_LABEL,
                       GAZE_POS_SURFACE_SUFFIX, GAZE_TIMESTAMP_FIELD,
                       NOUN_POS_LABEL, PART_OF_SPEECH_FIELD, ROUND_N,
                       SURFACE_COLUMNS, SURFACE_LIST, TIMES_FILE_SUFFIX,
                       TIMESTAMPS_FILE_SUFFIX, TRIAL_TIME_OFFSET, WORD_FIELD,
                       WORD_ID_FIELD, WORD_ONSET_FIELD)
from load_experiment import (create_experiment_outdir, get_experiment_id,
                             list_subject_files, load_config,
                             log_step_duration, parse_subject_ids,
                             subject_files_dict)
from utils import (get_continuous_indices, load_file_lines,
                   merge_dataframes_with_temp_transform, setdiff)

logger = logging.getLogger(__name__)


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


def get_per_subject_audio_and_time_files(audio_dir: str,
                                         times_dir: str,
                                         subject_id_regex: str,
                                         ) -> tuple[defaultdict, defaultdict, defaultdict]:
    audio_erp_files = subject_files_dict(
        dir=audio_dir,
        subject_regex=subject_id_regex,
        suffix=AUDIO_ERP_FILE_SUFFIX,
    )
    # Find subject times and timestamps files
    times_files = subject_files_dict(
        dir=times_dir,
        subject_regex=subject_id_regex,
        suffix=TIMES_FILE_SUFFIX,
    )
    timestamps_files = subject_files_dict(
        dir=times_dir,
        subject_regex=subject_id_regex,
        suffix=TIMESTAMPS_FILE_SUFFIX
    )

    # Ensure that the same subject IDs were found per file type
    try:
        assert set(audio_erp_files.keys()) == set(times_files.keys()) == set(timestamps_files.keys())
    except AssertionError as exc:
        raise ValueError("Unequal numbers of subject IDs found!") from exc

    # Ensure that the same numbers of files were found per subject
    for subject_id in audio_erp_files:
        try:
            assert len(audio_erp_files[subject_id]) == len(times_files[subject_id]) == len(timestamps_files[subject_id])
        except AssertionError as exc:
            raise ValueError(f"Unequal numbers of audio and/or time files found for subject ID={subject_id}") from exc

    return audio_erp_files, times_files, timestamps_files


def load_erp_file(erp_file: str) -> pd.DataFrame:
    """Load and preprocess ERP CSV file."""
    # Load ERP file data
    erp_file_data = pd.read_csv(erp_file)

    # Make sure second column is labeled as "time"
    if erp_file_data.columns[1] != WORD_ONSET_FIELD:
        renamed_columns = erp_file_data.columns.tolist()
        renamed_columns[1] = WORD_ONSET_FIELD
        erp_file_data.columns = renamed_columns

    return erp_file_data


def align_times_to_erp_word_timings(times: np.ndarray,
                                    erp_time_ids: dict,
                                    ) -> dict:
    """Align an array of timestamps with the nearest timestamp associated with a word ID.
    Returns a dictionary of time: aligned_time key-value pairs."""
    # Sort times for binary search
    sorted_times = np.sort(times)
    n_times = len(sorted_times)

    # For each new time, find nearest ERP time associated with a word ID
    # Initialize aligned word ID as 0, then replace with actual ID if aligned with a word timestamp
    word_aligned_times = {time_i: 0 for time_i in times}
    for time_i, word_id in erp_time_ids.items():
        # Find insertion position in sorted time list
        idx = np.searchsorted(sorted_times, time_i)

        # Determine neighbors
        time_before = sorted_times[idx - 1] if idx > 0 else None
        # NB: Because using insertion index, would be idx rather than idx + 1
        time_after = sorted_times[idx] if idx < n_times else None

        # Choose the nearer time
        if time_before is None:
            word_aligned_times[time_i] = erp_time_ids[time_after]
        elif time_after is None:
            word_aligned_times[time_i] = erp_time_ids[time_before]
        else:
            before_diff = abs(time_i - time_before)
            after_diff = abs(time_i - time_after)
            if before_diff == after_diff:
                logger.warning("Equal distance to previous and following ERP timestamps, defaulting to previous")
            aligned_time = time_before if before_diff <= after_diff else time_after
            word_aligned_times[aligned_time] = word_id

    return word_aligned_times


def filter_and_align_subject_gaze_data_with_audio(erp_file: str,
                                                  time_file: str,
                                                  timestamp_file: str,
                                                  raw_gaze_data: pd.DataFrame,
                                                  gaze_positions_subj: pd.DataFrame,
                                                  words_df: pd.DataFrame,
                                                  ) -> pd.DataFrame:
    # Load ERP file data
    erp_file_data = load_erp_file(erp_file)

    # Load times and timestamps files
    # NB: saved as CSV but actually just list of floats, one per line
    # timestamps file contains only 2 values (start and end)
    times, timestamps = map(load_file_lines, [time_file, timestamp_file])
    # Convert all times and timestamps to floats
    # Omit the first time entry, which is time=0
    times = np.array(times, dtype=float)[1:]
    timestamps = np.array(timestamps, dtype=float)

    # Get start and end time stamps and round to ROUND_N places
    start_timestamp = round(timestamps[0], ROUND_N)
    end_timestamp = round(timestamps[-1], ROUND_N)

    # Filter erp_file_data to only those entries whose gaze_timestamp is between the two timestamps
    filtered_gaze = raw_gaze_data[
        (raw_gaze_data[GAZE_TIMESTAMP_FIELD] >= start_timestamp) &
        (raw_gaze_data[GAZE_TIMESTAMP_FIELD] < end_timestamp)
    ].copy()
    # Add times array as new column "time" to filtered_gaze dataframe
    filtered_gaze[WORD_ONSET_FIELD] = times

    # Extract known ERP times and word IDs into dict mapping
    erp_times = np.array(erp_file_data[WORD_ONSET_FIELD], dtype=float)
    erp_word_ids = np.array(erp_file_data[WORD_ID_FIELD])
    erp_time_ids = dict(zip(erp_times, erp_word_ids))

    # Align each time to the nearest ERP time associated with a word ID
    word_aligned_times = align_times_to_erp_word_timings(times, erp_time_ids)
    # Add aligned word IDs to filtered_gaze dataframe
    filtered_gaze[WORD_ID_FIELD] = word_aligned_times.values()

    # Create temp copy of erp_file_data, renaming "time" to "audio_time"
    tmp_erp_file_data = erp_file_data.copy().rename(columns={WORD_ONSET_FIELD: "audio_time"})
    # Merge tmp_erp_file_data and filtered_gaze by word "id" column
    # Now there should be an "audio_time" column as well as "time" column
    filtered_gaze = filtered_gaze.merge(tmp_erp_file_data, on=WORD_ID_FIELD, how='left')

    # Lowercase text/word field of filtered_gaze
    filtered_gaze[WORD_FIELD] = filtered_gaze[WORD_FIELD].str.lower()

    # Add filtered_gaze to gaze_positions_s dataframe
    gaze_positions_subj = pd.concat([gaze_positions_subj, filtered_gaze], axis=0, ignore_index=True)

    # Add filtered_gaze to words_df
    words_df = pd.concat([words_df, erp_file_data], axis=0, ignore_index=True)

    return gaze_positions_subj, words_df


def add_trials_to_gaze_data(gaze_positions_subj: pd.DataFrame) -> pd.DataFrame:
    """Adds trials within time window to per-subject dataframe."""
    gaze_positions_subj["trial"] = pd.NA
    gaze_positions_subj["trial_time"] = pd.NA
    trial = 1
    noun_row_indices = gaze_positions_subj.index[
        gaze_positions_subj["condition"].isin(CONDITIONS) &
        (gaze_positions_subj[PART_OF_SPEECH_FIELD] == NOUN_POS_LABEL)
    ].tolist()

    # Iterate over noun row indices and add trial annotations for data points within trial time window
    progress_bar_n = len(noun_row_indices)
    with tqdm(total=progress_bar_n) as pbar:
        pbar.set_description("Adding trials...")
        for idx in noun_row_indices:
            row = gaze_positions_subj.loc[idx]
            pattern_id = row["pattern"]
            set_id = row["set"]
            time_at_idx = row[WORD_ONSET_FIELD]
            trial_start_time = time_at_idx - TRIAL_TIME_OFFSET
            trial_end_time = time_at_idx + TRIAL_TIME_OFFSET
            gaze_positions_subj.loc[idx, "trial"] = trial
            gaze_positions_subj.loc[idx, "trial_time"] = 0

            # Identify rows within trial time frame and within current block (coded by pattern and set columns)
            pre_indices_within_trial = gaze_positions_subj.index[
                (
                    (gaze_positions_subj["pattern"] == pattern_id) |
                    # Need to include NA as possible value for pattern and set columns and then apply secondary filter (see explanation below)
                    # Otherwise the condition and surface columns will not be set for those rows
                    (gaze_positions_subj["pattern"].isna())
                ) &
                (
                    (gaze_positions_subj["set"] == set_id) |
                    (gaze_positions_subj["set"].isna())
                ) &
                (gaze_positions_subj[WORD_ONSET_FIELD] > trial_start_time) &
                (gaze_positions_subj[WORD_ONSET_FIELD] < time_at_idx) &
                (abs(gaze_positions_subj[WORD_ONSET_FIELD] - time_at_idx) <= TRIAL_TIME_OFFSET)
            ].tolist()
            post_indices_within_trial = gaze_positions_subj.index[
                (
                    (gaze_positions_subj["pattern"] == pattern_id) |
                    (gaze_positions_subj["pattern"].isna())
                ) &
                (
                    (gaze_positions_subj["set"] == set_id) |
                    (gaze_positions_subj["set"].isna())
                ) &
                (gaze_positions_subj[WORD_ONSET_FIELD] < trial_end_time) &
                (gaze_positions_subj[WORD_ONSET_FIELD] > time_at_idx) &
                (abs(gaze_positions_subj[WORD_ONSET_FIELD] - time_at_idx) <= TRIAL_TIME_OFFSET)
            ].tolist()
            # Above filtering logic does not fully exclude indices from different blocks, which have NA values for "pattern" and "set" fieds
            # Solution is to take only the longest continuous set of indices around the central time point,
            # which would mimic forward and backward while loop behavior.
            # When a non-continuous index is found, it necessarily comes from different block
            # Exclude all indices which are non-contiguous to the central set
            pre_indices_within_trial = get_continuous_indices(idx, pre_indices_within_trial, direction="pre")
            post_indices_within_trial = get_continuous_indices(idx, post_indices_within_trial, direction="post")
            indices_to_set = pre_indices_within_trial + post_indices_within_trial

            # Set column values for data points within trial window
            gaze_positions_subj.loc[indices_to_set, "trial"] = trial
            gaze_positions_subj.loc[indices_to_set, "trial_time"] = gaze_positions_subj.loc[indices_to_set, WORD_ONSET_FIELD] - time_at_idx
            for col_name in SURFACE_COLUMNS + ["condition"]:
                gaze_positions_subj.loc[indices_to_set, col_name] = row[col_name]

            # Increment trial
            trial += 1

            # Update progress bar
            pbar.update(1)

    return gaze_positions_subj


def validate_surface_annotation(value: str | int) -> str:
    """Validates that a surface annotation value belongs to the list of surfaces.
    Corrects single-digit surfaces to corresponding double-digit surfaces, e.g. "2" -> "22"."""
    if pd.isna(value):
        return value
    value = str(int(value))
    if len(value) == 1:
        value = value * 2
    try:
        assert value in SURFACE_LIST
    except AssertionError as exc:
        raise ValueError(f"Invalid surface annotation: {value}") from exc
    return value


def add_surface_aoi_annotations(gaze_positions_subj: pd.DataFrame) -> pd.DataFrame:
    """Add annotations for surface areas of interest for a single subject."""

    logger.info("Annotating surface areas of interest...")
    surface_condition_indices = gaze_positions_subj.index[
        (gaze_positions_subj["trial_time"].notna()) &
        (gaze_positions_subj["condition"].isin(CONDITIONS)) &
        (gaze_positions_subj["surface"].notna())
    ].tolist()

    # Cast AOI columns to nullable Boolean type to avoid warnings when setting NaN values
    # NB: otherwise will cause an error in a future version of pandas
    for aoi_field in AOI_COLUMNS.keys():
        gaze_positions_subj[aoi_field] = gaze_positions_subj[aoi_field].astype('boolean')

    progress_bar_n = len(surface_condition_indices)
    with tqdm(total=progress_bar_n) as pbar:
        pbar.set_description("Annotating surface areas of interest...")
        for idx in surface_condition_indices:
            row = gaze_positions_subj.loc[idx]
            target, goal, otherTarget, fillerA, fillerB, competitor, otherCompetitor = map(
                validate_surface_annotation,
                [
                    row["surface"],
                    row["target_location"],
                    row["targetB_surface"],
                    row["fillerA_surface"],
                    row["fillerB_surface"],
                    row["surface_competitor"],
                    row["compB_surface"],
                ]
            )
            empty_surfaces = setdiff(SURFACE_LIST, {target, competitor, otherTarget, otherCompetitor, fillerA, fillerB})

            # Add area of interest (AOI) flags
            def set_aoi_flag(surface, aoi_field):
                if surface != ERROR_LABEL and pd.notna(surface):
                    gaze_positions_subj.at[idx, aoi_field] = row[surface]  # noqa: B023

            set_aoi_flag(target, 'aoi_target')
            set_aoi_flag(goal, 'aoi_goal')
            set_aoi_flag(otherTarget, 'aoi_otherTarget')
            set_aoi_flag(competitor, 'aoi_comp')
            set_aoi_flag(otherCompetitor, 'aoi_otherComp')
            set_aoi_flag(fillerA, 'aoi_fillerA')
            set_aoi_flag(fillerB, 'aoi_fillerB')

            # Mark whether any AOI surface is empty
            gaze_positions_subj.at[idx, "aoi_empty"] = any(
                pd.notna(surface) and gaze_positions_subj.at[idx, surface] is True
                for surface in empty_surfaces
            )

            # Update progress bar
            pbar.update(1)

    return gaze_positions_subj


def main(config: str | dict) -> dict:
    start_time = time.time()
    # Load experiment config
    if isinstance(config, str):
        config = load_config(config)
    experiment_id = get_experiment_id(config)

    # Retrieve paths to inputs
    input_dir = config["data"]["input"]["root"]
    audio_dir = config["data"]["input"]["audio_dir"]
    audio_indir = os.path.join(input_dir, audio_dir)
    gaze_dir = config["data"]["input"]["gaze_dir"]
    gaze_indir = os.path.join(input_dir, gaze_dir)
    gaze_pos_file = os.path.join(gaze_indir, "gaze_positions.csv")
    times_dir = config["data"]["input"]["times_dir"]
    times_indir = os.path.join(input_dir, times_dir)
    surface_dir = config["data"]["input"]["surfaces_dir"]
    surface_indir = os.path.join(input_dir, surface_dir)
    surface_files = list_subject_files(dir=surface_indir, subject_regex="^", suffix=GAZE_POS_SURFACE_SUFFIX)

    # Output paths
    output_dir = create_experiment_outdir(config, experiment_id)
    audio_outdir = os.path.join(output_dir, audio_dir)
    gaze_outdir = os.path.join(output_dir, gaze_dir)
    gaze_all_out = os.path.join(gaze_outdir, "gaze_positions_all_4analysis.csv")

    # Load surface and object position data
    logger.info("Loading surface fixation position data...")
    surface_pos_data = load_and_combine_surface_files(surface_files)
    # Round gaze_timestamp field in order to enable merge
    surface_pos_data[f"rounded_{GAZE_TIMESTAMP_FIELD}"] = round(surface_pos_data[GAZE_TIMESTAMP_FIELD], ROUND_N)

    # Load gaze file (columns of interest only)
    logger.info(f"Loading gaze data from {gaze_pos_file}")
    raw_gaze_data = pd.read_csv(
        gaze_pos_file,
        usecols=[
            GAZE_TIMESTAMP_FIELD,
            "world_index",
            "confidence",
            "norm_pos_x",
            "norm_pos_y",
            "base_data",
        ]
    )
    # Round gaze_timestamp field of raw_gaze_data to ROUND_N places
    raw_gaze_data[GAZE_TIMESTAMP_FIELD] = raw_gaze_data[GAZE_TIMESTAMP_FIELD].astype(float).round(ROUND_N)

    # Get selected subject IDs
    _, subject_id_regex = parse_subject_ids(config["experiment"]["subjects"])
    # Find per-subject audio ERP and time/timestamp files
    logger.info("Loading per-subject audio and timing files...")
    subj_audio_erp_dict, subj_times_dict, subj_timestamps_dict = get_per_subject_audio_and_time_files(
        audio_dir=audio_indir,
        times_dir=times_indir,
        subject_id_regex=subject_id_regex,
    )
    # Get subject IDs (should be identical for all 3 file types)
    subject_ids = sorted(list(subj_audio_erp_dict.keys()))
    logger.info(f"Processing {len(subject_ids)} subject ID(s): {', '.join(subject_ids)}")

    # Iterate over subjects and combine word data with gaze data
    logger.info("Loading word data and combining with gaze data...")
    gaze_positions_all = pd.DataFrame()
    for subject_id in subject_ids:
        logger.info(f"Processing subject '{subject_id}'...")

        # Create per-subject subject output directories
        for outdir_i in {audio_outdir, gaze_outdir}:
            subj_outdir_i = os.path.join(outdir_i, subject_id)
            os.makedirs(subj_outdir_i, exist_ok=True)
        # Designate per-subject output file paths
        word_outfile = os.path.join(audio_outdir, subject_id, "all_words_4analysis.csv")
        gaze_before_words_file = os.path.join(gaze_outdir, subject_id, "gaze_positions_before_words.csv")
        gaze_subj_out = os.path.join(gaze_outdir, subject_id, "gaze_positions_4analysis.csv")
        tmp_gaze_s = os.path.join(gaze_outdir, subject_id, "tmp_gaze_positions.csv")

        # Initialize empty dataframe to contain all processed gaze data per subject
        gaze_positions_subj = pd.DataFrame()
        words_df = pd.DataFrame()

        # Load word data and combine with gaze data
        audio_erp_files = subj_audio_erp_dict[subject_id]
        times_files = subj_times_dict[subject_id]
        timestamps_files = subj_timestamps_dict[subject_id]
        for erp_file, time_file, timestamp_file in zip(audio_erp_files, times_files, timestamps_files):
            logger.debug(f"ERP file: {os.path.basename(erp_file)}")
            logger.debug(f"Time file: {os.path.basename(time_file)}")
            logger.debug(f"Timestamp file: {os.path.basename(timestamp_file)}")
            gaze_positions_subj, words_df = filter_and_align_subject_gaze_data_with_audio(
                erp_file=erp_file,
                time_file=time_file,
                timestamp_file=timestamp_file,
                raw_gaze_data=raw_gaze_data,
                gaze_positions_subj=gaze_positions_subj,
                words_df=words_df,
            )

        # Merge gaze positions and surface positions by timestamp
        # Add another column with rounded timestamps in order to merge, floating point timestamps may not match exactly
        gaze_positions_subj = merge_dataframes_with_temp_transform(
            left_df=gaze_positions_subj,
            right_df=surface_pos_data,
            on=GAZE_TIMESTAMP_FIELD,
            how="left",
            transform=lambda x: round(x, ROUND_N),
            transform_left=True,
            # NB: not necessary, to transform right df
            # surface_pos_data is transformed only once in advance to avoid re-transforming in each subject loop iteration
            transform_right=False,
            temp_column_name=f"rounded_{GAZE_TIMESTAMP_FIELD}",
        )

        # Add subject ID to gaze_positions_subj dataframe and write CSV outfiles
        gaze_positions_subj["subj"] = subject_id
        gaze_positions_subj["subj"] = gaze_positions_subj["subj"].astype(str)
        gaze_positions_subj.to_csv(gaze_before_words_file, index=False)
        logger.info(f"Wrote subject {subject_id} data to {gaze_before_words_file}")
        words_df.to_csv(word_outfile, index=False)
        logger.info(f"Wrote word data to {word_outfile}")

        # Add trials
        gaze_positions_subj = add_trials_to_gaze_data(gaze_positions_subj)
        # Write CSV file with added trials
        gaze_positions_subj.to_csv(tmp_gaze_s, index=False)
        logger.info(f"Wrote data with trial annotations to {tmp_gaze_s}")

        # Initialize new area of interest (aoi) columns as False
        for new_column in AOI_COLUMNS.keys():
            gaze_positions_subj[new_column] = False
        # Set trackloss column to boolean value, whether confidence < DEFAULT_CONFIDENCE
        gaze_positions_subj["trackloss"] = gaze_positions_subj["confidence"] < DEFAULT_CONFIDENCE

        # Check if participants looked at a surface or not at a given time point
        gaze_positions_subj = add_surface_aoi_annotations(gaze_positions_subj)

        # Sort dataframe by gaze_timestamp
        gaze_positions_subj = gaze_positions_subj.sort_values(by=GAZE_TIMESTAMP_FIELD, ascending=True).reset_index(drop=True)

        # Write CSV file
        gaze_positions_subj.to_csv(gaze_subj_out, index=False)
        logger.info(f"Wrote per-subject gaze file (subject = {subject_id}) to {gaze_subj_out}")

        # Add per-subject trial data into running dataframe for all subjects
        gaze_positions_all = pd.concat([gaze_positions_all, gaze_positions_subj], axis=0, ignore_index=True)

    # Write full gaze positions CSV file for all subjects
    gaze_positions_all.to_csv(gaze_all_out, index=False)
    logger.info(f"Wrote full gaze file (all subjects) to {gaze_all_out}")

    # Log duration of this step in run config
    log_step_duration(config, start_time, step_id="Ca_preproc_et_data")

    return config


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Creates trials from the continuous gaze file.")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(args.config)
