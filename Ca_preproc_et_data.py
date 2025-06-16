import argparse
import logging
import os
import re
from collections import defaultdict

import numpy as np
import pandas as pd
from tqdm import tqdm

from constants import (AUDIO_ERP_FILE_SUFFIX, GAZE_POS_SURFACE_SUFFIX,
                       GAZE_TIMESTAMP_FIELD, ROUND_N, TIMES_FILE_SUFFIX,
                       TIMESTAMPS_FILE_SUFFIX, WORD_FIELD, WORD_ID_FIELD,
                       WORD_ONSET_FIELD)
from load_experiment import (list_subject_files, load_config,
                             load_object_positions_data, parse_subject_ids,
                             subject_files_dict)
from utils import load_file_lines

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


def main(config_path):
    # Load experiment config
    config = load_config(config_path)

    # Retrieve paths to inputs
    input_dir = config["data"]["input"]["root"]
    audio_dir = os.path.join(input_dir, config["data"]["input"]["audio_dir"])
    word_outfile = os.path.join(audio_dir, "all_words_4analysis.csv")
    gaze_dir = os.path.join(input_dir, config["data"]["input"]["gaze_dir"])
    gaze_pos_file = os.path.join(gaze_dir, "gaze_positions.csv")
    gaze_before_words_file = os.path.join(gaze_dir, "gaze_positions_before_words.csv")
    gaze_subj_out = os.path.join(gaze_dir, "gaze_positions_4analysis.csv")
    tmp_gaze_s = os.path.join(gaze_dir, "tmp_gaze_positions.csv")
    times_dir = os.path.join(input_dir, config["data"]["input"]["times_dir"])
    surface_dir = os.path.join(input_dir, config["data"]["input"]["surfaces_dir"])
    surface_files = list_subject_files(dir=surface_dir, subject_regex="^", suffix=GAZE_POS_SURFACE_SUFFIX)
    obj_pos_csv = os.path.join(input_dir, config["data"]["input"]["object_positions"])
    
    # Load surface and object position data
    logger.info("Loading surface fixation position data...")
    surface_pos_data = load_and_combine_surface_files(surface_files)
    logger.info("Loading object position data...")
    obj_pos_data = load_object_positions_data(obj_pos_csv)

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
        audio_dir=audio_dir,
        times_dir=times_dir,
        subject_id_regex=subject_id_regex,
    )
    # Get subject IDs (should be identical for all 3 file types)
    subject_ids = sorted(list(subj_audio_erp_dict.keys()))
    logger.info(f"Processing {len(subject_ids)} subject ID(s): {', '.join(subject_ids)}")

    # Iterate over subjects and combine word data with gaze data
    logger.info("Loading word data and combining with gaze data...")
    progress_bar_n = sum([len(subject_files) for _, subject_files in subj_audio_erp_dict.items()])
    with tqdm(total=progress_bar_n) as pbar:
        for subject_id in subject_ids:
            pbar.set_description(f"Processing subject '{subject_id}'...")
            # Initialize empty dataframe to contain all processed gaze data per subject
            gaze_positions_subj = pd.DataFrame()

            # Load word data and combine with gaze data
            audio_erp_files = subj_audio_erp_dict[subject_id]
            times_files = subj_times_dict[subject_id]
            timestamps_files = subj_timestamps_dict[subject_id]
            for erp_file, time_file, timestamp_file in zip(audio_erp_files, times_files, timestamps_files):
                logger.debug(f"ERP file: {os.path.basename(erp_file)}")
                logger.debug(f"Time file: {os.path.basename(time_file)}")
                logger.debug(f"Timestamp file: {os.path.basename(timestamp_file)}")
                # Load ERP file data
                erp_file_data = pd.read_csv(erp_file)
                # Make sure second column is labeled as "time"
                if erp_file_data.columns[1] != WORD_ONSET_FIELD:
                    renamed_columns = erp_file_data.columns.to_list()
                    renamed_columns[1] = WORD_ONSET_FIELD
                    erp_file_data.columns = renamed_columns
                
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
                erp_ids = np.array(erp_file_data[WORD_ID_FIELD])
                erp_time_ids = dict(zip(erp_times, erp_ids))

                # Sort ERP times for binary search
                sorted_erp_times = np.sort(erp_times)
                n_times = len(sorted_erp_times)

                # For each new time, find nearest ERP time associated with a word ID
                word_aligned_times = {}
                for time_i in times:
                    # Find insertion position in sorted ERP time list
                    idx = np.searchsorted(sorted_erp_times, time_i)

                    # Determine neighbors
                    time_before = sorted_erp_times[idx - 1] if idx > 0 else None
                    # NB: Because using insertion index, would be idx rather than idx + 1
                    time_after = sorted_erp_times[idx] if idx < n_times else None

                    # Choose the closer one
                    if time_before is None:
                        word_aligned_times[time_i] = erp_time_ids[time_after]
                    elif time_after is None:
                        word_aligned_times[time_i] = erp_time_ids[time_before]
                    else:
                        before_diff = abs(time_i - time_before)
                        after_diff = abs(time_i - time_after)
                        if before_diff == after_diff:
                            logger.warning("Equal distance to previous and following ERP timestamps, defaulting to previous")
                        word_aligned_times[time_i] = erp_time_ids[time_before if before_diff <= after_diff else time_after]
                
                # Add aligned word IDs to filtered_gaze dataframe
                filtered_gaze[WORD_ID_FIELD] = word_aligned_times.values()

                # Create temp copy of filtered_gaze, renaming "time" to "audio_time"
                tmp_filtered_gaze = filtered_gaze.copy().rename(columns={WORD_ONSET_FIELD: "audio_time"})
                # Merged temp filtered_gaze and erp_file_data by word "id" column
                # Now there should be an "audio_time" column as well as "time" column
                filtered_gaze = tmp_filtered_gaze.merge(erp_file_data, on=WORD_ID_FIELD, how='left')

                # Lowercase text/word field of filtered_gaze
                filtered_gaze[WORD_FIELD] = filtered_gaze[WORD_FIELD].str.lower()

                # Add filtered_gaze to gaze_positions_s dataframe
                gaze_positions_subj = pd.concat([gaze_positions_subj, filtered_gaze], axis=0, ignore_index=True)

                # Update progress bar
                pbar.update(1)

            # Merge gaze positions and surface positions by timestamp
            gaze_positions_subj = gaze_positions_subj.merge(surface_pos_data, on=GAZE_TIMESTAMP_FIELD, how='left')


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Creates trials from the continuous gaze file.")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(args.config)
