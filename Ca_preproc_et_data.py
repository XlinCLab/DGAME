import argparse
import glob
import logging
import os

import numpy as np
import pandas as pd

from constants import (AUDIO_ERP_FILE_SUFFIX, GAZE_POS_SURFACE_SUFFIX, ROUND_N,
                       TIMES_FILE_SUFFIX, TIMESTAMPS_FILE_SUFFIX, WORD_FIELD,
                       WORD_ID_FIELD, WORD_ONSET_FIELD)
from load_experiment import (list_subject_files, load_config,
                             load_object_positions_data, parse_subject_ids)
from utils import load_file_lines

logger = logging.getLogger(__name__)

def main(config_path):
    # Load experiment config
    config = load_config(config_path)

    # Get selected subject IDs
    _, subject_id_regex = parse_subject_ids(config["experiment"]["subjects"])

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
    obj_pos_data = load_object_positions_data(obj_pos_csv)

    # Find subject audio ERP files
    audio_erp_files = list_subject_files(dir=audio_dir, subject_regex=subject_id_regex, suffix=AUDIO_ERP_FILE_SUFFIX)
    # Find subject times and timestamps files
    times_files = list_subject_files(dir=times_dir, subject_regex=subject_id_regex, suffix=TIMES_FILE_SUFFIX)
    timestamps_files = list_subject_files(dir=times_dir, subject_regex=subject_id_regex, suffix=TIMESTAMPS_FILE_SUFFIX)

    # Load gaze file (columns of interest only)
    logger.info(f"Loading gaze data from {gaze_pos_file}")
    raw_gaze_data = pd.read_csv(
        gaze_pos_file,
        usecols=[
            "gaze_timestamp",
            "world_index",
            "confidence",
            "norm_pos_x",
            "norm_pos_y",
            "base_data",
        ]
    )
    # Round gaze_timestamp field of raw_gaze_data to ROUND_N places
    raw_gaze_data['gaze_timestamp'] = raw_gaze_data['gaze_timestamp'].astype(float).round(ROUND_N)

    # Load word data and combine with gaze data
    logger.info("Loading word data and combining with gaze data...")
    for erp_file, time_file, timestamp_file in zip(audio_erp_files, times_files, timestamps_files):
        logger.info(f"ERP file: {os.path.basename(erp_file)}")
        logger.info(f"Time file: {os.path.basename(time_file)}")
        logger.info(f"Timestamp file: {os.path.basename(timestamp_file)}")
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
            (raw_gaze_data["gaze_timestamp"] >= start_timestamp) &
            (raw_gaze_data["gaze_timestamp"] < end_timestamp)
        ].copy()
        # Add times array as new column "audio_time" to filtered_gaze dataframe
        filtered_gaze["audio_time"] = times

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

        # Merged filtered_gaze and erp_file_data by word "id" column
        filtered_gaze = filtered_gaze.merge(erp_file_data, on=WORD_ID_FIELD, how='left')

        # Lowercase text/word field of filtered_gaze
        filtered_gaze[WORD_FIELD] = filtered_gaze[WORD_FIELD].str.lower()


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Creates trials from the continuous gaze file.")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(args.config)
