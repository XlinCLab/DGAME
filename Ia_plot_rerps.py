import argparse
import os
import time

import pandas as pd

from constants import CHANNEL_COORDS_FILE, ERP_NOUN_FILE_SUFFIX, ERP_FIXATION_FILE_SUFFIX
from load_experiment import (create_experiment_outdir, get_experiment_id,
                             load_config, log_step_duration, parse_subject_ids,
                             subject_files_dict)
from utils import load_csv_list


def main(config: str | dict) -> dict:
    start_time = time.time()
    # Load experiment config
    if isinstance(config, str):
        config = load_config(config)
    experiment_id = get_experiment_id(config)

    # Get input/output paths
    input_dir = config["data"]["input"]["root"]
    eeg_dir = os.path.join(input_dir, config["data"]["input"]["eeg_dir"])
    output_dir = create_experiment_outdir(config, experiment_id)
    
    # Load channel coords file
    channel_coords_file = os.path.join(input_dir, CHANNEL_COORDS_FILE)  # TODO is this a static file across all subjects? where do we expect it to be found? eeg_dir?
    channel_coords = pd.read_csv(channel_coords_file, names=["channel", "lat", "sag", "z"])
    channel_coords["channel"] = channel_coords["channel"].astype(str)

    # Get selected subject IDs and per-subject files
    _, subject_id_regex = parse_subject_ids(config["experiment"]["subjects"])
    subject_erp_noun_files = subject_files_dict(
        dir=os.path.join(eeg_dir, "unfold_out", "results"),
        subject_regex=subject_id_regex,
        suffix=ERP_NOUN_FILE_SUFFIX,
    )
    subject_erp_fixation_files = subject_files_dict(
        dir=os.path.join(eeg_dir, "unfold_out", "results"),
        subject_regex=subject_id_regex,
        suffix=ERP_FIXATION_FILE_SUFFIX,
    )
    subject_ids = noun_subject_ids = sorted(list(subject_erp_noun_files.keys()))
    fixation_subject_ids = sorted(list(subject_erp_fixation_files.keys()))
    # Make sure the same subject IDs are found for both types
    assert fixation_subject_ids == noun_subject_ids

    # Iterate over subjects
    for subject_id in subject_ids:
        erp_noun_files = subject_erp_noun_files[subject_id]
        erp_fixation_files = subject_erp_fixation_files[subject_id]

        # Load noun data
        erp_noun_data = load_csv_list(erp_noun_files)
        erp_noun_data["channel"] = erp_noun_data["channel"].astype(str)

        # Load fixation data
        erp_fixation_data = load_csv_list(erp_fixation_files)
        erp_fixation_data["channel"] = erp_fixation_data["channel"].astype(str)

        # Merge channel coordinates into both dataframes
        erp_noun_data = erp_noun_data.merge(channel_coords, how="left", on="channel")
        erp_fixation_data = erp_fixation_data.merge(channel_coords, how="left", on="channel")

    # Log duration of this step in run config
    log_step_duration(config, start_time, step_id="Da_gaze_stats")

    return config


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Calculate and plot gaze statistics.")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(args.config)
