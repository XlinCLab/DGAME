import argparse
import logging
import os
import time
from datetime import timedelta

import pandas as pd
import rpy2.robjects as ro

from constants import RUN_CONFIG_KEY, WORD_END_FIELD, WORD_ONSET_FIELD
from load_experiment import (create_experiment_outdir, get_experiment_id,
                             load_config, parse_subject_ids, subject_dirs_dict)

logger = logging.getLogger(__name__)


def main(config: str | dict) -> dict:
    start_time = time.time()
    # Load experiment config
    if isinstance(config, str):
        config = load_config(config)
    experiment_id = get_experiment_id(config)

    # Get input/output paths
    output_dir = create_experiment_outdir(config, experiment_id)
    gaze_dir = config["data"]["input"]["gaze_dir"]
    gaze_outdir = os.path.join(output_dir, gaze_dir)

    # Get selected subject IDs and per-subject output directories
    _, subject_id_regex = parse_subject_ids(config["experiment"]["subjects"])
    subject_gaze_dirs_dict = subject_dirs_dict(root_dir=gaze_outdir, subject_regex=subject_id_regex)
    subject_ids = sorted(list(subject_gaze_dirs_dict.keys()))
    logger.info(f"Processing {len(subject_ids)} subject ID(s): {', '.join(subject_ids)}")

    # Iterate over subject directories
    for subject_id, subject_gaze_dirs in subject_gaze_dirs_dict.items():
        logger.info(f"Processing subject '{subject_id}'...")

        # Get per-subject gaze and fixation directories/files
        if len(subject_gaze_dirs) > 1:
            logger.warning(f">1 matching directory found for subject ID '{subject_id}'")
        subject_gaze_dir = subject_gaze_dirs[0]
        # TODO need to confirm whether this is meant to be a single file per subject or a single file with all subjects' data; current implementation assumes the former
        gaze_infile = os.path.join(subject_gaze_dir, "gaze_positions_4analysis.csv")

        # Load gaze_infile and drop all non-trial data points
        gaze_positions_all = pd.read_csv(gaze_infile)
        gaze_positions_all = gaze_positions_all[gaze_positions_all["condition"].notna()]
        gaze2analysis = (
            gaze_positions_all
            .loc[
                gaze_positions_all["trial_time"].notna() &
                (gaze_positions_all["trackloss"] is False)
            ]
            .drop_duplicates()
            .assign(duration=lambda df: df[WORD_END_FIELD] - df[WORD_ONSET_FIELD])
        )

        # Filter valid durations and then drop duration column
        gaze2analysis = gaze2analysis[
            gaze2analysis["duration"] >= 0
        ].drop(columns="duration")

    # Calculate duration of this step and add to run config
    end_time = time.time()
    duration = str(timedelta(seconds=int(end_time - start_time)))
    config[RUN_CONFIG_KEY]["duration"]["Da_gaze_stats"] = duration
    logger.info(f"Step Da completed successfully (duration: {duration}).")
    return config


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Process fixation files from pupil player and prepare them for further processing.")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(args.config)
