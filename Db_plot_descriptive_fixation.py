import argparse
import os
import time
from datetime import timedelta
import logging

from constants import RUN_CONFIG_KEY
from load_experiment import (create_experiment_outdir, get_experiment_id,
                             load_config, parse_subject_ids,
                             subject_files_dict)

logger = logging.getLogger(__name__)

def main(config: str | dict) -> dict:
    start_time = time.time()
    # Load experiment config
    if isinstance(config, str):
        config = load_config(config)
    experiment_id = get_experiment_id(config)


    # Calculate duration of this step and add to run config
    end_time = time.time()
    duration = str(timedelta(seconds=int(end_time - start_time)))
    config[RUN_CONFIG_KEY]["duration"]["Da_gaze_stats"] = duration
    logger.info(f"Step Da completed successfully (duration: {duration}).")
    return config


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Compute and plot descriptive fixation statistics.")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(args.config)
