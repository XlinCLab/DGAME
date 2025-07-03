import argparse
import logging
import os
import time

from load_experiment import (create_experiment_outdir, get_experiment_id,
                             load_config, log_step_duration, parse_subject_ids,
                             subject_files_dict)

logger = logging.getLogger(__name__)

def main(config: str | dict) -> dict:
    start_time = time.time()
    # Load experiment config
    if isinstance(config, str):
        config = load_config(config)
    experiment_id = get_experiment_id(config)


    # Log duration of this step in run config
    log_step_duration(config, start_time, step_id="Db_plot_descriptive_fixation")

    return config


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Compute and plot descriptive fixation statistics.")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(args.config)
