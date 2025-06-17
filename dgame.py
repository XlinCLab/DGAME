import argparse
import json
import logging
import os
import time
from datetime import timedelta

from B_prepare_words import main as step_b
from Ca_preproc_et_data import main as step_ca
from constants import RUN_CONFIG_KEY
from load_experiment import (create_experiment_outdir, dump_config,
                             get_experiment_id, load_config)

logger = logging.getLogger(__name__)


def run_dgame_analysis(config: dict) -> dict:
    """Run all component DGAME analysis steps."""
    # Run Step B: prepare words data
    config = step_b(config)
    # Run Step Ca: preproc ET data
    config = step_ca(config)
    return config


def main(config_path: str) -> dict:
    start_time = time.time()

    # Load config
    config = load_config(config_path)
    logger.info(json.dumps(config, indent=4))

    # Retrieve and set experiment ID and outdir
    experiment_id = get_experiment_id(config)
    experiment_outdir = create_experiment_outdir(config, experiment_id)
    config[RUN_CONFIG_KEY]["id"] = experiment_id
    config[RUN_CONFIG_KEY]["outdir"] = experiment_outdir

    # Run DGAME analysis
    config = run_dgame_analysis(config)

    # Add total duration to config output
    end_time = time.time()
    total_duration = str(timedelta(seconds=int(end_time - start_time)))
    config[RUN_CONFIG_KEY]["duration"]["total"] = total_duration

    # Write updated experiment config to experiment outdir
    config_outpath = os.path.join(experiment_outdir, "config.yml")
    dump_config(config, config_outpath)
    logger.info(f"Wrote experiment run config to {os.path.abspath(config_outpath)}")
    logger.info("Completed successfully.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Runs analysis of DGAME raw data.")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(args.config)
