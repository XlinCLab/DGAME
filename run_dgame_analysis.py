import argparse
import logging
import os

from dgame.dgame import DGAME

logger = logging.getLogger(__name__)

def main(config_path: str) -> dict:
    # Initialize DGAME experiment from config
    experiment = DGAME(config_path)

    # Run DGAME analysis
    experiment.run_analysis()

    # Log total experiment duration
    # Write updated experiment config to experiment outdir
    experiment.finish()
    logger.info("Completed successfully.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Runs analysis of DGAME raw data.")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(os.path.abspath(args.config))
