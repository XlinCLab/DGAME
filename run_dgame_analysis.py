import argparse
import logging
import os

from dgame.constants import DGAME_DEFAULT_CONFIG, REQUIRED_CONFIG_FIELDS
from dgame.dgame import DGAME
from experiment.experiment_gui import initialize_experiment_from_gui

logger = logging.getLogger(__name__)


def main(config: str = None) -> None:
    # Initialize from GUI based on default config if no config path provided
    if config is None:
        config = initialize_experiment_from_gui(
            config=DGAME_DEFAULT_CONFIG,
            required_fields=REQUIRED_CONFIG_FIELDS,
        )
    # Initialize DGAME experiment from config
    else:
        config = os.path.abspath(config)
    experiment = DGAME(config)

    # Run DGAME analysis
    experiment.run_analysis()

    # Log total experiment duration
    # Write updated experiment config to experiment outdir
    experiment.finish()
    logger.info("Completed successfully.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Runs analysis of DGAME raw data.")
    parser.add_argument('--config', default=None, help='Path to config.yml file')
    args = parser.parse_args()
    main(args.config)
