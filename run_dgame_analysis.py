import argparse
import logging
import os

from dgame.config import DGAME_DEFAULT_CONFIG, REQUIRED_CONFIG_FIELDS
from dgame.dgame import DGAME
from dgame.pipeline import FULL_DGAME_PIPELINE, STEP_REQUIREMENTS
from experiment.gui.main_window import run_config_gui

logger = logging.getLogger(__name__)


def main(config: str = None) -> None:
    # Open the experiment setup GUI, based on the default config, if no config path is provided
    # (the GUI runs the analysis itself, from within its window)
    if config is None:
        run_config_gui(
            config=DGAME_DEFAULT_CONFIG,
            required_fields=REQUIRED_CONFIG_FIELDS,
            pipeline_steps=FULL_DGAME_PIPELINE,
            step_requirements=STEP_REQUIREMENTS,
        )
        return
    # Initialize DGAME experiment from config
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
