import argparse
import os

from experiment.load_experiment import Experiment


def main(experiment: str | dict | Experiment) -> Experiment:

    # Initialize DGAME experiment from config
    if not isinstance(experiment, Experiment):
        from dgame.dgame import DGAME
        experiment = DGAME.from_input(experiment)
    logger = experiment.logger




    return experiment


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Extract Director Task instructions from word-based CSVs and derive abstract syntactic patterns.")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(os.path.abspath(args.config))
