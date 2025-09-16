import argparse
import logging
import os

from dgame.constants import STEP_JA_KEY
from dgame.J_lm_permute_and_plot_fixations_and_language import step_J_analysis
from experiment.load_experiment import Experiment

logger = logging.getLogger(__name__)
        

def main(experiment: str | dict | Experiment) -> Experiment:
    # Initialize DGAME experiment from config
    if not isinstance(experiment, Experiment):
        from dgame.dgame import DGAME
        experiment = DGAME.from_input(experiment)

    n_permutations = experiment.get_dgame_step_parameter(STEP_JA_KEY, "n_permutations")
    include_baseline = experiment.get_dgame_step_parameter(STEP_JA_KEY, "include_baseline")
    experiment = step_J_analysis(
        experiment,
        mode="FIX",
        n_permutations=n_permutations,
        include_baseline=include_baseline,
    )

    return experiment

if __name__ == "__main__":
    parser = argparse.ArgumentParser("Run permutation tests and plot fixation statistics.")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(os.path.abspath(args.config))
