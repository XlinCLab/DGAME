import argparse
import logging
import os

from dgame.constants import STEP_F_KEY
from experiment.load_experiment import Experiment
from experiment.test_subjects import subject_dirs_dict

logger = logging.getLogger(__name__)


def main(experiment: str | dict | Experiment) -> Experiment:
    # Initialize DGAME experiment from config
    if not isinstance(experiment, Experiment):
        from dgame.dgame import DGAME
        experiment = DGAME.from_input(experiment)

    # Get list of subject IDs and their corresponding XDF directory paths
    subject_xdf_dirs_dict = subject_dirs_dict(
        root_dir=experiment.xdf_indir,
        subject_regex=experiment.subject_id_regex,
    )
    subject_ids = list(subject_xdf_dirs_dict.keys())
    subject_xdf_dirs = list(subject_xdf_dirs_dict.values())
    assert all(len(subj_xdf_dir) == 1 for subj_xdf_dir in subject_xdf_dirs)
    subject_xdf_dirs = [subj_xdf_dir[0] for subj_xdf_dir in subject_xdf_dirs]

    # Run preproc_EEG step in MATLAB
    experiment.run_matlab_step(
        os.path.join(experiment.matlab_script_dir, f"{STEP_F_KEY}.m"),
        args=[
            subject_ids,
            subject_xdf_dirs,
            experiment.input_dir,
            experiment.matlab_root,
        ]
    )

    return experiment


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Preprocess EEG data.")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(os.path.abspath(args.config))
