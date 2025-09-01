import argparse
import logging
import os

from dgame.constants import SCRIPT_DIR, STEP_G_KEY
from experiment.load_experiment import Experiment
from experiment.test_subjects import subject_dirs_dict

logger = logging.getLogger(__name__)


def main(experiment: str | dict | Experiment) -> Experiment:
    # Initialize DGAME experiment from config
    if not isinstance(experiment, Experiment):
        from dgame.dgame import DGAME
        experiment = DGAME.from_input(experiment)

    # Get list of subject IDs and their corresponding EEG directory paths
    subject_eeg_dirs_dict = subject_dirs_dict(
        root_dir=experiment.eeg_outdir,
        subject_regex=experiment.subject_id_regex,
    )
    subject_ids = list(subject_eeg_dirs_dict.keys())
    subject_eeg_dirs = list(subject_eeg_dirs_dict.values())
    assert all(len(subj_eeg_dir) == 1 for subj_eeg_dir in subject_eeg_dirs)
    subject_eeg_dirs = [subj_eeg_dir[0] for subj_eeg_dir in subject_eeg_dirs]

    # Run deconvolution_ERPs step in MATLAB
    experiment.run_matlab_step(
        os.path.join(SCRIPT_DIR, f"{STEP_G_KEY}.m"),
        args=[
            subject_ids,
            subject_eeg_dirs,
            experiment.matlab_root,
        ]
    )

    return experiment


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Compute rERPs.")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(os.path.abspath(args.config))
