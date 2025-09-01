import argparse
import logging
import os

from dgame.constants import SCRIPT_DIR, STEP_G_KEY
from experiment.load_experiment import Experiment
from experiment.test_subjects import subject_dirs_dict
from utils.matlab_interface import run_matlab_script

logger = logging.getLogger(__name__)


def main(experiment: str | dict | Experiment) -> Experiment:
    # Initialize DGAME experiment from config
    if not isinstance(experiment, Experiment):
        from dgame.dgame import DGAME
        experiment = DGAME.from_input(experiment)

    # Get list of subject IDs and their corresponding XDF directory paths
    subject_eeg_dirs_dict = subject_dirs_dict(
        root_dir=experiment.eeg_outdir,
        subject_regex=experiment.subject_id_regex,
    )
    subject_ids = list(subject_eeg_dirs_dict.keys())
    subject_eeg_dirs = list(subject_eeg_dirs_dict.values())
    assert all(len(subj_eeg_dir) == 1 for subj_eeg_dir in subject_eeg_dirs)
    subject_eeg_dirs = [subj_eeg_dir[0] for subj_eeg_dir in subject_eeg_dirs]

    # Designate MATLAB logfile
    matlab_logdir = os.path.join(experiment.logdir, "MATLAB")
    os.makedirs(matlab_logdir, exist_ok=True)
    logfile = os.path.join(matlab_logdir, f"{STEP_G_KEY}.log")

    # Run preproc_EEG step in MATLAB
    run_matlab_script(
        os.path.join(SCRIPT_DIR, f"{STEP_G_KEY}.m"),
        args=[
            subject_ids,
            subject_eeg_dirs,
            experiment.input_dir,
            experiment.matlab_root,
        ],
        matlab_version="R2021a",  # R2021a required for MoBILAB dependency
        logfile=logfile,
    )

    return experiment


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Preprocess EEG data.")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(os.path.abspath(args.config))
