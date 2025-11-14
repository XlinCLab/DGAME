import argparse
import logging
import os

from dgame.constants import STEP_F_KEY
from experiment.load_experiment import Experiment

logger = logging.getLogger(__name__)


def main(experiment: str | dict | Experiment) -> Experiment:
    # Initialize DGAME experiment from config
    if not isinstance(experiment, Experiment):
        from dgame.dgame import DGAME
        experiment = DGAME.from_input(experiment)

    # Get flattened list of per-subject xdf directories as input to MATLAB script
    subject_xdf_dirs = [
        os.path.join(experiment.xdf_indir, subject_id)
        for subject_id in experiment.subject_ids
    ]

    # Create per-subject EEG and ICA outdirs
    for subj_id in experiment.subject_ids:
        os.makedirs(os.path.join(experiment.eeg_outdir, subj_id), exist_ok=True)
        os.makedirs(os.path.join(experiment.eeg_ica_outdir, subj_id), exist_ok=True)

    # Retrieve list of electrodes which were removed to fit eyetracking glasses
    removed_electrodes = experiment.get_dgame_step_parameter(STEP_F_KEY, "removed_electrodes")

    # Run preproc_EEG step in MATLAB
    experiment.run_matlab_step(
        os.path.join(experiment.matlab_script_dir, f"{STEP_F_KEY}.m"),
        args=[
            experiment.subject_ids,
            subject_xdf_dirs,
            experiment.input_dir,
            experiment.outdir,
            experiment.eeg_ica_outdir,
            experiment.matlab_root,
            experiment.dgame_version,
            removed_electrodes,
        ]
    )

    return experiment


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Preprocess EEG data.")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(os.path.abspath(args.config))
