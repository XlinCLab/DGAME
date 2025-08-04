import argparse
import logging
import os
import time

import matlab.engine

from dgame.constants import STEP_A_KEY
from experiment.load_experiment import Experiment
from experiment.test_subjects import subject_dirs_dict

logger = logging.getLogger(__name__)

def main(experiment: str | dict | Experiment):

    start_time = time.time()

    # Initialize DGAME experiment from config
    if not isinstance(experiment, Experiment):
        from dgame.dgame import DGAME
        experiment = DGAME.from_input(experiment)

    # Get list of subject IDs
    subject_xdf_dirs = subject_dirs_dict(
        root_dir=experiment.xdf_indir,
        subject_regex=experiment.subject_id_regex,
    )

    logger.info("Starting MATLAB engine...")
    eng = matlab.engine.start_matlab()
    # Add path to directory with MATLAB script (same as this directory) to MATLAB engine path
    eng.addpath(os.path.dirname(os.path.abspath(__file__)))
    eng.A_export_audio_and_et_times(
        subject_xdf_dirs, experiment.input_dir, experiment.matlab_root,
        nargout=0
    )
    eng.quit()
    logger.info("MATLAB engine ended.")

    # Log duration of this step in run config
    experiment.log_step_duration(start_time, step_id=STEP_A_KEY)

    return experiment


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Export audio data and times of the gaze data.")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(os.path.abspath(args.config))
