import argparse
import logging
import os

from dgame.constants import BLOCK_IDS, STEP_A_KEY
from experiment.input_validation import (OutputValidationError,
                                         assert_output_file_exists)
from experiment.load_experiment import Experiment

logger = logging.getLogger(__name__)


def validate_outputs(experiment, subject_ids: list) -> None:
    """Validate audio outputs from MATLAB script."""
    from dgame.dgame import validate_dgame_input
    experiment = validate_dgame_input(experiment)

    # Verify audio directories exist per subject and that expected files were created
    subject_audio_dirs_dict = experiment.get_subject_dirs_dict(experiment.audio_indir)

    # Make sure audio and xdf directories are found for exactly the same subjects
    audio_subject_ids = sorted(list(subject_audio_dirs_dict.keys()))
    if audio_subject_ids != subject_ids:
        missing_audio = [subject_id for subject_id in subject_ids if subject_id not in subject_audio_dirs_dict]
        if len(missing_audio) > 0:
            raise OutputValidationError(f"Audio directory missing for following subjects: {', '.join(missing_audio)}")

    # Verify that expected audio and times files were created
    for subject_id, subject_audio_dirs in subject_audio_dirs_dict.items():
        # Verify that there is only one audio directory per subject
        try:
            assert len(subject_audio_dirs) == 1
        except AssertionError as exc:
            raise OutputValidationError(f">1 audio directory found for subject <{subject_id}>") from exc
        subject_audio_dir = subject_audio_dirs[0]

        # Verify individual audio and time files
        subj_times_dir = os.path.join(experiment.times_dir, subject_id)
        for block in BLOCK_IDS:
            for condition_label in {"decke", "director"}:  # TODO could save these as a constant somewhere, since also referenced in MATLAB script
                audio_file = os.path.join(subject_audio_dir, f"{subject_id}_{condition_label}_{block}.wav")
                assert_output_file_exists(audio_file)

                # time files per subject per block
                timestamp_file = os.path.join(subj_times_dir, f"{subject_id}_timestamps_max-min_{block}.csv")
                times_file = os.path.join(subj_times_dir, f"{subject_id}_times_{block}.csv")
                assert_output_file_exists(timestamp_file)
                assert_output_file_exists(times_file)

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

    # Run export_audio_and_et_times step in MATLAB
    experiment.run_matlab_step(
        os.path.join(experiment.matlab_script_dir, f"{STEP_A_KEY}.m"),
        args=[
            experiment.subject_ids,
            subject_xdf_dirs,
            experiment.input_dir,
            experiment.times_dir,
            experiment.matlab_root,
            experiment.dgame_version,
        ]
    )

    # Validate outputs from MATLAB
    validate_outputs(experiment, experiment.subject_ids)

    return experiment


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Export audio data and times of the gaze data.")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(os.path.abspath(args.config))
