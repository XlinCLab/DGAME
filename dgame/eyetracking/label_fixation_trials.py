import argparse
import os

import pandas as pd

from dgame.constants import BLOCK_IDS
from dgame.eyetracking import COLUMN_DATA_TYPES, GAZE_TIMESTAMP_FIELD
from experiment.input_validation import (OutputValidationError,
                                         assert_output_file_exists)
from experiment.load_experiment import Experiment

ALPHANUMERIC_COLUMN_MAP = {
    '11': 'AA',
    '12': 'AB',
    '13': 'AC',
    '14': 'AD',
    '21': 'BA',
    '22': 'BB',
    '23': 'BC',
    '24': 'BD',
    '31': 'CA',
    '32': 'CB',
    '33': 'CC',
    '34': 'CD',
    '41': 'DA',
    '42': 'DB',
    '43': 'DC',
    '44': 'DD',
}


def validate_outputs(experiment, subject_ids: list) -> None:
    """Validate fixation-trial output files."""

    # Verify fixation directories exist per subject and that expected files were created
    subject_fixation_dirs_dict = experiment.get_subject_dirs_dict(experiment.fixations_outdir)

    # Make sure fixation directories were found for exactly the expected subjects
    fixation_subject_ids = sorted(list(subject_fixation_dirs_dict.keys()))
    if fixation_subject_ids != subject_ids:
        missing_fixations = [subject_id for subject_id in subject_ids if subject_id not in subject_fixation_dirs_dict]
        if len(missing_fixations) > 0:
            raise OutputValidationError(f"Fixations directory missing for following subjects: {', '.join(missing_fixations)}")

    # Verify that expected per-block fixation-trial files were created
    for subject_id, subject_fixation_dirs in subject_fixation_dirs_dict.items():
        # Verify that there is only one fixations directory per subject
        try:
            assert len(subject_fixation_dirs) == 1
        except AssertionError as exc:
            raise OutputValidationError(f">1 fixations directory found for subject <{subject_id}>") from exc
        subject_fixation_dir = subject_fixation_dirs[0]

        for block in BLOCK_IDS:
            fixation_trials_file = os.path.join(subject_fixation_dir, f"fixations_times_{block}_trials.csv")
            assert_output_file_exists(fixation_trials_file)


def main(experiment: str | dict | Experiment) -> Experiment:
    # Initialize DGAME experiment from config
    if not isinstance(experiment, Experiment):
        from dgame.dgame import DGAME
        experiment = DGAME.from_input(experiment)
    logger = experiment.logger

    # Process per subject
    subj_fixation_dirs_dict = experiment.get_subject_dirs_dict(experiment.fixations_outdir)
    for subject_id, subj_fixation_dir in subj_fixation_dirs_dict.items():
        if len(subj_fixation_dir) > 1:
            logger.warning(f">1 matching directory found for subject ID '{subject_id}'")
        subj_fixation_dir = subj_fixation_dir[0]
        fixation_file = os.path.join(subj_fixation_dir, "fixations_4analysis.csv")
        # Ensure subject IDs are read in as strings
        subj_fixation_data = pd.read_csv(fixation_file, dtype=COLUMN_DATA_TYPES)

        # Filter out NA trial_time or fixation_id
        subj_fixation_data = subj_fixation_data[
            subj_fixation_data["trial_time"].notna() &
            subj_fixation_data["fixation_id"].notna()
        ]

        # Rename columns 11, 12, 13 ... to AA, AB, AC ...
        subj_fixation_data = subj_fixation_data.rename(columns=ALPHANUMERIC_COLUMN_MAP)

        # Annotate fix_at column
        subj_fixation_data["fix_at"] = pd.NA
        subj_fixation_data.loc[subj_fixation_data["fixation_id"].notna(), "fix_at"] = "elsewhere"
        subj_fixation_data.loc[subj_fixation_data["aoi_target"] == True, "fix_at"] = "target"
        subj_fixation_data.loc[
            (subj_fixation_data["aoi_target"] == False) &  # noqa: E712
            (subj_fixation_data[["aoi_comp", "aoi_otherTarget", "aoi_otherComp", "aoi_fillerA", "aoi_fillerB"]] == True).any(axis=1),  # noqa: E712
            "fix_at"
        ] = "other"
        subj_fixation_data.loc[subj_fixation_data["aoi_goal"] == True, "fix_at"] = "goal"   # noqa: E712

        # Sort by gaze_timestamp field
        subj_fixation_data = subj_fixation_data.sort_values(by=GAZE_TIMESTAMP_FIELD)

        # Iterate over blocks and write filtered block data to CSV
        for block in BLOCK_IDS:
            block_data = subj_fixation_data[subj_fixation_data["block"] == block]
            fixation_outfile = os.path.join(subj_fixation_dir, f"fixations_times_{block}_trials.csv")
            block_data.to_csv(fixation_outfile, index=False)
        logger.info(f"Wrote subject <{subject_id}> block fixation files to {subj_fixation_dir}")

    # Validate outputs
    validate_outputs(experiment, experiment.subject_ids)

    return experiment


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Prepare fixations for inclusion in the EEG event structure")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(os.path.abspath(args.config))
