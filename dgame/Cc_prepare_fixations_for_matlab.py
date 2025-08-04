import argparse
import logging
import os
import time

import pandas as pd

from dgame.constants import (BLOCK_IDS, GAZE_TIMESTAMP_FIELD, STEP_CC_KEY,
                             WORD_ONSET_FIELD)
from experiment.load_experiment import Experiment
from experiment.test_subjects import subject_dirs_dict

logger = logging.getLogger(__name__)


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


def main(experiment: str | dict | Experiment) -> dict:
    start_time = time.time()

    # Initialize DGAME experiment from config
    if not isinstance(experiment, Experiment):
        from dgame.dgame import DGAME
        experiment = DGAME.from_input(experiment)

    # Process per subject
    subj_fixation_dirs_dict = subject_dirs_dict(
        root_dir=experiment.fixations_outdir,
        subject_regex=experiment.subject_id_regex
    )
    for subject_id, subj_fixation_dir in subj_fixation_dirs_dict.items():
        if len(subj_fixation_dir) > 1:
            logger.warning(f">1 matching directory found for subject ID '{subject_id}'")
        subj_fixation_dir = subj_fixation_dir[0]
        fixation_file = os.path.join(subj_fixation_dir, "fixations_4analysis.csv")
        subj_fixation_data = pd.read_csv(fixation_file)

        # Filter out NA trial_time or fixation_id
        subj_fixation_data = subj_fixation_data[
            subj_fixation_data["trial_time"].notna() &
            subj_fixation_data["fixation_id"].notna()
        ]

        # Rename columns 11, 12, 13 ... to AA, AB, AC ...
        subj_fixation_data = subj_fixation_data.rename(columns=ALPHANUMERIC_COLUMN_MAP)

        # Annotate block and fix_at columns
        subj_fixation_data["block"] = pd.NA
        subj_fixation_data["fix_at"] = pd.NA
        subj_fixation_data.loc[subj_fixation_data["fixation_id"].notna(), "fix_at"] = "elsewhere"
        subj_fixation_data.loc[
            (subj_fixation_data["aoi_target"] == False) &
            (subj_fixation_data[["aoi_comp", "aoi_otherTarget", "aoi_otherComp", "aoi_fillerA", "aoi_fillerB"]] == True).any(axis=1),
            "fix_at"
        ] = "other"
        subj_fixation_data.loc[subj_fixation_data["aoi_goal"] == True, "fix_at"] = "goal"

        pattern_id, set_id = 1, 1
        for idx in range(len(subj_fixation_data) - 1):
            current_time = subj_fixation_data.iloc[idx][WORD_ONSET_FIELD]
            next_time = subj_fixation_data.iloc[idx + 1][WORD_ONSET_FIELD]

            if pattern_id > 2:
                pattern_id = 1
                set_id += 1

            if current_time < next_time and pattern_id < 3:
                subj_fixation_data.iloc[idx, subj_fixation_data.columns.get_loc("block")] = f"{set_id}{pattern_id}"
            elif current_time > next_time and pattern_id < 3:
                subj_fixation_data.iloc[idx, subj_fixation_data.columns.get_loc("block")] = f"{set_id}{pattern_id}"
                pattern_id += 1

        # Sort by gaze_timestamp field
        subj_fixation_data = subj_fixation_data.sort_values(by=GAZE_TIMESTAMP_FIELD)

        # Iterate over blocks and write filtered block data to CSV
        for block in BLOCK_IDS:
            block_data = subj_fixation_data[subj_fixation_data["block"] == str(block)]
            fixation_outfile = os.path.join(subj_fixation_dir, f"fixations_times_{block}_trials.csv")
            block_data.to_csv(fixation_outfile, index=False)
        logger.info(f"Wrote subject <{subject_id}> block fixation files to {subj_fixation_dir}")

    # Log duration of this step in run config
    experiment.log_step_duration(start_time, step_id=STEP_CC_KEY)

    return experiment


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Prepare fixations for inclusion in EEG.event structure in MATLAB")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(os.path.abspath(args.config))
