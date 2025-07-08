import argparse
import logging
import os
import time

import pandas as pd

from constants import BLOCK_IDS, GAZE_TIMESTAMP_FIELD, WORD_ONSET_FIELD
from load_experiment import (create_experiment_outdir, get_experiment_id,
                             load_config, log_step_duration, parse_subject_ids,
                             subject_dirs_dict)

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


def main(config: str | dict) -> dict:
    start_time = time.time()
    # Load experiment config
    if isinstance(config, str):
        config = load_config(config)
    experiment_id = get_experiment_id(config)

    # Output paths
    output_dir = create_experiment_outdir(config, experiment_id)
    fixations_dir = config["data"]["input"]["fixations_dir"]
    fixations_outdir = os.path.join(output_dir, fixations_dir)

    # Get selected subject IDs
    _, subject_id_regex = parse_subject_ids(config["experiment"]["subjects"])
    subj_fixation_dirs_dict = subject_dirs_dict(root_dir=fixations_outdir, subject_regex=subject_id_regex)

    # Process per subject

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
    log_step_duration(config, start_time, step_id="Cc_prepare_fixations_for_matlab")

    return config


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Prepare fixations for inclusion in EEG.event structure in MATLAB")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(args.config)
