import argparse
import logging
import os
import time

import pandas as pd

from constants import AOI_COLUMNS, FIXATION_TIMES_TRIALS_SUFFIX
from load_experiment import (create_experiment_outdir, get_experiment_id,
                             list_subject_files, load_config,
                             log_step_duration, parse_subject_ids,
                             subject_dirs_dict)

logger = logging.getLogger(__name__)


def load_fixation_times_trials_files(subj_fixation_dirs_dict: dict) -> pd.DataFrame:
    """Loads fixation times trials files from selected subjects into a single dataframe."""
    fixation_times_trials_df = pd.DataFrame()
    for subject_id, subj_fixation_dirs in subj_fixation_dirs_dict.items():
        if len(subj_fixation_dirs) > 1:
            logger.warning(f">1 matching directory found for subject ID '{subject_id}'")
        subj_fixation_dir = subj_fixation_dirs[0]
        fixation_times_trials_files = list_subject_files(
            dir=subj_fixation_dir,
            subject_regex=r"^",
            suffix=FIXATION_TIMES_TRIALS_SUFFIX
        )
        for fixation_time_trial_file in fixation_times_trials_files:
            data = pd.read_csv(fixation_time_trial_file)
            fixation_times_trials_df = pd.concat([fixation_times_trials_df, data], axis=0, ignore_index=True)
    return fixation_times_trials_df


def determine_fixation_label(row: pd.Series) -> str:
    """Returns a fixation label for a given dataframe row given AOI surface values."""
    for aoi_col in AOI_COLUMNS:
        # NB: need to iterate over AOI columns in the specific order listed in AOI_COLUMNS dict
        if row[aoi_col] is True:
            # Simply remove the "aoi_" prefix, e.g. "aoi_comp" -> "comp"
            return aoi_col.replace("aoi_", "")
    # Elsewhere case for when fixation was not in any AOI
    return "elsewhere"


def main(config: str | dict) -> dict:
    start_time = time.time()
    # Load experiment config
    if isinstance(config, str):
        config = load_config(config)
    experiment_id = get_experiment_id(config)

    # Get input/output paths
    output_dir = create_experiment_outdir(config, experiment_id)
    fixations_dir = config["data"]["input"]["fixations_dir"]
    fixations_outdir = os.path.join(output_dir, fixations_dir)

    # Get selected subject IDs and per-subject fixation outdirs 
    _, subject_id_regex = parse_subject_ids(config["experiment"]["subjects"])
    subj_fixation_dirs_dict = subject_dirs_dict(root_dir=fixations_outdir, subject_regex=subject_id_regex)
    subject_ids = subj_fixation_dirs_dict.keys()

    # Iterate through subject IDs, retrieve relevant fixation_times_*_trials.csv files, and combine into single dataframe
    # (!) NB: these fixation_times_*_trials.csv are written in Step Cc, which is currently not reproduced
    logger.warning("********** (!) Using preset fixation_times_*_trials.csv files as Step Cc is not yet implemented **********")
    fixation_times_trials_df = load_fixation_times_trials_files(subj_fixation_dirs_dict)

    # Add fixation label to "fix_at" column
    fixation_times_trials_df["fix_at"] = fixation_times_trials_df.apply(determine_fixation_label, axis=1)
    # Convert to category type (equivalent to R's factor)
    fixation_times_trials_df["fix_at"] = fixation_times_trials_df["fix_at"].astype("category")

    # Iterate through subject IDs and get per-subject fixation data
    fixation_times_trials_df["subj"] = fixation_times_trials_df["subj"].astype(str)
    logger.warning("********** (!) Using hard-coded subject ID = 2 from preset input files **********")
    subject_ids = ['2'] # TODO remove this line once hard-coded/preset files are no longer used
    for subject_id in subject_ids:
        subj_fixation_data = fixation_times_trials_df[fixation_times_trials_df["subj"] == subject_id]


    # Log duration of this step in run config
    log_step_duration(config, start_time, step_id="Db_plot_descriptive_fixation")

    return config


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Compute and plot descriptive fixation statistics.")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(args.config)
