import argparse
import os
import re

import pandas as pd

from dgame.constants import (AUDIO_ERP_FILE_SUFFIX, DET_POS_LABEL,
                             NOUN_POS_LABEL, PART_OF_SPEECH_FIELD, PATTERN_IDS,
                             SET_IDS)
from dgame.eyetracking.utils import load_filtered_gaze_data
from experiment.load_experiment import Experiment
from utils.utils import list_matching_files


def main(experiment: str | dict | Experiment) -> Experiment:
    # Initialize DGAME experiment from config
    if not isinstance(experiment, Experiment):
        from dgame.dgame import DGAME
        experiment = DGAME.from_input(experiment)
    logger = experiment.logger

    # Get selected subject IDs and directories
    subject_audio_dirs = experiment.get_subject_dirs_dict(experiment.preproc_audio_indir)
    logger.info(f"Processing {len(experiment.subject_ids)} subject ID(s): {', '.join(experiment.subject_ids)}")

    # Load and filter gaze data (output of step Ca)
    gaze2analysis = load_filtered_gaze_data(experiment)

    # Iterate over subject directories
    for subject_id, subject_audio_dir in subject_audio_dirs.items():
        logger.info(f"Processing subject '{subject_id}'...")

        # Get per-subject audio ERP files
        if len(subject_audio_dir) > 1:
            logger.warning(f">1 audio directory found for subject {subject_id}")
        subject_audio_dir = subject_audio_dir[0]
        subj_audio_erp_files = list_matching_files(
            dir=subject_audio_dir,
            pattern=AUDIO_ERP_FILE_SUFFIX
        )

        # Designate and create per-subject audio outdir (if doesn't already exist)
        subj_audio_outdir = os.path.join(experiment.audio_outdir, subject_id)
        os.makedirs(subj_audio_outdir, exist_ok=True)

        # Iterate through subject ERP audio files
        trial_counter_nouns, trial_counter_determiners = 1, 1
        for word_infile in sorted(subj_audio_erp_files):
            # Designate name for outfile
            basename, _ = os.path.splitext(os.path.basename(word_infile))
            word_outfile = os.path.join(subj_audio_outdir, f"{basename}_trialtime.csv")

            # Parse pattern and set IDs
            # e.g. 03_words2erp_12 -> set_id = 1, pattern_id = 2
            block_id = re.search(r"(?<=_)(\d+)$", basename).group(1)
            set_id, pattern_id = map(int, list(block_id))
            assert set_id in SET_IDS
            assert pattern_id in PATTERN_IDS

            # Read infile data
            word_data = pd.read_csv(word_infile)

            # Re-initialize trial column as 0
            word_data["trial"] = 0

            # Iterate through dataframe and set trial for rows with nouns and determiners
            # Increment the trial counters after each row, continuously across files (do not reset counters after each file)
            for idx, row in word_data.iterrows():
                if row[PART_OF_SPEECH_FIELD] == NOUN_POS_LABEL:
                    word_data.loc[idx, "trial"] = trial_counter_nouns
                    trial_counter_nouns += 1
                elif row[PART_OF_SPEECH_FIELD] == DET_POS_LABEL:
                    word_data.loc[idx, "trial"] = trial_counter_determiners
                    trial_counter_determiners += 1

            # Filter subject data
            subj_data = gaze2analysis.loc[
                (gaze2analysis["subj"] == subject_id) &
                (gaze2analysis["aoi_target"] == True) &  # noqa: E712
                (gaze2analysis["set"] == set_id) &
                (gaze2analysis["pattern"] == pattern_id)
            ]

            # Select only "trial_time" and "trial" columns
            subj_data = subj_data[["trial_time", "trial"]]

            # Group by "trial" and take the mean of all other columns (which is now just the "trial_time" column)
            subj_data = subj_data.groupby("trial", as_index=False).mean()

            # Merge subj_data into word_data
            word_data = word_data.merge(subj_data, how="left")

            # Write outfile CSV
            # Encode NA values explicitly as "NA" so downstream CSV readers don't misparse empty fields
            word_data.to_csv(word_outfile, index=False, na_rep="NA")
            logger.info(f"Wrote subject {subject_id} word data to {word_outfile}")

    return experiment


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Compute per-trial gaze-to-target fixation times and merge into word/ERP data.")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(os.path.abspath(args.config))
