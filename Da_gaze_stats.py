import argparse
import logging
import os
import re
import time
from datetime import timedelta

import pandas as pd
import rpy2.robjects as ro

from constants import (AUDIO_ERP_FILE_SUFFIX, CONDITIONS, DET_POS_LABEL,
                       NOUN_POS_LABEL, PART_OF_SPEECH_FIELD, PATTERN_IDS,
                       ROUND_N, RUN_CONFIG_KEY, SET_IDS, WORD_END_FIELD,
                       WORD_ONSET_FIELD)
from load_experiment import (create_experiment_outdir, get_experiment_id,
                             load_config, parse_subject_ids,
                             subject_files_dict)

logger = logging.getLogger(__name__)


def main(config: str | dict) -> dict:
    start_time = time.time()
    # Load experiment config
    if isinstance(config, str):
        config = load_config(config)
    experiment_id = get_experiment_id(config)

    # Get input/output paths
    input_dir = config["data"]["input"]["root"]
    output_dir = create_experiment_outdir(config, experiment_id)
    gaze_dir = config["data"]["input"]["gaze_dir"]
    gaze_outdir = os.path.join(output_dir, gaze_dir)
    audio_dir = config["data"]["input"]["audio_dir"]
    audio_indir = os.path.join(input_dir, audio_dir)
    audio_outdir = os.path.join(output_dir, audio_dir)

    # Get selected subject IDs and per-subject audio ERP files
    _, subject_id_regex = parse_subject_ids(config["experiment"]["subjects"])
    audio_erp_files = subject_files_dict(
        dir=audio_indir,
        subject_regex=subject_id_regex,
        suffix=AUDIO_ERP_FILE_SUFFIX,
    )
    subject_ids = sorted(list(audio_erp_files.keys()))
    logger.info(f"Processing {len(subject_ids)} subject ID(s): {', '.join(subject_ids)}")

    # Get overall gaze input file (output from Ca script). # TODO confirm that this should use the aggregated file, not individual per subject
    gaze_infile = os.path.join(gaze_outdir, "gaze_positions_all_4analysis.csv")

    # Load gaze_infile and drop all non-trial data points
    gaze2analysis = pd.read_csv(gaze_infile)
    gaze2analysis = gaze2analysis.loc[
        gaze2analysis["condition"].notna() &
        gaze2analysis["trial_time"].notna() &
        (~gaze2analysis["trackloss"]) &
        gaze2analysis["subj"].isin(subject_ids)
    ].drop_duplicates()

    # Add column "duration": tmax - time (rounded to ROUND_N digits)
    gaze2analysis["duration"] = (gaze2analysis[WORD_END_FIELD] - gaze2analysis[WORD_ONSET_FIELD]).round(ROUND_N)

    # Filter valid durations and then drop duration column
    gaze2analysis = gaze2analysis[
        gaze2analysis["duration"] >= 0
    ].drop(columns="duration")

    # Iterate over subject directories
    all_words = pd.DataFrame()
    for subject_id, subj_audio_erp_files in audio_erp_files.items():
        logger.info(f"Processing subject '{subject_id}'...")
        
        # Designate and create per-subject audio outdir (if doesn't already exist)
        subj_audio_outdir = os.path.join(audio_outdir, subject_id)
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
                (gaze2analysis["aoi_target"] == True) &
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
            word_data.to_csv(word_outfile)
            logger.info(f"Wrote subject {subject_id} word data to {word_outfile}")

            # Add word_data to running all_words dataframe
            all_words = pd.concat([all_words, word_data], axis=0, ignore_index=True)

    # Compute median determiner onset time
    median_d_onset = (
        gaze2analysis[gaze2analysis[PART_OF_SPEECH_FIELD] == DET_POS_LABEL]
        .loc[:, [PART_OF_SPEECH_FIELD, 'trial_time']]
        .groupby(PART_OF_SPEECH_FIELD, as_index=False)
        .median()
        # convert from seconds to milliseconds 
        .loc[0, 'trial_time'] * 1000.
    )

    # Compute median noun offset
    filtered_nouns = gaze2analysis[
        gaze2analysis["condition"].isin(CONDITIONS) &
        (gaze2analysis[PART_OF_SPEECH_FIELD] == NOUN_POS_LABEL)
    ]
    filtered_nouns["duration"] = (filtered_nouns[WORD_END_FIELD] - filtered_nouns[WORD_ONSET_FIELD]).round(ROUND_N)
    median_noun_offset = (
        filtered_nouns.loc[:, [PART_OF_SPEECH_FIELD, "duration"]]
        .groupby(PART_OF_SPEECH_FIELD, as_index=False)
        .median()
        .loc[0, "duration"] * 1000   
    )

    # Calculate duration of this step and add to run config
    end_time = time.time()
    duration = str(timedelta(seconds=int(end_time - start_time)))
    config[RUN_CONFIG_KEY]["duration"]["Da_gaze_stats"] = duration
    logger.info(f"Step Da completed successfully (duration: {duration}).")
    return config


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Process fixation files from pupil player and prepare them for further processing.")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(args.config)
