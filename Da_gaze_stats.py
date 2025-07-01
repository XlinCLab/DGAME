import argparse
import logging
import os
import re
import time
from datetime import timedelta

import pandas as pd
from rpy2.robjects import StrVector
from rpy2.robjects.packages import importr

from constants import (AUDIO_ERP_FILE_SUFFIX, CONDITIONS, CONFLICT_LABEL,
                       DET_POS_LABEL, NO_CONFLICT_LABEL, NOUN_POS_LABEL,
                       PART_OF_SPEECH_FIELD, PATTERN_IDS, ROUND_N,
                       RUN_CONFIG_KEY, SET_IDS, WORD_END_FIELD,
                       WORD_ONSET_FIELD)
from load_experiment import (create_experiment_outdir, get_experiment_id,
                             load_config, parse_subject_ids,
                             subject_files_dict)
from r_utils import (convert_pandas2r_dataframe, convert_r2pandas_dataframe,
                     r_assign, r_install_packages, r_interface)

logger = logging.getLogger(__name__)

# Load and/or install R dependencies
r_install_packages([
    "dplyr",
    "eyetrackingR",
])
eyetrackingr = importr("eyetrackingR")


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

    # (!) temp hack to have multiple subjects: duplicate existing subj 02 data
    audio_erp_files['pseudo02'] = audio_erp_files['02']

    subject_ids = sorted(list(audio_erp_files.keys()))
    n_subjects = len(subject_ids)
    logger.info(f"Processing {len(subject_ids)} subject ID(s): {', '.join(subject_ids)}")
    # temp fix since using hardcoded input file:
    subject_ids.append(2)

    # Get overall gaze input file (output from Ca script). # TODO confirm that this should use the aggregated file, not individual per subject
    # gaze_infile = os.path.join(gaze_outdir, "gaze_positions_all_4analysis.csv")
    # TODO change back once I have figured out why the input files are different
    gaze_infile = "/Users/pgeorgis/Documents/projects/dgame2/dg2_analysis/gaze_positions/gaze_positions_4analysis.csv"

    # Load gaze_infile and drop all non-trial data points
    gaze2analysis = pd.read_csv(gaze_infile)

    # (!) temp hack to have multiple subjects: duplicate existing subj 02 data
    pseudodata = gaze2analysis.copy()
    pseudodata["subj"] = "pseudo02"
    gaze2analysis = pd.concat([gaze2analysis, pseudodata], axis=0, ignore_index=False)

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

    # DEBUG: numbers match through here if using same "gaze_positions_4analysis.csv" input file as R script
    # However, if using file output by Ca script, the numbers are different (60 rows difference in filtered result)
    # TODO identify why there is a discrepancy in the input file

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
            # temp fix because using hardcoded input file with subject_id = 2
            subject_id_hack = 2
            subj_data = gaze2analysis.loc[
                (gaze2analysis["subj"] == subject_id_hack) &
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

    # Convert trial_time seconds to milliseconds
    gaze2analysis["trial_time"] *= 1000

    # Convert gaze2analysis to R dataframe
    r_gaze2analysis = convert_pandas2r_dataframe(gaze2analysis)

    # Call make_eyetrackingr_data function from within R
    aoi_column_vector = StrVector([
        "aoi_target",
        "aoi_comp",
        "aoi_otherTarget",
        "aoi_otherComp",
        "aoi_fillerA",
        "aoi_fillerB",
        "aoi_goal",
    ])
    eyetrackingr_data = eyetrackingr.make_eyetrackingr_data(
        r_gaze2analysis,
        participant_column="subj",
        trial_column="trial",
        time_column="trial_time",
        trackloss_column="trackloss",
        aoi_columns=aoi_column_vector,
        treat_non_aoi_looks_as_missing=False
    )
    r_interface.assign("eyetrackingr_data", eyetrackingr_data)
    r_interface("eyetrackingr_data$condition <- as.factor(eyetrackingr_data$condition)")

    # Transform data for AOI analysis from within R
    response_time = eyetrackingr.make_time_sequence_data(
        eyetrackingr_data,
        time_bin_size=100,
        predictor_columns=StrVector(["condition"]),
        aois=aoi_column_vector,
        summarize_by=StrVector(["subj"])
    )
    r_interface.assign("response_time", response_time)
    response_time_pd = convert_r2pandas_dataframe(response_time)
    # Confirm that the number of subjects matches the number of subject IDs selected earlier
    assert n_subjects == len(response_time_pd['subj'].unique())

    # Statistical thresholds
    # Pick threshold t based on alpha = 0.05, two tailed
    threshold_t = r_assign("threshold_t", f"qt(p = 1 - .05/2, df = {n_subjects} - 1)")
    n_bins = r_assign("n_bins", "length(unique(response_time$TimeBin))")

    # Data for competitor analysis
    response_time_comp = eyetrackingr.make_time_sequence_data(
        eyetrackingr_data,
        time_bin_size=100,
        predictor_columns=StrVector(["condition"]),
        aois=StrVector([
            "aoi_comp",
            "aoi_otherTarget",
            "aoi_fillerA",
            "aoi_fillerB",
            "aoi_otherComp",
        ]),
        summarize_by=StrVector(["subj"])
    )
    response_time_comp = convert_r2pandas_dataframe(response_time_comp)
    response_time_comp["condition"] = response_time_comp["condition"].astype("category")
    response_time_comp["subj"] = response_time_comp["subj"].astype("category")
    response_time_comp["AOI"] = response_time_comp["AOI"].apply(
        lambda x: "aoi_AllOther" if x != "aoi_comp" else x
    )
    aoi_comp_df = response_time_comp[response_time_comp["AOI"] == "aoi_comp"]
    numeric_cols = response_time_comp.select_dtypes(include='number').columns
    other_aoi_df = (
        response_time_comp[response_time_comp["AOI"] == "aoi_AllOther"]
        .groupby(['subj', 'condition', 'TimeBin', 'Time', 'AOI'], observed=True, as_index=False)[numeric_cols]
        .mean()
    )
    response_time_comp = pd.concat([aoi_comp_df, other_aoi_df], axis=0, ignore_index=True)

    # Bootstrap for comp
    response_time_comp["aoi_fct"] = response_time_comp["AOI"]
    response_time_comp = response_time_comp[response_time_comp["condition"] == CONFLICT_LABEL]
    response_time_comp["AOI"] = "dummy"
    response_time_comp["aoi_fct"] = response_time_comp["aoi_fct"].astype("category")
    response_time_comp = response_time_comp[response_time_comp["aoi_fct"].isin({"aoi_comp", "aoi_AllOther"})]
    response_time_comp = convert_pandas2r_dataframe(response_time_comp)

    # Time cluster
    time_cluster_data_target = eyetrackingr.make_time_cluster_data(
        data=response_time,
        predictor_column="condition",
        aoi="aoi_target",
        test="t.test",
        paired=True,
        threshold=threshold_t,
    )
    cluster_analysis_target = eyetrackingr.analyze_time_clusters(
        time_cluster_data_target,
        samples=4000,
        within_subj=True,
        paired=True,
    )
    # Test for competitors
    time_cluster_data_comp = eyetrackingr.make_time_cluster_data(
        response_time_comp,
        predictor_column="aoi_fct",
        aoi="dummy",
        test="t.test",
        threshold=threshold_t,
        treatment_level="aoi_AllOther",
        paired=True
    )
    cluster_analysis_comp = eyetrackingr.analyze_time_clusters(
        time_cluster_data_comp,
        within_subj=True,
        paired=True,
        samples=4000,
    )
    # Test for goal
    time_cluster_data_goal = eyetrackingr.make_time_cluster_data(
        response_time,
        predictor_column="condition",
        aoi="aoi_goal",
        test="t.test",
        threshold=threshold_t,
        treatment_level=NO_CONFLICT_LABEL,
        paired=True,
    )
    cluster_analysis_goal = eyetrackingr.analyze_time_clusters(
        time_cluster_data_goal,
        within_subj=True,
        paired=True,
        samples=4000,
    )
    # TODO why is the result below named the same as an object above, and not used after this?
    time_cluster_data_comp = eyetrackingr.make_time_cluster_data(
        response_time_comp,
        predictor_column="timepoint",
        aoi="aoi_comp",
        test="t.test",
        paired=True,
        threshold=threshold_t,
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
