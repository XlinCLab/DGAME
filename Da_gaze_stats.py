import argparse
import logging
import os
import re
import time
from datetime import timedelta

import numpy as np
import pandas as pd
import rpy2.robjects as robjects
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
from r_utils import (RDataFrame, convert_pandas2r_dataframe,
                     convert_r2pandas_dataframe, r_eval, r_install_packages,
                     r_interface)
from utils import generate_variable_name

logger = logging.getLogger(__name__)

# Load and/or install R dependencies
r_install_packages([
    "dplyr",
    "eyetrackingR",
    "ggplot2",
    "pbapply",
])
eyetrackingr = importr("eyetrackingR")

# Source R script with custom plotting function
robjects.r["source"]("plot_gaze_proportions.R")
plot_ti1 = robjects.globalenv["plot_ti1"]
plot_ti3 = robjects.globalenv["plot_ti3"]


def compute_median_det_onset(df: pd.DataFrame) -> np.float64:
    """Compute median determiner onset time from a dataframe with POS and time annotations."""
    median_determiner_onset = (
        df[df[PART_OF_SPEECH_FIELD] == DET_POS_LABEL]
        .loc[:, [PART_OF_SPEECH_FIELD, 'trial_time']]
        .groupby(PART_OF_SPEECH_FIELD, as_index=False)
        .median()
        # convert from seconds to milliseconds
        .loc[0, 'trial_time'] * 1000.
    )
    return median_determiner_onset


def compute_median_noun_offset(df: pd.DataFrame) -> np.float64:
    """Compute the median noun offset from a dataframe with POS and duration annotation."""
    filtered_nouns = df[
        df["condition"].isin(CONDITIONS) &
        (df[PART_OF_SPEECH_FIELD] == NOUN_POS_LABEL)
    ]
    filtered_nouns["duration"] = (filtered_nouns[WORD_END_FIELD] - filtered_nouns[WORD_ONSET_FIELD]).round(ROUND_N)
    median_noun_offset = (
        filtered_nouns.loc[:, [PART_OF_SPEECH_FIELD, "duration"]]
        .groupby(PART_OF_SPEECH_FIELD, as_index=False)
        .median()
        .loc[0, "duration"] * 1000
    )
    return median_noun_offset


def r_postprocess_response_time_df(response_time_df: RDataFrame,
                                   aoi_label: str="AOI",
                                   aoi_comp_label: str="aoi_comp",
                                   aoi_other_label: str="aoi_AllOther",
                                   aoi_fct_label: str="aoi_fct",
                                   subject_label: str="subj",
                                   condition_label: str="condition",
                                   time_bin_label: str="TimeBin",
                                   time_label: str="Time",
                                   dummy_label: str="dummy",
                                   ) -> RDataFrame:
    """Perform postprocessing (in R) on a dataframe which has already been run through make_eyetrackingr_data."""
    temp_df_name = generate_variable_name()
    r_interface.assign(temp_df_name, response_time_df)
    df_processed = r_eval(
        f"{temp_df_name} %>% mutate({aoi_label} = case_when({aoi_label} != '{aoi_comp_label}' ~ '{aoi_other_label}', TRUE ~ {aoi_label}))",
        name="df_processed"
    )
    aoi_comp_df = r_eval(
        f"df_processed %>% filter({aoi_label} == '{aoi_comp_label}')",
        name="aoi_comp_df"
    )
    other_aoi_df = r_eval(
        f"df_processed %>% filter({aoi_label} == '{aoi_other_label}') %>% group_by({subject_label},{condition_label},{time_bin_label},{time_label},{aoi_label}) %>% summarise_all(mean)",
        name="other_aoi_df"
    )
    combined_df = r_eval(
        "rbind(aoi_comp_df, other_aoi_df)",
        name="combined_df"
    )

    # Bootstrap for comp
    r_interface(f"combined_df${aoi_fct_label} <- combined_df${aoi_label}")
    final_df = r_eval(
        f"combined_df %>% filter({condition_label} == '{CONFLICT_LABEL}')",
        name="final_df"
    )
    r_interface(f"final_df${aoi_label} <- '{dummy_label}'")
    r_interface(f"final_df${aoi_fct_label} <- as.factor(final_df${aoi_fct_label})")
    final_df = r_eval(
        f"final_df %>% filter({aoi_fct_label} == '{aoi_comp_label}' | {aoi_fct_label} == '{aoi_other_label}')",
        name="final_df"
    )

    return final_df


def run_time_cluster_analysis(response_time_df: RDataFrame,
                              response_time_comp_df: RDataFrame,
                              threshold_t: float,
                              ) -> dict:
    """Run time cluster analysis on gaze response time dataframes."""
    logger.info("Analyzing time clusters for target...")
    time_cluster_data_target = eyetrackingr.make_time_cluster_data(
        data=response_time_df,
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
    logger.info("Analyzing time clusters for competitors...")
    time_cluster_data_comp = eyetrackingr.make_time_cluster_data(
        response_time_comp_df,
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
    logger.info("Analyzing time clusters for goal...")
    time_cluster_data_goal = eyetrackingr.make_time_cluster_data(
        response_time_df,
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

    result_dict = {
        "target": (time_cluster_data_target, cluster_analysis_target),
        "comp": (time_cluster_data_comp, cluster_analysis_comp),
        "goal": (time_cluster_data_goal, cluster_analysis_goal),
    }
    return result_dict


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
    n_subjects = len(subject_ids)
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
    median_det_onset = compute_median_det_onset(gaze2analysis)

    # Compute median noun offset
    median_noun_offset = compute_median_noun_offset(gaze2analysis)

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
    threshold_t = r_eval(f"qt(p = 1 - .05/2, df = {n_subjects} - 1)")
    n_bins = r_eval("length(unique(response_time$TimeBin))")

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
    response_time_comp = r_postprocess_response_time_df(response_time_comp)

    # Time cluster analysis
    time_cluster_analysis_active = config.get("analysis", {}).get("time_cluster_analysis", False)
    if time_cluster_analysis_active:
        logger.info("Starting time cluster analysis...")
        cluster_analysis_results = run_time_cluster_analysis(
            response_time_df=response_time,
            response_time_comp_df=response_time_comp,
            threshold_t=threshold_t,
        )
        logger.info("Time cluster analysis finished")
        # TODO no further steps implemented here for what to do after running cluster analysis
    else:
        logger.info("Skipping time cluster analysis")

    # Plot results
    gaze_plot_outdir = os.path.join(gaze_outdir, "plots")
    os.makedirs(gaze_plot_outdir, exist_ok=True)
    plotti1_out = os.path.join(gaze_plot_outdir, "gaze_proportions_ti1.png") # TODO needs better name
    plot_ti1(response_time, float(median_det_onset), outfile=plotti1_out)
    logger.info(f"Plotted to {plotti1_out}")
    plotti3_out = os.path.join(gaze_plot_outdir, "gaze_proportions_ti3.png") # TODO needs better name
    plot_ti3(response_time_comp, float(median_noun_offset), outfile=plotti3_out)
    logger.info(f"Plotted to {plotti3_out}")

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
