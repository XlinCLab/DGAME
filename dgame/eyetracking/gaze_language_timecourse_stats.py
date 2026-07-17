import argparse
import logging
import os

import numpy as np
import pandas as pd
import rpy2.robjects as robjects
from rpy2.robjects import StrVector
from rpy2.robjects.packages import importr

from dgame.constants import (CONDITIONS, CONFLICT_LABEL, DET_POS_LABEL,
                             NO_CONFLICT_LABEL, NOUN_POS_LABEL,
                             PART_OF_SPEECH_FIELD, R_PLOT_SCRIPT_DIR, ROUND_N,
                             WORD_END_FIELD, WORD_ONSET_FIELD)
from dgame.eyetracking.utils import load_filtered_gaze_data
from dgame.pipeline import ET_GAZE_LANG_STATS_STEP
from experiment.load_experiment import Experiment
from utils.r_dependencies import R_DEPENDENCIES, r_install_packages
from utils.r_utils import (RDataFrame, convert_pandas2r_dataframe,
                           convert_r2pandas_dataframe, r_eval, r_interface)
from utils.utils import generate_variable_name

# Load and/or install R dependencies
r_install_packages(R_DEPENDENCIES)
eyetrackingr = importr("eyetrackingR")

# Source R script with custom plotting function
robjects.r["source"](os.path.join(R_PLOT_SCRIPT_DIR, "plot_gaze_proportions.R"))
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
    ].copy()
    filtered_nouns["duration"] = (filtered_nouns[WORD_END_FIELD] - filtered_nouns[WORD_ONSET_FIELD]).round(ROUND_N)
    median_noun_offset = (
        filtered_nouns.loc[:, [PART_OF_SPEECH_FIELD, "duration"]]
        .groupby(PART_OF_SPEECH_FIELD, as_index=False)
        .median()
        .loc[0, "duration"] * 1000
    )
    return median_noun_offset


def r_postprocess_response_time_df(response_time_df: RDataFrame,
                                   aoi_label: str = "AOI",
                                   aoi_comp_label: str = "aoi_comp",
                                   aoi_other_label: str = "aoi_AllOther",
                                   aoi_fct_label: str = "aoi_fct",
                                   subject_label: str = "subj",
                                   condition_label: str = "condition",
                                   time_bin_label: str = "TimeBin",
                                   time_label: str = "Time",
                                   dummy_label: str = "dummy",
                                   ) -> RDataFrame:
    """Perform postprocessing (in R) on a dataframe which has already been run through make_eyetrackingr_data."""
    temp_df_name = generate_variable_name()
    r_interface.assign(temp_df_name, response_time_df)
    df_processed = r_eval(  # noqa: F841
        f"{temp_df_name} %>% mutate({aoi_label} = case_when({aoi_label} != '{aoi_comp_label}' ~ '{aoi_other_label}', TRUE ~ {aoi_label}))",
        name="df_processed"
    )
    aoi_comp_df = r_eval(  # noqa: F841
        f"df_processed %>% filter({aoi_label} == '{aoi_comp_label}')",
        name="aoi_comp_df"
    )
    other_aoi_df = r_eval(  # noqa: F841
        f"df_processed %>% filter({aoi_label} == '{aoi_other_label}') %>% group_by({subject_label},{condition_label},{time_bin_label},{time_label},{aoi_label}) %>% summarise_all(mean)",  # noqa: E501
        name="other_aoi_df"
    )
    combined_df = r_eval(  # noqa: F841
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
                              logger: logging.Logger,
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


def main(experiment: str | dict | Experiment) -> Experiment:
    # Initialize DGAME experiment from config
    if not isinstance(experiment, Experiment):
        from dgame.dgame import DGAME
        experiment = DGAME.from_input(experiment)
    logger = experiment.logger

    n_subjects = len(experiment.subject_ids)
    logger.info(f"Processing {n_subjects} subject ID(s): {', '.join(experiment.subject_ids)}")

    # Load and filter gaze data (output of step Ca)
    gaze2analysis = load_filtered_gaze_data(experiment)

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
    # n_bins = r_eval("length(unique(response_time$TimeBin))")

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
    time_cluster_analysis_active = experiment.get_dgame_step_parameter(
        ET_GAZE_LANG_STATS_STEP, "time_cluster_analysis",
        default=False
    )
    if time_cluster_analysis_active:
        if len(experiment.subject_ids) >= 2:
            logger.info("Starting time cluster analysis...")
            cluster_analysis_results = run_time_cluster_analysis(
                response_time_df=response_time,
                response_time_comp_df=response_time_comp,
                threshold_t=threshold_t,
                logger=logger,
            )
            logger.info("Time cluster analysis finished")
            # TODO no further steps implemented here for what to do after running cluster analysis
        else:
            logger.info("Skipping time cluster analysis due to insufficient subjects (at least 2 required)")
    else:
        logger.info("Skipping time cluster analysis")

    # Plot results
    gaze_plot_outdir = os.path.join(experiment.gaze_outdir, "plots")
    os.makedirs(gaze_plot_outdir, exist_ok=True)
    plotti1_out = os.path.join(gaze_plot_outdir, "gaze_proportions_ti1.png")  # TODO needs better name
    plot_ti1(response_time, float(median_det_onset), outfile=plotti1_out)
    logger.info(f"Plotted to {plotti1_out}")
    plotti3_out = os.path.join(gaze_plot_outdir, "gaze_proportions_ti3.png")  # TODO needs better name
    plot_ti3(response_time_comp, float(median_noun_offset), outfile=plotti3_out)
    logger.info(f"Plotted to {plotti3_out}")

    return experiment


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Compute gaze x language time-course statistics and plots.")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(os.path.abspath(args.config))
