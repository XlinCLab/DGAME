import argparse
import glob
import logging
import os
from math import floor

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.regression.linear_model import RegressionResultsWrapper

from dgame.constants import CHANNEL_FIELD
from experiment.load_experiment import Experiment
from experiment.test_subjects import subject_dirs_dict
from utils.utils import load_csv_list

logger = logging.getLogger(__name__)


def load_unfold_out_fixation_data(per_subject_unfold_out_dirs: dict,
                                  channel_coords: pd.DataFrame,
                                  ) -> pd.DataFrame:
    """Load data from fixation CSV files within unfold_out EEG output directories for all subjects into a single dataframe."""
    all_subjects_fixation_data = pd.DataFrame()
    for subject, unfold_out_dir in per_subject_unfold_out_dirs.items():
        # Load all unfold_FIX CSV files into a single dataframe
        unfold_fix_files = [
            os.path.join(unfold_out_dir, filepath)
            for filepath in glob.glob(f"{subject}_*_unfold_FIX.csv", root_dir=unfold_out_dir)
        ]
        unfold_fix_data = load_csv_list(unfold_fix_files)
        unfold_fix_data[CHANNEL_FIELD] = unfold_fix_data[CHANNEL_FIELD].astype(str)

        # Merge with channel coordinates
        unfold_fix_data = unfold_fix_data.merge(channel_coords, how="left", on=CHANNEL_FIELD)

        # Compute laterality
        unfold_fix_data["laterality"] = np.where(
            unfold_fix_data["lat"] < 0, "left",
            np.where(unfold_fix_data["lat"] > 0, "right", "central")
        )

        # Compute saggitality
        sag_conditions = [
            (unfold_fix_data["sag"] > 0) & (unfold_fix_data["sag"] <= 0.0714),  # frontal
            (unfold_fix_data["sag"] > 0.0714),                                  # prefrontal
            (unfold_fix_data["sag"] < 0) & (unfold_fix_data["sag"] >= -0.0929), # posterior
            (unfold_fix_data["sag"] < -0.0929)                                  # occipital
            # elsewhere condition                                               # central
        ]
        sag_labels = ["frontal", "prefrontal", "posterior", "occipital"]
        unfold_fix_data["saggitality"] = np.select(sag_conditions, sag_labels, default="central")

        # Compute fixation time labels
        fix_time_conditions = [
            (unfold_fix_data["trial_time"] < -1),                                           # >1s_before_noun
            (unfold_fix_data["trial_time"] < 0) & (unfold_fix_data["trial_time"] >= -1),    # before_noun
            (unfold_fix_data["trial_time"] >= 0) & (unfold_fix_data["trial_time"] <= 0.5),  # during_noun
            (unfold_fix_data["trial_time"] > 0.471) & (unfold_fix_data["trial_time"] <= 1)  # after_noun
            # elsewhere condition                                                           # 1s_after_noun
        ]
        fix_time_labels = [
            ">1s_before_noun",
            "before_noun",
            "during_noun",
            "after_noun",
        ]
        unfold_fix_data["fix_time"] = np.select(fix_time_conditions, fix_time_labels, default=">1s_after_noun")

        # Add subject unfold fixation data into running dataframe for all subjects
        all_subjects_fixation_data = pd.concat([all_subjects_fixation_data, unfold_fix_data], axis=0, ignore_index=True)
    
    return all_subjects_fixation_data


def create_time_windows(data: pd.DataFrame,
                        window_size: int = 100,
                        ) -> pd.DataFrame:
    data["time_bin"] = data["time"].apply(lambda t: floor(t / window_size) * window_size)
    aggregated_data = (
        data
        .groupby(["time_bin", "laterality", "saggitality", "condition", "fix_time", "subject"], as_index=False)
        .agg(data_mean=("data", "mean"))
    )
    return aggregated_data


def summarize_stats_model(model: RegressionResultsWrapper) -> pd.DataFrame:
    """Return summary dataframe for a fitted regression model."""
    model_summary = pd.DataFrame({
        'predictor': model.params.index,
        'coef': model.params.values,
        'p_value': model.pvalues.values
    })
    return model_summary


def summarize_regression_by_time_bins(fixation_data_windows: pd.DataFrame) -> pd.DataFrame:
    """Fit a regression model for each time bin using all variables as predictors and collect model summaries into a dataframe."""
    model_formula = "data_mean ~ laterality * saggitality * condition * fix_time * baseline"
    unique_time_bins = fixation_data_windows["time_bin"].unique()
    time_bin_model_summaries = pd.DataFrame()
    for time_window in unique_time_bins:
        filtered_time_window_fixation = (
            fixation_data_windows
            .loc[(fixation_data_windows["time_bin"] == time_window)]
        )
        model = smf.ols(formula=model_formula, data=filtered_time_window_fixation).fit()
        n_predictors = model.model.exog.shape[-1]
        design_matrix_rank = np.linalg.matrix_rank(model.model.exog)
        if design_matrix_rank < n_predictors:
            logger.error(f"Model (time_window = {time_window}) is overparametrized: {n_predictors} predictors but only {design_matrix_rank} observations")
            continue
        model_summary = summarize_stats_model(model)
        time_bin_model_summaries = pd.concat([time_bin_model_summaries, model_summary], axis=0, ignore_index=True)

    return time_bin_model_summaries


def main(experiment: str | dict | Experiment) -> Experiment:
    # Initialize DGAME experiment from config
    if not isinstance(experiment, Experiment):
        from dgame.dgame import DGAME
        experiment = DGAME.from_input(experiment)

    # Load channel coords file
    channel_coords = experiment.load_channel_coords()

    # Find per-subject EEG output directories
    per_subject_eeg_outdirs = subject_dirs_dict(
        root_dir = experiment.eeg_outdir,
        subject_regex=experiment.subject_id_regex,
    )
    assert all(len(eeg_outdir) == 1 for eeg_outdir in per_subject_eeg_outdirs.values())
    # Find per-subject unfold_out directories within EEG outdirs
    per_subject_unfold_out_dirs = {
        subject_id: os.path.join(subject_eeg_outdirs[0], "unfold_out")
        for subject_id, subject_eeg_outdirs in per_subject_eeg_outdirs.items()
    }

    # Load fixation data for all subjects
    logger.info(f"Loading EEG fixation data for {len(per_subject_unfold_out_dirs)} subject(s)...")
    all_subjects_fixation_data = load_unfold_out_fixation_data(
        per_subject_unfold_out_dirs, channel_coords
    )
    
    # Compute baseline factor  # TODO baseline of what?
    # Summarize "data" field as "baseline" for every combination of subject/condition/laterality/saggitality/fix_time/fix_at
    # from a filtered subset of the data containing only time >= -250 but < 0
    # and where the fixation is on the target
    logger.info("Filtering fixation data to establish baseline...")
    baseline_data = (
        all_subjects_fixation_data
        # filter time >= -250 & time < 0
        .loc[(all_subjects_fixation_data["time"] >= -250) & (all_subjects_fixation_data["time"] < 0)]
        # select only these columns
        [["data", "subject", "condition", "laterality", "saggitality", "fix_time", "fix_at"]]
        # filter fix_at == "target"
        .loc[lambda df: df["fix_at"] == "target"]
    )
    baseline_aggregated = (
        baseline_data
        .groupby(["subject", "condition", "laterality", "saggitality", "fix_time", "fix_at"], as_index=False)
        .agg(baseline=("data", "mean"))
    )
    # Add time window bins and aggregate with mean
    logger.info("Aggregating fixation data by time window...")
    fixation_data_windows = create_time_windows(all_subjects_fixation_data)
    # Merge with aggregated_baseline data
    fixation_data_windows = fixation_data_windows.merge(
        baseline_aggregated,
        how="inner",
        on=["subject", "condition", "laterality", "saggitality", "fix_time"]
    )
    # Filter time windows >= 0 and <= 1000
    fixation_data_windows = (
        fixation_data_windows
        .loc[(fixation_data_windows["time_bin"] >= 0) & (fixation_data_windows["time_bin"] <= 1000)]
    )

    # Fit regression models for filtered data by time bin, using all predictors, and collect summaries
    logger.info("Fitting regression models on data per time window...")
    time_bin_model_summaries = summarize_regression_by_time_bins(fixation_data_windows)
    if len(time_bin_model_summaries) > 0:
        # Write model summaries to output csv
        regression_summary_outfile = os.path.join(
            experiment.fixations_outdir,
            "regression_results_fixations.csv",
        )
        time_bin_model_summaries.to_csv(regression_summary_outfile, index=False)
        logger.info(f"Wrote regression model summaries to {regression_summary_outfile}")
    else:
        logger.warning("No valid regression models were produced.")
    
    

    return experiment

if __name__ == "__main__":
    parser = argparse.ArgumentParser("Compute and plot fixation statistics.")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(os.path.abspath(args.config))
