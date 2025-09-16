import argparse
import glob
import logging
import os
from math import floor

import numpy as np
import pandas as pd
import rpy2.robjects as robjects
import statsmodels.formula.api as smf
from rpy2.rinterface_lib.embedded import RRuntimeError
from statsmodels.regression.linear_model import RegressionResultsWrapper
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm

from dgame.constants import (CHANNEL_FIELD, LATERAL_INPUT_FIELD,
                             LATERALITY_FIELD, R_PLOT_SCRIPT_DIR,
                             SAGGITAL_INPUT_FIELD, SAGGITALITY_FIELD,
                             STEP_JA_KEY)
from experiment.load_experiment import Experiment
from utils.r_utils import convert_pandas2r_dataframe
from utils.utils import load_csv_list

logger = logging.getLogger(__name__)


# Source R script with custom plotting function
robjects.r["source"](os.path.join(R_PLOT_SCRIPT_DIR, "plot_language_fixation_stats.R"))
create_language_fixation_plot = robjects.globalenv["create_language_fixations_plot"]


def annotate_laterality_and_saggitality(df: pd.DataFrame) -> pd.DataFrame:
        """Annotate pandas DataFrame with laterality and saggitality labels."""
        if SAGGITAL_INPUT_FIELD not in df.columns:
            raise ValueError(f"DataFrame missing '{SAGGITAL_INPUT_FIELD}' column")
        if LATERAL_INPUT_FIELD not in df.columns:
            raise ValueError(f"DataFrame missing '{LATERAL_INPUT_FIELD}' column")

        # Compute laterality
        df[LATERALITY_FIELD] = np.where(
            df[LATERAL_INPUT_FIELD] < 0, "left",
            np.where(df[LATERAL_INPUT_FIELD] > 0, "right", "central")
        )

        # Compute saggitality
        sag_conditions = [
            (df[SAGGITAL_INPUT_FIELD] > 0) & (df[SAGGITAL_INPUT_FIELD] <= 0.0714),  # frontal
            (df[SAGGITAL_INPUT_FIELD] > 0.0714),                                  # prefrontal
            (df[SAGGITAL_INPUT_FIELD] < 0) & (df[SAGGITAL_INPUT_FIELD] >= -0.0929), # posterior
            (df[SAGGITAL_INPUT_FIELD] < -0.0929)                                  # occipital
            # elsewhere condition                                               # central
        ]
        sag_labels = ["frontal", "prefrontal", "posterior", "occipital"]
        df[SAGGITALITY_FIELD] = np.select(sag_conditions, sag_labels, default="central")

        return df


def load_unfold_out_data(per_subject_unfold_out_dirs: dict,
                         channel_coords: pd.DataFrame,
                         mode: str,
                         ) -> pd.DataFrame:
    """Load data from fixation CSV files within unfold_out EEG output directories for all subjects into a single dataframe."""
    if mode not in {"FIX", "N"}:
        raise ValueError(f"Expected mode to be either 'FIX' or 'N', got '{mode}'")
    all_subjects_unfold_out_data = pd.DataFrame()
    for subject, unfold_out_dir in per_subject_unfold_out_dirs.items():
        # Load all unfold_FIX / unfold_N CSV files into a single dataframe
        unfold_out_files = [
            os.path.join(unfold_out_dir, filepath)
            for filepath in glob.glob(f"{subject}_*_unfold_{mode}.csv", root_dir=unfold_out_dir)
        ]
        unfold_out_data = load_csv_list(
            unfold_out_files,
            progress_bar_description=f"Loading {len(unfold_out_files)} EEG {mode} files for subject <{subject}> ..."
        )
        unfold_out_data[CHANNEL_FIELD] = unfold_out_data[CHANNEL_FIELD].astype(str)

        # Merge with channel coordinates
        unfold_out_data = unfold_out_data.merge(channel_coords, how="left", on=CHANNEL_FIELD)

        # Annotate laterality and saggitality
        logger.info("Annotating laterality and saggitality...")
        unfold_out_data = annotate_laterality_and_saggitality(unfold_out_data)

        # Annotate fixation time labels
        if mode == "FIX":
            logger.info("Annotating fixation times relative to noun...")
            fix_time_conditions = [
                (unfold_out_data["trial_time"] < -1),                                           # >1s_before_noun
                (unfold_out_data["trial_time"] < 0) & (unfold_out_data["trial_time"] >= -1),    # before_noun
                (unfold_out_data["trial_time"] >= 0) & (unfold_out_data["trial_time"] <= 0.5),  # during_noun
                (unfold_out_data["trial_time"] > 0.471) & (unfold_out_data["trial_time"] <= 1)  # after_noun
                # elsewhere condition                                                           # 1s_after_noun
            ]
            fix_time_labels = [
                ">1s_before_noun",
                "before_noun",
                "during_noun",
                "after_noun",
            ]
            unfold_out_data["fix_time"] = np.select(fix_time_conditions, fix_time_labels, default=">1s_after_noun")

        # Add subject unfold fixation data into running dataframe for all subjects
        all_subjects_unfold_out_data = pd.concat([all_subjects_unfold_out_data, unfold_out_data], axis=0, ignore_index=True)
    
    return all_subjects_unfold_out_data


def create_time_windows(data: pd.DataFrame,
                        mode: str,
                        window_size: int = 100,
                        ) -> pd.DataFrame:
    if mode not in {"FIX", "N"}:
        raise ValueError(f"Expected mode to be either 'FIX' or 'N', got '{mode}'")
    if mode == "FIX":
        groupby = ["time_bin", LATERALITY_FIELD, SAGGITALITY_FIELD, "condition", "fix_time", "subject"]
    elif mode == "N":
        groupby = ["time_bin", LATERALITY_FIELD, SAGGITALITY_FIELD, "condition", "mean_target_fixation", "subject"]
    data["time_bin"] = data["time"].apply(lambda t: floor(t / window_size) * window_size)
    aggregated_data = (
        data
        .groupby(groupby, as_index=False)
        .agg(data_mean=("data", "mean"))
    )
    return aggregated_data


def summarize_stats_model(model: RegressionResultsWrapper) -> pd.DataFrame:
    """Return summary dataframe for a fitted regression model."""
    model_summary = pd.DataFrame({
        'predictor': model.params.index,
        'coef': model.params.values,
        'p_value': model.pvalues.values,
        't_value': model.tvalues.values,
    })
    return model_summary


def summarize_regression_by_time_bins(fixation_data_windows: pd.DataFrame, mode: str) -> pd.DataFrame:
    """Fit a regression model for each time bin using all variables as predictors and collect model summaries into a dataframe."""
    if mode not in {"FIX", "N"}:
        raise ValueError(f"Expected mode to be either 'FIX' or 'N', got '{mode}'")
    if mode == "FIX":
        model_formula = f"data_mean ~ {LATERALITY_FIELD} * {SAGGITALITY_FIELD} * condition * fix_time * baseline"
    elif mode == "N":
        model_formula = f"data_mean ~ {LATERALITY_FIELD} * {SAGGITALITY_FIELD} * condition * mean_target_fixation * baseline"
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
        if design_matrix_rank <= n_predictors:
            logger.error(f"Model (time_window = {time_window}) is overparametrized: {n_predictors} predictors but only {design_matrix_rank} observations")
            continue
        model_summary = summarize_stats_model(model)
        time_bin_model_summaries = pd.concat([time_bin_model_summaries, model_summary], axis=0, ignore_index=True)

    return time_bin_model_summaries


def time_bin_permutation(fixation_data_windows: pd.DataFrame,
                         time_bin: int,
                         mode: str,
                         n_permutations: int = 2000,
                         include_baseline: bool = False,
                         ) -> pd.DataFrame:
    if mode not in {"FIX", "N"}:
        raise ValueError(f"Expected mode to be either 'FIX' or 'N', got '{mode}'")
    # Filter data for this time bin
    filtered_time_bin_data = (
        fixation_data_windows
        .loc[(fixation_data_windows["time_bin"] == time_bin)]
    )

    # Construct model formula
    if mode == "FIX":
        model_formula = f"data_mean ~ {LATERALITY_FIELD} * {SAGGITALITY_FIELD} * condition * fix_time"
    elif mode == "N":
        model_formula = f"data_mean ~ {LATERALITY_FIELD} * {SAGGITALITY_FIELD} * condition * mean_target_fixation"
    if include_baseline:
        model_formula += " * baseline"
    
    # Fit linear regression model and extract statistics
    model = smf.ols(formula=model_formula, data=filtered_time_bin_data).fit()
    model_summary = summarize_stats_model(model)
    n_predictors = model.model.exog.shape[-1]
    design_matrix_rank = np.linalg.matrix_rank(model.model.exog)
    if design_matrix_rank <= n_predictors:
        logger.error(f"Model (time_bin = {time_bin}) is overparametrized: {n_predictors} predictors but only {design_matrix_rank} observations")
        # return
    obs_tstats = np.array(model_summary["t_value"].to_list())
    coef_names = model_summary["predictor"].to_list()
    
    def permute_data_and_refit_model(df: pd.DataFrame, random_seed: int) -> list:
        # Make copy of dataframe and shuffle the rows of data_mean column
        df = df.copy()
        df["data_mean"] = df["data_mean"].sample(frac=1, random_state=random_seed).to_list()

        # Fit new model on shuffled data
        shuffled_model = smf.ols(formula=model_formula, data=df).fit()

        # Summarized new model
        shuffled_model_summary = summarize_stats_model(shuffled_model)

        # Retrieve t-statistics
        t_stats = shuffled_model_summary["t_value"].to_list()

        return t_stats

    # Permutation for t-statistics
    perm_tstats = np.array(
        [
            permute_data_and_refit_model(filtered_time_bin_data, random_seed=n)
            for n in range(n_permutations)
        ]
    # transpose the resulting matrix so that there are n_permutations rows of n_predictors columns
    ).transpose()

    # Compute empirical (raw) p-values
    assert perm_tstats.shape[0] == obs_tstats.shape[0]  # make sure the number of rows in perm_tstats matches the number of entries in obs_tstats
    p_vals = np.array([
        np.mean(np.abs(perm_tstats[k]) >= np.abs(obs_tstat))
        for k, obs_tstat in enumerate(obs_tstats)
    ])

    # Return summary dataframe of results
    results = pd.DataFrame({
        "time_bin": time_bin,
        "predictor": coef_names,
        "observed_tstat": obs_tstats,
        "permutation_p_value": p_vals,
    })

    return results


def fdr_adjust_pvals(p_values: np.ndarray, alpha: float = 0.05):
    """Run FDR correction (Benjamini-Hochberg) for an array of p-values from permutation testing."""
    _, pvals_corrected, _, _ = multipletests(
        p_values,
        alpha=alpha,
        method='fdr_bh',  # equivalent to R's method = "fdr" in p.adjust
    )
    return pvals_corrected


def block_permutation_test_tstats_fdr(time_windowed_data: pd.DataFrame,
                                      mode: str,
                                      n_permutations: int = 2000,
                                      include_baseline: bool = False,
                                      ) -> pd.DataFrame:
    if mode not in {"FIX", "N"}:
        raise ValueError(f"Expected mode to be either 'FIX' or 'N', got '{mode}'")

    time_bins = time_windowed_data["time_bin"].unique()

    # Run permutations and get raw p-values
    permutation_results = []
    with tqdm(time_bins, unit="bin") as pbar:
        for time_bin in pbar:
            # Update the progress bar to show the current bin
            pbar.set_description(f"Running block permutation test ({n_permutations} permutations) for time_bin={time_bin} ...")

            result = time_bin_permutation(
                time_windowed_data,
                mode=mode,
                time_bin=time_bin,
                n_permutations=n_permutations,
                include_baseline=include_baseline,
            )
            permutation_results.append(result)

    # Combine and apply FDR adjustment to p-values within each time_bin
    permutation_results = pd.concat(permutation_results, axis=0, ignore_index=True)
    permutation_results["fdr_q_value"] = (
        permutation_results
        .groupby("time_bin")["permutation_p_value"]
        .transform(lambda x: fdr_adjust_pvals(x.tolist()))
    )
    return permutation_results


def find_highest_order_significant_predictor_set(df: pd.DataFrame,
                                                 alpha: float = 0.05,
                                                 predictor_sep: str = ":",
                                                 ) -> pd.DataFrame:
    # Initialize highest_order column as False
    df["highest_order"] = False

    # Process each unique time bin
    unique_time_bins = df["time_bin"].unique()
    for time_bin in unique_time_bins:
        # Subset dataframe to only this time bin
        sub_df = df[df["time_bin"] == time_bin]

        # Only consider rows with significant p-values (< alpha) and exclude intercept
        significant_df = (
            sub_df
            .loc[(sub_df["fdr_q_value"] < alpha) & (sub_df["predictor"] != "Intercept")]
        )
        if len(significant_df) == 0:
            logger.info(f"No significant predictors (alpha = {alpha}) found for time bin = {time_bin}")
            continue
        else:
            significant_predictors = significant_df["predictor"].unique()
            n_significant_predictors = len(significant_predictors)
            logger.info(f"{n_significant_predictors} significant predictor(s) (alpha = {alpha}) found for time bin = {time_bin}")
        
        # Add column indicating number of predictor components
        # e.g. with ":" as predictor_sep 'laterality[T.right]:saggitality[T.posterior]:fix_time[T.>1s_before_noun]'
        significant_df["n_predictors"] = significant_df["predictor"].apply(lambda x: len(x.split(predictor_sep)))
        # Sort dataframe by ascending number of predictor components
        significant_df.sort_values(by=['n_predictors'], ascending=True)

        # For each significant predictor, check if a higher-order predictor exists that includes all its parts
        for idx_i, row_i in significant_df.iterrows():
            current_predictor = row_i["predictor"]
            current_predictor_parts = current_predictor.split(predictor_sep)

            found_superset = False
            for _, row_j in significant_df.iloc[idx_i + 1:].iterrows():
                other_predictor = row_j["predictor"]
                other_predictor_parts = other_predictor.split(predictor_sep)
                # Check if all parts of current predictor set are contained in the higher order predictor set
                if len(other_predictor_parts) > len(current_predictor_parts) and all(part in other_predictor_parts for part in current_predictor_parts):
                    found_superset = True
                    break
        
            # Mark as highest order if no superset was found
            row_i["highest_order"] = False if found_superset else True
        
        # Remove no-longer-needed n_predictors column
        significant_df = significant_df.drop("n_predictors", axis="columns")

        # Update original df for rows in this time_bin
        for idx, row in significant_df.iterrows():
            df_mask = (df["predictor"] == row["predictor"]) & (df["time_bin"] == time_bin)
            df.loc[df_mask, "highest_order"] = row["highest_order"]
        
    return df
        

def main(experiment: str | dict | Experiment) -> Experiment:
    # Initialize DGAME experiment from config
    if not isinstance(experiment, Experiment):
        from dgame.dgame import DGAME
        experiment = DGAME.from_input(experiment)
    if len(experiment.subject_ids) < 2:
        logger.warning(f"Fewer than 2 subjects; skipping analysis step {STEP_JA_KEY}")
        return experiment

    # Find per-subject EEG output directories
    per_subject_eeg_outdirs = experiment.get_subject_dirs_dict(experiment.eeg_outdir)
    assert all(len(eeg_outdir) == 1 for eeg_outdir in per_subject_eeg_outdirs.values())
    # Find per-subject unfold_out directories within EEG outdirs
    per_subject_unfold_out_dirs = {
        subject_id: os.path.join(subject_eeg_outdirs[0], "unfold_out")
        for subject_id, subject_eeg_outdirs in per_subject_eeg_outdirs.items()
    }

    # Load fixation data for all subjects
    logger.info(f"Loading EEG fixation data for {len(per_subject_unfold_out_dirs)} subject(s)...")
    all_subjects_fixation_data = load_unfold_out_data(
        per_subject_unfold_out_dirs,
        experiment.channel_coords,
        mode="FIX",
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
        [["data", "subject", "condition", LATERALITY_FIELD, SAGGITALITY_FIELD, "fix_time", "fix_at"]]
        # filter fix_at == "target"
        .loc[lambda df: df["fix_at"] == "target"]
    )
    baseline_aggregated = (
        baseline_data
        .groupby(["subject", "condition", LATERALITY_FIELD, SAGGITALITY_FIELD, "fix_time", "fix_at"], as_index=False)
        .agg(baseline=("data", "mean"))
    )
    # Add time window bins and aggregate with mean
    logger.info("Aggregating fixation data by time window...")
    fixation_data_windows = create_time_windows(all_subjects_fixation_data, mode="FIX")
    # Merge with aggregated_baseline data
    fixation_data_windows = fixation_data_windows.merge(
        baseline_aggregated,
        how="inner",
        on=["subject", "condition", LATERALITY_FIELD, SAGGITALITY_FIELD, "fix_time"]
    )
    # Filter time windows >= 0 and <= 1000
    fixation_data_windows = (
        fixation_data_windows
        .loc[(fixation_data_windows["time_bin"] >= 0) & (fixation_data_windows["time_bin"] <= 1000)]
    )

    # Fit regression models for filtered data by time bin, using all predictors, and collect summaries
    logger.info("Fitting regression models on data per time window...")
    time_bin_model_summaries = summarize_regression_by_time_bins(fixation_data_windows, mode="FIX")
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

    # Run block permutation test
    n_permutations = experiment.get_dgame_step_parameter(STEP_JA_KEY, "n_permutations")
    include_baseline = experiment.get_dgame_step_parameter(STEP_JA_KEY, "include_baseline")
    permutation_results = block_permutation_test_tstats_fdr(
        fixation_data_windows,
        mode="FIX",
        n_permutations=n_permutations,
        include_baseline=include_baseline,
    )
    # Write permutation results to output csv
    permutation_results_outfile = os.path.join(
        experiment.fixations_outdir,
        "permutation_results_fixations.csv"
    )
    permutation_results.to_csv(permutation_results_outfile, index=False)
    logger.info(f"Wrote permutation test results to {permutation_results_outfile}")

    # Find highest order significant permutation results
    alpha = experiment.get_dgame_step_parameter(STEP_JA_KEY, "alpha")
    significant_permutation_results = find_highest_order_significant_predictor_set(permutation_results, alpha=alpha)
    # NB: now has boolean column "highest_order" (True for highest order significant predictor set, else False)

    # Convert to R dataframe and plot in R
    significant_permutation_results_r = convert_pandas2r_dataframe(significant_permutation_results)
    fixation_plot_dir = os.path.join(experiment.fixations_outdir, "plots")
    fixation_plot_outfile = os.path.join(fixation_plot_dir, "fixation-timing.png")
    os.makedirs(fixation_plot_dir, exist_ok=True)
    try:
        create_language_fixation_plot(significant_permutation_results_r, outfile=fixation_plot_outfile)
        logger.info(f"Plotted fixations to {fixation_plot_outfile}")
    except RRuntimeError as exc:
        logger.error(f"Error plotting fixations:\n {exc}")

    return experiment

if __name__ == "__main__":
    parser = argparse.ArgumentParser("Run permutation tests and plot fixation statistics.")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(os.path.abspath(args.config))
