import argparse
import logging
import os

import rpy2.robjects as robjects
from rpy2.rinterface_lib.embedded import RRuntimeError

from dgame.constants import (LATERALITY_FIELD, R_PLOT_SCRIPT_DIR,
                             SAGGITALITY_FIELD, STEP_JB_KEY)
from dgame.Ja_lm_permute_and_plot_fixations import (
    block_permutation_test_tstats_fdr, create_time_windows,
    find_highest_order_significant_predictor_set, load_unfold_out_data,
    summarize_regression_by_time_bins)
from experiment.load_experiment import Experiment
from utils.r_utils import convert_pandas2r_dataframe

logger = logging.getLogger(__name__)


# Source R script with custom plotting function
robjects.r["source"](os.path.join(R_PLOT_SCRIPT_DIR, "plot_language_fixation_stats.R"))
create_language_fixation_plot = robjects.globalenv["create_language_fixations_plot"]
        

def main(experiment: str | dict | Experiment) -> Experiment:
    # Initialize DGAME experiment from config
    if not isinstance(experiment, Experiment):
        from dgame.dgame import DGAME
        experiment = DGAME.from_input(experiment)

    # Find per-subject EEG output directories
    per_subject_eeg_outdirs = experiment.get_subject_dirs_dict(experiment.eeg_outdir)
    assert all(len(eeg_outdir) == 1 for eeg_outdir in per_subject_eeg_outdirs.values())
    # Find per-subject unfold_out directories within EEG outdirs
    per_subject_unfold_out_dirs = {
        subject_id: os.path.join(subject_eeg_outdirs[0], "unfold_out")
        for subject_id, subject_eeg_outdirs in per_subject_eeg_outdirs.items()
    }

    # Load noun data for all subjects
    logger.info(f"Loading EEG language data for {len(per_subject_unfold_out_dirs)} subject(s)...")
    all_subjects_noun_data = load_unfold_out_data(
        per_subject_unfold_out_dirs,
        experiment.channel_coords,
        mode="N",
    )
    
    # Compute baseline factor  # TODO baseline of what?
    # Summarize "data" field as "baseline" for every combination of subject/condition/laterality/saggitality/mean_target_fixation
    # from a filtered subset of the data containing only time >= -250 but < 0
    logger.info("Filtering noun data to establish baseline...")
    baseline_data = (
        all_subjects_noun_data
        # filter time >= -250 & time < 0
        .loc[(all_subjects_noun_data["time"] >= -250) & (all_subjects_noun_data["time"] < 0)]
        # select only these columns
        [["data", "subject", "condition", LATERALITY_FIELD, SAGGITALITY_FIELD, "mean_target_fixation"]]
    )
    baseline_aggregated = (
        baseline_data
        .groupby(["subject", "condition", LATERALITY_FIELD, SAGGITALITY_FIELD, "mean_target_fixation"], as_index=False)
        .agg(baseline=("data", "mean"))
    )
    # Add time window bins and aggregate with mean
    logger.info("Aggregating noun data by time window...")
    noun_data_windows = create_time_windows(all_subjects_noun_data, mode="N")
    # Merge with aggregated_baseline data
    noun_data_windows = noun_data_windows.merge(
        baseline_aggregated,
        how="inner",
        on=["subject", "condition", LATERALITY_FIELD, SAGGITALITY_FIELD, "mean_target_fixation"]
    )
    # Filter time windows >= 0 and <= 1000
    noun_data_windows = (
        noun_data_windows
        .loc[(noun_data_windows["time_bin"] >= 0) & (noun_data_windows["time_bin"] <= 1000)]
    )

    # Fit regression models for filtered data by time bin, using all predictors, and collect summaries
    logger.info("Fitting regression models on data per time window...")
    time_bin_model_summaries = summarize_regression_by_time_bins(noun_data_windows, mode="N")
    if len(time_bin_model_summaries) > 0:
        # Write model summaries to output csv
        regression_summary_outfile = os.path.join(
            experiment.fixations_outdir,
            "regression_results_nouns.csv",
        )
        time_bin_model_summaries.to_csv(regression_summary_outfile, index=False)
        logger.info(f"Wrote regression model summaries to {regression_summary_outfile}")
    else:
        logger.warning("No valid regression models were produced.")

    # Run block permutation test
    n_permutations = experiment.get_dgame_step_parameter(STEP_JB_KEY, "n_permutations")
    include_baseline = experiment.get_dgame_step_parameter(STEP_JB_KEY, "include_baseline")
    permutation_results = block_permutation_test_tstats_fdr(
        noun_data_windows,
        mode="N",
        n_permutations=n_permutations,
        include_baseline=include_baseline,
    )
    # Write permutation results to output csv
    permutation_results_outfile = os.path.join(
        experiment.fixations_outdir,
        "permutation_results_nouns.csv"
    )
    permutation_results.to_csv(permutation_results_outfile, index=False)
    logger.info(f"Wrote permutation test results to {permutation_results_outfile}")

    # Find highest order significant permutation results
    alpha = experiment.get_dgame_step_parameter(STEP_JB_KEY, "alpha")
    significant_permutation_results = find_highest_order_significant_predictor_set(permutation_results, alpha=alpha)
    # NB: now has boolean column "highest_order" (True for highest order significant predictor set, else False)

    # Convert to R dataframe and plot in R
    significant_permutation_results_r = convert_pandas2r_dataframe(significant_permutation_results)
    fixation_plot_dir = os.path.join(experiment.fixations_outdir, "plots")
    language_plot_outfile = os.path.join(fixation_plot_dir, "language-stats.png")
    os.makedirs(fixation_plot_dir, exist_ok=True)
    try:
        create_language_fixation_plot(significant_permutation_results_r, outfile=language_plot_outfile)
        logger.info(f"Plotted language statistics to {language_plot_outfile}")
    except RRuntimeError as exc:
        logger.error(f"Error plotting language statistics:\n {exc}")

    return experiment

if __name__ == "__main__":
    parser = argparse.ArgumentParser("Run permutation tests and plot language statistics.")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(os.path.abspath(args.config))
