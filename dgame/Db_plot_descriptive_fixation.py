import argparse
import logging
import os

import numpy as np
import pandas as pd
import rpy2.robjects as robjects
from rpy2.rinterface_lib.embedded import RRuntimeError
from rpy2.robjects import FloatVector

from dgame.constants import (AOI_COLUMNS, FIXATION_TIMES_TRIALS_SUFFIX,
                             R_PLOT_SCRIPT_DIR, TRIAL_TIME_OFFSET)
from experiment.load_experiment import Experiment
from experiment.test_subjects import list_subject_files
from utils.r_utils import RDataFrame, convert_pandas2r_dataframe

logger = logging.getLogger(__name__)

# Source R script with custom plotting function
robjects.r["source"](os.path.join(R_PLOT_SCRIPT_DIR, "plot_histogram.R"))
plot_histogram = robjects.globalenv["plot_histogram"]


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
            # Ensure that the subj column is read as a string, e.g. '02' will be read in as an integer
            data = pd.read_csv(fixation_time_trial_file, dtype={"subj": object})
            # Ensure that the subj column has a single entry and matches subject_id
            assert len(data["subj"].unique()) == 1
            assert data["subj"].unique()[0] == subject_id
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


def compute_saccade_angles(df: pd.DataFrame) -> pd.DataFrame:
    """Compute radians and degrees of saccades and add to dataframe."""

    # Compute deltas
    df['dx'] = df['norm_pos_x'].diff()
    df['dy'] = df['norm_pos_y'].diff()

    # Compute angles in radians
    df["angles"] = np.arctan2(df["dy"], df["dx"])

    # Convert to degrees for easier interpretation
    df["angles_deg"] = df["angles"] * (360 / np.pi)

    return df


def get_per_subject_fixation_time_summary(df: pd.DataFrame, subject_ids: list) -> pd.DataFrame:
    per_subject_fixation_time_summary = []
    for subject_id in subject_ids:
        # Filter data per subject
        subj_fixation_data = df[df["subj"] == subject_id]

        # Extract trials
        trials = subj_fixation_data["trial"].unique().tolist()
        for trial in trials:

            # Filter data per trial
            trial_data = subj_fixation_data[subj_fixation_data["trial"] == trial]

            # Extract condition
            # NB: assumes all will be same condition
            assert len(trial_data["condition"].unique()) == 1
            condition = trial_data["condition"].unique().tolist()[0]

            # Filter target and competitor AOI data
            target_data = trial_data[trial_data["aoi_target"] == True]
            comp_data = trial_data[trial_data["aoi_comp"] == True]

            # Target summary
            target_summary = {
                "subj": subject_id,
                "condition": condition,
                "object": "target",
                "dur_min": target_data["duration"].min() if len(target_data) > 0 else np.nan,
                "dur_max": target_data["duration"].max() if len(target_data) > 0 else np.nan,
                "dur_median": target_data["duration"].median() if len(target_data) > 0 else np.nan,
                "tmin": np.nan,
                "tmax": np.nan,
                "tmedian": np.nan
            }
            per_subject_fixation_time_summary.append(target_summary)

            # Competitor summary
            comp_summary = {
                "subj": subject_id,
                "condition": condition,
                "object": "competitor",
                "dur_min": np.nan,
                "dur_max": np.nan,
                "dur_median": np.nan,
                "tmin": comp_data["trial_time"].min() if len(comp_data) > 0 else np.nan,
                "tmax": comp_data["trial_time"].max() if len(comp_data) > 0 else np.nan,
                "tmedian": comp_data["trial_time"].median() if len(comp_data) > 0 else np.nan
            }
            per_subject_fixation_time_summary.append(comp_summary)

    # Transform summaries list into dataframe
    per_subject_fixation_time_summary = pd.DataFrame(per_subject_fixation_time_summary)

    # Drop rows with all NA times
    per_subject_fixation_time_summary = per_subject_fixation_time_summary.dropna(
        how="all",
        subset=["dur_min", "dur_max", "dur_median", "tmin", "tmax", "tmedian"]
    )

    # Convert dataframe from wide to long format
    # TODO check if actually needed
    # long_format_fixation_time_summary = per_subject_fixation_time_summary.melt(
    #     id_vars=["subj", "condition", "object"],
    #     value_vars=["dur_min", "dur_max", "dur_median", "tmin", "tmax", "tmedian"],
    #     var_name="value",
    #     value_name="time"
    # )

    # Mean summary
    # TODO check if these are actually needed for anything; in R script they are computed and nothing is done with them
    # fixes_mean = per_subject_fixation_time_summary.groupby(["subj", "condition", "object"], as_index=False).mean(numeric_only=True)
    # fixes_mean_filtered = fixes_mean.drop(columns=['subj'])
    # fixes_grand_mean = fixes_mean_filtered.groupby('condition', as_index=False).median(numeric_only=True)

    return per_subject_fixation_time_summary


def preprocess_df_for_histogram(df: pd.DataFrame,
                                value_col: str,
                                scale: int | float = None
                                ) -> RDataFrame:
    """
    Prepares a DataFrame for histogram plotting by:
    - selecting relevant columns
    - filtering fix_at to relevant types
    - transforming 'condition' to 'pair' or 'singleton'
    - optionally scaling the value column (e.g., ms to s)
    - dropping NA values
    - converting preprocessed pandas.DataFrame to RDataFrame
    """
    cols = ["subj", value_col, "condition", "fix_at"]
    df_preproc = (
        df[cols]
        .loc[lambda d: d["fix_at"].isin({"comp", "target", "goal"})]
        .assign(
            condition=lambda d: np.where(d["condition"] == "conflict", "pair", "singleton"),
            **({value_col: lambda d: d[value_col] / scale} if scale else {})
        )
        .dropna()
    )

    # Convert to R Dataframe
    df_out = convert_pandas2r_dataframe(df_preproc)

    return df_out


def plot_histograms(data: pd.DataFrame, plot_outdir: str):
    """Preprocesses data and plots 4 histograms: duration, trial time, saccadic amplitude, angular distribution of saccades."""
    os.makedirs(plot_outdir, exist_ok=True)

    # Duration histogram
    duration_histogram_outfile = os.path.join(plot_outdir, "duration_histogram.png")
    try:
        plot_histogram(
            preprocess_df_for_histogram(
                data,
                value_col="duration",
                scale=1000
            ),
            x_var="duration",
            x_label="duration (s)",
            outfile=duration_histogram_outfile,
        )
        logger.info(f"Plotted duration histogram to {duration_histogram_outfile}")
    except RRuntimeError as exc:
        logger.error(f"Error plotting duration histogram:\n {exc}")

    # Trial time histogram
    trial_time_histogram_outfile = os.path.join(plot_outdir, "trial-time_histogram.png")
    try:
        plot_histogram(
            preprocess_df_for_histogram(data, value_col="trial_time"),
            x_var="trial_time",
            x_label="trial time (s)",
            x_limits=FloatVector([-TRIAL_TIME_OFFSET, TRIAL_TIME_OFFSET]),
            x_breaks=FloatVector([-TRIAL_TIME_OFFSET, 0, TRIAL_TIME_OFFSET]),
            outfile=trial_time_histogram_outfile,
        )
        logger.info(f"Plotted trial time histogram to {trial_time_histogram_outfile}")
    except RRuntimeError as exc:
        logger.error(f"Error plotting trial time histogram:\n {exc}")

    # Saccadic amplitude histogram
    saccadicAmpl_histogram_outfile = os.path.join(plot_outdir, "saccadic-amplitude_histogram.png")
    try:
        plot_histogram(
            preprocess_df_for_histogram(data, value_col="saccAmpl"),
            x_var="saccAmpl",
            x_label="normalized saccadic amplitude",
            histogram_position="identity",  # overlay histogram bars rather than positioning them side-by-side
            outfile=saccadicAmpl_histogram_outfile,
        )
        logger.info(f"Plotted saccadic amplitude histogram to {saccadicAmpl_histogram_outfile}")
    except RRuntimeError as exc:
        logger.error(f"Error plotting saccadic amplitude histogram:\n {exc}")

    # Plot angular distribution of saccades
    saccadicAngles_histogram_outfile = os.path.join(plot_outdir, "saccadic-angles_histogram.png")
    try:
        plot_histogram(
            preprocess_df_for_histogram(data, value_col="angles_deg"),
            x_var="angles_deg",
            x_label="Angle (degrees)",
            y_label="Frequency",
            title="Angular Distribution of Saccades",
            fill_label="condition",
            x_limits=FloatVector([0, 360]),
            x_breaks=FloatVector([0, 90, 180, 270]),
            x_minor_breaks=FloatVector(list(range(0, 360, 30))),
            binwidth=10,
            use_polar=True,
            outfile=saccadicAngles_histogram_outfile,
        )
        logger.info(f"Plotted saccadic angular distribution histogram to {saccadicAngles_histogram_outfile}")
    except RRuntimeError as exc:
        logger.error(f"Error plotting angular distribution of saccades:\n {exc}")


def main(experiment: str | dict | Experiment) -> Experiment:
    # Initialize DGAME experiment from config
    if not isinstance(experiment, Experiment):
        from dgame.dgame import DGAME
        experiment = DGAME.from_input(experiment)

    # Get selected subject IDs and per-subject fixation outdirs
    subj_fixation_dirs_dict = experiment.get_subject_dirs_dict(experiment.fixations_outdir)
    subject_ids = subj_fixation_dirs_dict.keys()

    # Iterate through subject IDs, retrieve relevant fixation_times_*_trials.csv files, and combine into single dataframe
    fixation_times_trials_df = load_fixation_times_trials_files(subj_fixation_dirs_dict)

    # Add fixation label to "fix_at" column
    fixation_times_trials_df["fix_at"] = fixation_times_trials_df.apply(determine_fixation_label, axis=1)
    # Convert to category type (equivalent to R's factor)
    fixation_times_trials_df["fix_at"] = fixation_times_trials_df["fix_at"].astype("category")

    # Iterate through subject IDs and get per-subject fixation data
    fixation_times_trials_df["subj"] = fixation_times_trials_df["subj"].astype(str)
    per_subject_fixation_time_summary = get_per_subject_fixation_time_summary(
        fixation_times_trials_df,
        subject_ids = subject_ids,
    )
    summary_outfile = os.path.join(experiment.fixations_outdir, "fixation_time_summary_per_subject.csv")
    per_subject_fixation_time_summary.to_csv(summary_outfile, index=False)
    logger.info(f"Wrote per-subject fixation time summary to {summary_outfile}")

    # Compute radians and degrees of saccades
    fixation_times_trials_df = compute_saccade_angles(fixation_times_trials_df)

    # Plot histograms
    plot_histograms(
        data=fixation_times_trials_df,
        plot_outdir=os.path.join(experiment.fixations_outdir, "plots"),
    )

    return experiment


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Compute and plot descriptive fixation statistics.")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(os.path.abspath(args.config))
