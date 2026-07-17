import os

import pandas as pd

from dgame.constants import (COLUMN_DATA_TYPES, ROUND_N, WORD_END_FIELD,
                             WORD_ONSET_FIELD)
from experiment.load_experiment import Experiment


def load_filtered_gaze_data(experiment: Experiment) -> pd.DataFrame:
    """Load the aggregate gaze positions file (output of step Ca) and filter to valid trial
    data points with non-negative word durations.

    Shared by Da_compute_trialtime and Da_gaze_language_timecourse_stats, which both derive
    their input from this same filtered dataframe but otherwise have no dependency on one another.
    """
    # TODO confirm that this should use the aggregated file, not individual per subject
    gaze_infile = os.path.join(experiment.gaze_outdir, "gaze_positions_all_4analysis.csv")

    # Ensure that the subj column is read as a string; e.g. '02' will be read in as an integer unless otherwise specified in dtype arg
    gaze_data = pd.read_csv(gaze_infile, dtype=COLUMN_DATA_TYPES)
    gaze_data = gaze_data.loc[
        gaze_data["condition"].notna() &
        gaze_data["trial_time"].notna() &
        (~gaze_data["trackloss"]) &
        gaze_data["subj"].isin(experiment.subject_ids)
    ].drop_duplicates()

    # Add column "duration": tmax - time (rounded to ROUND_N digits)
    gaze_data["duration"] = (gaze_data[WORD_END_FIELD] - gaze_data[WORD_ONSET_FIELD]).round(ROUND_N)

    # Filter valid durations and then drop duration column
    gaze_data = gaze_data[gaze_data["duration"] >= 0].drop(columns="duration")

    return gaze_data
