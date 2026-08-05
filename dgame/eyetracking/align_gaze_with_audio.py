import argparse
import os
import re
from collections import defaultdict

import numpy as np
import pandas as pd
from tqdm import tqdm

from dgame.constants import (CONDITIONS, DIRECTOR_LABEL, ROUND_N,
                             TRIAL_TIME_OFFSET)
from dgame.eyetracking import (AOI_COLUMNS, DEFAULT_CONFIDENCE, ERROR_LABEL,
                               GAZE_CONFIDENCE_FIELD, GAZE_TIMESTAMP_FIELD,
                               SURFACE_COLUMNS, SURFACE_LIST)
from dgame.eyetracking.utils import (load_and_combine_surface_files,
                                     load_gaze_data_from_xdf)
from dgame.paths import (GAZE_POS_SURFACE_SUFFIX, SYNC_REFERENCE_FILE_SUFFIX,
                         WORDS_ANNOTATED_FILE_SUFFIX)
from dgame.words import (NOUN_POS_LABEL, PART_OF_SPEECH_FIELD, WORD_END_FIELD,
                         WORD_FIELD, WORD_ID_FIELD, WORD_ONSET_FIELD)
from dgame.xdf import (AUDIO_STREAM, EYETRACKER_STREAM,
                       SYNC_REFERENCE_DRIFT_INTERCEPT_COLUMN,
                       SYNC_REFERENCE_DRIFT_SLOPE_COLUMN,
                       SYNC_REFERENCE_SYNCED_START_TIME_COLUMN)
from dgame.xdf.utils import (apply_clock_drift_correction,
                             convert_relative_time_to_synced,
                             load_stream_sync_reference)
from experiment.load_experiment import Experiment
from utils.utils import (get_continuous_indices, list_matching_files,
                         merge_dataframes_with_temp_transform, setdiff)


def get_per_subject_audio_and_stream_sync_files(experiment) -> tuple[defaultdict, dict]:
    from dgame.dgame import validate_dgame_input
    experiment = validate_dgame_input(experiment)

    # Get subject IDs and corresponding audio and times directories for each
    subject_audio_dirs = experiment.get_subject_dirs_dict(experiment.audio_outdir)
    subject_time_dirs = experiment.get_subject_dirs_dict(experiment.times_outdir)
    # Ensure that the same subject IDs were found per file type
    try:
        assert set(subject_audio_dirs.keys()) == set(subject_time_dirs.keys())
    except AssertionError as exc:
        raise ValueError("Unequal numbers of subject IDs found!") from exc

    audio_erp_files = {
        subject_id: list_matching_files(
            dir=subject_audio_dir[0],
            pattern=WORDS_ANNOTATED_FILE_SUFFIX,
        )
        for subject_id, subject_audio_dir in subject_audio_dirs.items()
    }
    # Find each subject's single consolidated stream sync-reference file (covers all blocks)
    sync_reference_files = {}
    for subject_id, subject_times_dir in subject_time_dirs.items():
        matches = list_matching_files(dir=subject_times_dir[0], pattern=SYNC_REFERENCE_FILE_SUFFIX)
        try:
            assert len(matches) == 1
        except AssertionError as exc:
            raise ValueError(f"Expected exactly 1 sync-reference file for subject ID={subject_id}, found {len(matches)}") from exc
        sync_reference_files[subject_id] = matches[0]

    return audio_erp_files, sync_reference_files


def load_erp_file(erp_file: str) -> pd.DataFrame:
    """Load and preprocess ERP CSV file."""
    # Load ERP file data
    erp_file_data = pd.read_csv(erp_file)

    # Make sure second column is labeled as "time"
    if erp_file_data.columns[1] != WORD_ONSET_FIELD:
        renamed_columns = erp_file_data.columns.tolist()
        renamed_columns[1] = WORD_ONSET_FIELD
        erp_file_data.columns = renamed_columns

    return erp_file_data


def align_times_to_erp_word_timings(times: np.ndarray,
                                    erp_time_ids: dict,
                                    ) -> dict:
    """Align an array of timestamps with the nearest timestamp associated with a word ID.
    Returns a dictionary of time: aligned_time key-value pairs."""
    # Sort times for binary search
    sorted_times = np.sort(times)
    n_times = len(sorted_times)

    # For each new time, find nearest ERP time associated with a word ID
    # Initialize aligned word ID as 0, then replace with actual ID if aligned with a word timestamp
    word_aligned_times = {time_i: 0 for time_i in times}
    for time_i, word_id in erp_time_ids.items():
        # Find insertion position in sorted time list
        idx = np.searchsorted(sorted_times, time_i)

        # Determine neighbors
        time_before = sorted_times[idx - 1] if idx > 0 else None
        # NB: Because using insertion index, would be idx rather than idx + 1
        time_after = sorted_times[idx] if idx < n_times else None

        # Choose the nearer time
        if time_before is None:
            word_aligned_times[time_after] = word_id
        elif time_after is None:
            word_aligned_times[time_before] = word_id
        else:
            before_diff = abs(time_i - time_before)
            after_diff = abs(time_i - time_after)
            # if before_diff == after_diff:
            #     logger.warning("Equal distance to previous and following ERP timestamps, defaulting to previous")
            aligned_time = time_before if before_diff <= after_diff else time_after
            word_aligned_times[aligned_time] = word_id

    return word_aligned_times


def align_subject_gaze_data_with_audio(erp_file: str,
                                       xdf_file: str,
                                       sync_reference: pd.DataFrame,
                                       block: int,
                                       gaze_positions_subj: pd.DataFrame,
                                       words_df: pd.DataFrame,
                                       ) -> pd.DataFrame:
    # Load ERP file data
    erp_file_data = load_erp_file(erp_file)

    # Load this block's gaze data directly from its own XDF eyetracker stream --
    # gaze_timestamp is on the stream's raw (un-synchronized) local clock
    gaze_data = load_gaze_data_from_xdf(xdf_file)

    # Look up this block's linear clock-drift correction (intercept and slope) per stream
    # (synced = raw + drift_intercept + drift_slope * raw)
    # so that audio-relative word onsets (erp_file_data) and eyetracker gaze times (gaze_data)
    # can be put on the same absolute LSL axis before being compared/matched
    audio_synced_start_time = sync_reference.loc[(block, AUDIO_STREAM), SYNC_REFERENCE_SYNCED_START_TIME_COLUMN]
    audio_drift_slope = sync_reference.loc[(block, AUDIO_STREAM), SYNC_REFERENCE_DRIFT_SLOPE_COLUMN]
    eyetracker_drift_intercept = sync_reference.loc[(block, EYETRACKER_STREAM), SYNC_REFERENCE_DRIFT_INTERCEPT_COLUMN]
    eyetracker_drift_slope = sync_reference.loc[(block, EYETRACKER_STREAM), SYNC_REFERENCE_DRIFT_SLOPE_COLUMN]

    # Convert eyetracker raw-clock gaze times to the shared absolute LSL axis
    # and add as new column "time" to gaze_data dataframe
    gaze_times = apply_clock_drift_correction(
        raw_time=gaze_data[GAZE_TIMESTAMP_FIELD],
        drift_intercept=eyetracker_drift_intercept,
        drift_slope=eyetracker_drift_slope,
    ).to_numpy()
    gaze_data[WORD_ONSET_FIELD] = gaze_times
    # Record which physical block these rows came from for every row
    gaze_data["block"] = block

    # Convert audio-relative word onset/end times to the same shared absolute LSL axis:
    # erp_file_data times are seconds elapsed since the audio stream's first sample,
    # thus are anchored directly on the stream's synced start time and scaled by its drift slope
    # NB: NOT shifted by the stream's raw start time, which for a regular-rate stream like audio
    # can differ from the raw time fit_drift_correction was actually fit against, since 
    # pyxdf dejitters regular-rate raw timestamps before applying its own clock drift correction
    erp_file_data[WORD_ONSET_FIELD] = convert_relative_time_to_synced(
        relative_time=erp_file_data[WORD_ONSET_FIELD].astype(float),
        synced_start_time=audio_synced_start_time,
        drift_slope=audio_drift_slope,
    )
    if WORD_END_FIELD in erp_file_data.columns:
        erp_file_data[WORD_END_FIELD] = convert_relative_time_to_synced(
            relative_time=erp_file_data[WORD_END_FIELD].astype(float),
            synced_start_time=audio_synced_start_time,
            drift_slope=audio_drift_slope,
        )

    # Extract known ERP times and word IDs into dict mapping
    erp_times = np.array(erp_file_data[WORD_ONSET_FIELD], dtype=float)
    erp_word_ids = np.array(erp_file_data[WORD_ID_FIELD])
    erp_time_ids = dict(zip(erp_times, erp_word_ids))

    # Align each time to the nearest ERP time associated with a word ID
    word_aligned_times = align_times_to_erp_word_timings(gaze_times, erp_time_ids)
    # Add aligned word IDs to gaze_data dataframe
    gaze_data[WORD_ID_FIELD] = word_aligned_times.values()

    # Create temp copy of erp_file_data, renaming "time" to "audio_time"
    tmp_erp_file_data = erp_file_data.copy().rename(columns={WORD_ONSET_FIELD: "audio_time"})
    # Merge tmp_erp_file_data and gaze_data by word "id" column
    # Now there should be an "audio_time" column as well as "time" column
    gaze_data = gaze_data.merge(tmp_erp_file_data, on=WORD_ID_FIELD, how='left')

    # Lowercase text/word field of gaze_data
    gaze_data[WORD_FIELD] = gaze_data[WORD_FIELD].str.lower()

    # Add gaze_data to gaze_positions_subj dataframe
    gaze_positions_subj = pd.concat([gaze_positions_subj, gaze_data], axis=0, ignore_index=True)

    # Add erp_file_data to words_df
    words_df = pd.concat([words_df, erp_file_data], axis=0, ignore_index=True)

    return gaze_positions_subj, words_df


def add_trials_to_gaze_data(gaze_positions_subj: pd.DataFrame) -> pd.DataFrame:
    """Adds trials within time window to per-subject dataframe."""
    gaze_positions_subj["trial"] = pd.NA
    gaze_positions_subj["trial_time"] = pd.NA
    trial = 1
    noun_row_indices = gaze_positions_subj.index[
        gaze_positions_subj["condition"].isin(CONDITIONS) &
        (gaze_positions_subj[PART_OF_SPEECH_FIELD] == NOUN_POS_LABEL)
    ].tolist()

    # Iterate over noun row indices and add trial annotations for data points within trial time window
    progress_bar_n = len(noun_row_indices)
    with tqdm(total=progress_bar_n) as pbar:
        pbar.set_description("Adding trials...")
        for idx in noun_row_indices:
            row = gaze_positions_subj.loc[idx]
            pattern_id = row["pattern"]
            set_id = row["set"]
            time_at_idx = row[WORD_ONSET_FIELD]
            trial_start_time = time_at_idx - TRIAL_TIME_OFFSET
            trial_end_time = time_at_idx + TRIAL_TIME_OFFSET
            gaze_positions_subj.loc[idx, "trial"] = trial
            gaze_positions_subj.loc[idx, "trial_time"] = 0

            # Identify rows within trial time frame and within current block (coded by pattern and set columns)
            pre_indices_within_trial = gaze_positions_subj.index[
                (
                    (gaze_positions_subj["pattern"] == pattern_id) |
                    # Need to include NA as possible value for pattern and set columns and then apply secondary filter (see explanation below)
                    # Otherwise the condition and surface columns will not be set for those rows
                    (gaze_positions_subj["pattern"].isna())
                ) &
                (
                    (gaze_positions_subj["set"] == set_id) |
                    (gaze_positions_subj["set"].isna())
                ) &
                (gaze_positions_subj[WORD_ONSET_FIELD] > trial_start_time) &
                (gaze_positions_subj[WORD_ONSET_FIELD] < time_at_idx) &
                (abs(gaze_positions_subj[WORD_ONSET_FIELD] - time_at_idx) <= TRIAL_TIME_OFFSET)
            ].tolist()
            post_indices_within_trial = gaze_positions_subj.index[
                (
                    (gaze_positions_subj["pattern"] == pattern_id) |
                    (gaze_positions_subj["pattern"].isna())
                ) &
                (
                    (gaze_positions_subj["set"] == set_id) |
                    (gaze_positions_subj["set"].isna())
                ) &
                (gaze_positions_subj[WORD_ONSET_FIELD] < trial_end_time) &
                (gaze_positions_subj[WORD_ONSET_FIELD] > time_at_idx) &
                (abs(gaze_positions_subj[WORD_ONSET_FIELD] - time_at_idx) <= TRIAL_TIME_OFFSET)
            ].tolist()
            # Above filtering logic does not fully exclude indices from different blocks, which have NA values for "pattern" and "set" fieds
            # Solution is to take only the longest continuous set of indices around the central time point,
            # which would mimic forward and backward while loop behavior.
            # When a non-continuous index is found, it necessarily comes from different block
            # Exclude all indices which are non-contiguous to the central set
            pre_indices_within_trial = get_continuous_indices(idx, pre_indices_within_trial, direction="pre")
            post_indices_within_trial = get_continuous_indices(idx, post_indices_within_trial, direction="post")
            indices_to_set = pre_indices_within_trial + post_indices_within_trial

            # Set column values for data points within trial window
            gaze_positions_subj.loc[indices_to_set, "trial"] = trial
            gaze_positions_subj.loc[indices_to_set, "trial_time"] = gaze_positions_subj.loc[indices_to_set, WORD_ONSET_FIELD] - time_at_idx
            for col_name in SURFACE_COLUMNS + ["condition"]:
                gaze_positions_subj.loc[indices_to_set, col_name] = row[col_name]

            # Increment trial
            trial += 1

            # Update progress bar
            pbar.update(1)

    return gaze_positions_subj


def validate_surface_annotation(value: str | int) -> str:
    """Validates that a surface annotation value belongs to the list of surfaces.
    Corrects single-digit surfaces to corresponding double-digit surfaces, e.g. "2" -> "22"."""
    if pd.isna(value):
        return value
    value = str(int(value))
    if len(value) == 1:
        value = value * 2
    try:
        assert value in SURFACE_LIST
    except AssertionError as exc:
        raise ValueError(f"Invalid surface annotation: {value}") from exc
    return value


def add_surface_aoi_annotations(gaze_positions_subj: pd.DataFrame) -> pd.DataFrame:
    """Add annotations for surface areas of interest for a single subject."""
    surface_condition_indices = gaze_positions_subj.index[
        (gaze_positions_subj["trial_time"].notna()) &
        (gaze_positions_subj["condition"].isin(CONDITIONS)) &
        (gaze_positions_subj["surface"].notna())
    ].tolist()

    # Cast AOI columns to nullable Boolean type to avoid warnings when setting NaN values
    # NB: otherwise will cause an error in a future version of pandas
    for aoi_field in AOI_COLUMNS.keys():
        gaze_positions_subj[aoi_field] = gaze_positions_subj[aoi_field].astype('boolean')

    progress_bar_n = len(surface_condition_indices)
    with tqdm(total=progress_bar_n) as pbar:
        pbar.set_description("Annotating surface areas of interest...")
        for idx in surface_condition_indices:
            row = gaze_positions_subj.loc[idx]
            target, goal, otherTarget, fillerA, fillerB, competitor, otherCompetitor = map(
                validate_surface_annotation,
                [
                    row["surface"],
                    row["target_location"],
                    row["targetB_surface"],
                    row["fillerA_surface"],
                    row["fillerB_surface"],
                    row["surface_competitor"],
                    row["compB_surface"],
                ]
            )
            empty_surfaces = setdiff(SURFACE_LIST, {target, competitor, otherTarget, otherCompetitor, fillerA, fillerB})

            # Add area of interest (AOI) flags
            def set_aoi_flag(surface, aoi_field):
                if surface != ERROR_LABEL and pd.notna(surface):
                    gaze_positions_subj.at[idx, aoi_field] = row[surface]  # noqa: B023

            set_aoi_flag(target, 'aoi_target')
            set_aoi_flag(goal, 'aoi_goal')
            set_aoi_flag(otherTarget, 'aoi_otherTarget')
            set_aoi_flag(competitor, 'aoi_comp')
            set_aoi_flag(otherCompetitor, 'aoi_otherComp')
            set_aoi_flag(fillerA, 'aoi_fillerA')
            set_aoi_flag(fillerB, 'aoi_fillerB')

            # Mark whether any AOI surface is empty
            gaze_positions_subj.at[idx, "aoi_empty"] = any(
                pd.notna(surface) and gaze_positions_subj.at[idx, surface] is True
                for surface in empty_surfaces
            )

            # Update progress bar
            pbar.update(1)

    return gaze_positions_subj


def main(experiment: str | dict | Experiment) -> Experiment:

    # Initialize DGAME experiment from config
    if not isinstance(experiment, Experiment):
        from dgame.dgame import DGAME
        experiment = DGAME.from_input(experiment)
    logger = experiment.logger

    # Find per-subject audio ERP and stream sync-reference files
    logger.info("Loading per-subject audio and timing files...")
    subj_audio_erp_dict, subj_sync_reference_dict = get_per_subject_audio_and_stream_sync_files(experiment)
    # Get subject IDs (should be identical for all file types)
    subject_ids = sorted(list(subj_audio_erp_dict.keys()))
    logger.info(f"Processing {len(subject_ids)} subject ID(s): {', '.join(subject_ids)}")

    # Iterate over subjects and combine word data with gaze data
    gaze_positions_all = pd.DataFrame()
    for subject_id in subject_ids:
        logger.info(f"Processing subject '{subject_id}'...")

        # Load surface and object position data
        logger.info("Loading surface fixation position data...")
        surface_files = list_matching_files(
            dir=os.path.join(experiment.surface_indir, subject_id),
            pattern=GAZE_POS_SURFACE_SUFFIX,
        )
        surface_pos_data = load_and_combine_surface_files(surface_files)
        # Round gaze_timestamp field in order to enable merge
        surface_pos_data[f"rounded_{GAZE_TIMESTAMP_FIELD}"] = round(surface_pos_data[GAZE_TIMESTAMP_FIELD], ROUND_N)

        # Create per-subject subject output directories
        for outdir_i in {experiment.audio_outdir, experiment.gaze_outdir}:
            subj_outdir_i = os.path.join(outdir_i, subject_id)
            os.makedirs(subj_outdir_i, exist_ok=True)
        # Designate per-subject output file paths
        word_outfile = os.path.join(experiment.audio_outdir, subject_id, "all_words_4analysis.csv")
        gaze_before_words_file = os.path.join(experiment.gaze_outdir, subject_id, "gaze_positions_before_words.csv")
        gaze_subj_out = os.path.join(experiment.gaze_outdir, subject_id, "gaze_positions_4analysis.csv")
        tmp_gaze_s = os.path.join(experiment.gaze_outdir, subject_id, "tmp_gaze_positions.csv")

        logger.info("Loading word data and combining with gaze data...")
        # Initialize empty dataframe to contain all processed gaze data per subject
        gaze_positions_subj = pd.DataFrame()
        words_df = pd.DataFrame()

        # Load this subject's consolidated stream sync-reference data (covers all blocks)
        sync_reference_file = subj_sync_reference_dict[subject_id]
        logger.debug(f"Sync reference file: {os.path.basename(sync_reference_file)}")
        sync_reference = load_stream_sync_reference(sync_reference_file)

        # Load word data and combine with gaze data
        audio_erp_files = subj_audio_erp_dict[subject_id]
        for erp_file in audio_erp_files:
            logger.debug(f"ERP file: {os.path.basename(erp_file)}")
            block = int(re.search(WORDS_ANNOTATED_FILE_SUFFIX, os.path.basename(erp_file)).group(1))
            xdf_file = experiment.get_xdf_file(subject_id, block, role=DIRECTOR_LABEL)
            gaze_positions_subj, words_df = align_subject_gaze_data_with_audio(
                erp_file=erp_file,
                xdf_file=xdf_file,
                sync_reference=sync_reference,
                block=block,
                gaze_positions_subj=gaze_positions_subj,
                words_df=words_df,
            )

        # Merge gaze positions and surface positions by timestamp
        # Add another column with rounded timestamps in order to merge, floating point timestamps may not match exactly
        gaze_positions_subj = merge_dataframes_with_temp_transform(
            left_df=gaze_positions_subj,
            right_df=surface_pos_data,
            on=GAZE_TIMESTAMP_FIELD,
            how="left",
            transform=lambda x: round(x, ROUND_N),
            transform_left=True,
            # NB: not necessary, to transform right df
            # surface_pos_data is transformed only once in advance to avoid re-transforming in each subject loop iteration
            transform_right=False,
            temp_column_name=f"rounded_{GAZE_TIMESTAMP_FIELD}",
        )

        # Add subject ID to gaze_positions_subj dataframe and write CSV outfiles
        gaze_positions_subj["subj"] = subject_id
        gaze_positions_subj["subj"] = gaze_positions_subj["subj"].astype(str)
        gaze_positions_subj.to_csv(gaze_before_words_file, index=False)
        logger.info(f"Wrote subject {subject_id} data to {gaze_before_words_file}")
        words_df.to_csv(word_outfile, index=False)
        logger.info(f"Wrote word data to {word_outfile}")

        # Add trials
        gaze_positions_subj = add_trials_to_gaze_data(gaze_positions_subj)
        # Write CSV file with added trials
        gaze_positions_subj.to_csv(tmp_gaze_s, index=False)
        logger.info(f"Wrote data with trial annotations to {tmp_gaze_s}")

        # Initialize new area of interest (aoi) columns as False
        for new_column in AOI_COLUMNS.keys():
            gaze_positions_subj[new_column] = False
        # Set trackloss column to boolean value, whether confidence < DEFAULT_CONFIDENCE
        gaze_positions_subj["trackloss"] = gaze_positions_subj[GAZE_CONFIDENCE_FIELD] < DEFAULT_CONFIDENCE

        # Check if participants looked at a surface or not at a given time point
        logger.info("Annotating surface areas of interest...")
        gaze_positions_subj = add_surface_aoi_annotations(gaze_positions_subj)

        # Sort dataframe by gaze_timestamp
        gaze_positions_subj = gaze_positions_subj.sort_values(by=GAZE_TIMESTAMP_FIELD, ascending=True).reset_index(drop=True)

        # Write CSV file
        gaze_positions_subj.to_csv(gaze_subj_out, index=False)
        logger.info(f"Wrote per-subject gaze file (subject = {subject_id}) to {gaze_subj_out}")

        # Add per-subject trial data into running dataframe for all subjects
        gaze_positions_all = pd.concat([gaze_positions_all, gaze_positions_subj], axis=0, ignore_index=True)

    # Write full gaze positions CSV file for all subjects
    gaze_all_out = os.path.join(experiment.gaze_outdir, "gaze_positions_all_4analysis.csv")
    gaze_positions_all.to_csv(gaze_all_out, index=False)
    logger.info(f"Wrote full gaze file (all subjects) to {gaze_all_out}")

    return experiment


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Creates trials from the continuous gaze file.")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(os.path.abspath(args.config))
