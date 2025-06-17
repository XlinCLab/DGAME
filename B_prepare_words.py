import argparse
import glob
import json
import logging
import os
import re
import time
from collections import defaultdict
from datetime import timedelta
from typing import Iterable

import pandas as pd
import requests
from tqdm import tqdm

from constants import (AUDIO_FILE_SUFFIX, CONFLICT_LABEL, CORPORA,
                       DEFAULT_CORPUS, DEFINITE_ARTICLES, DET_POS_LABEL,
                       FREQ_CLASS_FIELD, INPUT_LINE_ID_FIELD,
                       INPUT_WORD_ONSET_FIELD, NEXT_WORD_LABEL,
                       NO_CONFLICT_LABEL, NOUN_POS_LABEL, PREV_WORD_LABEL,
                       WORD_END_FIELD, WORD_FIELD, WORD_ID_FIELD,
                       WORD_ONSET_FIELD)
from load_experiment import (create_experiment_outdir, get_experiment_id,
                             list_subject_files, load_config,
                             load_object_positions_data, parse_subject_ids)
from utils import setdiff

logger = logging.getLogger(__name__)


def retrieve_word_data_from_corpus(wordlist: Iterable, corpus: str = DEFAULT_CORPUS) -> dict:
    """Retrieve word data including frequency class from a specified corpus."""
    # Validate corpus and map to URL
    if corpus not in CORPORA:
        raise ValueError(f"No URL available for corpus '{corpus}'. Supported corpora:\n{json.dumps(CORPORA, indent=4)}")
    url_le = CORPORA[corpus]

    word_corpus_data = {}
    for word in tqdm(wordlist, desc=f"Retrieving word data from '{corpus}' corpus..."):
        http_le = url_le + word
        response = requests.get(http_le, headers={"Accept": "application/json"})
        retrieved_data = response.json()

        if len(retrieved_data) == 2:
            word_data = {FREQ_CLASS_FIELD: None}
        else:
            word_data = retrieved_data
            assert FREQ_CLASS_FIELD in word_data

        word_corpus_data[word] = word_data
    return word_corpus_data


def parse_user_and_block_from_audio_file_name(audio_file: str) -> str:
    """Parses user ID and block ID from an audio file name, e.g. 02_words_11.csv"""
    audio_file_pattern = re.compile(rf"(\d+){AUDIO_FILE_SUFFIX}")
    matched = audio_file_pattern.match(os.path.basename(audio_file))
    if not matched:
        raise AssertionError(f"File name '{audio_file}' does not conform to expected format, e.g. 02_words_11.csv")
    user_id, block_id = matched.groups()
    return user_id, block_id


def preprocess_words_data(audio_infile: str,
                          corpus_data: dict,
                          objects: set,
                          fillers: set,
                          case_insensitive: bool = True,
                          skip_indices: set = None,
                          pattern_id: int = 1,
                          set_id: int = 1,
                          sep: str = ",") -> pd.DataFrame:
    # Load audio data CSV into dataframe
    # Explicitly specify columns of interest in case more are present
    audio_data = pd.read_csv(
        audio_infile,
        sep=sep,
        usecols=[INPUT_LINE_ID_FIELD, INPUT_WORD_ONSET_FIELD, WORD_FIELD, WORD_END_FIELD],
    )
    # Rename columns
    audio_data.columns = [WORD_ID_FIELD, WORD_ONSET_FIELD, WORD_FIELD, WORD_END_FIELD]
    # Drop duplicate entries
    audio_data.drop_duplicates(inplace=True)

    # Add column with corpus frequency class for words matching either target objects or filler words
    def retrieve_frequency_class(word):
        if isinstance(word, float) and str(word) == 'nan':
            return None
        if case_insensitive:
            word = word.title()
        if word in corpus_data:
            # Standardize to title casing for corpus query
            word = word.title()
            # Get corpus frequency
            return corpus_data[word][FREQ_CLASS_FIELD]
        return None
    audio_data[FREQ_CLASS_FIELD] = audio_data[WORD_FIELD].apply(retrieve_frequency_class)

    # Set missing frequency class entries to the 1 + maximum attested frequency class
    max_freq_class = audio_data[FREQ_CLASS_FIELD].max(skipna=True)
    audio_data[FREQ_CLASS_FIELD] = audio_data[FREQ_CLASS_FIELD].fillna(max_freq_class + 1)

    # Drop rows without any WORD_FIELD ("text") entry
    audio_data = audio_data.dropna(subset=[WORD_FIELD])

    # Initialize "pattern" and "set" columns
    audio_data["pattern"] = pattern_id
    audio_data["set"] = set_id

    # Iterate over words by line ID, skipping specified indices
    # If word matches one of target object words, check if preceded by definite article
    def idx_should_be_skipped(idx: int) -> bool:
        # List of indices
        if isinstance(skip_indices, list):
            return idx in skip_indices

        # Dict of one or more conditions
        if isinstance(skip_indices, dict):
            if "lte" in skip_indices and idx <= skip_indices["lte"]:
                return True
            if "gte" in skip_indices and idx >= skip_indices["gte"]:
                return True
            if "range" in skip_indices:
                start, end = skip_indices["range"]
                if start <= idx <= end:
                    return True

        # List of multiple condition dicts (e.g., multiple ranges)
        if isinstance(skip_indices, list) and all(isinstance(c, dict) for c in skip_indices):
            for cond in skip_indices:
                if "lte" in cond and idx <= cond["lte"]:
                    return True
                if "gte" in cond and idx >= cond["gte"]:
                    return True
                if "range" in cond:
                    start, end = cond["range"]
                    if start <= idx <= end:
                        return True

        return False

    words = audio_data[WORD_FIELD].to_list()
    conditions = [None] * len(words)
    condition_codes = [None] * len(words)
    pos = [None] * len(words)
    positions = [None] * len(words)
    counts = defaultdict(lambda: 0)
    for idx, word in enumerate(words):
        # line_id = int(line_ids[idx])  #  TODO use line IDs from raw file or index of post-filtered words?
        # if skip_indices is not None and idx_should_be_skipped(line_id):
        if skip_indices is not None and idx_should_be_skipped(idx):
            continue
        # Check if word matches either target objects or fillers
        if word in objects.union(fillers):
            # Check for preceding definite article
            if idx > 0 and words[idx - 1] in DEFINITE_ARTICLES:
                nback = 1
            else:
                nback = 2
            pos[idx] = NOUN_POS_LABEL
            pos[idx + 1] = NEXT_WORD_LABEL
            pos[idx + 2] = NEXT_WORD_LABEL
            # Update counts
            counts[word] += 1
            positions[idx - nback] = counts[word]
            positions[idx] = counts[word]
        # Check if word matches target objects or fillers and assign conditions/codes accordingly
        if word in objects:
            condition_codes[idx - nback] = 11
            condition_codes[idx] = 12
            conditions[idx - nback: idx + 1] = [CONFLICT_LABEL] * (nback + 1)
            pos[idx - 1] = DET_POS_LABEL
            pos[idx - 2] = PREV_WORD_LABEL
        elif word in fillers:
            condition_codes[idx - nback] = 21
            condition_codes[idx] = 22
            conditions[idx - nback: idx + 1] = [NO_CONFLICT_LABEL] * (nback + 1)
            pos[idx - nback] = DET_POS_LABEL
            pos[idx - (nback + 1)] = PREV_WORD_LABEL
            pos[idx - (nback + 2)] = PREV_WORD_LABEL
            if nback == 2:
                pos[idx - 1] = PREV_WORD_LABEL
    audio_data["condition"] = conditions
    audio_data["condition_code"] = condition_codes
    audio_data["pos"] = pos
    audio_data["position"] = positions

    return audio_data


def combine_words_and_obj_position_data(word_data: pd.DataFrame,
                                        object_positions: pd.DataFrame) -> pd.DataFrame:
    # Merge object position data
    combined_data = pd.merge(word_data, object_positions, how='left')

    # Iterate again through words
    for idx, row in combined_data.iterrows():
        if not pd.isna(row["surface"]):
            # Check if preceding word is definite article
            if idx > 0 and combined_data[WORD_FIELD][idx - 1] in DEFINITE_ARTICLES:
                nback = 1
            else:
                nback = 2
            combined_data.loc[idx - nback, "surface"] = combined_data.loc[idx, "surface"]
            combined_data.loc[idx - nback, "surface_competitor"] = combined_data.loc[idx, "surface_competitor"]
            combined_data.loc[idx - nback, "surface_end"] = combined_data.loc[idx, "surface_end"]

    # Add other object information to file
    # Get object position entries whose surface_competitor entry is non-NA
    # and take intersection with objects from audio data whose condition is CONFLICT_LABEL and POS == NOUN_POS_LABEL
    targets_lc = set(
        object_positions.loc[object_positions["surface_competitor"].notna(), WORD_FIELD].unique()
    ).intersection(
        combined_data.loc[(combined_data["condition"] == CONFLICT_LABEL) & (combined_data["pos"] == NOUN_POS_LABEL), WORD_FIELD].unique()
    )
    targets_lc = list(targets_lc)
    # Get object position entries whose surface_competitor entry is NA
    # and take intersection with objects from audio data whose condition is NO_CONFLICT_LABEL and POS == NOUN_POS_LABEL
    fillers_lc = set(
        object_positions.loc[object_positions["surface_competitor"].isna(), WORD_FIELD].unique()
    ).intersection(
        combined_data.loc[(combined_data["condition"] == NO_CONFLICT_LABEL) & (combined_data["pos"] == NOUN_POS_LABEL), WORD_FIELD].unique()
    )
    fillers_lc = list(fillers_lc)

    pattern = combined_data["pattern"].unique()[0]
    set_id = combined_data["set"].unique()[0]

    def get_surface(text: str, position: str, column: str) -> str:
        """Retrieves the surface column value for a dataframe value matching specified text, position, and pattern."""
        result = object_positions[
            (object_positions[WORD_FIELD] == text) &
            (object_positions["position"] == position) &
            (object_positions["pattern"] == pattern) &
            (object_positions["set"] == set_id)
        ][column]
        return result.iloc[0] if not result.empty else None

    where_is_targets = {target: get_surface(target, position=1, column="surface") for target in targets_lc}
    where_is_comps = {target: get_surface(target, position=1, column="surface_competitor") for target in targets_lc}
    where_is_fillers = {filler: get_surface(filler, position=1, column="surface") for filler in fillers_lc}
    # TODO confirm that where_will* are not needed/used
    # where_will_fillers = {filler: get_surface(filler, position=2, column='surface') for filler in fillers_lc}
    # where_will_targets = {target: get_surface(target, position=2, column='surface') for target in targets_lc}

    # Initialize new columns for surfaces/locations
    N = len(combined_data)
    targetA_surface = [pd.NA] * N
    targetB_surface = [pd.NA] * N
    fillerA_surface = [pd.NA] * N
    fillerB_surface = [pd.NA] * N
    compA_surface = [pd.NA] * N
    compB_surface = [pd.NA] * N
    target_location = [pd.NA] * N

    # Iterate again through words in combined dataframe
    for idx, row in combined_data.iterrows():
        word = row[WORD_FIELD]
        if pd.isna(row["position"]):
            continue

        # Check if preceding word is a definite article
        if combined_data[WORD_FIELD][idx - 1] in DEFINITE_ARTICLES:
            nback = 1
        else:
            nback = 2

        if word in targets_lc:
            other_target = other_comp = list(setdiff(targets_lc, {word}))[0]  # TODO why are both the same?  would there ever be >1 item in this result?
            target = comp = list(set(targets_lc).intersection({word}))[0]  # TODO why are both the same? would there ever be >1 item in this result?
            # Set values of new columns
            targetA_surface[idx - nback] = targetA_surface[idx] = where_is_targets[target]
            targetB_surface[idx - nback] = targetB_surface[idx] = where_is_targets[other_target]
            compA_surface[idx - nback] = compA_surface[idx] = where_is_comps[comp]
            compB_surface[idx - nback] = compB_surface[idx] = where_is_comps[other_comp]
            fillerA_surface[idx - nback] = fillerA_surface[idx] = where_is_fillers[fillers_lc[0]]  # TODO/NB: this assumes there are only ever exactly 2 fillers
            fillerB_surface[idx - nback] = fillerB_surface[idx] = where_is_fillers[fillers_lc[-1]]  # TODO/NB: this assumes there are only ever exactly 2 fillers

            # TODO why is this update necessary? seems to be just adding the same values again
            where_is_targets[word] = row["surface"]
            where_is_comps[word] = row["surface_competitor"]
        elif word in fillers_lc:
            other_filler = list(setdiff(fillers_lc, {word}))[0]  # TODO would there ever be >1 item in this result?
            current_filler = list(set(fillers_lc).intersection({word}))[0]  # TODO would there ever be >1 item in this result?
            targetA_surface[idx - nback] = targetA_surface[idx] = where_is_targets[targets_lc[0]]   # TODO/NB: this assumes there are only ever exactly 2 targets
            targetB_surface[idx - nback] = targetB_surface[idx] = where_is_targets[targets_lc[-1]]   # TODO/NB: this assumes there are only ever exactly 2 targets
            compA_surface[idx - nback] = compA_surface[idx] = where_is_comps[targets_lc[0]]
            compB_surface[idx - nback] = compB_surface[idx] = where_is_comps[targets_lc[-1]]
            fillerA_surface[idx - nback] = fillerA_surface[idx] = where_is_fillers[current_filler]
            fillerB_surface[idx - nback] = fillerB_surface[idx] = where_is_fillers[other_filler]

    # Add new surface/location columns to dataframe
    combined_data["targetA_surface"] = targetA_surface
    combined_data["targetB_surface"] = targetB_surface
    combined_data["fillerA_surface"] = fillerA_surface
    combined_data["fillerB_surface"] = fillerB_surface
    combined_data["compA_surface"] = compA_surface
    combined_data["compB_surface"] = compB_surface
    combined_data["target_location"] = target_location

    # Set goal/ending locations
    target1 = combined_data[(combined_data[WORD_FIELD] == targets_lc[0]) & (combined_data["pos"] == NOUN_POS_LABEL)]
    target1.loc[:, "target_location"] = target1["surface"].shift(-1)
    target1.loc[target1.index[-1], "target_location"] = target1["surface_end"].iloc[0]
    target2 = combined_data[(combined_data[WORD_FIELD] == targets_lc[-1]) & (combined_data["pos"] == NOUN_POS_LABEL)]
    target2.loc[:, "target_location"] = target2["surface"].shift(-1)
    target2.loc[target2.index[-1], "target_location"] = target2["surface_end"].iloc[0]
    filler1 = combined_data[(combined_data[WORD_FIELD] == fillers_lc[0]) & (combined_data["pos"] == NOUN_POS_LABEL)]
    filler1.loc[:, "target_location"] = filler1["surface"].shift(-1)
    filler1.loc[filler1.index[-1], "target_location"] = filler1["surface_end"].iloc[0]
    filler2 = combined_data[(combined_data[WORD_FIELD] == fillers_lc[-1]) & (combined_data["pos"] == NOUN_POS_LABEL)]
    filler2.loc[:, "target_location"] = filler2["surface"].shift(-1)
    filler2.loc[filler2.index[-1], "target_location"] = filler2["surface_end"].iloc[0]
    rest = combined_data[combined_data["pos"] != NOUN_POS_LABEL]

    # Concatenate filtered dataframes back together once end locations are added
    combined_data = pd.concat([target1, target2, filler1, filler2, rest], axis=0, ignore_index=True)

    # Sort dataframe by LINE_ID_FIELD ("id") column to ensure correct (original) order
    combined_data = combined_data.sort_values(by=WORD_ID_FIELD).reset_index(drop=True)

    # One final iteration through words
    for idx, row in combined_data.iterrows():
        word = row[WORD_FIELD]
        if not pd.isna(row["target_location"]):
            # Check if preceding word is a definite article
            if combined_data[WORD_FIELD][idx - 1] in DEFINITE_ARTICLES:
                nback = 1
            else:
                nback = 2
            combined_data.loc[combined_data.index[idx - nback], "target_location"] = row["target_location"]

    return combined_data


def main(config: str | dict):
    start_time = time.time()

    # Load experiment config
    if isinstance(config, str):
        config = load_config(config)

    # Retrieve or generate experiment ID
    experiment_id = get_experiment_id(config)

    # Get selected subject IDs
    _, subject_id_regex = parse_subject_ids(config["experiment"]["subjects"])

    # Find audio files
    input_dir = config["data"]["input"]["root"]
    audio_dir = config["data"]["input"]["audio_dir"]
    audio_indir = os.path.join(input_dir, audio_dir)
    audio_files = list_subject_files(dir=audio_indir, subject_regex=subject_id_regex, suffix=AUDIO_FILE_SUFFIX)

    # Designate and create output directory
    output_dir = create_experiment_outdir(config, experiment_id)

    # Initialize target object words and filler words
    # Standardize to title casing (NB: because German nouns are capitalized)
    case_insensitive = config["experiment"].get("case_insensitive", True)
    if case_insensitive:
        objects = set(obj.title() for obj in config["experiment"]["objects"])
        fillers = set(filler.title() for filler in config["experiment"]["fillers"])
    else:
        objects = set(config["experiment"]["objects"])
        objects = set(config["experiment"]["fillers"])

    # Fetch word frequency information from corpus for words of interest
    words_of_interest = objects.union(fillers)
    corpus_data = retrieve_word_data_from_corpus(words_of_interest)

    # Load object positions data
    obj_pos_csv = os.path.join(input_dir, config["data"]["input"]["object_positions"])
    obj_pos_data = load_object_positions_data(obj_pos_csv)

    # Process audio files
    last_subject = None
    new_subject = True
    pattern_id, set_id = 1, 1
    for audio_file in sorted(audio_files):
        logger.info(f"Processing file {audio_file}")
        user_id, block_id = parse_user_and_block_from_audio_file_name(audio_file)
        logger.info(f"User ID: {user_id}")
        logger.info(f"Block ID: {block_id}")
        if user_id == last_subject:
            new_subject = False
            last_subject = user_id
        basename = os.path.basename(audio_file)
        subj_audio_outdir = os.path.join(output_dir, audio_dir, user_id)
        os.makedirs(subj_audio_outdir, exist_ok=True)
        audio_outfile = os.path.join(subj_audio_outdir, re.sub(r"\.csv$", "analysis.csv", basename))
        skip_indices = config["experiment"]["skip_indices"].get(os.path.basename(audio_file))
        word_data = preprocess_words_data(
            audio_infile=audio_file,
            corpus_data=corpus_data,
            objects=objects,
            fillers=fillers,
            case_insensitive=case_insensitive,
            skip_indices=skip_indices,
            pattern_id=pattern_id,
            set_id=set_id,
        )
        combined_data = combine_words_and_obj_position_data(
            word_data=word_data,
            object_positions=obj_pos_data,
        )
        # Write output CSV
        combined_data.to_csv(audio_outfile, index=False)
        logger.info(f"Wrote CSV to {audio_outfile}")
        # Reset pattern and set IDs to 1 if new user
        if new_subject:
            pattern_id = 1
            set_id = 1
        # Otherwise increment/adjust pattern and set IDs for next file from same user
        else:
            if pattern_id == 2:
                pattern_id = 1
                set_id += 1
            else:
                pattern_id += 1

    # Calculate duration of this step and add to run config
    end_time = time.time()
    duration = str(timedelta(seconds=int(end_time - start_time)))
    if "run" in config:
        config["run"]["duration"]["B_prepare_words"] = duration
    logger.info(f"Step B completed successfully (duration: {duration}).")

    return config


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Preprocess audio transcript data and combine with object position data.")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(args.config)
