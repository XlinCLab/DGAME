import argparse
import glob
import json
import logging
import os
import pandas as pd
import re
import requests
from tqdm import tqdm

from collections import defaultdict

from constants import AUDIO_FILE_PATTERN, DEFINITE_ARTICLES
from utils import load_config, setdiff, create_timestamp

logger = logging.getLogger(__name__)


def retrieve_word_data_from_corpus(wordlist, corpus="deu_news_2012_3M"):
    # Validate corpus and map to URL
    if corpus == "deu_news_2012_3M":
        url_le = "http://api.wortschatz-leipzig.de/ws/words/deu_news_2012_3M/word/"
    else:
        raise ValueError(f"No URL available for corpus '{corpus}'")

    word_corpus_data = {}
    for word in tqdm(wordlist, desc=f"Retrieving word data from '{corpus}' corpus..."):
        http_le = url_le + word
        response = requests.get(http_le, headers={"Accept": "application/json"})
        retrieved_data = response.json()

        if len(retrieved_data) == 2:
            word_data = {
                "id": None,
                "word": word,
                "freq": None,
                "wordRank": None,
                "frequencyClass": None
            }
        else:
            word_data = retrieved_data

        word_corpus_data[word] = word_data
    return word_corpus_data


def find_subject_audio_files(audio_dir: str) -> list:
    """Find audio files matching the standard audio file pattern in a given audio directory."""
    audio_files = glob.glob(os.path.join(os.path.abspath(audio_dir), "*.csv"))
    # Filter by recognized pattern
    audio_files = [
        audio_file
        for audio_file in audio_files
        if AUDIO_FILE_PATTERN.match(os.path.basename(audio_file))
    ]
    return audio_files


def parse_user_and_block_from_audio_file_name(audio_file: str) -> str:
    """Parses user ID and block ID from an audio file name, e.g. 02_words_11.csv"""
    matched = AUDIO_FILE_PATTERN.match(os.path.basename(audio_file))
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
        usecols=["line", "tmin", "text", "tmax"],
    )
    # Rename columns
    audio_data.columns = ["id", "time", "text", "tmax"]
    # Drop duplicate entries
    audio_data.drop_duplicates(inplace=True)

    # Add column with corpus frequency class for words matching either target objects or filler words
    def retrieve_frequency_class(word):
        if type(word) is float and str(word) == 'nan':
            return None
        if case_insensitive:
            word = word.title()
        if word in corpus_data:
            # Standardize to title casing for corpus query
            word = word.title()
            # Get corpus frequency
            return corpus_data[word]["frequencyClass"]
        return None
    audio_data["frequencyClass"] = audio_data["text"].apply(retrieve_frequency_class)

    # Set missing frequency class entries to the 1 + maximum attested frequency class
    max_freq_class = audio_data["frequencyClass"].max(skipna=True)
    audio_data["frequencyClass"] = audio_data["frequencyClass"].fillna(max_freq_class + 1)

    # Drop rows without any "text" entry
    audio_data = audio_data.dropna(subset=["text"])

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

    words = audio_data["text"].to_list()
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
            pos[idx] = "N"
            pos[idx + 1] = "next"
            pos[idx + 2] = "next"
            # Update counts
            counts[word] += 1
            positions[idx - nback] = counts[word]
            positions[idx] = counts[word]
        # Check if word matches target objects or fillers and assign conditions/codes accordingly
        if word in objects:
            condition_codes[idx - nback] = 11
            condition_codes[idx] = 12
            conditions[idx - nback: idx + 1] = ["conflict"] * (nback + 1)
            pos[idx - 1] = "D"
            pos[idx - 2] = "prev"
        elif word in fillers:
            condition_codes[idx - nback] = 21
            condition_codes[idx] = 22
            conditions[idx - nback: idx + 1] = ["no_conflict"] * (nback + 1)
            pos[idx - nback] = "D"
            pos[idx - (nback + 1)] = "prev"
            pos[idx - (nback + 2)] = "prev"
            if nback == 2:
                pos[idx - 1] = "prev"
    audio_data["condition"] = conditions
    audio_data["condition_code"] = condition_codes
    audio_data["pos"] = pos
    audio_data["position"] = positions

    return audio_data


def combine_words_and_obj_position_data(word_data: pd.DataFrame,
                                        object_positions: pd.DataFrame):
    # Merge object position data
    combined_data = pd.merge(word_data, object_positions, how='left')

    # Iterate again through words
    for idx, row in combined_data.iterrows():
        if not pd.isna(row["surface"]):
            # Check if preceding word is definite article
            if idx > 0 and combined_data["text"][idx - 1] in DEFINITE_ARTICLES:
                nback = 1
            else:
                nback = 2
            combined_data.loc[idx - nback, "surface"] = combined_data.loc[idx, "surface"]
            combined_data.loc[idx - nback, "surface_competitor"] = combined_data.loc[idx, "surface_competitor"]
            combined_data.loc[idx - nback, "surface_end"] = combined_data.loc[idx, "surface_end"]

    # Add other object information to file
    # Get object position entries whose surface_competitor entry is non-NA
    # and take intersection with objects from audio data whose condition is "conflict" and POS == "N"
    targets_lc = set(
        object_positions.loc[object_positions["surface_competitor"].notna(), 'text'].unique()
    ).intersection(
        combined_data.loc[(combined_data["condition"] == "conflict") & (combined_data["pos"] == "N"), "text"].unique()
    )
    targets_lc = list(targets_lc)
    # Get object position entries whose surface_competitor entry is NA
    # and take intersection with objects from audio data whose condition is "no_conflict" and POS == "N"
    fillers_lc = set(
        object_positions.loc[object_positions["surface_competitor"].isna(), 'text'].unique()
    ).intersection(
        combined_data.loc[(combined_data["condition"] == "no_conflict") & (combined_data["pos"] == "N"), "text"].unique()
    )
    fillers_lc = list(fillers_lc)

    pattern = combined_data["pattern"].unique()[0]
    set_id = combined_data["set"].unique()[0]

    def get_surface(text, position, column):
        """Retrieves the surface column value for a dataframe value matching specified text, position, and pattern."""
        result = object_positions[
            (object_positions["text"] == text) &
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
        word = row["text"]
        if pd.isna(row["position"]):
            continue

        # Check if preceding word is a definite article
        if combined_data["text"][idx - 1] in DEFINITE_ARTICLES:
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
    target1 = combined_data[(combined_data["text"] == targets_lc[0]) & (combined_data["pos"] == "N")]
    target1.loc[:, "target_location"] = target1["surface"].shift(-1)
    target1.loc[target1.index[-1], "target_location"] = target1["surface_end"].iloc[0]
    target2 = combined_data[(combined_data["text"] == targets_lc[-1]) & (combined_data["pos"] == "N")]
    target2.loc[:, "target_location"] = target2["surface"].shift(-1)
    target2.loc[target2.index[-1], "target_location"] = target2["surface_end"].iloc[0]
    filler1 = combined_data[(combined_data["text"] == fillers_lc[0]) & (combined_data["pos"] == "N")]
    filler1.loc[:, "target_location"] = filler1["surface"].shift(-1)
    filler1.loc[filler1.index[-1], "target_location"] = filler1["surface_end"].iloc[0]
    filler2 = combined_data[(combined_data["text"] == fillers_lc[-1]) & (combined_data["pos"] == "N")]
    filler2.loc[:, "target_location"] = filler2["surface"].shift(-1)
    filler2.loc[filler2.index[-1], "target_location"] = filler2["surface_end"].iloc[0]
    rest = combined_data[combined_data["pos"] != "N"]

    # Concatenate filtered dataframes back together once end locations are added
    combined_data = pd.concat([target1, target2, filler1, filler2, rest], axis=0, ignore_index=True)

    # Sort dataframe by "id" column to ensure correct (original) order
    combined_data = combined_data.sort_values(by='id').reset_index(drop=True)

    # One final iteration through words
    for idx, row in combined_data.iterrows():
        word = row["text"]
        if not pd.isna(row["target_location"]):
            # Check if preceding word is a definite article
            if combined_data["text"][idx - 1] in DEFINITE_ARTICLES:
                nback = 1
            else:
                nback = 2
            combined_data.loc[combined_data.index[idx - nback], "target_location"] = row["target_location"]

    return combined_data


def load_object_positions_data(filepath, sep=','):
    """Load and preprocess CSV file with object positions data."""
    # Load from CSV file
    obj_pos_data = pd.read_csv(filepath, sep=sep)
    # Rename "object" column to "text" and change to title casing
    obj_pos_data.rename(columns={"object": "text"}, inplace=True)
    obj_pos_data["text"] = obj_pos_data["text"].apply(lambda x: x.title())
    # Drop condition column
    obj_pos_data = obj_pos_data.drop(["condition"], axis=1)
    return obj_pos_data


def main(config_path):
    # Get start time and start timestamp
    start_time, start_timestamp = create_timestamp()

    # Load experiment config
    config = load_config(config_path)
    logger.info(json.dumps(config, indent=4))

    # Use start timestamp as experiment ID if none specified,
    # otherwise combine the specified ID with the timestamp
    experiment_id = config["experiment"].get("id")
    if experiment_id is None or experiment_id.strip() == "":
        experiment_id = start_timestamp
    else:
        experiment_id = os.path.join(experiment_id, start_timestamp)

    # Find audio files
    input_dir = config["data"]["input"]["root"]
    audio_dir = os.path.join(input_dir, config["data"]["input"]["audio"])
    audio_files = find_subject_audio_files(audio_dir)

    # Designate and create output directory
    output_dir = config["experiment"].get("outdir")
    if output_dir is None or output_dir.strip() == "":
        output_dir = "out"
    output_dir = os.path.join(os.path.abspath(output_dir), experiment_id)
    os.makedirs(output_dir, exist_ok=True)

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
        audio_outfile = os.path.join(output_dir, re.sub(r"\.csv$", "analysis.csv", basename))
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Preprocess audio transcript data and combine with object position data.")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(args.config)
