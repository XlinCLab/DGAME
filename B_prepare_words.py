import argparse
import glob
import json
import logging
import os
import pandas as pd
import re
import requests
from tqdm import tqdm

from constants import AUDIO_FILE_PATTERN
from utils import load_config

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
                          audio_outfile: str,
                          corpus_data: dict,
                          objects: set,
                          fillers: set,
                          case_insensitive: bool = True,
                          skip_indices: set = None,
                          sep: str = ","):
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
        if type(word) == float and str(word) == 'nan':
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

    # Initialize "pattern" and "set" columns, both with value=1
    audio_data["pattern"] = 1
    audio_data["set"] = 1

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
    #line_ids = audio_data["id"].to_list()
    conditions = [None] * len(words)
    condition_codes = [None] * len(words)
    pos = [None] * len(words)
    for idx, word in enumerate(words):
        #line_id = int(line_ids[idx])
        #if skip_indices is not None and idx_should_be_skipped(line_id):
        if skip_indices is not None and idx_should_be_skipped(idx):
            continue
        # Check if word matches either target objects or fillers
        if word in objects.union(fillers):
            # Check for preceding definite article
            if idx > 0 and words[idx - 1] in {"die", "der"}: # TODO other definite articles? das, dem, den, ...?
                nback = 1
            else:
                nback = 2
            pos[idx] = "N"
            pos[idx + 1] = "next"
            pos[idx + 2] = "next"
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
            conditions[idx - nback: idx + 1] = ["non_conflict"] * (nback + 1)
            pos[idx - nback] = "D"
            pos[idx - (nback + 1)] = "prev"
            pos[idx - (nback + 2)] = "prev"
            if nback == 2:
                pos[idx - 1] = "prev"
        audio_data["condition"] = conditions
        audio_data["condition_code"] = condition_codes
        audio_data["pos"] = pos

    return audio_data


def main(config):
    # Load experiment config
    config = load_config("config.yml")
    logger.info(json.dumps(config, indent=4))

    # Find audio files
    experiment_root = config["data"]["input"]["root"]
    audio_dir = os.path.join(experiment_root, config["data"]["input"]["audio"])
    audio_files = find_subject_audio_files(audio_dir)
    
    # Initialize target object words and filler words
    # Standardize to title casing (NB: because German nouns are capitalized)
    case_insensitive = config["data"].get("case_insensitive", True)
    if case_insensitive:
        objects = set(obj.title() for obj in config["data"]["objects"])
        fillers = set(filler.title() for filler in config["data"]["fillers"])
    else:
        objects = set(config["data"]["objects"])
        objects = set(config["data"]["fillers"])

    # Fetch word frequency information from corpus for words of interest
    words_of_interest = objects.union(fillers)
    corpus_data = retrieve_word_data_from_corpus(words_of_interest)
    
    # Process audio files
    for audio_file in audio_files:
        logger.info(f"Processing file {audio_file}")
        user_id, block_id = parse_user_and_block_from_audio_file_name(audio_file)
        logger.info(f"User ID: {user_id}")
        logger.info(f"Block ID: {block_id}")
        audio_outfile = re.sub(r"\.csv$", "analysis.csv", audio_file)
        skip_indices = config["data"]["skip_indices"].get(os.path.basename(audio_file))
        preprocess_words_data(
            audio_infile=audio_file,
            audio_outfile=audio_outfile,
            corpus_data=corpus_data,
            objects=objects,
            fillers=fillers,
            skip_indices=skip_indices,
            case_insensitive=case_insensitive,
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser() # TODO add description
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(args.config)
