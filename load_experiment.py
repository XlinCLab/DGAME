import glob
import logging
import os
import re
from collections import defaultdict
from typing import Iterable

import pandas as pd
import yaml

from constants import OBJECT_FIELD, WORD_FIELD
from utils import create_timestamp

logger = logging.getLogger(__name__)


def recursively_inherit_dict_values(target: dict, source: dict) -> None:
    for key, value in source.items():
        if target is not None and key not in target:
            target[key] = value
        elif isinstance(value, dict):
            recursively_inherit_dict_values(target[key], value)


def to_absolute_path(path: str) -> str:
    if os.path.isabs(path):
        return path
    current_file_path = os.path.abspath(__file__)
    root_path = os.path.abspath(os.path.join(current_file_path, os.pardir))
    return os.path.join(root_path, path)


def load_config(config_path: str) -> dict:
    """Returns a dictionary containing parameters from a specified config.yml file

    Args:
        config_path (str): Path to config.yml file

    Returns:
        config: nested dictionary of parameter names and values
    """

    config_path = to_absolute_path(config_path)
    with open(config_path, 'r') as f:
        loaded_config: dict = yaml.safe_load(f)

    included_config_block: dict | None = loaded_config.get('config')
    if included_config_block is None:
        return loaded_config
    included_config_path: str | None = included_config_block.get('extend')
    if included_config_path is None:
        return loaded_config

    included_config = None
    if included_config_path is not None:
        included_config_path: str = to_absolute_path(included_config_path)
        included_config = load_config(included_config_path)

    recursively_inherit_dict_values(loaded_config, included_config)

    return loaded_config


def get_experiment_id(config: dict) -> str:
    """Retrieve experiment ID from config (if set) and combine with timestamp."""
    experiment_id = config["experiment"].get("id")
    _, timestamp = create_timestamp()
    if experiment_id is None or experiment_id.strip() == "":
        experiment_id = timestamp
    else:
        experiment_id = os.path.join(experiment_id, timestamp)
    return experiment_id


def create_experiment_outdir(config: dict, experiment_id: str = None) -> str:
    """Retrieve and create experiment output directory."""
    output_dir = config["experiment"].get("outdir")
    if output_dir is None or output_dir.strip() == "":
        output_dir = "out"
    experiment_id = get_experiment_id(config) if experiment_id is None else experiment_id
    output_dir = os.path.join(os.path.abspath(output_dir), experiment_id)
    os.makedirs(output_dir, exist_ok=True)
    return output_dir


def parse_subject_ids(subject_ids: Iterable | str | int | None) -> tuple[list, str]:
    if len(subject_ids) == 0 or subject_ids is None or isinstance(subject_ids, str) and subject_ids.strip() == "":
        subject_id_regex = r"*"
    elif isinstance(subject_ids, str):
        subject_id_regex = subject_ids.strip()
    elif isinstance(subject_ids, int):
        subject_id_regex = str(subject_ids)
    else:
        subject_id_regex = "|".join([str(s) for s in subject_ids])
    return subject_ids, subject_id_regex


def list_subject_files(dir: str,
                       subject_regex: str | None = None,
                       suffix: str | None = None) -> list:
    """Filters a directory and returns a list of files matching a subject ID and optional file suffix."""
    # Ensure that the base directory exists
    dir = to_absolute_path(dir)
    if not os.path.exists(dir):
        raise FileNotFoundError(f"Directory not found: {dir}")
    # Compile a regular expression matching subject(s) and relevant suffix
    suffix = r"*" if suffix is None else suffix
    path_regex = re.compile(f"{subject_regex}{suffix}$")
    # Filter files in directory by regex
    subject_files = [
        os.path.join(dir, filepath) for filepath in glob.glob("*", root_dir=dir)
        if path_regex.search(filepath)
    ]
    if len(subject_files) == 0:
        logger.warning(f"No matching subject files found in {dir}")
    return sorted(subject_files)


def subject_files_dict(dir: str,
                       subject_regex: str | None = None,
                       suffix: str | None = None) -> defaultdict:
    """Returns a dictionary of subject IDs and files per subject ID matching a given pattern within a directory."""
    subject_files = list_subject_files(
        dir=dir, subject_regex=subject_regex, suffix=suffix
    )
    suffix = r"*" if suffix is None else suffix
    files_per_subject = defaultdict(lambda: [])
    for subject_file in subject_files:
        subject_id = re.search(rf"{subject_regex}(?={suffix})", os.path.basename(subject_file))
        if not subject_id:
            logger.error(f"Could not extract subject ID from {subject_file}")
            continue
        subject_id = subject_id.group()
        files_per_subject[subject_id].append(subject_file)
    return files_per_subject
 

def load_object_positions_data(filepath: str, sep: str = ",") -> pd.DataFrame:
    """Load and preprocess CSV file with object positions data."""
    # Load from CSV file
    obj_pos_data = pd.read_csv(filepath, sep=sep)
    # Rename OBJECT_FIELD ("object") column to WORD_FIELD ("text") and change to title casing
    obj_pos_data.rename(columns={OBJECT_FIELD: WORD_FIELD}, inplace=True)
    obj_pos_data[WORD_FIELD] = obj_pos_data[WORD_FIELD].apply(lambda x: x.title())
    # Drop condition column
    obj_pos_data = obj_pos_data.drop(["condition"], axis=1)
    return obj_pos_data
