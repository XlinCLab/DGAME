import glob
import logging
import os
import re
from collections import defaultdict
from typing import Iterable

import pandas as pd
import yaml

from constants import OBJECT_FIELD, RUN_CONFIG_KEY, WORD_FIELD
from utils import convert_sets_to_lists, create_timestamp, get_git_commit_hash

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

    # TODO add default config from which to recursively inherit values
    # recursively_inherit_dict_values(loaded_config, included_config)

    return init_run_config(loaded_config)


def init_run_config(config: dict):
    """Initialize run config with software version, start time, and duration."""
    _, timestamp = create_timestamp()
    run_config = {
        "version": get_git_commit_hash(),
        "start_time": timestamp,
        "duration": {},
    }
    config[RUN_CONFIG_KEY] = run_config

    return config


def dump_config(config, config_outpath):
    """Dumps a YAML config to an output file."""
    with open(config_outpath, "w") as f:
        yaml.dump(convert_sets_to_lists(config), f)


def get_experiment_id(config: dict, add_timestamp: bool = False) -> str:
    """Retrieve experiment ID from config (if set) and optionally combine with timestamp."""
    # First try to retrieve the experiment_id from "run" section,
    # which is set by parent script when running multiple steps/modules of DGAME experiment
    if RUN_CONFIG_KEY in config and "id" in config[RUN_CONFIG_KEY]:
        return config[RUN_CONFIG_KEY]["id"]

    experiment_id = config["experiment"].get("id")
    _, timestamp = create_timestamp()
    if experiment_id is None or experiment_id.strip() == "":
        experiment_id = timestamp
    elif add_timestamp:
        experiment_id = os.path.join(experiment_id, timestamp)
    return experiment_id


def create_experiment_outdir(config: dict, experiment_id: str = None) -> str:
    """Retrieve and create experiment output directory."""
    # First try to retrieve the experiment_id from "run" section,
    # which is set by parent script when running multiple steps/modules of DGAME experiment
    if RUN_CONFIG_KEY in config and "outdir" in config[RUN_CONFIG_KEY]:
        return config[RUN_CONFIG_KEY]["outdir"]

    base_output_dir = config["experiment"].get("outdir")
    if base_output_dir is None or base_output_dir.strip() == "":
        base_output_dir = "out"
    experiment_id = get_experiment_id(config) if experiment_id is None else experiment_id
    output_dir = os.path.join(os.path.abspath(base_output_dir), experiment_id)
    if os.path.exists(output_dir):
        # If experiment outdir already exists, ask user to confirm
        # in order to avoid potentially overwriting previous results
        def ask_user_to_confirm_overwrite(output_dir):
            valid_answers = {"y", "yes", "n", "no"}
            overwrite_previous = None
            while overwrite_previous not in valid_answers:
                overwrite_previous = input(f"Output directory {output_dir} already exists. Overwrite previous results? [Y/N]")
                overwrite_previous = overwrite_previous.lower().strip()
            return overwrite_previous in {"y", "yes"}
        overwrite_previous = ask_user_to_confirm_overwrite(output_dir)
        if overwrite_previous is False:
            # If user specifies to NOT overwrite previous results, add timestamp
            # in order to disambiguate and keep previous results
            _, timestamp = create_timestamp()
            experiment_id = "_".join([experiment_id, timestamp])
            output_dir = os.path.join(os.path.abspath(base_output_dir), experiment_id)

    os.makedirs(output_dir, exist_ok=True)
    logger.info(f"Created experiment output directory: {output_dir}")
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


def subject_dirs_dict(root_dir: str,
                      subject_regex: str | None = None,
                      ) -> dict:
    """List subdirectories matching a subject regex."""
    subject_regex = r"*" if subject_regex is None else subject_regex
    subject_regex = re.compile(subject_regex)
    subject_subdirs = defaultdict(lambda: [])
    for path in glob.glob(os.path.join(root_dir, "*")):
        if os.path.isdir(path):
            subject_id = subject_regex.match(os.path.basename(path))
            if subject_id:
                subject_subdirs[subject_id.group()].append(path)
    return subject_subdirs


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
