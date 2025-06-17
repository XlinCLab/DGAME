import datetime
import os
import subprocess
from typing import Iterable


def get_git_commit_hash():
    """Get latest git commit hash of current repository."""
    try:
        commit_hash = subprocess.check_output(['git', 'rev-parse', "--short", 'HEAD']).strip().decode('utf-8')
        return commit_hash
    except subprocess.CalledProcessError as e:
        return ""


def create_timestamp():
    """Create time stamp with current date and time."""
    # Get the current date and time
    current_datetime = datetime.datetime.now()
    # Format the date and time as a string
    formatted_datetime = current_datetime.strftime("%Y-%m-%d_%H-%M-%S")
    return current_datetime, formatted_datetime


def setdiff(a: Iterable, b: Iterable) -> set:
    """Returns the difference between two sets."""
    if not isinstance(a, set):
        a = set(a)
    if not isinstance(b, set):
        b = set(b)
    return a.difference(b)


def convert_sets_to_lists(obj):
    """Recursively convert list objects to sets."""
    if isinstance(obj, dict):
        return {k: convert_sets_to_lists(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_sets_to_lists(i) for i in obj]
    elif isinstance(obj, set):
        return list(obj)
    else:
        return obj


def load_file_lines(filepath: str, **kwargs) -> list:
    """Loads a file into a list of lines."""
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")
    with open(filepath, "r", **kwargs) as f:
        lines = [line.strip() for line in f.readlines()]
    return lines
