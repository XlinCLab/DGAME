import datetime
import os
from typing import Iterable


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


def load_file_lines(filepath: str, **kwargs) -> list:
    """Loads a file into a list of lines."""
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")
    with open(filepath, "r", **kwargs) as f:
        lines = [line.strip() for line in f.readlines()]
    return lines
