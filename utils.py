import datetime
import os
import subprocess
from typing import Callable, Iterable

import pandas as pd


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


def merge_dataframes_with_temp_transform(left_df: pd.DataFrame,
                                        right_df: pd.DataFrame,
                                        on: str,
                                        how: str,
                                        transform: Callable,
                                        transform_left: bool = True,
                                        transform_right: bool = True,
                                        temp_column_name: str = None,
                                        **kwargs):
    """Perform a temporary transformation on a dataframe column,merge on the column with transformed values,
    then remove temporary column from the merged dataframe."""

    # Designate temporary column name and make sure it does not already exist in either dataframe
    temp_column_name = f'temp_{on}' if temp_column_name is None else temp_column_name
    try:
        assert temp_column_name not in left_df.columns or transform_left is False
        assert temp_column_name not in right_df.columns or transform_right is False
    except AssertionError:
        raise ValueError(f"Column '{temp_column_name}' already exists in dataframe")

    # Perform transform on one or both dataframes
    if not transform_left and not transform_right:
        raise ValueError("At least one of transform_left or transform_right must be True")
    if transform_left:
        left_df[temp_column_name] = left_df[on].apply(transform)
    elif temp_column_name not in left_df.columns:
        raise ValueError(f"transform_left = False and '{temp_column_name}' column is missing from left dataframe")
    if transform_right:
        right_df[temp_column_name] = right_df[on].apply(transform)
    elif temp_column_name not in right_df.columns:
        raise ValueError(f"transform_left = False and '{temp_column_name}' column is missing from right dataframe")

    # Merge on the transformed column
    merged_df = pd.merge(
        left=left_df,
        right=right_df,
        on=temp_column_name,
        how=how,
        **kwargs
    )

    # Drop the transformed column
    columns_to_drop = [temp_column_name]
    # Try to also drop the _x and/or _y columns corresponding to the untransformed target column, which appears as duplicate
    duplicate_x = f'{on}_x'
    duplicate_y = f'{on}_y'
    if duplicate_x in merged_df.columns:
        # If merge strategy is how = 'left', rename _x column instead to original target column name
        if how == "left":
            merged_df = merged_df.rename(columns={duplicate_x: on})
        else:
            columns_to_drop.append(duplicate_x)
    if duplicate_y in merged_df.columns:
        # If merge strategy is how = 'right', rename _y column instead to original target column name
        if how == "right":
            merged_df = merged_df.rename(columns={duplicate_y: on})
        else:
            columns_to_drop.append(duplicate_y)
    merged_df = merged_df.drop(columns=columns_to_drop)

    return merged_df


def idx_should_be_skipped(idx: int, skip_indices: list | dict) -> bool:
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


def get_continuous_indices(center_idx: int, candidate_indices: list, direction: str) -> list:
    """
    Get the longest continuous sequence of indices adjacent to center_idx.
    
    direction: "pre" for indices before center_idx, "post" for after.
    """
    if direction not in {"pre", "post"}:
        raise ValueError("direction must be 'pre' or 'post'")
    
    sorted_indices = sorted(candidate_indices, reverse=(direction == "pre"))
    continuous = []
    expected_idx = center_idx - 1 if direction == "pre" else center_idx + 1

    for idx in sorted_indices:
        if idx == expected_idx:
            continuous.append(idx)
            expected_idx += -1 if direction == "pre" else 1
        elif (direction == "pre" and idx < expected_idx) or (direction == "post" and idx > expected_idx):
            # Break once discontinuous index is reacheed
            break
    return continuous[::-1] if direction == "pre" else continuous
