import os

import pandas as pd


class InputValidationError(Exception):
    pass


class OutputValidationError(Exception):
    pass


def assert_input_file_exists(filepath) -> None:
    try:
        assert os.path.exists(filepath)
    except AssertionError as exc:
        raise InputValidationError(f"Expected input file {filepath} was not found") from exc


def assert_output_file_exists(filepath) -> None:
    try:
        assert os.path.exists(filepath)
    except AssertionError as exc:
        raise OutputValidationError(f"Expected output file {filepath} was not found") from exc


def ensure_columns_exist(df: pd.DataFrame,
                         cols: list[str],
                         label: str,
                         ) -> None:
    """Raise an InputValidationError if any column from a specified set of columns
    is missing from a given Pandas DataFrame."""
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise InputValidationError(f"Missing required column(s) in {label}: {', '.join(missing)}")
