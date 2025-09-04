import os


class InputValidationError(Exception):
    pass


class OutputValidationError(Exception):
    pass


def assert_input_file_exists(filepath) -> None:
    try:
        assert os.path.exists(filepath)
    except AssertionError as exc:
        raise InputValidationError(f"Expected input file {filepath} was not found")


def assert_output_file_exists(filepath) -> None:
    try:
        assert os.path.exists(filepath)
    except AssertionError as exc:
        raise OutputValidationError(f"Expected output file {filepath} was not found")
