import os

import pytest

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

# NB: XDF test data intentionally not stored in git due to privacy restrictions
SAMPLE_XDF_FILES = [
    os.path.join(REPO_ROOT, "data", "dgame-test-data", "recordings", "xdf", "02", "Director", "dgame2_02_Director_11.xdf"),
    os.path.join(REPO_ROOT, "data", "dgame-test-data", "recordings", "xdf", "27", "Director", "dgame2_27_Director_21.xdf"),
]

def _skip_if_missing(path: str) -> None:
    if not os.path.exists(path):
        pytest.skip(f"Sample XDF file not found: {path}")
