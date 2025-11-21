import re
import subprocess
from pathlib import Path

MINIMUM_R_VERSION = "4.4.0"

# Path to r_requirements.txt in DGAME root directory
R_REQUIREMENTS_FILE = (Path(__file__).parent.parent.parent / "r_requirements.txt").resolve()


class RInstallationError(Exception):
    pass


class RDependencyError(Exception):
    pass


def list_r_dependencies(r_requirements_file: str) -> list:
    """Read a file containing R package dependencies and return these as a list."""
    with open(r_requirements_file, "r") as f:
        packages = f.readlines()
    return [package.strip() for package in packages]


def get_r_version() -> str:
    """
    Returns the installed R version as a string, or raises RInstallationError if R is not installed.
    """
    try:
        # Run R --version and capture output
        result = subprocess.run(
            ["R", "--version"],
            capture_output=True,
            text=True,
            check=True
        )

        # e.g.:
        # "R version 4.5.0 (2025-04-11) -- "How About a Twenty-Six"
        first_line = result.stdout.splitlines()[0]
        # Parse version
        r_version = re.sub(r'R version (\d+\.\d(\.\d)?) .*', r'\1', first_line)
        return r_version

    except FileNotFoundError as exc:
        raise RInstallationError("R is not installed") from exc
    except subprocess.CalledProcessError as exc:
        raise subprocess.CalledProcessError("Unknown error occurred while checking R version") from exc


# Get list of R package dependencies
R_DEPENDENCIES = list_r_dependencies(R_REQUIREMENTS_FILE)
