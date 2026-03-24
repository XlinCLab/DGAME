import re
import subprocess
from pathlib import Path

import rpy2.robjects as robjects
from rpy2.robjects import StrVector

from utils.utils import load_file_lines

MINIMUM_R_VERSION = "4.4.0"
# Path to r_requirements.txt in DGAME root directory
R_REQUIREMENTS_FILE = (Path(__file__).parent.parent / "r_requirements.txt").resolve()
# Get list of R package dependencies
R_DEPENDENCIES = load_file_lines(R_REQUIREMENTS_FILE)

# Source R functions from dependencies.R and load and/or install dependencies
R_DEPENDENCY_FUNCTIONS = (Path(__file__).parent / "dependencies.R").resolve().as_posix()
robjects.r["source"](R_DEPENDENCY_FUNCTIONS)


class RInstallationError(Exception):
    pass


class RDependencyError(Exception):
    pass


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


def r_install_package(package: str) -> None:
    """Installs a single R package."""
    robjects.globalenv["install_if_needed"](package)


def r_install_packages(package_list: list) -> None:
    """Installs multiple R packages."""
    robjects.globalenv["install_packages_if_needed"](StrVector(package_list))
