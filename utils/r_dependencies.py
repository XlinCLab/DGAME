import re
import subprocess
from pathlib import Path

from utils.utils import load_file_lines

MINIMUM_R_VERSION = "4.4.0"
# Path to r_requirements.txt in DGAME root directory
R_REQUIREMENTS_FILE = (Path(__file__).parent.parent / "r_requirements.txt").resolve()
# Get list of R package dependencies
R_DEPENDENCIES = load_file_lines(R_REQUIREMENTS_FILE)

# Path to R functions sourced by dependencies.R, loaded lazily (see _load_r_dependency_functions):
# importing rpy2 embeds an R interpreter in-process and mutates LD_LIBRARY_PATH, which can break
# subprocess-based tools (e.g. Julia) if that happens before they've had a chance to run
R_DEPENDENCY_FUNCTIONS = (Path(__file__).parent / "dependencies.R").resolve().as_posix()
_r_dependency_functions_loaded = False


class RInstallationError(Exception):
    pass


class RDependencyError(Exception):
    pass


def _load_r_dependency_functions() -> None:
    """Source dependencies.R exactly once, the first time it's actually needed."""
    global _r_dependency_functions_loaded
    if not _r_dependency_functions_loaded:
        import rpy2.robjects as robjects
        robjects.r["source"](R_DEPENDENCY_FUNCTIONS)
        _r_dependency_functions_loaded = True


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
    import rpy2.robjects as robjects
    _load_r_dependency_functions()
    robjects.globalenv["install_if_needed"](package)


def r_install_packages(package_list: list) -> None:
    """Installs multiple R packages."""
    import rpy2.robjects as robjects
    from rpy2.robjects import StrVector
    _load_r_dependency_functions()
    robjects.globalenv["install_packages_if_needed"](StrVector(package_list))
