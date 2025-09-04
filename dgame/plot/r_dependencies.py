import re
import subprocess

MINIMUM_R_VERSION = "4.4.0"

R_DEPENDENCIES = [
    "dplyr",
    "eyetrackingR",
    "ggplot2",
    "pbapply",
    "stringr",
    "viridis",
]

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
        raise subprocess.CalledProcessError("Unknown error occurred while checking R version")
