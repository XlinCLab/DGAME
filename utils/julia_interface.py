import shutil
from pathlib import Path

from juliacall import Pkg as jlPkg

from utils.utils import load_file_lines

JULIA_DEPENDENCIES_FILE = (Path(__file__).parent.parent / "julia_requirements.txt").resolve().as_posix()
JULIA_DEPENDENCIES = load_file_lines(JULIA_DEPENDENCIES_FILE)


class JuliaInstallationError(Exception):
    pass


class JuliaDependencyError(Exception):
    pass


def ensure_julia_installed() -> None | JuliaInstallationError:
    # Check if Julia is installed and raise an error if not
    if shutil.which("julia") is None:
        raise JuliaInstallationError("Julia is not installed")


def setup_julia_environment(julia_dependencies: list = None,
                            julia_dir: str = "."
                            ) -> None | JuliaDependencyError:
    """Activate Julia environment and install dependencies."""
    # Make sure Julia is installed
    ensure_julia_installed()

    # Activate the environment
    jlPkg.activate(julia_dir)

    # Install dependencies in activated environment
    for pkg in julia_dependencies:
        try:
            jlPkg.add(pkg)
        except Exception as exc:
            raise JuliaDependencyError from exc
