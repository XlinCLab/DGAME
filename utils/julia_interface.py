import shutil
from pathlib import Path

from utils.utils import load_file_lines

JULIA_DEPENDENCIES_FILE = (Path(__file__).parent.parent / "julia_requirements.txt").resolve().as_posix()
JULIA_DEPENDENCIES = load_file_lines(JULIA_DEPENDENCIES_FILE)


class JuliaInstallationError(Exception):
    pass


class JuliaDependencyError(Exception):
    pass


def ensure_julia_installed() -> str | JuliaInstallationError:
    # Check if Julia is installed and raise an error if not
    julia_bin = shutil.which("julia")
    if julia_bin is None:
        raise JuliaInstallationError("Julia is not installed")
    return julia_bin


def setup_julia_environment(julia_dependencies: list = None,
                            julia_dir: str = "."
                            ) -> tuple | JuliaDependencyError:
    """Activate Julia environment and install dependencies."""
    # NB: these imports need to be within this function rather than top-level
    from juliacall import Main as jl
    from juliacall import Pkg as jlPkg

    # Ensure Julia is installed before proceeding
    ensure_julia_installed()

    # Activate Julia project directory
    jlPkg.activate(julia_dir)

    # Install Julia dependencies
    try:
        jlPkg.add(jl.Vector[jl.String](julia_dependencies))
    except Exception as exc:
        raise JuliaDependencyError from exc

    return jl
