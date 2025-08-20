import logging
import os
import platform
import re
import subprocess

logger = logging.getLogger(__name__)

MATLAB_VERSION_REGEX = re.compile(r'^[Rr]?20[12]\d[AaBb]$')
DEFAULT_MATLAB_VERSION = "R2021a"


def validate_matlab_version(version: str) -> str:
    """Check that MATLAB version name is valid and convert to standardized form."""
    if not MATLAB_VERSION_REGEX.match(version):
        raise ValueError(f"Invalid MATLAB version <{version}>")
    version = version.lower()
    version = re.sub(r"^[Rr]?(?=20)", "R", version)
    return version
    

def find_matlab_installation(version: str) -> str:
    """Find path to the installation of the specified version of MATLAB, depending on OS platform."""
    version = validate_matlab_version(version)
    system = platform.system()
    if system == "Darwin":  # Mac
        matlab_bin = f"/Applications/MATLAB_{version}.app/bin/matlab"
    elif system == "Linux":
        matlab_bin = f"/usr/local/MATLAB/{version}"
    elif system == "Windows":
        raise NotImplementedError("Not yet implemented for Windows OS")  # TODO
    else:
        raise ValueError(f"Unsupported OS '{system}'")
    if not os.path.exists(matlab_bin):
        raise FileNotFoundError(f"No MATLAB_{version} installation found at {matlab_bin}")
    return matlab_bin


def run_matlab_script(script_name: str,
                      args: list = None,
                      matlab_version: str = DEFAULT_MATLAB_VERSION,
                      ):
    """
    Run a MATLAB script from Python using the command line.
    
    Parameters
    ----------
    script_name : str
        Name of the MATLAB script or function (without .m extension).
    args : list, optional
        Arguments to pass to the script.
    matlab_version : str
        MATLAB installation version, e.g. '2021a' or 'R2021a'
    """
    matlab_version = validate_matlab_version(matlab_version)
    logger.info(f"Using MATLAB version: {matlab_version}")
    matlab_bin = find_matlab_installation(matlab_version)

    # Convert args to MATLAB-compatible syntax
    def to_matlab(arg):
        if isinstance(arg, list):
            return "{" + ",".join([to_matlab(a) for a in arg]) + "}"
        elif isinstance(arg, str):
            return f"'{arg}'"
        else:
            return str(arg)

    matlab_args = ",".join([to_matlab(a) for a in args]) if args else ""

    # Build MATLAB command: run script, then exit
    script_dir = os.path.dirname(script_name)
    if script_name.endswith(".m"):
        script_name = os.path.splitext(os.path.basename(script_name))[0]
    matlab_cmd = f"addpath('{script_dir}'); {script_name}({matlab_args}); exit"

    # Run MATLAB in no-GUI, no-desktop mode
    cmd = [
        matlab_bin,
        "-batch", matlab_cmd   # MATLAB R2019a+ supports -batch
    ]

    logger.info(f"Running MATLAB command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
