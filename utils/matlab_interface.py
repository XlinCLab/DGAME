import logging
import os
import platform
import re
import subprocess

logger = logging.getLogger(__name__)

MATLAB_VERSION_REGEX = re.compile(r'^[Rr]?20[12]\d[AaBb]$')
DEFAULT_MATLAB_VERSION = "R2021a"


class MATLABInstallationError(Exception):
    pass

class MATLABDependencyError(Exception):
    pass


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

    # Mac OS
    if system == "Darwin":
        matlab_bin = f"/Applications/MATLAB_{version}.app/bin/matlab"

    # Linux and Docker containers
    elif system == "Linux":
        matlab_bin = f"/usr/local/MATLAB/{version}"
        if not os.path.exists(matlab_bin):
            # Fallback: check `which matlab` and resolve symlink (e.g. for Docker containers)
            try:
                matlab_exe = subprocess.check_output(["which", "matlab"], text=True).strip()
            except subprocess.CalledProcessError:
                raise MATLABInstallationError(f"MATLAB {version} not found in PATH")

            if not os.path.exists(matlab_exe):
                raise MATLABInstallationError(f"MATLAB executable not found at {matlab_exe}")

            # Resolve symlink and remove /bin/matlab from end of path
            matlab_bin = os.path.dirname(os.path.dirname(os.path.realpath(matlab_exe)))

    # Windows
    elif system == "Windows":
        matlab_bin = fr"C:\Program Files\MATLAB\{version}\bin\matlab.exe"
        if not os.path.exists(matlab_bin):
            # MATLAB may also be installed in "Program Files (x86)" on 64-bit Windows; check there as fallback
            logger.warning(f"Could not find MATLAB {version} under C:\\Program Files\\MATLAB\\, searching instead in C:\\Program Files (x86)\\MATLAB")
            matlab_bin = fr"C:\Program Files (x86)\MATLAB\{version}\bin\matlab.exe"

    # Unknown OS
    else:
        raise OSError(f"Unsupported OS '{system}'")

    # Verify that the MATLAB binary path exists and return
    if not os.path.exists(matlab_bin):
        raise MATLABInstallationError(f"No MATLAB_{version} installation found at {matlab_bin}")
    return matlab_bin


def run_matlab_script(script_name: str,
                      args: list = None,
                      matlab_version: str = DEFAULT_MATLAB_VERSION,
                      logfile: str = None,
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
    logfile: str
        Optional logfile where MATLAB script output will be written.
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

    # Assemble command to run MATLAB in no-GUI, no-desktop mode
    cmd = [
        matlab_bin,
        "-nodisplay",
        "-nodesktop",
    ]

    # Optionally add MATLAB logfile and create its directory in case it does not already exist
    if logfile:
        logfile = os.path.abspath(logfile) if not os.path.isabs(logfile) else logfile
        logfile_dir = os.path.dirname(logfile)
        os.makedirs(logfile_dir, exist_ok=True)
        cmd.append(f"-logfile {logfile}")

    # Add further options depending on MATLAB version
    version_num = int(re.search(r'\d+', matlab_version).group())
    if version_num < 2025:  # -nosplash no longer supported as of version 2025
        cmd.append("-nosplash")
    if version_num >= 2019:
        cmd.append("-batch")  # -batch (non-interactive mode) supported for MATLAB R2019a and beyond

    # Add actual MATLAB command
    cmd.append(matlab_cmd)

    logger.info(f"Running MATLAB script: {script_name} ")
    logger.debug(f"Running MATLAB command: {' '.join(cmd)}")
    if logfile:
        logger.info(f"See MATLAB log file for detailed output and logs: {logfile}")
        with open(logfile, "w") as logfile:
            subprocess.run(cmd, check=True, stdout=logfile, stderr=logfile, text=True)
    else:
        subprocess.run(cmd, check=True)
