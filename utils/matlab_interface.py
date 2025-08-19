import logging
import os
import subprocess

logger = logging.getLogger(__name__)


def run_matlab_script(script_name: str,
                      args: list = None,
                      matlab_bin: str = "/Applications/MATLAB_R2021a.app/bin/matlab",
                      ):
    """
    Run a MATLAB script from Python using the command line.
    
    Parameters
    ----------
    script_name : str
        Name of the MATLAB script or function (without .m extension).
    args : list, optional
        Arguments to pass to the script.
    matlab_bin : str
        Path to the MATLAB binary.
    """
    # TODO add handling for matlab_bin, adjust per OS
    logger.info(f"Using MATLAB version: {matlab_bin}")

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
