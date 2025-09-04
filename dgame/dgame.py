import logging
import os
import time
from typing import Callable

import pandas as pd
from packaging.version import Version

from dgame.A_export_audio_and_et_times import main as step_a
from dgame.B_prepare_words import main as step_b
from dgame.Ca_preproc_et_data import main as step_ca
from dgame.Cb_preproc_fixations import main as step_cb
from dgame.Cc_prepare_fixations_for_matlab import main as step_cc
from dgame.constants import (CHANNEL_COORDS_FILE, CHANNEL_FIELD, OBJECT_FIELD,
                             SCRIPT_DIR, STEP_A_KEY, STEP_B_KEY, STEP_CA_KEY,
                             STEP_CB_KEY, STEP_CC_KEY, STEP_DA_KEY,
                             STEP_DB_KEY, STEP_F_KEY, STEP_G_KEY, STEP_H_KEY,
                             STEP_IA_KEY, STEP_JA_KEY, STEP_JB_KEY, WORD_FIELD)
from dgame.Da_gaze_stats import main as step_da
from dgame.Db_plot_descriptive_fixation import main as step_db
from dgame.F_preproc_EEG import main as step_f
from dgame.G_deconvolution_ERPs import main as step_g
from dgame.H_reconstruct_ERPs import main as step_h
from dgame.Ia_plot_rerps import main as step_ia
from dgame.Ja_lm_permute_and_plot_fixations import main as step_ja
from dgame.matlab_scripts.dependencies import (MATLAB_DEPENDENCIES,
                                               MATLAB_VERSION)
from dgame.plot.r_dependencies import (MINIMUM_R_VERSION, R_DEPENDENCIES,
                                       RDependencyError, RInstallationError,
                                       get_r_version)
from experiment.constants import PARAM_ENABLED_KEY
from experiment.load_experiment import Experiment
from utils.matlab_interface import (MATLABDependencyError,
                                    MATLABInstallationError,
                                    find_matlab_installation,
                                    run_matlab_script, validate_matlab_version)
from utils.r_utils import r_install_packages

logger = logging.getLogger(__name__)

DGAME_KEY = "dgame"

DGAME_ANALYSIS_STEPS = {
    STEP_A_KEY: step_a,
    STEP_B_KEY: step_b,
    STEP_CA_KEY: step_ca,
    STEP_CB_KEY: step_cb,
    STEP_CC_KEY: step_cc,
    STEP_DA_KEY: step_da,
    STEP_DB_KEY: step_db,
    STEP_F_KEY: step_f,
    STEP_G_KEY: step_g,
    STEP_H_KEY: step_h,
    STEP_IA_KEY: step_ia,
    STEP_JA_KEY: step_ja,
    # STEP_JB_KEY: step_jb,
}


class DGAME(Experiment):
    def __init__(self,
                 config_path: str,
                 matlab_version: str = MATLAB_VERSION,
                 minimum_r_version: str = MINIMUM_R_VERSION,
                 ):
        # Initialize Experiment from config
        super().__init__(config_path)

        # Set experiment data paths
        self.set_data_directories()

        # Configure compatible MATLAB version
        # Default version is MATLAB R2021a
        self.matlab_version = self.configure_matlab(matlab_version) 

        # Configure R version
        self.r_version = self.configure_r(minimum_r_version)

        # Load object and filler words of interest
        self.objects = self.load_target_words("objects")
        self.fillers = self.load_target_words("fillers")

    def set_data_directories(self) -> None:
        """Set paths to data input and output directories."""
        self.input_dir = os.path.abspath(self.config["data"]["input"]["root"])
        self.preproc_dir = os.path.join(self.input_dir, self.config["data"]["input"]["preproc_dir"])
        # Recording inputs
        self.recordings_dir = self.config["data"]["input"]["recordings_dir"]
        self.recordings_indir = os.path.join(self.input_dir, self.recordings_dir)
        # Audio
        self.audio_dir = self.config["data"]["input"]["audio_dir"]
        self.audio_indir = os.path.join(self.recordings_indir, self.audio_dir)
        self.audio_outdir = os.path.join(self.outdir, self.audio_dir)
        # EEG
        self.eeg_dir = self.config["data"]["input"]["eeg_dir"]
        self.eeg_outdir = os.path.join(self.preproc_dir, self.eeg_dir)
        # Fixations
        self.fixations_dir = self.config["data"]["input"]["fixations_dir"]
        self.fixations_outdir = os.path.join(self.outdir, self.fixations_dir)
        # Gaze
        self.gaze_dir = self.config["data"]["input"]["gaze_dir"]
        self.gaze_indir = os.path.join(self.preproc_dir, self.gaze_dir)
        self.gaze_outdir = os.path.join(self.outdir, self.gaze_dir)
        # Object positions
        self.object_pos_dir = self.config["data"]["input"]["object_positions"]
        self.object_pos_indir = os.path.join(self.preproc_dir, self.object_pos_dir)
        # Times
        self.times_dir = self.config["data"]["input"]["times_dir"]
        self.times_indir = os.path.join(self.preproc_dir, self.times_dir)
        # Surfaces
        self.surface_dir = self.config["data"]["input"]["surfaces_dir"]
        self.surface_indir = os.path.join(self.preproc_dir, self.surface_dir)
        # xdf
        self.xdf_dir = self.config["data"]["input"]["xdf_dir"]
        self.xdf_indir = os.path.join(self.recordings_indir, self.xdf_dir)

    def configure_r(self, minimum_r_version: str) -> str:
        """Validate that a compatible version of R is installed and install dependencies."""
        try:
            installed_r_version = get_r_version()
        except RInstallationError as exc:
            raise RInstallationError("R is required but is not installed") from exc
        if Version(installed_r_version) >= Version(minimum_r_version):
            logger.info(f"Running R version {installed_r_version}")
        else:
            raise RDependencyError(f"DGAME requires R version >= {minimum_r_version} but version {installed_r_version} is installed")
        # Ensure R dependencies are installed
        r_install_packages(R_DEPENDENCIES)
        return installed_r_version

    def configure_matlab(self, matlab_version: str) -> str:
        """Validate MATLAB version input, ensure that version is installed, and set up MATLAB directories."""
        # Validate the version number
        matlab_version = validate_matlab_version(matlab_version)
        # Ensure the correct version of MATLAB is installed before running any analyses
        # Result of this function is not needed here but it will raise an error if the installation is not found
        try:
            find_matlab_installation(matlab_version)
        except MATLABInstallationError as exc:
            raise MATLABInstallationError(
                f"MATLAB version {self.matlab_version} is required but not installed!"
            ) from exc
        logger.info(f"Running MATLAB version {matlab_version}")

        # MATLAB root directory, where dependencies/toolboxes are mounted
        self.matlab_root = os.path.abspath(self.config["analysis"]["matlab_root"])
        # Validate that all MATLAB dependencies can be found
        missing_matlab_dependencies = []
        for matlab_dep in MATLAB_DEPENDENCIES:
            full_matlab_dep_path = os.path.join(self.matlab_root, matlab_dep)
            if not os.path.exists(full_matlab_dep_path):
                dep_basename = os.path.basename(full_matlab_dep_path)
                logger.warning(f"Could not find MATLAB dependency <{dep_basename}> within specified MATLAB root {self.matlab_root}")
                missing_matlab_dependencies.append(dep_basename)
        if len(missing_matlab_dependencies) > 0:
            missing_str = ", ".join(missing_matlab_dependencies)
            raise MATLABDependencyError(f"One or more MATLAB dependencies are missing: {missing_str}")

        # Set path to MATLAB DGAME scripts
        self.matlab_script_dir = os.path.join(SCRIPT_DIR, "matlab_scripts")

        # Create directory for MATLAB logs
        self.matlab_logdir = os.path.join(self.logdir, "MATLAB")
        os.makedirs(self.matlab_logdir, exist_ok=True)
        return matlab_version

    def run_matlab_step(self,
                        script_path: str,
                        args: list = None,
                        ) -> None:
        """Run a MATLAB script using the configured MATLAB installation."""
        # Designate MATLAB log file
        script_basename, _ = os.path.splitext(os.path.basename(script_path))
        logfile = os.path.join(self.matlab_logdir, f"{script_basename}.log")

        # Run MATLAB script with specified arguments
        run_matlab_script(
            script_path,
            args=args,
            matlab_version=self.matlab_version,
            logfile=logfile,
        )

    def load_target_words(self, label: str) -> set:
        """Initialize target object words and filler words."""
        case_insensitive = self.get_parameter("case_insensitive", True)
        targets = self.get_parameter(label)
        if case_insensitive:
            # Standardize to title casing (NB: because German nouns are capitalized)
            targets = set(obj.title() for obj in targets)
        else:
            targets = set(targets)
        return targets

    @staticmethod
    def load_object_positions_data(filepath: str, sep: str = ",") -> pd.DataFrame:
        """Load and preprocess CSV file with object positions data."""
        # Load from CSV file
        obj_pos_data = pd.read_csv(filepath, sep=sep)

        # Rename OBJECT_FIELD ("object") column to WORD_FIELD ("text") and change to title casing
        obj_pos_data.rename(columns={OBJECT_FIELD: WORD_FIELD}, inplace=True)
        obj_pos_data[WORD_FIELD] = obj_pos_data[WORD_FIELD].apply(lambda x: x.title())

        # Drop condition column
        obj_pos_data = obj_pos_data.drop(["condition"], axis=1)

        return obj_pos_data
    
    def load_channel_coords(self, sep: str = ",") -> pd.DataFrame:
        """Load channel coordinates file."""
        # TODO is this a static file across all subjects? where do we expect it to be found? eeg_dir?
        channel_coords_file = os.path.join(self.input_dir, CHANNEL_COORDS_FILE)
        channel_coords = pd.read_csv(channel_coords_file, names=[CHANNEL_FIELD, "lat", "sag", "z"], sep=sep)
        channel_coords[CHANNEL_FIELD] = channel_coords[CHANNEL_FIELD].astype(str)
        return channel_coords
    
    def get_dgame_step_parameter(self, *parameter_keys: str, default=None):
        """Get a DGAME stage parameter from the experiment config."""
        return self.get_parameter(DGAME_KEY, *parameter_keys, default=default)
    
    def run_analysis_step(self, step_id: str, step_func: Callable) -> None:
        """Run a particular DGAME analysis step."""
        if self.get_dgame_step_parameter(step_id, PARAM_ENABLED_KEY):
            logger.info(f"Running analysis step {step_id} ...")
            start_time = time.time()
            step_func(self)
            self.log_step_duration(start_time, step_id=step_id)
        else:
            logger.info(f"Skipping analysis step {step_id}")

    def run_analysis(self) -> None:
        """Run all component DGAME analysis steps."""
        for step_id, step_func in DGAME_ANALYSIS_STEPS.items():
            self.run_analysis_step(step_id, step_func)


def validate_dgame_input(x) -> DGAME:
    """Validate that an input is a DGAME experiment instance."""
    if not isinstance(x, DGAME):
        raise TypeError(
            f"Expected experiment input to be a DGAME instance, instead found {type(x)}"
        )
    return x
