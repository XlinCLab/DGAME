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
from dgame.constants import (BLOCK_IDS, CHANNEL_COORDS_FILE, CHANNEL_FIELD,
                             DGAME_DEFAULT_CONFIG, OBJECT_FIELD,
                             OBJECT_POSITIONS_FILE, SCRIPT_DIR, STEP_A_KEY,
                             STEP_B_KEY, STEP_CA_KEY, STEP_CB_KEY, STEP_CC_KEY,
                             STEP_DA_KEY, STEP_DB_KEY, STEP_F_KEY, STEP_G_KEY,
                             STEP_H_KEY, STEP_IA_KEY, STEP_J_KEY, WORD_FIELD)
from dgame.Da_gaze_stats import main as step_da
from dgame.Db_plot_descriptive_fixation import main as step_db
from dgame.F_preproc_EEG import main as step_f
from dgame.G_deconvolution_ERPs import main as step_g
from dgame.H_reconstruct_ERPs import main as step_h
from dgame.Ia_plot_rerps import main as step_ia
from dgame.J_lm_permute_and_plot_fixations_and_language import main as step_j
from dgame.matlab_scripts.dependencies import (MATLAB_DEPENDENCIES,
                                               MATLAB_VERSION)
from dgame.plot.r_dependencies import (MINIMUM_R_VERSION, R_DEPENDENCIES,
                                       RDependencyError, RInstallationError,
                                       get_r_version)
from experiment.constants import PARAM_ENABLED_KEY
from experiment.input_validation import (InputValidationError,
                                         assert_input_file_exists)
from experiment.load_experiment import Experiment
from utils.matlab_interface import (MATLABDependencyError,
                                    MATLABInstallationError,
                                    find_matlab_installation,
                                    run_matlab_script, validate_matlab_version)
from utils.r_utils import r_install_packages

logger = logging.getLogger(__name__)


class DGAME(Experiment):
    def __init__(self,
                 config: str | dict,
                 matlab_version: str = MATLAB_VERSION,
                 minimum_r_version: str = MINIMUM_R_VERSION,
                 ):
        # Initialize Experiment from config
        super().__init__(config, default_config=DGAME_DEFAULT_CONFIG)

        # Set experiment data paths, validate input directory, and create output directories
        self.set_data_directories()
        self.validate_inputs()
        self.create_experiment_outdirs()

        # Configure compatible MATLAB version
        # Default version is MATLAB R2021a
        self.matlab_version = self.configure_matlab(matlab_version)

        # Configure R version
        self.r_version = self.configure_r(minimum_r_version)

        # Load EEG channel coordinates
        self.channel_coords = self.load_channel_coords()

        # Load object and filler words of interest
        self.objects = self.load_target_words("objects")
        self.fillers = self.load_target_words("fillers")

        # Initialize DGAME analysis steps
        self.analysis_steps = {
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
            STEP_J_KEY: step_j,
        }

    def set_data_directories(self) -> None:
        """Set paths to data input and output directories."""
        self.input_dir = os.path.abspath(self.config["data"]["input"]["root"])
        self.preproc_dir = os.path.join(self.input_dir, self.config["data"]["input"]["preproc_dir"])
        # Recording inputs
        self.recordings_dir = self.config["data"]["input"]["recordings_dir"]
        self.recordings_indir = os.path.join(self.input_dir, self.recordings_dir)
        # Audio input and output
        self.audio_dir = self.config["data"]["input"]["audio_dir"]
        self.audio_indir = os.path.join(self.recordings_indir, self.audio_dir)
        self.preproc_audio_indir = os.path.join(self.preproc_dir, self.audio_dir)
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

    def validate_inputs(self) -> None:
        """Validate that all required input directories and files exist."""

        # Get recordings/xdf directory paths per subject
        subject_xdf_dirs_dict = self.get_subject_dirs_dict(self.xdf_indir)
        subject_xdf_dir_list = []
        xdf_subject_ids = []
        for subject_id, subject_xdf_dirs in subject_xdf_dirs_dict.items():
            # Verify that there is only one xdf directory per subject
            try:
                assert len(subject_xdf_dirs) == 1
            except AssertionError as exc:
                raise InputValidationError(f">1 recordings/xdf directory found for subject <{subject_id}>") from exc
            subject_xdf_dir = subject_xdf_dirs[0]
            subject_xdf_dir_list.append(subject_xdf_dir)
            xdf_subject_ids.append(subject_id)

            # Verify that the xdf directory contains all required files
            subject_xdf_director_dir = os.path.join(subject_xdf_dir, "Director")
            for block in BLOCK_IDS:
                xdf_file = os.path.join(subject_xdf_director_dir, f"dgame2_{subject_id}_Director_{str(block)}.xdf")
                assert_input_file_exists(xdf_file)

        # Assert the found list of subject IDs matches the existing subject_ids attribute
        xdf_subject_ids.sort()
        if len(self.subject_ids) > 0 and sorted(self.subject_ids) != xdf_subject_ids:
            missing = ", ".join([subj_id for subj_id in self.subject_ids if subj_id not in xdf_subject_ids])
            raise InputValidationError(f"Subject ID(s) <{missing}> are missing from: {self.xdf_indir}")
        # Reassign subject IDs as xdf_subject_ids, if the result from parse_subject_ids was an empty list
        # (this would happen if no subject_ids were specified, in order to use data from all available subjects)
        elif len(self.subject_ids) == 0:
            self.subject_ids = xdf_subject_ids
            logger.info(f"Auto-identified {len(xdf_subject_ids)} subject(s) from xdf input files: {', '.join(xdf_subject_ids)}")

        # Ensure preproc/audio directory contains same subjects as recordings/xdf
        subj_preproc_audio_dirs_dict = self.get_subject_dirs_dict(self.preproc_audio_indir)
        audio_subj_ids = sorted(list(subj_preproc_audio_dirs_dict.keys()))
        if audio_subj_ids != xdf_subject_ids:
            missing_audio = [subject_id for subject_id in xdf_subject_ids if subject_id not in subj_preproc_audio_dirs_dict]
            missing_xdf = [subject_id for subject_id in audio_subj_ids if subject_id not in xdf_subject_ids]
            if len(missing_audio) > 0:
                raise InputValidationError(f"preproc/audio directory missing for following subjects: {', '.join(missing_audio)}")
            if len(missing_xdf) > 0:
                raise InputValidationError(f"recordings/xdf directory missing for following subjects: {', '.join(missing_audio)}")

        # Ensure other directories contain all expected files per subject
        for subject_id, subj_preproc_audio_dirs in subj_preproc_audio_dirs_dict.items():
            # Verify that there is only one preproc/audio directory per subject
            try:
                assert len(subject_xdf_dirs) == 1
            except AssertionError as exc:
                raise InputValidationError(f">1 preproc/audio directory found for subject <{subject_id}>") from exc
            subj_preproc_audio_dir = subj_preproc_audio_dirs[0]

            subj_times_dir = os.path.join(self.times_indir, subject_id)
            for block in BLOCK_IDS:
                # preproc/audio directory files per subject per block
                words_file = os.path.join(subj_preproc_audio_dir, f"{subject_id}_words_{block}.csv")
                assert_input_file_exists(words_file)
                words2erp_file = os.path.join(subj_preproc_audio_dir, f"{subject_id}_words2erp_{block}.csv")
                assert_input_file_exists(words2erp_file)

                # preproc/helper_files directory files per subject per block
                timestamp_file = os.path.join(subj_times_dir, f"{subject_id}_timestamps_max-min_{block}.csv")
                assert_input_file_exists(timestamp_file)
    
            # preproc/object_positions directory
            obj_positions_file = os.path.join(self.object_pos_indir, subject_id, OBJECT_POSITIONS_FILE)
            assert_input_file_exists(obj_positions_file)

    def create_experiment_outdirs(self):
        """Create all required per-subject output directories."""
        # Create directories per subject
        for subject_id in self.subject_ids:
            for base_dir in [
                # recordings/audio 
                self.audio_indir,
                # preproc/eeg
                self.eeg_outdir,
                # preproc/eyetracking/fixations
                self.fixations_outdir,
                # preproc/eyetracking/gaze_positions
                self.gaze_outdir,
            ]:
                subject_dir = os.path.join(base_dir, subject_id)
                os.makedirs(subject_dir, exist_ok=True)
            # preproc/eeg/{subject_id}/unfold_out
            unfold_out_dir = os.path.join(self.eeg_outdir, subject_id, "unfold_out")
            os.makedirs(unfold_out_dir, exist_ok=True)

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
        case_insensitive = self.get_dgame_step_parameter(STEP_B_KEY, "case_insensitive", default=True)
        targets = self.get_experiment_parameter(label)
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
        """Load EEG channel coordinates file."""
        channel_coords_file = os.path.join(CHANNEL_COORDS_FILE)
        channel_coords = pd.read_csv(channel_coords_file, names=[CHANNEL_FIELD, "lat", "sag", "z"], sep=sep)
        channel_coords[CHANNEL_FIELD] = channel_coords[CHANNEL_FIELD].astype(str)
        return channel_coords
    
    def get_dgame_step_parameter(self, *parameter_keys: str, default=None):
        """Get a DGAME stage parameter from the experiment config."""
        return self.get_analysis_parameter("steps", *parameter_keys, default=default)

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
        for step_id, step_func in self.analysis_steps.items():
            self.run_analysis_step(step_id, step_func)


def validate_dgame_input(x) -> DGAME:
    """Validate that an input is a DGAME experiment instance."""
    if not isinstance(x, DGAME):
        raise TypeError(
            f"Expected experiment input to be a DGAME instance, instead found {type(x)}"
        )
    return x
