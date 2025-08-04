import logging
import os

import pandas as pd

from dgame.A_export_audio_and_et_times import main as step_a
from dgame.B_prepare_words import main as step_b
from dgame.Ca_preproc_et_data import main as step_ca
from dgame.Cb_preproc_fixations import main as step_cb
from dgame.Cc_prepare_fixations_for_matlab import main as step_cc
from dgame.constants import (OBJECT_FIELD, STEP_A_KEY, STEP_B_KEY, STEP_CA_KEY,
                             STEP_CB_KEY, STEP_CC_KEY, STEP_DA_KEY,
                             STEP_DB_KEY, STEP_IA_KEY, WORD_FIELD)
from dgame.Da_gaze_stats import main as step_da
from dgame.Db_plot_descriptive_fixation import main as step_db
from dgame.Ia_plot_rerps import main as step_ia
from experiment.constants import PARAM_ENABLED_KEY
from experiment.load_experiment import Experiment

logger = logging.getLogger(__name__)

DGAME_KEY = "dgame"


class DGAME(Experiment):
    def __init__(self, config_path):
        # Initialize Experiment from config
        super().__init__(config_path)

        # Set experiment data paths
        self.set_data_directories()

        # Load object and filler words of interest
        self.objects = self.load_target_words("objects")
        self.fillers = self.load_target_words("fillers")

    def set_data_directories(self) -> None:
        """Set paths to data input and output directories."""
        self.input_dir = os.path.abspath(self.config["data"]["input"]["root"])
        self.preproc_dir = os.path.join(self.input_dir, self.config["data"]["input"]["preproc_dir"])
        # Audio
        self.audio_dir = self.config["data"]["input"]["audio_dir"]
        self.audio_indir = os.path.join(self.preproc_dir, self.audio_dir)
        self.audio_outdir = os.path.join(self.outdir, self.audio_dir)
        # EEG
        self.eeg_dir = self.config["data"]["input"]["eeg_dir"]
        self.eeg_indir = os.path.join(self.preproc_dir, self.eeg_dir)
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
        # Recordings
        self.recordings_dir = self.config["data"]["input"]["recordings_dir"]
        self.recordings_indir = os.path.join(self.input_dir, self.recordings_dir)
        # xdf
        self.xdf_dir = self.config["data"]["input"]["xdf_dir"]
        self.xdf_indir = os.path.join(self.recordings_indir, self.xdf_dir)
        # MATLAB root, where dependencies/toolboxes are mounted
        self.matlab_root = os.path.abspath(self.config["analysis"]["matlab_root"])

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
    
    def get_dgame_step_parameter(self, *parameter_keys: str, default=None):
        """Get a DGAME stage parameter from the experiment config."""
        return self.get_parameter(DGAME_KEY, *parameter_keys, default=default)

    def run_analysis(self) -> None:
        """Run all component DGAME analysis steps."""
        # Run Step A: export audio and ET times [via MATLAB]
        if self.get_dgame_step_parameter(STEP_A_KEY, PARAM_ENABLED_KEY):
            self = step_a(self)
        # Run Step B: prepare words data
        if self.get_dgame_step_parameter(STEP_B_KEY, PARAM_ENABLED_KEY):
            self = step_b(self)
        # Run Step Ca: preproc ET data
        if self.get_dgame_step_parameter(STEP_CA_KEY, PARAM_ENABLED_KEY):
            self = step_ca(self)
        # Run Step Cb: preproc fixations
        if self.get_dgame_step_parameter(STEP_CB_KEY, PARAM_ENABLED_KEY):
            self = step_cb(self)
        # Run Step Cc: prepare fixations for MATLAB
        if self.get_dgame_step_parameter(STEP_CC_KEY, PARAM_ENABLED_KEY):
            self = step_cc(self)
        # Run Step Da: calculate gaze statistics
        if self.get_dgame_step_parameter(STEP_DA_KEY, PARAM_ENABLED_KEY):
            self = step_da(self)
        # Run Step Db: plot descriptive fixation
        if self.get_dgame_step_parameter(STEP_DB_KEY, PARAM_ENABLED_KEY):
            self = step_db(self)
        # TODO Step E
        # TODO Step F
        # TODO Step G
        # TODO Step H
        if self.get_dgame_step_parameter(STEP_IA_KEY, PARAM_ENABLED_KEY):
            self = step_ia(self)
