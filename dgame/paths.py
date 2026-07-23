import os

# PATH TO DGAME SCRIPTS
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# FILE PATH PATTERNS
# Times files
TIMES_FILE_SUFFIX = r"_times_.*\.txt"
TIMESTAMPS_FILE_SUFFIX = r"_timestamps_.*\.txt"

# Words files
AUDIO_FILE_SUFFIX = r"_words_(\d+)\.csv"
WORDS_ANNOTATED_FILE_SUFFIX = r"_words_(\d+)_annotated\.csv"

# Eyetracking (fixation / gaze / object positions) files
FIXATIONS_FILE_SUFFIX = r"fixations_on_surface_(\d+)\.csv"
FIXATION_TIMES_TRIALS_SUFFIX = r"fixations_times_(\d+)_trials\.csv"
GAZE_POS_SURFACE_SUFFIX = r"gaze_positions_on_surface_\d+\.csv"
GAZE_POSITIONS_FILE = "gaze_positions.csv"
OBJECT_POSITIONS_FILE = "object_positions.csv"

# EEG / ERP files
ERP_NOUN_FILE_SUFFIX = r"_[\w\d]+_unfold_N\.csv"
ERP_FIXATION_FILE_SUFFIX = r"_[\w\d]+_unfold_FIX\.csv"
