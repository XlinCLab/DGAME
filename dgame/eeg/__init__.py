import os

# PATH TO THIS PACKAGE
EEG_DIR = os.path.dirname(os.path.abspath(__file__))

# FILE NAMES / SUFFIXES
CHANNEL_COORDS_FILE = os.path.join(EEG_DIR, "channel_positions.csv")   # constant across all experiments
HEAD_MONTAGE_FILE = os.path.join(EEG_DIR, "montage", "standard-10-5-cap385.elp")   # constant across all experiments

# INPUT DATA FIELDS
CHANNEL_FIELD = "channel"
SAGGITAL_INPUT_FIELD = "sag"
SAGGITALITY_FIELD = "saggitality"
LATERAL_INPUT_FIELD = "lat"
LATERALITY_FIELD = "laterality"
