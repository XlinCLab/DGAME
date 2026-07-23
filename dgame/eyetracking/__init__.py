import os

# PATH TO THIS PACKAGE
ET_DIR = os.path.dirname(os.path.abspath(__file__))

# R PLOTTING SCRIPTS
GAZE_PROPORTIONS_PLOT_SCRIPT = os.path.join(ET_DIR, "plot", "plot_gaze_proportions.R")
HISTOGRAM_PLOT_SCRIPT = os.path.join(ET_DIR, "plot", "plot_histogram.R")

# INPUT DATA FIELDS
GAZE_TIMESTAMP_FIELD = "gaze_timestamp"
FIXATION_ID_FIELD = "fixation_id"
SURFACE_COLUMNS = [
    "surface",
    "surface_end",
    "surface_competitor",
    "targetA_surface",
    "targetB_surface",
    "compA_surface",
    "compB_surface",
    "fillerA_surface",
    "fillerB_surface",
    "target_location",
]

# Mapping of AOI (area of interest) columns with corresponding lookup column
AOI_COLUMNS = {
    "aoi_target": "surface",
    "aoi_otherTarget": "surface_competitor",
    "aoi_comp": "targetA_surface",
    "aoi_otherComp": "targetB_surface",
    "aoi_fillerA": "fillerA_surface",
    "aoi_fillerB": "fillerB_surface",
    "aoi_empty": None,
    "aoi_other": None,
    "aoi_goal": None,
}

# DATA LABELS AND TYPES
FIXATION_LABEL = "fixation"
ERROR_LABEL = "Fehler"
SURFACE_LIST = [
    f"{first_digit}{second_digit}"
    for first_digit in range(1, 5)
    for second_digit in range(1, 5)
]
COLUMN_DATA_TYPES = {column: "boolean" for column in SURFACE_LIST + list(AOI_COLUMNS.keys())}
# Ensure that subject ID and condition always read in as strings
COLUMN_DATA_TYPES["subj"] = "string"
COLUMN_DATA_TYPES["condition"] = "string"

# NUMERICAL CONSTANTS
DEFAULT_CONFIDENCE = 0.6
