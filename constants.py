import logging
import re

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(name)s %(levelname)s: %(message)s')

# FILE SUFFIXES
AUDIO_FILE_SUFFIX = r"_words_(\d+)\.csv"
AUDIO_ERP_FILE_SUFFIX = r"_words2erp_(\d+)\.csv"
GAZE_POS_SURFACE_SUFFIX = r"gaze_positions_on_surface_\d+.csv"
TIMES_FILE_SUFFIX = r"_times_.*.csv"
TIMESTAMPS_FILE_SUFFIX = r"_timestamps_.*.csv"

# CONSTANT WORD CLASSES
DEFINITE_ARTICLES = {"die", "der"}

# CORPORA
CORPORA = {
    "deu_news_2012_3M": "http://api.wortschatz-leipzig.de/ws/words/deu_news_2012_3M/word/"
}
DEFAULT_CORPUS = "deu_news_2012_3M"

# INPUT DATA FIELDS (and, if relevant, what they should be renamed to)
# "line" -> "id"
INPUT_LINE_ID_FIELD = "line"
WORD_ID_FIELD = "id"
# "tmin" -> "time"
INPUT_WORD_ONSET_FIELD = "tmin"
WORD_ONSET_FIELD = "time"
WORD_END_FIELD = "tmax"
# "object" -> "text"
WORD_FIELD = "text"
OBJECT_FIELD = "object"
FREQ_CLASS_FIELD = "frequencyClass"
GAZE_TIMESTAMP_FIELD = "gaze_timestamp"
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
AOI_COLUMNS = [
    "aoi_target",
    "aoi_otherTarget",
    "aoi_comp",
    "aoi_otherComp",
    "aoi_fillerA",
    "aoi_fillerB",
    "aoi_empty",
    "aoi_other",
    "aoi_goal",
]

# DATA LABELS
DET_POS_LABEL = "D"  # determiner / definite article
NOUN_POS_LABEL = "N"
PREV_WORD_LABEL = "prev"
NEXT_WORD_LABEL = "next"
CONFLICT_LABEL = "conflict"  # label for objects where one is hidden for one participant
NO_CONFLICT_LABEL = "no_conflict"  # label for objects which are visible to both participants
CONDITIONS = {CONFLICT_LABEL, NO_CONFLICT_LABEL}
SURFACE_LIST = [
    f"{first_digit}{second_digit}"
    for first_digit in range(1, 5)
    for second_digit in range(1, 5)
]
ERROR_LABEL = "Fehler"

# NUMERICAL CONSTANTS
ROUND_N = 7
TRIAL_TIME_OFFSET = 3.5
DEFAULT_CONFIDENCE = 0.6

# SYSTEM CONSTANTS
RUN_CONFIG_KEY = "run"
