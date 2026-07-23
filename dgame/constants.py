import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(name)s %(levelname)s: %(message)s')

# EXPERIMENTAL DESIGN
BLOCK_IDS = [11, 12, 21, 22]
CONFLICT_LABEL = "conflict"  # label for objects where one is hidden for one participant
NO_CONFLICT_LABEL = "no_conflict"  # label for objects which are visible to both participants
CONDITIONS = {CONFLICT_LABEL, NO_CONFLICT_LABEL}

# PARTICIPANT CONDITION LABELS
DIRECTOR_LABEL = "director"
DECKE_LABEL = "decke"
PARTICIPANT_CONDITION_LABELS = {DIRECTOR_LABEL, DECKE_LABEL}

# NUMERICAL CONSTANTS
ROUND_N = 7
TRIAL_TIME_OFFSET = 3.5  # TODO change to 1.5 once issue has been fixed
