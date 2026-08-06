import re

# Default German spaCy model
# large (`de_core_news_lg`) model selected over newer transformer pipeline (`de_dep_news_trf`)
# as the latter does not include word vectors for frequency rank
DEFAULT_SPACY_MODEL = "de_core_news_lg"

# Universal POS tags from spaCy
DET_UPOS_TAG = "DET"
ADJ_UPOS_TAG = "ADJ"

# spaCy morphological feature/value marking a determiner as definite
DEFINITE_MORPH_FEATURE = "Definite"
DEFINITE_MORPH_VALUE = "Def"

# German disfluency/filler words (filled pauses)
DISFLUENCY_TAG = "[DISFLUENCY]"
DISFLUENCY_PATTERNS = {
    re.compile(pattern) for pattern in [
        r"^[äaeöu]+h*m+$",
        r"^[äaeöu]+h+m*$",
        r"^h+m+$",
        r"^m+h+m*$",
    ]
}

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
FREQ_RANK_FIELD = "frequencyRank"
PART_OF_SPEECH_FIELD = "pos"
SPACY_POS_FIELD = "spacy_pos"
SPACY_LEMMA_FIELD = "spacy_lemma"

# WORD DATA LABELS
DET_POS_LABEL = "D"  # determiner / definite article
NOUN_POS_LABEL = "N"
VERB_POS_LABEL = "VERB"
DIRECTION_WORD_LABEL = "DIR"
PREV_WORD_LABEL = "prev"
NEXT_WORD_LABEL = "next"
