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
PART_OF_SPEECH_FIELD = "pos"

# WORD DATA LABELS
DET_POS_LABEL = "D"  # determiner / definite article
NOUN_POS_LABEL = "N"
VERB_POS_LABEL = "VERB"
DIRECTION_WORD_LABEL = "DIR"
PREV_WORD_LABEL = "prev"
NEXT_WORD_LABEL = "next"
