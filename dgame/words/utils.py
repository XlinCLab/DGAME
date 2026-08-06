import pandas as pd
import spacy
from spacy.language import Language
from spacy.tokens import Doc, Token

from dgame.words import (DEFAULT_SPACY_MODEL, DET_POS_LABEL, NOUN_POS_LABEL,
                         PART_OF_SPEECH_FIELD)


def load_spacy_pipeline(model_name: str = DEFAULT_SPACY_MODEL) -> Language:
    """Load and initialize an NLP pipeline from a spaCy model."""
    return spacy.load(model_name)


def tag_pretokenized_words(words: list[str], nlp_pipeline: Language) -> Doc:
    """Run a spaCy pipeline over an already-tokenized sequence of words
    (e.g. ASR output), one word per input list entry.
    Builds the Doc directly from `words` and runs it through the pipeline's
    non-tokenizer components, rather than handing raw text to the pipeline's
    own tokenizer, so that the resulting Doc's tokens stay in strict 1:1
    index correspondence with the input list."""
    doc = Doc(nlp_pipeline.vocab, words=[str(word) for word in words])
    for _, component in nlp_pipeline.pipeline:
        doc = component(doc)
    return doc


def word_frequency_rank(token: Token) -> int | None:
    """Return a token's corpus frequency rank (lower = more frequent), or None if
    the word is out-of-vocabulary or the loaded spaCy model has no frequency/vector
    data (e.g. a model without word vectors)."""
    if token.is_oov or not token.has_vector:
        return None
    return token.rank


def assign_trial_numbers(word_data: pd.DataFrame,
                         trial_counter_nouns: int,
                         trial_counter_determiners: int,
                         ) -> tuple[pd.DataFrame, int, int]:
    """Assign sequential trial numbers to noun (N) and determiner (D) rows.

    Trial numbers are unique per subject, not reset per block, so the counters
    must be threaded through across a subject's files by the caller.
    """
    word_data = word_data.copy()
    word_data["trial"] = 0
    for idx, row in word_data.iterrows():
        if row[PART_OF_SPEECH_FIELD] == NOUN_POS_LABEL:
            word_data.loc[idx, "trial"] = trial_counter_nouns
            trial_counter_nouns += 1
        elif row[PART_OF_SPEECH_FIELD] == DET_POS_LABEL:
            word_data.loc[idx, "trial"] = trial_counter_determiners
            trial_counter_determiners += 1
    return word_data, trial_counter_nouns, trial_counter_determiners
