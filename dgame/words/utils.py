import logging
from re import Pattern

import pandas as pd
import spacy
from spacy.language import Language
from spacy.tokens import Doc, Token

from dgame.words import (ADJ_UPOS_TAG, DEFAULT_SPACY_MODEL,
                         DEFINITE_MORPH_FEATURE, DEFINITE_MORPH_VALUE,
                         DET_POS_LABEL, DET_UPOS_TAG, DISFLUENCY_PATTERNS,
                         NOUN_POS_LABEL, PART_OF_SPEECH_FIELD)

logger = logging.getLogger(__name__)


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


def tag_words_excluding_disfluencies(
        words: list[str],
        nlp_pipeline: Language,
        disfluency_patterns: set[Pattern] = DISFLUENCY_PATTERNS,
    ) -> tuple[Doc, list[int], dict]:
    """Tag `words` with spaCy after excluding disfluency/filler words (e.g. "äh", "ähm") from the input.

    Returns the resulting Doc, built only from the non-disfluency words, together with a
    list mapping each of its token indices back to the corresponding index in `words`
    (`doc[i]` corresponds to `words[mapping[i]]`)."""
    kept_indices = [
        i for i, word in enumerate(words)
        if not any(
            disfluency_pattern.match(word.strip().lower()) for disfluency_pattern in disfluency_patterns
        )
    ]
    kept_words = [words[i] for i in kept_indices]
    doc = tag_pretokenized_words(kept_words, nlp_pipeline)
    disfluencies = {i: word for i, word in enumerate(words) if i not in kept_indices}
    return doc, kept_indices, disfluencies


def word_frequency_rank(token: Token) -> int | None:
    """Return a token's corpus frequency rank (lower = more frequent), or None if
    the word is out-of-vocabulary or the loaded spaCy model has no frequency/vector
    data (e.g. a model without word vectors)."""
    if token.is_oov or not token.has_vector:
        return None
    return token.rank


def is_definite_determiner(token: Token) -> bool:
    """Return whether `token` is a definite-article determiner (e.g. "die"/"der"), as
    opposed to an indefinite article (e.g. "eine") or a demonstrative. The experiment's
    design always refers to target/filler objects with a definite article, never
    indefinite, so only these should count as a target/filler noun's determiner."""
    return token.pos_ == DET_UPOS_TAG and DEFINITE_MORPH_VALUE in token.morph.get(DEFINITE_MORPH_FEATURE)


def determiner_distance(doc: Doc, noun_idx: int, max_reasonable_distance: int = 5) -> int | None:
    """Return the token distance from a noun (at `noun_idx`) back to its (definite-article) determiner.
    Primarily uses the noun's dependency children, with fallback to a local, bounded backward
    walk over POS tags when the dependency attachment is missing or implausibly far away (> `max_reasonable_distance`).
    If neither approach finds a definite determiner, logs a warning and returns None: the
    experiment's design means a target/filler noun should always have a nearby definite
    determiner, so this is likely an ASR error or disfluency worth reviewing (e.g. via
    the `skip_indices` config)."""
    noun_token = doc[noun_idx]
    det_children = [child for child in noun_token.children if is_definite_determiner(child)]
    if det_children:
        # spaCy parse can attach more than one DET on unsegmented or disfluent text; take the nearest
        nearest = min(det_children, key=lambda child: noun_token.i - child.i)
        distance = noun_token.i - nearest.i
        if distance <= max_reasonable_distance:
            return distance

    # Fallback in case definite article determiner found within reasonable distance of target noun
    for distance in range(1, max_reasonable_distance + 1):
        idx = noun_idx - distance
        if idx < 0:
            break
        token = doc[idx]
        if is_definite_determiner(token):
            return distance
        if token.pos_ != ADJ_UPOS_TAG:
            break

    # Log warning and return None if still no definite article determiner is found
    logger.warning(
        f"No definite determiner found within {max_reasonable_distance} tokens before "
        f"'{noun_token.text}' (index {noun_idx}); skipping this occurrence."
    )
    return None


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
