from re import Pattern

import pandas as pd
import spacy
from spacy.language import Language
from spacy.tokens import Doc, Token

from dgame.words import (ADJ_UPOS_TAG, DEFAULT_SPACY_MODEL,
                         DEFINITE_MORPH_FEATURE, DEFINITE_MORPH_VALUE,
                         DET_POS_LABEL, DET_UPOS_TAG, DISFLUENCY_PATTERNS,
                         DISFLUENCY_TAG, NOUN_POS_LABEL, PART_OF_SPEECH_FIELD,
                         PUNCT_UPOS_TAG)


def load_spacy_pipeline(model_name: str = DEFAULT_SPACY_MODEL) -> Language:
    """Load and initialize an NLP pipeline from a spaCy model."""
    return spacy.load(model_name)


def tag_words_excluding_disfluencies(
        words: list[str],
        nlp_pipeline: Language,
        disfluency_patterns: set[Pattern] = DISFLUENCY_PATTERNS,
    ) -> tuple[Doc, list[int], dict[int, int], dict]:
    """Tag `words` with spaCy after excluding disfluency/filler words (e.g. "äh", "ähm") from the input.

    Returns:
    - the resulting Doc, tokenized from the non-disfluency words
    - `token_word_idx`: for each spaCy token index, the original index in `words` it belongs to
    - `word_to_doc_idx`: for each original word index (that wasn't filtered out), the spaCy
      token index of its representative (first non-punctuation, else first) token,
      whose POS/lemma/frequency actually describes the word, and the anchor used to search
      for a noun's determiner
    - `disfluencies`: filtered-out {index: word} pairs, for logging
    """
    kept_indices = [
        i for i, word in enumerate(words)
        if word.strip() != DISFLUENCY_TAG and not any(
            disfluency_pattern.match(word.strip().lower()) for disfluency_pattern in disfluency_patterns
        )
    ]
    kept_words = [str(words[i]) for i in kept_indices]
    disfluencies = {i: word for i, word in enumerate(words) if i not in kept_indices}

    # Track each kept word's character span in the joined text, to map spaCy's own tokens
    # (which may split a word, e.g. on attached punctuation) back to the word they came from
    char_ends = []
    char_pos = 0
    for word in kept_words:
        char_pos += len(word)
        char_ends.append(char_pos)
        char_pos += 1  # for the joining space
    doc = nlp_pipeline(" ".join(kept_words))

    token_word_idx = []
    word_ptr = 0
    for token in doc:
        while word_ptr < len(kept_words) - 1 and token.idx >= char_ends[word_ptr]:
            word_ptr += 1
        token_word_idx.append(kept_indices[word_ptr])

    word_to_doc_idx: dict[int, int] = {}
    for doc_idx, word_idx in enumerate(token_word_idx):
        if word_idx not in word_to_doc_idx:
            word_to_doc_idx[word_idx] = doc_idx
        elif doc[doc_idx].pos_ != PUNCT_UPOS_TAG and doc[word_to_doc_idx[word_idx]].pos_ == PUNCT_UPOS_TAG:
            word_to_doc_idx[word_idx] = doc_idx

    return doc, token_word_idx, word_to_doc_idx, disfluencies


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
    Returns None if no definite determiner is found."""
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
