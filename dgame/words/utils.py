import pandas as pd

from dgame.words import DET_POS_LABEL, NOUN_POS_LABEL, PART_OF_SPEECH_FIELD


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
