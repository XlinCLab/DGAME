import argparse
import os
import re
from collections import defaultdict

import pandas as pd
from spacy.language import Language

from dgame.constants import CONFLICT_LABEL, NO_CONFLICT_LABEL
from dgame.paths import AUDIO_FILE_SUFFIX, OBJECT_POSITIONS_FILE
from dgame.pipeline import WORDS_PREPROCESS_STEP
from dgame.words import (DEFAULT_SPACY_MODEL, DET_POS_LABEL, DET_UPOS_TAG,
                         FREQ_RANK_FIELD, INPUT_LINE_ID_FIELD,
                         INPUT_WORD_ONSET_FIELD, NEXT_WORD_LABEL,
                         NOUN_POS_LABEL, PART_OF_SPEECH_FIELD, PREV_WORD_LABEL,
                         SPACY_LEMMA_FIELD, SPACY_POS_FIELD, WORD_END_FIELD,
                         WORD_FIELD, WORD_ID_FIELD, WORD_ONSET_FIELD)
from dgame.words.utils import (assign_trial_numbers, load_spacy_pipeline,
                               tag_pretokenized_words, word_frequency_rank)
from experiment.load_experiment import Experiment
from utils.utils import idx_should_be_skipped, setdiff


def preprocess_words_data(audio_infile: str,
                          nlp_pipeline: Language,
                          objects: set,
                          fillers: set,
                          skip_indices: set = None,
                          pattern_id: int = 1,
                          set_id: int = 1,
                          sep: str = ",") -> pd.DataFrame:
    # Load audio data CSV into dataframe
    # Explicitly specify columns of interest in case more are present
    audio_data = pd.read_csv(
        audio_infile,
        sep=sep,
        usecols=[INPUT_LINE_ID_FIELD, INPUT_WORD_ONSET_FIELD, WORD_FIELD, WORD_END_FIELD],
    )
    # Rename columns
    audio_data.columns = [WORD_ID_FIELD, WORD_ONSET_FIELD, WORD_FIELD, WORD_END_FIELD]
    # Drop duplicate entries
    audio_data.drop_duplicates(inplace=True)

    # Drop rows without any WORD_FIELD ("text") entry before running the NLP pipeline,
    # and reset the index so it stays aligned with the (list-positional) spaCy tagging below
    audio_data = audio_data.dropna(subset=[WORD_FIELD]).reset_index(drop=True)

    # Tag the word sequence with spaCy
    words = audio_data[WORD_FIELD].astype(str).to_list()
    doc = tag_pretokenized_words(words, nlp_pipeline)
    spacy_pos_tags = [token.pos_ for token in doc]
    spacy_lemmas = [token.lemma_ for token in doc]
    audio_data[SPACY_POS_FIELD] = spacy_pos_tags
    audio_data[SPACY_LEMMA_FIELD] = spacy_lemmas
    audio_data[FREQ_RANK_FIELD] = [word_frequency_rank(token) for token in doc]

    # Set missing frequency rank entries (out-of-vocabulary words) to 1 + the maximum attested rank
    max_freq_rank = audio_data[FREQ_RANK_FIELD].max(skipna=True)
    audio_data[FREQ_RANK_FIELD] = audio_data[FREQ_RANK_FIELD].fillna(max_freq_rank + 1)

    # Initialize "pattern" and "set" columns
    audio_data["pattern"] = pattern_id
    audio_data["set"] = set_id

    # Iterate over words by line ID, skipping specified indices
    # If word's lemma matches one of target object words, check if preceded by a determiner
    conditions = [None] * len(words)
    condition_codes = [None] * len(words)
    pos = [None] * len(words)
    positions = [None] * len(words)
    counts = defaultdict(lambda: 0)
    for idx, lemma in enumerate(spacy_lemmas):
        # line_id = int(line_ids[idx])  #  TODO use line IDs from raw file or index of post-filtered words?
        # if skip_indices is not None and idx_should_be_skipped(line_id):
        if skip_indices is not None and idx_should_be_skipped(idx, skip_indices):
            continue
        # Check if word's lemma matches either target objects or fillers
        if lemma in objects.union(fillers):
            # Check for preceding determiner (spaCy POS tag)
            if idx > 0 and spacy_pos_tags[idx - 1] == DET_UPOS_TAG:
                nback = 1
            else:
                nback = 2
            pos[idx] = NOUN_POS_LABEL
            pos[idx + 1] = NEXT_WORD_LABEL
            pos[idx + 2] = NEXT_WORD_LABEL
            # Update counts
            counts[lemma] += 1
            positions[idx - nback] = counts[lemma]
            positions[idx] = counts[lemma]
        # Check if word's lemma matches target objects or fillers and assign conditions/codes accordingly
        if lemma in objects:
            condition_codes[idx - nback] = 11
            condition_codes[idx] = 12
            conditions[idx - nback: idx + 1] = [CONFLICT_LABEL] * (nback + 1)
            pos[idx - 1] = DET_POS_LABEL
            pos[idx - 2] = PREV_WORD_LABEL
        elif lemma in fillers:
            condition_codes[idx - nback] = 21
            condition_codes[idx] = 22
            conditions[idx - nback: idx + 1] = [NO_CONFLICT_LABEL] * (nback + 1)
            pos[idx - nback] = DET_POS_LABEL
            pos[idx - (nback + 1)] = PREV_WORD_LABEL
            pos[idx - (nback + 2)] = PREV_WORD_LABEL
            if nback == 2:
                pos[idx - 1] = PREV_WORD_LABEL
    audio_data["condition"] = conditions
    audio_data["condition_code"] = condition_codes
    audio_data[PART_OF_SPEECH_FIELD] = pos
    audio_data["position"] = positions

    return audio_data


def combine_words_and_obj_position_data(word_data: pd.DataFrame,
                                        object_positions: pd.DataFrame) -> pd.DataFrame:
    # Merge object position data
    combined_data = pd.merge(word_data, object_positions, how='left')
    # Ensure index is 0-based integer index (in order to use .loc in loop below)
    combined_data = combined_data.reset_index(drop=True)
    # Iterate again through words
    for idx, row in combined_data.iterrows():
        if not pd.isna(row["surface"]):
            # Check if preceding word is a determiner
            if idx > 0 and combined_data[SPACY_POS_FIELD][idx - 1] == DET_UPOS_TAG:
                nback = 1
            else:
                nback = 2
            for col in [
                "surface",
                "surface_competitor",
                "surface_end"
            ]:
                combined_data.loc[idx - nback, col] = row[col]

    # Add other object information to file
    # Get object position entries whose surface_competitor entry is non-NA
    # and take intersection with object lemmas from audio data whose condition is CONFLICT_LABEL and POS == NOUN_POS_LABEL
    target_words_from_positions = object_positions.loc[object_positions["surface_competitor"].notna(), WORD_FIELD].unique()
    target_words_from_audio = set(
        combined_data.loc[(combined_data["condition"] == CONFLICT_LABEL) & (combined_data[PART_OF_SPEECH_FIELD] == NOUN_POS_LABEL), SPACY_LEMMA_FIELD].unique()
    )
    targets_lc = [word for word in target_words_from_positions if word in target_words_from_audio]
    # Get object position entries whose surface_competitor entry is NA
    # and take intersection with object lemmas from audio data whose condition is NO_CONFLICT_LABEL and POS == NOUN_POS_LABEL
    filler_words_from_positions = object_positions.loc[object_positions["surface_competitor"].isna(), WORD_FIELD].unique()
    filler_words_from_audio = set(
        combined_data.loc[(combined_data["condition"] == NO_CONFLICT_LABEL) & (combined_data[PART_OF_SPEECH_FIELD] == NOUN_POS_LABEL), SPACY_LEMMA_FIELD].unique()
    )
    fillers_lc = [word for word in filler_words_from_positions if word in filler_words_from_audio]

    # There should be exactly 2 unique targets and fillers each
    assert len(fillers_lc) == 2
    assert len(targets_lc) == 2

    pattern = combined_data["pattern"].unique()[0]
    set_id = combined_data["set"].unique()[0]

    def get_surface(text: str, position: str, column: str) -> str:
        """Retrieves the surface column value for a dataframe value matching specified text, position, and pattern."""
        result = object_positions[
            (object_positions[WORD_FIELD] == text) &
            (object_positions["position"] == position) &
            (object_positions["pattern"] == pattern) &
            (object_positions["set"] == set_id)
        ][column]
        return result.iloc[0] if not result.empty else None

    where_is_targets = {target: get_surface(target, position=1, column="surface") for target in targets_lc}
    where_is_comps = {target: get_surface(target, position=1, column="surface_competitor") for target in targets_lc}
    where_is_fillers = {filler: get_surface(filler, position=1, column="surface") for filler in fillers_lc}
    # TODO confirm that where_will* are not needed/used
    # where_will_fillers = {filler: get_surface(filler, position=2, column='surface') for filler in fillers_lc}
    # where_will_targets = {target: get_surface(target, position=2, column='surface') for target in targets_lc}

    # Initialize new columns for surfaces/locations
    N = len(combined_data)
    targetA_surface = [pd.NA] * N
    targetB_surface = [pd.NA] * N
    fillerA_surface = [pd.NA] * N
    fillerB_surface = [pd.NA] * N
    compA_surface = [pd.NA] * N
    compB_surface = [pd.NA] * N
    target_location = [pd.NA] * N

    # Iterate again through words in combined dataframe
    for idx, row in combined_data.iterrows():
        lemma = row[SPACY_LEMMA_FIELD]
        if pd.isna(row["position"]):
            continue

        # Check if preceding word is a determiner
        if combined_data[SPACY_POS_FIELD][idx - 1] == DET_UPOS_TAG:
            nback = 1
        else:
            nback = 2

        if lemma in targets_lc:
            other_target = other_comp = list(setdiff(targets_lc, {lemma}))[0]
            target = comp = list(set(targets_lc).intersection({lemma}))[0]
            # Set values of new columns
            targetA_surface[idx - nback] = targetA_surface[idx] = where_is_targets[target]
            targetB_surface[idx - nback] = targetB_surface[idx] = where_is_targets[other_target]
            compA_surface[idx - nback] = compA_surface[idx] = where_is_comps[comp]
            compB_surface[idx - nback] = compB_surface[idx] = where_is_comps[other_comp]
            fillerA_surface[idx - nback] = fillerA_surface[idx] = where_is_fillers[fillers_lc[0]]
            fillerB_surface[idx - nback] = fillerB_surface[idx] = where_is_fillers[fillers_lc[-1]]

            # TODO why is this update necessary? seems to be just adding the same values again
            where_is_targets[lemma] = row["surface"]
            where_is_comps[lemma] = row["surface_competitor"]
        elif lemma in fillers_lc:
            other_filler = list(setdiff(fillers_lc, {lemma}))[0]
            current_filler = list(set(fillers_lc).intersection({lemma}))[0]
            targetA_surface[idx - nback] = targetA_surface[idx] = where_is_targets[targets_lc[0]]
            targetB_surface[idx - nback] = targetB_surface[idx] = where_is_targets[targets_lc[-1]]
            compA_surface[idx - nback] = compA_surface[idx] = where_is_comps[targets_lc[0]]
            compB_surface[idx - nback] = compB_surface[idx] = where_is_comps[targets_lc[-1]]
            fillerA_surface[idx - nback] = fillerA_surface[idx] = where_is_fillers[current_filler]
            fillerB_surface[idx - nback] = fillerB_surface[idx] = where_is_fillers[other_filler]

    # Add new surface/location columns to dataframe
    combined_data["targetA_surface"] = targetA_surface
    combined_data["targetB_surface"] = targetB_surface
    combined_data["fillerA_surface"] = fillerA_surface
    combined_data["fillerB_surface"] = fillerB_surface
    combined_data["compA_surface"] = compA_surface
    combined_data["compB_surface"] = compB_surface
    combined_data["target_location"] = target_location

    # Set goal/ending locations
    target1 = combined_data[(combined_data[SPACY_LEMMA_FIELD] == targets_lc[0]) & (combined_data[PART_OF_SPEECH_FIELD] == NOUN_POS_LABEL)]
    target1.loc[:, "target_location"] = target1["surface"].shift(-1)
    target1.loc[target1.index[-1], "target_location"] = target1["surface_end"].iloc[0]
    target2 = combined_data[(combined_data[SPACY_LEMMA_FIELD] == targets_lc[-1]) & (combined_data[PART_OF_SPEECH_FIELD] == NOUN_POS_LABEL)]
    target2.loc[:, "target_location"] = target2["surface"].shift(-1)
    target2.loc[target2.index[-1], "target_location"] = target2["surface_end"].iloc[0]
    filler1 = combined_data[(combined_data[SPACY_LEMMA_FIELD] == fillers_lc[0]) & (combined_data[PART_OF_SPEECH_FIELD] == NOUN_POS_LABEL)]
    filler1.loc[:, "target_location"] = filler1["surface"].shift(-1)
    filler1.loc[filler1.index[-1], "target_location"] = filler1["surface_end"].iloc[0]
    filler2 = combined_data[(combined_data[SPACY_LEMMA_FIELD] == fillers_lc[-1]) & (combined_data[PART_OF_SPEECH_FIELD] == NOUN_POS_LABEL)]
    filler2.loc[:, "target_location"] = filler2["surface"].shift(-1)
    filler2.loc[filler2.index[-1], "target_location"] = filler2["surface_end"].iloc[0]
    rest = combined_data[combined_data[PART_OF_SPEECH_FIELD] != NOUN_POS_LABEL]

    # Concatenate filtered dataframes back together once end locations are added
    combined_data = pd.concat([target1, target2, filler1, filler2, rest], axis=0, ignore_index=True)

    # Sort dataframe by LINE_ID_FIELD ("id") column to ensure correct (original) order
    combined_data = combined_data.sort_values(by=WORD_ID_FIELD).reset_index(drop=True)

    # One final iteration through words
    for idx, row in combined_data.iterrows():
        if not pd.isna(row["target_location"]):
            # Check if preceding word is a determiner
            if combined_data[SPACY_POS_FIELD][idx - 1] == DET_UPOS_TAG:
                nback = 1
            else:
                nback = 2
            combined_data.loc[idx - nback, "target_location"] = row["target_location"]

    return combined_data


def main(experiment: str | dict | Experiment) -> Experiment:

    # Initialize DGAME experiment from config
    if not isinstance(experiment, Experiment):
        from dgame.dgame import DGAME
        experiment = DGAME.from_input(experiment)
    logger = experiment.logger

    # Find audio transcript files
    per_subject_audio_files = experiment.get_subject_files_dict(
        dir=experiment.get_transcription_dir(),
        suffix=AUDIO_FILE_SUFFIX,
        recursive=True
    )

    # Retrieve object and filler words
    objects = experiment.objects
    fillers = experiment.fillers

    # Load spaCy NLP pipeline used to tag POS and frequency rank
    spacy_model = experiment.get_dgame_step_parameter(WORDS_PREPROCESS_STEP, "spacy_model", default=DEFAULT_SPACY_MODEL)
    logger.info(f"Loading spaCy model: {spacy_model}")
    nlp_pipeline = load_spacy_pipeline(spacy_model)

    # Process audio files
    skip_indices = experiment.get_dgame_step_parameter(WORDS_PREPROCESS_STEP, "skip_indices")
    for subject_id, audio_files in per_subject_audio_files.items():
        logger.info(f"Processing subject {subject_id}")
        # Reset pattern and set IDs to 1 for each new subject
        pattern_id, set_id = 1, 1
        # Reset trial-number counters for each new subject
        # (trial numbers are unique per subject, continuing across that subject's block files, not reset per block)
        trial_counter_nouns, trial_counter_determiners = 1, 1
        # Load object positions data
        obj_pos_csv = os.path.join(experiment.object_pos_indir, subject_id, OBJECT_POSITIONS_FILE)
        obj_pos_data = experiment.load_object_positions_data(obj_pos_csv)
        # Create subject's audio outdir
        subj_audio_outdir = os.path.join(experiment.outdir, experiment.audio_dir, subject_id)
        os.makedirs(subj_audio_outdir, exist_ok=True)
        for audio_file in sorted(audio_files):
            basename = os.path.basename(audio_file)
            block = re.search(AUDIO_FILE_SUFFIX, basename).group(1)
            audio_outfile = os.path.join(subj_audio_outdir, f"{subject_id}_words_{block}_annotated.csv")
            file_skip_indices = skip_indices.get(os.path.basename(audio_file))
            word_data = preprocess_words_data(
                audio_infile=audio_file,
                nlp_pipeline=nlp_pipeline,
                objects=objects,
                fillers=fillers,
                skip_indices=file_skip_indices,
                pattern_id=pattern_id,
                set_id=set_id,
            )
            combined_data = combine_words_and_obj_position_data(
                word_data=word_data,
                object_positions=obj_pos_data,
            )
            # Assign sequential trial numbers to noun/determiner rows
            combined_data, trial_counter_nouns, trial_counter_determiners = assign_trial_numbers(
                combined_data, trial_counter_nouns, trial_counter_determiners
            )
            # Write output CSV
            combined_data.to_csv(audio_outfile, index=False)
            logger.info(f"Wrote CSV to {audio_outfile}")

            # Increment/adjust pattern and set IDs for next file from same user
            if pattern_id == 2:
                pattern_id = 1
                set_id += 1
            else:
                pattern_id += 1

    return experiment


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Preprocess audio transcript data and combine with object position data.")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(os.path.abspath(args.config))
