import argparse
import os
import re

import pandas as pd

from dgame.constants import (AUDIO_ERP_TRIALTIME_FILE_SUFFIX, DET_POS_LABEL,
                             DIRECTION_WORD_LABEL, NOUN_POS_LABEL,
                             PART_OF_SPEECH_FIELD, VERB_POS_LABEL,
                             WORD_END_FIELD, WORD_FIELD, WORD_ONSET_FIELD)
from dgame.pipeline import STEP_DA_KEY, STEP_E_KEY
from experiment.load_experiment import Experiment


def assign_trial_context(df: pd.DataFrame, gap_threshold: float = 2.0) -> pd.DataFrame:
    """Assign trial membership to each row by extending D/N core trials to surrounding context.
    Preamble and postamble rows are absorbed greedily until a temporal gap exceeding
    gap_threshold seconds is encountered.

    Assumptions about trial segmentation:
    - Core rows are those tagged pos='D' or pos='N' with a nonzero trial number.
    - A trial consists of: (a) preamble = adjacent/word rows immediately
        before the D, (b) D-N core, (c) postamble = adjacent/word rows immediately
        after the N until the next D or until a large temporal gap.
    - A "large gap" is defined via gap_threshold (default: 2.0 s
        between consecutive words). This prevents temporally distant words from
        being incorrectly assigned to a trial.
    """
    df = df.sort_values(WORD_ONSET_FIELD).reset_index(drop=True)

    core_mask = (
        df[PART_OF_SPEECH_FIELD].isin([DET_POS_LABEL, NOUN_POS_LABEL])
        & df["trial"].notna()
        & (df["trial"] != 0)
    )

    # Save positional index as a column so groupby agg can reference it
    core_df = df[core_mask].copy()
    core_df["_idx"] = core_df.index
    trial_cores = (
        core_df.groupby("trial")
        .agg(core_start_idx=("_idx", "min"), core_end_idx=("_idx", "max"))
        .reset_index()
        .sort_values("core_start_idx")
    )

    trial_assigned = [pd.NA] * len(df)

    # Assign core rows
    for _, row in trial_cores.iterrows():
        tr = row["trial"]
        for i in range(int(row["core_start_idx"]), int(row["core_end_idx"]) + 1):
            trial_assigned[i] = tr

    n = len(df)

    # Preamble: walk backwards from each core start
    for _, row in trial_cores.iterrows():
        tr = row["trial"]
        j = int(row["core_start_idx"]) - 1
        while j >= 0 and pd.isna(trial_assigned[j]):
            gap = df.at[j + 1, WORD_ONSET_FIELD] - df.at[j, WORD_END_FIELD]
            if pd.isna(gap) or gap > gap_threshold:
                break
            trial_assigned[j] = tr
            j -= 1

    # Postamble: walk forwards from each core end
    for _, row in trial_cores.iterrows():
        tr = row["trial"]
        j = int(row["core_end_idx"]) + 1
        while j < n and pd.isna(trial_assigned[j]):
            gap = df.at[j, WORD_ONSET_FIELD] - df.at[j - 1, WORD_END_FIELD]
            if pd.isna(gap) or gap > gap_threshold:
                break
            trial_assigned[j] = tr
            j += 1

    df["trial_assigned"] = trial_assigned
    return df


def classify_slot(text: str,
                  pos: str,
                  direction_lemmas: set[str],
                  verb_lemmas: set[str]) -> str:
    """Map a word to its abstract slot label (determiner, noun, direction word, or verb)."""
    if pos == DET_POS_LABEL:
        return DET_POS_LABEL
    if pos == NOUN_POS_LABEL:
        return NOUN_POS_LABEL
    text_l = str(text).lower()
    if text_l in direction_lemmas:
        return DIRECTION_WORD_LABEL
    if verb_lemmas and text_l in verb_lemmas:
        return VERB_POS_LABEL
    return text_l


def aggregate_trials(df: pd.DataFrame,
                     direction_lemmas: set[str],
                     verb_lemmas: set[str]) -> pd.DataFrame:
    """Collapse word-level rows into one row per trial with wording and slot pattern."""
    df_assigned = df[df["trial_assigned"].notna()].sort_values(WORD_ONSET_FIELD)

    rows = []
    for (subject_id, source_file, trial), g in df_assigned.groupby(
        ["subject_id", "source_file", "trial_assigned"]
    ):
        wording = " ".join(str(t) for t in g[WORD_FIELD])
        wording = re.sub(r"\s+", " ", wording)
        pattern = " ".join(
            classify_slot(str(t), str(p), direction_lemmas, verb_lemmas)
            for t, p in zip(g[WORD_FIELD], g[PART_OF_SPEECH_FIELD])
        )
        pattern = re.sub(r"\s+", " ", pattern)
        condition_vals = g["condition"].dropna()
        rows.append({
            "subject_id": subject_id,
            "source_file": source_file,
            "trial": int(trial),
            "n_words": len(g),
            "t_start": g[WORD_ONSET_FIELD].min(),
            "t_end": g[WORD_END_FIELD].max(),
            "duration": g[WORD_END_FIELD].max() - g[WORD_ONSET_FIELD].min(),
            "condition": condition_vals.iloc[0] if len(condition_vals) > 0 else None,
            "wording": wording,
            "pattern": pattern,
            "has_D": (g[PART_OF_SPEECH_FIELD] == DET_POS_LABEL).any(),
            "has_N": (g[PART_OF_SPEECH_FIELD] == NOUN_POS_LABEL).any(),
        })

    return pd.DataFrame(rows)


def summarize_duration_stats(df: pd.DataFrame) -> dict:
    return {
        "n_trials": len(df),
        "mean_s": df["duration"].mean(),
        "median_s": df["duration"].median(),
        "sd_s": df["duration"].std(),
        "min_s": df["duration"].min(),
        "max_s": df["duration"].max(),
        "mean_words": df["n_words"].mean(),
    }


def main(experiment: str | dict | Experiment) -> Experiment:

    # Initialize DGAME experiment from config
    if not isinstance(experiment, Experiment):
        from dgame.dgame import DGAME
        experiment = DGAME.from_input(experiment)
    logger = experiment.logger

    # Load parameters from experiment config
    gap_threshold = experiment.get_dgame_step_parameter(STEP_E_KEY, "gap_threshold")
    n_examples = experiment.get_dgame_step_parameter(STEP_E_KEY, "n_examples")
    direction_lemmas = experiment.get_dgame_step_parameter(STEP_E_KEY, "direction_lemmas")
    direction_lemmas = set(dir_lemma.lower() for dir_lemma in direction_lemmas)
    # NB: verb_lemmas by default empty unless explicitly set
    verb_lemmas = experiment.get_dgame_step_parameter(STEP_E_KEY, "verb_lemmas")
    verb_lemmas = set(verb_lemma.lower() for verb_lemma in verb_lemmas)

    # Find words2erp trialtime files per subject written by step Da
    per_subject_erp_files = experiment.get_subject_files_dict(
        dir=experiment.audio_outdir,
        suffix=AUDIO_ERP_TRIALTIME_FILE_SUFFIX,
        recursive=True,
    )

    # Write all outputs to a shared instruction_patterns subdirectory under audio_outdir
    instruction_pattern_outdir = os.path.join(experiment.audio_outdir, "instruction_patterns")
    os.makedirs(instruction_pattern_outdir, exist_ok=True)

    # Process each subject's words2erp files
    all_trial_dfs = []
    for subject_id, erp_files in per_subject_erp_files.items():
        logger.info(f"Processing subject <{subject_id}>...")
        for erp_file in sorted(erp_files):
            logger.debug(f"Processing file: {os.path.basename(erp_file)}")
            df = pd.read_csv(erp_file)
            df["source_file"] = os.path.basename(erp_file)
            df["subject_id"] = subject_id
            df = assign_trial_context(df, gap_threshold=gap_threshold)
            trial_df = aggregate_trials(df, direction_lemmas=direction_lemmas, verb_lemmas=verb_lemmas)
            all_trial_dfs.append(trial_df)

    if not all_trial_dfs:
        logger.warning(
            f"No words2erp files found; skipping DGAME step <{STEP_E_KEY}>. Ordinarily these files are created by step <{STEP_DA_KEY}>."
        )
        return experiment

    all_trials = pd.concat(all_trial_dfs, ignore_index=True)

    # Flag and drop trials missing a D-N pair
    problematic = all_trials[~(all_trials["has_D"] & all_trials["has_N"])]
    if len(problematic) > 0:
        problematic_csv = os.path.join(instruction_pattern_outdir, "problematic_trials.csv")
        logger.warning(f"{len(problematic)} trial(s) without a complete D-N pair. See details in {problematic_csv}")
        problematic.to_csv(problematic_csv, index=False)
        all_trials = all_trials[all_trials["has_D"] & all_trials["has_N"]]

    # instructions_full.csv: one row per trial with full wording and slot pattern
    cols_full = ["subject_id", "source_file", "trial", "condition", "n_words", "duration", "wording", "pattern"]
    all_trials[cols_full].to_csv(os.path.join(instruction_pattern_outdir, "instructions_full.csv"), index=False)

    # patterns_abstract.csv: trial × abstract pattern only
    all_trials[["subject_id", "trial", "condition", "pattern"]].to_csv(
        os.path.join(instruction_pattern_outdir, "patterns_abstract.csv"), index=False
    )

    # pattern_counts.csv: frequency table of abstract patterns
    pattern_counts = (
        all_trials.groupby("pattern")
        .size()
        .reset_index(name="n")
        .sort_values(["n", "pattern"], ascending=[False, True])
        .reset_index(drop=True)
    )
    pattern_counts["rel_freq"] = pattern_counts["n"] / pattern_counts["n"].sum()
    pattern_counts.to_csv(os.path.join(instruction_pattern_outdir, "pattern_counts.csv"), index=False)

    # pattern_examples.csv: up to n_examples example wordings per pattern
    pattern_examples = (
        all_trials.groupby("pattern")
        .head(n_examples)
        .reset_index(drop=True)
        [["pattern", "wording", "subject_id", "trial"]]
        .sort_values(["pattern", "wording"])
        .reset_index(drop=True)
    )
    pattern_examples.to_csv(os.path.join(instruction_pattern_outdir, "pattern_examples.csv"), index=False)

    # Duration statistics
    duration_overall = pd.DataFrame([summarize_duration_stats(all_trials)])
    duration_overall.to_csv(os.path.join(instruction_pattern_outdir, "duration_overall.csv"), index=False)

    duration_by_subject = pd.DataFrame([
        {"subject_id": sid, **summarize_duration_stats(g)}
        for sid, g in all_trials.groupby("subject_id")
    ])
    duration_by_subject.to_csv(os.path.join(instruction_pattern_outdir, "duration_by_subject.csv"), index=False)

    trials_with_condition = all_trials.dropna(subset=["condition"])
    if len(trials_with_condition) > 0:
        duration_by_condition = pd.DataFrame([
            {"condition": cond, **summarize_duration_stats(g)}
            for cond, g in trials_with_condition.groupby("condition")
        ])
        duration_by_condition.to_csv(os.path.join(instruction_pattern_outdir, "duration_by_condition.csv"), index=False)

    logger.info(f"Processed {len(all_trials)} trials across {len(pattern_counts)} unique patterns.")
    logger.info(f"Output written to: {instruction_pattern_outdir}")

    return experiment


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Extract Director Task instructions from word-based CSVs and derive abstract syntactic patterns.")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(os.path.abspath(args.config))
