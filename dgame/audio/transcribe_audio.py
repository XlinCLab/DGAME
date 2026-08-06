import argparse
import gc
import os

import pandas as pd
import torch

from dgame.audio.run_asr import load_asr_model, transcribe_audio
from dgame.audio.utils import (create_text_grid, get_asr_results_df,
                               save_transcript)
from dgame.constants import BLOCK_IDS
from dgame.pipeline import TRANSCRIBE_AUDIO_STEP
from dgame.words import (INPUT_LINE_ID_FIELD, INPUT_WORD_ONSET_FIELD,
                         WORD_END_FIELD, WORD_FIELD)
from experiment.input_validation import assert_output_file_exists
from experiment.load_experiment import Experiment


def get_words_df(chunks: list[dict]) -> pd.DataFrame:
    """Convert Whisper word-level chunks into the "words" CSV format expected by
    step words.annotate_words: one row per word, with the recognized word text
    and its onset/offset timestamps."""
    words = [
        {
            INPUT_LINE_ID_FIELD: idx + 1,
            INPUT_WORD_ONSET_FIELD: chunk["timestamp"][0],
            WORD_FIELD: chunk["text"],
            WORD_END_FIELD: chunk["timestamp"][-1],
        }
        for idx, chunk in enumerate(chunks)
    ]
    return pd.DataFrame(words, columns=[INPUT_LINE_ID_FIELD, INPUT_WORD_ONSET_FIELD, WORD_FIELD, WORD_END_FIELD])


def validate_outputs(experiment, subject_ids: list) -> None:
    """Validate word transcript outputs."""
    director_label = experiment.director_label
    for subject_id in subject_ids:
        subject_words_dir = experiment.get_role_outdir(experiment.preproc_audio_indir, subject_id, director_label)
        for block in BLOCK_IDS:
            words_file = os.path.join(subject_words_dir, f"{subject_id}_words_{block}.csv")
            assert_output_file_exists(words_file)


def main(experiment: str | dict | Experiment) -> Experiment:
    # Initialize DGAME experiment from config
    if not isinstance(experiment, Experiment):
        from dgame.dgame import DGAME
        experiment = DGAME.from_input(experiment)
    logger = experiment.logger

    language = experiment.get_dgame_step_parameter(TRANSCRIBE_AUDIO_STEP, "language", default="de")
    model_size = experiment.get_dgame_step_parameter(TRANSCRIBE_AUDIO_STEP, "model_size", default="medium")
    device = experiment.get_dgame_step_parameter(TRANSCRIBE_AUDIO_STEP, "device", default=None)
    max_workers = experiment.get_dgame_step_parameter(TRANSCRIBE_AUDIO_STEP, "max_workers", default=None)

    model = load_asr_model(
        model_name=model_size,
        device=device,
        max_workers=max_workers,
    )

    # Transcribe the director's audio channel (the participant giving verbal instructions)
    # for each subject and block, writing one "words" CSV per subject per block
    failed_files = []
    for subject_id in experiment.subject_ids:
        subject_audio_outdir = experiment.get_role_outdir(
            experiment.audio_outdir,
            subject_id,
            experiment.director_label,
        )
        os.makedirs(subject_audio_outdir, exist_ok=True)
        for block in BLOCK_IDS:
            audio_basename = f"{subject_id}_{experiment.director_label}_{block}"
            audio_file = os.path.join(subject_audio_outdir, f"{audio_basename}.wav")
            try:
                asr_result = transcribe_audio(
                    file_path=audio_file,
                    model=model,
                    language=language,
                )
                chunks = asr_result["chunks"]

                # Write required "words" CSV, consumed by step words.annotate_words
                words_df = get_words_df(chunks)
                words_outfile = os.path.join(subject_audio_outdir, f"{subject_id}_words_{block}.csv")
                words_df.to_csv(words_outfile, index=False)
                logger.info(f"Wrote word transcript to {words_outfile}")

                # Write plain text transcript and Praat TextGrid alongside the source audio
                save_transcript(asr_result["text"], os.path.join(subject_audio_outdir, f"{audio_basename}.txt"))
                create_text_grid(get_asr_results_df(chunks), os.path.join(subject_audio_outdir, f"{audio_basename}.TextGrid"))
            except Exception as exc:
                # ASR can fail on individual audio content regardless of settings;
                # don't let one bad file abort transcription for every other subject/block.
                # Instead log it and continue
                logger.error(f"Transcription failed for {audio_file}: {exc}")
                failed_files.append(audio_file)
                continue
            finally:
                # Release this file's tensors before next file to ensure memory is reclaimed
                if "asr_result" in locals():
                    del asr_result
                if "words_df" in locals():
                    del words_df
                gc.collect()
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()
                elif torch.backends.mps.is_available():
                    torch.mps.empty_cache()

    if failed_files:
        logger.warning(
            f"Transcription failed for {len(failed_files)} file(s), skipped: "
            + ", ".join(failed_files)
        )

    # Validate outputs
    validate_outputs(experiment, experiment.subject_ids)

    return experiment


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Transcribe director audio with OpenAI Whisper ASR.")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(os.path.abspath(args.config))
