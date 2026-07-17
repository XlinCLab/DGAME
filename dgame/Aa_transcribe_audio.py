import argparse
import gc
import importlib.util
import logging
import multiprocessing
import os

# Must be set before torch initializes its CUDA allocator. Reduces OOM risk from allocator
# fragmentation (this ASR pipeline processes many long audio files in a single long-lived
# process, and the CUDA allocator's own troubleshooting message recommends this setting)
os.environ.setdefault("PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True")

import numpy as np
import pandas as pd
import textgrids as tg
import torch
from dotenv import load_dotenv
from huggingface_hub import login
from transformers import AutoModelForSpeechSeq2Seq, AutoProcessor, Pipeline, pipeline

from dgame.A_export_audio_and_et_times import get_role_outdir
from dgame.constants import BLOCK_IDS, DIRECTOR_LABEL, SCRIPT_DIR, STEP_AA_KEY
from experiment.input_validation import assert_output_file_exists
from experiment.load_experiment import Experiment

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(name)s %(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

MODEL_ID = "nyrahealth/CrisperWhisper"
MAX_CPU_WORKERS = multiprocessing.cpu_count()

# CrisperWhisper itself is not a pip-installable package: it is expected to be checked
# out as a sibling directory of this repo (i.e. "../CrisperWhisper" relative to the
# DGAME repo root). Its "utils" module is loaded directly from that file path (rather
# than via sys.path + a plain import) to avoid colliding with DGAME's own "utils" package,
# which is already present in sys.modules by the time this step runs.
CRISPERWHISPER_REPO_DIR = os.path.join(os.path.dirname(SCRIPT_DIR), "..", "CrisperWhisper")


def _load_adjust_pauses_function():
    """Load CrisperWhisper's `adjust_pauses_for_hf_pipeline_output` from its sibling repo checkout."""
    crisperwhisper_utils_path = os.path.join(CRISPERWHISPER_REPO_DIR, "utils.py")
    if not os.path.isfile(crisperwhisper_utils_path):
        raise ImportError(
            f"Could not find CrisperWhisper's utils.py at {crisperwhisper_utils_path}. "
            "Clone https://github.com/nyrahealth/CrisperWhisper as a sibling directory of this repo."
        )
    spec = importlib.util.spec_from_file_location("crisperwhisper_utils", crisperwhisper_utils_path)
    crisperwhisper_utils = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(crisperwhisper_utils)
    return crisperwhisper_utils.adjust_pauses_for_hf_pipeline_output


# Mapping of devices to float types
DEVICE_DTYPE_MAP = {
    "cuda:0": torch.float16,
    "mps": torch.float16,
    "cpu": torch.float32,
}
DEVICES = {"cuda", "mps", "cpu"}


def get_device(device: str = None) -> tuple[str, np.dtype]:
    """Select device and corresponding torch float type."""
    if device is None:
        if torch.cuda.is_available():
            device = "cuda:0"
        elif torch.backends.mps.is_available():
            device = "mps"
        else:
            device = "cpu"
    elif device.lower() not in DEVICES:
        raise ValueError(f"Unrecognized device '{device}'. Please specify one of {DEVICES}.")
    torch_dtype = DEVICE_DTYPE_MAP.get(device, torch.float32)
    return device, torch_dtype


def init_asr_pipeline(language: str,
                      repetition_penalty: float,
                      max_workers: int = MAX_CPU_WORKERS,
                      device: str = None,
                      batch_size: int = None,
                      num_beams: int = 1,
                      ) -> Pipeline:
    """Initialize a CrisperWhisper ASR pipeline."""
    device, torch_dtype = get_device(device)
    logger.info(f"Using device: <{device}>")
    if device == "cpu":
        n_threads = MAX_CPU_WORKERS if max_workers is None else min(MAX_CPU_WORKERS, max_workers)
        torch.set_num_threads(n_threads)
        logger.info(f"Worker threads: {n_threads}")

    model = AutoModelForSpeechSeq2Seq.from_pretrained(
        MODEL_ID,
        torch_dtype=torch_dtype,
        low_cpu_mem_usage=True,
        use_safetensors=True,
        attn_implementation="eager",
    )
    model.to(device)

    processor = AutoProcessor.from_pretrained(MODEL_ID)

    # batch_size=1 is more efficient on CPU/MPS; can be raised via the "batch_size" step
    # config parameter on GPUs with more headroom
    if batch_size is None:
        batch_size = 1

    logger.info("Initializing ASR pipeline...")
    asr_pipeline = pipeline(
        "automatic-speech-recognition",
        model=model,
        tokenizer=processor.tokenizer,
        feature_extractor=processor.feature_extractor,
        chunk_length_s=30,
        batch_size=batch_size,
        return_timestamps="word",
        torch_dtype=torch_dtype,
        device=device,
        generate_kwargs={
            "language": language,
            "task": "transcribe",
            "repetition_penalty": repetition_penalty,
            # CrisperWhisper's default generation_config uses num_beams=5: each beam keeps
            # its own growing KV-cache, which (combined with attn_implementation="eager",
            # required for word-level timestamp alignment) is the main driver of GPU memory
            # use here, not batch_size. Defaulting to greedy decoding (num_beams=1) cuts
            # that ~5x; raise it again via the "num_beams" step config parameter if your
            # GPU has enough headroom and you want beam search's usual quality benefit
            "num_beams": num_beams,
        },
    )
    return asr_pipeline


def transcribe_audio(file_path: str, asr_pipeline: Pipeline, adjust_pauses) -> dict:
    """Transcribe an audio file with an initialized ASR pipeline."""
    logger.info(f"Starting ASR transcription on input file {file_path} ...")
    hf_output = asr_pipeline(file_path)
    return adjust_pauses(hf_output)


def get_words_df(asr_result: dict) -> pd.DataFrame:
    """Convert an ASR result's word-level chunks into the "words" CSV format expected by
    step B_prepare_words: one row per word, with columns "line" (1-based word index),
    "tmin" (word onset), "text" (word), and "tmax" (word offset)."""
    words = [
        {
            "line": idx + 1,
            "tmin": chunk["timestamp"][0],
            "text": chunk["text"].strip(),
            "tmax": chunk["timestamp"][-1],
        }
        for idx, chunk in enumerate(asr_result["chunks"])
    ]
    return pd.DataFrame(words, columns=["line", "tmin", "text", "tmax"])


def save_transcript(text: str, outfile: str) -> None:
    """Save a plain text ASR transcript file."""
    with open(outfile, "w", encoding="utf-8") as f:
        f.write(text)
    logger.info(f"Saved transcript to: {outfile}")


def create_text_grid(words_df: pd.DataFrame, textgrid_outfile: str) -> None:
    """Create a Praat TextGrid from word-level ASR results, filling gaps between words with empty intervals."""
    xmin = float(words_df["tmin"].min())
    xmax = float(words_df["tmax"].max())

    text_grid = tg.TextGrid()
    text_grid.xmin = xmin
    text_grid.xmax = xmax

    tier = tg.Tier()
    tier.xmin = xmin
    tier.xmax = xmax

    last_end = xmin
    for _, row in words_df.iterrows():
        start, end, text = float(row["tmin"]), float(row["tmax"]), str(row["text"])
        # Skip words whose end is earlier than their start (ASR timestamp error)
        if end < start:
            logger.warning(f"Skipping incorrectly timed word <{text}>: start: {start} | end: {end}")
            continue
        # Fill gaps with empty intervals
        if start > last_end:
            tier.append(tg.Interval(text="", xmin=last_end, xmax=start))
        tier.append(tg.Interval(text=text, xmin=start, xmax=end))
        last_end = end

    text_grid["words"] = tier
    text_grid.write(textgrid_outfile)
    logger.info(f"Wrote Praat TextGrid to {textgrid_outfile}")


def validate_outputs(experiment, subject_ids: list) -> None:
    """Validate word transcript outputs."""
    for subject_id in subject_ids:
        subject_audio_dir = os.path.join(experiment.preproc_audio_indir, subject_id)
        for block in BLOCK_IDS:
            words_file = os.path.join(subject_audio_dir, f"{subject_id}_words_{block}.csv")
            assert_output_file_exists(words_file)


def main(experiment: str | dict | Experiment) -> Experiment:
    # Initialize DGAME experiment from config
    if not isinstance(experiment, Experiment):
        from dgame.dgame import DGAME
        experiment = DGAME.from_input(experiment)
    logger = experiment.logger

    # Load Hugging Face access token from environment variable, if set
    load_dotenv()
    hf_access_token = os.getenv("HF_ACCESS_TOKEN")
    if hf_access_token is not None:
        login(hf_access_token)

    adjust_pauses = _load_adjust_pauses_function()

    language = experiment.get_dgame_step_parameter(STEP_AA_KEY, "language", default="de")
    repetition_penalty = experiment.get_dgame_step_parameter(STEP_AA_KEY, "repetition_penalty", default=1.1)
    max_workers = experiment.get_dgame_step_parameter(STEP_AA_KEY, "max_workers", default=MAX_CPU_WORKERS)
    device = experiment.get_dgame_step_parameter(STEP_AA_KEY, "device", default=None)
    batch_size = experiment.get_dgame_step_parameter(STEP_AA_KEY, "batch_size", default=None)
    num_beams = experiment.get_dgame_step_parameter(STEP_AA_KEY, "num_beams", default=1)

    asr_pipeline = init_asr_pipeline(
        language=language,
        repetition_penalty=repetition_penalty,
        max_workers=max_workers,
        device=device,
        batch_size=batch_size,
        num_beams=num_beams,
    )

    # Transcribe the director's audio channel (the participant giving verbal instructions)
    # for each subject and block, writing one "words" CSV per subject per block
    failed_files = []
    for subject_id in experiment.subject_ids:
        subject_words_outdir = os.path.join(experiment.preproc_audio_indir, subject_id)
        os.makedirs(subject_words_outdir, exist_ok=True)
        director_audio_dir = get_role_outdir(
            experiment.audio_outdir, subject_id, DIRECTOR_LABEL, experiment.dgame_version
        )
        for block in BLOCK_IDS:
            audio_basename = f"{subject_id}_{DIRECTOR_LABEL}_{block}"
            audio_file = os.path.join(director_audio_dir, f"{audio_basename}.wav")
            try:
                asr_result = transcribe_audio(audio_file, asr_pipeline, adjust_pauses)
                words_df = get_words_df(asr_result)

                # Write required "words" CSV, consumed by step B_prepare_words
                words_outfile = os.path.join(subject_words_outdir, f"{subject_id}_words_{block}.csv")
                words_df.to_csv(words_outfile, index=False)
                logger.info(f"Wrote word transcript to {words_outfile}")

                # Write plain text transcript and Praat TextGrid alongside the source audio, for reference/QA
                save_transcript(asr_result["text"], os.path.join(director_audio_dir, f"{audio_basename}.txt"))
                create_text_grid(words_df, os.path.join(director_audio_dir, f"{audio_basename}.TextGrid"))
            except Exception as exc:
                # ASR (word-timestamp DTW alignment specifically) can fail on individual
                # audio content regardless of settings (observed independent of num_beams).
                # Don't let one bad file abort transcription for every other subject/block -
                # log it and keep going; validate_outputs() below will surface exactly which
                # files are still missing once the run finishes.
                logger.error(f"Transcription failed for {audio_file}: {exc}")
                failed_files.append(audio_file)
                continue
            finally:
                # Release this file's GPU tensors (cross-attention weights etc.) before the
                # next one: this is a single long-lived process transcribing many files
                # back-to-back, and neither Python's refcounting nor torch's caching allocator
                # reliably reclaims that memory on their own between iterations
                if "asr_result" in locals():
                    del asr_result
                if "words_df" in locals():
                    del words_df
                gc.collect()
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()

    if failed_files:
        logger.warning(
            f"Transcription failed for {len(failed_files)} file(s), skipped: "
            + ", ".join(failed_files)
        )

    # Validate outputs
    validate_outputs(experiment, experiment.subject_ids)

    return experiment


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Transcribe director audio with CrisperWhisper ASR.")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(os.path.abspath(args.config))
