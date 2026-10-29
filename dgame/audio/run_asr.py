import argparse
import logging
import os

import torch
import whisper

from dgame.audio import MAX_CPU_WORKERS
from dgame.audio.utils import (create_text_grid, get_asr_results_df,
                               get_device, save_transcript, write_asr_df)

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(name)s %(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


def load_asr_model(model_name: str, device: str = None, max_workers: int = None) -> whisper.Whisper:
    """Load a Whisper ASR model."""
    device = get_device(device)
    logger.info(f"Using device: <{device}>")
    if device == "cpu":
        n_threads = MAX_CPU_WORKERS if max_workers is None else min(MAX_CPU_WORKERS, max_workers)
        torch.set_num_threads(n_threads)
        logger.info(f"Worker threads: {n_threads}")

    logger.info(f"Loading Whisper model '{model_name}'...")
    return whisper.load_model(model_name, device=device)


def transcribe_audio(file_path: str, model: whisper.Whisper, language: str) -> dict:
    """Transcribe an audio file with a loaded Whisper model, including word-level timestamps."""
    logger.info(f"Starting ASR transcription on input file {file_path} ...")
    result = model.transcribe(
        file_path,
        language=language,
        task="transcribe",
        word_timestamps=True,
        # # Avoids aggressive context smoothing that may reduce disfluencies
        # condition_on_previous_text=False,
    )
    chunks = [
        {"text": word["word"].strip(), "timestamp": (word["start"], word["end"])}
        for segment in result["segments"]
        for word in segment["words"]
    ]
    return {"text": result["text"], "chunks": chunks}


def main(audio_files: list,
         outdir: str = None,
         model_name: str = "medium",
         language: str = "de",
         max_workers: int = MAX_CPU_WORKERS,
         device: str = None,
         ) -> None:

    # Designate output directory
    if outdir is not None:
        outdir = os.path.abspath(outdir)
    else:
        outdir = os.path.abspath("./asr_out")
    os.makedirs(outdir, exist_ok=True)
    logger.info(f"ASR output directory: {outdir}")

    # Load Whisper model with specified parameters
    model = load_asr_model(
        model_name=model_name,
        device=device,
        max_workers=max_workers,
    )

    # Transcribe each audio file and write transcript, ASR dataframe, and Praat TextGrid to output directory
    n_audios = len(audio_files)
    for audio_path in audio_files:
        result = transcribe_audio(
            file_path=audio_path,
            model=model,
            language=language,
        )
        text = result["text"]
        outfile = os.path.join(
            outdir,
            os.path.basename(audio_path).replace(".wav", ".txt"),
        )
        # Write transcript text file
        save_transcript(text=text, outfile=outfile)

        # Extract word time stamps and write complete ASR dataframe as CSV
        time_chunks = result["chunks"]
        outfile_times = outfile.replace(".txt", ".csv")
        asr_df = get_asr_results_df(time_chunks)
        write_asr_df(asr_df, outfile=outfile_times)

        # Create and write Praat TextGrid
        outfile_textgrid = outfile.replace(".txt", ".TextGrid")
        create_text_grid(asr_df, outfile_textgrid)

    logger.info(f"Completed transcription of {n_audios} audio files.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Transcribe audio with OpenAI Whisper ASR.")
    parser.add_argument(
        "--audio",
        type=str,
        nargs="+",
        required=True,
        help="Path to one or more input .wav files"
    )
    parser.add_argument(
        "--outdir",
        type=str,
        default=None,
        help="Optional path to output directory where transcripts will be written"
    )
    parser.add_argument(
        "--model",
        type=str,
        default="medium",
        help="Whisper model size (tiny, base, small, medium, large)"
    )
    parser.add_argument(
        "--language",
        type=str,
        default="de",
        help="Language code (default: de for German)"
    )
    parser.add_argument(
        "--max-workers",
        type=int,
        default=MAX_CPU_WORKERS,
        help="Maximum number of parallel worker threads for CPU"
    )
    parser.add_argument(
        "--device",
        type=str,
        default=None,
        help="Device (cpu, cuda, mps) on which to run ASR"
    )

    args = parser.parse_args()
    main(
        audio_files=args.audio,
        outdir=args.outdir,
        model_name=args.model,
        language=args.language,
        max_workers=args.max_workers,
        device=args.device,
    )
