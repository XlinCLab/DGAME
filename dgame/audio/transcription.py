import logging
import os

import pandas as pd
import textgrids as tg

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(name)s %(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


def save_transcript(text: str, outfile: str) -> None:
    """Save a plain text ASR transcript file."""
    outdir = os.path.abspath(os.path.dirname(outfile))
    os.makedirs(outdir, exist_ok=True)
    with open(outfile, "w", encoding="utf-8") as f:
        f.write(text)
    logger.info(f"Saved transcript to: {outfile}")


def get_asr_results_df(time_chunks: list) -> pd.DataFrame:
    """
    Change format of ASR results object from:
    {'text': str, 'timestamp': tuple(start, end)}
    to:
    {'text': str, 'start': float, 'end': float, 'duration': float}
    and then convert to Pandas DataFrame.
    """
    time_chunks = [
        {
            'text': chunk['text'],
            'start': chunk['timestamp'][0],
            'end': chunk['timestamp'][-1],
            'duration': chunk['timestamp'][-1] - chunk['timestamp'][0],
        }
        for chunk in time_chunks
    ]
    return pd.DataFrame(time_chunks)


def write_asr_df(asr_df: pd.DataFrame, outfile: str) -> None:
    """Write a Pandas DataFrame of ASR results to CSV."""
    outdir = os.path.abspath(os.path.dirname(outfile))
    os.makedirs(outdir, exist_ok=True)
    asr_df.to_csv(outfile, index=False)
    logger.info(f"Saved full ASR results to: {outfile}")


def create_text_grid(asr_df: pd.DataFrame,
                     textgrid_outfile: str | None = None,
                     ) -> tg.TextGrid:
    """Create a Praat TextGrid from a Pandas DataFrame of ASR results
    (recognized words with start and end times.)"""

    # Overall time range: start from 0.0 so the TextGrid aligns with the audio file
    xmin = 0.0
    xmax = float(asr_df["end"].max())

    # Create TextGrid and Tier with correct bounds
    text_grid = tg.TextGrid()
    text_grid.xmin = xmin
    text_grid.xmax = xmax

    tier = tg.Tier()
    tier.xmin = xmin
    tier.xmax = xmax

    last_end = xmin

    for _, row in asr_df.iterrows():
        start = float(row["start"])
        end = float(row["end"])
        text = str(row["text"])

        # Skip segments whose end is earlier than the start
        # ASR timestamp error
        if end < start:
            logger.warning(
                f"Skipping incorrectly timed segment <{text}>: start: {start} | end: {end}"
            )
            continue

        # Skip words that overlap the previous interval
        if start < last_end - 1e-6:
            logger.warning(
                f"Skipping overlapping word <{text}> at {start:.3f}-{end:.3f}s "
                f"(previous interval ended at {last_end:.3f}s)"
            )
            continue

        # Fill gaps with empty intervals
        if start > last_end:
            tier.append(
                tg.Interval(
                    text="",
                    xmin=last_end,
                    xmax=start,
                )
            )

        tier.append(
            tg.Interval(
                text=text,
                xmin=start,
                xmax=end,
            )
        )

        last_end = end

    text_grid["words"] = tier

    if textgrid_outfile:
        text_grid.write(textgrid_outfile)
        logger.info(f"Wrote Praat TextGrid to {textgrid_outfile}")

    return text_grid
