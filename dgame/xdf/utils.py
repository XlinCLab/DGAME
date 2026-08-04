import numpy as np
import pandas as pd

from dgame.xdf import SYNC_REFERENCE_STREAM_COLUMN
from dgame.xdf.xdf_stream import XDFStream


def load_stream_sync_reference(sync_reference_file: str) -> pd.DataFrame:
    """Load a stream sync reference CSV, indexed by stream name,
    so callers can look up a given stream's LSL-synchronized
    start/end time by name rather than by position."""
    return pd.read_csv(sync_reference_file, index_col=SYNC_REFERENCE_STREAM_COLUMN)


def extract_audio_stream_channels(audio_stream: XDFStream) -> list[np.ndarray, float]:
    """Extract and separately normalize audio channel samples from a single audio stream."""
    # Extract samples and sampling rate
    samples = audio_stream.time_series.astype(np.float32)
    fs = audio_stream.nominal_srate

    # Ensure shape is (samples, channels)
    # NB: audio streams typically have no per-channel label metadata, so
    # XDFStream.time_series can't disambiguate orientation on its own here
    if samples.shape[0] < samples.shape[1]:
        # likely (channels, samples) -> transpose
        samples = samples.T

    # Extract and normalize each channel independently
    channels = []
    for ch in range(samples.shape[1]):
        channel = samples[:, ch]
        max_val = np.max(np.abs(channel))
        if max_val > 0:
            channel = channel / max_val
        channel_int16 = (channel * 32767).astype(np.int16)
        channels.append(channel_int16)

    return channels, fs
