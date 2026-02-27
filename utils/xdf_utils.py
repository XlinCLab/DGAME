import numpy as np
from pyxdf import load_xdf

STREAM_TIMESTAMPS_LABEL = "time_stamps"


def get_xdf_stream(
        stream_label: str,
        xdf_file: str = None,
        xdf_data: list = None,
        verbose: bool = False,
        **kwargs
        ) -> dict:
    """Fetches a specific stream by its label from an XDF file (optionally preloaded)."""
    if xdf_data is None:
        assert xdf_file is not None, "xdf_file argument is required if no xdf_data argument is provided"
        xdf_data, _ = load_xdf(xdf_file, verbose=verbose, **kwargs)
    stream_idx = None
    for idx, stream in enumerate(xdf_data):
        stream_name = stream['info']['name'][0]
        if stream_name == stream_label:
            stream_idx = idx
            break
    if stream_idx is None:
        raise ValueError(f"No <{stream_label}> stream found in {xdf_file}")
    stream = xdf_data[stream_idx]
    return stream


def extract_audio_stream_channels(audio_stream: dict) -> list[np.ndarray, float]:
    """Extract and separately normalize audio channel samples from a single audio stream."""
    # Extract samples and sampling rate
    samples = np.array(audio_stream['time_series'], dtype=np.float32)
    fs = float(audio_stream['info']['nominal_srate'][0])

    # Ensure shape is (samples, channels)
    if samples.ndim == 1:
        samples = samples[:, None]  # mono -> (N,1)
    elif samples.shape[0] < samples.shape[1]:
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


def get_relative_times_from_stream(
        stream: dict,
        round_n: int = None
        ) -> np.array:
    ts = np.array(stream[STREAM_TIMESTAMPS_LABEL], dtype=np.float64)

    if len(ts) == 0:
        return ts

    rel = ts - ts[0]

    if round_n is not None:
        rel = np.round(rel, decimals=round_n)

    return rel
