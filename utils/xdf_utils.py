import json
import os

import numpy as np
from pyxdf import load_xdf

STREAM_TIMESTAMPS_LABEL = "time_stamps"
# Stream type (not name) under which the EEG stream is registered in the XDF file
EEG_STREAM_TYPE = "EEG"

STREAM_ORIGINS_FILENAME_TEMPLATE = "{subject_id}_stream_origins_{block}.json"


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


def get_xdf_stream_by_type(
        stream_type: str,
        xdf_file: str = None,
        xdf_data: list = None,
        verbose: bool = False,
        **kwargs
        ) -> dict:
    """Fetches a specific stream by its type from an XDF file (optionally preloaded)."""
    if xdf_data is None:
        assert xdf_file is not None, "xdf_file argument is required if no xdf_data argument is provided"
        xdf_data, _ = load_xdf(xdf_file, verbose=verbose, **kwargs)
    stream_idx = None
    for idx, stream in enumerate(xdf_data):
        stream_type_val = stream.get('info', {}).get('type', "")
        if isinstance(stream_type_val, list):
            stream_type_val = stream_type_val[0] if stream_type_val else ""
        if str(stream_type_val).lower() == str(stream_type).lower():
            stream_idx = idx
            break
    if stream_idx is None:
        raise ValueError(f"No <{stream_type}> stream found in {xdf_file}")
    stream = xdf_data[stream_idx]
    return stream


def extract_stream_labels(stream: dict) -> list[str]:
    """Extract channel labels from a stream's metadata."""
    desc = stream.get("info", {}).get("desc", [])
    if isinstance(desc, list) and len(desc) > 0:
        desc = desc[0]
    channels = desc.get("channels", {}) if isinstance(desc, dict) else {}
    if isinstance(channels, list) and len(channels) > 0:
        channels = channels[0]
    channels = channels.get("channel", []) if isinstance(channels, dict) else []
    labels = []
    for ch in channels:
        if isinstance(ch, dict) and "label" in ch:
            label = ch["label"]
            if isinstance(label, list):
                label = label[0] if label else ""
            labels.append(str(label))
        else:
            labels.append("")
    return labels


def extract_eeg_stream_samples(eeg_stream: dict) -> tuple[np.ndarray, float, list[str], np.ndarray]:
    """Extract EEG samples as (channels, samples), sampling rate, labels, and per-sample timestamps.

    The timestamps are returned (rather than discarded, as previously) so that callers can
    anchor the EEG recording to the same LSL-synchronized clock used for the other streams
    (audio, eyetracker) in the same XDF file, instead of implicitly assuming the EEG stream's
    first sample is time zero.
    """
    samples = np.array(eeg_stream['time_series'], dtype=np.float64)
    if samples.ndim == 1:
        samples = samples[:, None]

    srate = eeg_stream.get('info', {}).get('nominal_srate', 0)
    if isinstance(srate, list):
        srate = srate[0] if srate else 0
    srate = float(srate)
    labels = extract_stream_labels(eeg_stream)
    if len(labels) == samples.shape[1] and len(labels) != samples.shape[0]:
        samples = samples.T
    time_stamps = np.array(eeg_stream[STREAM_TIMESTAMPS_LABEL], dtype=np.float64)
    return samples, srate, labels, time_stamps


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


def get_stream_first_timestamp(stream: dict) -> float:
    """Return a stream's first LSL clock-synchronized timestamp.

    This is a stream's "origin" on the clock shared by every stream in the same XDF file
    (when loaded with `synchronize_clocks=True`, which is pyxdf's default). Subtracting one
    stream's first timestamp from another's yields the real-world offset between when each
    device started streaming, which is exactly the information that is lost if a stream's
    timestamps are instead only ever expressed relative to its own first sample.
    """
    return float(stream[STREAM_TIMESTAMPS_LABEL][0])


def get_times_relative_to_reference(
        stream: dict,
        reference_timestamp: float,
        round_n: int = None,
        ) -> np.array:
    """Compute a stream's per-sample times relative to a shared reference timestamp
    (typically another stream's first timestamp, from `get_stream_first_timestamp`).

    Expressing two different streams' times relative to the *same* reference point (rather
    than each relative to its own start) keeps their timelines directly comparable/addable
    downstream, which is what actually makes use of the LSL clock synchronization performed
    when the XDF file is loaded.
    """
    ts = np.array(stream[STREAM_TIMESTAMPS_LABEL], dtype=np.float64)

    if len(ts) == 0:
        return ts

    rel = ts - reference_timestamp

    if round_n is not None:
        rel = np.round(rel, decimals=round_n)

    return rel


def stream_origins_filepath(times_outdir: str, subject_id: str, block) -> str:
    """Path convention for the per-subject/block JSON file recording each XDF stream's first
    synchronized timestamp for a given recording block.

    Step A (which reads the audio/eyetracker streams) and step F (which reads the EEG stream)
    are separate process invocations that each independently reload and clock-synchronize the
    same XDF file. This file is how they agree on one shared time origin for that file instead
    of each silently assuming its own stream started recording at the same instant as the others.
    """
    return os.path.join(
        times_outdir,
        str(subject_id),
        STREAM_ORIGINS_FILENAME_TEMPLATE.format(subject_id=subject_id, block=block),
    )


def write_stream_origins(filepath: str, **first_timestamps: float) -> None:
    """Persist first (LSL-synchronized) timestamps for one or more XDF streams to a JSON file."""
    with open(filepath, "w") as f:
        json.dump(first_timestamps, f, indent=2)


def read_stream_origins(filepath: str) -> dict:
    """Load first (LSL-synchronized) timestamps previously written by `write_stream_origins`."""
    with open(filepath) as f:
        return json.load(f)
