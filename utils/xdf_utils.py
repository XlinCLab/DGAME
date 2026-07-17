import numpy as np
from pyxdf import load_xdf

STREAM_TIMESTAMPS_LABEL = "time_stamps"


def fill_stream_gaps(stream: dict) -> dict:
    """Return a copy of a stream with zero-valued samples inserted between its jitter-removal
    segments (see pyxdf's "Segments and clock-segments differ" warning), so that the time
    series is contiguous and regularly sampled at the stream's nominal rate. Adds an 'is_gap'
    boolean array (same length as the filled time series) marking which samples were
    synthetically inserted, so they can be identified/filtered out downstream.

    No-op (aside from adding an all-False 'is_gap' array) if the stream has an irregular
    sample rate, has no gaps between segments, or has no samples at all.
    """
    samples = np.asarray(stream["time_series"])
    timestamps = np.asarray(stream["time_stamps"], dtype=np.float64)
    segments = stream.get("info", {}).get("segments", [])
    srate = stream.get("info", {}).get("nominal_srate", 0)
    if isinstance(srate, list):
        srate = srate[0] if srate else 0
    srate = float(srate)

    filled_stream = dict(stream)
    if srate <= 0 or len(segments) <= 1 or len(timestamps) == 0:
        filled_stream["is_gap"] = np.zeros(len(timestamps), dtype=bool)
        return filled_stream

    period = 1.0 / srate
    sample_chunks, timestamp_chunks, gap_chunks = [], [], []

    for i, (start, end) in enumerate(segments):
        sample_chunks.append(samples[start:end + 1])
        timestamp_chunks.append(timestamps[start:end + 1])
        gap_chunks.append(np.zeros(end + 1 - start, dtype=bool))

        if i < len(segments) - 1:
            next_start = segments[i + 1][0]
            gap_duration = timestamps[next_start] - timestamps[end]
            n_missing = int(round(gap_duration / period)) - 1
            if n_missing > 0:
                gap_shape = (n_missing,) + samples.shape[1:]
                sample_chunks.append(np.zeros(gap_shape, dtype=samples.dtype))
                timestamp_chunks.append(timestamps[end] + period * np.arange(1, n_missing + 1))
                gap_chunks.append(np.ones(n_missing, dtype=bool))

    filled_stream["time_series"] = np.concatenate(sample_chunks, axis=0)
    filled_stream["time_stamps"] = np.concatenate(timestamp_chunks)
    filled_stream["is_gap"] = np.concatenate(gap_chunks)
    return filled_stream


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


def get_stream_hostname(stream: dict) -> str:
    """Extract the hostname of a stream's source recording machine from its metadata."""
    hostname = stream.get("info", {}).get("hostname", [""])
    if isinstance(hostname, list):
        hostname = hostname[0] if hostname else ""
    return str(hostname)


def filter_streams_by_hostname(streams: list[dict], hostname: str | list[str]) -> list[dict]:
    """Filter a list of streams down to those recorded on one of one or more given host
    machines (a role may have had more than one hostname, e.g. if its recording rig was
    replaced partway through data collection)."""
    hostnames = [hostname] if isinstance(hostname, str) else hostname
    allowed = {str(h).lower() for h in hostnames}
    return [stream for stream in streams if get_stream_hostname(stream).lower() in allowed]


def get_xdf_streams_by_name(
        stream_label: str,
        xdf_file: str = None,
        xdf_data: list = None,
        verbose: bool = False,
        **kwargs
        ) -> list[dict]:
    """Fetches all streams matching a specific name/label from an XDF file (optionally preloaded)."""
    if xdf_data is None:
        assert xdf_file is not None, "xdf_file argument is required if no xdf_data argument is provided"
        xdf_data, _ = load_xdf(xdf_file, verbose=verbose, **kwargs)
    return [stream for stream in xdf_data if stream['info']['name'][0] == stream_label]


def get_xdf_streams_by_type(
        stream_type: str,
        xdf_file: str = None,
        xdf_data: list = None,
        exclude_name_substring: str = None,
        verbose: bool = False,
        **kwargs
        ) -> list[dict]:
    """Fetches all streams matching a specific type from an XDF file (optionally preloaded),
    optionally excluding streams whose name contains a given substring (e.g. to exclude
    impedance-check streams which otherwise share the same type as the real data stream)."""
    if xdf_data is None:
        assert xdf_file is not None, "xdf_file argument is required if no xdf_data argument is provided"
        xdf_data, _ = load_xdf(xdf_file, verbose=verbose, **kwargs)
    matches = []
    for stream in xdf_data:
        stream_type_val = stream.get('info', {}).get('type', "")
        if isinstance(stream_type_val, list):
            stream_type_val = stream_type_val[0] if stream_type_val else ""
        if str(stream_type_val).lower() != str(stream_type).lower():
            continue
        if exclude_name_substring is not None:
            stream_name = stream.get('info', {}).get('name', [""])
            stream_name = stream_name[0] if isinstance(stream_name, list) and stream_name else stream_name
            if exclude_name_substring.lower() in str(stream_name).lower():
                continue
        matches.append(stream)
    return matches


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


def extract_eeg_stream_samples(eeg_stream: dict) -> tuple[np.ndarray, float, list[str]]:
    """Extract EEG samples as (channels, samples), sampling rate, and labels."""
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
    return samples, srate, labels


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
