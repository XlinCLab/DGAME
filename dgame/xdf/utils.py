import numpy as np
import pandas as pd

from dgame.xdf import SYNC_REFERENCE_BLOCK_COLUMN, SYNC_REFERENCE_STREAM_COLUMN
from dgame.xdf.xdf_stream import XDFStream


def load_stream_sync_reference(sync_reference_file: str) -> pd.DataFrame:
    """Load a subject's stream sync reference CSV, indexed by (block, stream),
    so callers can look up a given stream's raw/LSL-synced start/end time
    for a given block by name."""
    return pd.read_csv(sync_reference_file, index_col=[SYNC_REFERENCE_BLOCK_COLUMN, SYNC_REFERENCE_STREAM_COLUMN])


def apply_clock_drift_correction(raw_time, drift_intercept: float, drift_slope: float):
    """Convert a raw (un-synchronized) timestamp or an array of raw timestamps
    onto the shared LSL-synchronized axis using a linear clock-drift correction
    with intercept and slope fit by XDFStream.fit_drift_correction:
    `synced = raw + drift_intercept + drift_slope * raw`

    Only valid for timestamps already on the same raw-clock domain as the stream's own
    clock_times/clock_values (e.g. a stream's own time_stamps array from a synchronize_clocks=False load).
    For regular-rate streams, pyxdf's default dejitter_timestamps=True regularizes raw
    per-sample timestamps *before* this correction is applied internally, so an XDF footer's
    first_timestamp/last_timestamp (which preserve the true, non-dejittered raw values) are NOT
    on this domain and must not be passed here.
    Use convert_relative_time_to_synced instead for times measured relative to a stream's own start
    (e.g. seconds into an exported WAV file).
    """
    return raw_time + drift_intercept + drift_slope * raw_time


def convert_relative_time_to_synced(relative_time, synced_start_time: float, drift_slope: float):
    """Convert a time measured in seconds elapsed since a stream's first sample (e.g. seconds
    into a WAV file exported from an audio stream) onto the shared LSL-synchronized axis:
    `synced = synced_start_time + relative_time * (1 + drift_slope)`

    Anchors directly on the stream's own already-synchronized start time rather than its raw start time,
    since for regular-rate streams the raw start time preserved in the XDF footer can differ
    from the raw start time pyxdf's clock-drift correction was actually fit against.
    """
    return synced_start_time + relative_time * (1 + drift_slope)


def filter_streams_by_hostname(streams: list[XDFStream], hostname: str | list[str]) -> list[XDFStream]:
    """Filter a list of streams down to those recorded on one of one or more given host
    machines. Used to disambiguate same-named streams in a DGAME3 recording, which
    multiplexes both dyad members' streams (recorded on separate rigs) into one XDF file.
    A role may have more than one hostname if its recording rig was replaced partway
    through data collection."""
    hostnames = [hostname] if isinstance(hostname, str) else hostname
    allowed = {str(h).lower() for h in hostnames}
    return [stream for stream in streams if stream.hostname.lower() in allowed]


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
