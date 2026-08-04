import numpy as np
from pyxdf import load_xdf

# Raw pyxdf stream dict keys
STREAM_INFO = "info"
STREAM_TIME_SERIES = "time_series"
STREAM_TIME_STAMPS = "time_stamps"
STREAM_CLOCK_TIMES = "clock_times"
STREAM_CLOCK_VALUES = "clock_values"
INFO_NAME = "name"
INFO_TYPE = "type"
INFO_HOSTNAME = "hostname"
INFO_NOMINAL_SRATE = "nominal_srate"
INFO_DESC = "desc"
INFO_CHANNELS = "channels"
INFO_CHANNEL = "channel"
INFO_LABEL = "label"
INFO_UNIT = "unit"


class XDFStream:
    """Wraps a single stream dict as returned by pyxdf.load_xdf(), exposing its metadata,
    timing, and sample data as attributes/methods instead of raw dict indexing."""

    def __init__(self, stream: dict):
        self._stream = stream

    def get_metadata(self, field: str, result_type=None):
        """Extract a field from the stream's info metadata."""
        result = self._stream.get(STREAM_INFO, {}).get(field, [""])
        if isinstance(result, list):
            result = result[0] if result else ""
        if result_type is not None:
            return result_type(result)
        return result

    @property
    def name(self) -> str:
        return self.get_metadata(INFO_NAME, str)

    @property
    def type(self) -> str:
        return self.get_metadata(INFO_TYPE, str)

    @property
    def hostname(self) -> str:
        return self.get_metadata(INFO_HOSTNAME, str)

    @property
    def nominal_srate(self) -> float:
        return self.get_metadata(INFO_NOMINAL_SRATE, float)

    @property
    def is_regular(self) -> bool:
        """Whether this stream has a fixed nominal sampling rate, as opposed to an
        irregular/event stream (e.g. markers) sampled asynchronously."""
        return self.nominal_srate > 0

    @property
    def time_stamps(self) -> np.ndarray:
        """This stream's per-sample timestamps, exactly as loaded (LSL-synchronized
        absolute seconds if the parent XDFFile was loaded with synchronize_clocks=True,
        this stream's own raw local-clock seconds otherwise)."""
        return np.asarray(self._stream.get(STREAM_TIME_STAMPS, []), dtype=np.float64)

    @property
    def start_time(self) -> float:
        """This stream's first sample's timestamp (see `time_stamps` for which clock
        domain this is in, depending on how the parent XDFFile was loaded)."""
        timestamps = self.time_stamps
        if len(timestamps) == 0:
            raise ValueError(f"Stream {self.name!r} has no samples")
        return float(timestamps[0])

    @property
    def n_samples(self) -> int:
        return len(self.time_stamps)

    @property
    def duration(self) -> float:
        """Time span in seconds from the first to the last sample."""
        timestamps = self.time_stamps
        if len(timestamps) < 2:
            return 0.0
        return float(timestamps[-1] - timestamps[0])

    @property
    def clock_times(self) -> np.ndarray:
        """Local times (on this stream's own local clock) at which a clock-offset calibration
        measurement was taken during recording (periodic LSL clock-sync pings)."""
        return np.asarray(self._stream.get(STREAM_CLOCK_TIMES, []), dtype=np.float64)

    @property
    def clock_values(self) -> np.ndarray:
        """Measured clock offset at each of clock_times: the correction pyxdf applies
        (when loaded with synchronize_clocks=True) to map this stream's local clock onto
        the file's shared clock."""
        return np.asarray(self._stream.get(STREAM_CLOCK_VALUES, []), dtype=np.float64)

    def _raw_channel_field(self, field: str) -> list[str]:
        """
        Per-channel string field extraction from info.desc.channels.channel[]
        (e.g. "label" or "unit"), with no fallback -- used internally both 
        to disambiguate time_series orientation and to look up channel labels/units.
        Returns [] if the metadata is absent or malformed (e.g. irregular
        streams typically don't declare per-channel metadata at all).
        """
        desc = self._stream.get(STREAM_INFO, {}).get(INFO_DESC, [])
        if isinstance(desc, list) and len(desc) > 0:
            desc = desc[0]
        channels = desc.get(INFO_CHANNELS, {}) if isinstance(desc, dict) else {}
        if isinstance(channels, list) and len(channels) > 0:
            channels = channels[0]
        channels = channels.get(INFO_CHANNEL, []) if isinstance(channels, dict) else []

        values = []
        for channel in channels:
            value = channel.get(field) if isinstance(channel, dict) else None
            if isinstance(value, list):
                value = value[0] if value else ""
            values.append(str(value) if value else "")
        return values

    def _raw_channel_labels(self) -> list[str]:
        return self._raw_channel_field(INFO_LABEL)

    @property
    def time_series(self) -> np.ndarray:
        """
        Sample data as a (n_samples, n_channels) array.
        pyxdf normally already returns this shape;
        the channel-label count (when available) is used to catch the rare
        case where it comes back transposed.
        """
        samples = np.asarray(self._stream.get(STREAM_TIME_SERIES, []))
        if samples.ndim == 1:
            samples = samples[:, None]
        labels = self._raw_channel_labels()
        if len(labels) and len(labels) == samples.shape[0] and len(labels) != samples.shape[1]:
            samples = samples.T
        return samples

    @property
    def n_channels(self) -> int:
        time_series = self.time_series
        return time_series.shape[1] if time_series.ndim == 2 else 1

    @property
    def channel_labels(self) -> list[str]:
        """Per-channel labels, falling back to "Ch 0", "Ch 1", ... if the stream's
        metadata doesn't declare a usable label for every channel."""
        labels = self._raw_channel_labels()
        n_channels = self.n_channels
        if len(labels) == n_channels and all(labels):
            return labels
        return [f"Ch {i}" for i in range(n_channels)]

    @property
    def channel_units(self) -> list[str]:
        """
        Per-channel unit strings (e.g. "microvolts", "kohms").
        Returns empty strings for channels without a declared unit.
        """
        units = self._raw_channel_field(INFO_UNIT)
        n_channels = self.n_channels
        if len(units) == n_channels:
            return units
        return [""] * n_channels

    def __repr__(self):
        return f"XDFStream(name={self.name!r}, type={self.type!r}, hostname={self.hostname!r})"


class XDFFile:
    """Wraps the streams loaded from a single .xdf file via pyxdf.load_xdf()."""

    def __init__(self, path: str, **load_kwargs):
        self.path = path
        raw_streams, self.header = load_xdf(path, **load_kwargs)
        self.streams = [XDFStream(s) for s in raw_streams]

    def streams_by_type(self, stream_type: str) -> list[XDFStream]:
        """Return all streams matching a given type (case-insensitive)."""
        return [s for s in self.streams if s.type.lower() == str(stream_type).lower()]

    def stream_by_type(self, stream_type: str) -> XDFStream:
        """Return the first stream matching a given type (case-insensitive)."""
        matches = self.streams_by_type(stream_type)
        if not matches:
            raise ValueError(f"No <{stream_type}> stream found in {self.path}")
        return matches[0]

    def stream_by_name(self, name: str) -> XDFStream:
        """Return the stream with the given exact name."""
        for stream in self.streams:
            if stream.name == name:
                return stream
        raise ValueError(f"No <{name}> stream found in {self.path}")

    def __repr__(self):
        return f"XDFFile(path={self.path!r}, n_streams={len(self.streams)})"
