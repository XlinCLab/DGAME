import numpy as np
from pyxdf import load_xdf
from pyxdf.pyxdf import _robust_fit

from dgame.xdf import (FOOTER_INFO_FIRST_TIMESTAMP, FOOTER_INFO_LAST_TIMESTAMP,
                       INFO_CHANNEL, INFO_CHANNELS, INFO_CLOCK_SEGMENTS,
                       INFO_DESC, INFO_HOSTNAME, INFO_LABEL, INFO_NAME,
                       INFO_NOMINAL_SRATE, INFO_TYPE, INFO_UNIT,
                       STREAM_CLOCK_TIMES, STREAM_CLOCK_VALUES, STREAM_FOOTER,
                       STREAM_INFO, STREAM_TIME_SERIES, STREAM_TIME_STAMPS,
                       WINSOR_THRESHOLD)
from experiment.input_validation import InputValidationError


class ClockResetError(InputValidationError):
    """Raised when a stream's clock-offset calibration data shows more than one
    segment, i.e. pyxdf detected a mid-recording clock discontinuity
    (e.g. if the recording device was restarted during recording, which is not expected).
    Fitting a single linear drift correction across a reset would be wrong,
    so this is treated as an error rather than handled."""


class XDFStream:
    """Wraps a single stream dict as returned by pyxdf.load_xdf(), exposing its metadata,
    timing, and sample data as attributes/methods instead of raw dict indexing."""

    def __init__(self, stream: dict, clocks_synced: bool = True):
        self._stream = stream
        self._clocks_synced = clocks_synced

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
    def start_time(self) -> float:
        """This stream's first sample's timestamp (see `time_stamps` for which clock
        domain this is in, depending on how the parent XDFFile was loaded)."""
        timestamps = self.time_stamps
        if len(timestamps) == 0:
            raise ValueError(f"Stream {self.name!r} has no samples")
        return float(timestamps[0])

    @property
    def end_time(self) -> float:
        """This stream's last sample's timestamp (see `time_stamps` for which clock
        domain this is in, depending on how the parent XDFFile was loaded)."""
        timestamps = self.time_stamps
        if len(timestamps) == 0:
            raise ValueError(f"Stream {self.name!r} has no samples")
        return float(timestamps[-1])

    def get_footer_info(self, field: str) -> float:
        """Read field from the XDF footer chunk.
        The footer is parsed unconditionally and is never touched by clock synchronization,
        so timestamps listed here are always the stream's raw (un-synchronized) boundary timestamps,
        regardless of whether the parent XDFFile was loaded with synchronize_clocks=True."""
        footer_info = self._stream.get(STREAM_FOOTER, {}).get(STREAM_INFO, {})
        value = footer_info.get(field)
        if not value:
            raise ValueError(f"Stream {self.name!r} has no {field!r} in its XDF footer chunk")
        if isinstance(value, list):
            value = value[0]
        return float(value)

    @property
    def raw_start_time(self) -> float:
        """This stream's raw (un-synchronized, own local clock) first-sample timestamp,
        read from the XDF footer chunk."""
        return self.get_footer_info(FOOTER_INFO_FIRST_TIMESTAMP)

    @property
    def raw_end_time(self) -> float:
        """This stream's raw (un-synchronized, own local clock) last-sample timestamp,
        read from the XDF footer chunk."""
        return self.get_footer_info(FOOTER_INFO_LAST_TIMESTAMP)

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

    @property
    def clock_segments(self) -> list[tuple[int, int]]:
        """Segment boundaries (start, end) index pairs of clock_times/clock_values,
        identified by pyxdf's clock-reset detection. 
        A single segment means no clock discontinuity was detected during the recording.
        Only meaningful when the parent XDFFile was loaded with synchronize_clocks=True
        (empty otherwise, since reset detection only runs as part of that step)."""
        segments = self._stream.get(STREAM_INFO, {}).get(INFO_CLOCK_SEGMENTS, [])
        return [tuple(segment) for segment in segments]

    def fit_drift_correction(self) -> tuple[float, float]:
        """Fit this stream's linear clock-offset correction with pyxdf's Huber-robust 
        regression (pyxdf.pyxdf._robust_fit) applied to clock_times/clock_values,
        reproducing the exact preprocessing pyxdf.pyxdf._clock_sync does,
        in order to make the drift correction intercept and slope available
        (otherwise discarded and not exposed in pyxdf).

        Raises ClockResetError if this stream's clock-offset data has more than one
        clock_segments entry (a detected mid-recording clock reset), since a single
        linear fit across a reset would be wrong. Also requires the XDFFile to have
        been loaded with synchronize_clocks=True (so that clock_segments is populated)."""
        if not self._clocks_synced:
            raise ValueError(
                f"Stream {self.name!r} was loaded with synchronize_clocks=False; "
                "clock_segments is only populated when synchronize_clocks=True."
            )
        segments = self.clock_segments
        if len(segments) != 1:
            raise ClockResetError(
                f"Stream {self.name!r} has {len(segments)} clock segments (expected "
                "exactly 1), meaning that a mid-recording clock reset was detected, "
                "which drift correction does not currently handle."
            )
        clock_times = self.clock_times
        clock_values = self.clock_values
        if len(clock_times) < 2:
            raise ValueError(
                f"Stream {self.name!r} has fewer than 2 clock-offset measurements; "
                "cannot fit a drift correction."
            )
        # Reproduces pyxdf.pyxdf._clock_sync's own preprocessing exactly
        design_matrix = np.column_stack([
            np.ones(len(clock_times)),
            clock_times / WINSOR_THRESHOLD,
        ])
        target = clock_values / WINSOR_THRESHOLD
        intercept, slope = _robust_fit(design_matrix, target)
        intercept *= WINSOR_THRESHOLD
        return float(intercept), float(slope)

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
    def relative_times(self) -> np.array:
        """Return XDF stream's timestamps relative to its start time (first timestamp)."""
        timestamps = self.time_stamps
        start_time = self.start_time
        relative_timestamps = timestamps - start_time
        return relative_timestamps

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
        self._clocks_synced = load_kwargs.get("synchronize_clocks", True)  # NB: pyxdf's default for `synchronize_clocks` is True
        self.streams = [XDFStream(s, clocks_synced=self._clocks_synced) for s in raw_streams]

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
