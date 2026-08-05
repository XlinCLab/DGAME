import inspect

from pyxdf import load_xdf

# pyxdf stream dict keys
STREAM_INFO = "info"
STREAM_TIME_SERIES = "time_series"
STREAM_TIME_STAMPS = "time_stamps"
STREAM_CLOCK_TIMES = "clock_times"
STREAM_CLOCK_VALUES = "clock_values"
STREAM_FOOTER = "footer"
FOOTER_INFO_FIRST_TIMESTAMP = "first_timestamp"
FOOTER_INFO_LAST_TIMESTAMP = "last_timestamp"
INFO_NAME = "name"
INFO_TYPE = "type"
INFO_HOSTNAME = "hostname"
INFO_NOMINAL_SRATE = "nominal_srate"
INFO_DESC = "desc"
INFO_CHANNELS = "channels"
INFO_CHANNEL = "channel"
INFO_LABEL = "label"
INFO_UNIT = "unit"
# Segment boundaries (index pairs into clock_times/clock_values) that pyxdf's own
# clock-reset detection identified; only populated when loaded with synchronize_clocks=True
INFO_CLOCK_SEGMENTS = "clock_segments"
# pyxdf scales clock_times/clock_values by winsor_threshold before fitting drift correction
# set value directly from load_xdf's own default so this can never drift out of sync with pyxdf
WINSOR_THRESHOLD = inspect.signature(load_xdf).parameters["winsor_threshold"].default

# Default DGAME XDF stream labels
AUDIO_STREAM = "audio"
EYETRACKER_STREAM = "pupil_capture"

# Stream sync-reference file columns
SYNC_REFERENCE_BLOCK_COLUMN = "block"
SYNC_REFERENCE_STREAM_COLUMN = "stream"
# Raw (un-synchronized, stream's local clock) start/end times
SYNC_REFERENCE_RAW_START_TIME_COLUMN = "raw_start_time"
SYNC_REFERENCE_RAW_END_TIME_COLUMN = "raw_end_time"
# LSL-synchronized (shared, cross-stream-comparable) start/end times
SYNC_REFERENCE_SYNCED_START_TIME_COLUMN = "synced_start_time"
SYNC_REFERENCE_SYNCED_END_TIME_COLUMN = "synced_end_time"
# Linear clock-drift correction: synced = raw + drift_intercept + drift_slope * raw
SYNC_REFERENCE_DRIFT_INTERCEPT_COLUMN = "drift_intercept"
SYNC_REFERENCE_DRIFT_SLOPE_COLUMN = "drift_slope"
