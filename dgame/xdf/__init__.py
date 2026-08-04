# pyxdf stream dict keys
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
