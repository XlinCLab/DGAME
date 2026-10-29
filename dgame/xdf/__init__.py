import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(name)s %(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

# Default DGAME XDF stream labels
AUDIO_STREAM = "audio"
EYETRACKER_STREAM = "pupil_capture"

# Stream sync-reference file columns
SYNC_REFERENCE_BLOCK_COLUMN = "block"
SYNC_REFERENCE_STREAM_COLUMN = "stream"
SYNC_REFERENCE_ROLE_COLUMN = "role"  # participant role (e.g. Director/Matcher) that a stream belongs to, if relevant
# Raw (un-synchronized, stream's local clock) start/end times
SYNC_REFERENCE_RAW_START_TIME_COLUMN = "raw_start_time"
SYNC_REFERENCE_RAW_END_TIME_COLUMN = "raw_end_time"
# LSL-synchronized (shared, cross-stream-comparable) start/end times
SYNC_REFERENCE_SYNCED_START_TIME_COLUMN = "synced_start_time"
SYNC_REFERENCE_SYNCED_END_TIME_COLUMN = "synced_end_time"
# Linear clock-drift correction: synced = raw + drift_intercept + drift_slope * raw
SYNC_REFERENCE_DRIFT_INTERCEPT_COLUMN = "drift_intercept"
SYNC_REFERENCE_DRIFT_SLOPE_COLUMN = "drift_slope"
