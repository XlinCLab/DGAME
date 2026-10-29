import pandas as pd

from dgame.xdf import SYNC_REFERENCE_BLOCK_COLUMN, SYNC_REFERENCE_STREAM_COLUMN


def load_stream_sync_reference(sync_reference_file: str) -> pd.DataFrame:
    """Load a subject's stream sync reference CSV, indexed by (block, stream),
    so callers can look up a given stream's raw/LSL-synced start/end time
    for a given block by name."""
    return pd.read_csv(sync_reference_file, index_col=[SYNC_REFERENCE_BLOCK_COLUMN, SYNC_REFERENCE_STREAM_COLUMN])
