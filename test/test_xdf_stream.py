from test import SAMPLE_XDF_FILES, _skip_if_missing

import numpy as np
import pytest
import pyxdf

from dgame.xdf import (_EXPECTED_WINSOR_THRESHOLD, AUDIO_STREAM,
                       EYETRACKER_STREAM, INFO_NAME, STREAM_INFO,
                       STREAM_TIME_STAMPS, WINSOR_THRESHOLD)
from dgame.xdf.utils import apply_clock_drift_correction
from dgame.xdf.xdf_stream import XDFFile


def test_winsor_default_unchanged():
    """Verify that pyxdf's `winsor_default` (extracted from pyxdf.load_xdf() parameter defaults)
    has not changed internally; must match _EXPECTED_WINSOR_THRESHOLD constant."""
    assert WINSOR_THRESHOLD == _EXPECTED_WINSOR_THRESHOLD


@pytest.mark.parametrize("xdf_file", SAMPLE_XDF_FILES)
def test_start_time_matches_raw_pyxdf(xdf_file):
    """XDFStream.start_time (via XDFFile) must match pyxdf's own synchronized time_stamps[0] exactly."""
    _skip_if_missing(xdf_file)

    # Independent oracle: call pyxdf directly, bypassing XDFStream/XDFFile entirely
    oracle_streams, _ = pyxdf.load_xdf(xdf_file, synchronize_clocks=True, verbose=False)
    oracle_start_times = {}
    for stream in oracle_streams:
        name = stream[STREAM_INFO][INFO_NAME][0]
        oracle_start_times[name] = float(stream[STREAM_TIME_STAMPS][0])

    xdf = XDFFile(xdf_file, synchronize_clocks=True, verbose=False)
    audio_stream = xdf.stream_by_name(AUDIO_STREAM)
    eyetracker_stream = xdf.stream_by_name(EYETRACKER_STREAM)
    eeg_stream = xdf.stream_by_type("EEG")

    assert audio_stream.start_time == oracle_start_times[AUDIO_STREAM]
    assert eyetracker_stream.start_time == oracle_start_times[EYETRACKER_STREAM]
    assert eeg_stream.start_time == oracle_start_times[eeg_stream.name]


@pytest.mark.parametrize("xdf_file", SAMPLE_XDF_FILES)
def test_inter_host_stream_offset_is_non_trivial(xdf_file):
    """Regression guard: audio, eyetracker, and EEG streams in `SAMPLE_XDF_FILES` come
    from two different host computers, so their LSL-synchronized start times must
    differ by a measurable amount. If this ever collapses to ~0, it means
    clock synchronization silently stopped being applied (e.g. synchronize_clocks
    reverted to False, or streams got re-zeroed independently)."""
    _skip_if_missing(xdf_file)

    xdf = XDFFile(xdf_file, synchronize_clocks=True, verbose=False)
    audio_stream = xdf.stream_by_name(AUDIO_STREAM)
    audio_start = audio_stream.start_time
    eyetracker_stream = xdf.stream_by_name(EYETRACKER_STREAM)
    eyetracker_start = eyetracker_stream.start_time
    eeg_stream = xdf.stream_by_type("EEG")
    eeg_start = eeg_stream.start_time

    assert audio_stream.hostname != eeg_stream.hostname
    assert audio_stream.hostname != eyetracker_stream.hostname
    assert abs(audio_start - eeg_start) > 1e-3
    assert abs(audio_start - eyetracker_start) > 1e-3
    assert abs(eyetracker_start - eeg_start) > 1e-3


@pytest.mark.parametrize("xdf_file", SAMPLE_XDF_FILES)
def test_raw_and_synced_clocks_are_in_different_domains(xdf_file):
    """Audio/webcam host and the eyetracker/EEG host are on unrelated raw clock
    epochs; only synchronize_clocks=True should bring them into a common frame.
    This guards against ever comparing raw timestamps across hosts by mistake."""
    _skip_if_missing(xdf_file)

    xdf_raw = XDFFile(xdf_file, synchronize_clocks=False, verbose=False)
    raw_audio_stream = xdf_raw.stream_by_name(AUDIO_STREAM)
    audio_raw_start = raw_audio_stream.start_time
    raw_eyetracker_stream = xdf_raw.stream_by_name(EYETRACKER_STREAM)
    eyetracker_raw_start = raw_eyetracker_stream.start_time
    assert raw_audio_stream.hostname != raw_eyetracker_stream.hostname
    # Raw clocks from the two hosts differ by many orders of magnitude more than
    # any real inter-stream recording-start offset could ever be
    assert abs(audio_raw_start - eyetracker_raw_start) > 1000000000

    # Check that reloading with synchronize_clocks=True brings them into alignment
    xdf_synced = XDFFile(xdf_file, synchronize_clocks=True, verbose=False)
    synced_audio_stream = xdf_synced.stream_by_name(AUDIO_STREAM)
    audio_synced_start = synced_audio_stream.start_time
    synced_eyetracker_stream = xdf_synced.stream_by_name(EYETRACKER_STREAM)
    eyetracker_synced_start = synced_eyetracker_stream.start_time
    assert abs(audio_synced_start - eyetracker_synced_start) < 1
    assert abs(eyetracker_synced_start - eyetracker_raw_start) > 1000000000


@pytest.mark.parametrize("xdf_file", SAMPLE_XDF_FILES)
@pytest.mark.parametrize("stream_name", [AUDIO_STREAM, EYETRACKER_STREAM])
def test_fit_drift_correction_matches_pyxdf_own_synced_time_stamps(xdf_file, stream_name):
    """XDFStream.fit_drift_correction() must reproduce pyxdf's own applied clock
    correction almost exactly (to floating-point noise), since it reuses pyxdf's own
    _robust_fit on the same clock_times/clock_values."""
    _skip_if_missing(xdf_file)

    xdf_synced = XDFFile(xdf_file, synchronize_clocks=True, verbose=False)
    stream = xdf_synced.stream_by_name(stream_name)
    intercept, slope = stream.fit_drift_correction()

    # Independent oracle: pyxdf's own raw (un-synchronized) time_stamps for this stream
    raw_streams, _ = pyxdf.load_xdf(xdf_file, synchronize_clocks=False, verbose=False)
    raw_stream = next(s for s in raw_streams if s[STREAM_INFO][INFO_NAME][0] == stream_name)
    raw_time_stamps = np.asarray(raw_stream[STREAM_TIME_STAMPS], dtype=np.float64)

    # predicted_synced = raw_time_stamps + intercept + slope * raw_time_stamps
    predicted_synced = apply_clock_drift_correction(
        raw_time=raw_time_stamps,
        drift_intercept=intercept,
        drift_slope=slope,
    )
    actual_synced = stream.time_stamps

    assert predicted_synced.shape == actual_synced.shape
    # Sub-microsecond agreement with pyxdf's own already-applied correction
    assert np.max(np.abs(predicted_synced - actual_synced)) < 1e-5


@pytest.mark.parametrize("xdf_file", SAMPLE_XDF_FILES)
@pytest.mark.parametrize("stream_name", [AUDIO_STREAM, EYETRACKER_STREAM])
def test_clock_segments_is_single_segment(xdf_file, stream_name):
    """Drift correction assumes exactly one clock segment (no mid-recording clock reset).
    This should hold for all current sample recordings; if it ever doesn't,
    XDFStream.fit_drift_correction() must raise ClockResetError rather than silently
    fitting one line across a discontinuity."""
    _skip_if_missing(xdf_file)

    xdf_synced = XDFFile(xdf_file, synchronize_clocks=True, verbose=False)
    stream = xdf_synced.stream_by_name(stream_name)
    assert len(stream.clock_segments) == 1
