from test import SAMPLE_XDF_FILES, _skip_if_missing

import numpy as np
import pytest
import pyxdf

from dgame.xdf import AUDIO_STREAM, INFO_NAME, STREAM_INFO, STREAM_TIME_STAMPS
from dgame.xdf.utils import (apply_clock_drift_correction,
                             convert_relative_time_to_synced)
from dgame.xdf.xdf_stream import XDFFile


@pytest.mark.parametrize("raw_start,intercept,slope,relative_time", [
    (1000.0, 0.5, 1e-6, 10.0),
    (1_669_975_520.208033, -55.0944633, 3.299117231754625e-08, 3.00712),
    (0.0, 0.0, 0.0, 100.0),
    (500.0, -2.5, -5e-7, 0.0),
    (12345.6789, 1.23456, 2e-8, -50.0),
])
def test_convert_relative_time_to_synced_matches_apply_clock_drift_correction(raw_start, intercept, slope, relative_time):
    """convert_relative_time_to_synced(t, synced_start, slope) must exactly agree with
    apply_clock_drift_correction(raw_start + t, intercept, slope) whenever synced_start is
    itself apply_clock_drift_correction(raw_start, intercept, slope)."""
    synced_start = apply_clock_drift_correction(raw_start, intercept, slope)
    expected = apply_clock_drift_correction(raw_start + relative_time, intercept, slope)
    actual = convert_relative_time_to_synced(relative_time, synced_start, slope)
    assert actual == pytest.approx(expected, abs=1e-6)


def test_convert_relative_time_to_synced_at_zero_returns_synced_start_time_exactly():
    """The boundary case relative_time=0 (the stream's very first sample) must return
    synced_start_time exactly, with no floating-point drift introduced by the multiply."""
    synced_start_time = 1669975520.1383622
    assert convert_relative_time_to_synced(0.0, synced_start_time, drift_slope=3.3e-08) == synced_start_time


@pytest.mark.parametrize("xdf_file", SAMPLE_XDF_FILES)
def test_convert_relative_time_to_synced_matches_pyxdf_for_regular_rate_stream(xdf_file):
    """For a regular-rate stream (e.g. audio) and given only the stream's `synced_start_time`
    and `drift_slope`, `convert_relative_time_to_synced` must reproduce pyxdf's own synced
    `time_stamps` at arbitrary points across the recording, when `relative_time` is computed from
    pyxdf's own *dejittered* raw time_stamps."""
    _skip_if_missing(xdf_file)

    xdf_synced = XDFFile(xdf_file, synchronize_clocks=True, verbose=False)
    audio_stream = xdf_synced.stream_by_name(AUDIO_STREAM)
    assert audio_stream.is_regular  # dejittering only applies to regular-rate streams
    _, drift_slope = audio_stream.fit_drift_correction()
    synced_start_time = audio_stream.start_time

    # Independent oracle: pyxdf's own dejittered raw time_stamps and its own synced time_stamps
    raw_streams, _ = pyxdf.load_xdf(xdf_file, synchronize_clocks=False, verbose=False)
    raw_stream = next(s for s in raw_streams if s[STREAM_INFO][INFO_NAME][0] == AUDIO_STREAM)
    raw_time_stamps = np.asarray(raw_stream[STREAM_TIME_STAMPS], dtype=np.float64)
    actual_synced_time_stamps = audio_stream.time_stamps

    sample_indices = np.linspace(0, len(raw_time_stamps) - 1, num=5, dtype=int)
    for i in sample_indices:
        relative_time = raw_time_stamps[i] - raw_time_stamps[0]
        predicted = convert_relative_time_to_synced(relative_time, synced_start_time, drift_slope)
        assert predicted == pytest.approx(actual_synced_time_stamps[i], abs=1e-5)


@pytest.mark.parametrize("xdf_file", SAMPLE_XDF_FILES)
def test_convert_relative_time_to_synced_corrects_naive_footer_anchored_divergence(xdf_file):
    """For a real word onset (audio-relative seconds into an exported WAV file), naively anchoring on
    the footer's raw_start_time and running the full apply_clock_drift_correction formula
    diverges measurably from pyxdf's own synced value, while convert_relative_time_to_synced
    (anchored on the trusted synced_start_time) stays within microseconds of it."""
    _skip_if_missing(xdf_file)

    xdf_synced = XDFFile(xdf_file, synchronize_clocks=True, verbose=False)
    stream = xdf_synced.stream_by_name(AUDIO_STREAM)
    intercept, slope = stream.fit_drift_correction()
    synced_start_time = stream.start_time
    footer_raw_start_time = stream.raw_start_time

    relative_time = 3.00712  # an arbitrary WAV-relative word onset, seconds into the audio file
    true_synced_time = synced_start_time + relative_time * (1 + slope)

    raw_footer_synced_result = apply_clock_drift_correction(footer_raw_start_time + relative_time, intercept, slope)
    relative_synced_result = convert_relative_time_to_synced(relative_time, synced_start_time, slope)

    # Verify that the approach using raw start timestamp from XDF footer is off from the true synced time by a measurable amount
    assert abs(raw_footer_synced_result - true_synced_time) > 1e-3
    # Verify that relative timestamp conversion approach reproduces the true synced time to floating-point precision
    assert relative_synced_result == pytest.approx(true_synced_time, abs=1e-9)
