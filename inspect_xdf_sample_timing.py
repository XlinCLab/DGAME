import sys

import numpy as np
from pyxdf import load_xdf

from utils.xdf_utils import get_stream_hostname

if len(sys.argv) != 2:
    print("Usage: python inspect_xdf_sample_timing.py /path/to/file.xdf")
    sys.exit(1)

xdf_file = sys.argv[1]
streams, _header = load_xdf(xdf_file, synchronize_clocks=True, dejitter_timestamps=True, verbose=False)

for stream in streams:
    info = stream["info"]
    stype = info.get("type", [""])[0]
    if str(stype).lower() != "eeg":
        continue

    name = info.get("name", [""])[0]
    hostname = get_stream_hostname(stream)
    nominal_srate = info.get("nominal_srate", ["0"])[0]
    nominal_srate = float(nominal_srate)
    effective_srate = info.get("effective_srate", None)

    timestamps = np.asarray(stream["time_stamps"], dtype=np.float64)
    if len(timestamps) < 2:
        print(f"name={name!r} hostname={hostname!r}: fewer than 2 samples, skipping")
        continue

    intervals = np.diff(timestamps)
    nominal_period = 1.0 / nominal_srate if nominal_srate > 0 else None

    print(f"\n=== name={name!r} hostname={hostname!r} stream_id={info.get('stream_id')} ===")
    print(f"nominal_srate={nominal_srate} Hz  (period={nominal_period})")
    print(f"effective_srate (pyxdf)={effective_srate} Hz")
    print(f"n_samples={len(timestamps)}  duration={timestamps[-1] - timestamps[0]:.3f}s")
    print(f"interval stats: mean={intervals.mean():.6f}s  median={np.median(intervals):.6f}s  "
          f"std={intervals.std():.6f}s  min={intervals.min():.6f}s  max={intervals.max():.6f}s")

    if nominal_period is not None:
        # Classify each inter-sample interval relative to the nominal period:
        # roughly-nominal (device sampling on schedule) vs. an outlier gap (likely dropped sample(s))
        outlier_mask = intervals > (1.5 * nominal_period)
        n_outliers = int(outlier_mask.sum())
        print(f"intervals > 1.5x nominal period: {n_outliers} / {len(intervals)}")
        if n_outliers > 0:
            outlier_gaps = intervals[outlier_mask]
            implied_missing = np.round(outlier_gaps / nominal_period).astype(int) - 1
            print(f"  outlier gap sizes (s): min={outlier_gaps.min():.4f} max={outlier_gaps.max():.4f} "
                  f"total={outlier_gaps.sum():.3f}s")
            print(f"  implied missing samples per outlier: min={implied_missing.min()} "
                  f"max={implied_missing.max()} total={implied_missing.sum()}")
        # How much of the nominal-vs-effective rate gap is explained by these discrete outliers,
        # vs. a broad/uniform deviation across all (non-outlier) intervals
        non_outlier_intervals = intervals[~outlier_mask]
        if len(non_outlier_intervals) > 0:
            non_outlier_mean_rate = 1.0 / non_outlier_intervals.mean()
            print(f"  mean rate excluding outlier gaps: {non_outlier_mean_rate:.4f} Hz "
                  f"(vs nominal {nominal_srate} Hz)")
