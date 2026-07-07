import argparse
import os

import numpy as np
import yaml
from pyxdf import load_xdf

from dgame.constants import BLOCK_IDS, IMPEDANCE_STREAM_NAME_SUBSTRING
from utils.xdf_utils import (filter_streams_by_hostname, get_stream_hostname,
                             get_xdf_streams_by_type)


def load_rig_hostnames(config_path: str) -> dict:
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)
    return config["experiment"]["rig_hostnames"]


def summarize_stream_gaps(stream: dict) -> dict:
    """Summarize sample-timing quality for a single EEG stream."""
    timestamps = np.asarray(stream["time_stamps"], dtype=np.float64)
    nominal_srate = stream.get("info", {}).get("nominal_srate", ["0"])
    nominal_srate = float(nominal_srate[0] if isinstance(nominal_srate, list) else nominal_srate)
    n_samples = len(timestamps)

    if n_samples < 2 or nominal_srate <= 0:
        return {
            "n_samples": n_samples,
            "duration": 0.0,
            "effective_srate": 0.0,
            "n_gaps": 0,
            "total_missing": 0,
            "max_gap_missing": 0,
        }

    duration = timestamps[-1] - timestamps[0]
    effective_srate = (n_samples - 1) / duration if duration > 0 else 0.0
    nominal_period = 1.0 / nominal_srate

    intervals = np.diff(timestamps)
    outlier_mask = intervals > (1.5 * nominal_period)
    outlier_gaps = intervals[outlier_mask]
    implied_missing = np.round(outlier_gaps / nominal_period).astype(int) - 1
    implied_missing = implied_missing[implied_missing > 0]

    return {
        "n_samples": n_samples,
        "duration": duration,
        "effective_srate": effective_srate,
        "nominal_srate": nominal_srate,
        "n_gaps": len(implied_missing),
        "total_missing": int(implied_missing.sum()) if len(implied_missing) else 0,
        "max_gap_missing": int(implied_missing.max()) if len(implied_missing) else 0,
    }


def classify(summary: dict, expected_n_samples: int = None) -> str:
    if summary["n_samples"] < 2:
        return "NO DATA"
    if expected_n_samples is not None and summary["n_samples"] < 0.5 * expected_n_samples:
        return "CORRUPT / TOO SHORT"
    # Sustained rate mismatch with no discrete gaps: every interval is a bit off nominal,
    # rather than a few isolated dropouts -- classify() would otherwise call this "ok"
    if summary["n_gaps"] == 0:
        rate_ratio = summary["effective_srate"] / summary["nominal_srate"]
        if abs(1 - rate_ratio) > 0.02:
            return "SUSTAINED RATE MISMATCH (no discrete gaps, likely real amplifier clock deviation)"
        return "ok"
    missing_fraction = summary["total_missing"] / summary["n_samples"]
    if missing_fraction < 0.005:
        return "ok (negligible gaps)"
    if summary["n_gaps"] <= 3 and missing_fraction < 0.15:
        return "salvageable (few concentrated gaps)"
    return "SEVERE (many/large gaps)"


def main(xdf_root: str, config_path: str, dgame_version: str = "3") -> None:
    rig_hostnames = load_rig_hostnames(config_path)
    roles = list(rig_hostnames.keys())

    dyad_ids = sorted(
        d for d in os.listdir(xdf_root)
        if os.path.isdir(os.path.join(xdf_root, d))
    )

    header = f"{'dyad':>6} {'block':>6} {'role':>10} {'n_samples':>10} {'dur(s)':>8} {'eff_Hz':>8} {'n_gaps':>7} {'missing':>8}  verdict"
    print(header)
    print("-" * len(header))

    for dyad_id in dyad_ids:
        for block in sorted(BLOCK_IDS):
            xdf_file = os.path.join(
                xdf_root, dyad_id, "Director", f"dgame{dgame_version}_{dyad_id}_Director_{block}.xdf"
            )
            if not os.path.exists(xdf_file):
                continue

            try:
                streams, _header = load_xdf(xdf_file, synchronize_clocks=True, verbose=False)
            except Exception as exc:
                print(f"{dyad_id:>6} {block:>6}  FAILED TO LOAD: {exc}")
                continue

            eeg_streams = get_xdf_streams_by_type(
                "EEG", xdf_data=streams, exclude_name_substring=IMPEDANCE_STREAM_NAME_SUBSTRING
            )

            role_summaries = {}
            for role in roles:
                matches = filter_streams_by_hostname(eeg_streams, rig_hostnames[role])
                if len(matches) != 1:
                    hosts_found = [get_stream_hostname(s) for s in eeg_streams]
                    print(f"{dyad_id:>6} {block:>6} {role:>10}  UNEXPECTED STREAM COUNT ({len(matches)}); "
                          f"hosts in file: {hosts_found}")
                    continue
                role_summaries[role] = summarize_stream_gaps(matches[0])

            if not role_summaries:
                continue

            expected_n_samples = max(s["n_samples"] for s in role_summaries.values())
            for role, summary in role_summaries.items():
                verdict = classify(summary, expected_n_samples=expected_n_samples)
                print(
                    f"{dyad_id:>6} {block:>6} {role:>10} {summary['n_samples']:>10} "
                    f"{summary['duration']:>8.1f} {summary['effective_srate']:>8.1f} "
                    f"{summary['n_gaps']:>7} {summary['total_missing']:>8}  {verdict}"
                )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Batch-check EEG sample-timing quality (gaps, effective rate) across all dyads/blocks."
    )
    parser.add_argument("xdf_root", help="Path to recordings/xdf directory (contains one subdir per dyad)")
    parser.add_argument(
        "--config", default="config/dgame3_defaults.yml",
        help="Path to a config yml containing experiment.rig_hostnames (default: config/dgame3_defaults.yml)"
    )
    parser.add_argument("--dgame-version", default="3")
    args = parser.parse_args()
    main(args.xdf_root, args.config, args.dgame_version)
