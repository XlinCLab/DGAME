"""
gaze_subsampling_sensitivity.py
================================================================================
Subsampling sensitivity analysis on eyetrackingR-aggregated data.

Input: one CSV per analysis, each the output of eyetrackingR's
make_time_cluster_data() — i.e., the per-(subj, condition, TimeBin) Prop table
that the cluster permutation test runs on. Schema:

    subj, condition, TimeBin, SamplesInAOI, SamplesTotal, AOI,
    Elog, Weights, Prop, LogitAdjusted, ArcSin, Time, ot1..ot7,
    ClusterNeg, ClusterPos

By taking these as input we avoid re-implementing eyetrackingR's preprocessing
and exactly match what Da_gaze_stats.r feeds into analyze_time_clusters.

Mirrors analyze_time_clusters(..., within_subj = TRUE, paired = TRUE):
  1. Threshold: qt(0.975, df = N - 1) from the (sub)sample size.
  2. Per-bin one-sample paired t-test on Prop[level_a] - Prop[level_b],
     dropping NaN per bin (na.rm).
  3. Sign-stable cluster runs (ClusterPos vs ClusterNeg kept separate).
  4. Sign-flip permutation null distribution of signed extreme cluster sums.
  5. Two-tailed p-value: p_left + p_right against the signed null.

The three analyses match Da_gaze_stats.r:
  - target:     condition='conflict' minus condition='no_conflict' on aoi_target
  - goal:       condition='no_conflict' minus condition='conflict' on aoi_goal
                (matches treatment_level='no_conflict' in R)
  - competitor: condition='conflict' only, predictor was relabeled aoi_fct,
                so the input CSV for this analysis must use aoi_fct as the
                grouping column (values 'aoi_comp' vs 'aoi_AllOther').
                Pass --pred_col_competitor aoi_fct.

Usage (full sample only):
    python gaze_subsampling_sensitivity.py \\
        --target_csv  time_cluster_data_target.csv \\
        --goal_csv    time_cluster_data_goal.csv \\
        --comp_csv    time_cluster_data_comp.csv \\
        --n_min 26 --n_max 26 --n_draws 1 \\
        --output_dir gaze_sensitivity_results

Usage (full sensitivity sweep):
    python gaze_subsampling_sensitivity.py \\
        --target_csv  time_cluster_data_target.csv \\
        --goal_csv    time_cluster_data_goal.csv \\
        --comp_csv    time_cluster_data_comp.csv \\
        --n_min 15 --n_max 26 --n_draws 10 \\
        --n_permutations 4000 \\
        --output_dir gaze_sensitivity_results
================================================================================
"""

import argparse
import itertools
import os
from math import comb

import numpy as np
import pandas as pd
from scipy import stats
from tqdm import tqdm

# =============================================================================
# Constants
# =============================================================================

ALPHA = 0.05
N_PERMUTATIONS = 4000


# =============================================================================
# Cluster permutation test (matches eyetrackingR analyze_time_clusters)
# =============================================================================

def find_clusters(tstats, threshold):
    """Sign-stable contiguous runs above ±threshold (eyetrackingR ClusterPos/Neg).

    Returns list of dicts with start_idx, end_idx, cluster_sum, n_bins, direction.
    """
    clusters = []
    pos = tstats >= threshold
    i = 0
    while i < len(pos):
        if pos[i]:
            j = i + 1
            while j < len(pos) and pos[j]:
                j += 1
            clusters.append({"start_idx": i, "end_idx": j - 1,
                             "cluster_sum": float(np.sum(tstats[i:j])),
                             "n_bins": j - i, "direction": "Positive"})
            i = j
        else:
            i += 1
    neg = tstats <= -threshold
    i = 0
    while i < len(neg):
        if neg[i]:
            j = i + 1
            while j < len(neg) and neg[j]:
                j += 1
            clusters.append({"start_idx": i, "end_idx": j - 1,
                             "cluster_sum": float(np.sum(tstats[i:j])),
                             "n_bins": j - i, "direction": "Negative"})
            i = j
        else:
            i += 1
    return clusters


def _max_cluster_signed_sum_batch(perm_tstats, threshold):
    """Signed cluster sum with the largest |sum| per permutation row.

    Mirrors eyetrackingR (within-subj branch): ss[which.max(abs(ss))].
    """
    n_perm = perm_tstats.shape[0]
    out = np.zeros(n_perm)
    for pi in range(n_perm):
        row = perm_tstats[pi]
        best_abs, best_signed = 0.0, 0.0
        for sign_mask in (row >= threshold, row <= -threshold):
            if sign_mask.any():
                padded = np.concatenate([[False], sign_mask, [False]])
                d = np.diff(padded.astype(np.int8))
                starts = np.where(d == 1)[0]
                ends = np.where(d == -1)[0]
                for s, e in zip(starts, ends):
                    cs = float(np.sum(row[s:e]))
                    if abs(cs) > best_abs:
                        best_abs = abs(cs)
                        best_signed = cs
        out[pi] = best_signed
    return out


def _ttest_per_bin_nan(diff_matrix):
    """One-sample paired t-test per bin with na.rm.

    Returns (tstats,) of shape (n_bins,). Bins with <2 valid subjects: t=0.
    """
    n_sub, n_bins = diff_matrix.shape
    tstats = np.zeros(n_bins)
    for b in range(n_bins):
        col = diff_matrix[:, b]
        v = col[~np.isnan(col)]
        if len(v) < 2:
            continue
        sd = v.std(ddof=1)
        if sd > 0:
            tstats[b] = v.mean() / (sd / np.sqrt(len(v)))
    return tstats


def cluster_permutation_test(diff_matrix, time_ms, n_permutations, rng):
    """Paired cluster-based permutation test.

    diff_matrix: (n_sub, n_bins), Prop[a] - Prop[b] per subject/bin (NaN OK).
    time_ms:     (n_bins,) bin times in ms (eyetrackingR's Time column).
    Threshold:   qt(0.975, df = n_sub - 1) from full subsample N.
    p-value:     two-tailed against signed null of extreme cluster sums.
    """
    n_sub, _ = diff_matrix.shape
    if n_sub < 2:
        return _empty_cluster_result()

    threshold_t = float(stats.t.ppf(1 - ALPHA / 2, df=n_sub - 1))

    obs_tstats = _ttest_per_bin_nan(diff_matrix)
    obs_clusters = find_clusters(obs_tstats, threshold_t)

    if not obs_clusters:
        return {**_empty_cluster_result(),
                "p_value": 1.0, "max_cluster_sum": 0.0,
                "max_abs_cluster_sum": 0.0, "threshold_t": threshold_t}

    max_cluster = max(obs_clusters, key=lambda c: abs(c["cluster_sum"]))
    obs_signed = float(max_cluster["cluster_sum"])
    obs_max = abs(obs_signed)

    # Sign-flip permutation, NaN-aware
    signs = rng.integers(0, 2, size=(n_permutations, n_sub)) * 2 - 1
    valid_mask = ~np.isnan(diff_matrix)
    diff_filled = np.where(valid_mask, diff_matrix, 0.0)

    perm_diff = diff_filled[np.newaxis, :, :] * signs[:, :, np.newaxis]
    n_valid = valid_mask.sum(axis=0).astype(float)
    sum_x = perm_diff.sum(axis=1)
    sum_x2 = (perm_diff ** 2).sum(axis=1)

    with np.errstate(invalid="ignore", divide="ignore"):
        means = sum_x / n_valid[np.newaxis, :]
        var = (sum_x2 - n_valid[np.newaxis, :] * means ** 2) \
              / np.maximum(n_valid - 1, 1)
        var = np.where(var < 0, 0.0, var)
        sds = np.sqrt(var)
        se = sds / np.sqrt(np.maximum(n_valid, 1))[np.newaxis, :]
        perm_tstats = np.where(se > 0, means / se, 0.0)

    insufficient = n_valid < 2
    if insufficient.any():
        perm_tstats[:, insufficient] = 0.0

    perm_signed = _max_cluster_signed_sum_batch(perm_tstats, threshold_t)
    p_left = float(np.mean(perm_signed <= -obs_max))
    p_right = float(np.mean(perm_signed >= obs_max))
    p_value = p_left + p_right

    return {
        "n_clusters": len(obs_clusters),
        "max_cluster_sum": obs_signed,
        "max_abs_cluster_sum": obs_max,
        "p_value": p_value,
        "significant": p_value < ALPHA,
        "cluster_start_ms": float(time_ms[max_cluster["start_idx"]]),
        "cluster_end_ms": float(time_ms[max_cluster["end_idx"]]),
        "cluster_n_bins": max_cluster["n_bins"],
        "cluster_direction": max_cluster["direction"],
        "threshold_t": threshold_t,
    }


def _empty_cluster_result():
    return {
        "n_clusters": 0,
        "max_cluster_sum": np.nan,
        "max_abs_cluster_sum": np.nan,
        "p_value": np.nan,
        "significant": False,
        "cluster_start_ms": np.nan,
        "cluster_end_ms": np.nan,
        "cluster_n_bins": 0,
        "cluster_direction": "",
        "threshold_t": np.nan,
    }


# =============================================================================
# Build a diff matrix from an eyetrackingR-aggregated CSV
# =============================================================================

def _load_time_cluster_csv(path):
    """Load an eyetrackingR time_cluster_data CSV.

    Tolerates the unnamed first column (R rownames) that read.csv produces.
    Returns (df, time_bin_size_ms): the bin size is recovered from
    Time / TimeBin (eyetrackingR multiplies by time_bin_size, so the ratio
    is the bin size in original units).
    """
    df = pd.read_csv(path, index_col=0)
    required = {"subj", "condition", "TimeBin", "Prop", "Time"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path}: missing required columns: {missing}. "
                         f"Got: {list(df.columns)}")
    # Recover bin size in ms (Time = TimeBin * time_bin_size in eyetrackingR)
    nz = df[df["TimeBin"] != 0]
    if len(nz) > 0:
        bin_size = int(round(float(nz["Time"].iloc[0])
                              / float(nz["TimeBin"].iloc[0])))
    else:
        bin_size = 100
    return df, bin_size


def build_diff_matrix(df, subjects, pred_col, level_a, level_b):
    """Pivot pred_col, compute Prop[level_a] - Prop[level_b] per (subj, TimeBin).

    NaN cells (subject without data for one level in a bin) are preserved —
    eyetrackingR's per-bin t.test drops NA per bin (na.rm).
    Returns (diff_matrix, time_ms) or (None, None) if data is insufficient.
    """
    sub_df = df[df["subj"].isin(subjects)]
    if pred_col not in sub_df.columns:
        return None, None
    try:
        pivot = sub_df.pivot_table(
            index=["subj", "TimeBin"], columns=pred_col,
            values="Prop", aggfunc="mean"
        )
    except Exception:
        return None, None
    if level_a not in pivot.columns or level_b not in pivot.columns:
        return None, None

    diff = pivot[level_a] - pivot[level_b]
    diff_wide = diff.unstack(level="TimeBin")
    bins = np.array(sorted(diff_wide.columns.tolist()))
    subj_ordered = sorted(diff_wide.index.tolist())
    diff_matrix = diff_wide.loc[subj_ordered, bins].values

    # Recover Time (ms) from the pivot — match a (subj, TimeBin) row to
    # its Time value in the source frame
    time_lookup = (sub_df.drop_duplicates(subset=["TimeBin"])
                   .set_index("TimeBin")["Time"])
    time_ms = np.array([float(time_lookup.loc[b]) for b in bins])

    return diff_matrix, time_ms


# =============================================================================
# Subsample generation
# =============================================================================

def generate_subsamples(n_total, n_target, n_draws, rng):
    """Up to n_draws distinct index subsamples of size n_target."""
    n_possible = comb(n_total, n_target)
    if n_possible <= n_draws:
        return [list(c) for c in itertools.combinations(range(n_total), n_target)]
    seen, result = set(), []
    while len(result) < n_draws:
        pick = tuple(sorted(rng.choice(n_total, size=n_target,
                                       replace=False).tolist()))
        if pick not in seen:
            seen.add(pick)
            result.append(list(pick))
    return result


# =============================================================================
# Run all three analyses on a subject subset
# =============================================================================

def run_analysis_on_subset(target_df, goal_df, comp_df,
                           pred_col_target, pred_col_goal, pred_col_comp,
                           comp_levels,
                           subjects, n_permutations, rng):
    """Run target, goal, and competitor cluster permutation tests."""
    results = {}

    # Target: conflict - no_conflict on aoi_target
    if target_df is not None:
        dm, t = build_diff_matrix(target_df, subjects, pred_col_target,
                                  "conflict", "no_conflict")
        if dm is not None and dm.shape[0] >= 2:
            results["target"] = cluster_permutation_test(
                dm, t, n_permutations, rng)
        else:
            results["target"] = _empty_cluster_result()

    # Goal: no_conflict - conflict on aoi_goal  (treatment_level='no_conflict')
    if goal_df is not None:
        dm, t = build_diff_matrix(goal_df, subjects, pred_col_goal,
                                  "no_conflict", "conflict")
        if dm is not None and dm.shape[0] >= 2:
            results["goal"] = cluster_permutation_test(
                dm, t, n_permutations, rng)
        else:
            results["goal"] = _empty_cluster_result()

    # Competitor: aoi_comp - aoi_AllOther on conflict only
    if comp_df is not None:
        dm, t = build_diff_matrix(comp_df, subjects, pred_col_comp,
                                  comp_levels[0], comp_levels[1])
        if dm is not None and dm.shape[0] >= 2:
            results["competitor"] = cluster_permutation_test(
                dm, t, n_permutations, rng)
        else:
            results["competitor"] = _empty_cluster_result()

    return results


# =============================================================================
# Stability metrics
# =============================================================================

def compute_stability_metrics(long_df, full_results, output_dir):
    analyses = long_df["analysis"].unique()
    summary_rows = []
    for analysis in analyses:
        sub = long_df[long_df["analysis"] == analysis]
        for N, grp in sub.groupby("N"):
            row = {
                "analysis": analysis, "N": N, "n_draws": len(grp),
                "significance_rate": grp["significant"].mean(),
                "max_abs_cluster_sum_mean": grp["max_abs_cluster_sum"].mean(),
                "max_abs_cluster_sum_sd": grp["max_abs_cluster_sum"].std(),
                "max_abs_cluster_sum_median": grp["max_abs_cluster_sum"].median(),
                "max_abs_cluster_sum_iqr_low":
                    grp["max_abs_cluster_sum"].quantile(0.25),
                "max_abs_cluster_sum_iqr_high":
                    grp["max_abs_cluster_sum"].quantile(0.75),
                "n_clusters_mean": grp["n_clusters"].mean(),
                "cluster_start_ms_mean": grp["cluster_start_ms"].mean(),
                "cluster_end_ms_mean": grp["cluster_end_ms"].mean(),
                "cluster_n_bins_mean": grp["cluster_n_bins"].mean(),
            }
            if full_results and analysis in full_results:
                full = full_results[analysis]
                orig = full.get("max_abs_cluster_sum", np.nan)
                row["original_max_abs_cluster_sum"] = orig
                row["original_significant"] = full.get("significant", False)
                row["original_in_iqr"] = bool(
                    (not np.isnan(orig)) and
                    row["max_abs_cluster_sum_iqr_low"] <= orig
                    <= row["max_abs_cluster_sum_iqr_high"]
                )
            summary_rows.append(row)

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(os.path.join(output_dir, "gaze_stability_summary.csv"),
                   index=False)

    threshold_rows = []
    for analysis in analyses:
        sub = summary[summary["analysis"] == analysis]
        crossed = sub[sub["significance_rate"] >= 0.95]
        thr = int(crossed["N"].min()) if not crossed.empty else None
        threshold_rows.append({"analysis": analysis, "threshold_N_95": thr})
    pd.DataFrame(threshold_rows).to_csv(
        os.path.join(output_dir, "gaze_sample_size_threshold.csv"),
        index=False)
    return summary


# =============================================================================
# Top-level driver
# =============================================================================

def run_sensitivity(target_csv, goal_csv, comp_csv,
                    pred_col_target, pred_col_goal, pred_col_comp,
                    comp_levels,
                    n_min, n_max, n_draws, n_permutations,
                    output_dir, seed):

    rng = np.random.default_rng(seed)
    os.makedirs(output_dir, exist_ok=True)

    print("Loading aggregated CSVs...")
    target_df = goal_df = comp_df = None
    if target_csv:
        target_df, _ = _load_time_cluster_csv(target_csv)
        print(f"  target:     {len(target_df):>6} rows from {target_csv}")
    if goal_csv:
        goal_df, _ = _load_time_cluster_csv(goal_csv)
        print(f"  goal:       {len(goal_df):>6} rows from {goal_csv}")
    if comp_csv:
        comp_df, _ = _load_time_cluster_csv(comp_csv)
        print(f"  competitor: {len(comp_df):>6} rows from {comp_csv}")

    # Subjects: union across loaded frames
    sub_sets = [set(d["subj"].unique()) for d in (target_df, goal_df, comp_df)
                if d is not None]
    if not sub_sets:
        raise ValueError("No input CSVs provided.")
    all_subjects = sorted(set.union(*sub_sets))
    n_total = len(all_subjects)
    print(f"  subjects ({n_total}): {all_subjects}")

    # Full-sample reference (cached)
    full_csv = os.path.join(output_dir, "gaze_full_sample_results.csv")
    if os.path.exists(full_csv):
        print("Full-sample reference already exists, loading...")
        full_results = {row["analysis"]: dict(row)
                        for _, row in pd.read_csv(full_csv).iterrows()}
    else:
        print(f"Running full-sample reference ({n_total} subjects, "
              f"{n_permutations} permutations)...")
        full_results = run_analysis_on_subset(
            target_df, goal_df, comp_df,
            pred_col_target, pred_col_goal, pred_col_comp, comp_levels,
            all_subjects, n_permutations, rng)
        pd.DataFrame([{"analysis": k, **v}
                       for k, v in full_results.items()]).to_csv(
            full_csv, index=False)

    print("\nFull-sample results:")
    for analysis, res in full_results.items():
        start, end = res.get("cluster_start_ms"), res.get("cluster_end_ms")
        cluster_str = (f"{int(start)}-{int(end)} ms"
                       if not (isinstance(start, float) and np.isnan(start))
                       else "no cluster")
        print(f"  {analysis:12s}  significant={res['significant']},  "
              f"p={res['p_value']:.4f},  cluster={cluster_str}")

    # Subsampling loop
    sensitivity_csv = os.path.join(output_dir, "gaze_sensitivity_results.csv")
    indices_csv = os.path.join(output_dir, "gaze_subsample_indices.csv")
    first_write = True

    for N in range(n_min, min(n_max, n_total) + 1):
        # At N == n_total there is exactly one possible "subsample"
        # (the full sample itself). Re-running the permutation analysis
        # would only differ from the cached full-sample reference by
        # Monte-Carlo noise. Reusing `full_results` here guarantees that
        # the N=n_total row of gaze_sensitivity_results matches the values
        # reported in the manuscript.
        if N == n_total:
            print(f"\n=== N = {N}: reusing full-sample reference ===")
            picked = list(all_subjects)
            write_mode = "w" if first_write else "a"
            pd.DataFrame([{
                "N": N, "draw": 0,
                "subject_ids": ";".join(str(s) for s in picked),
            }]).to_csv(indices_csv, mode=write_mode,
                       header=first_write, index=False)
            rows = [{"N": N, "draw": 0, "analysis": k, **v}
                    for k, v in full_results.items()]
            pd.DataFrame(rows).to_csv(
                sensitivity_csv, mode=write_mode,
                header=first_write, index=False)
            first_write = False
            continue

        subsamples = generate_subsamples(n_total, N, n_draws, rng)
        n_actual = len(subsamples)
        print(f"\n=== N = {N}: {n_actual} draws ===")

        with tqdm(subsamples, desc=f"N={N}", unit="draw") as pbar:
            for d_idx, idx in enumerate(pbar):
                picked = [all_subjects[i] for i in idx]
                write_mode = "w" if first_write else "a"

                pd.DataFrame([{
                    "N": N, "draw": d_idx,
                    "subject_ids": ";".join(str(s) for s in picked),
                }]).to_csv(indices_csv, mode=write_mode,
                           header=first_write, index=False)

                sub_results = run_analysis_on_subset(
                    target_df, goal_df, comp_df,
                    pred_col_target, pred_col_goal, pred_col_comp,
                    comp_levels,
                    picked, n_permutations, rng)

                rows = [{"N": N, "draw": d_idx, "analysis": k, **v}
                        for k, v in sub_results.items()]
                pd.DataFrame(rows).to_csv(
                    sensitivity_csv, mode=write_mode,
                    header=first_write, index=False)
                first_write = False

    print("\nComputing stability metrics...")
    long_df = pd.read_csv(sensitivity_csv)
    compute_stability_metrics(long_df, full_results, output_dir)
    print(f"\nDone. Results in {output_dir}/")


def main():
    p = argparse.ArgumentParser(
        description="Subsampling sensitivity on eyetrackingR-aggregated data.")
    p.add_argument("--target_csv",
                   help="time_cluster_data CSV for target analysis "
                        "(aoi_target, predictor=condition)")
    p.add_argument("--goal_csv",
                   help="time_cluster_data CSV for goal analysis "
                        "(aoi_goal, predictor=condition)")
    p.add_argument("--comp_csv",
                   help="time_cluster_data CSV for competitor analysis "
                        "(aoi_fct: aoi_comp vs aoi_AllOther, conflict only)")
    p.add_argument("--pred_col_target", default="condition",
                   help="Predictor column in target CSV (default: condition)")
    p.add_argument("--pred_col_goal", default="condition",
                   help="Predictor column in goal CSV (default: condition)")
    p.add_argument("--pred_col_competitor", default="aoi_fct",
                   help="Predictor column in competitor CSV "
                        "(default: aoi_fct)")
    p.add_argument("--comp_level_a", default="aoi_comp",
                   help="Competitor analysis: level A "
                        "(default: aoi_comp)")
    p.add_argument("--comp_level_b", default="aoi_AllOther",
                   help="Competitor analysis: level B "
                        "(default: aoi_AllOther)")
    p.add_argument("--n_min", type=int, default=15)
    p.add_argument("--n_max", type=int, default=26)
    p.add_argument("--n_draws", type=int, default=10)
    p.add_argument("--n_permutations", type=int, default=N_PERMUTATIONS)
    p.add_argument("--output_dir", default="gaze_sensitivity_results")
    p.add_argument("--seed", type=int, default=42)
    args = p.parse_args()

    if not (args.target_csv or args.goal_csv or args.comp_csv):
        p.error("Provide at least one of --target_csv / --goal_csv / "
                "--comp_csv")

    run_sensitivity(
        target_csv=args.target_csv, goal_csv=args.goal_csv,
        comp_csv=args.comp_csv,
        pred_col_target=args.pred_col_target,
        pred_col_goal=args.pred_col_goal,
        pred_col_comp=args.pred_col_competitor,
        comp_levels=(args.comp_level_a, args.comp_level_b),
        n_min=args.n_min, n_max=args.n_max, n_draws=args.n_draws,
        n_permutations=args.n_permutations,
        output_dir=args.output_dir, seed=args.seed,
    )


if __name__ == "__main__":
    main()
