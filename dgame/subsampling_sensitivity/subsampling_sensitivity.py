"""
subsampling_sensitivity.py
================================================================================
Subsampling sensitivity analysis for the FRP paper.

Standalone script — no dependency on the DGAME Experiment class or
config.yml. Reads the reconstructed ERP CSVs and channel positions
directly from disk.

The statistical analysis matches the R scripts
(Ja_lm_permute_and_plot_fixations.r / Jb_lm_permute_and_plot_language.r)
exactly:
  - FIX mode: data filtered to fix_at=="target" before everything;
    formula: data ~ laterality * saggitality * condition * trial_time * baseline
  - N mode: mean_target_fixation rounded to 3 decimals;
    formula: data ~ laterality * saggitality * condition * mean_target_fixation * baseline
  - Baseline computed with channel in groupby, then averaged over channels
  - Time bins 0-900ms; N mode permutation on bins [4:8] (300-700) only
  - DV column is 'data' (not 'data_mean'); permutation shuffles 'data'
  - FDR per time_bin + global FDR (fdr_q_value_global used for significance)

Usage:
    python subsampling_sensitivity.py \
        --data_dir /path/to/unfold_out_csvs \
        --channel_positions channel_positions.csv \
        --mode both \
        --n_min 15 --n_max 26 --n_draws 10 \
        --output_dir sensitivity_results
================================================================================
"""

import argparse
import glob
import itertools
import os
import time as _time
from concurrent.futures import ThreadPoolExecutor, as_completed
from math import comb, floor

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from patsy import dmatrices as _patsy_dmatrices
from tqdm import tqdm

# =============================================================================
# Constants (matching the original J script)
# =============================================================================

ALPHA = 0.05
N_PERMUTATIONS = 2000

# Column names
CHANNEL_FIELD = "channel"
LATERALITY_FIELD = "laterality"
SAGGITALITY_FIELD = "saggitality"
# X/Y columns from channel_positions.csv after merge
LATERAL_INPUT_FIELD = "X"
SAGGITAL_INPUT_FIELD = "Y"


# =============================================================================
# Data loading
# =============================================================================

def load_channel_positions(csv_path):
    """Load channel_positions.csv (no header: label, X, Y, Z)."""
    return pd.read_csv(
        csv_path, header=None,
        names=[CHANNEL_FIELD, LATERAL_INPUT_FIELD, SAGGITAL_INPUT_FIELD, "Z"])


def _process_subject(args):
    """Worker: process one subject. Top-level so threads can call it cleanly."""
    (sid, csvs, mode, channel_coords,
     baseline_groupby, baseline_join_on,
     final_cols, final_groupby) = args

    subj_frames = []
    skipped = []
    for fpath in csvs:
        fname = os.path.basename(fpath)
        try:
            try:
                df = pd.read_csv(fpath, engine="c")
            except UnicodeDecodeError:
                df = pd.read_csv(fpath, engine="c", encoding="unicode_escape")
            df[CHANNEL_FIELD] = df[CHANNEL_FIELD].astype(str)
            subj_frames.append(df)
        except (OSError, pd.errors.ParserError, UnicodeDecodeError) as e:
            skipped.append((fname, str(e)))

    if not subj_frames:
        return sid, None, skipped, 0

    subj_data = pd.concat(subj_frames, axis=0, ignore_index=True)
    subj_data = subj_data.merge(channel_coords, how="left", on=CHANNEL_FIELD)
    subj_data = annotate_laterality_and_saggitality(subj_data)

    if mode == "FIX":
        subj_data = subj_data[subj_data["fix_at"] == "target"].copy()
        subj_data = subj_data.drop(columns=["fix_at"], errors="ignore")

    if mode == "N":
        subj_data["mean_target_fixation"] = subj_data["mean_target_fixation"].round(3)

    baseline_cols = ["data"] + baseline_groupby
    baseline_data = subj_data.loc[
        (subj_data["time"] >= -250) & (subj_data["time"] < 0)][baseline_cols]
    baseline_agg = (
        baseline_data
        .groupby(baseline_groupby, as_index=False)
        .agg(baseline=("data", "mean"))
    )

    subj_data["time_bin"] = subj_data["time"].apply(lambda t: floor(t / 100) * 100)
    subj_data = subj_data.merge(baseline_agg, how="inner", on=baseline_join_on)
    subj_data = subj_data.loc[
        (subj_data["time_bin"] >= 0) & (subj_data["time_bin"] <= 900)]

    subj_agg = (
        subj_data[final_cols]
        .groupby(final_groupby, as_index=False)
        .mean()
    )
    return sid, subj_agg, skipped, len(subj_agg)


def load_and_preprocess(data_dir, mode, channel_coords, n_workers=None):
    """Load, annotate, and aggregate data in parallel, one worker per subject.

    Processing order matches the R scripts exactly:
      1. Load subject CSVs, merge channel coords, annotate lat/sag
      2. FIX mode: filter entire dataset to fix_at == "target", drop fix_at
      3. Compute baseline WITH channel in the groupby
      4. Add time_bin, inner_join baseline (WITH channel)
      5. N mode: round mean_target_fixation to 3 decimals
      6. Filter time_bin >= 0 & <= 900
      7. Drop channel column
      8. group_by (without channel) + mean over both data and baseline

    Returns a DataFrame with columns:
        time_bin, laterality, saggitality, condition, subject,
        trial_time (FIX) or mean_target_fixation (N), data, baseline
    """
    # Find subject subdirectories
    subdirs = sorted([
        d for d in glob.glob(os.path.join(data_dir, "*"))
        if os.path.isdir(d)
    ])
    flat_files = sorted(glob.glob(os.path.join(data_dir, f"*_unfold_{mode}.csv")))

    if not subdirs and not flat_files:
        raise FileNotFoundError(
            f"No subject subdirectories or *_unfold_{mode}.csv files "
            f"found in {data_dir}")

    # Build mapping: subject_id -> list of CSV files
    subject_files = {}
    if subdirs:
        for subdir in subdirs:
            sid = os.path.basename(subdir)
            csvs = sorted(glob.glob(
                os.path.join(subdir, f"*_unfold_{mode}.csv")))
            if not csvs:
                # Pipeline writes CSVs one level deeper under unfold_out/
                csvs = sorted(glob.glob(
                    os.path.join(subdir, "unfold_out", f"*_unfold_{mode}.csv")))
            if csvs:
                subject_files[sid] = csvs
    else:
        for fpath in flat_files:
            fname = os.path.basename(fpath)
            sid = fname.split("_")[0]
            subject_files.setdefault(sid, []).append(fpath)

    print(f"  Found {sum(len(v) for v in subject_files.values())} CSV files "
          f"across {len(subject_files)} subjects for mode={mode}")

    # Column definitions matching R exactly
    if mode == "FIX":
        # R: group_by(subject, condition, laterality, saggitality, channel, trial_time)
        baseline_groupby = ["subject", "condition", LATERALITY_FIELD,
                            SAGGITALITY_FIELD, CHANNEL_FIELD, "trial_time"]
        baseline_join_on = ["subject", "condition", CHANNEL_FIELD,
                            LATERALITY_FIELD, SAGGITALITY_FIELD, "trial_time"]
        # R: select(time_bin, laterality, saggitality, condition, trial_time, subject, data, baseline)
        final_cols = ["time_bin", LATERALITY_FIELD, SAGGITALITY_FIELD,
                      "condition", "trial_time", "subject", "data", "baseline"]
        # R: group_by(time_bin, laterality, saggitality, condition, trial_time, subject)
        final_groupby = ["time_bin", LATERALITY_FIELD, SAGGITALITY_FIELD,
                         "condition", "trial_time", "subject"]
    elif mode == "N":
        # R: group_by(subject, condition, channel, laterality, saggitality, mean_target_fixation)
        baseline_groupby = ["subject", "condition", CHANNEL_FIELD,
                            LATERALITY_FIELD, SAGGITALITY_FIELD,
                            "mean_target_fixation"]
        baseline_join_on = ["subject", "condition", LATERALITY_FIELD,
                            SAGGITALITY_FIELD, CHANNEL_FIELD,
                            "mean_target_fixation"]
        # R: select(time_bin, laterality, saggitality, condition, mean_target_fixation, subject, baseline, data)
        final_cols = ["time_bin", LATERALITY_FIELD, SAGGITALITY_FIELD,
                      "condition", "mean_target_fixation", "subject",
                      "baseline", "data"]
        # R: group_by(time_bin, laterality, saggitality, condition, mean_target_fixation, subject)
        final_groupby = ["time_bin", LATERALITY_FIELD, SAGGITALITY_FIELD,
                         "condition", "mean_target_fixation", "subject"]
    else:
        raise ValueError(f"mode must be 'FIX' or 'N', got '{mode}'")

    aggregated_frames = []
    skipped_files = []
    _t0 = _time.perf_counter()

    worker_args = [
        (sid, csvs, mode, channel_coords,
         baseline_groupby, baseline_join_on,
         final_cols, final_groupby)
        for sid, csvs in subject_files.items()
    ]
    n_subjects = len(worker_args)
    workers = min(n_workers or os.cpu_count() or 1, n_subjects)

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(_process_subject, a): a[0] for a in worker_args}
        with tqdm(as_completed(futures), total=n_subjects,
                  unit="subj", desc="Loading subjects") as pbar:
            for fut in pbar:
                sid = futures[fut]
                try:
                    sid_out, subj_agg, skipped, n_rows = fut.result()
                except Exception as exc:
                    print(f"\n  ERROR subject {sid}: {exc}")
                    continue
                for fname, err in skipped:
                    print(f"\n  FAILED {fname}: {err}")
                    skipped_files.append(fname)
                if subj_agg is None:
                    print(f"\n  WARNING: no data for subject {sid_out}, skipping")
                    continue
                aggregated_frames.append(subj_agg)
                pbar.set_postfix({"last": sid_out, "rows": n_rows})

    if skipped_files:
        print(f"\n  WARNING: {len(skipped_files)} file(s) skipped:")
        for s in skipped_files:
            print(f"    - {s}")

    if not aggregated_frames:
        raise RuntimeError("No data loaded successfully.")

    all_data = pd.concat(aggregated_frames, axis=0, ignore_index=True)
    print(f"  Total: {len(all_data)} aggregated rows "
          f"from {len(aggregated_frames)} subjects "
          f"in {_time.perf_counter() - _t0:.1f}s ({workers} workers)")

    return all_data


# =============================================================================
# Annotation (identical to original J script)
# =============================================================================

def annotate_laterality_and_saggitality(df):
    """Annotate DataFrame with laterality and saggitality labels.
    Identical logic to the original J script."""

    df[LATERALITY_FIELD] = np.where(
        df[LATERAL_INPUT_FIELD] < 0, "left",
        np.where(df[LATERAL_INPUT_FIELD] > 0, "right", "central")
    )

    sag_conditions = [
        (df[SAGGITAL_INPUT_FIELD] > 0) & (df[SAGGITAL_INPUT_FIELD] <= 0.0714),
        (df[SAGGITAL_INPUT_FIELD] > 0.0714),
        (df[SAGGITAL_INPUT_FIELD] < 0) & (df[SAGGITAL_INPUT_FIELD] >= -0.0929),
        (df[SAGGITAL_INPUT_FIELD] < -0.0929),
    ]
    sag_labels = ["frontal", "prefrontal", "posterior", "occipital"]
    df[SAGGITALITY_FIELD] = np.select(sag_conditions, sag_labels, default="central")

    return df


# =============================================================================
# Statistics (matching the R scripts)
# =============================================================================

def summarize_stats_model(model):
    """Extract predictor names, coefficients, and t-values from a
    fitted statsmodels OLS result."""
    return pd.DataFrame({
        "predictor": model.params.index.tolist(),
        "coefficient": model.params.values,
        "std_err": model.bse.values,
        "t_value": model.tvalues.values,
        "p_value": model.pvalues.values,
    })


def fdr_adjust_pvals(p_values):
    """Benjamini-Hochberg FDR adjustment. Returns q-values as a list."""
    p = np.array(p_values)
    n = len(p)
    if n == 0:
        return []
    order = np.argsort(p)
    ranked = p[order]
    rank = np.arange(1, n + 1)
    q = np.minimum(1.0, ranked * n / rank)
    for i in range(n - 2, -1, -1):
        q[i] = min(q[i], q[i + 1])
    q_out = np.empty(n)
    q_out[order] = q
    return q_out.tolist()


def time_bin_permutation(filtered_data, mode, time_bin,
                         n_permutations, include_baseline):
    """Permutation test for one time bin.

    Matches the R scripts exactly:
      - FIX formula: data ~ laterality * saggitality * condition * trial_time * baseline
      - N formula:   data ~ laterality * saggitality * condition * mean_target_fixation * baseline
      - DV is 'data' (not 'data_mean')
      - Shuffles 'data' column (R: g_perm$data <- sample(g_perm$data))

    Optimised: design matrix X is built once via patsy; QR decomposition is
    reused across all permutations so the entire permutation loop collapses
    into a single batched matrix solve.
    """
    if mode == "FIX":
        model_formula = (f"data ~ {LATERALITY_FIELD} * {SAGGITALITY_FIELD}"
                         f" * condition * trial_time")
    elif mode == "N":
        model_formula = (f"data ~ {LATERALITY_FIELD} * {SAGGITALITY_FIELD}"
                         f" * condition * mean_target_fixation")
    if include_baseline:
        model_formula += " * baseline"

    # Build design matrix ONCE — X stays identical across all permutations
    y_df, X_df = _patsy_dmatrices(model_formula, filtered_data,
                                   return_type="dataframe")
    y = np.asarray(y_df).ravel()
    X = np.asarray(X_df)
    coef_names = list(X_df.columns)
    n_obs, n_params = X.shape

    # QR decomposition once; precompute diagonal of (X'X)^{-1} for std errors
    Q, R = np.linalg.qr(X)
    R_inv = np.linalg.solve(R, np.eye(n_params))
    diag_XtX_inv = np.einsum("ij,ij->i", R_inv, R_inv)

    def _tstats_batch(Y):
        """Y: (n_draws, n_obs) — returns t-stats (n_params, n_draws)."""
        QtY  = Q.T @ Y.T                                         # (p, n_draws)
        beta = np.linalg.solve(R, QtY)                          # (p, n_draws)
        resid = Y.T - X @ beta                                   # (n, n_draws)
        sigma2 = np.sum(resid ** 2, axis=0) / (n_obs - n_params) # (n_draws,)
        var_b = diag_XtX_inv[:, None] * sigma2[None, :]          # (p, n_draws)
        return beta / np.sqrt(var_b)

    # Observed t-stats
    obs_tstats = _tstats_batch(y[None, :]).ravel()

    # All permutations at once — vectorised, no Python loop
    perm_Y = np.tile(y, (n_permutations, 1))
    np.random.default_rng(seed=0).permuted(perm_Y, axis=1, out=perm_Y)
    perm_tstats = _tstats_batch(perm_Y)  # (n_params, n_permutations)

    p_vals = np.mean(np.abs(perm_tstats) >= np.abs(obs_tstats)[:, None], axis=1)

    return pd.DataFrame({
        "time_bin": time_bin,
        "predictor": coef_names,
        "observed_tstat": obs_tstats,
        "permutation_p_value": p_vals,
    })


def block_permutation_test(time_windowed_data, mode,
                           n_permutations, include_baseline):
    """Run permutation tests across time bins with FDR.

    Matches R time bin selection:
      - FIX: all time bins
      - N: time_bins[4:8] only (i.e. 300, 400, 500, 600, 700)

    FDR is applied per time_bin (R: group_by(time_bin) %>%
    mutate(fdr_q_value = p.adjust(..., "fdr"))), plus a global FDR
    across all bins (R: fdr_q_value_global).
    """
    all_bins = sorted(time_windowed_data["time_bin"].unique())

    if mode == "N":
        # R: unique(dat_n_windows$time_bin)[4:8] — 1-indexed, so Python [3:8]
        time_bins = all_bins[3:8] if len(all_bins) >= 8 else all_bins
    else:
        # FIX: all bins
        time_bins = all_bins

    results = []
    with tqdm(time_bins, unit="bin") as pbar:
        for time_bin in pbar:
            pbar.set_description(f"Permutation test (time_bin={time_bin})")
            filtered = time_windowed_data[
                time_windowed_data["time_bin"] == time_bin]
            result = time_bin_permutation(
                filtered, mode, time_bin, n_permutations, include_baseline)
            results.append(result)

    combined = pd.concat(results, axis=0, ignore_index=True)

    # FDR per time_bin (R: group_by(time_bin) %>% mutate(fdr_q_value = p.adjust(...)))
    combined["fdr_q_value"] = (
        combined
        .groupby("time_bin")["permutation_p_value"]
        .transform(lambda x: fdr_adjust_pvals(x.tolist()))
    )
    # Global FDR across all bins (R: fdr_q_value_global)
    combined["fdr_q_value_global"] = fdr_adjust_pvals(
        combined["permutation_p_value"].tolist())

    return combined


# =============================================================================
# Run analysis on a subject subset
# =============================================================================

def run_analysis_on_subset(preprocessed_data, subject_ids_subset, mode,
                           n_permutations=N_PERMUTATIONS,
                           include_baseline=True):
    """Run permutation + FDR on a subset of subjects.

    The input is ALREADY time-windowed, baseline-merged, and filtered
    (output of load_and_preprocess). This function just filters by
    subject and runs the statistical tests.
    """
    data = preprocessed_data[
        preprocessed_data["subject"].isin(subject_ids_subset)
    ].copy()

    results = block_permutation_test(
        data, mode, n_permutations, include_baseline)

    return results


# =============================================================================
# Subsample generation
# =============================================================================

def generate_subsamples(n_total, n_target, n_draws, rng):
    """Generate up to n_draws distinct subsamples of size n_target."""
    n_possible = comb(n_total, n_target)
    if n_possible <= n_draws:
        return [list(c) for c in itertools.combinations(range(n_total), n_target)]
    seen = set()
    result = []
    while len(result) < n_draws:
        pick = tuple(sorted(rng.choice(n_total, size=n_target, replace=False)))
        if pick not in seen:
            seen.add(pick)
            result.append(list(pick))
    return result


# =============================================================================
# Stability metrics
# =============================================================================

def compute_stability_metrics(long_df, full_results, output_dir, mode_label):
    """Aggregate subsampling results into stability summary."""

    # Use fdr_q_value_global for significance (matching R:
    # filter(fdr_q_value_global <= 0.05))
    long_df["significant"] = long_df["fdr_q_value_global"] < ALPHA

    summary = (
        long_df
        .groupby(["predictor", "time_bin", "N"])
        .agg(
            significance_rate=("significant", "mean"),
            tstat_mean=("observed_tstat", "mean"),
            tstat_sd=("observed_tstat", "std"),
            tstat_median=("observed_tstat", "median"),
            tstat_iqr_low=("observed_tstat", lambda x: np.percentile(x, 25)),
            tstat_iqr_high=("observed_tstat", lambda x: np.percentile(x, 75)),
            n_draws=("observed_tstat", "size"),
        )
        .reset_index()
    )

    # Attach original full-sample t-stats
    orig = (
        full_results[["predictor", "time_bin", "observed_tstat"]]
        .rename(columns={"observed_tstat": "original_tstat"})
    )
    summary = summary.merge(orig, on=["predictor", "time_bin"], how="left")
    summary["original_in_iqr"] = (
        (summary["original_tstat"] >= summary["tstat_iqr_low"]) &
        (summary["original_tstat"] <= summary["tstat_iqr_high"])
    )

    summary.to_csv(
        os.path.join(output_dir, f"stability_summary_{mode_label}.csv"),
        index=False)

    # Sample-size threshold
    threshold_rows = []
    for (pred, tbin), grp in summary.groupby(["predictor", "time_bin"]):
        crossed = grp[grp["significance_rate"] >= 0.95]
        thr = int(crossed["N"].min()) if not crossed.empty else None
        threshold_rows.append({
            "predictor": pred, "time_bin": tbin, "threshold_N_95": thr})
    pd.DataFrame(threshold_rows).to_csv(
        os.path.join(output_dir, f"sample_size_threshold_{mode_label}.csv"),
        index=False)


# =============================================================================
# Main
# =============================================================================

def run_sensitivity(data_dir, channel_positions_path, mode,
                    n_min=15, n_max=26, n_draws=10,
                    n_permutations=N_PERMUTATIONS,
                    include_baseline=True,
                    output_dir="sensitivity_results",
                    seed=42):
    """Top-level driver."""

    rng = np.random.default_rng(seed)
    os.makedirs(output_dir, exist_ok=True)
    mode_label = "fixations" if mode == "FIX" else "language"

    # Load channel positions
    channel_coords = load_channel_positions(channel_positions_path)
    print(f"  Channel positions: {len(channel_coords)} entries")

    # Load and preprocess per-subject (saves RAM)
    print(f"Loading and preprocessing data ({mode_label})...")
    all_data = load_and_preprocess(data_dir, mode, channel_coords)
    all_subject_ids = sorted(all_data["subject"].unique())
    n_total = len(all_subject_ids)
    print(f"  {n_total} subjects, {len(all_data)} aggregated rows")

    # Full-sample reference
    print("Running full-sample reference...")
    full_csv = os.path.join(output_dir, f"original_full_sample_{mode_label}.csv")
    if os.path.exists(full_csv):
        print("Full-sample reference already exists, skipping...")
        full_results = pd.read_csv(full_csv)
    else:
        print("Running full-sample reference...")
        full_results = run_analysis_on_subset(
            all_data, all_subject_ids, mode, n_permutations, include_baseline)
        full_results.to_csv(full_csv, index=False)

    # Subsampling loop — results written to disk after each draw to avoid RAM buildup
    sensitivity_csv = os.path.join(output_dir, f"sensitivity_results_{mode_label}.csv")
    indices_csv = os.path.join(output_dir, f"subsample_indices_{mode_label}.csv")
    first_write = True

    for N in range(n_min, min(n_max, n_total) + 1):
        # At N == n_total there is exactly one possible "subsample"
        # (the full sample itself). Re-running the permutation analysis
        # on it would produce a result that differs from the cached
        # full-sample reference only by Monte-Carlo noise in the
        # permutation null. Reusing `full_results` here guarantees that
        # the N=n_total row in sensitivity_results matches the values
        # reported in the manuscript and avoids spurious sig-rate jumps.
        if N == n_total:
            print(f"\n=== N = {N}: reusing full-sample reference ===")
            picked_ids = list(all_subject_ids)
            write_mode = "w" if first_write else "a"
            pd.DataFrame([{
                "N": N, "draw": 0,
                "subject_ids": ";".join(str(s) for s in picked_ids),
            }]).to_csv(indices_csv, mode=write_mode,
                       header=first_write, index=False)
            result = full_results.copy()
            result["N"] = N
            result["draw"] = 0
            result.to_csv(sensitivity_csv, mode=write_mode,
                          header=first_write, index=False)
            first_write = False
            continue

        subsamples = generate_subsamples(n_total, N, n_draws, rng)
        n_actual = len(subsamples)
        print(f"\n=== N = {N}: {n_actual} draws ===")

        for d_idx, idx in enumerate(subsamples):
            picked_ids = [all_subject_ids[i] for i in idx]
            write_mode = "w" if first_write else "a"
            pd.DataFrame([{
                "N": N, "draw": d_idx,
                "subject_ids": ";".join(str(s) for s in picked_ids),
            }]).to_csv(indices_csv, mode=write_mode, header=first_write, index=False)

            print(f"  draw {d_idx+1}/{n_actual} (N={N})...")
            result = run_analysis_on_subset(
                all_data, picked_ids, mode, n_permutations, include_baseline)
            result["N"] = N
            result["draw"] = d_idx
            result.to_csv(sensitivity_csv, mode=write_mode, header=first_write, index=False)
            first_write = False

    # Stability metrics — read back from disk instead of keeping all results in RAM
    long_df = pd.read_csv(sensitivity_csv)
    compute_stability_metrics(long_df, full_results, output_dir, mode_label)

    print(f"\nDone. Results in {output_dir}")


def main():
    parser = argparse.ArgumentParser(
        description="Subsampling sensitivity analysis for the FRP paper")
    parser.add_argument("--data_dir", required=True,
                        help="Directory with *_unfold_FIX.csv / "
                             "*_unfold_N.csv files (searched recursively)")
    parser.add_argument("--channel_positions", required=True,
                        help="Path to channel_positions.csv (label,X,Y,Z)")
    parser.add_argument("--mode", choices=["FIX", "N", "both"], default="both")
    parser.add_argument("--n_min", type=int, default=15)
    parser.add_argument("--n_max", type=int, default=26)
    parser.add_argument("--n_draws", type=int, default=10)
    parser.add_argument("--n_permutations", type=int, default=N_PERMUTATIONS)
    parser.add_argument("--include_baseline", action="store_true", default=True)
    parser.add_argument("--no_baseline", dest="include_baseline",
                        action="store_false")
    parser.add_argument("--output_dir", default="sensitivity_results")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    modes = ["FIX", "N"] if args.mode == "both" else [args.mode]

    for mode in modes:
        print(f"\n{'='*60}")
        print(f"  Sensitivity analysis: mode = {mode}")
        print(f"{'='*60}")
        run_sensitivity(
            data_dir=args.data_dir,
            channel_positions_path=args.channel_positions,
            mode=mode,
            n_min=args.n_min,
            n_max=args.n_max,
            n_draws=args.n_draws,
            n_permutations=args.n_permutations,
            include_baseline=args.include_baseline,
            output_dir=args.output_dir,
            seed=args.seed,
        )


if __name__ == "__main__":
    main()