"""
plot_sensitivity_main_figures.py
================================

Standalone plotting script for the three main sensitivity-analysis figures of
the BRM revision. Each figure has 6 subplots arranged 2 rows × 3 columns.

  fig_sensitivity_gaze.pdf
    Top row    — significance rate × N (target / competitor / goal)
                 alpha = .05
    Bottom row — cluster-sum stability × N (target / competitor / goal)
                 median ± IQR over draws, full-sample reference

  fig_sensitivity_frp.pdf  (and analogously fig_sensitivity_language.pdf)
    Top row    — significance rate × N for three central interactions:
                 condition × slope,
                 condition × slope × saggitality[occipital],
                 condition × slope × saggitality[prefrontal]
                 (slope = trial_time for FRP, mean_target_fixation for LRP)
                 one line per time bin, alpha = .05
    Bottom row — t-stat stability × N for the same three predictors
                 one line per time bin, median ± IQR over draws,
                 per-bin full-sample t-stat as reference

In all subplots:
  - Three horizontal stability reference lines at .80, .90, .95 (sig-rate row)
  - One vertical reference line at N_max (full sample)

Inputs (CLI args):
  --gaze_csv             raw per-draw gaze results
                         (gaze_sensitivity_results.csv)
  --gaze_full_csv        full-sample gaze results
                         (gaze_full_sample_results.csv)
  --frp_csv              long per-draw FRP results
                         (sensitivity_results_fixations.csv)
  --frp_full_csv         full-sample FRP results
                         (original_full_sample_fixations.csv); optional
  --language_csv         long per-draw LRP results
                         (sensitivity_results_language.csv)
  --language_full_csv    full-sample LRP results
                         (original_full_sample_language.csv); optional
  --output_dir           where to write the three PDFs

Notes:
  - For the ERP plots, significance is determined by fdr_q_value_global < alpha
    (matching the original analysis pipeline).
  - For the gaze plot, significance is permutation_p_value < alpha.
  - The script does NOT recompute thresholds; it only visualises sig-rates
    and t-stat distributions across subsamples.
"""

from __future__ import annotations

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------

ALPHAS = [0.05]
ALPHA_COLOR = {0.05: "#c0392b"}
ALPHA_LABEL = {0.05: r"$\alpha = .05$"}

STAB_LEVELS = [0.80, 0.90, 0.95]
STAB_LINESTYLE = {0.80: (0, (5, 5)), 0.90: (0, (3, 2)), 0.95: "solid"}
STAB_COLOR = "#888888"

FULL_SAMPLE_COLOR = "#1f77b4"
FULL_SAMPLE_LS = (0, (4, 2))

TSTAT_COLOR = "#2c3e50"
TSTAT_FS_COLOR = "#c0392b"

# ---------------------------------------------------------------------------
# Predictor selections
# ---------------------------------------------------------------------------

GAZE_ANALYSES = [("target", "Target AOI"),
                 ("competitor", "Competitor AOI"),
                 ("goal", "Goal AOI")]

FRP_PREDICTORS = [
    ("condition[T.no_conflict]:trial_time",
     "condition × trial_time"),
    ("saggitality[T.occipital]:condition[T.no_conflict]:trial_time",
     "occipital × condition × trial_time"),
    ("saggitality[T.prefrontal]:condition[T.no_conflict]:trial_time",
     "prefrontal × condition × trial_time"),
]

LANGUAGE_PREDICTORS = [
    ("condition[T.no_conflict]:mean_target_fixation",
     "condition × mean_target_fixation"),
    ("saggitality[T.occipital]:condition[T.no_conflict]:mean_target_fixation",
     "occipital × condition × mean_target_fixation"),
    ("saggitality[T.prefrontal]:condition[T.no_conflict]:mean_target_fixation",
     "prefrontal × condition × mean_target_fixation"),
]


# ---------------------------------------------------------------------------
# Subplot helpers
# ---------------------------------------------------------------------------

def add_stability_lines(ax):
    for stab in STAB_LEVELS:
        ax.axhline(stab, color=STAB_COLOR, lw=0.7,
                   ls=STAB_LINESTYLE[stab], alpha=0.7, zorder=0)


def add_full_sample_line(ax, n_max):
    ax.axvline(n_max, color=FULL_SAMPLE_COLOR, lw=1.0,
               ls=FULL_SAMPLE_LS, alpha=0.6, zorder=0)


def style_sigrate_axis(ax, n_min, n_max, title=None,
                       ylabel="Significance rate"):
    ax.set_ylim(-0.02, 1.05)
    ax.set_xlim(n_min - 0.4, n_max + 0.4)
    ax.set_xticks(np.arange(n_min, n_max + 1, 3))
    ax.set_xlabel("Sample size N", fontsize=9)
    ax.set_ylabel(ylabel, fontsize=9)
    if title is not None:
        ax.set_title(title, fontsize=9)
    ax.tick_params(axis="both", labelsize=8)
    ax.grid(True, axis="y", alpha=0.25, lw=0.4)


def style_value_axis(ax, n_min, n_max, title=None, ylabel="t-statistic"):
    ax.set_xlim(n_min - 0.4, n_max + 0.4)
    ax.set_xticks(np.arange(n_min, n_max + 1, 3))
    ax.set_xlabel("Sample size N", fontsize=9)
    ax.set_ylabel(ylabel, fontsize=9)
    if title is not None:
        ax.set_title(title, fontsize=9)
    ax.tick_params(axis="both", labelsize=8)
    ax.grid(True, axis="y", alpha=0.25, lw=0.4)
    ax.axhline(0, color="#bbbbbb", lw=0.5, zorder=0)


def reference_legend_handles(include_alphas=True, alpha_set=ALPHAS,
                             include_full_sample=True,
                             include_iqr=False, include_minmax=False):
    handles = []
    if include_alphas:
        for a in alpha_set:
            handles.append(Line2D([0], [0], color=ALPHA_COLOR[a], lw=1.6,
                                  label=ALPHA_LABEL[a]))
    for stab in STAB_LEVELS:
        handles.append(Line2D([0], [0], color=STAB_COLOR, lw=0.9,
                              ls=STAB_LINESTYLE[stab],
                              label=f"{int(stab*100)}% stability"))
    if include_iqr:
        handles.append(Patch(facecolor="#888888", alpha=0.25, label="IQR"))
    if include_minmax:
        handles.append(Line2D([0], [0], color="#888888", lw=0.7,
                              ls=(0, (1, 1)), label="min / max"))
    #if include_full_sample:
    #    handles.append(Line2D([0], [0], color=FULL_SAMPLE_COLOR, lw=1.0,
    #                          ls=FULL_SAMPLE_LS,
    #                          label="Full sample (N max)"))
    return handles


def time_bin_colors(bins, cmap_name="viridis"):
    """Build a colour mapping for the time bins, ordered ascending."""
    bins = sorted(bins)
    cmap = plt.get_cmap(cmap_name)
    if len(bins) == 1:
        return {bins[0]: cmap(0.5)}
    return {b: cmap(i / (len(bins) - 1)) for i, b in enumerate(bins)}


def time_bin_legend_handles(bins, colors):
    return [Line2D([0], [0], color=colors[b], lw=1.6,
                   label=f"{int(b)} ms") for b in sorted(bins)]


# ---------------------------------------------------------------------------
# IO helpers — minimise memory by reading only the columns we need
# ---------------------------------------------------------------------------

def _read_erp_long_for_sigrate(path, predictors):
    df = pd.read_csv(path,
                     usecols=["predictor", "time_bin", "N",
                              "fdr_q_value_global"])
    return df[df["predictor"].isin(predictors)].copy()


def _read_erp_long_for_tstat(path, predictors):
    df = pd.read_csv(path,
                     usecols=["predictor", "time_bin", "N", "observed_tstat"])
    return df[df["predictor"].isin(predictors)].copy()


# ---------------------------------------------------------------------------
# Aggregations
# ---------------------------------------------------------------------------

def compute_gaze_sig_rates(df_raw, alphas=ALPHAS):
    frames = []
    for a in alphas:
        rate = (df_raw
                .assign(sig=(df_raw["p_value"] < a).astype(float))
                .groupby(["analysis", "N"], as_index=False)
                .agg(sig_rate=("sig", "mean")))
        rate["alpha"] = a
        frames.append(rate)
    return pd.concat(frames, ignore_index=True)


def compute_gaze_cluster_bands(df_raw):
    return (df_raw.groupby(["analysis", "N"], as_index=False)
            .agg(median=("max_abs_cluster_sum", "median"),
                 q25=("max_abs_cluster_sum",
                      lambda x: np.percentile(x, 25)),
                 q75=("max_abs_cluster_sum",
                      lambda x: np.percentile(x, 75))))


def compute_gaze_extent_bands(df_raw):
    """Cluster start/end times per (analysis, N), restricted to draws that
    actually had a detected cluster. Returns median ± IQR bands for both
    boundaries."""
    sig = df_raw[df_raw["max_abs_cluster_sum"] > 0].copy()
    if sig.empty:
        return pd.DataFrame(columns=[
            "analysis", "N",
            "start_median", "start_q25", "start_q75",
            "end_median", "end_q25", "end_q75"])
    return (sig.groupby(["analysis", "N"], as_index=False)
            .agg(start_median=("cluster_start_ms", "median"),
                 start_q25=("cluster_start_ms",
                            lambda x: np.percentile(x, 25)),
                 start_q75=("cluster_start_ms",
                            lambda x: np.percentile(x, 75)),
                 end_median=("cluster_end_ms", "median"),
                 end_q25=("cluster_end_ms",
                          lambda x: np.percentile(x, 25)),
                 end_q75=("cluster_end_ms",
                          lambda x: np.percentile(x, 75))))


def compute_erp_sig_rates(df_long, alphas=ALPHAS):
    if df_long.empty:
        return pd.DataFrame(columns=[
            "predictor", "time_bin", "N", "sig_rate", "alpha"])
    frames = []
    for a in alphas:
        rate = (df_long
                .assign(sig=(df_long["fdr_q_value_global"] < a).astype(float))
                .groupby(["predictor", "time_bin", "N"], as_index=False)
                .agg(sig_rate=("sig", "mean")))
        rate["alpha"] = a
        frames.append(rate)
    return pd.concat(frames, ignore_index=True)


def compute_erp_tstat_bands(df_long):
    """t-stat distribution per (predictor, time_bin, N), aggregating over
    draws. The long ERP CSV does not carry a per-draw identifier, so the
    set of t-stats per (predictor, time_bin, N) cell is the draw distribution
    at that cell."""
    if df_long.empty:
        return pd.DataFrame(columns=[
            "predictor", "time_bin", "N", "median", "q25", "q75"])
    return (df_long.groupby(["predictor", "time_bin", "N"], as_index=False)
            .agg(median=("observed_tstat", "median"),
                 q25=("observed_tstat",
                      lambda x: np.percentile(x, 25)),
                 q75=("observed_tstat",
                      lambda x: np.percentile(x, 75))))


# ---------------------------------------------------------------------------
# Plotters
# ---------------------------------------------------------------------------

def plot_gaze(gaze_csv, gaze_full_csv, output_path):
    df = pd.read_csv(gaze_csv)
    if "p_value" not in df.columns:
        raise ValueError(f"{gaze_csv} is missing 'p_value'.")
    
    df = df[df["N"] < df["N"].max()]
    sig = compute_gaze_sig_rates(df, ALPHAS)
    bands = compute_gaze_cluster_bands(df)
    extent = compute_gaze_extent_bands(df)
    full = pd.read_csv(gaze_full_csv).set_index("analysis")
    n_min, n_max = int(sig["N"].min()), int(sig["N"].max())

    fig, axes = plt.subplots(3, 3, figsize=(11, 9.5), sharex=True)

    # Top row — sig rate (alpha = .05 only)
    for col, (key, title) in enumerate(GAZE_ANALYSES):
        ax = axes[0, col]
        add_stability_lines(ax)
        sub = sig[(sig["analysis"] == key) & (sig["alpha"] == 0.05)]
        sub = sub.sort_values("N")
        ax.plot(sub["N"], sub["sig_rate"],
                color=ALPHA_COLOR[0.05], lw=1.7, marker="o", ms=3.2)
        style_sigrate_axis(ax, n_min, n_max, title=title)

    # Middle row — cluster-sum stability (median ± IQR over draws)
    for col, (key, _) in enumerate(GAZE_ANALYSES):
        ax = axes[1, col]
        sub = bands[bands["analysis"] == key].sort_values("N")
        if not sub.empty:
            ax.fill_between(sub["N"], sub["q25"], sub["q75"],
                            color=TSTAT_COLOR, alpha=0.18, lw=0)
            ax.plot(sub["N"], sub["median"], color=TSTAT_COLOR, lw=1.8,
                    marker="o", ms=3.2)
        if key in full.index:
            ref = float(full.loc[key, "max_abs_cluster_sum"])
            ax.axhline(ref, color=TSTAT_FS_COLOR, lw=1.0, ls="--",
                       alpha=0.8)
        style_value_axis(ax, n_min, n_max, title=None,
                         ylabel="Max |cluster sum (t)|")

    # Bottom row — detected cluster time window (start/end ms over N)
    for col, (key, _) in enumerate(GAZE_ANALYSES):
        ax = axes[2, col]
        sub = extent[extent["analysis"] == key].sort_values("N")
        if not sub.empty:
            # Median window as filled band between start and end medians
            ax.fill_between(sub["N"], sub["start_median"], sub["end_median"],
                            color=TSTAT_COLOR, alpha=0.12, lw=0)
            # IQR shading at start (lower edge) and end (upper edge)
            ax.fill_between(sub["N"], sub["start_q25"], sub["start_q75"],
                            color=TSTAT_COLOR, alpha=0.25, lw=0)
            ax.fill_between(sub["N"], sub["end_q25"], sub["end_q75"],
                            color=TSTAT_COLOR, alpha=0.25, lw=0)
            # Median lines for both boundaries
            ax.plot(sub["N"], sub["start_median"], color=TSTAT_COLOR, lw=1.5,
                    marker="o", ms=2.8)
            ax.plot(sub["N"], sub["end_median"], color=TSTAT_COLOR, lw=1.5,
                    marker="s", ms=2.8)
        if key in full.index:
            fs_start = float(full.loc[key, "cluster_start_ms"])
            fs_end = float(full.loc[key, "cluster_end_ms"])
            if not (np.isnan(fs_start) or np.isnan(fs_end)):
                ax.axhline(fs_start, color=TSTAT_FS_COLOR, lw=1.0, ls="--",
                           alpha=0.8)
                ax.axhline(fs_end, color=TSTAT_FS_COLOR, lw=1.0, ls="--",
                           alpha=0.8)
        style_value_axis(ax, n_min, n_max, title=None, ylabel="Time (ms)")

    handles = reference_legend_handles(include_alphas=True, alpha_set=ALPHAS,
                                       include_full_sample=False,
                                       include_iqr=True)
    handles.append(Line2D([0], [0], color=TSTAT_FS_COLOR, lw=1.0, ls="--",
                          label="Full-sample reference"))
    handles.append(Line2D([0], [0], color=TSTAT_COLOR, lw=1.5,
                          marker="o", ms=4, label="Cluster start (median)"))
    handles.append(Line2D([0], [0], color=TSTAT_COLOR, lw=1.5,
                          marker="s", ms=4, label="Cluster end (median)"))
    fig.legend(handles=handles, loc="center right",
               bbox_to_anchor=(1.0, 0.5), fontsize=7.5,
               frameon=True, borderaxespad=0.5)

    fig.suptitle(
        "Gaze sensitivity — significance rate (top), "
        "cluster-sum stability (middle), "
        "detected cluster time window (bottom); "
        "median ± IQR over draws", fontsize=10.5)
    fig.tight_layout(rect=(0, 0, 0.84, 0.96))
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {output_path}")


def plot_erp(long_csv, full_csv, predictors_with_titles, output_path,
             title_prefix):
    pred_keys = [k for k, _ in predictors_with_titles]

    sig_long = _read_erp_long_for_sigrate(long_csv, pred_keys)
    sig_long = sig_long[sig_long["N"] < sig_long["N"].max()]
    sig = compute_erp_sig_rates(sig_long, ALPHAS)
    if sig.empty:
        raise RuntimeError(
            f"No matching predictors found in {long_csv}. "
            f"Looked for: {pred_keys}")

    tstat_long = _read_erp_long_for_tstat(long_csv, pred_keys)
    tstat_long = tstat_long[tstat_long["N"] < tstat_long["N"].max()]
    bands = compute_erp_tstat_bands(tstat_long)

    # Per-bin full-sample t-stat (one reference value per time bin)
    full_ref = None  # dict: predictor -> {time_bin: t}
    if full_csv and os.path.exists(full_csv):
        try:
            full_df = pd.read_csv(
                full_csv,
                usecols=["predictor", "time_bin", "observed_tstat"])
            full_ref = {}
            for pred in pred_keys:
                rows = full_df[full_df["predictor"] == pred]
                full_ref[pred] = dict(zip(rows["time_bin"],
                                          rows["observed_tstat"]))
        except Exception as e:
            print(f"  warning: could not read {full_csv}: {e}")

    n_min, n_max = int(sig["N"].min()), int(sig["N"].max())
    bins = sorted(sig["time_bin"].unique())
    colors = time_bin_colors(bins)

    fig, axes = plt.subplots(2, 3, figsize=(11, 6.5), sharex=True)

    # Top row — sig rate, one line per time bin (alpha = .05)
    sig_05 = sig[sig["alpha"] == 0.05]
    for col, (key, title) in enumerate(predictors_with_titles):
        ax = axes[0, col]
        add_stability_lines(ax)
       # add_full_sample_line(ax, n_max)
        sub_pred = sig_05[sig_05["predictor"] == key]
        for b in bins:
            sb = sub_pred[sub_pred["time_bin"] == b].sort_values("N")
            if sb.empty:
                continue
            ax.plot(sb["N"], sb["sig_rate"],
                    color=colors[b], lw=1.4, marker="o", ms=2.6)
        style_sigrate_axis(ax, n_min, n_max, title=title)

    # Bottom row — t-stat stability, one line per time bin
    # (median over draws + IQR shade)
    for col, (key, _) in enumerate(predictors_with_titles):
        ax = axes[1, col]
        #add_full_sample_line(ax, n_max)
        sub_pred = bands[bands["predictor"] == key]
        for b in bins:
            sb = sub_pred[sub_pred["time_bin"] == b].sort_values("N")
            if sb.empty:
                continue
            color = colors[b]
            ax.fill_between(sb["N"], sb["q25"], sb["q75"],
                            color=color, alpha=0.15, lw=0)
            ax.plot(sb["N"], sb["median"], color=color, lw=1.4,
                    marker="o", ms=2.6)
            # per-bin full-sample reference
            if full_ref and key in full_ref and b in full_ref[key]:
                ax.axhline(float(full_ref[key][b]),
                           color=color, lw=0.7, ls="--", alpha=0.55)
        style_value_axis(ax, n_min, n_max, title=None, ylabel="t-statistic")

    # Legend: time-bin colours first, then references
    bin_handles = time_bin_legend_handles(bins, colors)
    ref_handles = reference_legend_handles(
        include_alphas=False, include_full_sample=True, include_iqr=True)
   # if full_ref:
   #     ref_handles.append(Line2D([0], [0], color="#888888", lw=0.7, ls="--",
   #                               label="Full-sample t-stat"))
    fig.legend(handles=bin_handles + ref_handles, loc="center right",
               bbox_to_anchor=(1.0, 0.5), fontsize=7.5,
               frameon=True, borderaxespad=0.5)

    fig.suptitle(
        f"{title_prefix} — significance rate (top) and t-stat stability "
        f"(bottom, median ± IQR over draws), one line per time bin "
        f"(α = .05)", fontsize=10.5)
    fig.tight_layout(rect=(0, 0, 0.82, 0.95))
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {output_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gaze_csv", required=True)
    ap.add_argument("--gaze_full_csv", required=True)
    ap.add_argument("--frp_csv", required=True)
    ap.add_argument("--frp_full_csv", default=None)
    ap.add_argument("--language_csv", required=True)
    ap.add_argument("--language_full_csv", default=None)
    ap.add_argument("--output_dir", required=True)
    args = ap.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    print("Plotting gaze sensitivity figure …")
    plot_gaze(args.gaze_csv, args.gaze_full_csv,
              os.path.join(args.output_dir, "fig_sensitivity_gaze.pdf"))

    print("Plotting FRP sensitivity figure …")
    plot_erp(args.frp_csv, args.frp_full_csv,
             FRP_PREDICTORS,
             os.path.join(args.output_dir, "fig_sensitivity_frp.pdf"),
             title_prefix="Fixation-related potentials")

    print("Plotting language sensitivity figure …")
    plot_erp(args.language_csv, args.language_full_csv,
             LANGUAGE_PREDICTORS,
             os.path.join(args.output_dir, "fig_sensitivity_language.pdf"),
             title_prefix="Language-related potentials")

    print("Done.")


if __name__ == "__main__":
    main()