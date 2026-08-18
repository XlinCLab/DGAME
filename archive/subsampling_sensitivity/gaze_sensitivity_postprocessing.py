"""
gaze_sensitivity_postprocessing.py
================================================================================
Post-processing and plotting for gaze_subsampling_sensitivity.py output.

Sweeps alpha thresholds {0.20, 0.10, 0.05, 0.01} and stability thresholds
{0.80, 0.90, 0.95} over gaze_sensitivity_results.csv — no need to re-run
the permutation step.

Produces four figure types:
  1. Significance-rate curves  — sig rate vs N, lines per alpha, subplots
                                  per analysis (target / competitor / goal)
  2. Cluster-sum stability     — max abs cluster sum median ± IQR vs N,
                                  with full-sample reference, subplots per analysis
  3. Cluster temporal extent   — detected cluster start/end time vs N,
                                  median ± IQR, subplots per analysis
  4. Threshold-N heatmap       — rows = analysis, cols = alpha, panels per
                                  stability level

Summary CSV with the full threshold table is also saved.

Usage:
    python gaze_sensitivity_postprocessing.py \
        --results_dir gaze_sensitivity_results \
        --output_dir gaze_sensitivity_plots
================================================================================
"""

import argparse
import os
import warnings
import subprocess
import sys


def _ensure_packages():
    required = {
        "numpy": "numpy",
        "pandas": "pandas",
        "matplotlib": "matplotlib",
        "seaborn": "seaborn",
    }
    for import_name, pip_name in required.items():
        try:
            __import__(import_name)
        except ImportError:
            print(f"Installing {pip_name}...")
            subprocess.check_call(
                [sys.executable, "-m", "pip", "install", pip_name, "-q"])


_ensure_packages()

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns

warnings.filterwarnings("ignore", category=FutureWarning)


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

ALPHAS = [0.20, 0.10, 0.05, 0.01]
STABILITY_THRESHOLDS = [0.80, 0.90, 0.95]
ANALYSES = ["target", "competitor", "goal"]

_ALPHA_PALETTE = sns.color_palette("rocket_r", n_colors=len(ALPHAS) + 1)[1:]
ALPHA_COLORS = dict(zip(ALPHAS, _ALPHA_PALETTE))
ALPHA_LABELS = {a: f"α = {a}" for a in ALPHAS}

STAB_LINESTYLES = {0.80: (0, (5, 5)), 0.90: (0, (3, 2)), 0.95: "solid"}

ANALYSIS_COLORS = {
    "target": "#4878d0",
    "competitor": "#e24a33",
    "goal": "#6acc65",
}
ANALYSIS_LABELS = {
    "target": "Target AOI",
    "competitor": "Competitor AOI",
    "goal": "Goal AOI",
}


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def load_sensitivity_results(results_dir):
    path = os.path.join(results_dir, "gaze_sensitivity_results.csv")
    if not os.path.exists(path):
        raise FileNotFoundError(f"Not found: {path}")
    df = pd.read_csv(path)
    print(f"  Loaded {len(df):,} rows from {os.path.basename(path)}")
    # Keep only the three analyses of interest (guard against stale rows)
    df = df[df["analysis"].isin(ANALYSES)].reset_index(drop=True)
    return df


def load_full_results(results_dir):
    path = os.path.join(results_dir, "gaze_full_sample_results.csv")
    if not os.path.exists(path):
        print("  Full-sample reference not found, skipping reference lines")
        return None
    return pd.read_csv(path).set_index("analysis")


# ---------------------------------------------------------------------------
# Core computations
# ---------------------------------------------------------------------------

def compute_significance_rates(df):
    """Significance rate per (analysis, N) for each alpha (re-threshold p_value)."""
    frames = []
    for alpha in ALPHAS:
        rate = (
            df.assign(significant=(df["p_value"] < alpha).astype(float))
            .groupby(["analysis", "N"], as_index=False)
            .agg(significance_rate=("significant", "mean"),
                 n_draws=("significant", "size"))
        )
        rate["alpha"] = alpha
        frames.append(rate)
    return pd.concat(frames, ignore_index=True)


def compute_cluster_sum_bands(df):
    """max_abs_cluster_sum: median, IQR, mean, SD per (analysis, N)."""
    return (
        df.groupby(["analysis", "N"], as_index=False)
        .agg(
            cs_median=("max_abs_cluster_sum", "median"),
            cs_mean=("max_abs_cluster_sum", "mean"),
            cs_sd=("max_abs_cluster_sum", "std"),
            cs_q25=("max_abs_cluster_sum", lambda x: np.percentile(x, 25)),
            cs_q75=("max_abs_cluster_sum", lambda x: np.percentile(x, 75)),
            n_draws=("max_abs_cluster_sum", "size"),
        )
    )


def compute_extent_bands(df):
    """Cluster start/end ms: median and IQR per (analysis, N), draws with clusters only."""
    sig = df[df["max_abs_cluster_sum"] > 0].copy()
    if sig.empty:
        return pd.DataFrame()
    return (
        sig.groupby(["analysis", "N"], as_index=False)
        .agg(
            start_median=("cluster_start_ms", "median"),
            start_q25=("cluster_start_ms", lambda x: np.percentile(x, 25)),
            start_q75=("cluster_start_ms", lambda x: np.percentile(x, 75)),
            end_median=("cluster_end_ms", "median"),
            end_q25=("cluster_end_ms", lambda x: np.percentile(x, 25)),
            end_q75=("cluster_end_ms", lambda x: np.percentile(x, 75)),
            n_with_cluster=("cluster_start_ms", "count"),
        )
    )


def compute_threshold_n(sig_rates):
    """Minimum N to reach stability_threshold significance rate per
    (analysis, alpha, stability_threshold)."""
    rows = []
    for stab in STABILITY_THRESHOLDS:
        for (analysis, alpha), grp in sig_rates.groupby(["analysis", "alpha"]):
            passed = grp[grp["significance_rate"] >= stab]
            thr = int(passed["N"].min()) if not passed.empty else None
            rows.append({
                "analysis": analysis,
                "alpha": alpha,
                "stability_threshold": stab,
                "threshold_N": thr,
            })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Figure 1 — Significance-rate curves
# ---------------------------------------------------------------------------

def plot_sig_rate_curves(sig_rates, output_dir):
    analyses = [a for a in ANALYSES if a in sig_rates["analysis"].unique()]
    n = len(analyses)

    fig, axes = plt.subplots(1, n, figsize=(4.2 * n, 3.8), squeeze=False, sharey=True)
    fig.suptitle("Gaze sensitivity — significance rate", fontsize=11, y=1.02)

    for col_i, analysis in enumerate(analyses):
        ax = axes[0][col_i]
        adat = sig_rates[sig_rates["analysis"] == analysis]

        for alpha in ALPHAS:
            d = adat[adat["alpha"] == alpha].sort_values("N")
            ax.plot(d["N"], d["significance_rate"],
                    color=ALPHA_COLORS[alpha], marker="o", ms=3.5,
                    lw=1.5, label=ALPHA_LABELS[alpha])

        for stab, ls in STAB_LINESTYLES.items():
            ax.axhline(stab, color="#888", lw=0.8, ls=ls, alpha=0.75,
                       label=f"{int(stab*100)}% stability" if col_i == 0 else "")
            # Mark full sample (N_max) with a vertical reference line
            n_max = int(adat["N"].max())
            ax.axvline(n_max, color="#c43c3c", lw=1.0, ls=(0, (4, 2)),
                   alpha=0.7, zorder=1)

        ax.set_title(ANALYSIS_LABELS.get(analysis, analysis), fontsize=10, pad=4)
        ax.set_ylim(-0.05, 1.05)
        ax.yaxis.set_major_locator(mticker.MultipleLocator(0.25))
        ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True, nbins=5))
        ax.tick_params(labelsize=8)
        ax.set_xlabel("Sample size N", fontsize=9)
        if col_i == 0:
            ax.set_ylabel("Significance rate", fontsize=9)

    # Legend on last axis
    handles_a = [
        plt.Line2D([0], [0], color=ALPHA_COLORS[a], marker="o", ms=4,
                   lw=1.5, label=ALPHA_LABELS[a])
        for a in ALPHAS
    ]
    handles_s = [
        plt.Line2D([0], [0], color="#888", lw=0.8, ls=ls,
                   label=f"{int(s*100)}% stability")
        for s, ls in STAB_LINESTYLES.items()
    ]
    handles_full = [
        plt.Line2D([0], [0], color="#c43c3c", lw=1.0, ls=(0, (4, 2)),
                   label="Full sample (N_max)")
    ]
    axes[0][-1].legend(handles=handles_a + handles_s + handles_full,
                       fontsize=7.5, frameon=True,
                       loc="lower right")
    
    fig.tight_layout()
    fname = os.path.join(output_dir, "fig1_sig_rate_gaze.pdf")
    fig.savefig(fname, bbox_inches="tight")
    plt.close(fig)
    print(f"    saved {os.path.basename(fname)}")


# ---------------------------------------------------------------------------
# Figure 2 — Cluster-sum stability bands
# ---------------------------------------------------------------------------

def plot_cluster_sum_bands(cs_bands, full_results, output_dir):
    analyses = [a for a in ANALYSES if a in cs_bands["analysis"].unique()]
    n = len(analyses)

    fig, axes = plt.subplots(1, n, figsize=(4.2 * n, 3.8), squeeze=False)
    fig.suptitle("Gaze sensitivity — cluster sum stability", fontsize=11, y=1.02)

    for col_i, analysis in enumerate(analyses):
        ax = axes[0][col_i]
        d = cs_bands[cs_bands["analysis"] == analysis].sort_values("N")
        color = ANALYSIS_COLORS.get(analysis, "#4878d0")

        ax.fill_between(d["N"], d["cs_q25"], d["cs_q75"],
                        color=color, alpha=0.25)
        ax.plot(d["N"], d["cs_median"],
                color=color, lw=1.6, marker="o", ms=3.5,
                label="Median ± IQR")

        if full_results is not None and analysis in full_results.index:
            ref = full_results.loc[analysis, "max_abs_cluster_sum"]
            if pd.notna(ref):
                ax.axhline(float(ref), color="#e24a33", lw=1.2, ls="--",
                           label="Full sample")

        ax.axhline(0, color="#bbb", lw=0.6)
        ax.set_title(ANALYSIS_LABELS.get(analysis, analysis), fontsize=10, pad=4)
        ax.set_xlabel("Sample size N", fontsize=9)
        ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True, nbins=5))
        ax.tick_params(labelsize=8)
        if col_i == 0:
            ax.set_ylabel("Max abs cluster sum (t)", fontsize=9)

        handles = [
            plt.Line2D([0], [0], color=color, lw=1.6, marker="o", ms=4,
                       label="Median ± IQR (subsamples)"),
            plt.Line2D([0], [0], color="#e24a33", lw=1.2, ls="--",
                       label="Full sample"),
        ]
        ax.legend(handles=handles, fontsize=7.5, frameon=True, loc="lower right")

    fig.tight_layout()
    fname = os.path.join(output_dir, "fig2_cluster_sum_gaze.pdf")
    fig.savefig(fname, bbox_inches="tight")
    plt.close(fig)
    print(f"    saved {os.path.basename(fname)}")


# ---------------------------------------------------------------------------
# Figure 3 — Cluster temporal extent
# ---------------------------------------------------------------------------

def plot_cluster_extent(extent_bands, full_results, output_dir):
    if extent_bands.empty:
        print("    No cluster extent data, skipping Fig 3")
        return

    analyses = [a for a in ANALYSES if a in extent_bands["analysis"].unique()]
    n = len(analyses)

    fig, axes = plt.subplots(1, n, figsize=(4.2 * n, 3.8), squeeze=False)
    fig.suptitle("Gaze sensitivity — detected cluster time window", fontsize=11, y=1.02)

    for col_i, analysis in enumerate(analyses):
        ax = axes[0][col_i]
        d = extent_bands[extent_bands["analysis"] == analysis].sort_values("N")
        color = ANALYSIS_COLORS.get(analysis, "#4878d0")

        if d.empty:
            ax.set_title(ANALYSIS_LABELS.get(analysis, analysis), fontsize=10)
            ax.text(0.5, 0.5, "No clusters found", transform=ax.transAxes,
                    ha="center", va="center", fontsize=9, color="#888")
            continue

        # Shaded band: median start to median end
        ax.fill_between(d["N"], d["start_median"], d["end_median"],
                        color=color, alpha=0.20, label="Median window")
        # IQR for start and end separately
        ax.fill_between(d["N"], d["start_q25"], d["start_q75"],
                        color=color, alpha=0.35)
        ax.fill_between(d["N"], d["end_q25"], d["end_q75"],
                        color=color, alpha=0.35)
        ax.plot(d["N"], d["start_median"], color=color, lw=1.4,
                marker="o", ms=3, ls="-")
        ax.plot(d["N"], d["end_median"], color=color, lw=1.4,
                marker="s", ms=3, ls="-")

        # Full-sample reference window
        if full_results is not None and analysis in full_results.index:
            row = full_results.loc[analysis]
            fs_start = row.get("cluster_start_ms", np.nan)
            fs_end = row.get("cluster_end_ms", np.nan)
            if pd.notna(fs_start) and pd.notna(fs_end):
                ax.axhline(float(fs_start), color="#e24a33", lw=1.0, ls="--",
                           label=f"Full sample: {int(fs_start)}–{int(fs_end)} ms")
                ax.axhline(float(fs_end), color="#e24a33", lw=1.0, ls="--")

        ax.axhline(0, color="#bbb", lw=0.6)
        ax.set_title(ANALYSIS_LABELS.get(analysis, analysis), fontsize=10, pad=4)
        ax.set_xlabel("Sample size N", fontsize=9)
        ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True, nbins=5))
        ax.tick_params(labelsize=8)
        if col_i == 0:
            ax.set_ylabel("Time (ms)", fontsize=9)
        ax.legend(fontsize=7.5, frameon=True, loc="best")

    fig.tight_layout()
    fname = os.path.join(output_dir, "fig3_cluster_extent_gaze.pdf")
    fig.savefig(fname, bbox_inches="tight")
    plt.close(fig)
    print(f"    saved {os.path.basename(fname)}")


# ---------------------------------------------------------------------------
# Figure 4 — Threshold-N heatmap
# ---------------------------------------------------------------------------

def plot_threshold_heatmap(thresh_df, output_dir):
    n_stab = len(STABILITY_THRESHOLDS)
    analyses = [a for a in ANALYSES if a in thresh_df["analysis"].unique()]

    fig, axes = plt.subplots(
        1, n_stab,
        figsize=(3.8 * n_stab, 0.55 * len(analyses) + 2.2),
        squeeze=False,
    )
    fig.suptitle("Gaze sensitivity — minimum N for stability",
                 fontsize=11, y=1.02)

    all_Ns = thresh_df["threshold_N"].dropna()
    vmin = int(all_Ns.min()) if not all_Ns.empty else 0
    vmax = int(all_Ns.max()) if not all_Ns.empty else 26
    cmap = sns.color_palette("YlOrRd", as_cmap=True)

    for col_i, stab in enumerate(STABILITY_THRESHOLDS):
        ax = axes[0][col_i]
        panel = thresh_df[thresh_df["stability_threshold"] == stab]

        pivot = (
            panel.pivot_table(
                index="analysis", columns="alpha",
                values="threshold_N", aggfunc="mean")
            .reindex(analyses)
            .reindex(columns=ALPHAS)
        )
        annot = pivot.map(lambda v: str(int(v)) if pd.notna(v) else "—")

        sns.heatmap(
            pivot, ax=ax,
            cmap=cmap, vmin=vmin, vmax=vmax,
            annot=annot, fmt="s", annot_kws={"size": 9},
            linewidths=0.4, linecolor="#ddd",
            cbar=(col_i == n_stab - 1),
            cbar_kws={"label": "Threshold N", "shrink": 0.8},
            yticklabels=[ANALYSIS_LABELS.get(a, a) for a in analyses],
        )
        ax.set_title(f"{int(stab*100)}% stability", fontsize=10, pad=6)
        ax.set_xlabel("α threshold", fontsize=8)
        ax.set_ylabel("Analysis" if col_i == 0 else "", fontsize=8)
        ax.set_xticklabels([str(a) for a in ALPHAS], fontsize=8)
        ax.set_yticklabels(ax.get_yticklabels(), fontsize=9, rotation=0)

    fig.tight_layout()
    fname = os.path.join(output_dir, "fig4_threshold_heatmap_gaze.pdf")
    fig.savefig(fname, bbox_inches="tight")
    plt.close(fig)
    print(f"    saved {os.path.basename(fname)}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run_postprocessing(results_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    print(f"\n{'='*60}")
    print(f"  Gaze sensitivity post-processing")
    print(f"{'='*60}")

    df = load_sensitivity_results(results_dir)
    full_results = load_full_results(results_dir)

    print("  Computing significance rates …")
    sig_rates = compute_significance_rates(df)

    print("  Computing cluster sum bands …")
    cs_bands = compute_cluster_sum_bands(df)

    print("  Computing cluster extent bands …")
    extent_bands = compute_extent_bands(df)

    print("  Computing threshold-N table …")
    thresh_df = compute_threshold_n(sig_rates)

    # Save summary CSVs
    sig_rates.to_csv(os.path.join(output_dir, "gaze_sig_rates.csv"), index=False)
    thresh_df.to_csv(os.path.join(output_dir, "gaze_threshold_n.csv"), index=False)
    cs_bands.to_csv(os.path.join(output_dir, "gaze_cluster_sum_bands.csv"), index=False)
    if not extent_bands.empty:
        extent_bands.to_csv(os.path.join(output_dir, "gaze_cluster_extent_bands.csv"), index=False)
    print(f"  Saved summary CSVs to {output_dir}/")

    print("  Plotting Fig 1: significance-rate curves …")
    plot_sig_rate_curves(sig_rates, output_dir)

    print("  Plotting Fig 2: cluster sum stability …")
    plot_cluster_sum_bands(cs_bands, full_results, output_dir)

    print("  Plotting Fig 3: cluster temporal extent …")
    plot_cluster_extent(extent_bands, full_results, output_dir)

    print("  Plotting Fig 4: threshold-N heatmap …")
    plot_threshold_heatmap(thresh_df, output_dir)

    print(f"\n  Done. All output in {output_dir}/")


def main():
    parser = argparse.ArgumentParser(
        description="Post-process gaze subsampling sensitivity results")
    parser.add_argument("--results_dir", default="gaze_sensitivity_results",
                        help="Directory written by gaze_subsampling_sensitivity.py")
    parser.add_argument("--output_dir", default="gaze_sensitivity_plots")
    args = parser.parse_args()

    run_postprocessing(
        results_dir=args.results_dir,
        output_dir=args.output_dir,
    )


if __name__ == "__main__":
    main()
