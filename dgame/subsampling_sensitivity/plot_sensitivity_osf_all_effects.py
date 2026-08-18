"""
plot_sensitivity_osf_all_effects.py
=====================================

Multi-page PDF for OSF — one page per predictor.

Predictor filter (applied to both FRP and LRP):
  keep  if 'condition[T.no_conflict]' in name
           OR 'trial_time' / 'mean_target_fixation' in name
  drop  if name == 'Intercept' or 'baseline' in name
  => keeps all condition and slope effects incl. topo interactions,
     drops pure topography predictors and baseline nuisance terms.

Each page = one predictor, 2 subplots (shared x-axis):
  Top    — significance rate × N (one line per time bin, alpha = .05)
  Bottom — t-statistic median +/- IQR × N (one line per time bin;
           dashed line = full-sample reference per time bin)

Outputs:
  <output_dir>/fig_osf_all_effects_frp.pdf       (~45 pages)
  <output_dir>/fig_osf_all_effects_language.pdf   (~45 pages)

Inputs (CLI):
  --frp_csv           sensitivity_results_fixations.csv
  --frp_full_csv      original_full_sample_fixations.csv   (optional)
  --language_csv      sensitivity_results_language.csv
  --language_full_csv original_full_sample_language.csv    (optional)
  --output_dir        directory for output PDFs
"""

from __future__ import annotations

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from plot_sensitivity_main_figures import (
    ALPHAS,
    STAB_LEVELS, STAB_LINESTYLE, STAB_COLOR,
    add_stability_lines,
    style_sigrate_axis, style_value_axis,
    reference_legend_handles,
    time_bin_colors, time_bin_legend_handles,
)


# ---------------------------------------------------------------------------
# Predictor filter
# ---------------------------------------------------------------------------

def _keep_frp(pred: str) -> bool:
    if pred == "Intercept":
        return False
    if "baseline" in pred:
        return False
    return "condition[T.no_conflict]" in pred or "trial_time" in pred


def _keep_lrp(pred: str) -> bool:
    if pred == "Intercept":
        return False
    if "baseline" in pred:
        return False
    return "condition[T.no_conflict]" in pred or "mean_target_fixation" in pred


# ---------------------------------------------------------------------------
# Core plotting
# ---------------------------------------------------------------------------

def plot_all_effects(long_csv: str, full_csv: str | None,
                     keep_fn, output_path: str, title_prefix: str) -> None:

    df_sig = pd.read_csv(long_csv,
                         usecols=["predictor", "time_bin", "N",
                                  "fdr_q_value_global"])
    df_sig = df_sig[df_sig["predictor"].apply(keep_fn)].copy()
    df_sig = df_sig[df_sig["N"] < df_sig["N"].max()]

    df_t = pd.read_csv(long_csv,
                       usecols=["predictor", "time_bin", "N",
                                "observed_tstat"])
    df_t = df_t[df_t["predictor"].apply(keep_fn)].copy()
    df_t = df_t[df_t["N"] < df_t["N"].max()]

    full_ref: dict[str, dict] = {}
    if full_csv and os.path.exists(full_csv):
        try:
            fdf = pd.read_csv(full_csv,
                              usecols=["predictor", "time_bin",
                                       "observed_tstat"])
            for pred in df_sig["predictor"].unique():
                rows = fdf[fdf["predictor"] == pred]
                if not rows.empty:
                    full_ref[pred] = dict(zip(rows["time_bin"],
                                             rows["observed_tstat"]))
        except Exception as e:
            print(f"  warning: could not read {full_csv}: {e}")

    predictors = [p for p in df_sig["predictor"].unique() if keep_fn(p)]
    n_min = int(df_sig["N"].min())
    n_max = int(df_sig["N"].max())
    bins = sorted(df_sig["time_bin"].unique())
    colors = time_bin_colors(bins)

    bin_handles = time_bin_legend_handles(bins, colors)
    ref_handles = reference_legend_handles(
        include_alphas=False, include_full_sample=False, include_iqr=True)
    all_handles = bin_handles + ref_handles

    print(f"  {len(predictors)} predictors → {output_path}")

    with PdfPages(output_path) as pdf:
        for pred in predictors:
            fig, axes = plt.subplots(2, 1, figsize=(7, 5.5), sharex=True)

            # --- top: significance rate ---
            ax = axes[0]
            sub = df_sig[df_sig["predictor"] == pred].copy()
            sub["sig"] = (sub["fdr_q_value_global"] < 0.05).astype(float)
            rate = (sub.groupby(["time_bin", "N"], as_index=False)
                    .agg(sig_rate=("sig", "mean")))

            add_stability_lines(ax)
            for b in bins:
                sb = rate[rate["time_bin"] == b].sort_values("N")
                if sb.empty:
                    continue
                ax.plot(sb["N"], sb["sig_rate"],
                        color=colors[b], lw=1.4, marker="o", ms=2.6)
            style_sigrate_axis(ax, n_min, n_max)

            # --- bottom: t-stat median ± IQR ---
            ax2 = axes[1]
            subt = df_t[df_t["predictor"] == pred]
            bands = (subt.groupby(["time_bin", "N"], as_index=False)
                     .agg(median=("observed_tstat", "median"),
                          q25=("observed_tstat",
                               lambda x: np.percentile(x, 25)),
                          q75=("observed_tstat",
                               lambda x: np.percentile(x, 75))))

            for b in bins:
                sb = bands[bands["time_bin"] == b].sort_values("N")
                if sb.empty:
                    continue
                c = colors[b]
                ax2.fill_between(sb["N"], sb["q25"], sb["q75"],
                                 color=c, alpha=0.15, lw=0)
                ax2.plot(sb["N"], sb["median"],
                         color=c, lw=1.4, marker="o", ms=2.6)
                if pred in full_ref and b in full_ref[pred]:
                    ax2.axhline(float(full_ref[pred][b]),
                                color=c, lw=0.7, ls="--", alpha=0.55)
            style_value_axis(ax2, n_min, n_max)

            fig.legend(handles=all_handles, loc="center right",
                       bbox_to_anchor=(1.0, 0.5), fontsize=7,
                       frameon=True, borderaxespad=0.5)
            fig.suptitle(f"{title_prefix}\n{pred}", fontsize=9)
            fig.tight_layout(rect=(0, 0, 0.82, 0.94))
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

    print(f"  done.")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--frp_csv", required=True)
    ap.add_argument("--frp_full_csv", default=None)
    ap.add_argument("--language_csv", required=True)
    ap.add_argument("--language_full_csv", default=None)
    ap.add_argument("--output_dir", required=True)
    args = ap.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    print("Plotting OSF all-effects FRP …")
    plot_all_effects(
        args.frp_csv, args.frp_full_csv,
        _keep_frp,
        os.path.join(args.output_dir, "fig_osf_all_effects_frp.pdf"),
        title_prefix="Fixation-related potentials")

    print("Plotting OSF all-effects LRP …")
    plot_all_effects(
        args.language_csv, args.language_full_csv,
        _keep_lrp,
        os.path.join(args.output_dir, "fig_osf_all_effects_language.pdf"),
        title_prefix="Language-related potentials")

    print("Done.")


if __name__ == "__main__":
    main()
