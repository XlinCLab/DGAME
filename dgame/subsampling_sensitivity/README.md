# Subsampling sensitivity analysis — code and run order

Code accompanying the subsampling sensitivity analysis in *Combined eye tracking
and electroencephalography during referential selection in dyadic interaction*
(Brilmayer, Georgis & Schumacher; manuscript BR-Org-25-839).

The analysis evaluates how reliably the reported effects are detected as a
function of sample size, by repeatedly drawing subsamples of N participants and
re-running the original permutation tests on each draw. The Python permutation
code reproduces the original R analyses exactly (Pearson r = 1.0 between Python
and R t-statistics; max absolute deviation 8 × 10⁻¹²).

## Prerequisites

Run the main DGAME pipeline first; it produces the inputs for the Python steps
below:

| Upstream step                                                   | Produces                                              |
|----------------------------------------------------------------|-------------------------------------------------------|
| `H_reconstruct_ERPs` (MATLAB/Unfold) + `J_lm_permute_and_plot` | per-subject reconstructed ERP CSVs (`*_unfold_{FIX,N}.csv`) in each subject's `unfold_out/` directory |
| `Da_gaze_stats` with `time_cluster_analysis: true` (see note)  | `time_cluster_data_{target,goal,comp}.csv` in the gaze output directory |

> **Note:** The `time_cluster_analysis` sub-step of `Da_gaze_stats` is **disabled by default**. To enable it, add the following to your `config.yml`:
> ```yml
> analysis:
>   steps:
>     Da_gaze_stats:
>       time_cluster_analysis: true
> ```

Python dependencies: `numpy`, `pandas`, `scipy`, `patsy`, `matplotlib`,
`statsmodels`.

## Run order

The whole pipeline is wrapped in **`run_sensitivity_pipeline.sh`**. Before running it, edit **`sensitivity_config.sh`** to set `STUDYROOT` (path to the DGAME project root), `DGAME_OUTDIR` (the `outdir` value from your `config.yml`), and optionally the analysis parameters. All input and output paths are then derived automatically. Run with:
```bash
bash run_sensitivity_pipeline.sh
```
Steps individually:

| # | Script | Input | Output | Depends on |
|---|--------|-------|--------|------------|
| 1 | `subsampling_sensitivity.py --mode both` | reconstructed ERP CSVs + `channel_positions.csv` | `sensitivity_results_{fixations,language}.csv`, `original_full_sample_{fixations,language}.csv`, `stability_summary_*`, `sample_size_threshold_*` | — |
| 2 | `gaze_subsampling_sensitivity.py` | `time_cluster_data_{target,goal,comp}.csv` | `gaze_sensitivity_results.csv`, `gaze_full_sample_results.csv`, `gaze_stability_summary.csv`, `gaze_sample_size_threshold.csv` | — |
| 3 | `gaze_sensitivity_postprocessing.py` | output of step 2 | gaze figures (`fig1–4_*_gaze.pdf`) + α/threshold sweep tables | 2 |
| 4 | `plot_sensitivity_main_figures.py` | output of steps 1 **and** 2 | `fig_sensitivity_{gaze,frp,language}.pdf` (manuscript Figs 14–16) | 1 + 2 |
| 5 | `plot_sensitivity_osf_all_effects.py` | output of step 1 | `fig_osf_all_effects_{frp,language}.pdf` (OSF supplement, one page per predictor) | 1 |

Steps 1 and 2 are independent and may run in either order or in parallel.
Steps 3 and 4 only post-process / plot; they do not re-run any permutation test.

The α-level sweep (0.01 / 0.05 / 0.10 / 0.20) and the stability thresholds
(80 % / 90 % / 95 %) reported in the paper are computed in the post-processing /
plotting steps (3 and 4) from the per-draw significance results — the
permutation engines (steps 1–2) write one significance decision per draw.

Note on the asymmetry between the two tracks: there is a dedicated
post-processing script for the gaze analysis (step 3) but not for the ERP
analysis. This is a difference in how the work is split between scripts, not in
what is computed. For the ERPs, the stability metrics are produced by the engine
itself (`compute_stability_metrics` in `subsampling_sensitivity.py`, written to
`stability_summary_*`), and the α-level sweep is performed on the fly in the
shared plotting step (4) via `fdr_q_value_global < alpha`. The gaze cluster test
additionally yields quantities that have no ERP counterpart (signed cluster sum,
detected cluster start/end in ms), so the gaze sweep and its extra diagnostic
figures were collected in a separate step (3).

## Parameters used in the paper

| Parameter | Value | Note |
|-----------|-------|------|
| Sample sizes N | 15 – 26 | N = 26 is the full sample (reference point) |
| Draws per N | **200** | the script **default is 10**; pass `--n_draws 200` |
| Permutations | 2000 (EEG) / 4000 (gaze) | `--n_permutations` |
| RNG seed | 42 | `--seed` |


