#!/usr/bin/env bash
# =============================================================================
# run_sensitivity_pipeline.sh
# -----------------------------------------------------------------------------
# Master runner for the subsampling sensitivity analysis of the BRM paper
# "Combined eye tracking and electroencephalography during referential
#  selection in dyadic interaction" (Brilmayer, Georgis & Schumacher).
#
# Runs the four Python steps in the correct order and produces the three
# main sensitivity figures (manuscript Figs 14-16) plus the gaze-specific
# diagnostic figures.
#
# The PYTHON steps below depend on inputs produced UPSTREAM by the R/MATLAB
# pipeline (not run here):
#   * EEG : Unfold deconvolution + Ja_lm_permute_and_plot_fixations.r /
#           Jb_lm_permute_and_plot_language.r
#           -> per-subject reconstructed ERP CSVs + channel_positions.csv
#   * Gaze: Da_gaze_stats.r  (eyetrackingR make_time_cluster_data)
#           -> time_cluster_data_{target,goal,comp}.csv
#
# Usage:  bash run_sensitivity_pipeline.sh
# =============================================================================

set -euo pipefail

# =============================================================================
# CONFIG  ---  EDIT THESE PATHS/PARAMETERS, everything else is relative to root
# =============================================================================
STUDYROOT="./"                                   # project root (as in the R/Matlab scripts)
SCRIPTDIR="${STUDYROOT}"                          # where the .py scripts live
PYTHON="python3"                                  # or 'python' / a venv interpreter

# ---- INPUTS (produced by the R/Unfold pipeline) ----------------------------
ERP_CSV_DIR="${STUDYROOT}derivatives/erp_csvs"                    # reconstructed single-subject ERP CSVs (Unfold output)
CHANLOCS="${STUDYROOT}derivatives/channel_positions.csv"         # channel x/y positions
GAZE_TARGET_CSV="${STUDYROOT}derivatives/time_cluster_data_target.csv"
GAZE_GOAL_CSV="${STUDYROOT}derivatives/time_cluster_data_goal.csv"
GAZE_COMP_CSV="${STUDYROOT}derivatives/time_cluster_data_comp.csv"

# ---- OUTPUTS ---------------------------------------------------------------
EEG_OUT="${STUDYROOT}sensitivity_results"         # EEG (FRP + language) results
GAZE_OUT="${STUDYROOT}gaze_sensitivity_results"   # gaze results
GAZE_PLOTS="${STUDYROOT}gaze_sensitivity_plots"   # gaze-specific figures
FIG_OUT="${STUDYROOT}figures"                     # the three main manuscript figures

# ---- PARAMETERS (values used in the paper) ---------------------------------
N_MIN=15            # smallest subsample size
N_MAX=26            # full sample (N=26 row = full-sample reference)
N_DRAWS=200         # draws per N  (NOTE: script default is 10 -- the paper uses 200)
SEED=42             # RNG seed for reproducibility
EEG_PERM=2000       # permutations, EEG second-order test
GAZE_PERM=4000      # permutations, gaze cluster test
# =============================================================================

echo "================================================================"
echo " Sensitivity pipeline  |  studyroot: ${STUDYROOT}"
echo " N=${N_MIN}-${N_MAX}, draws=${N_DRAWS}, seed=${SEED}"
echo "================================================================"

# -----------------------------------------------------------------------------
# STEP 1 — EEG / ERP sensitivity  (FRP + language, both modes)
#          independent of the gaze step
# -----------------------------------------------------------------------------
echo; echo ">>> [1/4] EEG sensitivity (fixations + language) ..."
"${PYTHON}" "${SCRIPTDIR}subsampling_sensitivity.py" \
    --data_dir "${ERP_CSV_DIR}" \
    --channel_positions "${CHANLOCS}" \
    --mode both \
    --n_min "${N_MIN}" --n_max "${N_MAX}" --n_draws "${N_DRAWS}" \
    --n_permutations "${EEG_PERM}" \
    --seed "${SEED}" \
    --output_dir "${EEG_OUT}"

# -----------------------------------------------------------------------------
# STEP 2 — Gaze sensitivity  (target / goal / competitor)
#          independent of step 1 (could run in parallel)
# -----------------------------------------------------------------------------
echo; echo ">>> [2/4] Gaze sensitivity (target / goal / competitor) ..."
"${PYTHON}" "${SCRIPTDIR}gaze_subsampling_sensitivity.py" \
    --target_csv "${GAZE_TARGET_CSV}" \
    --goal_csv "${GAZE_GOAL_CSV}" \
    --comp_csv "${GAZE_COMP_CSV}" \
    --n_min "${N_MIN}" --n_max "${N_MAX}" --n_draws "${N_DRAWS}" \
    --n_permutations "${GAZE_PERM}" \
    --seed "${SEED}" \
    --output_dir "${GAZE_OUT}"

# -----------------------------------------------------------------------------
# STEP 3 — Gaze post-processing  (alpha {0.20,0.10,0.05,0.01} x stability
#          {0.80,0.90,0.95} sweep + gaze-specific figures); needs STEP 2
# -----------------------------------------------------------------------------
echo; echo ">>> [3/4] Gaze post-processing (alpha/threshold sweep + figures) ..."
"${PYTHON}" "${SCRIPTDIR}gaze_sensitivity_postprocessing.py" \
    --results_dir "${GAZE_OUT}" \
    --output_dir "${GAZE_PLOTS}"

# -----------------------------------------------------------------------------
# STEP 4 — Main manuscript figures (Figs 14-16); needs STEP 1 AND STEP 2
# -----------------------------------------------------------------------------
echo; echo ">>> [4/4] Main sensitivity figures (manuscript Figs 14-16) ..."
"${PYTHON}" "${SCRIPTDIR}plot_sensitivity_main_figures.py" \
    --gaze_csv      "${GAZE_OUT}/gaze_sensitivity_results.csv" \
    --gaze_full_csv "${GAZE_OUT}/gaze_full_sample_results.csv" \
    --frp_csv       "${EEG_OUT}/sensitivity_results_fixations.csv" \
    --frp_full_csv  "${EEG_OUT}/original_full_sample_fixations.csv" \
    --language_csv      "${EEG_OUT}/sensitivity_results_language.csv" \
    --language_full_csv "${EEG_OUT}/original_full_sample_language.csv" \
    --output_dir "${FIG_OUT}"

echo; echo "================================================================"
echo " Done. Main figures in: ${FIG_OUT}"
echo "       Gaze figures in: ${GAZE_PLOTS}"
echo "================================================================"


