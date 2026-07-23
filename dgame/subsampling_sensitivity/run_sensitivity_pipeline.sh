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
# The PYTHON steps below depend on inputs produced UPSTREAM by the main
# DGAME pipeline:
#   * EEG : H_reconstruct_ERPs  ->  per-subject *_unfold_{FIX,N}.csv files
#                                    under {eeg_outdir}/{subject}/unfold_out/
#   * Gaze: Da_gaze_stats (with time_cluster_analysis: true)
#                            ->  time_cluster_data_{target,goal,comp}.csv
#                                under {gaze_outdir}/
#
# Usage:
#   1. Edit sensitivity_config.sh to set STUDYROOT, DGAME_OUTDIR, and
#      analysis parameters.
#   2. Run:  bash run_sensitivity_pipeline.sh
# =============================================================================

set -euo pipefail

# Locate this script's own directory so sensitivity_config.sh is always found
# regardless of the working directory from which this script is invoked.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Load user configuration (edit sensitivity_config.sh, not this file)
# shellcheck source=sensitivity_config.sh
source "${SCRIPT_DIR}/sensitivity_config.sh"

# =============================================================================
# Derived paths — based on the default pipeline output directory structure.
# If you changed eeg_dir or gaze_dir in your config.yml from their defaults
# ("eeg/" and "eyetracking/gaze_positions/"), adjust the two lines below.
# =============================================================================
SCRIPTDIR="${STUDYROOT}/dgame/subsampling_sensitivity"
CHANLOCS="${STUDYROOT}/dgame/channel_positions.csv"

EEG_OUTDIR="${DGAME_OUTDIR}/eeg"
GAZE_OUTDIR="${DGAME_OUTDIR}/eyetracking/gaze_positions"

ERP_CSV_DIR="${EEG_OUTDIR}"
GAZE_TARGET_CSV="${GAZE_OUTDIR}/time_cluster_data_target.csv"
GAZE_GOAL_CSV="${GAZE_OUTDIR}/time_cluster_data_goal.csv"
GAZE_COMP_CSV="${GAZE_OUTDIR}/time_cluster_data_comp.csv"

EEG_OUT="${DGAME_OUTDIR}/sensitivity_results"
GAZE_OUT="${DGAME_OUTDIR}/gaze_sensitivity_results"
GAZE_PLOTS="${DGAME_OUTDIR}/gaze_sensitivity_plots"
FIG_OUT="${DGAME_OUTDIR}/sensitivity_figures"

echo "================================================================"
echo " Sensitivity pipeline  |  outdir: ${DGAME_OUTDIR}"
echo " N=${N_MIN}-${N_MAX}, draws=${N_DRAWS}, seed=${SEED}"
echo "================================================================"

# -----------------------------------------------------------------------------
# STEP 1 — EEG / ERP sensitivity  (FRP + language, both modes)
#          independent of the gaze step
# -----------------------------------------------------------------------------
echo; echo ">>> [1/4] EEG sensitivity (fixations + language) ..."
"${PYTHON}" "${SCRIPTDIR}/subsampling_sensitivity.py" \
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
"${PYTHON}" "${SCRIPTDIR}/gaze_subsampling_sensitivity.py" \
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
"${PYTHON}" "${SCRIPTDIR}/gaze_sensitivity_postprocessing.py" \
    --results_dir "${GAZE_OUT}" \
    --output_dir "${GAZE_PLOTS}"

# -----------------------------------------------------------------------------
# STEP 4 — Main manuscript figures (Figs 14-16); needs STEP 1 AND STEP 2
# -----------------------------------------------------------------------------
echo; echo ">>> [4/4] Main sensitivity figures (manuscript Figs 14-16) ..."
"${PYTHON}" "${SCRIPTDIR}/plot_sensitivity_main_figures.py" \
    --gaze_csv      "${GAZE_OUT}/gaze_sensitivity_results.csv" \
    --gaze_full_csv "${GAZE_OUT}/gaze_full_sample_results.csv" \
    --frp_csv       "${EEG_OUT}/sensitivity_results_fixations.csv" \
    --frp_full_csv  "${EEG_OUT}/original_full_sample_fixations.csv" \
    --language_csv      "${EEG_OUT}/sensitivity_results_language.csv" \
    --language_full_csv "${EEG_OUT}/original_full_sample_language.csv" \
    --output_dir "${FIG_OUT}"

# -----------------------------------------------------------------------------
# STEP 5 — OSF supplement figures (one page per predictor); needs STEP 1 only
# -----------------------------------------------------------------------------
echo; echo ">>> [5/5] OSF all-effects figures (FRP + language) ..."
"${PYTHON}" "${SCRIPTDIR}/plot_sensitivity_osf_all_effects.py" \
    --frp_csv           "${EEG_OUT}/sensitivity_results_fixations.csv" \
    --frp_full_csv      "${EEG_OUT}/original_full_sample_fixations.csv" \
    --language_csv      "${EEG_OUT}/sensitivity_results_language.csv" \
    --language_full_csv "${EEG_OUT}/original_full_sample_language.csv" \
    --output_dir "${FIG_OUT}"

echo; echo "================================================================"
echo " Done. Main figures in: ${FIG_OUT}"
echo "       Gaze figures in: ${GAZE_PLOTS}"
echo "================================================================"
