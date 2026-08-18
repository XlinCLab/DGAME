#!/usr/bin/env bash
# =============================================================================
# sensitivity_config.sh
# -----------------------------------------------------------------------------
# Configuration for the subsampling sensitivity pipeline.
# Edit this file, then run:  bash run_sensitivity_pipeline.sh
# =============================================================================

# Path to the DGAME project root (the directory containing dgame/, config/, etc.)
STUDYROOT="/absolute/path/to/dgame2py"

# Path to the experiment pipeline output directory
# (the value of 'outdir' in your experiment config.yml)
DGAME_OUTDIR="/absolute/path/to/your/experiment/output"

# Python interpreter — use the same environment as the main pipeline
PYTHON="${STUDYROOT}/venv/bin/python"

# ---------------------------------------------------------------------------
# Analysis parameters (values used in the published paper)
# ---------------------------------------------------------------------------
N_MIN=15        # smallest subsample size
N_MAX=26        # full sample size (reference point)
N_DRAWS=200     # draws per N  (script default is 10; paper used 200)
SEED=42         # RNG seed for reproducibility
EEG_PERM=2000   # permutations for EEG second-order test
GAZE_PERM=4000  # permutations for gaze cluster test
