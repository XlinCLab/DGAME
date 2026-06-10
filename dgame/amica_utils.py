"""
amica_utils.py — run the AMICA binary from Python and inject the results into an
MNE ICA object that is compatible with mne-icalabel's label_components().

The AMICA binary interface (documented in runamica15.m) is entirely file-based:
  1. EEG data is written as float32 in column-major (Fortran) order.
  2. A plain-text `input.param` file lists all parameters.
  3. The binary is called: amica15{mac|c} <outdir>/input.param
  4. Output binary files (W, S, mean, …) are read back as float64.

MNE ICA attribute mapping from AMICA outputs (W, S, mean):
  pre_whitener_          = ones((n_chan, 1))         — cancels MNE's z-score step
  pca_mean_              = mean_data                 — AMICA's channel means
  pca_components_        = S[0:n_pcs, :]             — sphering / whitening rows
  pca_explained_variance_= ones(n_pcs)               — no-op (not used in inference)
  unmixing_matrix_       = W[:, :, 0]                — ICA unmixing (post-sphering)
  mixing_matrix_         = pinv(S_pcs).T @ pinv(W@S_pcs)  — for correct topographies

This gives:
  get_sources()  → W @ S_pcs @ (data − mean)        ✓ matches EEGLAB convention
  get_components() → pinv(W @ S_pcs)                ✓ correct scalp topographies
"""

import logging
import os
import platform
import subprocess
import tempfile
from typing import Optional

import mne
import numpy as np

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Binary selection
# ---------------------------------------------------------------------------

def _amica_binary(amica_plugin_dir: str) -> str:
    """Return the path to the appropriate AMICA binary for this platform."""
    system = platform.system()
    if system == "Darwin":
        name = "amica15mac"
    elif system == "Linux":
        name = "amica15c"
    elif system == "Windows":
        name = "amica15mkl.exe"
    else:
        raise RuntimeError(f"Unsupported platform for AMICA binary: {system}")
    path = os.path.join(amica_plugin_dir, name)
    if not os.path.exists(path):
        raise FileNotFoundError(f"AMICA binary not found: {path}")
    os.chmod(path, 0o755)
    return path


# ---------------------------------------------------------------------------
# Data / param file I/O
# ---------------------------------------------------------------------------

def _write_data_file(data: np.ndarray, path: str) -> None:
    """Write [n_chan × n_timepoints] as float32 in Fortran/column-major order.

    This matches MATLAB's floatwrite() convention used by runamica15.m.
    Column-major layout means: all channels at t=0, then all channels at t=1, …
    """
    data.astype(np.float32).flatten("F").tofile(path)


def _write_param_file(
    param_path: str,
    data_path: str,
    out_dir: str,
    n_chan: int,
    n_timepoints: int,
    n_pcs: int,
    max_iter: int = 2000,
    max_threads: int = 4,
    do_reject: int = 1,
    num_rej: int = 10,
    rej_sig: float = 2.5,
    rej_int: int = 1,
    rej_start: int = 2,
) -> None:
    """Write the AMICA input.param file matching runamica15.m defaults.

    All parameters and their values are taken directly from runamica15.m so
    that the Python run is bit-for-bit equivalent to the MATLAB pipeline.
    Key points:
      - do_opt_block 0  : fixed block size (MATLAB default); adaptive blocks
                          (do_opt_block 1) triggers different early-stopping
                          behaviour and causes premature convergence.
      - do_newton 1     : Newton-step acceleration from iter newt_start=50;
                          keeps the optimisation stable for the full max_iter.
      - use_min_dll / use_grad_norm : explicit log-likelihood and gradient-norm
                          convergence guards (both are tight but in practice
                          not triggered before max_iter in the MATLAB runs).
    """
    lines = [
        f"files {data_path}",
        f"outdir {out_dir}{os.sep}",
        # --- block / iteration settings ---
        "block_size 128",
        "do_opt_block 0",         # MATLAB default: fixed block size
        "blk_min 256",
        "blk_step 256",
        "blk_max 1024",
        "num_models 1",
        f"max_threads {max_threads}",
        f"max_iter {max_iter}",
        "num_samples 1",
        f"data_dim {n_chan}",
        f"field_dim {n_timepoints}",
        "field_blocksize 1",
        # --- history / sharing (off by default in runamica15.m) ---
        "do_history 0",
        "histstep 10",
        "share_comps 0",
        "share_start 100",
        "comp_thresh 0.990000",
        "share_iter 100",
        # --- learning rate schedule ---
        "lrate 0.050000",
        "minlrate 1.000000e-08",
        "mineig 1.000000e-12",
        "lratefact 0.500000",
        # --- rho (shape) update ---
        "rholrate 0.050000",
        "rho0 1.500000",
        "minrho 1.000000",
        "maxrho 2.000000",
        "rholratefact 0.500000",
        # --- kurtosis update ---
        "kurt_start 3",
        "num_kurt 5",
        "kurt_int 1",
        # --- Newton update (enabled in MATLAB; stabilises long runs) ---
        "do_newton 1",
        "newt_start 50",
        "newt_ramp 10",
        "newtrate 1.000000",
        # --- data rejection ---
        f"do_reject {do_reject}",
        f"numrej {num_rej}",
        f"rejsig {rej_sig:.6f}",
        f"rejstart {rej_start}",
        f"rejint {rej_int}",
        # --- convergence guards ---
        "decwindow 1",
        "max_decs 3",
        "use_min_dll 1",
        "min_dll 1.000000e-09",
        "use_grad_norm 1",
        "min_grad_norm 1.000000e-07",
        # --- initialisation flags (all off) ---
        "fix_init 0",
        "load_rej 0",
        "load_W 0",
        "load_c 0",
        "load_gm 0",
        "load_alpha 0",
        "load_mu 0",
        "load_beta 0",
        "load_rho 0",
        "load_comp_list 0",
        # --- update flags (all on) ---
        "update_A 1",
        "update_c 1",
        "update_gm 1",
        "update_alpha 1",
        "update_mu 1",
        "update_beta 1",
        # --- sigma bounds ---
        "invsigmax 100.000000",
        "invsigmin 0.000000",
        "do_rho 1",
        # --- PCA / sphere / scaling ---
        f"pcakeep {n_pcs}",
        "pcadb 30.000000",
        "do_mean 1",
        "do_sphere 1",
        "doPCA 1",
        "doscaling 1",
        "scalestep 1",
        # --- PDF type (0 = extended infomax generalized Gaussian) ---
        "pdftype 0",
        "num_mix_comps 3",
        # --- output ---
        "byte_size 4",
        "writestep 20",
        "write_nd 0",
        "write_LLt 1",
    ]
    with open(param_path, "w") as fid:
        fid.write("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# Output reading  (mirrors loadmodout15.m)
# ---------------------------------------------------------------------------

def _read_amica_output(out_dir: str) -> dict:
    """Read the binary output files written by the AMICA binary.

    Returns a dict with keys: W, S, mean, n_pcs, n_chan, LL
    """
    results = {}

    # --- model mixture weights (gm) → num_models ---
    gm_path = os.path.join(out_dir, "gm")
    if os.path.exists(gm_path):
        gm = np.fromfile(gm_path, dtype=np.float64)
        num_models = len(gm)
    else:
        num_models = 1

    # --- W: unmixing matrix [n_pcs × n_pcs × num_models] (Fortran order) ---
    W_path = os.path.join(out_dir, "W")
    if not os.path.exists(W_path):
        raise FileNotFoundError(f"AMICA output W not found in {out_dir}")
    W_flat = np.fromfile(W_path, dtype=np.float64)
    n_pcs = round((len(W_flat) / num_models) ** 0.5)
    W = W_flat.reshape(n_pcs, n_pcs, num_models, order="F")
    results["W"] = W
    results["n_pcs"] = n_pcs

    # --- mean: channel means [n_chan] ---
    mean_path = os.path.join(out_dir, "mean")
    if os.path.exists(mean_path):
        mean_data = np.fromfile(mean_path, dtype=np.float64)
        n_chan = len(mean_data)
    else:
        logger.warning("AMICA output 'mean' not found; setting to zero.")
        mean_data = None
        n_chan = n_pcs  # fallback (pcakeep == n_chan case)
    results["mean"] = mean_data
    results["n_chan"] = n_chan

    # --- S: sphering matrix [n_chan × n_chan] (Fortran order) ---
    S_path = os.path.join(out_dir, "S")
    if os.path.exists(S_path):
        S = np.fromfile(S_path, dtype=np.float64).reshape(n_chan, n_chan, order="F")
    else:
        logger.warning("AMICA output 'S' not found; using identity.")
        S = np.eye(n_chan)
    results["S"] = S

    # --- LL: log-likelihood per iteration ---
    LL_path = os.path.join(out_dir, "LL")
    if os.path.exists(LL_path):
        results["LL"] = np.fromfile(LL_path, dtype=np.float64)
    else:
        results["LL"] = np.array([])

    return results


# ---------------------------------------------------------------------------
# MNE ICA construction from AMICA outputs
# ---------------------------------------------------------------------------

def _build_mne_ica_from_amica(
    amica_results: dict,
    raw: mne.io.BaseRaw,
    n_components: int,
) -> mne.preprocessing.ICA:
    """Construct an MNE ICA object populated with AMICA's decomposition.

    MNE's _transform() computes:
        unmixing_matrix_ @ pca_components_ @ (data / pre_whitener_ − pca_mean_)

    We set:
        pre_whitener_           = ones((n_chan, 1))   → cancels z-score step
        pca_mean_               = AMICA mean          → centres the data
        pca_components_         = S[0:n_pcs, :]       → sphering rows
        unmixing_matrix_        = W[:, :, 0]          → ICA unmixing

    This gives: W @ S_pcs @ (data − mean)  ✓

    For get_components(), which returns (mixing_matrix_.T @ pca_components_).T:
    We need mixing_matrix_ such that the result equals pinv(W @ S_pcs).

    Derivation:
        (M.T @ S_pcs).T = A  where A = pinv(W @ S_pcs), shape [n_chan × n_pcs]
        ⟹  M.T = A.T @ pinv(S_pcs)
        ⟹  M   = pinv(S_pcs).T @ A
    """
    W = amica_results["W"][:, :, 0]      # [n_pcs × n_pcs]
    S = amica_results["S"]               # [n_chan × n_chan]
    mean_data = amica_results["mean"]    # [n_chan]
    n_chan = amica_results["n_chan"]
    n_pcs = amica_results["n_pcs"]

    S_pcs = S[:n_pcs, :]                 # [n_pcs × n_chan]

    if mean_data is None:
        mean_data = np.zeros(n_chan)

    # Mixing matrix (scalp topographies):  A = pinv(W @ S_pcs)
    A = np.linalg.pinv(W @ S_pcs)       # [n_chan × n_pcs]

    # Construct ica object (unfitted — we set all attributes manually)
    ica = mne.preprocessing.ICA(
        n_components=n_components,
        method="infomax",
        fit_params={"extended": True},
        random_state=97,
    )

    # Populate attributes that MNE uses during inference
    ica.n_components_ = n_pcs
    ica.ch_names = list(raw.ch_names)
    ica.info = mne.pick_info(raw.info, mne.pick_types(raw.info, eeg=True, exclude=[]))

    # Pre-whitener: set to ones so MNE's z-score step is a no-op
    ica.pre_whitener_ = np.ones((n_chan, 1))

    # PCA attributes
    ica.pca_mean_ = mean_data                    # [n_chan]
    ica.pca_components_ = S_pcs                  # [n_pcs × n_chan]
    ica.pca_explained_variance_ = np.ones(n_pcs) # not used in inference; safe default

    # ICA unmixing
    ica.unmixing_matrix_ = W                     # [n_pcs × n_pcs]

    # Mixing matrix derived so get_components() returns correct topographies
    S_pcs_pinv = np.linalg.pinv(S_pcs)          # [n_chan × n_pcs]
    ica.mixing_matrix_ = S_pcs_pinv.T @ A        # [n_pcs × n_pcs]

    # Required for apply_ica / exclude list handling
    ica.exclude = []
    ica.labels_ = {}
    ica.n_iter_ = len(amica_results.get("LL", []))

    # Mark as fitted on raw data so ica.save() / apply_ica() don't raise
    ica.current_fit = "raw"
    ica.reject_ = None

    return ica


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def run_amica(
    raw: mne.io.BaseRaw,
    n_pcs: int,
    amica_plugin_dir: str,
    work_dir: Optional[str] = None,
    max_iter: int = 2000,
    max_threads: int = 4,
    log_prefix: str = "",
) -> mne.preprocessing.ICA:
    """Run the AMICA binary on raw EEG data and return a fitted MNE ICA object.

    Parameters
    ----------
    raw : mne.io.BaseRaw
        Preprocessed, average-referenced EEG data (all bad channels already dropped).
    n_pcs : int
        Number of PCA components to retain (= data rank; matches MATLAB's pcakeep).
    amica_plugin_dir : str
        Directory containing the AMICA binary (e.g. …/eeglab2025.1.0/plugins/amica).
    work_dir : str, optional
        Directory for temporary data / output files.  A subdirectory ``amica_run``
        is created inside it.  Defaults to a system temp directory.
    max_iter : int
        Maximum AMICA iterations (default 2000, matches MATLAB pipeline).
    max_threads : int
        Number of threads for the AMICA binary.
    log_prefix : str
        Optional string prepended to log messages (e.g. "Subject 02: ").

    Returns
    -------
    mne.preprocessing.ICA
        ICA object with AMICA weights loaded, compatible with label_components().
    """
    binary = _amica_binary(amica_plugin_dir)
    data = raw.get_data()          # [n_chan × n_timepoints], in V (MNE convention)
    n_chan, n_timepoints = data.shape

    if work_dir is None:
        work_dir = tempfile.mkdtemp(prefix="amica_")
    out_dir = os.path.join(work_dir, "amica_out")
    os.makedirs(out_dir, exist_ok=True)

    data_path = os.path.join(work_dir, "amica_data.fdt")
    param_path = os.path.join(out_dir, "input.param")

    # Write data file
    logger.info(f"{log_prefix}Writing AMICA data file ({n_chan} ch × {n_timepoints} samples)...")
    _write_data_file(data, data_path)

    # Write param file
    _write_param_file(
        param_path=param_path,
        data_path=data_path,
        out_dir=out_dir,
        n_chan=n_chan,
        n_timepoints=n_timepoints,
        n_pcs=n_pcs,
        max_iter=max_iter,
        max_threads=max_threads,
    )

    # Run AMICA binary
    logger.info(
        f"{log_prefix}Running AMICA binary: {os.path.basename(binary)} "
        f"(n_pcs={n_pcs}, max_iter={max_iter}, max_threads={max_threads})"
    )
    result = subprocess.run(
        [binary, param_path],
        capture_output=False,   # let AMICA print to stdout (it logs iteration progress)
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"AMICA binary exited with code {result.returncode}. "
            f"Check output above for details."
        )

    # Read output
    logger.info(f"{log_prefix}Reading AMICA output from {out_dir}...")
    amica_results = _read_amica_output(out_dir)
    n_pcs_out = amica_results["n_pcs"]
    n_iters = len(amica_results.get("LL", []))
    logger.info(
        f"{log_prefix}AMICA finished: {n_pcs_out} components, "
        f"{n_iters} iterations completed."
    )

    # Build MNE ICA object
    ica = _build_mne_ica_from_amica(amica_results, raw, n_components=n_pcs_out)
    return ica
