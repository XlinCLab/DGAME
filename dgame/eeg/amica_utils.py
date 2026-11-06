"""
Fit AMICA ICA with pAMICA and return an MNE ICA object that is compatible with mne-icalabel's label_components().

pAMICA's own defaults differ from EEGLAB's `runamica15.m`. To keep the decomposition comparable
with EEGLAB and with the AMICA binary previously used by DGAME, the `runamica15.m` settings are
passed explicitly as constants.
"""

import contextlib
import logging
import os
from collections.abc import Generator

import mne
import torch
from pamica.mne_compat import AMICAICA

logger = logging.getLogger(__name__)

# Number of mixture components per source (AMICA `num_mix_comps`)
AMICA_N_MIX = 3

# Fit settings of EEGLAB's runamica15.m (names are pAMICA's; EEGLAB's names in comments)
AMICA_FIT_PARAMS = {
    "lrate": 0.05,           # lrate
    "minlrate": 1e-8,        # minlrate
    "lratefact": 0.5,        # lratefact
    "rholratefact": 0.5,     # rholratefact
    "maxdecs": 3,            # max_decs
    "do_newton": True,       # do_newton
    "newt_start": 50,        # newt_start
    "newt_ramp": 10,         # newt_ramp
    "newtrate": 1.0,         # newtrate
    "do_reject": True,       # do_reject: iterative rejection of unlikely samples
    "maxrej": 10,            # numrej
    "rejsig": 2.5,           # rejsig
    "rejstart": 2,           # rejstart
    "rejint": 1,             # rejint
    "invsigmax": 100.0,      # invsigmax
    "invsigmin": 0.0,        # invsigmin
    "mineig": 1e-12,         # mineig
    "use_min_dll": True,     # use_min_dll
    "min_dll": 1e-9,         # min_dll
    "use_grad_norm": True,   # use_grad_norm
    "min_nd": 1e-7,          # min_grad_norm
}


@contextlib.contextmanager
def _torch_threads(n_threads: int | None) -> Generator[None]:
    """Temporarily set the number of CPU threads used by PyTorch."""
    if n_threads is None:
        yield
        return
    previous = torch.get_num_threads()
    torch.set_num_threads(n_threads)
    try:
        yield
    finally:
        torch.set_num_threads(previous)


@contextlib.contextmanager
def _log_to_file(path: str | None) -> Generator[None]:
    """Write pAMICA's log messages (e.g. sample rejections) and this module's to a file."""
    if path is None:
        yield
        return
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    handler = logging.FileHandler(path, mode="w")
    handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s: %(message)s"))
    targets = [logging.getLogger("pamica"), logger]
    previous_levels = [t.level for t in targets]
    for target in targets:
        target.addHandler(handler)
        if target.getEffectiveLevel() > logging.INFO:
            target.setLevel(logging.INFO)
    try:
        yield
    finally:
        for target, level in zip(targets, previous_levels):
            target.removeHandler(handler)
            target.setLevel(level)
        handler.close()


def run_amica(
    raw: mne.io.BaseRaw,
    n_pcs: int,
    max_iter: int = 2000,
    max_threads: int | None = None,
    device: str | None = None,
    random_state: int | None = None,
    log_prefix: str = "",
    amica_log_path: str | None = None,
) -> mne.preprocessing.ICA:
    """Run AMICA (pAMICA) on raw EEG data and return a fitted MNE ICA object.

    Parameters
    ----------
    raw : mne.io.BaseRaw
        Preprocessed, average-referenced EEG data (all bad channels already dropped).
    n_pcs : int
        Number of PCA components to retain (= data rank; matches MATLAB's pcakeep).
    max_iter : int
        Maximum AMICA iterations (default 2000, matches MATLAB pipeline).
    max_threads : int, optional
        Number of CPU threads for PyTorch during the fit. None keeps PyTorch's current setting.
    device : str, optional
        PyTorch device (e.g. "cpu", "cuda"). None selects automatically; on Apple silicon the
        float64 fit runs on the CPU, as the MPS device does not support float64.
    random_state : int, optional
        Seed for the AMICA initialization, for reproducible decompositions.
    log_prefix : str
        Optional string prepended to log messages (e.g. "Subject 02: ").
    amica_log_path : str, optional
        File path to write the AMICA fit log (sample rejections, log-likelihood trajectory).

    Returns
    -------
    mne.preprocessing.ICA
        ICA object with the AMICA decomposition, compatible with label_components().
    """
    n_chan = len(mne.pick_types(raw.info, eeg=True, exclude="bads"))
    logger.info(
        f"{log_prefix}Running AMICA (pAMICA): {n_chan} channels x {raw.n_times} samples, "
        f"n_pcs={n_pcs}, max_iter={max_iter}, max_threads={max_threads}, device={device or 'auto'}"
    )
    amica = AMICAICA(
        n_models=1,
        n_mix=AMICA_N_MIX,
        random_state=random_state,
        device=device,
        verbose=False,
    )
    with _log_to_file(amica_log_path), _torch_threads(max_threads):
        amica.fit(raw, picks="eeg", max_iter=max_iter, pcakeep=n_pcs, **AMICA_FIT_PARAMS)
        ll_history = amica.amica_.ll_history_
        logger.info(f"AMICA stop reason: {amica.stop_reason_}; {len(ll_history)} iterations")
        for i in range(0, len(ll_history), 10):
            logger.info(f"iteration {i + 1}: LL = {float(ll_history[i]):.6f}")
        if len(ll_history) > 0:
            logger.info(f"iteration {len(ll_history)}: LL = {float(ll_history[-1]):.6f}")

    if not amica.converged_:
        raise RuntimeError(
            f"{log_prefix}AMICA fit was degenerate (stop reason: {amica.stop_reason_})."
            + (f" See {amica_log_path}." if amica_log_path else "")
        )
    logger.info(
        f"{log_prefix}AMICA finished: {amica.n_components_} components, "
        f"{len(ll_history)} iterations completed ({amica.stop_reason_})."
    )
    return amica.to_mne_ica()
