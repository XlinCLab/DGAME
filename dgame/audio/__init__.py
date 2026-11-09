import logging
import multiprocessing
import os
import platform
import sys

import torch

# Whisper uses Numba for its word-level timestamp alignment, and Numba's default OpenMP threading layer is unsafe in this process:
# PyTorch and R (embedded via rpy2) each bring their own copy of the OpenMP runtime, and the layer can end
# up mixing them, which crashes the whole process with a segmentation fault.
# Numba's "workqueue" threading layer does not depend on OpenMP.
NUMBA_THREADING_LAYER = "workqueue"
os.environ.setdefault("NUMBA_THREADING_LAYER", NUMBA_THREADING_LAYER)
if "numba" in sys.modules:  # already imported by then (the environment variable is read on import)
    sys.modules["numba"].config.THREADING_LAYER = os.environ["NUMBA_THREADING_LAYER"]

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(name)s %(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

DEVICES = {"cuda", "cpu", "mps"}
MAX_CPU_WORKERS = multiprocessing.cpu_count()


def _patch_mps_dtw() -> None:
    """
    whisper.timing.dtw casts its cost matrix to float64 before moving it off
    the device (`x.double().cpu().numpy()`), which raises a hard TypeError on
    MPS tensors since Apple's MPS backend does not support float64 at all.
    This only affects word-level timestamp alignment; move to CPU first
    (a no-op for already-CPU/CUDA tensors) so word_timestamps=True works on MPS.
    """
    from whisper import timing

    original_dtw = timing.dtw

    def patched_dtw(x: torch.Tensor):
        if not x.is_cuda:
            x = x.cpu()
        return original_dtw(x)

    timing.dtw = patched_dtw


# Run patch on macOS only
if platform.system() == "Darwin":
    logger.warning("Patching whisper.timing.dtw to support MPS device for macOS")
    _patch_mps_dtw()
