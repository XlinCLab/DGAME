import logging
import multiprocessing
import platform

import torch

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
