import torch

from dgame.audio import DEVICES


def get_device(device: str = None) -> str:
    """Select device on which to run Whisper."""
    if device is None:
        if torch.cuda.is_available():
            device = "cuda:0"
        elif torch.backends.mps.is_available():
            device = "mps"
        else:
            device = "cpu"
    elif device.lower() not in DEVICES:
        raise ValueError(f"Unrecognized device '{device}'. Please specify one of {DEVICES}.")
    return device
