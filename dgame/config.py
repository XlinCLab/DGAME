import os
from pathlib import Path

from dgame.paths import SCRIPT_DIR
from utils.run_config import load_config

CONFIG_DIR = os.path.join(Path(SCRIPT_DIR).parent.absolute(), "config")

# Fallback default config, used only to render the initial experiment setup GUI
DGAME_DEFAULT_CONFIG = load_config(
    os.path.join(
        Path(SCRIPT_DIR).parent.absolute(),
        "config",
        "dgame_defaults.yml",
    )
)

# Per-version default config file, keyed by dgame_version
SUPPORTED_DGAME_VERSIONS = {
    "2",
    "3",
}
DGAME_DEFAULT_CONFIG_FILES = {
    str(v): os.path.join(CONFIG_DIR, f"dgame{v}_defaults.yml")
    for v in SUPPORTED_DGAME_VERSIONS
}


# REQUIRED CONFIG FIELDS
REQUIRED_CONFIG_FIELDS = [
    r"^data.input\.*",
    r"^experiment\.objects",
    r"^experiment\.fillers",
]
