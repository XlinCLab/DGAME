import os
from pathlib import Path

from dgame.paths import SCRIPT_DIR
from utils.run_config import load_config

# DGAME DEFAULT CONFIGS
DGAME_DEFAULT_CONFIG = load_config(
    os.path.join(
        Path(SCRIPT_DIR).parent.absolute(),
        "config",
        "dgame2_defaults.yml",
    )
)

# REQUIRED CONFIG FIELDS
REQUIRED_CONFIG_FIELDS = [
    r"^data.input\.*",
    r"^experiment\.objects",
    r"^experiment\.fillers",
]
