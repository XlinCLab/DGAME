import os

SUPPORTED_MATLAB_VERSIONS = [
    "R2021a",
    "R2021b",
    "R2022a",
    "R2022b",
    "R2023a",
    "R2023b",
    "R2024a",
    "R2024b",
]
LATEST_MATLAB_VERSION = sorted(SUPPORTED_MATLAB_VERSIONS)[-1]

EEGLAB_DEPENDENCY = "eeglab2021.1"
EEGLAB_PLUGIN_PATH = os.path.join(EEGLAB_DEPENDENCY, "plugins")

MATLAB_DEPENDENCIES = [
    # EEGLAB 2021.1
    EEGLAB_DEPENDENCY,

    # Plugins to EEGLAB
    # Amica
    os.path.join(EEGLAB_PLUGIN_PATH, "amica"),
    # CleanLine
    os.path.join(EEGLAB_PLUGIN_PATH, "cleanline"),
    # Dipfit
    os.path.join(EEGLAB_PLUGIN_PATH, "dipfit"),
    # unfold
    os.path.join(EEGLAB_PLUGIN_PATH, "unfold"),
    # xdf-EEGLAB
    os.path.join(EEGLAB_PLUGIN_PATH, "xdf-EEGLAB"),

]
