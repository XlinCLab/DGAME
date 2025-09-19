import os

MATLAB_VERSION = "R2021a"  # version R2021a required for MoBILAB dependency

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
    # ERPLAB
    os.path.join(EEGLAB_PLUGIN_PATH, "erplab"),
    # Dipfit
    os.path.join(EEGLAB_PLUGIN_PATH, "dipfit"),
    # MoBILAB
    os.path.join(EEGLAB_PLUGIN_PATH, "mobilab"),
    # xdf-EEGLAB
    os.path.join(EEGLAB_PLUGIN_PATH, "xdf-EEGLAB"),

]
