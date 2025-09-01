import os

MATLAB_VERSION = "R2021a"  # version R2021a required for MoBILAB dependency

EEGLAB_DEPENDENCY = "eeglab2021.1"

MATLAB_DEPENDENCIES = [
    # MoBILAB
    "mobilab",

    # EEGLAB 2021.1
    EEGLAB_DEPENDENCY,

    # Plugins to EEGLAB
    # Amica
    os.path.join(EEGLAB_DEPENDENCY, "plugins", "amica"),
    # Dipfit
    os.path.join(EEGLAB_DEPENDENCY, "plugins", "dipfit5.3"),

]
