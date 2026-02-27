# DGAME
## Table of Contents
* [Introduction](#introduction)
* [Setup](#setup)
    * [Setup using Docker](#setup-using-docker)
    * [Manual setup without Docker](#manual-setup-without-docker)
* [Running DGAME](#running-dgame)
    * [Experiment config](#experiment-configyml-file)
    * [Running analysis](#running-analysis)

## Introduction
This repo contains data processing and analysis scripts for "DGAME" experiments run in the University of Cologne's "XlinC" experimental linguistics lab, involving combined eye tracking and electroencephalography during referential selection in dyadic interaction.

## Setup
Running DGAME requires:
- [Python 3.11](https://www.python.org/downloads/release/python-3110/)
- [R](https://www.r-project.org/) (version 4.4.0 or higher)
- [MATLAB](https://www.mathworks.com/help/install/ug/install-products-with-internet-connection.html) (version R2024b)
    - Toolboxes:
        - [EEGLAB 2021.1](https://sccn.ucsd.edu/eeglab/download/daily/eeglab2021.1.zip)
            - Plugins:
                - [amica](https://sccn.ucsd.edu/~jason/amica_web.html)
                - [CleanLine](https://github.com/sccn/cleanline)
                - [dipfit](https://eeglab.org/plugins/dipfit/) (included by default in EEGLAB 2021.1)
                - [unfold](https://www.unfoldtoolbox.org/) (NB: Requires [forked branch](https://github.com/XlinCLab/unfold/tree/DGAME2) with slight adjustments)
                - [xdf-EEGLAB](https://github.com/xdf-modules/xdf-EEGLAB/)

Follow the instructions below to set up the environment and install dependencies.

### Setup using Docker
Instead of manually installing these dependencies, we provide a Dockerfile which builds an isolated Docker image and container including all required analysis tools and dependencies.

#### Prerequisites
- Create a [MathWorks account](https://www.mathworks.com/mwaccount/account/create?uri=) and license (for running MATLAB)
- Download [Docker](https://docs.docker.com/get-started/get-docker/) to your machine.
    - NB: If using MacOS or Windows, you must explicitly open the Docker or Docker Desktop application before building or invoking Docker.

#### Docker container setup
- Copy or move your input data into this `DGAME` repo's `data/` directory, which will be mounted inside the container.
    - NB: This step is optional during setup and can instead be done later. However, the Docker container will not have access to your device's file system outside of this `DGAME` directory, so input data would eventually need to be copied here.
- Run the setup script, which will build the Docker image and launch the Docker container:
```
./setup.sh
```
Note that it may take some time for all dependencies to be installed the first time, upwards of 30 minutes, depending on your machine.
The same setup script can subsequently be reused to run and enter the Docker container, which should be immediate once its image has been built.

#### MATLAB authentication
Once you have entered the container, run the following command to launch MATLAB and authenticate with your MathWorks account credentials:
```
matlab -nodesktop -nosplash
```
Once you have authenticated, you can exit the MATLAB shell and follow further instructions below to begin running a DGAME experiment analysis.
```
exit
```

#### (Optional) Verify MATLAB toolboxes in license
Optionally, you can verify whether the required toolboxes are included in your MathWorks license. Open a MATLAB shell:
```
matlab -nodesktop -nosplash
```
Run the following in the MATLAB shell:
```
license('test', 'Distrib_Computing_Toolbox')
license('test', 'Image_Toolbox')
license('test', 'Optimization_Toolbox')
license('test', 'Signal_Toolbox')
license('test', 'Statistics_Toolbox')
```
If the toolbox is installed, you should see, e.g.:
```
>> license('test', 'Distrib_Computing_Toolbox')

ans =

     1
```


If any of the above return `ans = 0` instead of `ans = 1`, then the relevant toolbox is NOT included in your license and DGAME analyses may not run as expected.


### Manual setup without Docker
Ensure that the versions of Python, R, and MATLAB specified above are installed on your machine.

#### Python environment
Run the following commands to set up the Python virtual environment for DGAME:
```
# Create Python virtual environment
python3.11 -m venv venv

# Activate virtual environment
source venv/bin/activate

# Install requirements
venv/bin/pip install --upgrade pip
venv/bin/pip install -r requirements.txt
```

#### MATLAB toolboxes and plugins
Run the following script to install the required MATLAB toolboxes and plugins:
```
./install_matlab_plugins.sh
```
The MATLAB toolbox and plugin dependencies will be installed into a new `matlab` directory within this repo.

(!) Note that if you install or have previously installed these MATLAB dependencies into some location other than `./matlab` within this directory, this must be specified in your DGAME experiment `config.yml` file.

## Running DGAME
(!) Please follow above instructions to set up the DGAME analysis environment before proceeding to these steps.

### Input data
The expected input data directory structure for DGAME is the following, given an analysis of two test subjects with participant IDs `testsubject01` and `testsubject02`:
```
dgame_data_root/
    preproc/
        audio/
            testsubject01/
                testsubject01_words_11.csv
                testsubject01_words_12.csv
                testsubject01_words_21.csv
                testsubject01_words_22.csv
                testsubject01_words2erp_11.csv
                testsubject01_words2erp_12.csv
                testsubject01_words2erp_21.csv
                testsubject01_words2erp_22.csv
            testsubject02/
                testsubject02_words_11.csv
                testsubject02_words_12.csv
                testsubject02_words_21.csv
                testsubject02_words_22.csv
                testsubject02_words2erp_11.csv
                testsubject02_words2erp_12.csv
                testsubject02_words2erp_21.csv
                testsubject02_words2erp_22.csv
        eyetracking/
            fixations/
                testsubject01/
                    fixations_times_11_trials.csv
                    fixations_times_12_trials.csv
                    fixations_times_21_trials.csv
                    fixations_times_22_trials.csv
                testsubject02/
                    fixations_times_11_trials.csv
                    fixations_times_12_trials.csv
                    fixations_times_21_trials.csv
                    fixations_times_22_trials.csv
            gaze_positions/
                testsubject01/
                    gaze_positions.csv
                testsubject02/
                    gaze_positions.csv
            surfaces/
                testsubject01/
                    fixations_on_surface_11.csv
                    fixations_on_surface_12.csv
                    fixations_on_surface_13.csv
                    fixations_on_surface_14.csv
                    ...
                    gaze_positions_on_surface_14.csv
                    gaze_positions_on_surface_14.csv
                    gaze_positions_on_surface_14.csv
                    gaze_positions_on_surface_14.csv
                    ...
                testsubject02/
                    fixations_on_surface_11.csv
                    fixations_on_surface_12.csv
                    fixations_on_surface_13.csv
                    fixations_on_surface_14.csv
                    ...
                    gaze_positions_on_surface_14.csv
                    gaze_positions_on_surface_14.csv
                    gaze_positions_on_surface_14.csv
                    gaze_positions_on_surface_14.csv
                    ...
        object_positions/
            testsubject01/
                object_positions.csv
            testsubject02/
                object_positions.csv
    recordings/
        xdf/
            testsubject01/
            testsubject02/
```

If the input data directory does not match the expected structure or the expected input files are not found, a validation error may occur.

### Experiment `config.yml` file
To run DGAME data postprocessing and/or analysis, create a `config.yml` file which defines several experimental parameters, including:
- Location of your input data
- Desired location for output files
- Experimental subjects to process
- DGAME setup parameters (e.g. target and filler words)
- Desired processing/analysis steps

A full config file with default values as a template is saved at:
[`config/dgame2_defaults.yml`](config/dgame2_defaults.yml)

Note that any parameters left unspecified in your experimental config file will be automatically inherited from this default config. Therefore, only experimental parameters which either have no default values (e.g. target and filler words) or which differ from the defaults need to be specified explicitly in your config file.
A sample config file specifying only those minimally required parameters is saved at:
 [`config/sample_config.yml`](config/sample_config.yml)

Specific DGAME analysis steps can be enabled or disabled by setting the relevant `enabled` parameter to either `true` or `false`, e.g. adding the following block to your experimental config file would disable step `A_export_audio_and_et_times` while keeping other analysis steps enabled.
```
analysis:
  steps:
    A_export_audio_and_et_times:
      enabled: false
```
(!) Please note that certain steps depend upon outputs of earlier processing steps, and therefore would not work as expected in isolation.

### Running analysis
Once your environment and config file are set up, simply run the following from the repository root:
```
python3 run_dgame_analysis.py --config /path/to/your/config.yml
```

