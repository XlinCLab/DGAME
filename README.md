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
- [Julia](https://julialang.org/)
- [R](https://www.r-project.org/) (version 4.4.0 or higher)

For EEG preprocessing with AMICA ICA (optional):
- [AMICA](https://sccn.ucsd.edu/~jason/amica_web.html) standalone binary (see [AMICA installation](#amica-ica-binary-optional) below)

Follow the instructions below to set up the environment and install dependencies.

### Setup using Docker
Instead of manually installing these dependencies, we provide a Dockerfile which builds an isolated Docker image and container including all required analysis tools and dependencies.

#### Prerequisites
- Download [Docker](https://docs.docker.com/get-started/get-docker/) to your machine.
    - NB: If using MacOS or Windows, you must explicitly open the Docker or Docker Desktop application before building or invoking Docker.

#### Docker container setup
- Copy or move your input data into this `DGAME` repo's `data/` directory, which will be mounted inside the container.
    - NB: This step is optional during setup and can instead be done later. However, the Docker container will not have access to your device's file system outside of this `DGAME` directory, so input data would eventually need to be copied here.
- Run the setup script, which will build the Docker image and launch the Docker container:
```
./setup.sh
```
The setup script performs the following steps in order:
1. Installs Julia on the host machine (skipped if already installed).
2. Installs the AMICA ICA binary into `./plugins/amica/` on the host (skipped if already installed). This directory is mounted inside the container and accessible at runtime if `ica.method: amica` is set in your config.
3. Builds the Docker image (Python 3.11, R 4.4.x, and Julia are installed inside the image).
4. Starts the Docker container with this repository mounted at `/app`.

Note that it may take some time (upwards of 30 minutes) for all dependencies to be installed the first time, depending on your machine.
The same setup script can subsequently be reused to run and enter the Docker container, which should be immediate once its image has been built.

On first DGAME run inside the container, Julia package dependencies will be installed automatically; this adds a one-time delay of a few minutes.

### Manual setup without Docker
Ensure that the versions of Python, R, and Julia specified above are installed on your machine.

#### Python environment
Run the following commands to set up the Python virtual environment for DGAME:
```bash
# Create Python virtual environment
python3.11 -m venv venv

# Activate virtual environment
source venv/bin/activate

# Install requirements
venv/bin/pip install --upgrade pip
venv/bin/pip install -r requirements.txt
```

#### Julia installation
Run the following script to install Julia, if not previously installed on your machine:
```bash
./install_julia.sh
```

Julia package dependencies will be automatically installed upon running the `DGAME` code. 

#### AMICA ICA binary (optional)
AMICA is an optional ICA algorithm for EEG preprocessing. It is only required if you set `ica.method: amica` in your experiment config. Run the following script to install the AMICA binary:
```bash
./install_amica.sh
```
The binary will be installed into `./plugins/amica/`. If you install it elsewhere, specify the path via `analysis.dependencies.amica.dir` in your config.

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

Set `experiment.dgame_version` to the DGAME experiment version you are processing (currently supported: `2` and `3`). A full config file with default values as a template for each version is saved at:
[`config/dgame2_defaults.yml`](config/dgame2_defaults.yml) and [`config/dgame3_defaults.yml`](config/dgame3_defaults.yml), respectively.

Note that any parameters left unspecified in your experimental config file will be automatically inherited from the relevant version's default config. Therefore, only experimental parameters which either have no default values (e.g. target and filler words) or which differ from the defaults need to be specified explicitly in your config file. If `experiment.dgame_version` itself is left unspecified, it defaults to `2`.
A sample config file specifying only those minimally required parameters is saved at:
 [`config/sample_config.yml`](config/sample_config.yml)

Specific DGAME analysis steps can be enabled or disabled by setting the relevant `enabled` parameter to either `true` or `false`, e.g. adding the following block to your experimental config file would disable step `A_export_audio_and_et_times` while keeping other analysis steps enabled.
```yml
analysis:
  steps:
    A_export_audio_and_et_times:
      enabled: false
```
(!) Please note that certain steps depend upon outputs of earlier processing steps, and therefore would not work as expected in isolation.

### Running analysis
Once your environment and config file are set up, simply run the following from the repository root:
```bash
python3 run_dgame_analysis.py --config /path/to/your/config.yml
```

