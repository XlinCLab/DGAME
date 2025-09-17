# DGAME
## Introduction
This repo contains data processing and analysis scripts for "DGAME" experiments run in the University of Cologne's "XlinC" experimental linguistics lab, involving combined eye tracking and electroencephalography during referential selection in dyadic interaction.

## Setup
Running DGAME requires:
- Python 3.11
- R version 4.4.0 or higher
- MATLAB R2021a
    - MATLAB plug-ins:
        - EEGLAB 2021.1
            - [amica](https://sccn.ucsd.edu/~jason/amica_web.html)
            - dipfit5.3
        - MoBILAB

### Setup using Docker:
Instead of manually installing these dependencies, we provide a Dockerfile which builds an isolated container including all required analysis tools and dependencies.

- Download [Docker](https://docs.docker.com/get-started/get-docker/) to your machine.
    - NB: If using MacOS or Windows, you must explicitly open the Docker or Docker Desktop application before building or invoking Docker.
- Copy or move your input data into this `DGAME` repo's `data/` directory, which will be mounted inside the container. 
- Run the setup script, which will build and enter the Docker container:
```
./setup.sh
```

The same script can then subsequently be used to invoke and enter the already built Docker container.

### Setup without Docker:
Ensure that the versions of Python, R, and MATLAB (including plugins) specified above are installed on your machine. Then run the following commands to set up the Python virtual environment for DGAME:

```
# Create Python virtual environment
python3.11 -m venv venv

# Activate virtual environment
source venv/bin/activate

# Install requirements
venv/bin/pip install --upgrade pip
venv/bin/pip install -r requirements.txt
```

### Final Notes

Follow the instructions below to run DGAME analyses. Note that the first time the DGAME scripts are called from within the Docker container (and, depending on what R packages you already have installed, potentially also if not using Docker) it may take some time to install required R libraries first.
