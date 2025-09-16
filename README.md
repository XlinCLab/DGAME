# DGAME
## Introduction
This repo contains data processing and analysis scripts for "DGAME" experiments run in the University of Cologne's "XlinC" experimental linguistics lab, involving combined eye tracking and electroencephalography during referential selection in dyadic interaction.

## Setup
Running DGAME requires:
- Python 3.8 or higher
- R version 4.4.0 or higher
- MATLAB R2021a
    - MATLAB plug-ins:
        - EEGLAB 2021.1
            - [amica](https://sccn.ucsd.edu/~jason/amica_web.html)
            - dipfit5.3
        - MoBILAB

Instead of manually installing these dependencies, we provide a Dockerfile which builds an isolated container including all required analysis tools and dependencies.

Setup using Docker:
- Download [Docker](https://docs.docker.com/get-started/get-docker/) to your machine.
    - NB: If using MacOS or Windows, you must explicitly open the Docker or Docker Desktop application before building or invoking Docker.
- Copy or move your input data into this `DGAME` repo's `data/` directory, which will be mounted inside the container. 
- Run the setup script, which will build and enter the Docker container:
```
./setup.sh
```

Then follow the instructions below to run DGAME analyses.
