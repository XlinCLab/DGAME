# DGAME
## Introduction
This repo contains data processing and analysis scripts for "DGAME" experiments run in the University of Cologne's "XlinC" experimental linguistics lab, involving combined eye tracking and electroencephalography during referential selection in dyadic interaction.

## Setup
Running DGAME requires:
- [Python 3.11](https://www.python.org/downloads/release/python-3110/)
- [R](https://www.r-project.org/) (version 4.4.0 or higher)
- [MATLAB](https://www.mathworks.com/help/install/ug/install-products-with-internet-connection.html) (version R2021a)
    - Toolboxes:
        - [EEGLAB 2021.1](https://sccn.ucsd.edu/eeglab/download/daily/eeglab2021.1.zip)
            - Plugins:
                - [amica](https://sccn.ucsd.edu/~jason/amica_web.html)
                - [CleanLine](https://github.com/sccn/cleanline)
                - [dipfit](https://eeglab.org/plugins/dipfit/) (included by default in EEGLAB 2021.1)
                - [ERPLAB](https://erpinfo.org/erplab)
                - [MoBILAB](https://github.com/sccn/mobilab)
                - [unfold](https://www.unfoldtoolbox.org/)
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
