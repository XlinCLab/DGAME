# Start from MATLAB container as base image
FROM mathworks/matlab:r2021a

# Add Linux timestamp as source for reproducibility and to find previously created Docker image 
ARG SOURCE_DATE_EPOCH=1758040333

# Switch to root for package installs
USER root

ENV DEBIAN_FRONTEND=noninteractive

# Install Python 3.11, R, and essential tools
RUN apt-get update && apt-get upgrade -y && \
    apt-get install -y --no-install-recommends \
        git build-essential software-properties-common curl gnupg \
        # System-level dependencies for installing certain R packages
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libfontconfig1-dev \
        libfreetype6-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libtiff5-dev \
        libcairo2-dev \
        libpng-dev \
        libtiff5-dev \
        libjpeg-dev \
        libwebp-dev \
        pkg-config && \
    # Add deadsnakes PPA for Python 3.11 (otherwise MATLAB base image with Ubuntu 20.04 supports only until Python 3.8)
    add-apt-repository ppa:deadsnakes/ppa -y && \
    apt-get update && \
    apt-get install -y --no-install-recommends \
        python3.11 python3.11-venv python3.11-dev python3.11-distutils python3-pip \
        # Also install python3.11-tk (tkinter) for Python GUI, which otherwise is not included in image by default 
        python3.11-tk && \
    # Symlinks so "python" and "python3" point to Python 3.11
    ln -sf /usr/bin/python3.11 /usr/bin/python3 && \
    ln -sf /usr/bin/python3.11 /usr/bin/python && \
    # Add CRAN repo and key for R
    curl -fsSL https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc \
        | gpg --dearmor -o /usr/share/keyrings/cran-archive-keyring.gpg && \
    echo "deb [signed-by=/usr/share/keyrings/cran-archive-keyring.gpg] https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/" \
        | tee /etc/apt/sources.list.d/cran.list > /dev/null && \
    # Install R
    apt-get update && \
    apt-get install -y --no-install-recommends r-base r-base-dev && \
    # Clean up
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Set default working directory where repo is mounted
WORKDIR /app

# Copy requirements.txt files for Python and R into container
COPY ./requirements.txt /app/requirements.txt
COPY ./r_requirements.txt /app/r_requirements.txt

# Pre-install R packages listed in r_requirements.txt
RUN R -e "packages <- readLines('/app/r_requirements.txt'); \
    install.packages(packages, repos='https://cloud.r-project.org', INSTALL_opts=c('--no-test-load','--no-multiarch'))"

# Create venv and install dependencies
RUN python3 -m venv /opt/venv && \
    /opt/venv/bin/pip install --upgrade pip && \
    /opt/venv/bin/pip install -r /app/requirements.txt

# Make the venv default for interactive shells
RUN echo "source /opt/venv/bin/activate" >> /etc/bash.bashrc

# Install MATLAB toolboxes
ENV TOOLBOX_DIR=/opt

# Download EEGLAB toolbox
ENV EEGLAB_VERSION=eeglab2021.1
ENV EEGLAB_DIR=${TOOLBOX_DIR}/${EEGLAB_VERSION}
RUN apt-get update && apt-get install -y wget unzip && \
    wget -O /tmp/eeglab.zip https://sccn.ucsd.edu/eeglab/download/daily/$EEGLAB_VERSION.zip && \
    unzip /tmp/eeglab.zip -d $TOOLBOX_DIR && \
    rm /tmp/eeglab.zip

# Download amica plugin for EEGLAB
ENV EEGLAB_PLUGINS_DIR=${EEGLAB_DIR}/plugins
ENV AMICA_DIR=${EEGLAB_PLUGINS_DIR}/amica
WORKDIR $AMICA_DIR
RUN wget -O /tmp/amica.zip https://sccn.ucsd.edu/~jason/amica1.5.zip && \
    unzip /tmp/amica.zip -d $AMICA_DIR && \
    rm /tmp/amica.zip

# Clone MoBILAB toolbox
WORKDIR $TOOLBOX_DIR
RUN git clone https://github.com/sccn/mobilab.git mobilab

# Set working directory to /app
WORKDIR /app
