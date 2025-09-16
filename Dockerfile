# Start from MATLAB container as base image
FROM mathworks/matlab:r2021a

# Add Linux timestamp as source for reproducibility and to find previously created Docker image 
ARG SOURCE_DATE_EPOCH=1758040333

# Switch to root for package installs
USER root

ENV DEBIAN_FRONTEND=noninteractive

# Update and install Python (3.8), R, and essential tools
RUN apt-get update && apt-get upgrade -y && \
    apt-get install -y --no-install-recommends \
        git build-essential \
        python3 python3-pip python3-venv python3-dev \
        software-properties-common curl gnupg && \
    # Add symlink from python to python3
    ln -sf /usr/bin/python3 /usr/bin/python && \
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

# Copy only requirements.txt first to leverage Docker caching
COPY ./requirements.txt /app/requirements.txt

# Create venv and install dependencies
RUN python3 -m venv /opt/venv && \
    /opt/venv/bin/pip install --upgrade pip && \
    /opt/venv/bin/pip install -r /app/requirements.txt

# Make the venv default for interactive shells
RUN echo "source /opt/venv/bin/activate" >> /etc/bash.bashrc

# Switch back to MATLAB default user
USER matlab
