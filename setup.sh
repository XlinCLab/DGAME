#!/bin/bash
# setup.sh
set -euo pipefail

# Detect OS and architecture
OS=$(uname -s | tr '[:upper:]' '[:lower:]')
ARCH=$(uname -m)

#echo "Detected OS: $OS"
#echo "Detected Arch: $ARCH"

# Default values
PLATFORM=""

# Adjust platform based on architecture
if [[ "$ARCH" == "arm64" || "$ARCH" == "aarch64" ]]; then
    #echo "Running on ARM (e.g., Apple Silicon). Forcing amd64 platform."
    PLATFORM="--platform=linux/amd64"
fi

# Build Docker image
docker build $PLATFORM -t dgame .

# Run the container with current repo mounted at /app
docker run -it $PLATFORM \
    -v "$(pwd)":/app \
    --name dgame_container \
    dgame /bin/bash
