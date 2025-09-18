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
IMAGE_NAME="dgame"
CONTAINER_NAME="dgame_container"

# Adjust platform based on architecture
if [[ "$ARCH" == "arm64" || "$ARCH" == "aarch64" ]]; then
    #echo "Running on ARM (e.g., Apple Silicon). Forcing amd64 platform."
    PLATFORM="--platform=linux/amd64"
fi

# Check if Docker image "dgame" exists; if not, build it
if [[ -z $(docker images -q $IMAGE_NAME) ]]; then
    echo "Building Docker image '$IMAGE_NAME'..."
    docker build $PLATFORM -t $IMAGE_NAME .
fi

# Check container status
if docker ps -a --format '{{.Names}}' | grep -Eq "^${CONTAINER_NAME}\$"; then
    if docker ps --format '{{.Names}}' | grep -Eq "^${CONTAINER_NAME}\$"; then
        echo "Attaching to running container '$CONTAINER_NAME'..."
        docker exec -it $CONTAINER_NAME /bin/bash
    else
        echo "Starting existing container '$CONTAINER_NAME'..."
        docker start -ai $CONTAINER_NAME
    fi
else
    # Run the container with current repo mounted at /app
    echo "Creating and running new container '$CONTAINER_NAME'..."
    docker run -it $PLATFORM \
        -v "$(pwd)":/app \
        --name $CONTAINER_NAME \
        $IMAGE_NAME /bin/bash
fi
