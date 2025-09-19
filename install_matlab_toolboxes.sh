#!/usr/bin/env bash
set -euo pipefail

# Detect OS
OS="$(uname -s)"

# Install dependencies on Linux/macOS if needed
install_command() {
    local cmd="$1"
    if ! command -v $cmd &> /dev/null; then
        echo "$cmd not found. Installing..."
        case "$OS" in
            Linux*) sudo apt-get update && sudo apt-get install -y $cmd ;;
            Darwin*) brew install $cmd ;;
            *) echo "Please install $cmd manually." ;;
        esac
    fi
}

install_command wget
install_command unzip
install_command git

# Default installation directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_TOOLBOX_DIR="${SCRIPT_DIR}/matlab"

# Use first argument as TOOLBOX_DIR if provided, otherwise default
TOOLBOX_DIR="${1:-$DEFAULT_TOOLBOX_DIR}"

echo "Installing toolboxes to: $TOOLBOX_DIR"

# EEGLAB settings
EEGLAB_VERSION="eeglab2021.1"
EEGLAB_DIR="${TOOLBOX_DIR}/${EEGLAB_VERSION}"

# Plugin settings
EEGLAB_PLUGINS_DIR="${EEGLAB_DIR}/plugins"
AMICA_DIR="${EEGLAB_PLUGINS_DIR}/amica"
MOBILAB_DIR="${EEGLAB_PLUGINS_DIR}/mobilab"

# Create toolbox and plugin directories
mkdir -p "$TOOLBOX_DIR" "$EEGLAB_PLUGINS_DIR" "$AMICA_DIR"

# Download EEGLAB
echo "Downloading EEGLAB ${EEGLAB_VERSION}..."
wget -O "$TOOLBOX_DIR/eeglab.zip" "https://sccn.ucsd.edu/eeglab/download/daily/${EEGLAB_VERSION}.zip"
unzip -q "$TOOLBOX_DIR/eeglab.zip" -d "$TOOLBOX_DIR"
rm "$TOOLBOX_DIR/eeglab.zip"

# Download AMICA plugin
echo "Downloading AMICA plugin..."
wget -O "$TOOLBOX_DIR/amica.zip" "https://sccn.ucsd.edu/~jason/amica1.5.zip"
unzip -q "$TOOLBOX_DIR/amica.zip" -d "$AMICA_DIR"
rm "$TOOLBOX_DIR/amica.zip"

# Clone MoBILAB toolbox
echo "Cloning MoBILAB toolbox..."
cd "$EEGLAB_PLUGINS_DIR"
if [ ! -d "mobilab" ]; then
    git clone https://github.com/sccn/mobilab.git mobilab
else
    echo "MoBILAB already exists, skipping clone."
fi

echo "All MATLAB toolboxes installed successfully!"
