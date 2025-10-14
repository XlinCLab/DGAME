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
ERPLAB_DIR="${EEGLAB_PLUGINS_DIR}/erplab"

# Create toolbox directory
mkdir -p "$TOOLBOX_DIR"

# -----------------------------
# EEGLAB
# -----------------------------
if [ -d "$EEGLAB_DIR" ]; then
    echo "EEGLAB already installed at $EEGLAB_DIR, skipping..."
else
    mkdir -p "$EEGLAB_PLUGINS_DIR"
    echo "Installing EEGLAB ${EEGLAB_VERSION}..."
    wget -O "$TOOLBOX_DIR/eeglab.zip" "https://sccn.ucsd.edu/eeglab/download/daily/${EEGLAB_VERSION}.zip"
    unzip -q "$TOOLBOX_DIR/eeglab.zip" -d "$TOOLBOX_DIR"
    rm "$TOOLBOX_DIR/eeglab.zip"
    echo "✅️ Successfully installed $EEGLAB_VERSION"
fi

# -----------------------------
# ERPLAB
# -----------------------------
if [ -d "$ERPLAB_DIR" ]; then
    echo "ERPLAB plugin already installed at $ERPLAB_DIR, skipping..."
else
    echo "Installing ERPLAB plugin..."
    mkdir -p "$ERPLAB_DIR"
    wget -O "$TOOLBOX_DIR/erplab.zip" "https://github.com/ucdavis/erplab/releases/download/12.00/erplab12.01.zip"
    unzip -q "$TOOLBOX_DIR/erplab.zip" -d "$ERPLAB_DIR"
    rm "$TOOLBOX_DIR/erplab.zip"
    echo "✅️ Successfully installed ERPLAB"
fi

# -----------------------------
# AMICA
# -----------------------------
if [ -d "$AMICA_DIR" ]; then
    echo "AMICA plugin already installed at $AMICA_DIR, skipping..."
else
    echo "Installing AMICA plugin..."
    mkdir -p "$AMICA_DIR"
    wget -O "$TOOLBOX_DIR/amica.zip" "https://sccn.ucsd.edu/~jason/amica1.5.zip"
    unzip -q "$TOOLBOX_DIR/amica.zip" -d "$AMICA_DIR"
    rm "$TOOLBOX_DIR/amica.zip"
    echo "✅️ Successfully installed AMICA"
fi

# -----------------------------
# CleanLine
# -----------------------------
if [ -d "$EEGLAB_PLUGINS_DIR/cleanline" ]; then
    echo "CleanLine plugin already installed at $EEGLAB_PLUGINS_DIR/cleanline, skipping..."
else
    echo "Installing CleanLine plugin..."
    git -C "$EEGLAB_PLUGINS_DIR" clone https://github.com/sccn/cleanline
    echo "✅️ Successfully installed CleanLine"
fi

# -----------------------------
# MoBILAB
# -----------------------------
if [ -d "$EEGLAB_PLUGINS_DIR/mobilab" ]; then
    echo "MoBILAB plugin already installed at $EEGLAB_PLUGINS_DIR/mobilab, skipping..."
else
    echo "Installing MoBILAB toolbox as EEGLAB plugin..."
    git -C "$EEGLAB_PLUGINS_DIR" clone https://github.com/sccn/mobilab.git
    echo "✅️ Successfully installed MoBILAB"
fi

# -----------------------------
# unfold
# -----------------------------
if [ -d "$EEGLAB_PLUGINS_DIR/unfold" ]; then
    echo "unfold plugin already installed at $EEGLAB_PLUGINS_DIR/unfold, skipping..."
else
    echo "Installing unfold plugin..."
    git -C "$EEGLAB_PLUGINS_DIR" clone https://github.com/XlinCLab/unfold.git
    cd "$EEGLAB_PLUGINS_DIR/unfold"
    # Initialize and link submodules
    git submodule update --init --recursive --remote
    # Check out DGAME2 branch with adjustments
    git switch DGAME2
    git pull -v
    echo "✅️ Successfully installed unfold"
fi

# -----------------------------
# xdf-EEGLAB
# -----------------------------
if [ -d "$EEGLAB_PLUGINS_DIR/xdf-EEGLAB" ]; then
    echo "xdf-EEGLAB plugin already installed at $EEGLAB_PLUGINS_DIR/xdf-EEGLAB, skipping..."
else
    echo "Installing xdf-EEGLAB plugin..."
    git -C "$EEGLAB_PLUGINS_DIR" clone https://github.com/xdf-modules/xdf-EEGLAB.git
    cd "$EEGLAB_PLUGINS_DIR/xdf-EEGLAB"
    # Initialize and link submodules
    git submodule update --init --recursive
    echo "✅️ Successfully installed xdf-EEGLAB"
fi

echo "✅️ All MATLAB toolboxes and plugins installed successfully!"
