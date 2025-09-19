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

# Create toolbox and plugin directories
mkdir -p "$TOOLBOX_DIR" "$EEGLAB_PLUGINS_DIR" "$AMICA_DIR"

# Download EEGLAB
echo "Installing EEGLAB ${EEGLAB_VERSION}..."
wget -O "$TOOLBOX_DIR/eeglab.zip" "https://sccn.ucsd.edu/eeglab/download/daily/${EEGLAB_VERSION}.zip"
unzip -q "$TOOLBOX_DIR/eeglab.zip" -d "$TOOLBOX_DIR"
rm "$TOOLBOX_DIR/eeglab.zip"
echo "✅️ Successfully installed $EEGLAB_VERSION"

# Download ERPLAB as plugin to EEGLAB
echo "Installing ERPLAB plugin..."
wget -O "$TOOLBOX_DIR/erplab.zip" "https://github.com/ucdavis/erplab/releases/download/12.00/erplab12.01.zip"
unzip -q "$TOOLBOX_DIR/erplab.zip" -d "$ERPLAB_DIR"
rm "$TOOLBOX_DIR/erplab.zip"
echo "✅️ Successfully installed ERPLAB"

# Download AMICA plugin
echo "Installing AMICA plugin..."
wget -O "$TOOLBOX_DIR/amica.zip" "https://sccn.ucsd.edu/~jason/amica1.5.zip"
unzip -q "$TOOLBOX_DIR/amica.zip" -d "$AMICA_DIR"
rm "$TOOLBOX_DIR/amica.zip"
echo "✅️ Successfully installed amica"

# Clone CleanLine as plugin to EEGLAB
echo "Installing CleanLine plugin..."
cd "$EEGLAB_PLUGINS_DIR"
git clone https://github.com/sccn/cleanline
echo "✅️ Successfully installed CleanLine"

# Clone MoBILAB as plugin to EEGLAB
echo "Installing MoBILAB toolbox as EEGLAB plugin..."
cd "$EEGLAB_PLUGINS_DIR"
git clone https://github.com/sccn/mobilab.git
echo "✅️ Successfully installed MoBILAB"

# Download xdf-EEGLAB plugin
# Ensure xdf submodule is initiated and linked
cd "$EEGLAB_PLUGINS_DIR"
echo "Installing xdf-EEGLAB plugin..."
git clone https://github.com/xdf-modules/xdf-EEGLAB.git
cd "$EEGLAB_PLUGINS_DIR"/xdf-EEGLAB
git submodule update --init --recursive
echo "✅️ Successfully installed xdf-EEGLAB"

echo "✅️ All MATLAB toolboxes and plugins installed successfully!"
