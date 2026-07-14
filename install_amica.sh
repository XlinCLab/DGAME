#!/usr/bin/env bash
set -euo pipefail

# Install the AMICA standalone binary.
# Usage: ./install_amica.sh [DEST_DIR]
# Default destination: ./plugins/amica

OS="$(uname -s)"

install_command() {
    local cmd="$1"
    if ! command -v "$cmd" &> /dev/null; then
        echo "$cmd not found. Installing..."
        case "$OS" in
            Linux*) sudo apt-get update && sudo apt-get install -y "$cmd" ;;
            Darwin*) brew install "$cmd" ;;
            *) echo "Please install $cmd manually." ;;
        esac
    fi
}

install_command wget
install_command unzip

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AMICA_DIR="${1:-${SCRIPT_DIR}/plugins/amica}"

if [ -d "$AMICA_DIR" ] && [ -n "$(ls -A "$AMICA_DIR" 2>/dev/null)" ]; then
    echo "AMICA already installed at $AMICA_DIR, skipping..."
    exit 0
fi

echo "Installing AMICA to $AMICA_DIR..."
mkdir -p "$AMICA_DIR"
TMP_ZIP="$(mktemp).zip"
wget -O "$TMP_ZIP" "https://sccn.ucsd.edu/~jason/amica1.5.zip"
unzip -q "$TMP_ZIP" -d "$AMICA_DIR"
rm "$TMP_ZIP"
echo "Successfully installed AMICA to $AMICA_DIR"
