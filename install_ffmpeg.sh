#!/bin/bash
set -euo pipefail

# Install ffmpeg, required by openai-whisper for audio decoding

if command -v ffmpeg >/dev/null 2>&1; then
    echo "ffmpeg is already installed."
    exit 0
fi

echo "ffmpeg not found. Attempting installation..."

OS="$(uname -s)"

case "$OS" in
    Linux)
        if command -v apt >/dev/null 2>&1; then
            echo "Installing ffmpeg via apt..."
            sudo apt update && sudo apt install -y ffmpeg

        elif command -v pacman >/dev/null 2>&1; then
            echo "Installing ffmpeg via pacman..."
            sudo pacman -S --noconfirm ffmpeg

        else
            echo "Unsupported Linux package manager. Please install ffmpeg manually."
            exit 1
        fi
        ;;

    Darwin)
        if command -v brew >/dev/null 2>&1; then
            echo "Installing ffmpeg via Homebrew..."
            brew install ffmpeg
        else
            echo "Homebrew not found. Install it from https://brew.sh/ and rerun this script."
            exit 1
        fi
        ;;

    MINGW*|MSYS*|CYGWIN*)
        if command -v choco >/dev/null 2>&1; then
            echo "Installing ffmpeg via Chocolatey..."
            choco install ffmpeg -y

        elif command -v scoop >/dev/null 2>&1; then
            echo "Installing ffmpeg via Scoop..."
            scoop install ffmpeg

        else
            echo "No Windows package manager found (Chocolatey or Scoop)."
            echo "Install ffmpeg manually from https://ffmpeg.org/download.html"
            exit 1
        fi
        ;;

    *)
        echo "Unsupported OS: $OS"
        exit 1
        ;;
esac

echo "ffmpeg installation step completed."
