#!/usr/bin/env bash

set -e

REPO_DIR="$(cd "$(dirname "$0")" && pwd)"
VENV_DIR="$REPO_DIR/.venv"
CRISPERWHISPER_DIR="$REPO_DIR/../CrisperWhisper"

echo "Setting up Python virtual environment..."
if [ -d "$VENV_DIR" ]; then
    echo "Virtual environment already exists at $VENV_DIR, reusing it."
else
    python3 -m venv "$VENV_DIR"
fi
source "$VENV_DIR/bin/activate"

pip install --upgrade pip

# Pin torchvision/torchaudio to the same torch release as requirements.txt (rather than
# leaving them unpinned) so pip resolves a mutually compatible set in one step. Otherwise,
# installing unpinned torch/torchvision/torchaudio here and then separately installing
# requirements.txt's exact torch==X.Y.Z pin later leaves an orphaned, mismatched torchvision.
TORCH_VERSION=$(grep -oP '^torch==\K[0-9.]+' "$REPO_DIR/requirements.txt")

# Install PyTorch: use ROCm build on AMD GPU systems, otherwise default (CPU/CUDA)
ROCM_VERSION=$(dpkg -l 2>/dev/null | grep 'rocm-smi ' | grep -oP '\d+\.\d+(?=\.\d+-)' | head -1)
if [ -n "$ROCM_VERSION" ]; then
    echo "ROCm detected (version ${ROCM_VERSION}); installing PyTorch with ROCm ${ROCM_VERSION} support..."
    # ROCm 5.7 tops out at torch 2.3.1; ROCm 6.0+ supports torch 2.4+. This intentionally
    # does NOT match TORCH_VERSION (requirements.txt targets a newer, non-ROCm build) -
    # the plain "torch" entry is excluded from the requirements.txt install below so this
    # ROCm-specific build doesn't get silently overwritten.
    if [[ "$ROCM_VERSION" == 5.* ]]; then
        pip install torch==2.3.1 torchvision torchaudio --index-url "https://download.pytorch.org/whl/rocm${ROCM_VERSION}"
    else
        pip install torch==2.4.0 torchvision torchaudio --index-url "https://download.pytorch.org/whl/rocm${ROCM_VERSION}"
    fi
else
    echo "No ROCm detected; installing default PyTorch (CPU/CUDA)..."
    pip install "torch==${TORCH_VERSION}" torchvision torchaudio
fi

# Install remaining requirements, excluding the "torch" line: its version/build (CPU, CUDA,
# or ROCm-specific) was already selected above and must not be silently overwritten
pip install -r <(grep -v '^torch==' "$REPO_DIR/requirements.txt")

# Some AMD GPUs (e.g. gfx1032 / RX 6600) are not included in prebuilt PyTorch ROCm
# wheels, which target gfx1030. Setting HSA_OVERRIDE_GFX_VERSION=10.3.0 tells the
# ROCm runtime to use gfx1030 kernels, which are compatible enough to work.
GPU_ARCH=$(rocminfo 2>/dev/null | grep -oP 'gfx\d+' | head -1)
if [ "$GPU_ARCH" = "gfx1032" ]; then
    ACTIVATE="$VENV_DIR/bin/activate"
    if ! grep -q "HSA_OVERRIDE_GFX_VERSION" "$ACTIVATE"; then
        echo "" >> "$ACTIVATE"
        echo "# ROCm gfx1032 workaround: use gfx1030 kernels from PyTorch wheel" >> "$ACTIVATE"
        echo "export HSA_OVERRIDE_GFX_VERSION=10.3.0" >> "$ACTIVATE"
    fi
    echo "gfx1032 GPU detected; added HSA_OVERRIDE_GFX_VERSION=10.3.0 to $ACTIVATE"
    deactivate
    source "$ACTIVATE"
fi

echo "Downloading spaCy German transformer model..."
python3 -m spacy download de_dep_news_trf

echo "Checking for ffmpeg..."

if command -v ffmpeg >/dev/null 2>&1; then
    echo "ffmpeg is already installed."
else
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
fi

# CrisperWhisper is not pip-installable; check it out as a sibling directory of this repo
# (i.e. "../CrisperWhisper" relative to DGAME) rather than nesting it inside DGAME
echo "Setting up CrisperWhisper..."
if [ -d "$CRISPERWHISPER_DIR" ]; then
    echo "CrisperWhisper already present at $CRISPERWHISPER_DIR, skipping clone."
else
    echo "Cloning CrisperWhisper into $CRISPERWHISPER_DIR..."
    git clone https://github.com/nyrahealth/CrisperWhisper "$CRISPERWHISPER_DIR"
fi

echo "Installing requirements for CrisperWhisper..."
pip install -r "$CRISPERWHISPER_DIR/requirements.txt"

echo "Checking for git-lfs..."

if command -v git-lfs >/dev/null 2>&1; then
    echo "git-lfs is already installed."
else
    echo "git-lfs not found. Attempting installation..."

    OS="$(uname -s)"

    case "$OS" in
        Linux)
            if command -v apt >/dev/null 2>&1; then
                echo "Installing git-lfs via apt..."
                sudo apt update && sudo apt install -y git-lfs

            elif command -v pacman >/dev/null 2>&1; then
                echo "Installing git-lfs via pacman..."
                sudo pacman -S --noconfirm git-lfs

            else
                echo "Unsupported Linux package manager. Please install git-lfs manually."
                exit 1
            fi
            ;;

        Darwin)
            if command -v brew >/dev/null 2>&1; then
                echo "Installing git-lfs via Homebrew..."
                brew install git-lfs
            else
                echo "Homebrew not found. Install it from https://brew.sh/ and rerun this script."
                exit 1
            fi
            ;;

        MINGW*|MSYS*|CYGWIN*)
            if command -v choco >/dev/null 2>&1; then
                echo "Installing git-lfs via Chocolatey..."
                choco install git-lfs -y

            elif command -v scoop >/dev/null 2>&1; then
                echo "Installing git-lfs via Scoop..."
                scoop install git-lfs

            else
                echo "No Windows package manager found (Chocolatey or Scoop)."
                echo "Install git-lfs manually from https://git-lfs.com/"
                exit 1
            fi
            ;;

        *)
            echo "Unsupported OS: $OS"
            exit 1
            ;;
    esac

    echo "git-lfs installation step completed."
fi

# Download files in git LFS
echo "Downloading CrisperWhisper's files in Git Large File Storage..."
git -C "$CRISPERWHISPER_DIR" lfs install
git -C "$CRISPERWHISPER_DIR" lfs fetch --all
git -C "$CRISPERWHISPER_DIR" lfs pull

echo "Setup complete."
