#!/bin/bash
set -euo pipefail

if command -v julia >/dev/null 2>&1; then
    echo "Julia is already installed."
    exit 0
fi

echo "Julia not found. Installing..."

OS="$(uname -s)"

if [[ "$OS" == "Linux" || "$OS" == "Darwin" ]]; then
    # macOS or Linux
    curl -fsSL https://install.julialang.org | sh

    # Reload shell configuration so that Julia is found in PATH
    echo "Reloading shell configuration..."

    if [[ -n "$BASH_VERSION" ]]; then
        [[ -f "$HOME/.bashrc" ]] && . "$HOME/.bashrc"
        [[ -f "$HOME/.bash_profile" ]] && . "$HOME/.bash_profile"
    elif [[ -n "$ZSH_VERSION" ]]; then
        [[ -f "$HOME/.zshrc" ]] && . "$HOME/.zshrc"
    fi

elif [[ "$OS" == MINGW* || "$OS" == CYGWIN* || "$OS" == MSYS* ]]; then
    # Windows (Git Bash / MSYS / Cygwin)
    winget install --name Julia --id 9NJNWW8PVKMN -e -s msstore

else
    echo "Unsupported OS: $OS"
    return 1
fi

# Re-check installation
if command -v julia >/dev/null 2>&1; then
    echo "Julia installed successfully."
    exit 0
else
    echo "Julia installation may have failed or requires a new shell."
    exit 1
fi
