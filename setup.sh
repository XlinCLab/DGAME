#!/usr/bin/env bash

# Create Python virtual environment
python3.11 -m venv venv

# Activate virtual environment
source venv/bin/activate

# Install requirements
venv/bin/pip install --upgrade pip
venv/bin/pip install -r requirements.txt
