#!/bin/bash
# Devcontainer setup script

set -e

echo "Installing additional dependencies..."

# Install f2c for Fortran to C conversion (needed for ARPACK)
sudo apt-get update
sudo apt-get install -y f2c libf2c2-dev

# Install Python numpy (for testing/verification)
pip install numpy

# Install Node.js dependencies
pnpm install

echo "Setup complete!"
