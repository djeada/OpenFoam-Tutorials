#!/bin/bash
# =========================================================================
#  Heated Pipe Flow — Docker runner
#
#  Runs the full case inside the cfdengine/openfoam Docker image.
#
#  Usage:
#    ./run_docker.sh
# =========================================================================
set -euo pipefail
cd "${0%/*}" || exit 1

IMAGE="cfdengine/openfoam"

echo "========================================"
echo "  Heated Pipe Flow — Docker Runner"
echo "========================================"

echo "Cleaning previous run..."
docker run --rm -v "$PWD":/case -w /case "$IMAGE" \
    bash -lc 'source /opt/openfoam6/etc/bashrc && ./Allclean'

echo "Running OpenFOAM case in Docker..."
docker run --rm -v "$PWD":/case -w /case "$IMAGE" \
    bash -lc 'source /opt/openfoam6/etc/bashrc && ./Allrun'

echo "Docker run complete."
