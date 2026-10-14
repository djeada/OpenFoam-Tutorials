#!/bin/bash
# Run the dam break case inside the cfdengine/openfoam Docker image
# (OpenFOAM 6). The case directory is mounted at /case.
#
# Usage: ./run_docker.sh
set -euo pipefail
cd "${0%/*}" || exit 1

IMAGE="cfdengine/openfoam"

echo "Cleaning previous run..."
docker run --rm -v "$PWD":/case -w /case "$IMAGE" \
    bash -lc 'source /opt/openfoam6/etc/bashrc && ./Allclean'

echo "Running case in Docker..."
docker run --rm -v "$PWD":/case -w /case "$IMAGE" \
    bash -lc 'source /opt/openfoam6/etc/bashrc && ./Allrun'

echo "Verifying results..."
python3 verify.py
