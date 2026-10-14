#!/bin/bash
# =========================================================================
#  Mixing T-Junction - Docker runner
#
#  Cleans the case, runs Allrun inside the cfdengine/openfoam image
#  (OpenFOAM 6, case directory mounted at /case), then checks the results
#  on the host with verify.py (exit code 0 only if all checks pass).
#
#  Usage:  ./run_docker.sh
# =========================================================================
set -euo pipefail
cd "${0%/*}" || exit 1

IMAGE="cfdengine/openfoam"
SETUP='source /opt/openfoam6/etc/bashrc'

echo "Cleaning previous run..."
docker run --rm -v "$PWD":/case -w /case "$IMAGE" \
    bash -lc "$SETUP && ./Allclean"

echo "Running case in Docker..."
docker run --rm -v "$PWD":/case -w /case "$IMAGE" \
    bash -lc "$SETUP && ./Allrun"

echo "Verifying results (python3 verify.py)..."
python3 verify.py
