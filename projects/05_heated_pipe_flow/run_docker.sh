#!/bin/bash
# =========================================================================
#  Heated Pipe Flow - Docker runner
#
#  Cleans the case, then runs Allrun inside the cfdengine/openfoam image
#  (OpenFOAM 6) with this directory mounted at /case.
#
#  Usage:  ./run_docker.sh
# =========================================================================
set -euo pipefail
cd "${0%/*}" || exit 1

IMAGE="cfdengine/openfoam"

docker run --rm -v "$PWD":/case -w /case "$IMAGE" \
    bash -lc 'source /opt/openfoam6/etc/bashrc && ./Allclean && ./Allrun'

# Check the results (reads log.* and the time directories; no Docker needed)
python3 verify.py
