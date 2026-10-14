#!/bin/bash
# =========================================================================
#  Wind Turbine Blade (MRF) - Docker runner
#
#  Cleans and runs the case inside the cfdengine/openfoam image
#  (OpenFOAM 6). The case directory is mounted at /case.
#
#  Usage:  ./run_docker.sh
# =========================================================================
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

IMAGE="cfdengine/openfoam"
SRC="source /opt/openfoam6/etc/bashrc"

docker run --rm -v "$PWD":/case -w /case "$IMAGE" bash -lc "$SRC && ./Allclean"
docker run --rm -v "$PWD":/case -w /case "$IMAGE" bash -lc "$SRC && ./Allrun"

# Check the results (reads logs and fields only; needs python3 on the host)
python3 verify.py
