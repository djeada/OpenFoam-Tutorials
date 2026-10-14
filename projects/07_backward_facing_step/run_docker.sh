#!/bin/bash
# =========================================================================
#  Backward-facing step - Docker runner
#
#  Cleans the case, runs Allrun in the cfdengine/openfoam image
#  (OpenFOAM 6, case directory mounted at /case) and checks the results
#  with verify.py on the host.
#
#  Usage:  ./run_docker.sh
# =========================================================================
set -euo pipefail
cd "${0%/*}" || exit 1

IMAGE="cfdengine/openfoam"
# Run as the calling user so generated files are not owned by root.
DOCKER_RUN=(docker run --rm --user "$(id -u):$(id -g)" -v "$PWD":/case -w /case "$IMAGE")

echo "Cleaning previous run..."
"${DOCKER_RUN[@]}" bash -lc 'source /opt/openfoam6/etc/bashrc && ./Allclean'

echo "Running OpenFOAM case in Docker..."
"${DOCKER_RUN[@]}" bash -lc 'source /opt/openfoam6/etc/bashrc && ./Allrun'

echo "Verifying results..."
python3 verify.py
