#!/bin/bash
# =========================================================================
#  Ahmed body aerodynamics: Docker runner
#
#  1. cleans the case, 2. generates the STL on the host (python3),
#  3. runs Allrun inside the cfdengine/openfoam image (OpenFOAM 6),
#  4. checks the results with verify.py.
#
#  Usage:
#    ./run_docker.sh               25 degree rear slant (default)
#    ./run_docker.sh --slant 35    35 degree rear slant
#    ./run_docker.sh 35            same as above
# =========================================================================
set -euo pipefail
cd "${0%/*}" || exit 1

IMAGE="cfdengine/openfoam"
CONTAINER="ahmed_body_$$"
SLANT_ANGLE="${1:-25}"

if [ "${1:-}" = "--slant" ]; then
    SLANT_ANGLE="${2:?Usage: $0 [--slant ANGLE]}"
fi

# Stop the container if this script is interrupted or killed (e.g. by
# `timeout`). The container runs in the background and is waited for, so that
# bash can run the trap immediately instead of after the container ends.
trap 'docker kill "$CONTAINER" >/dev/null 2>&1 || true; exit 143' INT TERM

foam() {
    docker run --rm --name "$CONTAINER" -v "$PWD":/case -w /case "$IMAGE" \
        bash -lc "source /opt/openfoam6/etc/bashrc && $1" &
    wait $!
}

echo "Ahmed body, slant angle ${SLANT_ANGLE} deg"

echo "Cleaning previous run..."
foam ./Allclean

echo "Generating STL..."
python3 generate_ahmed_body.py --slant-angle "$SLANT_ANGLE"

echo "Running OpenFOAM case in Docker..."
foam ./Allrun

echo "Docker run complete. Verifying results..."
python3 verify.py --slant-angle "$SLANT_ANGLE"
