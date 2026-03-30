#!/bin/bash
# =========================================================================
#  Ahmed Body Aerodynamics — Docker runner
#
#  Generates the Ahmed body STL on the host with Python, then runs
#  the parallel OpenFOAM case inside the cfdengine/openfoam Docker image.
#
#  Usage:
#    ./run_docker.sh               25° slant angle (default)
#    ./run_docker.sh --slant 35    35° slant angle
# =========================================================================
set -euo pipefail
cd "${0%/*}" || exit 1

IMAGE="cfdengine/openfoam"
SLANT_ANGLE="${1:-25}"

# Accept --slant N syntax
if [ "${1:-}" = "--slant" ]; then
    SLANT_ANGLE="${2:?Usage: $0 [--slant ANGLE]}"
fi

echo "========================================"
echo "  Ahmed Body — Docker Runner"
echo "  Slant angle: ${SLANT_ANGLE}°"
echo "========================================"

# ── Clean ─────────────────────────────────────────────────────────────────
echo "Cleaning previous run..."
docker run --rm -v "$PWD":/case -w /case "$IMAGE" \
    bash -lc 'source /opt/openfoam6/etc/bashrc && ./Allclean'

# ── Generate STL on host ──────────────────────────────────────────────────
echo "Generating Ahmed body STL (${SLANT_ANGLE}° slant)..."
python3 generate_ahmed_body.py --slant-angle "$SLANT_ANGLE"

# ── Run OpenFOAM in Docker ───────────────────────────────────────────────
echo "Running OpenFOAM case in Docker..."
docker run --rm -v "$PWD":/case -w /case "$IMAGE" \
    bash -lc 'source /opt/openfoam6/etc/bashrc && ./Allrun'

echo "Docker run complete."
