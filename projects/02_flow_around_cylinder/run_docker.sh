#!/bin/bash
# =========================================================================
#  2D Flow Around Cylinder — Docker runner
#
#  Runs the full case inside the cfdengine/openfoam Docker image.
#  Gmsh and Python run on the host; OpenFOAM runs in the container.
#
#  Usage:
#    ./run_docker.sh
# =========================================================================
set -eu

cd "${0%/*}" || exit 1

IMAGE="cfdengine/openfoam"

echo "========================================"
echo "  2D Cylinder — Docker Runner"
echo "========================================"

# ── Clean ─────────────────────────────────────────────────────────────────
echo "Cleaning previous run..."
docker run --rm -v "$PWD":/case -w /case "$IMAGE" \
    bash -lc 'source /opt/openfoam6/etc/bashrc && ./Allclean'

# ── Generate geometry & mesh on host ──────────────────────────────────────
echo "Generating Gmsh geometry..."
python3 generate_cylinder_geo.py --output constant/cylinder.geo

echo "Meshing with Gmsh..."
gmsh -3 -format msh2 constant/cylinder.geo -o constant/cylinder.msh

# ── Run OpenFOAM in Docker ───────────────────────────────────────────────
echo "Running OpenFOAM case in Docker..."
docker run --rm -v "$PWD":/case -w /case "$IMAGE" \
    bash -lc 'source /opt/openfoam6/etc/bashrc && ./Allrun'

echo "Docker run complete."
