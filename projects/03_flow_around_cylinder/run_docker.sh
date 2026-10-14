#!/bin/bash
# =========================================================================
#  2D Flow Around Cylinder — Docker runner
#
#  Python and Gmsh run on the host (the image has neither Gmsh nor a
#  recent Python); OpenFOAM 6 runs in the cfdengine/openfoam container.
#
#  Requires: docker, python3, gmsh
#  Ends with `python3 verify.py` (PASS/FAIL checks of mesh, solver, Cd/Cl/St).
#
#  Usage:
#    ./run_docker.sh
# =========================================================================
set -eu

cd "${0%/*}" || exit 1

IMAGE="cfdengine/openfoam"

for tool in docker python3 gmsh; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "Error: $tool is required but not found in PATH" >&2
        exit 1
    fi
done

# Gmsh 4.12.x (Ubuntu 24.04 apt package) segfaults in the 2D quad mesher on this
# geometry; 4.11.1, 4.13.1 and 4.15.2 mesh it fine.
case "$(gmsh --version 2>&1 | head -n1)" in
    4.12.*)
        echo "Error: Gmsh 4.12.x crashes on this geometry (segmentation fault)." >&2
        echo "Use another version, for example:" >&2
        echo "  python3 -m venv ~/gmsh-venv && ~/gmsh-venv/bin/pip install gmsh==4.15.2" >&2
        echo "  PATH=~/gmsh-venv/bin:\$PATH ./run_docker.sh" >&2
        exit 1
        ;;
esac

# Run as the calling user so generated files are not owned by root.
docker_foam() {
    docker run --rm --user "$(id -u):$(id -g)" -v "$PWD":/case -w /case "$IMAGE" \
        bash -lc "source /opt/openfoam6/etc/bashrc && $1"
}

echo "========================================"
echo "  2D Cylinder — Docker Runner"
echo "========================================"

# ── Clean ─────────────────────────────────────────────────────────────────
echo "Cleaning previous run..."
docker_foam ./Allclean

# ── Generate geometry & mesh on host ──────────────────────────────────────
echo "Generating Gmsh geometry..."
python3 generate_cylinder_geo.py --output constant/cylinder.geo

echo "Meshing with Gmsh..."
gmsh -3 -format msh2 constant/cylinder.geo -o constant/cylinder.msh

# ── Run OpenFOAM in Docker ───────────────────────────────────────────────
echo "Running OpenFOAM case in Docker..."
docker_foam ./Allrun

echo "Docker run complete. Checking results..."

# ── Verify results (host, reads logs and postProcessing only) ─────────────
python3 verify.py
