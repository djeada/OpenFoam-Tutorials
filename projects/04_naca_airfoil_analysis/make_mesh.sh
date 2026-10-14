#!/bin/bash
# Generate the Gmsh geometry and mesh for the steady NACA 0012 case.
# Needs python3 and gmsh on PATH. An existing constant/airfoil_benchmark.msh is reused.
set -eu

cd "${0%/*}" || exit 1

GEO=constant/airfoil_benchmark.geo
MSH=constant/airfoil_benchmark.msh

# --angle-deg -4 rotates the airfoil clockwise about the quarter chord, so the
# horizontal freestream sees a +4 degree angle of attack.
echo "Generating $GEO ..."
python3 mesh_generation_scripts/generate_naca0012_benchmark_geo.py \
    --output "$GEO" \
    --angle-deg -4 \
    --chord-points 101 \
    --domain-xmin -15 \
    --domain-xmax 20 \
    --domain-yhalf 15 \
    --mesh-size-airfoil 0.012 \
    --first-layer 0.004 \
    --bl-thickness 0.05 \
    --bl-ratio 1.12 \
    --near-size 0.03 \
    --near-dist 1.5 \
    --wake-size 0.04 \
    --wake-xmax 10 \
    --wake-yhalf 3 \
    --split-farfield

if [ -f "$MSH" ]; then
    echo "Reusing existing $MSH"
elif command -v gmsh >/dev/null 2>&1; then
    echo "Meshing with Gmsh (this can take several minutes) ..."
    gmsh -3 -format msh2 "$GEO" -o "$MSH"
else
    echo "gmsh not found and $MSH does not exist." >&2
    echo "Run ./make_mesh.sh on a machine with gmsh, or use ./run_docker.sh." >&2
    exit 1
fi
