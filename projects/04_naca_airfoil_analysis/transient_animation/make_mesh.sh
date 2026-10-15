#!/bin/bash
# Generate the Gmsh geometry and mesh for the transient NACA 0012 case.
# Needs python3 and gmsh on PATH. An existing constant/airfoil_animation.msh is reused.
set -eu

cd "${0%/*}" || exit 1

GEO=constant/airfoil_animation.geo
MSH=constant/airfoil_animation.msh
GENERATOR=../mesh_generation_scripts/generate_naca0012_benchmark_geo.py

# --angle-deg -15 rotates the airfoil clockwise about the quarter chord, so the
# horizontal freestream sees a +15 degree angle of attack.
# The mesh is deliberately moderate (a few tens of thousands of cells): an earlier,
# finer setting (121 chord points, airfoil size 0.008, wake size 0.025, far field 0.8)
# made Gmsh's 2D recombination (Blossom, about 650 000 edges) run for more than 25
# minutes without finishing, and would have made the flow solution take hours.
echo "Generating $GEO ..."
python3 "$GENERATOR" \
    --output "$GEO" \
    --angle-deg -15 \
    --chord-points 101 \
    --domain-xmin -15 \
    --domain-xmax 25 \
    --domain-yhalf 15 \
    --mesh-size-airfoil 0.012 \
    --mesh-size-farfield 1.5 \
    --first-layer 0.003 \
    --bl-thickness 0.05 \
    --bl-ratio 1.12 \
    --near-size 0.03 \
    --near-dist 1.5 \
    --wake-size 0.04 \
    --wake-xmax 12 \
    --wake-yhalf 3 \
    --split-farfield

if [ -f "$MSH" ]; then
    echo "Reusing existing $MSH"
elif command -v gmsh >/dev/null 2>&1; then
    echo "Meshing with Gmsh (this can take several minutes) ..."
    gmsh -3 -format msh2 "$GEO" -o "$MSH"
else
    echo "gmsh not found and $MSH does not exist." >&2
    echo "Run ./make_mesh.sh on a machine with gmsh, or use ../run_docker.sh transient." >&2
    exit 1
fi
