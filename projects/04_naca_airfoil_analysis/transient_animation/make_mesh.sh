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
echo "Generating $GEO ..."
python3 "$GENERATOR" \
    --output "$GEO" \
    --angle-deg -15 \
    --chord-points 121 \
    --domain-xmin -15 \
    --domain-xmax 25 \
    --domain-yhalf 15 \
    --mesh-size-airfoil 0.008 \
    --mesh-size-farfield 0.8 \
    --first-layer 0.002 \
    --bl-thickness 0.05 \
    --bl-ratio 1.12 \
    --near-size 0.018 \
    --near-dist 2.0 \
    --wake-size 0.025 \
    --wake-xmax 15 \
    --wake-yhalf 4 \
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
