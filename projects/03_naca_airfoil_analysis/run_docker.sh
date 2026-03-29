#!/bin/bash
set -eu

cd "${0%/*}" || exit 1

docker run --rm \
    -v "$PWD":/case \
    -w /case \
    cfdengine/openfoam \
    bash -lc 'source /opt/openfoam6/etc/bashrc && ./Allclean'

python3 mesh_generation_scripts/generate_naca0012_benchmark_geo.py \
    --output constant/airfoil_benchmark.geo \
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

gmsh -3 -format msh2 constant/airfoil_benchmark.geo -o constant/airfoil_benchmark.msh

docker run --rm \
    -v "$PWD":/case \
    -w /case \
    cfdengine/openfoam \
    bash -lc 'source /opt/openfoam6/etc/bashrc && ./Allrun'
