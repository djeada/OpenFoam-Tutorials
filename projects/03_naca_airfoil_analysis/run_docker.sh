#!/bin/bash
# =========================================================================
#  NACA 0012 — unified Docker runner
#
#  Usage:
#    ./run_docker.sh              Run the steady-state case (default)
#    ./run_docker.sh steady       Run the steady-state RANS case
#    ./run_docker.sh transient    Run the transient animation case
#    ./run_docker.sh all          Run both cases sequentially
# =========================================================================
set -eu

cd "${0%/*}" || exit 1

CASE="${1:-steady}"
IMAGE="cfdengine/openfoam"

usage() {
    cat <<EOF
Usage: $(basename "$0") [CASE]

Cases:
  steady      Steady-state RANS with simpleFoam (default)
  transient   Transient laminar animation with pimpleFoam
  all         Run both cases sequentially
EOF
}

run_steady() {
    echo "========================================"
    echo "  NACA 0012 — Steady-State (simpleFoam)"
    echo "========================================"

    docker run --rm -v "$PWD":/case -w /case "$IMAGE" \
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

    docker run --rm -v "$PWD":/case -w /case "$IMAGE" \
        bash -lc 'source /opt/openfoam6/etc/bashrc && ./Allrun'

    echo "Steady-state case complete."
}

run_transient() {
    echo "========================================"
    echo "  NACA 0012 — Transient (pimpleFoam)"
    echo "========================================"

    docker run --rm -v "$PWD/transient_animation":/case -w /case "$IMAGE" \
        bash -lc 'source /opt/openfoam6/etc/bashrc && ./Allclean'

    python3 mesh_generation_scripts/generate_naca0012_benchmark_geo.py \
        --output transient_animation/constant/airfoil_animation.geo \
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

    gmsh -3 -format msh2 \
        transient_animation/constant/airfoil_animation.geo \
        -o transient_animation/constant/airfoil_animation.msh

    docker run --rm -v "$PWD/transient_animation":/case -w /case "$IMAGE" \
        bash -lc 'source /opt/openfoam6/etc/bashrc && ./Allrun'

    echo "Transient case complete."
}

case "$CASE" in
    steady)     run_steady ;;
    transient)  run_transient ;;
    all)        run_steady; run_transient ;;
    -h|--help|help) usage ;;
    *)
        echo "Unknown case: $CASE" >&2
        usage
        exit 1
        ;;
esac
