#!/bin/bash
# =========================================================================
#  NACA 0012 -- Docker runner
#
#  Python and Gmsh run on the host (the OpenFOAM image has no Gmsh);
#  OpenFOAM 6 runs in the cfdengine/openfoam container.
#
#  Requires: docker, python3, gmsh
#
#  Usage:
#    ./run_docker.sh              Steady-state RANS case (default), then verify.py
#    ./run_docker.sh steady       same
#    ./run_docker.sh transient    Transient animation case (pimpleFoam)
#    ./run_docker.sh all          Both cases, one after the other
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

check_tools() {
    for tool in docker python3 gmsh; do
        if ! command -v "$tool" >/dev/null 2>&1; then
            echo "Error: $tool is required but not found in PATH" >&2
            exit 1
        fi
    done
}

# docker_foam WORKDIR COMMAND: run COMMAND in the OpenFOAM container with the
# project mounted at /case, as the calling user so no root-owned files appear.
docker_foam() {
    docker run --rm --user "$(id -u):$(id -g)" -v "$PWD":/case -w "/case/$1" "$IMAGE" \
        bash -lc "source /opt/openfoam6/etc/bashrc && $2"
}

run_steady() {
    echo "========================================"
    echo "  NACA 0012 -- Steady-State (simpleFoam)"
    echo "========================================"

    docker_foam . ./Allclean
    ./make_mesh.sh
    docker_foam . ./Allrun

    echo "Steady-state case complete."
    python3 verify.py
}

run_transient() {
    echo "========================================"
    echo "  NACA 0012 -- Transient (pimpleFoam)"
    echo "========================================"

    docker_foam transient_animation ./Allclean
    ./transient_animation/make_mesh.sh
    docker_foam transient_animation ./Allrun

    echo "Transient case complete."
}

case "$CASE" in
    steady)     check_tools; run_steady ;;
    transient)  check_tools; run_transient ;;
    all)        check_tools; run_steady; run_transient ;;
    -h|--help|help) usage ;;
    *)
        echo "Unknown case: $CASE" >&2
        usage
        exit 1
        ;;
esac
