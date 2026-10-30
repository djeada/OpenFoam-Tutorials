#!/bin/bash
# =========================================================================
#  Vertical-axis wind turbine (rotating AMI mesh) - Docker runner
#
#  Python and Gmsh run on the host (the image has no Gmsh); OpenFOAM 6 runs
#  in the cfdengine/openfoam container with the case mounted at /case.
#
#  Requires: docker, python3, gmsh (tested with 4.12.1)
#  Ends with `python3 verify.py` (PASS/FAIL checks of mesh, AMI, torque, power).
#
#  Usage:   ./run_docker.sh             (all cores, at most 8)
#           NPROCS=4 ./run_docker.sh
#  Run time: about 70 minutes on 8 cores.
# =========================================================================
set -euo pipefail
cd "${0%/*}" || exit 1

IMAGE="cfdengine/openfoam"

for tool in docker python3 gmsh; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "Error: $tool is required but not found in PATH" >&2
        exit 1
    fi
done

# Run as the calling user so generated files are not owned by root.
docker_foam() {
    docker run --rm --user "$(id -u):$(id -g)" -e NPROCS="${NPROCS:-}" \
        -v "$PWD":/case -w /case "$IMAGE" \
        bash -lc "source /opt/openfoam6/etc/bashrc && $1"
}

echo "Cleaning previous run..."
docker_foam ./Allclean

echo "Meshing with Gmsh on the host..."
python3 generate_vawt_geo.py --output constant/vawt.geo
gmsh -3 -format msh2 constant/vawt.geo -o constant/vawt.msh > log.gmsh 2>&1

echo "Running OpenFOAM in Docker..."
docker_foam ./Allrun

echo "Verifying results..."
python3 verify.py
