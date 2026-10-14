#!/bin/bash
# =========================================================================
#  Fixed Boat Hull — Docker runner
#
#  Runs Allclean and Allrun inside the cfdengine/openfoam image, then checks
#  the results with verify.py on the host.
#
#  The repository root is mounted (not just this directory) because Allrun
#  reads ../../boat_hull.stl. The container runs as the calling user, so no
#  root-owned files are left in the repository.
#
#  Usage:
#    ./run_docker.sh
# =========================================================================
set -euo pipefail
cd "${0%/*}" || exit 1

IMAGE="cfdengine/openfoam"
REPO_ROOT="$(cd ../.. && pwd)"
PROJECT="projects/11_boat_hull_fixed"

echo "========================================"
echo "  Fixed Boat Hull — Docker Runner"
echo "========================================"

echo "Cleaning previous run..."
docker run --rm -u "$(id -u):$(id -g)" -v "$REPO_ROOT":/repo -w /repo/"$PROJECT" "$IMAGE" \
    bash -lc 'source /opt/openfoam6/etc/bashrc && ./Allclean'

echo "Running OpenFOAM case in Docker..."
docker run --rm -u "$(id -u):$(id -g)" -v "$REPO_ROOT":/repo -w /repo/"$PROJECT" "$IMAGE" \
    bash -lc 'source /opt/openfoam6/etc/bashrc && ./Allrun'

echo "Docker run complete. Checking results..."
python3 verify.py
