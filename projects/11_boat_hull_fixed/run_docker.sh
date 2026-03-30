#!/bin/bash
# =========================================================================
#  Fixed Boat Hull — Docker runner
#
#  The Allrun script references ../../boat_hull.stl, so we mount the
#  repository root to preserve relative paths inside the container.
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
docker run --rm -v "$REPO_ROOT":/repo -w /repo/"$PROJECT" "$IMAGE" \
    bash -lc 'source /opt/openfoam6/etc/bashrc && ./Allclean'

echo "Running OpenFOAM case in Docker..."
docker run --rm -v "$REPO_ROOT":/repo -w /repo/"$PROJECT" "$IMAGE" \
    bash -lc 'source /opt/openfoam6/etc/bashrc && ./Allrun'

echo "Docker run complete."
