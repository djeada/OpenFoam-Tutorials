#!/bin/bash
# =========================================================================
#  Floating Boat Hull - Docker runner (image cfdengine/openfoam, OpenFOAM 6)
#
#  Allrun reads ../../boat_hull.stl, so the repository root is mounted and
#  the working directory is the project inside it.
#
#  Usage:  ./run_docker.sh
#  The case needs 1-2 hours or more for the full 3 s; stop it with Ctrl-C
#  or set "stopAt writeNow;" in system/controlDict (runTimeModifiable).
# =========================================================================
set -euo pipefail
cd "${0%/*}" || exit 1

IMAGE="cfdengine/openfoam"
REPO_ROOT="$(cd ../.. && pwd)"
PROJECT="projects/12_boat_hull_floating"

echo "Cleaning previous run..."
docker run --rm -v "$REPO_ROOT":/repo -w /repo/"$PROJECT" "$IMAGE" \
    bash -lc 'source /opt/openfoam6/etc/bashrc && ./Allclean'

echo "Running OpenFOAM case in Docker..."
docker run --rm -v "$REPO_ROOT":/repo -w /repo/"$PROJECT" "$IMAGE" \
    bash -lc 'source /opt/openfoam6/etc/bashrc && ./Allrun'

echo "Docker run complete. Verifying results..."
python3 verify.py
