#!/bin/bash
set -euo pipefail
cd "${0%/*}" || exit 1

echo "=== Cleaning case ==="
docker run --rm -v "$PWD":/case -w /case cfdengine/openfoam \
    bash -lc 'source /opt/openfoam6/etc/bashrc && ./Allclean'

echo "=== Running case ==="
docker run --rm -v "$PWD":/case -w /case cfdengine/openfoam \
    bash -lc 'source /opt/openfoam6/etc/bashrc && ./Allrun'
