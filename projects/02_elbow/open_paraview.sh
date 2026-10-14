#!/bin/bash
set -euo pipefail
cd "${0%/*}" || exit 1

if ! command -v paraview >/dev/null 2>&1; then
    echo "Error: paraview not found in PATH" >&2
    exit 1
fi

if [ ! -f log.icoFoam ]; then
    echo "Error: no results found, run ./run_docker.sh (or ./Allrun) first" >&2
    exit 1
fi

touch elbow.foam
paraview --script paraview_elbow.py
