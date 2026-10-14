#!/bin/bash
# Open the case in ParaView (creates an empty dam_break.foam marker file).
# Run ./Allrun or ./run_docker.sh first.
set -euo pipefail
cd "${0%/*}" || exit 1

if ! command -v paraview &>/dev/null; then
    echo "Error: paraview is not installed or not in PATH" >&2
    exit 1
fi

touch dam_break.foam
paraview dam_break.foam
