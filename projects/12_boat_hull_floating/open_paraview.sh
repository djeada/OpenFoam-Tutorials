#!/bin/bash
# =========================================================================
#  Open ParaView for the floating boat hull simulation
# =========================================================================
set -euo pipefail
cd "${0%/*}" || exit 1

touch boat_hull_floating.foam

if ! command -v paraview &>/dev/null; then
    echo "Error: paraview is not installed or not in PATH" >&2
    exit 1
fi

paraview boat_hull_floating.foam
