#!/bin/bash
# =========================================================================
#  Open ParaView with the cylinder vortex-shedding visualization
# =========================================================================
set -euo pipefail
cd "${0%/*}" || exit 1

touch cylinder.foam

if ! command -v paraview &>/dev/null; then
    echo "Error: paraview is not installed or not in PATH" >&2
    exit 1
fi

paraview --script paraview_cylinder.py
