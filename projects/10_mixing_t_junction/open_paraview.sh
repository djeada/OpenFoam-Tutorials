#!/bin/bash
# =========================================================================
#  Open ParaView for the T-junction mixing simulation
# =========================================================================
set -euo pipefail
cd "${0%/*}" || exit 1

touch t_junction.foam

if ! command -v paraview &>/dev/null; then
    echo "Error: paraview is not installed or not in PATH" >&2
    exit 1
fi

paraview t_junction.foam
