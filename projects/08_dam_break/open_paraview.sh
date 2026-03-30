#!/bin/bash
# =========================================================================
#  Open ParaView for the dam break simulation
# =========================================================================
set -euo pipefail
cd "${0%/*}" || exit 1

touch dam_break.foam

if ! command -v paraview &>/dev/null; then
    echo "Error: paraview is not installed or not in PATH" >&2
    exit 1
fi

paraview dam_break.foam
