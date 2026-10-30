#!/bin/bash
# =========================================================================
#  Open the results in ParaView (run the case first). Use the Play button:
#  the rotor turns, because the mesh points are written at every time.
# =========================================================================
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

if ! command -v paraview &>/dev/null; then
    echo "Error: paraview is not installed or not in PATH" >&2
    exit 1
fi

touch wind_turbine.foam
paraview wind_turbine.foam
