#!/bin/bash
# =========================================================================
#  Open ParaView for the wind turbine MRF simulation
# =========================================================================
set -euo pipefail
cd "${0%/*}" || exit 1

touch wind_turbine.foam

if ! command -v paraview &>/dev/null; then
    echo "Error: paraview is not installed or not in PATH" >&2
    exit 1
fi

paraview wind_turbine.foam
