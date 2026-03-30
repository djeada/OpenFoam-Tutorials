#!/bin/bash
# =========================================================================
#  Open ParaView for the Ahmed body aerodynamics case
# =========================================================================
set -euo pipefail
cd "${0%/*}" || exit 1

touch ahmed_body.foam

if ! command -v paraview &>/dev/null; then
    echo "Error: paraview is not installed or not in PATH" >&2
    exit 1
fi

paraview ahmed_body.foam
