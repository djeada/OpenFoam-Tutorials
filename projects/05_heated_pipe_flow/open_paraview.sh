#!/bin/bash
# =========================================================================
#  Open ParaView for the heated pipe flow case
# =========================================================================
set -euo pipefail
cd "${0%/*}" || exit 1

touch heated_pipe.foam

if ! command -v paraview &>/dev/null; then
    echo "Error: paraview is not installed or not in PATH" >&2
    exit 1
fi

paraview heated_pipe.foam
