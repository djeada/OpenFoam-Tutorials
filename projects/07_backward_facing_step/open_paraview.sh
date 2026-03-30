#!/bin/bash
# =========================================================================
#  Open ParaView for the backward-facing step case
# =========================================================================
set -euo pipefail
cd "${0%/*}" || exit 1

touch backward_step.foam

if ! command -v paraview &>/dev/null; then
    echo "Error: paraview is not installed or not in PATH" >&2
    exit 1
fi

paraview backward_step.foam
