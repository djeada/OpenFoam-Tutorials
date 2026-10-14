#!/bin/bash
# =========================================================================
#  Open the backward-facing step case in ParaView
#  (creates the empty marker file backward_step.foam, removed by Allclean)
# =========================================================================
set -euo pipefail
cd "${0%/*}" || exit 1

if ! command -v paraview &>/dev/null; then
    echo "Error: paraview is not installed or not in PATH" >&2
    exit 1
fi

touch backward_step.foam
paraview backward_step.foam
