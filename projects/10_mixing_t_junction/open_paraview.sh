#!/bin/bash
# =========================================================================
#  Open the T-junction case in ParaView (creates an empty t_junction.foam
#  file as the reader entry point; Allclean removes it).
# =========================================================================
set -euo pipefail
cd "${0%/*}" || exit 1

if ! command -v paraview &>/dev/null; then
    echo "Error: paraview is not installed or not in PATH" >&2
    exit 1
fi

touch t_junction.foam
paraview t_junction.foam
