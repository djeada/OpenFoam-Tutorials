#!/bin/bash
# =========================================================================
#  Open ParaView with the prepared view: centre slice (y = 0), hull surface
#  and streamlines. The view is built by paraview_streamlines.py.
# =========================================================================
set -euo pipefail
cd "${0%/*}" || exit 1

touch boat_hull_fixed.foam

if ! command -v paraview &>/dev/null; then
    echo "Error: paraview is not installed or not in PATH" >&2
    exit 1
fi

paraview --script=paraview_streamlines.py
