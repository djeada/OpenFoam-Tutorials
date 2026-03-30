#!/bin/bash
# =========================================================================
#  Open ParaView for the fixed boat hull case
#
#  Launches ParaView with the streamlines visualization script.
#  For a plain view without streamlines, run:
#    paraview boat_hull_fixed.foam
# =========================================================================
set -euo pipefail
cd "${0%/*}" || exit 1

touch boat_hull_fixed.foam

if ! command -v paraview &>/dev/null; then
    echo "Error: paraview is not installed or not in PATH" >&2
    exit 1
fi

paraview --script paraview_streamlines.py
