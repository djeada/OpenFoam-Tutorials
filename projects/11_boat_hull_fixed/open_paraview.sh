#!/bin/bash
# =========================================================================
#  Open the fixed boat hull results in ParaView
#
#  For the prepared view (centre slice, hull, streamlines) use
#  ./open_paraview_streamlines.sh instead.
# =========================================================================
set -euo pipefail
cd "${0%/*}" || exit 1

touch boat_hull_fixed.foam

if ! command -v paraview &>/dev/null; then
    echo "Error: paraview is not installed or not in PATH" >&2
    exit 1
fi

paraview boat_hull_fixed.foam
