#!/usr/bin/env bash
# Mesh on the host, matching the Gmsh workflow in projects 03 and 04.
set -euo pipefail
cd "$(dirname "$0")"
command -v gmsh >/dev/null || { echo "Install Gmsh 4.x (tested with 4.12 and 4.15) and put gmsh on PATH." >&2; exit 1; }
python3 prepare_cases.py "$@"
for scenario in wide_open tilted; do
    case_dir="cases/$scenario"
    echo "Meshing $scenario (see $case_dir/log.gmsh)"
    gmsh "$case_dir/constant/window.geo" -3 -format msh2 \
        -o "$case_dir/constant/window.msh" > "$case_dir/log.gmsh" 2>&1
done
