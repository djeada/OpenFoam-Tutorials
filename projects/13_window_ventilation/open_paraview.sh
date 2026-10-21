#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
for scenario in wide_open tilted; do
    [[ -d "cases/$scenario/constant/polyMesh" ]] || { echo "Run the simulations first." >&2; exit 1; }
    touch "cases/$scenario/$scenario.foam"
done
echo "Add cases/tilted/tilted.foam with File > Open for the comparison."
exec paraview --data="cases/wide_open/wide_open.foam"
