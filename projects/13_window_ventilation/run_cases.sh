#!/usr/bin/env bash
# Requires a sourced OpenFOAM Foundation 6 environment and prepared meshes.
set -euo pipefail
cd "$(dirname "$0")"
if [[ "${WM_PROJECT_VERSION:-}" != 6 ]]; then
    echo "Source OpenFOAM Foundation 6 (WM_PROJECT_VERSION must be 6)." >&2
    exit 1
fi
run() {
    local case_dir=$1 app=$2
    shift 2
    echo "$case_dir: $app"
    if ! "$app" -case "$case_dir" "$@" > "$case_dir/log.$app" 2>&1; then
        tail -n 40 "$case_dir/log.$app" >&2
        return 1
    fi
}
for scenario in wide_open tilted; do
    case_dir="cases/$scenario"
    [[ -f "$case_dir/constant/window.msh" ]] || { echo "Missing mesh; run mesh_cases.sh." >&2; exit 1; }
    run "$case_dir" gmshToFoam "$case_dir/constant/window.msh"
    run "$case_dir" changeDictionary
    run "$case_dir" checkMesh
    grep -q 'Mesh OK\.' "$case_dir/log.checkMesh" || { echo "checkMesh rejected $scenario" >&2; exit 1; }
    run "$case_dir" topoSet
    run "$case_dir" setFields
    run "$case_dir" buoyantBoussinesqPimpleFoam
    touch "$case_dir/$scenario.foam"
done
