#!/bin/bash
set -eu

SCRIPT_DIR="$(CDPATH= cd -- "$(dirname "$0")" && pwd)"
CASE_DIR="$(CDPATH= cd -- "$SCRIPT_DIR/.." && pwd)"

T=0.12
C=1.0
N_POINTS=200
THICKNESS=0.1
VTK_FILENAME="${SCRIPT_DIR}/naca0012.vtk"
STL_FILENAME="${SCRIPT_DIR}/naca0012.stl"
REPAIRED_STL="${SCRIPT_DIR}/naca0012_repaired.stl"
MESHED_STL="${SCRIPT_DIR}/naca0012_meshed.stl"
CASE_STL="${CASE_DIR}/constant/triSurface/airfoil.stl"

check_status() {
    if [ $? -ne 0 ]; then
        echo "Error: $1 failed."
        exit 1
    fi
}

echo "Generating NACA 0012 geometry in ${SCRIPT_DIR}"
python3 "${SCRIPT_DIR}/generate_naca_0012_airfoil.py" \
    --t "$T" \
    --c "$C" \
    --n_points "$N_POINTS" \
    --thickness "$THICKNESS" \
    --vtk_filename "$VTK_FILENAME" \
    --stl_filename "$STL_FILENAME"
check_status "geometry generation"

echo "Repairing STL in Blender"
blender --factory-startup --background \
    --python "${SCRIPT_DIR}/repair_naca_0012_airfoil.py" -- \
    --input_stl "$STL_FILENAME" \
    --repaired_stl "$REPAIRED_STL"
check_status "STL repair"

echo "Triangulating STL in Blender"
blender --factory-startup --background \
    --python "${SCRIPT_DIR}/mesh_naca_0012_airfoil.py" -- \
    --input_stl "$REPAIRED_STL" \
    --output_stl "$MESHED_STL"
check_status "STL meshing"

cp "$MESHED_STL" "$CASE_STL"
echo "Updated case geometry: ${CASE_STL}"
