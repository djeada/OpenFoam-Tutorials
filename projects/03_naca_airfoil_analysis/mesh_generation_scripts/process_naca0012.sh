#!/bin/bash

# Script name
SCRIPT_NAME="NACA 0012 Airfoil Processing Script"

# Default parameters
T=0.12
C=1.0
N_POINTS=100
THICKNESS=0.1
VTK_FILENAME="naca0012.vtk"
STL_FILENAME="naca0012.stl"
REPAIRED_STL="naca0012_repaired.stl"
MESHED_STL="naca0012_meshed.stl"

# Helper function to check if previous command was successful
check_status() {
    if [ $? -ne 0 ]; then
        echo "Error: $1 failed. Exiting script."
        exit 1
    else
        echo "$1 completed successfully."
    fi
}

echo "Starting $SCRIPT_NAME"

# Step 1: Generate NACA 0012 Airfoil
echo "Step 1: Generating NACA 0012 Airfoil"
python generate_naca_0012_airfoil.py --t $T --c $C --n_points $N_POINTS --thickness $THICKNESS --vtk_filename $VTK_FILENAME --stl_filename $STL_FILENAME
check_status "NACA 0012 Airfoil Generation"

# Validate STL file creation
if [ ! -f $STL_FILENAME ]; then
    echo "Error: STL file $STL_FILENAME not found. Exiting script."
    exit 1
fi

# Step 2: Repair NACA 0012 Airfoil
echo "Step 2: Repairing NACA 0012 Airfoil"
blender --background --python repair_naca_0012_airfoil.py -- --input_stl $STL_FILENAME --repaired_stl $REPAIRED_STL
check_status "NACA 0012 Airfoil Repair"

# Validate repaired STL file creation
if [ ! -f $REPAIRED_STL ]; then
    echo "Error: Repaired STL file $REPAIRED_STL not found. Exiting script."
    exit 1
fi

# Step 3: Mesh NACA 0012 Airfoil
echo "Step 3: Meshing NACA 0012 Airfoil"
blender --background --python mesh_naca_0012_airfoil.py -- --input_stl $REPAIRED_STL --output_stl $MESHED_STL
check_status "NACA 0012 Airfoil Meshing"

# Validate meshed STL file creation
if [ ! -f $MESHED_STL ]; then
    echo "Error: Meshed STL file $MESHED_STL not found. Exiting script."
    exit 1
fi

echo "All steps completed successfully. Meshed STL file saved as $MESHED_STL."
