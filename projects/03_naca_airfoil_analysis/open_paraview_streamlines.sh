#!/bin/bash
set -euo pipefail

cd "${0%/*}" || exit
touch airfoil.foam
paraview --script=paraview_airfoil.py -- "$@"
