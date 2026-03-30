#!/bin/bash
set -euo pipefail
cd "${0%/*}" || exit 1

touch elbow.foam
# check paraview is available
paraview --script paraview_elbow.py
