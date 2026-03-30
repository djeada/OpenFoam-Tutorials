#!/bin/bash
set -euo pipefail

cd "${0%/*}" || exit
paraview --script=paraview_streamlines.py
