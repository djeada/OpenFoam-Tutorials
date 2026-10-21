#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
command -v docker >/dev/null || { echo "Docker is required." >&2; exit 1; }
./mesh_cases.sh "$@"
docker run --rm --user "$(id -u):$(id -g)" -v "$PWD":/case -w /case \
    cfdengine/openfoam bash -lc 'source /opt/openfoam6/etc/bashrc && ./run_cases.sh'
python3 verify.py
python3 analyze.py
