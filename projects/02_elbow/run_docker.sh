#!/bin/bash
set -eu
cd "${0%/*}" || exit 1

docker run --rm -v "$PWD":/case -w /case cfdengine/openfoam bash -lc 'source /opt/openfoam6/etc/bashrc && ./allclean'
docker run --rm -v "$PWD":/case -w /case cfdengine/openfoam bash -lc 'source /opt/openfoam6/etc/bashrc && ./allrun'
