#!/usr/bin/env bash

THIS_DIR=$(dirname -- "${BASH_SOURCE[0]}")
THIS_SCRIPT=$(basename -- ${BASH_SOURCE[0]})
THIS_APP_AND_VERSION=$(cd -- "$THIS_DIR" &> /dev/null && pwd | rev | cut -d'/' -f2,3 | rev )

source "${PATHOGEN_APPLICATIONS}"

track_usage "$THIS_APP_AND_VERSION" "$THIS_SCRIPT" "${@:1:25}"

nextflow run ${THIS_DIR}/main.nf ${@:1}
