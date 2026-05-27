#!/usr/bin/env bash

THIS_DIR=$(dirname -- "${BASH_SOURCE[0]}")
THIS_SCRIPT=$(basename -- ${BASH_SOURCE[0]})
THIS_APP_AND_VERSION=$(cd -- "$THIS_DIR" &> /dev/null && pwd | rev | cut -d'/' -f2,3 | rev )

if [[ ! -s "$PATHOGEN_APPLICATIONS" ]]; then
    echo "Environment variable PATHOGEN_APPLICATIONS is undefined or refers to a file that does not exist" 1>&2
    exit 255
fi

# comment below stops shellcheck complaining that it doesn't know where this file is
# shellcheck disable=SC1090
source "${PATHOGEN_APPLICATIONS}"

track_usage "$THIS_APP_AND_VERSION" "$THIS_SCRIPT" "${@:1:25}"

nextflow run ${THIS_DIR}/main.nf "${@:1}"
