#!/usr/bin/env bash

# exit when any command fails
set -e

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# scenarios to test
declare -a scenarios=(
    "vanKeken1997_case1a"
    "Crameri2012_case2"
    "continental_rift"
)

# Run scenarios
for scenario in "${scenarios[@]}"
do
    echo -e "/n==> Running scenario: ${scenario} <===/n"

    cp -r "${SCRIPT_DIR}/data/${scenario}/input" "${SCRIPT_DIR}/data/${scenario}/output"

    cd "${SCRIPT_DIR}/data/${scenario}/output" && \
        MANDYOC=${MANDYOC} bash run.sh
done

# Run tests
cd "${SCRIPT_DIR}" && \
    pytest -v testing_results.py

# Clean up
for scenario in "${scenarios[@]}"
do
    rm -rf "${SCRIPT_DIR}/data/${scenario}/output"
done
