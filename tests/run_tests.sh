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

    # Remove old output
    rm -rf "${SCRIPT_DIR}/data/${scenario}/output"

    # Copy input files
    cp -r "${SCRIPT_DIR}/data/${scenario}/input" "${SCRIPT_DIR}/data/${scenario}/output"

    # Run the simulation
    cd "${SCRIPT_DIR}/data/${scenario}/output" && \
        MANDYOC=${MANDYOC} bash run.sh
done

# Run tests
cd "${SCRIPT_DIR}" && \
    pytest -v testing_results.py
