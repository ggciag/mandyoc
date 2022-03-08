#!/usr/bin/env bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# scenarios to test
declare -a scenarios=(
    "vanKeken1997_case1a"
    "Crameri2012_case2"
    "continental_rift"
)

for scenario in "${scenarios[@]}"
do
    echo -e "\nRunning ${scenario}\n"

    cp -r "${SCRIPT_DIR}/data/${scenario}/input" "${SCRIPT_DIR}/data/${scenario}/output"

    cd "${SCRIPT_DIR}/data/${scenario}/output" && \
        bash run.sh &&

    cd "${SCRIPT_DIR}/data/${scenario}/" && \
        pytest -v testing_results.py
done
