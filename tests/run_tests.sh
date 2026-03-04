#!/usr/bin/env bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# scenarios to test
declare -a scenarios=(
    "vanKeken1997_case1a"
    "Crameri2012_case2"
    "continental_rift"
    "punch"
    "Schmalholz_2011"
    "two_layers_with_variable_conductivity"
)

# Delete old test summary
rm -rf "${SCRIPT_DIR}/test_summary.txt"

# Run scenarios
for scenario in "${scenarios[@]}"
do
    echo -e "\n==> Running simulation: ${scenario} <==\n"

    # Remove old output
    rm -rf "${SCRIPT_DIR}/data/${scenario}/output"

    # Copy input files
    cp -r "${SCRIPT_DIR}/data/${scenario}/input" "${SCRIPT_DIR}/data/${scenario}/output"

    # Run the simulation
    cd "${SCRIPT_DIR}/data/${scenario}/output" && \
        MANDYOC=${MANDYOC} bash run.sh
done

echo -e "\n==> Verifying results <==\n"

cd "${SCRIPT_DIR}" && \
    pytest -v testing_results.py --tb=line -rf | tee test_summary.txt
