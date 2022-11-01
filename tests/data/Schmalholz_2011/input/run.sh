#!/bin/sh

# exit when any command fails
set -e

${MPIEXEC} -n 2 ${MANDYOC} -seed 0 -strain_seed 0.0
