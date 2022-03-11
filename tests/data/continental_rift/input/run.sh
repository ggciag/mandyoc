#!/bin/sh

# exit when any command fails
set -e

${MPIEXEC} -n 2 ${MANDYOC} -seed 2 -strain_seed 1.0
