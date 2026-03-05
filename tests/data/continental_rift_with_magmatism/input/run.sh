#!/bin/sh

# exit when any command fails
set -e

${MPIEXEC} -n 2 ${MANDYOC} -seed 0,2 -strain_seed 0.0,1.0
