#!/bin/sh

# exit when any command fails
set -e

${MPIEXEC} -n 1 ${MANDYOC}
