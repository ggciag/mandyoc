#!/bin/sh

# exit when any command fails
set -e

${MPIEXEC} -n 2 ${MANDYOC}
