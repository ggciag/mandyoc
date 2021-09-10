#!/bin/sh

MPI_PATH=$HOME/petsc/arch-label-optimized/bin
MANDYOC_PATH=../..
NUMBER_OF_CORES=2

$MPI_PATH/mpirun -n $NUMBER_OF_CORES $MANDYOC_PATH/mandyoc