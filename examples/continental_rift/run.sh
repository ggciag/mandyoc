#!/bin/sh

MANDYOC_PATH = ../..
MPI_PATH=${HOME}/${PETSC_DIR}/${PETSC_PATH}/bin/
NUMBER_OF_CORES=2

MPI_PATH/mpirun -n $NUMBER_OF_CORES $MANDYOC_PATH/mandyoc -seed 0,2 -strain_seed 0.0,1.0 
