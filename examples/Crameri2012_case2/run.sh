#!/bin/sh

MANDYOC_PATH = ../..
MPI_PATH=$HOME/petsc/arch-label-optimized/bin/

MPI_PATH/mpirun -n 4 $MANDYOC_PATH/mandyoc
