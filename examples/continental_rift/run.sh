#!/bin/sh

MANDYOC_PATH = ../..
MPI_PATH=$HOME/petsc/arch-label-optimized/bin/

MPI_PATH/mpirun -n 8 $MANDYOC_PATH/mandyoc -seed 0,2 -strain_seed 0.0,1.0 
