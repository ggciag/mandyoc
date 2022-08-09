/Users/victorsacek/Documents/petsc/arch-label-optimized/bin/mpirun -n 8 /Users/victorsacek/Documents/gits/mandyoc/bin/mandyoc \
-veloc_ksp_type fgmres \
-veloc_pc_type mg \
-veloc_pc_mg_galerkin \
-veloc_pc_mg_levels 4 \
-veloc_mg_levels_ksp_max_it 9 \
 -seed 0,2 -strain_seed 0.0,1.0 

