#include <petscdmda.h>
#include <petscdmshell.h>
#include <petscdmswarm.h>
#include <petsc/private/dmimpl.h>
#include <petscmath.h>
#include <petscviewer.h>

#include <petscksp.h>

PetscErrorCode magmatic_extrusive(PetscReal magma_volume,PetscReal x_center);
PetscErrorCode sp_update_surface_swarm_particles_properties(PetscInt active_layer);


extern DM da_Thermal;
extern DM dms;
extern DM da_Veloc;

extern DM dms_s;
extern PetscInt dms_s_ppe;

extern Vec Phi;

extern double dx_const;
extern double dz_const;

extern PetscReal x_magma_high;
extern PetscReal z_magma_high;

extern PetscInt magmatic_layer;

extern PetscReal previous_magmatic_volume;

PetscErrorCode calc_magmatic_extraction(){

    PetscErrorCode ierr;

    PetscReal total_magmatism;
    PetscReal magmatism_volume;

    VecSum(Phi,&total_magmatism);
    total_magmatism*=dx_const*dz_const;

    magmatism_volume = total_magmatism - previous_magmatic_volume;

    previous_magmatic_volume = total_magmatism;

    PetscPrintf(PETSC_COMM_WORLD,"Volume of magmatism created: %lg m2\n",magmatism_volume);



    if (magmatism_volume>0.0){

        printf("x,z: %lf %lf",x_magma_high,z_magma_high);

        PetscScalar global_x_magma_high;
        PetscScalar global_z_magma_high;
        PetscScalar local_x = x_magma_high;
        PetscScalar local_z = z_magma_high;

        PetscMPIInt rank, max_rank;
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

        /* Step 1: find global maximum of z */
        MPI_Allreduce(&local_z, &global_z_magma_high, 1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD);

        /* Step 2: find which rank owns the max */
        PetscMPIInt candidate = (local_z == global_z_magma_high) ? rank : -1;

        MPI_Allreduce(&candidate, &max_rank, 1, MPI_INT, MPI_MAX, PETSC_COMM_WORLD);

        /* Step 3: broadcast x from max_rank to all processes */
        MPI_Bcast(&local_x, 1, MPI_DOUBLE, max_rank, PETSC_COMM_WORLD);

        /* Step 4: store result */
        global_x_magma_high = local_x;

        x_magma_high = global_x_magma_high;

        PetscPrintf(PETSC_COMM_WORLD,"x_magma_high: %lf\n",x_magma_high);
        PetscPrintf(PETSC_COMM_WORLD,"z_magma_high: %lf\n",z_magma_high);

        magmatic_extrusive(magmatism_volume, local_x);

        ierr = sp_update_surface_swarm_particles_properties(magmatic_layer); CHKERRQ(ierr);


    }

}


PetscErrorCode magmatic_extrusive(PetscReal magma_volume,PetscReal x_center)
{
    PetscErrorCode ierr;
    PetscMPIInt rank;

    PetscFunctionBeginUser;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    PetscInt n;
    PetscInt j;
    PetscInt si;
    PetscInt bs;

    PetscInt nlocal;
    PetscInt seq_surface_size;

    Vec global_surface;
    Vec seq_surface;
    VecScatter ctx;
    PetscReal *seq_array = NULL;
    PetscReal *seq_array_copy = NULL;
    PetscReal *seq_array_aux = NULL;
    PetscReal *array = NULL;

    ierr = DMSwarmCreateGlobalVectorFromField(dms_s, DMSwarmPICField_coor, &global_surface); CHKERRQ(ierr);
    ierr = VecScatterCreateToZero(global_surface, &ctx, &seq_surface); CHKERRQ(ierr);
    ierr = VecScatterBegin(ctx, global_surface, seq_surface, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterEnd(ctx, global_surface, seq_surface, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterDestroy(&ctx); CHKERRQ(ierr);
    ierr = DMSwarmDestroyGlobalVectorFromField(dms_s, DMSwarmPICField_coor, &global_surface); CHKERRQ(ierr);

    ierr = VecGetSize(seq_surface, &seq_surface_size); CHKERRQ(ierr);
    ierr = VecGetArray(seq_surface, &seq_array); CHKERRQ(ierr);

    n = seq_surface_size/2;


    if (!rank) {

        PetscReal diff_h;
        PetscReal dx_sed = seq_array[2*1]-seq_array[2*0];
        
        PetscReal Delta_x;

        //magmatic

        PetscReal sigma = 10000.0; // !!! make it free for the user

        //maximum thickness of the magmatism: Area/(sigma*sqrt(2*pi))
        PetscReal magma_amplitude = magma_volume/(sigma*sqrt(2*3.14159));

        for (j=0; j<n; j++){
            Delta_x = j*dx_sed - x_center;
            diff_h = magma_amplitude*exp(-Delta_x*Delta_x/(2*sigma*sigma));

            seq_array[2*j+1] += diff_h;
        }

        
    }

    ierr = MPI_Bcast(&seq_surface_size, 1, MPI_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);

    if (rank) {
        ierr = PetscCalloc1(seq_surface_size, &seq_array_copy); CHKERRQ(ierr);
        ierr = MPI_Bcast(seq_array_copy, seq_surface_size, MPIU_SCALAR, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
        seq_array_aux = seq_array_copy;
    } else {
        ierr = MPI_Bcast(seq_array, seq_surface_size, MPIU_SCALAR, 0, PETSC_COMM_WORLD);
        seq_array_aux = seq_array;
    }

    ierr = DMDAGetCorners(da_Veloc, &si, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);
    ierr = DMSwarmGetLocalSize(dms_s, &nlocal); CHKERRQ(ierr);
    ierr = DMSwarmGetField(dms_s, DMSwarmPICField_coor, &bs, NULL, (void**)&array); CHKERRQ(ierr);

    for (j = 0; j < nlocal; j++) {
        array[2*j] = seq_array_aux[si*dms_s_ppe*2+2*j];
        array[2*j+1] = seq_array_aux[si*dms_s_ppe*2+2*j+1];
    }

    ierr = DMSwarmRestoreField(dms_s, DMSwarmPICField_coor, &bs, NULL, (void **)&array); CHKERRQ(ierr);
    ierr = VecRestoreArray(seq_surface, &seq_array); CHKERRQ(ierr);
    ierr = VecDestroy(&seq_surface); CHKERRQ(ierr);
    if (rank) {
        ierr = PetscFree(seq_array_copy); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}