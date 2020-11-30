#include <petscdmda.h>
#include <petscdmswarm.h>
#include <petscmath.h>

extern long Nx;
extern double Lx;
extern int n_interfaces;
extern double dx_const;
extern double dz_const;
extern double depth;
extern double tempo;
extern DM dms;
extern DM da_Veloc;

extern PetscScalar *interfaces;
extern PetscScalar *inter_rho;
extern PetscScalar *inter_geoq;
extern PetscScalar *inter_H;
extern PetscScalar *inter_A;
extern PetscScalar *inter_n;
extern PetscScalar *inter_Q;
extern PetscScalar *inter_V;

extern Vec sp_surface_global;
extern Vec sp_surface_global_n;
extern Vec sp_surface_coords_global;
extern Vec sp_top_surface_global;
extern Vec sp_bot_surface_global;
extern Vec sp_mean_surface_global;
extern Vec sp_top_surface_global_n;
extern Vec sp_bot_surface_global_n;
extern PetscReal sp_dt;
extern PetscReal sp_eval_time;
extern PetscReal sp_last_eval_time;
extern PetscReal sp_c;
extern PetscScalar *topo_var_time;
extern PetscScalar *topo_var_rate;
extern long sp_n_profiles;

extern Vec sp_surface_global_aux;

extern PetscScalar *global_surface_array_helper;
extern PetscScalar *global_surface_array_helper_aux;

extern PetscInt sp_mode;
extern PetscScalar sp_d_c;
extern PetscBool a2l;

extern PetscScalar K_fluvial;
extern PetscScalar sea_level;

PetscErrorCode sp_create_surface_vec();
PetscErrorCode sp_interpolate_surface_particles_to_vec();
PetscErrorCode evaluate_surface_processes();
PetscErrorCode update_particles_properties();
PetscErrorCode sp_write_surface_vec(PetscInt i);
PetscErrorCode sp_destroy();
PetscErrorCode sp_topo_var(PetscReal dt, PetscInt size);
PetscErrorCode sp_diffusion(PetscReal dt, PetscInt size);
PetscErrorCode sp_fluvial(PetscReal dt, PetscInt size);
PetscErrorCode sp_fluvial2(PetscReal dt, PetscInt size);
PetscErrorCode DMDAGetElementCorners(DM da,PetscInt *sx,PetscInt *sz,PetscInt *mx,PetscInt *mz);


PetscErrorCode sp_create_surface_vec()
{
    PetscErrorCode ierr;
    PetscMPIInt rank;
    PetscMPIInt size;
    PetscInt i;
    PetscInt si;
    PetscInt sk;
    PetscInt milocal;
    PetscInt mklocal;
    PetscInt local_size;
    PetscInt global_size;
    PetscInt low;
    PetscInt high;
    PetscReal x;
    PetscReal dx;
    PetscReal sp_depth;


    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    ierr = DMDAGetCorners(da_Veloc, &si, &sk, NULL, &milocal, &mklocal, NULL); CHKERRQ(ierr);

    local_size = 0;

    for (i = si; i < si+milocal; i++) {
        sp_depth = interfaces[i + Nx*(n_interfaces-1)];

        if (sp_depth > sk*dz_const-depth && sp_depth < (sk+mklocal)*dz_const-depth) {
            local_size = milocal;
            break;
        }
    }

    ierr = VecCreate(PETSC_COMM_WORLD, &sp_surface_global); CHKERRQ(ierr);
    ierr = VecSetSizes(sp_surface_global, local_size, PETSC_DECIDE); CHKERRQ(ierr);
    ierr = VecSetType(sp_surface_global, VECMPI); CHKERRQ(ierr);

    ierr = VecDuplicate(sp_surface_global, &sp_top_surface_global_n); CHKERRQ(ierr);
    ierr = VecDuplicate(sp_surface_global, &sp_bot_surface_global_n); CHKERRQ(ierr);
    ierr = VecDuplicate(sp_surface_global, &sp_surface_coords_global); CHKERRQ(ierr);
    ierr = VecDuplicate(sp_surface_global, &sp_top_surface_global); CHKERRQ(ierr);
    ierr = VecDuplicate(sp_surface_global, &sp_bot_surface_global); CHKERRQ(ierr);
    ierr = VecDuplicate(sp_surface_global, &sp_mean_surface_global); CHKERRQ(ierr);

    ierr = VecDuplicate(sp_surface_global, &sp_surface_global_aux); CHKERRQ(ierr);


    ierr = VecGetSize(sp_surface_global, &global_size); CHKERRQ(ierr);

    dx = Lx/(global_size-1);
    ierr = VecGetOwnershipRange(sp_surface_coords_global, &low, &high); CHKERRQ(ierr);

    for (i = low; i < high; i++) {
        x = i*dx;
        ierr = VecSetValue(sp_surface_coords_global, i, x, INSERT_VALUES); CHKERRQ(ierr);
    }

    ierr = VecAssemblyBegin(sp_surface_coords_global); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(sp_surface_coords_global); CHKERRQ(ierr);

    ierr = VecZeroEntries(sp_surface_global); CHKERRQ(ierr);
    ierr = VecZeroEntries(sp_top_surface_global); CHKERRQ(ierr);
    ierr = VecZeroEntries(sp_bot_surface_global); CHKERRQ(ierr);
    ierr = VecZeroEntries(sp_mean_surface_global); CHKERRQ(ierr);

    // Allocate helper arrays
    PetscCalloc1(Nx, &global_surface_array_helper);
    PetscCalloc1(Nx, &global_surface_array_helper_aux);

    ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD, "rank=%d local_size=%-3d global_size=%-3d\n", rank, local_size, global_size); CHKERRQ(ierr);
    ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

PetscErrorCode sp_interpolate_surface_particles_to_vec()
{
    PetscErrorCode ierr;
    PetscMPIInt rank;
    PetscScalar bot_threshold = -1.0e99;
    PetscReal *sp_top_surface_global_n_values;
    PetscReal *sp_bot_surface_global_n_values;
    PetscReal *sp_top_surface_global_values;
    PetscReal *sp_bot_surface_global_values;
    PetscReal *sp_mean_surface_global_values;
    PetscInt nlocal;
    PetscInt bs;
    PetscReal *pcoords;
    PetscInt *layer;
    PetscInt sex;
    PetscInt sey;
    PetscInt mx;
    PetscInt my;
    PetscInt p;
    PetscReal px;
    PetscReal py;
    PetscInt i;
    PetscInt sp_surface_local_size;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    ierr = VecZeroEntries(sp_top_surface_global_n); CHKERRQ(ierr);
    ierr = VecZeroEntries(sp_bot_surface_global_n); CHKERRQ(ierr);
    ierr = VecZeroEntries(sp_top_surface_global); CHKERRQ(ierr);
    ierr = VecSet(sp_bot_surface_global, bot_threshold); CHKERRQ(ierr);
    ierr = VecZeroEntries(sp_mean_surface_global); CHKERRQ(ierr);

    ierr = VecGetLocalSize(sp_top_surface_global, &sp_surface_local_size);
    ierr = VecGetArray(sp_top_surface_global_n, &sp_top_surface_global_n_values); CHKERRQ(ierr);
    ierr = VecGetArray(sp_bot_surface_global_n, &sp_bot_surface_global_n_values); CHKERRQ(ierr);
    ierr = VecGetArray(sp_top_surface_global, &sp_top_surface_global_values); CHKERRQ(ierr);
    ierr = VecGetArray(sp_bot_surface_global, &sp_bot_surface_global_values); CHKERRQ(ierr);
    ierr = VecGetArray(sp_mean_surface_global, &sp_mean_surface_global_values); CHKERRQ(ierr);

	ierr = DMSwarmGetLocalSize(dms, &nlocal);CHKERRQ(ierr);
    ierr = DMSwarmGetField(dms, DMSwarmPICField_coor, &bs,NULL, (void**)&pcoords); CHKERRQ(ierr);
    ierr = DMSwarmGetField(dms, "layer", &bs, NULL, (void**)&layer); CHKERRQ(ierr);

    ierr = DMDAGetElementCorners(da_Veloc, &sex, &sey, &mx, &my); CHKERRQ(ierr);

    for (p = 0; sp_surface_local_size && p < nlocal; p++) {
        px = pcoords[2*p];
		py = pcoords[2*p+1];

        i = (PetscInt)(px/dx_const);
        if (i < 0 || i > sex+mx) {ierr = PetscPrintf(PETSC_COMM_SELF, "[sp_interpolate_surface_particles_to_vec] rank=%d estranho i=%d (px=%.6e, p=%d)\n", rank, i, px, p); CHKERRQ(ierr); exit(1);}

        if (layer[p] == n_interfaces) {
            if (py < sp_top_surface_global_values[i-sex]) {
                sp_top_surface_global_values[i-sex] = py;
                sp_top_surface_global_n_values[i-sex] = 1;
            }
        } else {
            if (py > sp_bot_surface_global_values[i-sex]) {
                sp_bot_surface_global_values[i-sex] = py;
                sp_bot_surface_global_n_values[i-sex] = 1;
            }
        }
    }

    MPI_Barrier(PETSC_COMM_WORLD);

    for (i = 0; i < sp_surface_local_size; i++) {
        if (i > 0 && (int)sp_top_surface_global_values[i] == 0) {
            sp_top_surface_global_values[i] = sp_top_surface_global_values[i-1];
            sp_top_surface_global_n_values[i] = 1;
        }

        if (i > 0 && sp_bot_surface_global_values[i] < -1.0e90) {
            sp_bot_surface_global_values[i] = sp_bot_surface_global_values[i-1];
            sp_bot_surface_global_n_values[i] = 1;
        }

        if ((int)sp_top_surface_global_n_values[i] == 0) {
            ierr = PetscPrintf(PETSC_COMM_SELF, "ERROR [sp_interpolate_surface_particles_to_vec] rank=%d sp_top_surface_global_n_values i=%d\n", rank, i); CHKERRQ(ierr);
            ierr = PetscPrintf(PETSC_COMM_SELF, "ERROR [sp_interpolate_surface_particles_to_vec] rank=%d sp_top_surface_global_values=%e\n", rank, sp_top_surface_global_values[i]); CHKERRQ(ierr);
            exit(1);
        }

        if ((int)sp_bot_surface_global_n_values[i] == 0) {
            ierr = PetscPrintf(PETSC_COMM_SELF, "ERROR [sp_interpolate_surface_particles_to_vec] rank=%d sp_bot_surface_global_n_values i=%d\n", rank, i); CHKERRQ(ierr);
            ierr = PetscPrintf(PETSC_COMM_SELF, "ERROR [sp_interpolate_surface_particles_to_vec] rank=%d sp_bot_surface_global_values=%e\n", rank, sp_bot_surface_global_values[i]); CHKERRQ(ierr);
            exit(1);
        }

        sp_mean_surface_global_values[i] = (sp_top_surface_global_values[i] + sp_bot_surface_global_values[i]) / 2.0;
    }

    MPI_Barrier(PETSC_COMM_WORLD);

    ierr = VecRestoreArray(sp_top_surface_global_n, &sp_top_surface_global_n_values); CHKERRQ(ierr);
    ierr = VecRestoreArray(sp_bot_surface_global_n, &sp_bot_surface_global_n_values); CHKERRQ(ierr);
    ierr = VecRestoreArray(sp_top_surface_global, &sp_top_surface_global_values); CHKERRQ(ierr);
    ierr = VecRestoreArray(sp_bot_surface_global, &sp_bot_surface_global_values); CHKERRQ(ierr);
    ierr = VecRestoreArray(sp_mean_surface_global, &sp_mean_surface_global_values); CHKERRQ(ierr);

    ierr = DMSwarmRestoreField(dms, "layer", &bs, NULL, (void**)&layer); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(dms, DMSwarmPICField_coor, &bs, NULL, (void**)&pcoords); CHKERRQ(ierr);

    ierr = VecCopy(sp_mean_surface_global, sp_surface_global); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode evaluate_surface_processes()
{
    PetscErrorCode ierr;
    PetscMPIInt rank;
    PetscInt j;
    PetscInt low;
    PetscInt high;
    PetscInt size;
    PetscReal dt;
    PetscReal *y;
    PetscReal *y_aux;
    PetscReal dt_sp;
    PetscReal dt_aux;

    VecScatter ctx;
    Vec seq_surface;
    PetscReal *seq_y;


    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    ierr = VecScatterCreateToZero(sp_surface_global, &ctx, &seq_surface); CHKERRQ(ierr);
    ierr = VecScatterBegin(ctx, sp_surface_global, seq_surface, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterEnd(ctx, sp_surface_global, seq_surface, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterDestroy(&ctx); CHKERRQ(ierr);

    ierr = VecCopy(sp_surface_global, sp_surface_global_aux); CHKERRQ(ierr);

    ierr = VecGetArray(sp_surface_global, &y);
    ierr = VecGetOwnershipRange(sp_surface_global, &low, &high); CHKERRQ(ierr);
    ierr = VecGetSize(sp_surface_global, &size); CHKERRQ(ierr);

    ierr = VecGetArray(sp_surface_global_aux, &y_aux);

    if (!rank) {
        ierr = PetscPrintf(PETSC_COMM_SELF, "[rank %d] Evaluating surface processes\n", rank); CHKERRQ(ierr);

        ierr = VecGetArray(seq_surface, &seq_y);
        for (j = 0; j < size; j++) {
            global_surface_array_helper[j] = seq_y[j];
            global_surface_array_helper_aux[j] = seq_y[j];
        }
        ierr = VecRestoreArray(seq_surface, &seq_y); CHKERRQ(ierr);

        dt = tempo-sp_last_eval_time;

        if (1 == sp_mode) {
            sp_topo_var(dt, size);
        } else if (2 == sp_mode) {
            sp_diffusion(dt, size);
        } else if (3 == sp_mode) {
            sp_fluvial(dt, size);
        } else if (4 == sp_mode) {
            sp_fluvial2(dt, size);
        } else {
            ierr = PetscPrintf(PETSC_COMM_WORLD,"ERROR evaluate_surface_processes sp_mode undetermined\n"); CHKERRQ(ierr);
            exit(1);
        }

    } else {
        ierr = PetscPrintf(PETSC_COMM_SELF, "[rank %d] Wating for surface processes evaluation\n", rank); CHKERRQ(ierr);
    }

    ierr = VecDestroy(&seq_surface); CHKERRQ(ierr);

    MPI_Barrier(PETSC_COMM_WORLD);

    MPI_Bcast(global_surface_array_helper, size, MPIU_SCALAR, 0, PETSC_COMM_WORLD);
    MPI_Bcast(global_surface_array_helper_aux, size, MPIU_SCALAR, 0, PETSC_COMM_WORLD);

    MPI_Barrier(PETSC_COMM_WORLD);

    ierr = PetscPrintf(PETSC_COMM_SELF, "[rank %d] Adjusting local part of global vectors from global seq arrays\n", rank); CHKERRQ(ierr);
    for (j = low; j < high; j++) {
        y[j-low] = global_surface_array_helper[j];
        y_aux[j-low] = global_surface_array_helper_aux[j];
    }

    // PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
    // MPI_Barrier(PETSC_COMM_WORLD);

    ierr = VecRestoreArray(sp_surface_global, &y); CHKERRQ(ierr);
    ierr = VecRestoreArray(sp_surface_global_aux, &y_aux); CHKERRQ(ierr);

    MPI_Barrier(PETSC_COMM_WORLD);

    update_particles_properties();

    PetscFunctionReturn(0);
}

PetscErrorCode update_particles_properties()
{
    PetscErrorCode ierr;
    PetscMPIInt rank;
    PetscInt i;
    PetscInt p;
    PetscInt bs;
    PetscInt mx;
    PetscInt my;
    PetscInt sex;
    PetscInt sey;
    PetscInt nlocal;
    PetscInt sp_size;
    PetscInt sp_size_local;
	PetscInt *iarray;
	PetscInt *carray;
    PetscInt *layer;
    PetscReal px;
    PetscReal py;
    PetscReal ys;
    PetscReal rx;
    PetscReal sp_dx;
    PetscReal *y;
    PetscReal *pcoords;
	PetscReal *geoq_fac;
	PetscReal *strain_fac;

    PetscReal ys_aux;
    PetscReal *y_aux;
    PetscReal dh;
    PetscReal epsilon;


    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    ierr = PetscPrintf(PETSC_COMM_WORLD, "[update_particles_properties] (started) [rank %d]\n", rank); CHKERRQ(ierr);

	ierr = DMSwarmGetLocalSize(dms, &nlocal);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms, DMSwarmPICField_coor, &bs, NULL, (void**)&pcoords);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms, "cont", &bs, NULL, (void**)&carray);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms, "itag", &bs, NULL, (void**)&iarray);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms, "layer", &bs, NULL, (void**)&layer);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms, "geoq_fac", &bs, NULL, (void**)&geoq_fac);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms, "strain_fac", &bs, NULL, (void**)&strain_fac);CHKERRQ(ierr);

    ierr = VecGetArray(sp_surface_global, &y); CHKERRQ(ierr);
    ierr = VecGetSize(sp_surface_global, &sp_size); CHKERRQ(ierr);
    ierr = VecGetLocalSize(sp_surface_global, &sp_size_local); CHKERRQ(ierr);
    sp_dx = Lx/(sp_size-1);

    ierr = DMDAGetElementCorners(da_Veloc, &sex, &sey, &mx, &my); CHKERRQ(ierr);

    ierr = VecGetArray(sp_surface_global_aux, &y_aux); CHKERRQ(ierr);


    // ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD, "(debug) [update_particles_properties][rank %d] nlocal %d | sp_size %d | sp_size_local | %d | sex %d | mx %d\n", rank, nlocal, sp_size, sp_size_local, sex, mx); CHKERRQ(ierr);
    // ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT); CHKERRQ(ierr);

    for (p = 0; sp_size_local && p < nlocal; p++) {
        px = pcoords[2*p];
		py = pcoords[2*p+1];

        i = (PetscInt)(px/dx_const);
        if (i < 0 || i > sex+mx) {ierr = PetscPrintf(PETSC_COMM_SELF, "(error) [sp_interpolate_surface_particles_to_vec] rank=%d estranho i=%d (px=%.6e, p=%d)\n", rank, i, px, p); CHKERRQ(ierr); exit(1);}

        // rx = (px-i*sp_dx)/sp_dx;
        rx = 0.5;

        if (rx < 0 || rx > 1) {
            ierr = PetscPrintf(PETSC_COMM_SELF, "(error) [update_particles_properties] rx=%f px=%e i=%d sp_dx=%e dx_const=%e sp_size=%d\n", rx, px, i, sp_dx, dx_const, sp_size); CHKERRQ(ierr);
            exit(1);
        }

        if (i == sex+mx-1) {
            i--;
        }

        ys = y[i-sex] * (1.0 - rx) + y[i-sex+1] * rx;

        ys_aux = y_aux[i-sex] * (1.0 - rx) + y_aux[i-sex+1] * rx;

        dh = ys - ys_aux;

        epsilon = 1.00e-10;

        // (A2L) air particle bellow the surface => assign land properties
        if (a2l == PETSC_TRUE && layer[p] == n_interfaces) {

            if (py < ys_aux) {
                // set local surface to a level below the air particle
                ys_aux = py - epsilon;
            }

            if (py < (ys_aux + dh)) {
                ierr = PetscPrintf(PETSC_COMM_SELF, "(debug) [update_particles_properties] [rank %d] AIR  => LAND : p %d | layer %d=>%d | py %.8e < ys %.8e | x(%.8e < %.8e < %.8e)\n", rank, p, layer[p], n_interfaces-1, py, ys, i*sp_dx, px, (i+1)*sp_dx); CHKERRQ(ierr);
                // ierr = PetscPrintf(PETSC_COMM_SELF, "(debug) [update_particles_properties] [rank %d] A2L i %3d | ie %3d | rx %.2f | %e %e %e %e %d\n", rank, i, i-sex, rx, y[i-sex], y[i-sex] * (1.0 - rx) + y[i-sex+1] * rx, y[i-sex+1], py, layer[p]); CHKERRQ(ierr);

                layer[p] = n_interfaces - 1;
                geoq_fac[p] = inter_geoq[n_interfaces - 1];
                strain_fac[p] = 0.0;
                // carray[p] = 0;
                // iarray[p] = 0;
            }
        }

        // (L2A) land particle above the surface => assign air properties
        if (layer[p] != n_interfaces) {

            if (py > ys_aux) {
                // set local surface to a level above the land particle
                ys_aux = py + epsilon;
            }

            if (py > (ys_aux + dh)) {
                ierr = PetscPrintf(PETSC_COMM_SELF, "(debug) [update_particles_properties] [rank %d] LAND => AIR  : p %d | layer %d=>%d | py %.8e > ys %.8e | x(%.8e < %.8e < %.8e)\n", rank, p, layer[p], n_interfaces, py, ys, i*sp_dx, px, (i+1)*sp_dx); CHKERRQ(ierr);
                // ierr = PetscPrintf(PETSC_COMM_SELF, "(debug) [update_particles_properties] [rank %d] L2A i %3d | ie %3d | rx %.2f | %e %e %e %e %d\n", rank, i, i-sex, rx, y[i-sex], y[i-sex] * (1.0 - rx) + y[i-sex+1] * rx, y[i-sex+1], py, layer[p]); CHKERRQ(ierr);

                layer[p] = n_interfaces;
                geoq_fac[p] = inter_geoq[n_interfaces];
                strain_fac[p] = 0.0;
                // carray[p] = 0;
                // iarray[p] = 0;
            }
        }
    }

	ierr = DMSwarmRestoreField(dms, "cont", &bs, NULL, (void**)&carray);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms, "itag", &bs, NULL, (void**)&iarray);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms, "layer", &bs, NULL, (void**)&layer);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms, "geoq_fac", &bs, NULL, (void**)&geoq_fac);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms, "strain_fac", &bs, NULL, (void**)&strain_fac);CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(dms, DMSwarmPICField_coor, &bs, NULL, (void**)&pcoords); CHKERRQ(ierr);
    ierr = DMSwarmMigrate(dms, PETSC_TRUE); CHKERRQ(ierr);
    ierr = VecRestoreArray(sp_surface_global, &y); CHKERRQ(ierr);
    pcoords = NULL;

    ierr = PetscPrintf(PETSC_COMM_WORLD, "[update_particles_properties] (completed) [rank %d]\n", rank); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode sp_write_surface_vec(PetscInt i)
{
    PetscErrorCode ierr;
    PetscMPIInt rank;
    PetscViewer viewer;
    char filename[100];

	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    sprintf(filename, "sp_surface_global_%010d.txt", i); CHKERRQ(ierr);
    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer); CHKERRQ(ierr);
    ierr = VecView(sp_surface_global, viewer); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

    MPI_Barrier(PETSC_COMM_WORLD);

    PetscFunctionReturn(0);
}

PetscErrorCode sp_destroy()
{
    PetscErrorCode ierr;

    ierr = VecDestroy(&sp_surface_global); CHKERRQ(ierr);
    ierr = VecDestroy(&sp_surface_global_n); CHKERRQ(ierr);
    ierr = VecDestroy(&sp_surface_coords_global); CHKERRQ(ierr);
    ierr = VecDestroy(&sp_top_surface_global); CHKERRQ(ierr);
    ierr = VecDestroy(&sp_bot_surface_global); CHKERRQ(ierr);
    ierr = VecDestroy(&sp_mean_surface_global); CHKERRQ(ierr);
    ierr = VecDestroy(&sp_top_surface_global_n); CHKERRQ(ierr);
    ierr = VecDestroy(&sp_bot_surface_global_n); CHKERRQ(ierr);
    ierr = VecDestroy(&sp_surface_global_aux); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode sp_topo_var(PetscReal dt, PetscInt size)
{
    PetscInt t;
    PetscInt j;


    for (t = 0; t < sp_n_profiles-1; t++) {
        if ((sp_last_eval_time > topo_var_time[t] || fabs(sp_last_eval_time-topo_var_time[t]) < 1.0e-13) && (sp_last_eval_time < topo_var_time[t+1] || fabs(sp_last_eval_time - topo_var_time[t+1]) < 1.0e-13)) {
            break;
        }
    }

    for (j = 0; j < size; j++) {
        global_surface_array_helper[j] = global_surface_array_helper[j] - topo_var_rate[j + t*Nx]*dt;
    }

    PetscFunctionReturn(0);
}

PetscErrorCode sp_diffusion(PetscReal dt, PetscInt size)
{
    PetscErrorCode ierr;
    PetscInt t;
    PetscInt j;
    PetscInt max_steps;
    PetscReal r;
    PetscReal sp_dx;
    PetscReal sp_dt;
    PetscScalar *sp_y_aux;


    ierr = PetscCalloc1(size, &sp_y_aux); CHKERRQ(ierr);
    ierr = PetscArraycpy(sp_y_aux, global_surface_array_helper, size); CHKERRQ(ierr);

    sp_dx = Lx/(size-1);

    sp_dt = sp_dx*sp_dx/sp_d_c/4.0;
    max_steps = (int)(PetscFloorReal(dt/sp_dt)) + 1;
    if (max_steps < 10) {
        max_steps = 10;
    }
    sp_dt = dt/max_steps;

    r = sp_d_c*sp_dt/(sp_dx*sp_dx);

    ierr = PetscPrintf(PETSC_COMM_SELF, "[sp_diffusion] dt=%e\n", dt); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "[sp_diffusion] size=%d\n", size); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "[sp_diffusion] sp_dx=%e\n", sp_dx); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "[sp_diffusion] sp_dt=%e\n", sp_dt); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "[sp_diffusion] max_steps=%d\n", max_steps); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "[sp_diffusion] sp_d_c=%e\n", sp_d_c); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "[sp_diffusion] r=%e\n", r); CHKERRQ(ierr);
    // for (j = 0; j < size; j++) {
    //     ierr = PetscPrintf(PETSC_COMM_SELF, "[sp_diffusion] sp_y_aux=%e\n", sp_y_aux[j]); CHKERRQ(ierr);
    // }


    global_surface_array_helper[0] = global_surface_array_helper[1];
    global_surface_array_helper[size-1] = global_surface_array_helper[size-2];

    for (t = 0; t < max_steps; t++) {
        for (j = 1; j < size-1; j++) {
            global_surface_array_helper[j] = sp_y_aux[j] + r*(sp_y_aux[j+1] - 2.0*sp_y_aux[j] + sp_y_aux[j-1]);
        }

        for (j = 0; j < size; j++) {
            sp_y_aux[j] = global_surface_array_helper[j];
        }
    }

    PetscFunctionReturn(0);
}



PetscErrorCode sp_fluvial(PetscReal dt, PetscInt size)
{
    PetscErrorCode ierr;
    PetscInt t;
    PetscInt j,cont;
    PetscInt max_steps;

    PetscReal sp_dx;
    PetscReal sp_dt;
    PetscScalar *h,*h_aux,*q,*fc;
    PetscInt *hi;
    PetscInt sum;
    
    PetscReal K = K_fluvial;//2.0E-7;


    ierr = PetscCalloc1(size, &h); CHKERRQ(ierr);
    ierr = PetscArraycpy(h, global_surface_array_helper, size); CHKERRQ(ierr);

    sp_dx = Lx/(size-1);

    //sp_dt = sp_dx*sp_dx/sp_d_c/4.0;
    //max_steps = (int)(PetscFloorReal(dt/sp_dt)) + 1;
    //if (max_steps < 10) {
    //    max_steps = 10;
    //}
    sp_dt = 5000.0;
    max_steps=(int)(dt/sp_dt);
    sp_dt = dt/max_steps;

    //r = sp_d_c*sp_dt/(sp_dx*sp_dx);

    ierr = PetscPrintf(PETSC_COMM_SELF, "[sp_diffusion] dt=%e\n", dt); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "[sp_diffusion] size=%d\n", size); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "[sp_diffusion] sp_dx=%e\n", sp_dx); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "[sp_diffusion] sp_dt=%e\n", sp_dt); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "[sp_diffusion] max_steps=%d\n", max_steps); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "[sp_diffusion] sp_d_c=%e\n", sp_d_c); CHKERRQ(ierr);
    //ierr = PetscPrintf(PETSC_COMM_SELF, "[sp_diffusion] r=%e\n", r); CHKERRQ(ierr);
    // for (j = 0; j < size; j++) {
    //     ierr = PetscPrintf(PETSC_COMM_SELF, "[sp_diffusion] sp_y_aux=%e\n", sp_y_aux[j]); CHKERRQ(ierr);
    // }

    PetscReal mean_h=0.0;
    for (j = 0,cont=0; j < size/4; j++) {
        mean_h+=h[j];
        cont++;
    }
    mean_h/=cont;

    PetscReal hsl = mean_h+sea_level;


    ierr = PetscCalloc1(size, &h_aux); CHKERRQ(ierr);
    ierr = PetscCalloc1(size, &q); CHKERRQ(ierr);
    ierr = PetscCalloc1(size, &fc); CHKERRQ(ierr);
    ierr = PetscCalloc1(size, &hi); CHKERRQ(ierr);

    for (j=0;j<size;j++) h_aux[j]=h[j];

    for (t=0; t < max_steps; t++) {
        for (j=0;j<size;j++) q[j] = 2000.0;    
        for (j=0;j<size;j++) {
            fc[j]=1;
            if (h[j]<hsl) fc[j]=0;
        }
        fc[0]=0;
        fc[size-1]=0;
        for (j=0;j<size;j++) {if (h[j]<hsl) h[j]=hsl;}
        for (j=0;j<size;j++) hi[j]=j;
        
        for (sum=0,j=0;j<size;j++) sum+=fc[j];
        while (sum>0){
            for (j=0;j<size;j++) {if (fc[j]==1) fc[j]=2;}
            for (j=1;j<size-1;j++){
                if (fc[j]>0){
                    if (h[j]>h[j-1]){
                        if (h[j+1]<h[j-1]){
                            hi[j]=j+1;
                            if (fc[j+1]>1) fc[j+1]=1;
                        }
                        else{
                            hi[j]=j-1;
                            if (fc[j-1]>1) fc[j-1]=1;
                        }
                    }
                    else if (h[j]>h[j+1]){
                        hi[j]=j+1;
                        if (fc[j+1]>1) fc[j+1]=1;
                    }
                }
            }
            for (j=1;j<size-1;j++){
                if (fc[j]==2){
                    h_aux[j]=h[j]-q[j]*K*sp_dt*(h[j]-h[hi[j]])/sp_dx;
                    q[hi[j]]+=q[j];
                    fc[j]=0;
                    if (h_aux[j]<hsl) h_aux[j]=hsl;
                }
            }
            for (sum=0.0,j=0;j<size;j++) sum+=fc[j];
            printf("sum = %d\n",sum);
        }
        for (j=0;j<size;j++) h[j]=h_aux[j];

    }

    for (j = 0; j < size; j++) {
        global_surface_array_helper[j] = h[j];
    }

    PetscFree(h);
    PetscFree(h_aux);
    PetscFree(fc);
    PetscFree(q);
    PetscFree(hi);
    

   

    PetscFunctionReturn(0);
}




PetscErrorCode sp_fluvial2(PetscReal dt, PetscInt size)
{
    PetscErrorCode ierr;
    PetscInt t;
    PetscInt j,jj,cont,imax;
    PetscInt max_steps;

    PetscReal sp_dx;
    PetscReal sp_dt;
    PetscScalar *h,*h_aux,*q;
    PetscInt *hi,*fc;
    PetscInt ci,cf,out;
    PetscInt sum;
    PetscReal hmax;
    
    PetscReal K = K_fluvial;//2.0E-7;


    ierr = PetscCalloc1(size, &h); CHKERRQ(ierr);
    ierr = PetscArraycpy(h, global_surface_array_helper, size); CHKERRQ(ierr);

    sp_dx = Lx/(size-1);

    //sp_dt = sp_dx*sp_dx/sp_d_c/4.0;
    //max_steps = (int)(PetscFloorReal(dt/sp_dt)) + 1;
    //if (max_steps < 10) {
    //    max_steps = 10;
    //}
    sp_dt = 5000.0;
    max_steps=(int)(dt/sp_dt);
    sp_dt = dt/max_steps;

    //r = sp_d_c*sp_dt/(sp_dx*sp_dx);

    ierr = PetscPrintf(PETSC_COMM_SELF, "[sp_diffusion] dt=%e\n", dt); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "[sp_diffusion] size=%d\n", size); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "[sp_diffusion] sp_dx=%e\n", sp_dx); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "[sp_diffusion] sp_dt=%e\n", sp_dt); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "[sp_diffusion] max_steps=%d\n", max_steps); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "[sp_diffusion] sp_d_c=%e\n", sp_d_c); CHKERRQ(ierr);
    //ierr = PetscPrintf(PETSC_COMM_SELF, "[sp_diffusion] r=%e\n", r); CHKERRQ(ierr);
    // for (j = 0; j < size; j++) {
    //     ierr = PetscPrintf(PETSC_COMM_SELF, "[sp_diffusion] sp_y_aux=%e\n", sp_y_aux[j]); CHKERRQ(ierr);
    // }

    PetscReal mean_h=0.0;
    for (j = 0,cont=0; j < size/4; j++) {
        mean_h+=h[j];
        cont++;
    }
    mean_h/=cont;

    PetscReal hsl = mean_h+sea_level;


    ierr = PetscCalloc1(size, &h_aux); CHKERRQ(ierr);
    ierr = PetscCalloc1(size, &q); CHKERRQ(ierr);
    ierr = PetscCalloc1(size, &fc); CHKERRQ(ierr);
    ierr = PetscCalloc1(size, &hi); CHKERRQ(ierr);

    for (j=0;j<size;j++) h_aux[j]=h[j];

    for (t=0; t < max_steps; t++) {
        for (j=0;j<size;j++) q[j] = sp_dx;    
        for (j=0;j<size;j++) {
            fc[j]=1;
            if (h[j]<hsl) fc[j]=0;
        }
        for (j=0;j<size;j++){
            if (j*sp_dx<40.0E3 || j*sp_dx>Lx-40.0E3) fc[j]=0;
        }
        fc[0]=0;
        fc[size-1]=0;
        for (j=0;j<size;j++) {if (h[j]<hsl) h[j]=hsl;}
        for (j=0;j<size;j++) hi[j]=j;


        ci=1;
        cf=1;
        out = 0;
        for (j=1;j<size;j++){
            if (out==0 && fc[j]==1){ 
                ci=j;
                cf=j;
                out=1;
            }
            if (out==1 && fc[j]==1){
                cf=j;
            }
            if (out==1 && fc[j]==0){
                imax=ci;
                hmax=h[imax];
                for (jj=ci;jj<cf+1;jj++){//ii in range(ci,cf+1):
                    if (hmax<h[jj]){
                        imax=jj;
                        hmax=h[imax];
                    }
                }
                for (jj=ci;jj<imax;jj++){//ii in range(ci,imax):
                    hi[jj]=jj-1;
                }
                for (jj=imax;jj<cf+1;jj++){//ii in range(imax,cf+1):
                    hi[jj]=jj+1;
                }
                out=0;
            }
        }


        for (sum=0,j=0;j<size;j++) sum+=fc[j];
        printf("sum = %d\n",sum);
        while (sum>0){
            for (j=0;j<size;j++) {if (fc[j]==1) fc[j]=2;}
            for (j=1;j<size-1;j++){
                if (fc[j]>0){
                    if (fc[hi[j]]>0){
                        fc[hi[j]]=1;
                    }
                }
            }
            for (j=1;j<size-1;j++){
                if (fc[j]==2){
                    if (h[j]>h[hi[j]]) h_aux[j]=h[j]-q[j]*K*sp_dt*(h[j]-h[hi[j]])/sp_dx;
                    q[hi[j]]+=q[j];
                    fc[j]=0;
                    if (h_aux[j]<hsl) h_aux[j]=hsl;
                }
            }
            for (sum=0,j=0;j<size;j++) sum+=fc[j];
            printf("sum = %d\n",sum);
        }
        for (j=0;j<size;j++) h[j]=h_aux[j];

    }

    for (j = 0; j < size; j++) {
        global_surface_array_helper[j] = h[j];
    }

    PetscFree(h);
    PetscFree(h_aux);
    PetscFree(fc);
    PetscFree(q);
    PetscFree(hi);
    



    /*
    global_surface_array_helper[0] = global_surface_array_helper[1];
    global_surface_array_helper[size-1] = global_surface_array_helper[size-2];

    for (t = 0; t < max_steps; t++) {
        for (j = 1; j < size-1; j++) {
            global_surface_array_helper[j] = sp_y_aux[j] + r*(sp_y_aux[j+1] - 2.0*sp_y_aux[j] + sp_y_aux[j-1]);
        }

        for (j = 0; j < size; j++) {
            sp_y_aux[j] = global_surface_array_helper[j];
        }
    }
    */

   

    PetscFunctionReturn(0);
}

