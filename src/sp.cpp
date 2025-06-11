#include <petscdmda.h>
#include <petscdmshell.h>
#include <petscdmswarm.h>
#include <petsc/private/dmimpl.h>
#include <petscmath.h>
#include <petscviewer.h>

#include "sp_mode.h"

extern long Nx;
extern double Lx;
extern int n_interfaces;
extern double dx_const;
extern double dz_const;
extern double depth;
extern DM dms;
extern DM da_Veloc;

extern double seg_per_ano;

extern PetscScalar *interfaces;
extern PetscScalar *inter_rho;
extern PetscScalar *inter_geoq;

extern SP_Mode sp_mode;
extern PetscScalar sp_d_c;
extern PetscBool a2l;
extern PetscScalar sea_level;
extern PetscReal Ks;
extern PetscReal lambda_s;

extern DM dmcell_s;
extern DM dms_s;
extern PetscInt buffer_s;
extern PetscInt dms_s_ppe;
extern PetscReal epsilon_x;
extern Vec local_V;
extern Vec Veloc_weight;

typedef struct {
	PetscScalar u;
	PetscScalar w;
} Stokes2d;

PetscErrorCode sp_create_surface_swarm_2d();
PetscErrorCode sp_move_surface_swarm(PetscInt dimensions, PetscReal dt);
PetscErrorCode sp_move_surface_swarm_advection_2d(PetscReal dt);
PetscReal sp_evaluate_adjusted_mean_elevation_with_sea_level();
PetscErrorCode sp_surface_swarm_interpolation();
PetscErrorCode sp_evaluate_surface_processes(PetscInt dimensions, PetscReal dt);
PetscErrorCode sp_evaluate_surface_processes_2d(PetscReal dt);
PetscErrorCode sp_evaluate_surface_processes_2d_diffusion(PetscReal dt);
PetscErrorCode sp_evaluate_surface_processes_2d_sedimentation_only(PetscReal dt);
PetscErrorCode sp_evaluate_surface_processes_2d_diffusion_sedimentation_only(PetscReal dt);
PetscErrorCode sp_update_surface_swarm_particles_properties();
PetscErrorCode DMLocatePoints_DMDARegular_2d(DM dm,Vec pos,DMPointLocationType ltype, PetscSF cellSF);
PetscErrorCode DMGetNeighbors_DMDARegular_2d(DM dm,PetscInt *nneighbors,const PetscMPIInt **neighbors);
PetscErrorCode sp_view_2d(DM dm, const char prefix[]);
PetscErrorCode sp_destroy();
PetscErrorCode DMDAGetElementCorners_2d(DM da,PetscInt *sx,PetscInt *sz,PetscInt *mx,PetscInt *mz);
PetscReal linear_interpolation(PetscReal rx, PetscReal rz,PetscScalar V0, PetscScalar V1, PetscScalar V2, PetscScalar V3);
PetscInt get_i(PetscReal cx);
PetscInt get_k(PetscReal cz);
PetscReal get_rx(PetscReal cx, PetscInt i);
PetscReal get_rz(PetscReal cz, PetscInt k);


PetscErrorCode sp_create_surface_swarm_2d()
{
    PetscErrorCode ierr;
    PetscMPIInt rank;

    PetscFunctionBeginUser;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    ierr = DMShellCreate(PETSC_COMM_WORLD, &dmcell_s); CHKERRQ(ierr);
    ierr = DMSetApplicationContext(dmcell_s, (void*)da_Veloc); CHKERRQ(ierr);
    dmcell_s->ops->locatepoints = DMLocatePoints_DMDARegular_2d;
    dmcell_s->ops->getneighbors = DMGetNeighbors_DMDARegular_2d;

    ierr = DMCreate(PETSC_COMM_WORLD, &dms_s); CHKERRQ(ierr);
    ierr = DMSetType(dms_s, DMSWARM); CHKERRQ(ierr);
    ierr = DMSetDimension(dms_s, 2); CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)dms_s, "surface"); CHKERRQ(ierr);

    ierr = DMSwarmSetType(dms_s, DMSWARM_PIC); CHKERRQ(ierr);
    ierr = DMSwarmSetCellDM(dms_s, dmcell_s); CHKERRQ(ierr);

    /* fields */
    ierr = DMSwarmRegisterPetscDatatypeField(dms_s, "itag", 1, PETSC_INT); CHKERRQ(ierr);
    ierr = DMSwarmFinalizeFieldRegister(dms_s); CHKERRQ(ierr);

    {
        PetscInt p;
        PetscInt nlocal;
        PetscInt si;
        PetscInt sj;
        PetscInt milocal;
        PetscInt mjlocal;
        PetscInt bs;
        PetscInt cnt;
        PetscInt *iarray;
        PetscReal *array;
        PetscReal px;
        PetscReal pz;
        PetscReal rx;

        PetscReal z_top;
        PetscReal z_bottom;

        ierr = DMDAGetCorners(da_Veloc, &si, &sj, NULL, &milocal, &mjlocal, NULL); CHKERRQ(ierr);

        buffer_s = milocal * dms_s_ppe / 2;
        ierr = DMSwarmSetLocalSizes(dms_s, milocal * dms_s_ppe, buffer_s); CHKERRQ(ierr);
        ierr = DMSwarmGetLocalSize(dms_s, &nlocal); CHKERRQ(ierr);

        ierr = DMSwarmGetField(dms_s, DMSwarmPICField_coor, &bs, NULL, (void **)&array); CHKERRQ(ierr);

        z_top = -depth + (sj + mjlocal - 1) * dz_const;
        z_bottom = -depth + sj * dz_const;

        cnt = 0;

        for (p = si; p < si+milocal; p++) {
            for (PetscInt i = 0; i < dms_s_ppe; i++) {
                rx = (i * dx_const)/(dms_s_ppe);
                px = p * dx_const + rx;

                if (p < milocal - 1) {
                    pz = (1.0 - (rx / dx_const)) * (interfaces[(n_interfaces-1) * Nx + p]) + (rx / dx_const) * (interfaces[(n_interfaces-1) * Nx + p + 1]);
                } else {
                    pz = interfaces[(n_interfaces-1) * Nx + p];
                }

                if ((px >= 0) && (px <= Lx) && (pz >= z_bottom) && (pz <= z_top)) {
                    array[bs * cnt + 0] = px;
                    array[bs * cnt + 1] = pz;

                    if (px == 0) {
                        array[bs * cnt + 0] += epsilon_x;
                    }

                    if (px == Lx) {
                        array[bs * cnt + 0] -= epsilon_x;
                    }

                    cnt++;
                }
            }
        }

        ierr = DMSwarmRestoreField(dms_s, DMSwarmPICField_coor, &bs, NULL, (void **)&array); CHKERRQ(ierr);
        ierr = DMSwarmSetLocalSizes(dms_s, cnt, buffer_s); CHKERRQ(ierr);

        ierr = DMSwarmGetLocalSize(dms_s, &nlocal); CHKERRQ(ierr);

        ierr = DMSwarmGetField(dms_s, "itag", &bs, NULL, (void **)&iarray); CHKERRQ(ierr);
        ierr = DMSwarmGetField(dms_s, DMSwarmPICField_coor, &bs, NULL, (void **)&array); CHKERRQ(ierr); //
        for (p = 0; p < nlocal; p++) {
            iarray[p] = si * dms_s_ppe + p;
        }
        ierr = DMSwarmRestoreField(dms_s, "itag", &bs, NULL, (void **)&iarray); CHKERRQ(ierr);
        ierr = DMSwarmRestoreField(dms_s, DMSwarmPICField_coor, &bs, NULL, (void **)&array); CHKERRQ(ierr); //

        ierr = DMView(dms_s, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
        ierr = sp_view_2d(dms_s, "inital_surface"); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}

PetscErrorCode sp_move_surface_swarm(PetscInt dimensions, PetscReal dt)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    if (dimensions == 2) {
        ierr = sp_move_surface_swarm_advection_2d(dt); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}

PetscErrorCode sp_move_surface_swarm_advection_2d(PetscReal dt)
{
    PetscErrorCode ierr;
    PetscMPIInt rank;

    PetscFunctionBeginUser;

	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    // 2D
	// Velocity
	Stokes2d **VV;

    ierr = VecZeroEntries(local_V); CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(da_Veloc, Veloc_weight, INSERT_VALUES, local_V); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(da_Veloc, Veloc_weight, INSERT_VALUES, local_V); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_Veloc, local_V, &VV); CHKERRQ(ierr);

    PetscInt p;
    PetscInt i;
    PetscInt k;
    PetscReal cx;
    PetscReal cz;
    PetscReal vx;
    PetscReal vz;
    PetscReal rx;
    PetscReal rz;
    PetscInt nlocal;
    PetscInt bs;
    PetscReal *array;

    ierr = DMSwarmGetLocalSize(dms_s, &nlocal); CHKERRQ(ierr);
    ierr = DMSwarmGetField(dms_s, DMSwarmPICField_coor, &bs, NULL, (void**)&array); CHKERRQ(ierr);

    for (p = 0; p < nlocal; p++) {

        cx = array[2 * p];
        cz = array[2 * p + 1];

        if (cx >= Lx) {
            printf("[SP] surface horizontal advetion - outside: cx=%.4e>=%.4e\n\n", cx, Lx);
            cx = Lx - epsilon_x;
        }

        if (cx <= 0.0) {
            printf("[SP] surface horizontal advetion - outside: cx=%.4e<=0.0\n\n", cx);
            cx = epsilon_x;
        }

        if (cz >= 0.0){
            printf("[SP] surface horizontal advetion - outside: cz=%.4e>=0.0\n\n", cz);
            cz = -epsilon_x;
        }

        if (cz <= -depth){
            printf("[SP] surface horizontal advetion - outside: cz=%.4e<=-%.4e\n\n", cz, depth);
            cz = -depth + epsilon_x;
        }

        i = get_i(cx);
        k = get_k(cz);
        rx = get_rx(cx, i);
        rz = get_rz(cz, k);
        vx = linear_interpolation(rx, rz, VV[k][i].u, VV[k][i+1].u, VV[k+1][i].u, VV[k+1][i+1].u);
        vz = linear_interpolation(rx, rz, VV[k][i].w, VV[k][i+1].w, VV[k+1][i].w, VV[k+1][i+1].w);
        array[2 * p] += dt * vx;
        array[2 * p + 1] += dt * vz;
    }

    ierr = DMSwarmRestoreField(dms_s, DMSwarmPICField_coor, &bs, NULL, (void**)&array); CHKERRQ(ierr);

    ierr = DMSwarmMigrate(dms_s, PETSC_TRUE); CHKERRQ(ierr);

    ierr = DMDAVecRestoreArray(da_Veloc, local_V, &VV); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscReal sp_evaluate_adjusted_mean_elevation_with_sea_level()
{
    PetscErrorCode ierr;
    PetscMPIInt rank;

    PetscReal mean_h = 0.0;
    PetscInt j;
    PetscInt cont;
    PetscReal hsl; // mean elevation plus sea level height

    Vec global_surface;
    Vec seq_surface;
    VecScatter ctx;
    PetscInt seq_surface_size;
    PetscReal *seq_array;


    PetscFunctionBeginUser;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    ierr = DMSwarmCreateGlobalVectorFromField(dms_s, DMSwarmPICField_coor, &global_surface); CHKERRQ(ierr);
    ierr = VecScatterCreateToAll(global_surface, &ctx, &seq_surface); CHKERRQ(ierr);
    ierr = VecScatterBegin(ctx, global_surface, seq_surface, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterEnd(ctx, global_surface, seq_surface, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterDestroy(&ctx); CHKERRQ(ierr);
    ierr = DMSwarmDestroyGlobalVectorFromField(dms_s, DMSwarmPICField_coor, &global_surface); CHKERRQ(ierr);

    ierr = VecGetSize(seq_surface, &seq_surface_size); CHKERRQ(ierr);
    ierr = VecGetArray(seq_surface, &seq_array); CHKERRQ(ierr);

    mean_h = 0.0;

    for (j = 0, cont = 0; j < seq_surface_size / 2 / 4; j++) {
        mean_h += seq_array[2 * j + 1];
        cont++;
    }

    mean_h /= cont;

    hsl = mean_h + sea_level;

    ierr = VecRestoreArray(seq_surface, &seq_array); CHKERRQ(ierr);
    ierr = VecDestroy(&seq_surface); CHKERRQ(ierr);

    return hsl;
}

PetscErrorCode sp_surface_swarm_interpolation()
{
    PetscMPIInt rank;
	PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    PetscInt si;
    PetscInt sj;
    PetscInt milocal;
    PetscInt mjlocal;
    PetscReal rx;
    PetscReal px;

    PetscInt nlocal;
    PetscInt bs;
    PetscInt *iarray;
    PetscReal *array;

    PetscReal *seq_array;

    PetscReal cx_aux;
    PetscReal cz_aux;
    PetscReal cz;
    PetscReal cx_l;
    PetscReal cz_l;
    PetscReal cx_r;
    PetscReal cz_r;

    PetscReal delta_x_l;
    PetscReal delta_x_r;
    PetscBool p_aux_l;
    PetscBool p_aux_r;

    PetscInt i;
    PetscInt j;
    PetscInt p;

    Vec global_surface;
    Vec seq_surface;
    VecScatter ctx;
    PetscInt seq_surface_size;


    ierr = DMDAGetCorners(da_Veloc, &si, &sj, NULL, &milocal, &mjlocal, NULL); CHKERRQ(ierr);

    ierr = DMSwarmGetLocalSize(dms_s, &nlocal); CHKERRQ(ierr);

    ierr = DMSwarmCreateGlobalVectorFromField(dms_s, DMSwarmPICField_coor, &global_surface); CHKERRQ(ierr);
    ierr = VecScatterCreateToAll(global_surface, &ctx, &seq_surface); CHKERRQ(ierr);
    ierr = VecScatterBegin(ctx, global_surface, seq_surface, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterEnd(ctx, global_surface, seq_surface, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterDestroy(&ctx); CHKERRQ(ierr);
    ierr = DMSwarmDestroyGlobalVectorFromField(dms_s, DMSwarmPICField_coor, &global_surface); CHKERRQ(ierr);

    ierr = VecGetSize(seq_surface, &seq_surface_size); CHKERRQ(ierr);
    ierr = VecGetArray(seq_surface, &seq_array); CHKERRQ(ierr);

    PetscInt dms_s_size;
    ierr = DMSwarmGetSize(dms_s, &dms_s_size); CHKERRQ(ierr);

    PetscReal *interpolated_surface;
    PetscInt interpolated_size = 2 * dms_s_ppe * (Nx - 1) + 1;

    if (nlocal) {
        ierr = PetscCalloc1(interpolated_size, &interpolated_surface); CHKERRQ(ierr);
    }

    PetscInt cnt = 0;

    for (i = si; nlocal && i < si+milocal; i++) {
        for (j = 0; j < dms_s_ppe; j++) {
            rx = (j * dx_const)/(dms_s_ppe);
            px = i * dx_const + rx;

            if (px > Lx) {
                break;
            }

            delta_x_l = Lx + epsilon_x;
            delta_x_r = Lx + epsilon_x;
            cx_l = seq_array[0];
            cz_l = seq_array[1];
            cx_r = seq_array[2*(dms_s_size-1)];
            cz_r = seq_array[2*(dms_s_size-1)+1];

            p_aux_l = PETSC_FALSE;
            p_aux_r = PETSC_FALSE;

            for (p = 0; p < dms_s_size; p++) {
                cx_aux = seq_array[2*p];
                cz_aux = seq_array[2*p+1];

                if (cx_aux < px && (px - cx_aux) < delta_x_l) {
                    cx_l = cx_aux;
                    cz_l = cz_aux;
                    delta_x_l = px - cx_aux;
                    p_aux_l = PETSC_TRUE;
                }

                if (cx_aux > px && (cx_aux - px) < delta_x_r) {
                    cx_r = cx_aux;
                    cz_r = cz_aux;
                    delta_x_r = cx_aux - px;
                    p_aux_r = PETSC_TRUE;
                }
            }

            if (p_aux_l == PETSC_TRUE && p_aux_r == PETSC_FALSE) {
                cx_r = px;
                cz_r = cz_l;
            } else if (p_aux_l == PETSC_FALSE && p_aux_r == PETSC_TRUE) {
                cx_l = px;
                cz_l = cz_r;
            }

            cz = (1.0 - (px - cx_l)/(cx_r - cx_l)) * cz_l + (px - cx_l)/(cx_r - cx_l) * cz_r;

            interpolated_surface[2*cnt] = px;
            interpolated_surface[2*cnt+1] = cz;

            cnt++;
        }
    }

    ierr = VecRestoreArray(seq_surface, &seq_array); CHKERRQ(ierr);
    ierr = VecDestroy(&seq_surface); CHKERRQ(ierr);

    for (i = 0; i < nlocal; i++) {
        ierr = DMSwarmRemovePoint(dms_s); CHKERRQ(ierr);
    }
    ierr = DMSwarmAddNPoints(dms_s, cnt);

    ierr = DMSwarmGetField(dms_s, DMSwarmPICField_coor, &bs, NULL, (void**)&array); CHKERRQ(ierr);
    ierr = DMSwarmGetField(dms_s, "itag", &bs, NULL, (void **)&iarray); CHKERRQ(ierr);
    for (p = 0; p < cnt; p++) {
        iarray[p] = si * dms_s_ppe + p;
        array[2*p] = interpolated_surface[2*p];
        array[2*p+1] = interpolated_surface[2*p+1];
    }
    ierr = DMSwarmRestoreField(dms_s, DMSwarmPICField_coor, &bs, NULL, (void **)&array); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(dms_s, "itag", &bs, NULL, (void **)&iarray); CHKERRQ(ierr);

    if (nlocal) {
        ierr = PetscFree(interpolated_surface); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}

PetscErrorCode sp_evaluate_surface_processes(PetscInt dimensions, PetscReal dt)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    if (dimensions == 2) {
        ierr = sp_evaluate_surface_processes_2d(dt); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}

PetscErrorCode sp_evaluate_surface_processes_2d(PetscReal dt)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    if (sp_mode == SP_DIFFUSION) {
        ierr = sp_evaluate_surface_processes_2d_diffusion(dt); CHKERRQ(ierr);
    }
    else if (sp_mode == SP_SEDIMENTATION_ONLY) {
        ierr = sp_evaluate_surface_processes_2d_sedimentation_only(dt); CHKERRQ(ierr);
    }
    else if (sp_mode == SP_DIFFUSION_SEDIMENTATION_ONLY) {
        ierr = sp_evaluate_surface_processes_2d_diffusion_sedimentation_only(dt); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}

PetscErrorCode sp_evaluate_surface_processes_2d_diffusion(PetscReal dt)
{
    PetscErrorCode ierr;
    PetscMPIInt rank;

    Vec global_surface;
    Vec seq_surface;
    VecScatter ctx;
    PetscInt seq_surface_size;
    PetscInt p;
    PetscInt bs;
    PetscReal *seq_array;
    PetscReal sp_dt;
    PetscInt t;
    PetscInt j;
    PetscInt max_steps;
    PetscReal r;
    PetscReal sp_dx;
    PetscScalar *z;
    PetscInt si;
    PetscInt nlocal;
    PetscReal *array;
    PetscInt n;

    PetscFunctionBeginUser;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

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
        ierr = PetscCalloc1(n, &z); CHKERRQ(ierr);

        sp_dx = Lx/(n-1);

        sp_dt = sp_dx*sp_dx/sp_d_c/4.0;
        max_steps = (int)(PetscFloorReal(dt/sp_dt)) + 1;
        if (max_steps < 10) {
            max_steps = 10;
        }
        sp_dt = dt/max_steps;

        r = sp_d_c*sp_dt/(sp_dx*sp_dx);

        seq_array[1] = seq_array[3];
        seq_array[2*n-1] = seq_array[2*n-3];

        t = 0;
        do {
            for (j = 0; j < n; j++) {
                z[j] = seq_array[2*j+1];
            }

            // TODO: exclude left/right borders?
            for (j = 1; j < n-1; j++) {
                seq_array[2*j+1] = z[j] + r*(z[j+1] - 2.0*z[j] + z[j-1]);
            }

            t++;
        } while (t < max_steps);
    }

    // TODO: change to scatter all?

    ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr);

    ierr = MPI_Bcast(&seq_surface_size, 1, MPI_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);

    if (rank) {
        ierr = PetscCalloc1(seq_surface_size, &seq_array); CHKERRQ(ierr);
    }
    ierr = MPI_Bcast(seq_array, seq_surface_size, MPIU_SCALAR, 0, PETSC_COMM_WORLD);

    ierr = DMDAGetCorners(da_Veloc, &si, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);
    ierr = DMSwarmGetLocalSize(dms_s, &nlocal); CHKERRQ(ierr);
    ierr = DMSwarmGetField(dms_s, DMSwarmPICField_coor, &bs, NULL, (void**)&array); CHKERRQ(ierr);

    for (p = 0; p < nlocal; p++) {
        array[2*p] = seq_array[si*dms_s_ppe*2+2*p];
        array[2*p+1] = seq_array[si*dms_s_ppe*2+2*p+1];
    }

    ierr = DMSwarmRestoreField(dms_s, DMSwarmPICField_coor, &bs, NULL, (void **)&array); CHKERRQ(ierr);
    ierr = VecRestoreArray(seq_surface, &seq_array); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode sp_evaluate_surface_processes_2d_sedimentation_only(PetscReal dt)
{
    PetscErrorCode ierr;
    PetscMPIInt rank;

    PetscReal hsl;
    PetscInt p;
    PetscInt bs;
    PetscInt si;
    PetscInt nlocal;
    PetscReal *array;

    PetscFunctionBeginUser;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    hsl = sp_evaluate_adjusted_mean_elevation_with_sea_level();

    ierr = DMDAGetCorners(da_Veloc, &si, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);
    ierr = DMSwarmGetLocalSize(dms_s, &nlocal); CHKERRQ(ierr);
    ierr = DMSwarmGetField(dms_s, DMSwarmPICField_coor, &bs, NULL, (void**)&array); CHKERRQ(ierr);

    for (p = 0; p < nlocal; p++) {
        if (array[2*p+1] < hsl) {
            array[2*p+1] = hsl;
            PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[SP sedimentation_only] [%d] updated surface [%.3e | %d] = %.3e\n", rank, array[2*p], p, hsl);
        }
    }

    PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);

    ierr = DMSwarmRestoreField(dms_s, DMSwarmPICField_coor, &bs, NULL, (void **)&array); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode sp_evaluate_surface_processes_2d_diffusion_sedimentation_only(PetscReal dt)
{
    PetscErrorCode ierr;
    PetscMPIInt rank;

    PetscFunctionBeginUser;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    Ks = Ks / seg_per_ano; // convert to [m2 s-1]

    PetscReal sea_level;
    PetscReal *hw;
    PetscReal *Kcell;
    PetscReal *Kint;
    PetscReal *flux;
    PetscReal dx_s;

    PetscInt n;
    PetscInt i;
    PetscInt j;
    PetscInt si;
    PetscInt bs;
    PetscInt nlocal;
    PetscInt seq_surface_size;

    Vec global_surface;
    Vec seq_surface;
    VecScatter ctx;
    PetscReal *seq_array;
    PetscReal *array;

    ierr = DMSwarmCreateGlobalVectorFromField(dms_s, DMSwarmPICField_coor, &global_surface); CHKERRQ(ierr);
    ierr = VecScatterCreateToZero(global_surface, &ctx, &seq_surface); CHKERRQ(ierr);
    ierr = VecScatterBegin(ctx, global_surface, seq_surface, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterEnd(ctx, global_surface, seq_surface, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterDestroy(&ctx); CHKERRQ(ierr);
    ierr = DMSwarmDestroyGlobalVectorFromField(dms_s, DMSwarmPICField_coor, &global_surface); CHKERRQ(ierr);

    ierr = VecGetSize(seq_surface, &seq_surface_size); CHKERRQ(ierr);
    ierr = VecGetArray(seq_surface, &seq_array); CHKERRQ(ierr);

    n = seq_surface_size/2;

    sea_level = sp_evaluate_adjusted_mean_elevation_with_sea_level();

    PetscInt nt = (PetscInt)ceil((4.0 * Ks * dt) / (dx_s * dx_s));
    PetscReal dt_s = dt / nt;

    if (!rank) {
        ierr = PetscCalloc1(n, &hw); CHKERRQ(ierr);
        ierr = PetscCalloc1(n, &Kcell); CHKERRQ(ierr);
        ierr = PetscCalloc1(n, &Kint); CHKERRQ(ierr);
        ierr = PetscCalloc1(n, &flux); CHKERRQ(ierr);

        dx_s = seq_array[2];

        for (i = 0; i < nt; i++) {

            for (j = 0; j < n; j++) {
                hw[j] = PetscMax(0.0, sea_level - seq_array[2*j+1]);
                Kcell[j] = Ks * PetscExpReal(-lambda_s * hw[j]);
            }

            for (j = 1; j < n; j++) {
                Kint[j] = 0.5 * (Kcell[j-1] + Kcell[j]);
                flux[j] = -Kint[j] * (seq_array[2*j+1] - seq_array[2*(j-1)+1]) / dx_s;
            }

            flux[0] = 0.0;
            flux[n-1] = 0.0;

            for (j = 0; j < n; j++) {
                seq_array[2*j+1] += (dt_s / dx_s) * (flux[j+1] - flux[j]);
            }
        }
    }


    ierr = MPI_Bcast(&seq_surface_size, 1, MPI_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);

    if (rank) {
        ierr = PetscCalloc1(seq_surface_size, &seq_array); CHKERRQ(ierr);
    }
    ierr = MPI_Bcast(seq_array, seq_surface_size, MPIU_SCALAR, 0, PETSC_COMM_WORLD);

    ierr = DMDAGetCorners(da_Veloc, &si, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);
    ierr = DMSwarmGetLocalSize(dms_s, &nlocal); CHKERRQ(ierr);
    ierr = DMSwarmGetField(dms_s, DMSwarmPICField_coor, &bs, NULL, (void**)&array); CHKERRQ(ierr);

    for (j = 0; j < nlocal; j++) {
        array[2*j] = seq_array[si*dms_s_ppe*2+2*j];
        array[2*j+1] = seq_array[si*dms_s_ppe*2+2*j+1];
    }

    ierr = DMSwarmRestoreField(dms_s, DMSwarmPICField_coor, &bs, NULL, (void **)&array); CHKERRQ(ierr);
    ierr = VecRestoreArray(seq_surface, &seq_array); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode sp_update_surface_swarm_particles_properties()
{
    PetscMPIInt rank;
	PetscErrorCode ierr;

    PetscInt i;
    PetscInt p;
    PetscInt bs;
    PetscInt mx;
    PetscInt sex;
    PetscInt nlocal;
    PetscInt *layer;
    PetscReal px;
    PetscReal py;
    PetscReal rx;
    PetscReal sp_dx;
    PetscReal *pcoords;
	PetscReal *geoq_fac;
	PetscReal *strain_fac;

    Vec global_surface;
    Vec seq_surface;
    VecScatter ctx;
    PetscInt global_surface_size;
    PetscInt seq_surface_size;
    PetscReal *seq_array;

    PetscReal surface;
    PetscReal _sl;
    PetscReal _sr;

    PetscFunctionBeginUser;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

	ierr = DMSwarmGetLocalSize(dms, &nlocal); CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms, DMSwarmPICField_coor, &bs, NULL, (void**)&pcoords); CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms, "layer", &bs, NULL, (void**)&layer); CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms, "geoq_fac", &bs, NULL, (void**)&geoq_fac); CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms, "strain_fac", &bs, NULL, (void**)&strain_fac); CHKERRQ(ierr);

    ierr = DMSwarmCreateGlobalVectorFromField(dms_s, DMSwarmPICField_coor, &global_surface); CHKERRQ(ierr);
    ierr = VecScatterCreateToAll(global_surface, &ctx, &seq_surface); CHKERRQ(ierr);
    ierr = VecScatterBegin(ctx, global_surface, seq_surface, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterEnd(ctx, global_surface, seq_surface, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterDestroy(&ctx); CHKERRQ(ierr);
    ierr = DMSwarmDestroyGlobalVectorFromField(dms_s, DMSwarmPICField_coor, &global_surface); CHKERRQ(ierr);

    ierr = DMSwarmGetSize(dms_s, &global_surface_size); CHKERRQ(ierr);
    ierr = DMSwarmGetLocalSize(dms_s, &seq_surface_size); CHKERRQ(ierr);
    ierr = VecGetArray(seq_surface, &seq_array); CHKERRQ(ierr);

    ierr = DMDAGetElementCorners_2d(da_Veloc, &sex, NULL, &mx, NULL); CHKERRQ(ierr);

    PetscInt si;
    PetscInt milocal;
    ierr = DMDAGetCorners(da_Veloc, &si, NULL, NULL, &milocal, NULL, NULL); CHKERRQ(ierr);

    sp_dx = Lx/((global_surface_size-1));

    for (p = 0; seq_surface_size && p < nlocal; p++) {
        px = pcoords[2*p];
		py = pcoords[2*p+1];

        i = (PetscInt)(px/sp_dx);

        rx = 0.5;

        _sl = seq_array[2*i+1];
        _sr = seq_array[2*i+3];

        surface = _sl * (1.0 - rx) + _sr * rx;

        // (A2L) air particle bellow the surface => assign land properties
        if (a2l == PETSC_TRUE && layer[p] == n_interfaces) {
            if (py < surface) {
                layer[p] = n_interfaces - 1;
                geoq_fac[p] = inter_geoq[n_interfaces - 1];
                strain_fac[p] = 0.0;

                PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[%d] particle updated - a2l - sedimentation | px=%.3e py=%.3e surface=%.3e\n", rank, px, py, surface);
            }
        }

        // (L2A) land particle above the surface => assign air properties
        if (layer[p] != n_interfaces) { //  && sp_mode != SP_SEDIMENTATION_ONLY) {
            if (py > surface) {
                layer[p] = n_interfaces;
                geoq_fac[p] = inter_geoq[n_interfaces];
                strain_fac[p] = 0.0;
            }
        }
    }

    PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);

	ierr = DMSwarmRestoreField(dms, "layer", &bs, NULL, (void**)&layer);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms, "geoq_fac", &bs, NULL, (void**)&geoq_fac);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms, "strain_fac", &bs, NULL, (void**)&strain_fac);CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(dms, DMSwarmPICField_coor, &bs, NULL, (void**)&pcoords); CHKERRQ(ierr);
    ierr = VecRestoreArray(seq_surface, &seq_array); CHKERRQ(ierr);
    pcoords = NULL;

    ierr = DMSwarmMigrate(dms, PETSC_TRUE); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode sp_view_2d(DM dm, const char prefix[])
{
    PetscMPIInt rank;
	PetscErrorCode ierr;
    PetscInt p;
	FILE *fp;
	char name[PETSC_MAX_PATH_LEN];

    Vec global_surface;
    Vec seq_surface;
    VecScatter ctx;
    PetscInt seq_surface_size;
    PetscReal *seq_array;
    PetscInt n;

    PetscFunctionBeginUser;

	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

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
        PetscSNPrintf(name, PETSC_MAX_PATH_LEN-1, "%s.txt", prefix);
        fp = fopen(name, "w");

        if (fp == NULL) {
            SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Cannot open file %s", name);
        }

        for (p = 0; p < n; p++) {
            fprintf(fp, "%+1.5e %+1.5e\n", seq_array[2*p], seq_array[2*p+1]);
        }

        fclose(fp);
    }

    ierr = VecRestoreArray(seq_surface, &seq_array); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}


PetscErrorCode sp_destroy()
{
	PetscErrorCode ierr;

    PetscFunctionBeginUser;

	ierr = DMDestroy(&dmcell_s);CHKERRQ(ierr);
	ierr = DMDestroy(&dms_s);CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
