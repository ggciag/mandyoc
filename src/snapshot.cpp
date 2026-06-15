#include <petscviewerhdf5.h>
#include <petscvec.h>
#include <petscdmswarm.h>
#include <petscdmda.h>
#include <petscis.h>


typedef struct {
    PetscBool is_snapshot;
    PetscInt  max_snapshots;
} SnapshotOptions;

typedef enum {
    PARTICLE_OUTPUT_FULL,
    PARTICLE_OUTPUT_FILTERED
} ParticleOutputMode;

PetscErrorCode update_snapshot_log(
    const char *new_snapshot,
    PetscInt max_snapshots
)
{
    PetscMPIInt rank;

    FILE *fp;

    char *entries[1024];
    PetscInt count = 0;
    PetscInt i;

    char line[PETSC_MAX_PATH_LEN];

    PetscFunctionBeginUser;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    if (rank != 0) {
        PetscFunctionReturn(0);
    }

    fp = fopen("snapshots.log", "r");

    if (fp) {
        while (fgets(line, sizeof(line), fp)) {
            line[strcspn(line, "\n")] = '\0';

            entries[count] = strdup(line);

            count++;
        }

        fclose(fp);
    }

    entries[count] = strdup(new_snapshot);
    count++;

    while (count > max_snapshots) {
        remove(entries[0]);
        free(entries[0]);

        for (i = 1; i < count; i++) {
            entries[i-1] = entries[i];
        }

        count--;
    }

    fp = fopen("snapshots.log.tmp", "w");

    if (!fp) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_WRITE, "Cannot open snapshots.log.tmp");
    }

    for (i = 0; i < count; i++) {
        fprintf(fp, "%s\n", entries[i]);
        free(entries[i]);
    }

    fclose(fp);

    remove("snapshots.log");
    rename("snapshots.log.tmp", "snapshots.log");

    PetscFunctionReturn(0);
}

PetscErrorCode get_processor_partitioning(
    DM da,
    PetscInt *Px,
    PetscInt *Py,
    PetscInt *Pz
)
{
    PetscErrorCode ierr;
    PetscInt dim, M, N, P;
    PetscInt m, n, p;

    const PetscInt *lx, *ly, *lz;

    PetscFunctionBegin;

    // try direct retrieval first
    ierr = DMDAGetInfo(da, &dim, &M, &N, &P, &m, &n, &p, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);

    if (dim == 2) {
        if (m > 0 && n > 0) {
            *Px = m;
            *Pz = n;
            *Py = PETSC_DECIDE;
            PetscFunctionReturn(0);
        }
    } else {
        if (m > 0 && n > 0 && p > 0) {
            *Px = m;
            *Py = n;
            *Pz = p;
            PetscFunctionReturn(0);
        }
    }

    // manually compute partitioning from ownership ranges
    ierr = DMDAGetOwnershipRanges(da, &lx, &ly, &lz); CHKERRQ(ierr);

    PetscInt sum, i, j, k;

    // count Px
    sum = 0;
    for (i = 0; sum < M; i++) {
        sum += lx[i];
    }

    // count Py
    sum = 0;
    for (j = 0; sum < N; j++) {
        sum += ly[j];
    }

    // count Pz
    sum = 0;
    for (k = 0; sum < P; k++) {
        sum += lz[k];
    }

    if (dim == 2) {
        *Px = i;
        *Py = PETSC_DECIDE;
        *Pz = j;
    } else {
        *Px = i;
        *Py = j;
        *Pz = k;
    }

    PetscFunctionReturn(0);
}

PetscErrorCode select_particles(
    DM dms,
    PetscBool plot_sediment,
    int n_interfaces,
    ParticleOutputMode mode,
    PetscInt **selected,
    PetscInt *nselected
)
{
    PetscErrorCode ierr;

    PetscInt npoints;
    PetscInt *itag;
    PetscInt *layer;

    PetscInt *idx;
    PetscInt count = 0;
    PetscInt p;

    PetscFunctionBeginUser;

    ierr = DMSwarmGetLocalSize(dms, &npoints); CHKERRQ(ierr);
    ierr = DMSwarmGetField(dms, "itag", NULL, NULL, (void**)&itag); CHKERRQ(ierr);
    ierr = DMSwarmGetField(dms, "layer", NULL, NULL, (void**)&layer); CHKERRQ(ierr);

    ierr = PetscMalloc1(npoints, &idx); CHKERRQ(ierr);

    if (mode == PARTICLE_OUTPUT_FULL) {
        for (p = 0; p < npoints; p++) {
            idx[count] = p;
            count++;
        }
    } else if (mode == PARTICLE_OUTPUT_FILTERED) {
        for (p = 0; p < npoints; p++) {
            if (itag[p] > 9999 || (plot_sediment && layer[p] == n_interfaces - 1)) {
                idx[count] = p;
                count++;
            }
        }
    }

    ierr = DMSwarmRestoreField(dms, "itag", NULL, NULL, (void**)&itag); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(dms, "layer", NULL, NULL, (void**)&layer); CHKERRQ(ierr);

    *selected = idx;
    *nselected = count;

    PetscFunctionReturn(0);
}

PetscErrorCode save_swarm_litho(
    DM dms,
    PetscViewer viewer,
    double dx_const,
    double dz_const
)
{
    PetscErrorCode ierr;

    PetscInt npoints;
    PetscInt bs;
    PetscInt p;
    PetscInt *layer;
    PetscScalar *array;
    PetscScalar *litho;
    PetscScalar *vec_array;
    Vec vec;

    PetscFunctionBeginUser;

    ierr = DMSwarmGetLocalSize(dms, &npoints); CHKERRQ(ierr);
    ierr = DMSwarmGetField(dms, DMSwarmPICField_coor, &bs, NULL, (void**)&array); CHKERRQ(ierr);
    ierr = DMSwarmGetField(dms, "layer", NULL, NULL, (void**)&layer); CHKERRQ(ierr);

    ierr = PetscMalloc1(npoints * 3, &litho); CHKERRQ(ierr);

    for (p = 0; p < npoints; p++) {
        litho[p*3 + 0] = (PetscScalar)((PetscInt)(array[p*bs + 0] * 5.0 / dx_const));
        litho[p*3 + 1] = (PetscScalar)((PetscInt)(array[p*bs + 1] * 5.0 / dz_const));
        litho[p*3 + 2] = (PetscScalar)layer[p];
    }

    ierr = VecCreateMPI(PETSC_COMM_WORLD, npoints * 3, PETSC_DECIDE, &vec); CHKERRQ(ierr);

    ierr = VecGetArray(vec, &vec_array); CHKERRQ(ierr);
    ierr = PetscMemcpy(vec_array, litho, npoints * 3 * sizeof(PetscScalar)); CHKERRQ(ierr);
    ierr = VecRestoreArray(vec, &vec_array); CHKERRQ(ierr);

    ierr = PetscObjectSetName((PetscObject)vec, "litho"); CHKERRQ(ierr);
    ierr = VecView(vec, viewer); CHKERRQ(ierr);

    ierr = VecDestroy(&vec); CHKERRQ(ierr);
    ierr = PetscFree(litho); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(dms, "layer", NULL, NULL, (void**)&layer); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(dms, DMSwarmPICField_coor, &bs, NULL, (void**)&array); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode save_particle_counts(
    DM dms,
    PetscViewer viewer
)
{
    PetscErrorCode ierr;
    PetscMPIInt rank, size;

    PetscInt nlocal;
    PetscInt *counts;
    PetscScalar *array;
    Vec vec;
    PetscInt i;
    PetscInt local_size;

    PetscFunctionBeginUser;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    ierr = DMSwarmGetLocalSize(dms,&nlocal); CHKERRQ(ierr);

    ierr = PetscMalloc1(size, &counts); CHKERRQ(ierr);
    MPI_Allgather(&nlocal, 1, MPIU_INT, counts, 1, MPIU_INT, PETSC_COMM_WORLD);

    local_size = (rank == 0) ? size : 0;
    ierr = VecCreateMPI(PETSC_COMM_WORLD, local_size, size, &vec); CHKERRQ(ierr);

    if (rank == 0) {
        ierr = VecGetArray(vec, &array); CHKERRQ(ierr);

        for (i = 0; i < size; i++) {
            array[i] = (PetscScalar)counts[i];
        }

        ierr = VecRestoreArray(vec, &array); CHKERRQ(ierr);
    }

    ierr = PetscObjectSetName((PetscObject)vec, "npoints_local"); CHKERRQ(ierr);
    ierr = VecView(vec, viewer); CHKERRQ(ierr);

    ierr = VecDestroy(&vec); CHKERRQ(ierr);
    ierr = PetscFree(counts); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode save_swarm_field(
    DM dms,
    PetscViewer viewer,
    const char *fieldname,
    PetscInt *selected,
    PetscInt nselected
)
{
    PetscErrorCode ierr;

    PetscScalar *buffer;
    PetscScalar *field_array;
    PetscScalar *vec_array;
    Vec vec;

    PetscInt bs;
    PetscInt i, j;

    PetscFunctionBeginUser;

    ierr = DMSwarmGetField(dms, fieldname, &bs, NULL, (void**)&field_array); CHKERRQ(ierr);

    ierr = PetscMalloc1(nselected * bs, &buffer); CHKERRQ(ierr);

    for (i = 0; i < nselected; i++) {
        PetscInt p = selected[i];

        for (j = 0; j < bs; j++) {
            buffer[i*bs + j] = field_array[p*bs + j];
        }
    }

    ierr = VecCreateMPI(PETSC_COMM_WORLD, nselected * bs, PETSC_DECIDE, &vec); CHKERRQ(ierr);

    ierr = VecGetArray(vec, &vec_array); CHKERRQ(ierr);
    ierr = PetscMemcpy(vec_array, buffer, nselected * bs * sizeof(PetscScalar)); CHKERRQ(ierr);
    ierr = VecRestoreArray(vec, &vec_array); CHKERRQ(ierr);

    ierr = PetscObjectSetName((PetscObject)vec, fieldname); CHKERRQ(ierr);
    ierr = VecView(vec, viewer); CHKERRQ(ierr);

    ierr = VecDestroy(&vec); CHKERRQ(ierr);
    ierr = PetscFree(buffer); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(dms, fieldname, &bs, NULL, (void**)&field_array); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode save_swarm_int_field(
    DM dms,
    PetscViewer viewer,
    const char *fieldname,
    PetscInt *selected,
    PetscInt nselected
)
{
    PetscErrorCode ierr;

    PetscScalar *buffer;
    PetscInt *field_array;
    PetscScalar *vec_array;
    Vec vec;

    PetscInt i;

    PetscFunctionBeginUser;

    ierr = DMSwarmGetField(dms, fieldname, NULL, NULL, (void**)&field_array); CHKERRQ(ierr);

    ierr = PetscMalloc1(nselected, &buffer); CHKERRQ(ierr);

    for (i = 0; i < nselected; i++) {
        buffer[i] = (PetscScalar)field_array[selected[i]];
    }

    ierr = VecCreateMPI(PETSC_COMM_WORLD, nselected, PETSC_DECIDE, &vec); CHKERRQ(ierr);

    ierr = VecGetArray(vec, &vec_array); CHKERRQ(ierr);
    ierr = PetscMemcpy(vec_array, buffer, nselected * sizeof(PetscScalar)); CHKERRQ(ierr);
    ierr = VecRestoreArray(vec, &vec_array); CHKERRQ(ierr);

    ierr = PetscObjectSetName((PetscObject)vec, fieldname); CHKERRQ(ierr);
    ierr = VecView(vec, viewer); CHKERRQ(ierr);

    ierr = VecDestroy(&vec); CHKERRQ(ierr);
    ierr = PetscFree(buffer); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(dms, fieldname, NULL, NULL, (void**)&field_array); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode save_particles_to_snapshot(
    DM dms,
    PetscViewer viewer,
    PetscBool magmatism_flag,
    PetscBool export_lithology,
    PetscBool plot_sediment,
    PetscBool is_snapshot,
    int n_interfaces,
    double dx_const,
    double dz_const,
    ParticleOutputMode mode
)
{
    PetscErrorCode ierr;
    PetscMPIInt rank, size;

    PetscInt *selected;
    PetscInt nselected;

    PetscFunctionBeginUser;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    ierr = select_particles(dms, plot_sediment, n_interfaces, mode, &selected, &nselected); CHKERRQ(ierr);

    // -- write fields
    ierr = PetscViewerHDF5PushGroup(viewer, "/particle_fields"); CHKERRQ(ierr);

    // -- particles metadata
    if (mode == PARTICLE_OUTPUT_FULL) {
        ierr = save_particle_counts(dms, viewer); CHKERRQ(ierr);
    }

    // -- fields
    ierr = save_swarm_int_field(dms, viewer, "itag", selected, nselected); CHKERRQ(ierr);
    ierr = save_swarm_int_field(dms, viewer, "layer", selected, nselected); CHKERRQ(ierr);
    ierr = save_swarm_int_field(dms, viewer, "cont", selected, nselected); CHKERRQ(ierr);

    ierr = save_swarm_field(dms, viewer, DMSwarmPICField_coor, selected, nselected); CHKERRQ(ierr);
    ierr = save_swarm_field(dms, viewer, "geoq_fac", selected, nselected); CHKERRQ(ierr);
    ierr = save_swarm_field(dms, viewer, "strain_fac", selected, nselected); CHKERRQ(ierr);
    ierr = save_swarm_field(dms, viewer, "strain_rate_fac", selected, nselected); CHKERRQ(ierr);

    // -- optional fields
    if (magmatism_flag) {
        ierr = save_swarm_field(dms, viewer, "X", selected, nselected); CHKERRQ(ierr);
        ierr = save_swarm_field(dms, viewer, "Phi", selected, nselected); CHKERRQ(ierr);
        ierr = save_swarm_field(dms, viewer, "dPhi", selected, nselected); CHKERRQ(ierr);
    }

    if (export_lithology && is_snapshot == PETSC_FALSE) {
        ierr = save_swarm_litho(dms, viewer, dx_const, dz_const); CHKERRQ(ierr);
    }

    ierr = PetscViewerHDF5PopGroup(viewer); CHKERRQ(ierr);

    ierr = PetscFree(selected); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode save_surface_particles_to_snapshot(
    DM dms_s,
    PetscViewer viewer
)
{
    Vec surface;
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = DMSwarmCreateGlobalVectorFromField(dms_s, DMSwarmPICField_coor, &surface); CHKERRQ(ierr);

    ierr = PetscObjectSetName((PetscObject)surface, "surface"); CHKERRQ(ierr);
    ierr = VecView(surface, viewer); CHKERRQ(ierr);

    ierr = DMSwarmDestroyGlobalVectorFromField(dms_s, DMSwarmPICField_coor, &surface); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode write_hdf5(
    int step,
    double time,
    double dt,
    long nx,
    long nz,
    double lx,
    double lz,
    PetscInt Px,
    PetscInt Pz,
    PetscInt cont_sediment_layer,
    PetscInt active_sediment_layer,
    PetscInt cont_sedimentation_rate,
    PetscReal sedimentation_rate,
    PetscInt cont_bl_level,
    PetscInt variable_baselevel,
    PetscInt cont_var_bcv,
    IS is_lx,
    IS is_lz,
    Vec velocity,
    Vec temperature,
    Vec pressure,
    Vec viscosity,
    Vec density,
    Vec heat,
    Vec strain,
    Vec strain_rate,
    Vec thermal_diffusivity,
    Vec X_depletion,
    Vec Phi,
    Vec dPhi,
    DM dms,
    DM dms_s,
    PetscBool magmatism_flag,
    PetscBool sp_surface_tracking,
    PetscBool export_lithology,
    PetscBool plot_sediment,
    int n_interfaces,
    ParticleOutputMode mode,
    const SnapshotOptions *opts
)
{
    PetscErrorCode ierr;
    PetscMPIInt rank;
    PetscViewerFormat format = PETSC_VIEWER_HDF5_PETSC;
    PetscViewer viewer;

    PetscFunctionBeginUser;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    char filename[PETSC_MAX_PATH_LEN];
    char tmp_filename[PETSC_MAX_PATH_LEN];

    if (opts->is_snapshot) {
        ierr = PetscSNPrintf(filename, PETSC_MAX_PATH_LEN-1, "snapshot_step_%06d_t_%.1fMyr.h5", step, time/1.0e6); CHKERRQ(ierr);
    } else {
        ierr = PetscSNPrintf(filename, PETSC_MAX_PATH_LEN-1, "output_step_%06d.h5", step); CHKERRQ(ierr);
    }
    ierr = PetscSNPrintf(tmp_filename, PETSC_MAX_PATH_LEN-1, "%s.tmp", filename); CHKERRQ(ierr);

    if (rank == 0) {
        remove(tmp_filename);
    }

    // to do: check if file exists

    MPI_Barrier(PETSC_COMM_WORLD);

    ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD, tmp_filename, FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);
    ierr = PetscViewerPushFormat(viewer, format); CHKERRQ(ierr);
    ierr = PetscViewerHDF5SetCollective(viewer, PETSC_TRUE); CHKERRQ(ierr);

    // -- Simulation metadata
    ierr = PetscViewerHDF5PushGroup(viewer, "/metadata"); CHKERRQ(ierr);

    // -- state
    ierr = PetscViewerHDF5PushGroup(viewer, "simulation"); CHKERRQ(ierr);
    ierr = PetscViewerHDF5WriteAttribute(viewer, NULL, "step", PETSC_INT, &step); CHKERRQ(ierr);
    ierr = PetscViewerHDF5WriteAttribute(viewer, NULL, "time", PETSC_REAL, &time); CHKERRQ(ierr);
    ierr = PetscViewerHDF5WriteAttribute(viewer, NULL, "dt", PETSC_REAL, &dt); CHKERRQ(ierr);
    if (opts->is_snapshot) {
        ierr = PetscViewerHDF5WriteAttribute(viewer, NULL, "cont_sediment_layer", PETSC_INT, &cont_sediment_layer); CHKERRQ(ierr);
        ierr = PetscViewerHDF5WriteAttribute(viewer, NULL, "active_sediment_layer", PETSC_INT, &active_sediment_layer); CHKERRQ(ierr);
        ierr = PetscViewerHDF5WriteAttribute(viewer, NULL, "cont_sedimentation_rate", PETSC_INT, &cont_sedimentation_rate); CHKERRQ(ierr);
        ierr = PetscViewerHDF5WriteAttribute(viewer, NULL, "sedimentation_rate", PETSC_REAL, &sedimentation_rate); CHKERRQ(ierr);
        ierr = PetscViewerHDF5WriteAttribute(viewer, NULL, "cont_bl_level", PETSC_INT, &cont_bl_level); CHKERRQ(ierr);
        ierr = PetscViewerHDF5WriteAttribute(viewer, NULL, "variable_baselevel", PETSC_INT, &variable_baselevel); CHKERRQ(ierr);
        ierr = PetscViewerHDF5WriteAttribute(viewer, NULL, "cont_var_bcv", PETSC_INT, &cont_var_bcv); CHKERRQ(ierr);
    }
    ierr = PetscViewerHDF5PopGroup(viewer); CHKERRQ(ierr);
    //

    // -- mesh
    ierr = PetscViewerHDF5PushGroup(viewer, "mesh"); CHKERRQ(ierr);
    ierr = PetscViewerHDF5WriteAttribute(viewer, NULL, "nx", PETSC_INT, &nx); CHKERRQ(ierr);
    ierr = PetscViewerHDF5WriteAttribute(viewer, NULL, "nz", PETSC_INT, &nz); CHKERRQ(ierr);
    ierr = PetscViewerHDF5WriteAttribute(viewer, NULL, "lx", PETSC_REAL, &lx); CHKERRQ(ierr);
    ierr = PetscViewerHDF5WriteAttribute(viewer, NULL, "lz", PETSC_REAL, &lz); CHKERRQ(ierr);
    ierr = PetscViewerHDF5PopGroup(viewer); CHKERRQ(ierr);
    //

    // -- processor layout
    if (opts->is_snapshot) {
        ierr = PetscViewerHDF5PushGroup(viewer, "processor_layout"); CHKERRQ(ierr);
        ierr = PetscViewerHDF5WriteAttribute(viewer, NULL, "Px", PETSC_INT, &Px); CHKERRQ(ierr);
        ierr = PetscViewerHDF5WriteAttribute(viewer, NULL, "Pz", PETSC_INT, &Pz); CHKERRQ(ierr);
        ierr = PetscViewerHDF5PopGroup(viewer); CHKERRQ(ierr);

        ierr = PetscViewerHDF5PushGroup(viewer, "dm_ownership"); CHKERRQ(ierr);
        {
            Vec vec_lx;
            Vec vec_lz;

            const PetscInt *lx_idx;
            const PetscInt *lz_idx;

            PetscInt n_lx;
            PetscInt n_lz;

            PetscInt i;
            PetscInt local_size;
            PetscScalar *array;

            ierr = ISGetLocalSize(is_lx, &n_lx); CHKERRQ(ierr);
            ierr = ISGetLocalSize(is_lz, &n_lz); CHKERRQ(ierr);

            ierr = ISGetIndices(is_lx, &lx_idx); CHKERRQ(ierr);
            ierr = ISGetIndices(is_lz, &lz_idx); CHKERRQ(ierr);

            local_size = (rank == 0) ? n_lx : 0;
            ierr = VecCreateMPI(PETSC_COMM_WORLD, local_size, n_lx, &vec_lx); CHKERRQ(ierr);

            if (rank == 0) {
                ierr = VecGetArray(vec_lx, &array); CHKERRQ(ierr);

                for (i = 0; i < n_lx; i++) {
                    array[i] = (PetscScalar)lx_idx[i];
                }

                ierr = VecRestoreArray(vec_lx, &array); CHKERRQ(ierr);
            }

            ierr = PetscObjectSetName((PetscObject)vec_lx, "lx"); CHKERRQ(ierr);
            ierr = VecView(vec_lx, viewer); CHKERRQ(ierr);

            ierr = VecDestroy(&vec_lx); CHKERRQ(ierr);

            local_size = (rank == 0) ? n_lz : 0;
            ierr = VecCreateMPI(PETSC_COMM_WORLD, local_size, n_lz, &vec_lz); CHKERRQ(ierr);

            if (rank == 0) {
                ierr = VecGetArray(vec_lz, &array); CHKERRQ(ierr);

                for (i = 0; i < n_lz; i++) {
                    array[i] = (PetscScalar)lz_idx[i];
                }

                ierr = VecRestoreArray(vec_lz, &array); CHKERRQ(ierr);
            }

            ierr = PetscObjectSetName((PetscObject)vec_lz,"lz"); CHKERRQ(ierr);
            ierr = VecView(vec_lz, viewer); CHKERRQ(ierr);

            ierr = VecDestroy(&vec_lz); CHKERRQ(ierr);

            ierr = ISRestoreIndices(is_lx, &lx_idx); CHKERRQ(ierr);
            ierr = ISRestoreIndices(is_lz, &lz_idx); CHKERRQ(ierr);
        }
        ierr = PetscViewerHDF5PopGroup(viewer); CHKERRQ(ierr);
    }

    ierr = PetscViewerHDF5PopGroup(viewer); CHKERRQ(ierr);
    // (metadata)

    // -- fields
    ierr = PetscViewerHDF5PushGroup(viewer, "/fields"); CHKERRQ(ierr);

    #define SAVE_VEC(v,name) \
        ierr = PetscObjectSetName((PetscObject)v, name); CHKERRQ(ierr); \
        ierr = VecView(v, viewer); CHKERRQ(ierr);

    SAVE_VEC(velocity, "velocity");
    SAVE_VEC(temperature, "temperature");
    SAVE_VEC(pressure, "pressure");
    SAVE_VEC(viscosity, "viscosity");
    SAVE_VEC(density, "density");
    SAVE_VEC(heat, "heat");
    SAVE_VEC(strain, "strain");
    SAVE_VEC(strain_rate, "strain_rate");
    SAVE_VEC(thermal_diffusivity, "thermal_diffusivity");

    if (magmatism_flag) {
        SAVE_VEC(X_depletion, "X_depletion");
        SAVE_VEC(Phi, "Phi");
        SAVE_VEC(dPhi, "dPhi");
    }

    ierr = PetscViewerHDF5PopGroup(viewer); CHKERRQ(ierr); CHKERRQ(ierr);
    //

    // -- particles
    ierr = save_particles_to_snapshot(dms, viewer, magmatism_flag, export_lithology, plot_sediment, opts->is_snapshot, n_interfaces, lx/(nx-1), -lz/(nz-1), mode); CHKERRQ(ierr);

    if (sp_surface_tracking) {
        ierr = save_surface_particles_to_snapshot(dms_s, viewer); CHKERRQ(ierr);
    }
    //

    ierr = PetscViewerPopFormat(viewer); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);


    MPI_Barrier(PETSC_COMM_WORLD);

    if (rank == 0) {
        int r = rename(tmp_filename, filename);
        if (r != 0) {
            SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_WRITE, "Failed to rename snapshot temp file");
        }
    }

    MPI_Barrier(PETSC_COMM_WORLD);


    if (opts->is_snapshot) {
        ierr = update_snapshot_log(filename, opts->max_snapshots); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}

PetscErrorCode save_snapshot(
    int step,
    double time,
    double dt,
    long nx,
    long nz,
    double lx,
    double lz,
    PetscInt Px,
    PetscInt Pz,
    PetscInt cont_sediment_layer,
    PetscInt active_sediment_layer,
    PetscInt cont_sedimentation_rate,
    PetscReal sedimentation_rate,
    PetscInt cont_bl_level,
    PetscInt variable_baselevel,
    PetscInt cont_var_bcv,
    IS is_lx,
    IS is_lz,
    Vec velocity,
    Vec temperature,
    Vec pressure,
    Vec viscosity,
    Vec density,
    Vec heat,
    Vec strain,
    Vec strain_rate,
    Vec thermal_diffusivity,
    Vec X_depletion,
    Vec Phi,
    Vec dPhi,
    DM dms,
    DM dms_s,
    PetscBool magmatism_flag,
    PetscBool sp_surface_tracking,
    PetscBool export_lithology,
    PetscBool plot_sediment,
    int n_interfaces,
    PetscInt max_snapshots
)
{
    SnapshotOptions opts;

    opts.is_snapshot = PETSC_TRUE;
    opts.max_snapshots = max_snapshots;

    return write_hdf5(
        step,
        time,
        dt,
        nx,
        nz,
        lx,
        lz,
        Px,
        Pz,
        cont_sediment_layer,
        active_sediment_layer,
        cont_sedimentation_rate,
        sedimentation_rate,
        cont_bl_level,
        variable_baselevel,
        cont_var_bcv,
        is_lx,
        is_lz,
        velocity,
        temperature,
        pressure,
        viscosity,
        density,
        heat,
        strain,
        strain_rate,
        thermal_diffusivity,
        X_depletion,
        Phi,
        dPhi,
        dms,
        dms_s,
        magmatism_flag,
        sp_surface_tracking,
        export_lithology,
        PETSC_FALSE,
        0,
        PARTICLE_OUTPUT_FULL,
        &opts
    );
}

PetscErrorCode save_hdf5(
    int step,
    double time,
    double dt,
    long nx,
    long nz,
    double lx,
    double lz,
    PetscInt Px,
    PetscInt Pz,
    IS is_lx,
    IS is_lz,
    Vec velocity,
    Vec temperature,
    Vec pressure,
    Vec viscosity,
    Vec density,
    Vec heat,
    Vec strain,
    Vec strain_rate,
    Vec thermal_diffusivity,
    Vec X_depletion,
    Vec Phi,
    Vec dPhi,
    DM dms,
    DM dms_s,
    PetscBool magmatism_flag,
    PetscBool sp_surface_tracking,
    PetscBool export_lithology,
    PetscBool plot_sediment,
    int n_interfaces
)
{
    SnapshotOptions opts;

    opts.is_snapshot = PETSC_FALSE;
    opts.max_snapshots = 0;

    return write_hdf5(
        step,
        time,
        dt,
        nx,
        nz,
        lx,
        lz,
        Px,
        Pz,
        0,
        0,
        0,
        0.0,
        0,
        0,
        0,
        is_lx,
        is_lz,
        velocity,
        temperature,
        pressure,
        viscosity,
        density,
        heat,
        strain,
        strain_rate,
        thermal_diffusivity,
        X_depletion,
        Phi,
        dPhi,
        dms,
        dms_s,
        magmatism_flag,
        sp_surface_tracking,
        export_lithology,
        plot_sediment,
        n_interfaces,
        PARTICLE_OUTPUT_FILTERED,
        &opts
    );
}

PetscErrorCode load_snapshot_metadata(
    const char *filename,
    int *step,
    double *time,
    double *dt,
    long *nx,
    long *nz,
    double *lx,
    double *lz,
    PetscInt *Px,
    PetscInt *Pz,
    PetscInt *cont_sediment_layer,
    PetscInt *active_sediment_layer,
    PetscInt *cont_sedimentation_rate,
    PetscReal *sedimentation_rate,
    PetscInt *cont_bl_level,
    PetscInt *variable_baselevel,
    PetscInt *cont_var_bcv,
    PetscInt **dm_lx,
    PetscInt **dm_lz
)
{
    PetscErrorCode ierr;
    PetscViewer viewer;
    PetscViewerFormat format = PETSC_VIEWER_HDF5_PETSC;

    PetscInt step_aux;
    PetscReal time_aux, dt_aux;
    PetscInt nx_aux;
    PetscInt nz_aux;
    PetscReal lx_aux;
    PetscReal lz_aux;
    PetscInt Px_aux;
    PetscInt Pz_aux;
    PetscInt cont_sediment_layer_aux;
    PetscInt active_sediment_layer_aux;
    PetscInt cont_sedimentation_rate_aux;
    PetscReal sedimentation_rate_aux;
    PetscInt cont_bl_level_aux;
    PetscInt variable_baselevel_aux;
    PetscInt cont_var_bcv_aux;

    PetscFunctionBeginUser;

    ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD, filename, FILE_MODE_READ, &viewer); CHKERRQ(ierr);
    ierr = PetscViewerPushFormat(viewer, format); CHKERRQ(ierr);

    // -- simulation metadata
    ierr = PetscViewerHDF5PushGroup(viewer, "/metadata/simulation"); CHKERRQ(ierr);

    ierr = PetscViewerHDF5ReadAttribute(viewer, NULL, "step", PETSC_INT, NULL, &step_aux); CHKERRQ(ierr);
    ierr = PetscViewerHDF5ReadAttribute(viewer, NULL, "time", PETSC_REAL, NULL, &time_aux); CHKERRQ(ierr);
    ierr = PetscViewerHDF5ReadAttribute(viewer, NULL, "dt", PETSC_REAL, NULL, &dt_aux);   CHKERRQ(ierr);
    ierr = PetscViewerHDF5ReadAttribute(viewer, NULL, "cont_sediment_layer", PETSC_INT, NULL, &cont_sediment_layer_aux); CHKERRQ(ierr);
    ierr = PetscViewerHDF5ReadAttribute(viewer, NULL, "active_sediment_layer", PETSC_INT, NULL, &active_sediment_layer_aux); CHKERRQ(ierr);
    ierr = PetscViewerHDF5ReadAttribute(viewer, NULL, "cont_sedimentation_rate", PETSC_INT, NULL, &cont_sedimentation_rate_aux); CHKERRQ(ierr);
    ierr = PetscViewerHDF5ReadAttribute(viewer, NULL, "sedimentation_rate", PETSC_REAL, NULL, &sedimentation_rate_aux); CHKERRQ(ierr);
    ierr = PetscViewerHDF5ReadAttribute(viewer, NULL, "cont_bl_level", PETSC_INT, NULL, &cont_bl_level_aux); CHKERRQ(ierr);
    ierr = PetscViewerHDF5ReadAttribute(viewer, NULL, "variable_baselevel", PETSC_INT, NULL, &variable_baselevel_aux); CHKERRQ(ierr);
    ierr = PetscViewerHDF5ReadAttribute(viewer, NULL, "cont_var_bcv", PETSC_INT, NULL, &cont_var_bcv_aux); CHKERRQ(ierr);
    ierr = PetscViewerHDF5PopGroup(viewer); CHKERRQ(ierr);

    // -- Mesh metadata
    ierr = PetscViewerHDF5PushGroup(viewer, "/metadata/mesh"); CHKERRQ(ierr);

    ierr = PetscViewerHDF5ReadAttribute(viewer, NULL, "nx", PETSC_INT, NULL, &nx_aux); CHKERRQ(ierr);
    ierr = PetscViewerHDF5ReadAttribute(viewer, NULL, "nz", PETSC_INT, NULL, &nz_aux); CHKERRQ(ierr);
    ierr = PetscViewerHDF5ReadAttribute(viewer, NULL, "lx", PETSC_REAL, NULL, &lx_aux); CHKERRQ(ierr);
    ierr = PetscViewerHDF5ReadAttribute(viewer, NULL, "lz", PETSC_REAL, NULL, &lz_aux); CHKERRQ(ierr);

    ierr = PetscViewerHDF5PopGroup(viewer); CHKERRQ(ierr);

    // -- processor layout
    ierr = PetscViewerHDF5PushGroup(viewer, "/metadata/processor_layout"); CHKERRQ(ierr);

    ierr = PetscViewerHDF5ReadAttribute(viewer, NULL, "Px", PETSC_INT, NULL, &Px_aux); CHKERRQ(ierr);
    ierr = PetscViewerHDF5ReadAttribute(viewer, NULL, "Pz", PETSC_INT, NULL, &Pz_aux); CHKERRQ(ierr);

    ierr = PetscViewerHDF5PopGroup(viewer); CHKERRQ(ierr);


    // -- dm ownership
    ierr = PetscViewerHDF5PushGroup(viewer, "/metadata/dm_ownership"); CHKERRQ(ierr);

    Vec vec_lx;
    Vec vec_lz;

    const PetscScalar *array;

    PetscInt n_lx;
    PetscInt n_lz;
    PetscInt i;

    ierr = VecCreate(PETSC_COMM_SELF, &vec_lx); CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)vec_lx, "lx"); CHKERRQ(ierr);
    ierr = VecLoad(vec_lx, viewer); CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_SELF, &vec_lz); CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)vec_lz, "lz"); CHKERRQ(ierr);
    ierr = VecLoad(vec_lz, viewer); CHKERRQ(ierr);

    ierr = VecGetSize(vec_lx, &n_lx); CHKERRQ(ierr);
    ierr = VecGetSize(vec_lz, &n_lz); CHKERRQ(ierr);

    PetscInt *dm_lx_aux;
    PetscInt *dm_lz_aux;

    PetscMalloc1(n_lx, &dm_lx_aux);
    PetscMalloc1(n_lz, &dm_lz_aux);

    ierr = VecGetArrayRead(vec_lx, &array); CHKERRQ(ierr);
    for (i = 0; i < n_lx; i++) {
        dm_lx_aux[i] = (PetscInt)array[i];
    }
    ierr = VecRestoreArrayRead(vec_lx, &array); CHKERRQ(ierr);
    ierr = VecDestroy(&vec_lx); CHKERRQ(ierr);

    ierr = VecGetArrayRead(vec_lz, &array); CHKERRQ(ierr);
    for (i = 0; i < n_lz; i++) {
        dm_lz_aux[i] = (PetscInt)array[i];
    }
    ierr = VecRestoreArrayRead(vec_lz, &array); CHKERRQ(ierr);
    ierr = VecDestroy(&vec_lz); CHKERRQ(ierr);

    *dm_lx = dm_lx_aux;
    *dm_lz = dm_lz_aux;

    ierr = PetscViewerHDF5PopGroup(viewer); CHKERRQ(ierr);

    ierr = PetscViewerPopFormat(viewer); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

    *step = (int)step_aux;
    *time = (double)time_aux;
    *dt = (double)dt_aux;
    *nx = (long)nx_aux;
    *nz = (long)nz_aux;
    *lx = (double)lx_aux;
    *lz = (double)lz_aux;
    *Px = (PetscInt)Px_aux;
    *Pz = (PetscInt)Pz_aux;
    *cont_sediment_layer = (PetscInt)cont_sediment_layer_aux;
    *active_sediment_layer = (PetscInt)active_sediment_layer_aux;
    *cont_sedimentation_rate = (PetscInt)cont_sedimentation_rate_aux;
    *sedimentation_rate = (PetscReal)sedimentation_rate_aux;
    *cont_bl_level = (PetscInt)cont_bl_level_aux;
    *variable_baselevel = (PetscInt)variable_baselevel_aux;
    *cont_var_bcv = (PetscInt)cont_var_bcv_aux;

    PetscFunctionReturn(0);
}

PetscErrorCode load_snapshot_fields(
    const char *filename,
    Vec velocity,
    Vec temperature,
    Vec pressure,
    Vec viscosity,
    Vec density,
    Vec heat,
    Vec strain,
    Vec strain_rate,
    Vec thermal_diffusivity,
    Vec X_depletion,
    Vec Phi,
    Vec dPhi,
    PetscBool magmatism_flag
)
{
    PetscErrorCode ierr;
    PetscViewer viewer;
    PetscViewerFormat format = PETSC_VIEWER_HDF5_PETSC;

    ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD, filename, FILE_MODE_READ, &viewer); CHKERRQ(ierr);
    ierr = PetscViewerPushFormat(viewer, format); CHKERRQ(ierr);

    ierr = PetscViewerHDF5PushGroup(viewer, "/fields"); CHKERRQ(ierr);

    #define LOAD_VEC(v,name) \
        ierr = PetscObjectSetName((PetscObject)v, name); CHKERRQ(ierr); \
        ierr = VecLoad(v, viewer); CHKERRQ(ierr);

    LOAD_VEC(velocity, "velocity");
    LOAD_VEC(temperature, "temperature");
    LOAD_VEC(pressure, "pressure");
    LOAD_VEC(viscosity, "viscosity");
    LOAD_VEC(density, "density");
    LOAD_VEC(heat, "heat");
    LOAD_VEC(strain, "strain");
    LOAD_VEC(strain_rate, "strain_rate");
    LOAD_VEC(thermal_diffusivity, "thermal_diffusivity");

    if (magmatism_flag) {
        LOAD_VEC(X_depletion, "X_depletion");
        LOAD_VEC(Phi, "Phi");
        LOAD_VEC(dPhi, "dPhi");
    }

    ierr = PetscViewerHDF5PopGroup(viewer); CHKERRQ(ierr);

    ierr = PetscViewerPopFormat(viewer); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode load_swarm_local_sizes(const char *filename, PetscInt *nlocal)
{
    PetscErrorCode ierr;
    PetscMPIInt rank;
    PetscViewer viewer;
    Vec vec;
    const PetscScalar *array;
    PetscInt low, high;

    PetscFunctionBeginUser;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD, filename, FILE_MODE_READ, &viewer); CHKERRQ(ierr);

    ierr = PetscViewerHDF5PushGroup(viewer, "/particle_fields"); CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_WORLD, &vec); CHKERRQ(ierr);

    ierr = PetscObjectSetName((PetscObject)vec, "npoints_local"); CHKERRQ(ierr);

    ierr = VecLoad(vec, viewer); CHKERRQ(ierr);

    ierr = VecGetOwnershipRange(vec, &low, &high); CHKERRQ(ierr);

    ierr = VecGetArrayRead(vec, &array); CHKERRQ(ierr);

    if (high > low) {
        *nlocal = (PetscInt)array[0];
    } else {
        *nlocal = 0;
    }

    ierr = VecRestoreArrayRead(vec, &array); CHKERRQ(ierr);

    ierr = VecDestroy(&vec); CHKERRQ(ierr);

    ierr = PetscViewerHDF5PopGroup(viewer); CHKERRQ(ierr);

    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode load_swarm_int_field(DM dms, PetscViewer viewer, const char *fieldname)
{
    Vec vec;
    PetscInt *field_array;
    const PetscScalar *vec_array;
    PetscInt nlocal, i;
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = DMSwarmGetLocalSize(dms, &nlocal); CHKERRQ(ierr);

    ierr = DMSwarmGetField(dms, fieldname, NULL, NULL, (void**)&field_array); CHKERRQ(ierr);

    ierr = VecCreateMPI(PETSC_COMM_WORLD, nlocal, PETSC_DECIDE, &vec); CHKERRQ(ierr);

    ierr = PetscObjectSetName((PetscObject)vec, fieldname); CHKERRQ(ierr);

    ierr = VecLoad(vec, viewer); CHKERRQ(ierr);

    ierr = VecGetArrayRead(vec, &vec_array); CHKERRQ(ierr);

    for (i = 0; i < nlocal; i++) {
        field_array[i] = (PetscInt)vec_array[i];
    }

    ierr = VecRestoreArrayRead(vec, &vec_array); CHKERRQ(ierr);

    ierr = VecDestroy(&vec); CHKERRQ(ierr);

    ierr = DMSwarmRestoreField(dms, fieldname, NULL, NULL, (void**)&field_array); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode load_swarm_field(DM dms, PetscViewer viewer, const char *fieldname)
{
    Vec vec;
    PetscScalar *field_array;
    const PetscScalar *vec_array;

    PetscInt nlocal, bs;
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = DMSwarmGetLocalSize(dms, &nlocal);CHKERRQ(ierr);

    ierr = DMSwarmGetField(dms, fieldname, &bs, NULL, (void**)&field_array); CHKERRQ(ierr);

    ierr = VecCreateMPI(PETSC_COMM_WORLD, nlocal*bs, PETSC_DECIDE, &vec); CHKERRQ(ierr);

    ierr = PetscObjectSetName((PetscObject)vec, fieldname); CHKERRQ(ierr);

    ierr = VecLoad(vec, viewer);CHKERRQ(ierr);

    ierr = VecGetArrayRead(vec, &vec_array);CHKERRQ(ierr);

    PetscMemcpy(field_array, vec_array, nlocal * bs * sizeof(PetscScalar));

    ierr = VecRestoreArrayRead(vec, &vec_array);CHKERRQ(ierr);

    ierr = VecDestroy(&vec);CHKERRQ(ierr);

    ierr = DMSwarmRestoreField(dms, fieldname, &bs, NULL, (void**)&field_array); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode load_particles_from_snapshot(const char *filename, DM dms, PetscBool magmatism_flag)
{
    PetscViewer viewer;
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD, filename, FILE_MODE_READ, &viewer); CHKERRQ(ierr);

    ierr = PetscViewerHDF5PushGroup(viewer, "/particle_fields"); CHKERRQ(ierr);

    ierr = load_swarm_int_field(dms, viewer, "itag"); CHKERRQ(ierr);
    ierr = load_swarm_int_field(dms, viewer, "layer"); CHKERRQ(ierr);
    ierr = load_swarm_int_field(dms, viewer, "cont"); CHKERRQ(ierr);

    ierr = load_swarm_field(dms, viewer, DMSwarmPICField_coor); CHKERRQ(ierr);
    ierr = load_swarm_field(dms, viewer, "geoq_fac"); CHKERRQ(ierr);
    ierr = load_swarm_field(dms, viewer, "strain_fac"); CHKERRQ(ierr);
    ierr = load_swarm_field(dms, viewer, "strain_rate_fac"); CHKERRQ(ierr);

    if (magmatism_flag) {
        ierr = load_swarm_field(dms, viewer, "X"); CHKERRQ(ierr);
        ierr = load_swarm_field(dms, viewer, "Phi"); CHKERRQ(ierr);
        ierr = load_swarm_field(dms, viewer, "dPhi"); CHKERRQ(ierr);
    }

    ierr = PetscViewerHDF5PopGroup(viewer); CHKERRQ(ierr);

    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
