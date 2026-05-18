#include <petscviewerhdf5.h>
#include <petscvec.h>
#include <petscdmswarm.h>
#include <petscdmda.h>
#include <petscis.h>


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
        sum += ly[i];
    }

    // count Pz
    sum = 0;
    for (k = 0; sum < P; k++) {
        sum += lz[i];
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

PetscErrorCode save_swarm_field(
    DM dms,
    PetscViewer viewer,
    const char *fieldname
)
{
    Vec vec;
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = DMSwarmCreateGlobalVectorFromField(dms, fieldname, &vec); CHKERRQ(ierr);

    ierr = PetscObjectSetName((PetscObject)vec, fieldname); CHKERRQ(ierr);
    ierr = VecView(vec, viewer); CHKERRQ(ierr);

    ierr = DMSwarmDestroyGlobalVectorFromField(dms, fieldname, &vec); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode save_swarm_int_field(
    DM dms,
    PetscViewer viewer,
    const char *fieldname
)
{
    Vec vec;
    const PetscInt *iarray;
    PetscScalar *varray;
    PetscInt nlocal, i;
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = DMSwarmGetLocalSize(dms, &nlocal); CHKERRQ(ierr);

    ierr = DMSwarmGetField(dms, fieldname, NULL, NULL, (void**)&iarray); CHKERRQ(ierr);

    // Create standard Vec (PETSc scalar)
    ierr = VecCreateMPI(PETSC_COMM_WORLD, nlocal, PETSC_DECIDE, &vec); CHKERRQ(ierr);

    ierr = VecGetArray(vec, &varray); CHKERRQ(ierr);

    // Convert int to scalar
    for (i = 0; i < nlocal; i++) {
        varray[i] = (PetscScalar)iarray[i];
    }

    ierr = VecRestoreArray(vec, &varray); CHKERRQ(ierr);

    ierr = PetscObjectSetName((PetscObject)vec, fieldname); CHKERRQ(ierr);
    ierr = VecView(vec, viewer); CHKERRQ(ierr);

    ierr = VecDestroy(&vec); CHKERRQ(ierr);

    ierr = DMSwarmRestoreField(dms, fieldname, NULL, NULL, (void**)&iarray); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode save_particles_to_snapshot(DM dms, PetscViewer viewer, PetscBool magmatism_flag)
{
    Vec all_local;
    PetscScalar *all_nlocal_array_aux;
    PetscInt *all_nlocal_array;
    PetscInt i, npoints_local;
    PetscMPIInt rank, size;
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    // -- write fields
    ierr = PetscViewerHDF5PushGroup(viewer, "/particle_fields"); CHKERRQ(ierr);

    // -- particles metadata
    ierr = DMSwarmGetLocalSize(dms, &npoints_local); CHKERRQ(ierr);
    PetscMalloc1(size, &all_nlocal_array);
    PetscMalloc1(size, &all_nlocal_array_aux);
    MPI_Allgather(&npoints_local, 1, MPIU_INT, all_nlocal_array, 1, MPIU_INT, PETSC_COMM_WORLD);

    for (i = 0; i < size; i++) {
        all_nlocal_array_aux[i] = (PetscScalar)all_nlocal_array[i];
    }

    ierr = VecCreateSeqWithArray(PETSC_COMM_SELF, 1, size, all_nlocal_array_aux, &all_local); CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)all_local, "npoints_local"); CHKERRQ(ierr);

    ierr = VecView(all_local, viewer); CHKERRQ(ierr);

    ierr = VecDestroy(&all_local); CHKERRQ(ierr);
    ierr = PetscFree(all_nlocal_array); CHKERRQ(ierr);
    ierr = PetscFree(all_nlocal_array_aux); CHKERRQ(ierr);

    // -- fields
    ierr = save_swarm_int_field(dms, viewer, "itag"); CHKERRQ(ierr);
    ierr = save_swarm_int_field(dms, viewer, "layer"); CHKERRQ(ierr);
    ierr = save_swarm_int_field(dms, viewer, "cont"); CHKERRQ(ierr);

    ierr = save_swarm_field(dms, viewer, DMSwarmPICField_coor); CHKERRQ(ierr);
    ierr = save_swarm_field(dms, viewer, "geoq_fac"); CHKERRQ(ierr);
    ierr = save_swarm_field(dms, viewer, "strain_fac"); CHKERRQ(ierr);
    ierr = save_swarm_field(dms, viewer, "strain_rate_fac"); CHKERRQ(ierr);

    // -- optional fields
    if (magmatism_flag) {
        ierr = save_swarm_field(dms, viewer, "X"); CHKERRQ(ierr);
        ierr = save_swarm_field(dms, viewer, "Phi"); CHKERRQ(ierr);
        ierr = save_swarm_field(dms, viewer, "dPhi"); CHKERRQ(ierr);
    }

    ierr = PetscViewerHDF5PopGroup(viewer); CHKERRQ(ierr);

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
    PetscBool sp_surface_tracking
)
{
    PetscErrorCode ierr;
    PetscViewerFormat format = PETSC_VIEWER_HDF5_PETSC;
    PetscViewer viewer;

    PetscFunctionBeginUser;

    char filename[PETSC_MAX_PATH_LEN];

    ierr = PetscSNPrintf(filename, PETSC_MAX_PATH_LEN-1, "snapshot_%d.h5", step); CHKERRQ(ierr);

    // to do: check if file exists

    ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD, filename, FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);
    ierr = PetscViewerPushFormat(viewer, format); CHKERRQ(ierr);

    // -- Simulation metadata
    ierr = PetscViewerHDF5PushGroup(viewer, "/metadata"); CHKERRQ(ierr);

    // -- state
    ierr = PetscViewerHDF5PushGroup(viewer, "simulation"); CHKERRQ(ierr);
    ierr = PetscViewerHDF5WriteAttribute(viewer, NULL, "step", PETSC_INT, &step); CHKERRQ(ierr);
    ierr = PetscViewerHDF5WriteAttribute(viewer, NULL, "time", PETSC_REAL, &time); CHKERRQ(ierr);
    ierr = PetscViewerHDF5WriteAttribute(viewer, NULL, "dt", PETSC_REAL, &dt); CHKERRQ(ierr);
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
    ierr = PetscViewerHDF5PushGroup(viewer, "processor_layout"); CHKERRQ(ierr);
    ierr = PetscViewerHDF5WriteAttribute(viewer, NULL, "Px", PETSC_INT, &Px); CHKERRQ(ierr);
    ierr = PetscViewerHDF5WriteAttribute(viewer, NULL, "Pz", PETSC_INT, &Pz); CHKERRQ(ierr);
    ierr = PetscViewerHDF5PopGroup(viewer); CHKERRQ(ierr);

    // -- dm ownership
    ierr = PetscViewerHDF5PushGroup(viewer, "dm_ownership"); CHKERRQ(ierr);

    ierr = PetscObjectSetName((PetscObject)is_lx, "lx"); CHKERRQ(ierr);
    ierr = ISView(is_lx, viewer); CHKERRQ(ierr);

    ierr = PetscObjectSetName((PetscObject)is_lz, "lz"); CHKERRQ(ierr);
    ierr = ISView(is_lz, viewer); CHKERRQ(ierr);

    ierr = PetscViewerHDF5PopGroup(viewer); CHKERRQ(ierr);

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
    ierr = save_particles_to_snapshot(dms, viewer, magmatism_flag); CHKERRQ(ierr);

    if (sp_surface_tracking) {
        ierr = save_surface_particles_to_snapshot(dms_s, viewer); CHKERRQ(ierr);
    }
    //

    // to do: add litho; active state for bc velocity, sed layers, sed rate, basal level

    ierr = PetscViewerPopFormat(viewer); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);

    PetscFunctionReturn(0);
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
    const PetscInt *dm_lx,
    const PetscInt *dm_lz
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

    PetscFunctionBeginUser;

    ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD, filename, FILE_MODE_READ, &viewer); CHKERRQ(ierr);
    ierr = PetscViewerPushFormat(viewer, format); CHKERRQ(ierr);

    // -- simulation metadata
    ierr = PetscViewerHDF5PushGroup(viewer, "/metadata/simulation"); CHKERRQ(ierr);

    ierr = PetscViewerHDF5ReadAttribute(viewer, NULL, "step", PETSC_INT, NULL, &step_aux); CHKERRQ(ierr);
    ierr = PetscViewerHDF5ReadAttribute(viewer, NULL, "time", PETSC_REAL, NULL, &time_aux); CHKERRQ(ierr);
    ierr = PetscViewerHDF5ReadAttribute(viewer, NULL, "dt", PETSC_REAL, NULL, &dt_aux);   CHKERRQ(ierr);

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
    IS is_lx, is_lz;

    ierr = PetscViewerHDF5PushGroup(viewer, "/metadata/dm_ownership"); CHKERRQ(ierr);

    ISCreate(PETSC_COMM_SELF, &is_lx);
    PetscObjectSetName((PetscObject)is_lx, "lx");
    ISLoad(is_lx, viewer);

    ISCreate(PETSC_COMM_SELF, &is_lz);
    PetscObjectSetName((PetscObject)is_lz, "lz");
    ISLoad(is_lz, viewer);

    ierr = PetscViewerHDF5PopGroup(viewer); CHKERRQ(ierr);

    ISGetIndices(is_lx, &dm_lx);
    ISGetIndices(is_lz, &dm_lz);

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
    PetscViewer viewer;
    Vec vec;
    const PetscScalar *array;
    PetscMPIInt rank;
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD, filename, FILE_MODE_READ, &viewer); CHKERRQ(ierr);

    ierr = PetscViewerHDF5PushGroup(viewer, "/particle_fields"); CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_WORLD, &vec); CHKERRQ(ierr);

    ierr = PetscObjectSetName((PetscObject)vec, "npoints_local"); CHKERRQ(ierr);

    ierr = VecLoad(vec, viewer); CHKERRQ(ierr);

    ierr = VecGetArrayRead(vec, &array); CHKERRQ(ierr);

    *nlocal = (PetscInt)array[0];

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
