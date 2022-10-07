#include <petscsys.h>


extern PetscInt Px;
extern PetscInt Pz;
extern int n_interfaces;
extern PetscInt interfaces_from_ascii;
extern PetscInt *seed_layer;
extern PetscInt seed_layer_size;
extern PetscBool seed_layer_set;
extern PetscReal *strain_seed_layer;
extern PetscInt strain_seed_layer_size;
extern PetscBool strain_seed_layer_set;
extern PetscBool strain_seed_constant;
extern PetscBool strain_seed_constant_set;
extern double h_air;;
extern PetscInt pressure_in_rheol;
extern PetscBool sp_surface_processes;
extern PetscBool sp_surface_tracking;
extern PetscInt sp_mode;
extern PetscBool set_sp_d_c;
extern PetscScalar sp_d_c;
extern PetscScalar K_fluvial;
extern PetscScalar sea_level;


PetscErrorCode load_topo_var(int rank);


PetscErrorCode parse_options(int rank)
{
	PetscErrorCode ierr;

	ierr = PetscOptionsGetInt(NULL , NULL, "-Px", &Px, NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL , NULL, "-Pz", &Pz, NULL); CHKERRQ(ierr);

	if (n_interfaces > 0 && interfaces_from_ascii == 1) {
		ierr = PetscCalloc1(n_interfaces, &seed_layer); CHKERRQ(ierr);
		seed_layer_size = n_interfaces + 1;
		ierr = PetscOptionsGetIntArray(NULL, NULL, "-seed", seed_layer, &seed_layer_size, &seed_layer_set); CHKERRQ(ierr);

		ierr = PetscCalloc1(n_interfaces, &strain_seed_layer); CHKERRQ(ierr);
		strain_seed_layer_size = n_interfaces + 1;
		ierr = PetscOptionsGetRealArray(NULL, NULL, "-strain_seed", strain_seed_layer, &strain_seed_layer_size, &strain_seed_layer_set); CHKERRQ(ierr);
		if (strain_seed_layer_set == PETSC_TRUE && seed_layer_set == PETSC_FALSE) {
			PetscPrintf(PETSC_COMM_WORLD, "Specify the seed layer with the flag -seed (required by -strain_seed)\n");
			exit(1);
		}
		if (strain_seed_layer_set == PETSC_TRUE && seed_layer_set == PETSC_TRUE && seed_layer_size != strain_seed_layer_size) {
			PetscPrintf(PETSC_COMM_WORLD, "Specify the same number of values in the list for flags -seed and -strain_seed\n");
			exit(1);
		}
		if (strain_seed_layer_set == PETSC_FALSE && seed_layer_set == PETSC_TRUE) {
			PetscPrintf(PETSC_COMM_WORLD, "Using default value '0.5' for -strain_seed (for all seed layers)\n");
			for (int k = 0; k < seed_layer_size; k++) {
				strain_seed_layer[k] = 0.5;
			}
		}
		PetscPrintf(PETSC_COMM_WORLD, "Number of seed layers: %d\n", seed_layer_size);
		for (int k = 0; k < seed_layer_size; k++) {
			PetscPrintf(PETSC_COMM_WORLD, "seed layer: %d - strain: %lf\n", seed_layer[k], strain_seed_layer[k]);
		}
		PetscPrintf(PETSC_COMM_WORLD, "\n");

		ierr = PetscOptionsGetBool(NULL, NULL, "-strain_seed_constant", &strain_seed_constant, &strain_seed_constant_set); CHKERRQ(ierr);
		if (strain_seed_constant_set == PETSC_TRUE && seed_layer_set == PETSC_FALSE) {
			PetscPrintf(PETSC_COMM_WORLD, "Specify the seed layer with the flags -seed and -strain_seed (required by -strain_seed_constant)\n");
			exit(1);
		}
	}

	h_air = -1.0;
	ierr = PetscOptionsGetReal(NULL, NULL, "-h_air", &h_air, NULL); CHKERRQ(ierr);
	if (pressure_in_rheol == 0 && h_air < 0.0) {
		PetscPrintf(PETSC_COMM_WORLD, "Specify the thickness of the air layer with the flag -h_air\n");
		PetscPrintf(PETSC_COMM_WORLD, "(you adopted depth dependent rheology: pressure_in_rheol = 0)\n");
		exit(1);
	} else {
		h_air = 0.0;
	}

	if (sp_surface_processes && sp_surface_tracking && sp_mode == 1) {
		load_topo_var(rank);
	}

	if (sp_mode == 2 && PETSC_FALSE == set_sp_d_c) {
		ierr = PetscPrintf(PETSC_COMM_WORLD, "-sp_mode 2 (diffusion) using default value: sp_d_c %e\n", sp_d_c); CHKERRQ(ierr);
	} else if (sp_mode == 2) {
		ierr = PetscPrintf(PETSC_COMM_WORLD, "-sp_mode 2 (diffusion) using custom value: sp_d_c %e\n", sp_d_c); CHKERRQ(ierr);
	} else if (sp_mode == 3) {
		ierr = PetscPrintf(PETSC_COMM_WORLD, "-sp_mode 3 (fluvial erosion) using K_fluvial: %e and sea_level %e\n", K_fluvial, sea_level); CHKERRQ(ierr);
	} else if (sp_mode == 4) {
		ierr = PetscPrintf(PETSC_COMM_WORLD, "-sp_mode 4 (fluvial erosion mode 2) using K_fluvial: %e and sea_level %e\n", K_fluvial, sea_level); CHKERRQ(ierr);
	}

    PetscFunctionReturn(0);
}
