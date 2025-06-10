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
extern double h_air;

PetscInt *seed_layer_aux;
PetscReal *strain_seed_layer_aux;


PetscErrorCode parse_options(int rank)
{
	PetscErrorCode ierr;

	PetscFunctionBeginUser;

	ierr = PetscOptionsGetInt(NULL , NULL, "-Px", &Px, NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL , NULL, "-Pz", &Pz, NULL); CHKERRQ(ierr);

	ierr = PetscCalloc1(100000, &seed_layer_aux); CHKERRQ(ierr);
	seed_layer_size = 100000;
	ierr = PetscOptionsGetIntArray(NULL, NULL, "-seed", seed_layer_aux, &seed_layer_size, &seed_layer_set); CHKERRQ(ierr);

	ierr = PetscCalloc1(100000, &strain_seed_layer_aux); CHKERRQ(ierr);
	strain_seed_layer_size = 100000;
	ierr = PetscOptionsGetRealArray(NULL, NULL, "-strain_seed", strain_seed_layer_aux, &strain_seed_layer_size, &strain_seed_layer_set); CHKERRQ(ierr);

	ierr = PetscCalloc1(seed_layer_size, &seed_layer); CHKERRQ(ierr);
	ierr = PetscCalloc1(strain_seed_layer_size, &strain_seed_layer); CHKERRQ(ierr);

	for (int k = 0; k < seed_layer_size; k++) {
		seed_layer[k] = seed_layer_aux[k];
		strain_seed_layer[k] = strain_seed_layer_aux[k];
	}

	PetscFree(seed_layer_aux);
	PetscFree(strain_seed_layer_aux);

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

	ierr = PetscOptionsGetBool(NULL, NULL, "-strain_seed_constant", &strain_seed_constant, &strain_seed_constant_set); CHKERRQ(ierr);
	if (strain_seed_constant_set == PETSC_TRUE && seed_layer_set == PETSC_FALSE) {
		PetscPrintf(PETSC_COMM_WORLD, "Specify the seed layer with the flags -seed and -strain_seed (required by -strain_seed_constant)\n");
		exit(1);
	}

	h_air = -1.0;
	ierr = PetscOptionsGetReal(NULL, NULL, "-h_air", &h_air, NULL); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
