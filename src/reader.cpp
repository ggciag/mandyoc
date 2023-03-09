#include <petscksp.h>

// Prototypes
int check_a_b(char tkn_w[], char tkn_v[], const char str_a[], const char str_b[]);
PetscBool check_a_b_bool(char tkn_w[], char tkn_v[], const char str_a[], const char str_b[]);
//void ErrorInterfaces(int rank, const char fname[], int flag);
void ErrorInterfaces();

// Parameter file variables
extern int dimensions;
extern long Nx, Ny, Nz;
extern long layers;
extern double Lx, Ly, depth;
extern int ContMult;
extern long stepMAX;
extern double timeMAX;
extern double dt_MAX;
extern long print_step;
extern double visco_r;
extern double visc_MAX;
extern double visc_MIN;
extern int geoq_on;
extern double escala_viscosidade;
extern double veloc_superf;
extern double RHOM;
extern double alpha_exp_thermo;
extern double kappa;
extern double gravity;
extern double Delta_T;
extern int n_interfaces;
extern double H_per_mass;
extern double c_heat_capacity;
extern int T_initial_cond;
extern int rheol;
extern int bcv_top_normal;
extern int bcv_top_slip;
extern int bcv_bot_normal;
extern int bcv_bot_slip;
extern int bcv_left_normal;
extern int bcv_left_slip;
extern int bcv_right_normal;
extern int bcv_right_slip;
extern int bcT_top;
extern int bcT_bot;
extern int bcT_left;
extern int bcT_right;
extern PetscInt WITH_NON_LINEAR;
extern PetscInt WITH_ADIABATIC_H;
extern PetscInt WITH_RADIOGENIC_H;
extern PetscInt PLASTICITY;
extern PetscReal denok_min;
extern PetscInt particles_per_ele;
extern PetscReal theta_FSSA;
extern PetscInt direct_solver;
extern PetscInt visc_const_per_element;
extern PetscReal sub_division_time_step;
extern PetscInt visc_harmonic_mean;
extern PetscInt pressure_in_rheol;
extern PetscReal particles_perturb_factor;
extern PetscInt variable_bcv;
extern PetscInt interfaces_from_ascii;
extern PetscReal rtol;
extern PetscInt temper_extern;
extern PetscInt veloc_extern;
extern PetscInt bcv_extern;
extern PetscInt binary_output;
extern PetscInt sticky_blanket_air;
extern PetscInt multi_velocity;
extern PetscInt sp_mode;
extern PetscInt precipitation_profile;
extern PetscInt climate_change;
extern PetscInt free_surface_stab;
extern PetscInt print_step_files;
extern PetscInt RK4;
extern PetscReal Xi_min;
extern PetscReal random_initial_strain;
extern PetscInt checkered;
extern PetscReal pressure_const;
extern PetscInt initial_dynamic_range;
extern PetscInt periodic_boundary;
extern PetscInt nx_ppe;
extern PetscInt ny_ppe;
extern PetscInt nz_ppe;
extern PetscInt initial_print_step;
extern PetscReal initial_print_max_time;
extern PetscScalar K_fluvial;
extern PetscScalar m_fluvial;
extern PetscScalar sea_level;
extern PetscScalar basal_heat;
extern PetscBool sp_surface_tracking;
extern PetscBool sp_surface_processes;
extern PetscReal sp_dt;
extern PetscReal sp_d_c;
extern PetscBool set_sp_dt;
extern PetscBool set_sp_d_c;
extern PetscInt high_kappa_in_asthenosphere;
extern PetscBool plot_sediment;
extern PetscBool a2l;

// Removed from parameter file
extern double H_lito;
extern double beta_max;
extern double ramp_begin;
extern double ramp_end;
extern double h_air;

// Parameters from the interfaces.txt
extern PetscScalar *interfaces;
extern PetscScalar *inter_rho;
extern PetscScalar *inter_geoq;
extern PetscScalar *inter_H;
extern PetscScalar *inter_A;
extern PetscScalar *inter_n;
extern PetscScalar *inter_Q;
extern PetscScalar *inter_V;
extern PetscScalar weakening_min;
extern PetscScalar weakening_max;
extern PetscScalar *weakening_seed;
extern PetscScalar *cohesion_min;
extern PetscScalar *cohesion_max;
extern PetscScalar *friction_angle_min;
extern PetscScalar *friction_angle_max;

extern PetscScalar *mv_time;
extern PetscInt n_mv;

extern PetscScalar *var_bcv_time;
extern PetscScalar *var_bcv_scale;
extern PetscInt n_var_bcv;

extern PetscInt variable_climate;

extern PetscScalar *var_climate_time;
extern PetscScalar *var_climate_scale;
extern PetscInt n_var_climate;

extern PetscInt climate_change;

char str[100];

extern long sp_n_profiles;
extern PetscScalar *topo_var_time;
extern PetscScalar *topo_var_rate;
extern PetscScalar *global_surface_array_helper;
extern PetscScalar *global_surface_array_helper_aux;

extern PetscInt non_dim;

extern PetscReal h0_scaled;
extern PetscReal visc0_scaled;
extern PetscReal g0_scaled;
extern PetscReal rho0_scaled;

extern PetscReal time0_scaled;
extern PetscReal veloc0_scaled;
extern PetscReal kappa0_scaled;
extern PetscReal temperature0_scaled;

extern PetscReal advection_scaled;
extern PetscReal diffusion_scaled;
extern PetscReal radiogenic_scaled;
extern PetscReal adiabatic_scaled;

extern PetscReal pressure0_scaled;
extern PetscReal strain_rate0_scaled;

extern PetscReal air_threshold_density;

PetscErrorCode load_topo_var(int rank);

// Reads input ASCII files
PetscErrorCode reader(int rank, const char fName[]){
	int nline;
	int size = 1024;
	void *nread;
	char *tkn_w, *tkn_v;
	char line[size];
	// char *mandatoryParam = ["nx", "nz"];

	if (rank==0){
		nline = 0;
		FILE *f_parameters;

		f_parameters = fopen(fName,"r");
		if (f_parameters == NULL) {PetscPrintf(PETSC_COMM_WORLD, "ERROR. The <%s> file was NOT FOUND.\n", fName); exit(1);}
		while (!feof(f_parameters))
		{
			// Increment line number nline
			nline += 1;

			// Read each line from the parameters f_parameters, trim undesireble characters
			// and split line in thk_w[] and thk_v[].
			nread = fgets(line, size, f_parameters);
			if ((nread == NULL) || (line[0] == '\n') || (line[0] == '#')) continue;
			tkn_w 	= strtok(line, " \t=\n");
			if (tkn_w[0] == '#') continue;
			tkn_v 	= strtok(NULL, " \t=#\n");
			// PetscPrintf(PETSC_COMM_WORLD, "(%s)(%s)\n", tkn_w, tkn_v);

			// Store and check every parameter.
			// Parameters with values
			if (strcmp(tkn_w, "dimensions") == 0) {dimensions = atoi(tkn_v);}
			else if (strcmp(tkn_w, "nx") == 0) {Nx = atol(tkn_v);}
			else if (strcmp(tkn_w, "ny") == 0) {Ny = atol(tkn_v);}
			else if (strcmp(tkn_w, "nz") == 0) {Nz = atol(tkn_v);}
			else if (strcmp(tkn_w, "lx") == 0) {Lx = atof(tkn_v);}
			else if (strcmp(tkn_w, "ly") == 0) {Ly = atof(tkn_v);}
			else if (strcmp(tkn_w, "lz") == 0) {depth = atof(tkn_v);}
			else if (strcmp(tkn_w, "multigrid") == 0) {ContMult = atoi(tkn_v);}
			else if (strcmp(tkn_w, "step_max") == 0) {stepMAX = atol(tkn_v);}
			else if (strcmp(tkn_w, "step_print") == 0) {print_step = atol(tkn_v);}
			else if (strcmp(tkn_w, "time_max") == 0) {timeMAX = atof(tkn_v);}
			else if (strcmp(tkn_w, "dt_max") == 0) {dt_MAX = atof(tkn_v);}
			else if (strcmp(tkn_w, "viscosity_reference") == 0) {visco_r = atof(tkn_v);}
			else if (strcmp(tkn_w, "viscosity_max") == 0) {visc_MAX = atof(tkn_v);}
			else if (strcmp(tkn_w, "viscosity_min") == 0) {visc_MIN = atof(tkn_v);}
			else if (strcmp(tkn_w, "n_interfaces") == 0) {n_interfaces = atoi(tkn_v);}
			else if (strcmp(tkn_w, "geoq_fac") == 0) {escala_viscosidade = atof(tkn_v);}
			else if (strcmp(tkn_w, "surface_velocity") == 0) {veloc_superf = atof(tkn_v);}
			else if (strcmp(tkn_w, "temperature_difference") == 0) {Delta_T = atof(tkn_v);}
			else if (strcmp(tkn_w, "thermal_expansion_coefficient") == 0) {alpha_exp_thermo = atof(tkn_v);}
			else if (strcmp(tkn_w, "thermal_diffusivity_coefficient") == 0) {kappa = atof(tkn_v);}
			else if (strcmp(tkn_w, "gravity_acceleration") == 0) {gravity = atof(tkn_v);}
			else if (strcmp(tkn_w, "density_mantle") == 0) {RHOM = atof(tkn_v);}
			else if (strcmp(tkn_w, "external_heat") == 0) {H_per_mass = atof(tkn_v);}
			else if (strcmp(tkn_w, "heat_capacity") == 0) {c_heat_capacity = atof(tkn_v);}
			else if (strcmp(tkn_w, "rheology_model") == 0) {rheol = atoi(tkn_v);}
			else if (strcmp(tkn_w, "T_initial") == 0) {T_initial_cond = atoi(tkn_v);}
			else if (strcmp(tkn_w, "denok") == 0) {denok_min = atof(tkn_v);}
			else if (strcmp(tkn_w, "particles_per_element") == 0) {particles_per_ele = atol(tkn_v);}
			else if (strcmp(tkn_w, "theta_FSSA") == 0) {theta_FSSA = atof(tkn_v);}
			else if (strcmp(tkn_w, "sub_division_time_step") == 0) {sub_division_time_step = atof(tkn_v);}
			else if (strcmp(tkn_w, "particles_perturb_factor") == 0) {particles_perturb_factor = atof(tkn_v);}
			else if (strcmp(tkn_w, "rtol") == 0) {rtol = atof(tkn_v);}
			else if (strcmp(tkn_w, "sp_mode") == 0) {sp_mode = atoi(tkn_v);}
			else if (strcmp(tkn_w, "Xi_min") == 0) {Xi_min = atof(tkn_v);}
			else if (strcmp(tkn_w, "random_initial_strain") == 0) {random_initial_strain = atof(tkn_v);}
			else if (strcmp(tkn_w, "pressure_const") == 0) {pressure_const = atof(tkn_v);}
			else if (strcmp(tkn_w, "particles_per_element_x") == 0) {nx_ppe = atoi(tkn_v);}
			else if (strcmp(tkn_w, "particles_per_element_y") == 0) {ny_ppe = atoi(tkn_v);}
			else if (strcmp(tkn_w, "particles_per_element_z") == 0) {nz_ppe = atoi(tkn_v);}
			else if (strcmp(tkn_w, "initial_print_step") == 0) {initial_print_step = atoi(tkn_v);}
			else if (strcmp(tkn_w, "initial_print_max_time") == 0) {initial_print_max_time = atof(tkn_v);}
			else if (strcmp(tkn_w, "K_fluvial") == 0) {K_fluvial = atof(tkn_v);}
			else if (strcmp(tkn_w, "m_fluvial") == 0) {m_fluvial = atof(tkn_v);}
			else if (strcmp(tkn_w, "sea_level") == 0) {sea_level = atof(tkn_v);}
			else if (strcmp(tkn_w, "basal_heat") == 0) {basal_heat = atof(tkn_v);}
			else if (strcmp(tkn_w, "sp_dt") == 0) {sp_dt = atof(tkn_v);}
			else if (strcmp(tkn_w, "sp_d_c") == 0) {sp_d_c = atof(tkn_v);}
			else if (strcmp(tkn_w, "weakening_min") == 0) {weakening_min = atof(tkn_v);}
			else if (strcmp(tkn_w, "weakening_max") == 0) {weakening_max = atof(tkn_v);}

			// Boolean parameters
			else if (strcmp(tkn_w, "geoq") == 0) {geoq_on = check_a_b(tkn_w, tkn_v, "on", "off");}
			else if (strcmp(tkn_w, "non_linear_method") == 0) {WITH_NON_LINEAR = check_a_b(tkn_w, tkn_v, "on", "off");}
			else if (strcmp(tkn_w, "adiabatic_component") == 0) {WITH_ADIABATIC_H = check_a_b(tkn_w, tkn_v, "on", "off");}
			else if (strcmp(tkn_w, "radiogenic_component") == 0) {WITH_RADIOGENIC_H = check_a_b(tkn_w, tkn_v, "on", "off");}
			else if (strcmp(tkn_w, "plasticity") == 0) {PLASTICITY = check_a_b(tkn_w, tkn_v, "on", "off");}
			else if (strcmp(tkn_w, "top_normal_velocity") == 0) {bcv_top_normal = check_a_b(tkn_w, tkn_v, "fixed", "free");}
			else if (strcmp(tkn_w, "bot_normal_velocity") == 0) {bcv_bot_normal = check_a_b(tkn_w, tkn_v, "fixed", "free");}
			else if (strcmp(tkn_w, "left_normal_velocity") == 0) {bcv_left_normal = check_a_b(tkn_w, tkn_v, "fixed", "free");}
			else if (strcmp(tkn_w, "right_normal_velocity") == 0) {bcv_right_normal = check_a_b(tkn_w, tkn_v, "fixed", "free");}
			else if (strcmp(tkn_w, "top_tangential_velocity") == 0) {bcv_top_slip = check_a_b(tkn_w, tkn_v, "fixed", "free");}
			else if (strcmp(tkn_w, "bot_tangential_velocity") == 0) {bcv_bot_slip = check_a_b(tkn_w, tkn_v, "fixed", "free");}
			else if (strcmp(tkn_w, "left_tangential_velocity") == 0) {bcv_left_slip = check_a_b(tkn_w, tkn_v, "fixed", "free");}
			else if (strcmp(tkn_w, "right_tangential_velocity") == 0) {bcv_right_slip = check_a_b(tkn_w, tkn_v, "fixed", "free");}
			else if (strcmp(tkn_w, "top_temperature") == 0) {bcT_top = check_a_b(tkn_w, tkn_v, "fixed", "free");}
			else if (strcmp(tkn_w, "bot_temperature") == 0) {bcT_bot = check_a_b(tkn_w, tkn_v, "fixed", "free");}
			else if (strcmp(tkn_w, "left_temperature") == 0) {bcT_left = check_a_b(tkn_w, tkn_v, "fixed", "free");}
			else if (strcmp(tkn_w, "right_temperature") == 0) {bcT_right = check_a_b(tkn_w, tkn_v, "fixed", "free");}
			else if (strcmp(tkn_w, "interfaces_from_ascii") == 0) {interfaces_from_ascii = check_a_b(tkn_w, tkn_v, "True", "False");}
			else if (strcmp(tkn_w, "temperature_from_ascii") == 0) {temper_extern = check_a_b(tkn_w, tkn_v, "True", "False");}
			else if (strcmp(tkn_w, "velocity_from_ascii") == 0) {veloc_extern = check_a_b(tkn_w, tkn_v, "True", "False");}
			else if (strcmp(tkn_w, "solver") == 0) {direct_solver = check_a_b(tkn_w, tkn_v, "direct", "iterative");}
			else if (strcmp(tkn_w, "viscosity_per_element") == 0) {visc_const_per_element = check_a_b(tkn_w, tkn_v, "constant", "variable");}
			else if (strcmp(tkn_w, "viscosity_mean_method") == 0) {visc_harmonic_mean = check_a_b(tkn_w, tkn_v, "harmonic", "arithmetic");}
			else if (strcmp(tkn_w, "viscosity_dependence") == 0) {pressure_in_rheol = check_a_b(tkn_w, tkn_v, "pressure", "depth");}
			else if (strcmp(tkn_w, "variable_bcv") == 0) {variable_bcv = check_a_b(tkn_w, tkn_v, "True", "False");}
			else if (strcmp(tkn_w, "bcv_extern") == 0) {bcv_extern = check_a_b(tkn_w, tkn_v, "True", "False");}
			else if (strcmp(tkn_w, "binary_output") == 0) {binary_output = check_a_b(tkn_w, tkn_v, "True", "False");}
			else if (strcmp(tkn_w, "sticky_blanket_air") == 0) {sticky_blanket_air = check_a_b(tkn_w, tkn_v, "True", "False");}
			else if (strcmp(tkn_w, "multi_velocity") == 0) {multi_velocity = check_a_b(tkn_w, tkn_v, "True", "False");}
			else if (strcmp(tkn_w, "precipitation_profile_from_ascii") == 0) {precipitation_profile = check_a_b(tkn_w, tkn_v, "True", "False");}
			else if (strcmp(tkn_w, "climate_change_from_ascii") == 0) {climate_change = check_a_b(tkn_w, tkn_v, "True", "False");}
			else if (strcmp(tkn_w, "free_surface_stab") == 0) {free_surface_stab = check_a_b(tkn_w, tkn_v, "True", "False");}
			else if (strcmp(tkn_w, "print_step_files") == 0) {print_step_files = check_a_b(tkn_w, tkn_v, "True", "False");}
			else if (strcmp(tkn_w, "RK4") == 0) {RK4 = check_a_b(tkn_w, tkn_v, "Runge-Kutta", "Euler");}
			else if (strcmp(tkn_w, "checkered") == 0) {checkered = check_a_b(tkn_w, tkn_v, "True", "False");}
			else if (strcmp(tkn_w, "initial_dynamic_range") == 0) {initial_dynamic_range = check_a_b(tkn_w, tkn_v, "True", "False");}
			else if (strcmp(tkn_w, "periodic_boundary") == 0) {periodic_boundary = check_a_b(tkn_w, tkn_v, "True", "False");}
			else if (strcmp(tkn_w, "sp_surface_tracking") == 0) {sp_surface_tracking = check_a_b_bool(tkn_w, tkn_v, "True", "False");}
			else if (strcmp(tkn_w, "sp_surface_processes") == 0) {sp_surface_processes = check_a_b_bool(tkn_w, tkn_v, "True", "False");}
			else if (strcmp(tkn_w, "plot_sediment") == 0) {plot_sediment = check_a_b_bool(tkn_w, tkn_v, "True", "False");}
			else if (strcmp(tkn_w, "a2l") == 0) {a2l = check_a_b_bool(tkn_w, tkn_v, "True", "False");}

			else if (strcmp(tkn_w, "high_kappa_in_asthenosphere") == 0) {high_kappa_in_asthenosphere = check_a_b(tkn_w, tkn_v, "True", "False");}

			else if (strcmp(tkn_w, "nondimensionalization") == 0) {non_dim = check_a_b(tkn_w, tkn_v, "True", "False");}
			/*else if (strcmp(tkn_w, "h0_scaled") == 0) {h0_scaled = atof(tkn_v);}
			else if (strcmp(tkn_w, "visc0_scaled") == 0) {visc0_scaled = atof(tkn_v);}
			else if (strcmp(tkn_w, "g0_scaled") == 0) {g0_scaled = atof(tkn_v);}
			else if (strcmp(tkn_w, "rho0_scaled") == 0) {rho0_scaled = atof(tkn_v);}

			else if (strcmp(tkn_w, "diffusivity0_scaled") == 0) {kappa0_scaled = atof(tkn_v);}
			else if (strcmp(tkn_w, "temperature0_scaled") == 0) {temperature0_scaled = atof(tkn_v);}*/
			


			// Else
			else
			{
				fprintf(stderr, "Error. Unrecognized keyword <%s> on line <%d> in the parameter file.\n", tkn_w, nline);
				exit(1);
			}
		}
		fclose(f_parameters);

		// check/set dimensions
		if (dimensions != 2 && dimensions != 3) {
			fprintf(stderr, "Error. \"dimensions\" keyword must be 2 (for 2D mode) or 3 (for 3D mode).\n");
			exit(1);
		}
		if (dimensions == 3 && Ny < 2) {
			fprintf(stderr, "Error. \"Ny\" keyword must be a least 2 for 3D mode.\n");
			exit(1);
		}

		// Parameters values exceptions and pre-processing
		if (particles_perturb_factor > 1.0) particles_perturb_factor = 1.0;
		if (particles_perturb_factor < 0.0) particles_perturb_factor = 0.0;
		layers = Nz;
		if (nx_ppe < 0) nx_ppe = 0;
		if (ny_ppe < 0) ny_ppe = 0;
		if (nz_ppe < 0) nz_ppe = 0;
		if (initial_print_step < 0) initial_print_step = 0;
		if (sp_dt > 0) set_sp_dt = PETSC_TRUE;
		if (sp_d_c > 0) set_sp_d_c = PETSC_TRUE;

		/*
		fscanf(f_parameters,"%s",str);
		if (strcmp (str,"H_lito") == 0) fscanf(f_parameters,"%lf",&H_lito);
		else {printf("H_lito error\n"); exit(1);}

		fscanf(f_parameters,"%s",str);
		if (strcmp (str,"h_air") == 0) fscanf(f_parameters,"%lf",&h_air);
		else {printf("h_air error\n"); exit(1);}

		fscanf(f_parameters,"%s",str);
		if (strcmp (str,"beta_max") == 0) fscanf(f_parameters,"%lf",&beta_max);
		else {printf("beta_max error\n"); exit(1);}

		fscanf(f_parameters,"%s",str);
		if (strcmp (str,"ramp_begin") == 0) fscanf(f_parameters,"%lf",&ramp_begin);
		else {printf("ramp_begin error\n"); exit(1);}

		fscanf(f_parameters,"%s",str);
		if (strcmp (str,"ramp_end") == 0) fscanf(f_parameters,"%lf",&ramp_end);
		else {printf("ramp_end error\n"); exit(1);}
		*/

	}

	// Broadcast every parameter from the parameter file.
	MPI_Bcast(&dimensions,1,MPI_INT,0,PETSC_COMM_WORLD);
	MPI_Bcast(&Nx,1,MPI_LONG,0,PETSC_COMM_WORLD);
	MPI_Bcast(&Ny,1,MPI_LONG,0,PETSC_COMM_WORLD);
	MPI_Bcast(&Nz,1,MPI_LONG,0,PETSC_COMM_WORLD);
	MPI_Bcast(&Lx,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
	MPI_Bcast(&Ly,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
	MPI_Bcast(&depth,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
	MPI_Bcast(&ContMult,1,MPI_INT,0,PETSC_COMM_WORLD);
	MPI_Bcast(&stepMAX,1,MPI_LONG,0,PETSC_COMM_WORLD);
	MPI_Bcast(&timeMAX,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
	MPI_Bcast(&dt_MAX,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
	MPI_Bcast(&print_step,1,MPI_LONG,0,PETSC_COMM_WORLD);
	MPI_Bcast(&visco_r,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
	MPI_Bcast(&visc_MAX,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
	MPI_Bcast(&visc_MIN,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
	MPI_Bcast(&geoq_on,1,MPI_INT,0,PETSC_COMM_WORLD);
	MPI_Bcast(&escala_viscosidade,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
	MPI_Bcast(&veloc_superf,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
	MPI_Bcast(&Delta_T,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
	MPI_Bcast(&alpha_exp_thermo,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
	MPI_Bcast(&kappa,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
	MPI_Bcast(&gravity,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
	MPI_Bcast(&RHOM,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
	MPI_Bcast(&H_per_mass,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
	MPI_Bcast(&c_heat_capacity,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
	MPI_Bcast(&WITH_NON_LINEAR,1,MPI_INT,0,PETSC_COMM_WORLD);
	MPI_Bcast(&WITH_ADIABATIC_H,1,MPI_INT,0,PETSC_COMM_WORLD);
	MPI_Bcast(&WITH_RADIOGENIC_H,1,MPI_INT,0,PETSC_COMM_WORLD);
	MPI_Bcast(&PLASTICITY,1,MPI_INT,0,PETSC_COMM_WORLD);
	MPI_Bcast(&bcv_top_normal,1,MPI_INT,0,PETSC_COMM_WORLD);
	MPI_Bcast(&bcv_top_slip,1,MPI_INT,0,PETSC_COMM_WORLD);
	MPI_Bcast(&bcv_bot_normal,1,MPI_INT,0,PETSC_COMM_WORLD);
	MPI_Bcast(&bcv_bot_slip,1,MPI_INT,0,PETSC_COMM_WORLD);
	MPI_Bcast(&bcv_left_normal,1,MPI_INT,0,PETSC_COMM_WORLD);
	MPI_Bcast(&bcv_left_slip,1,MPI_INT,0,PETSC_COMM_WORLD);
	MPI_Bcast(&bcv_right_normal,1,MPI_INT,0,PETSC_COMM_WORLD);
	MPI_Bcast(&bcv_right_slip,1,MPI_INT,0,PETSC_COMM_WORLD);
	MPI_Bcast(&bcT_top,1,MPI_INT,0,PETSC_COMM_WORLD);
	MPI_Bcast(&bcT_bot,1,MPI_INT,0,PETSC_COMM_WORLD);
	MPI_Bcast(&bcT_left,1,MPI_INT,0,PETSC_COMM_WORLD);
	MPI_Bcast(&bcT_right,1,MPI_INT,0,PETSC_COMM_WORLD);
	MPI_Bcast(&rheol,1,MPI_INT,0,PETSC_COMM_WORLD);
	MPI_Bcast(&T_initial_cond,1,MPI_INT,0,PETSC_COMM_WORLD);
	MPI_Bcast(&denok_min,1,MPIU_REAL,0,PETSC_COMM_WORLD);
	MPI_Bcast(&particles_per_ele,1,MPI_LONG,0,PETSC_COMM_WORLD);
	MPI_Bcast(&theta_FSSA,1,MPIU_REAL,0,PETSC_COMM_WORLD);
	MPI_Bcast(&sub_division_time_step,1,MPIU_REAL,0,PETSC_COMM_WORLD);
	MPI_Bcast(&particles_perturb_factor,1,MPIU_REAL,0,PETSC_COMM_WORLD);
	MPI_Bcast(&direct_solver,1,MPI_INT,0,PETSC_COMM_WORLD);
	MPI_Bcast(&visc_const_per_element,1,MPI_INT,0,PETSC_COMM_WORLD);
	MPI_Bcast(&visc_harmonic_mean,1,MPI_INT,0,PETSC_COMM_WORLD);
	MPI_Bcast(&pressure_in_rheol,1,MPI_INT,0,PETSC_COMM_WORLD);
	MPI_Bcast(&variable_bcv,1,MPI_INT,0,PETSC_COMM_WORLD);
	MPI_Bcast(&interfaces_from_ascii,1,MPI_INT,0,PETSC_COMM_WORLD);
	MPI_Bcast(&rtol,1,MPIU_REAL,0,PETSC_COMM_WORLD);
	MPI_Bcast(&temper_extern,1,MPI_INT,0,PETSC_COMM_WORLD);
	MPI_Bcast(&veloc_extern,1,MPI_INT,0,PETSC_COMM_WORLD);
	MPI_Bcast(&bcv_extern,1,MPI_INT,0,PETSC_COMM_WORLD);
	MPI_Bcast(&binary_output,1,MPI_INT,0,PETSC_COMM_WORLD);
	MPI_Bcast(&sticky_blanket_air,1,MPI_INT,0,PETSC_COMM_WORLD);
	MPI_Bcast(&multi_velocity,1,MPI_INT,0,PETSC_COMM_WORLD);
	MPI_Bcast(&sp_mode,1,MPI_INT,0,PETSC_COMM_WORLD);
	MPI_Bcast(&precipitation_profile,1,MPI_INT,0,PETSC_COMM_WORLD);
	MPI_Bcast(&climate_change,1,MPI_INT,0,PETSC_COMM_WORLD);
	MPI_Bcast(&free_surface_stab,1,MPI_INT,0,PETSC_COMM_WORLD);
	MPI_Bcast(&print_step_files,1,MPI_INT,0,PETSC_COMM_WORLD);
	MPI_Bcast(&RK4,1,MPI_INT,0,PETSC_COMM_WORLD);
	MPI_Bcast(&Xi_min,1,MPIU_REAL,0,PETSC_COMM_WORLD);
	MPI_Bcast(&random_initial_strain,1,MPIU_REAL,0,PETSC_COMM_WORLD);
	MPI_Bcast(&checkered,1,MPI_INT,0,PETSC_COMM_WORLD);
	MPI_Bcast(&pressure_const,1,MPIU_REAL,0,PETSC_COMM_WORLD);
	MPI_Bcast(&initial_dynamic_range,1,MPI_INT,0,PETSC_COMM_WORLD);
	MPI_Bcast(&periodic_boundary,1,MPI_INT,0,PETSC_COMM_WORLD);
	MPI_Bcast(&nx_ppe,1,MPI_INT,0,PETSC_COMM_WORLD);
	MPI_Bcast(&ny_ppe,1,MPI_INT,0,PETSC_COMM_WORLD);
	MPI_Bcast(&nz_ppe,1,MPI_INT,0,PETSC_COMM_WORLD);
	MPI_Bcast(&initial_print_step,1,MPI_INT,0,PETSC_COMM_WORLD);
	MPI_Bcast(&initial_print_max_time,1,MPIU_REAL,0,PETSC_COMM_WORLD);
	MPI_Bcast(&high_kappa_in_asthenosphere,1,MPI_INT,0,PETSC_COMM_WORLD);
	MPI_Bcast(&n_interfaces,1,MPI_INT,0,PETSC_COMM_WORLD); // Broadcast after interfaces.txt
	MPI_Bcast(&K_fluvial,1,MPI_REAL,0,PETSC_COMM_WORLD);
	MPI_Bcast(&m_fluvial,1,MPI_REAL,0,PETSC_COMM_WORLD);
	MPI_Bcast(&sea_level,1,MPI_REAL,0,PETSC_COMM_WORLD);
	MPI_Bcast(&basal_heat,1,MPI_REAL,0,PETSC_COMM_WORLD);
	MPI_Bcast(&sp_surface_tracking,1,MPI_C_BOOL,0,PETSC_COMM_WORLD);
	MPI_Bcast(&sp_surface_processes,1,MPI_C_BOOL,0,PETSC_COMM_WORLD);
	MPI_Bcast(&sp_dt,1,MPIU_REAL,0,PETSC_COMM_WORLD);
	MPI_Bcast(&sp_d_c,1,MPIU_REAL,0,PETSC_COMM_WORLD);
	MPI_Bcast(&set_sp_dt,1,MPI_C_BOOL,0,PETSC_COMM_WORLD);
	MPI_Bcast(&set_sp_d_c,1,MPI_C_BOOL,0,PETSC_COMM_WORLD);
	MPI_Bcast(&plot_sediment,1,MPI_C_BOOL,0,PETSC_COMM_WORLD);
	MPI_Bcast(&a2l,1,MPI_C_BOOL,0,PETSC_COMM_WORLD);
	MPI_Bcast(&weakening_min,1,MPI_REAL,0,PETSC_COMM_WORLD);
	MPI_Bcast(&weakening_max,1,MPI_REAL,0,PETSC_COMM_WORLD);

	MPI_Bcast(&non_dim,1,MPI_INT,0,PETSC_COMM_WORLD);

	if (non_dim==1){
		h0_scaled = depth;
		depth /= h0_scaled;
		Lx /= h0_scaled;
		if (dimensions==3) Ly /= h0_scaled;

		visc0_scaled = visc_MAX;
		visc_MAX /= visc0_scaled;
		visc_MIN /= visc0_scaled;
		visco_r /= visc0_scaled;

		g0_scaled = gravity;
		gravity /= g0_scaled;

		rho0_scaled = RHOM;
		RHOM /= rho0_scaled;

		kappa0_scaled = kappa;
		kappa /= kappa0_scaled;

	}

	if (rank==0){
		FILE *f_dim;
		f_dim = fopen("dimensional_convertion.txt","w");
		fprintf(f_dim,"depth %lg\n",h0_scaled);
		fprintf(f_dim,"viscosity %lg\n",visc0_scaled);
		fprintf(f_dim,"gravity %lg\n",g0_scaled);
		fprintf(f_dim,"density %lg\n",rho0_scaled);
		fprintf(f_dim,"diffusivity %lg\n",kappa0_scaled);
		fclose(f_dim);
	}

	/*MPI_Bcast(&h0_scaled,1,MPIU_REAL,0,PETSC_COMM_WORLD);
	MPI_Bcast(&visc0_scaled,1,MPIU_REAL,0,PETSC_COMM_WORLD);
	MPI_Bcast(&g0_scaled,1,MPIU_REAL,0,PETSC_COMM_WORLD);
	MPI_Bcast(&rho0_scaled,1,MPIU_REAL,0,PETSC_COMM_WORLD);

	MPI_Bcast(&kappa0_scaled,1,MPIU_REAL,0,PETSC_COMM_WORLD);
	MPI_Bcast(&temperature0_scaled,1,MPIU_REAL,0,PETSC_COMM_WORLD);*/

	pressure0_scaled = rho0_scaled*g0_scaled*h0_scaled;
	strain_rate0_scaled = rho0_scaled*g0_scaled*h0_scaled/visc0_scaled;
	veloc0_scaled = rho0_scaled*g0_scaled*h0_scaled*h0_scaled/visc0_scaled;

	time0_scaled = h0_scaled/veloc0_scaled;
	//if (time0_scaled!=1.0 || veloc0_scaled!=1) non_dim=1;

	advection_scaled = veloc0_scaled*time0_scaled/h0_scaled;
	diffusion_scaled = kappa0_scaled*time0_scaled/h0_scaled/h0_scaled;
	radiogenic_scaled = time0_scaled/temperature0_scaled;
	adiabatic_scaled = time0_scaled*g0_scaled*veloc0_scaled;


	
	air_threshold_density = 100.0/rho0_scaled;



	//MPI_Bcast(&,1,,0,PETSC_COMM_WORLD);

	/*MPI_Bcast(&H_lito,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
	MPI_Bcast(&h_air,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
	MPI_Bcast(&beta_max,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
	MPI_Bcast(&ramp_begin,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
	MPI_Bcast(&ramp_end,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);*/

	if (rank==0){
		if (dimensions == 2) {
			printf("Mesh size:   %ld %ld\n",Nx,Nz);
			printf("Domain size  %lf %lf\n\n",Lx,depth);
		} else {
			printf("Mesh size:   %ld %ld %ld\n",Nx,Ny,Nz);
			printf("Domain size  %lf %lf %lf\n\n",Lx,Ly,depth);
		}
	}

	// Interfaces
	if (n_interfaces>0 && interfaces_from_ascii==1) {
		if (dimensions == 2) {
			PetscCalloc1(Nx*n_interfaces,&interfaces);
		} else {
			PetscCalloc1(Nx*Ny*n_interfaces,&interfaces);
		}
	}
	PetscCalloc1(n_interfaces+1,&inter_geoq);
	PetscCalloc1(n_interfaces+1,&inter_rho);
	PetscCalloc1(n_interfaces+1,&inter_H);
	PetscCalloc1(n_interfaces+1,&inter_A);
	PetscCalloc1(n_interfaces+1,&inter_n);
	PetscCalloc1(n_interfaces+1,&inter_Q);
	PetscCalloc1(n_interfaces+1,&inter_V);
	PetscCalloc1(n_interfaces+1,&weakening_seed);
	PetscCalloc1(n_interfaces+1,&cohesion_min);
	PetscCalloc1(n_interfaces+1,&cohesion_max);
	PetscCalloc1(n_interfaces+1,&friction_angle_min);
	PetscCalloc1(n_interfaces+1,&friction_angle_max);


	// Read interfaces.txt
	if ((interfaces_from_ascii==1) && (rank==0)){
		nline = 0;
//		int nunits;
//		int nheader = 7;
//		char cp_line_1[size], cp_line_2[size];
		FILE *f_interfaces;

		f_interfaces = fopen("interfaces.txt","r");
		if (f_interfaces==NULL) {PetscPrintf(PETSC_COMM_WORLD, "ERROR. The <interfaces.txt> file was NOT FOUND.\n"); exit(1);}

//		while(!feof(f_interfaces))
//		{
//			nunits = 0;
//			// Read each line of the interfaces.txt file, trim undesireble
//			// characters.
//			nread = fgets(line, size, f_interfaces);
//			if (nline<=nheader) strcpy(cp_line_1, line);
//			else strcpy(cp_line_2, line);
//			if ((nread == NULL) || (line[0] == '\n') || (line[0] == '#')) continue;
//			tkn_w 	= strtok(line, " \t\n");
//			while ((tkn_w != NULL) && (strcmp(tkn_w, "#") != 0) && (tkn_w[0] != '#'))
//			{
//				PetscPrintf(PETSC_COMM_WORLD, "(%s)\n", tkn_w);
//				tkn_w = strtok(NULL, " \t\n");
//				nunits += 1;
//			}
//			PetscPrintf(PETSC_COMM_WORLD, "nunits: (%d)\n", nunits);
//			if (nline == 0)
//			{
//				if (n_interfaces < 0) {n_interfaces = nunits - 1;} // Allocate variables based on interfaces.txt file
//				PetscCalloc1(Nx * n_interfaces, &interfaces);
//				PetscCalloc1(n_interfaces + 1, &inter_geoq);
//				PetscCalloc1(n_interfaces + 1, &inter_rho);
//				PetscCalloc1(n_interfaces + 1, &inter_H);
//				PetscCalloc1(n_interfaces + 1, &inter_A);
//				PetscCalloc1(n_interfaces + 1, &inter_n);
//				PetscCalloc1(n_interfaces + 1, &inter_Q);
//				PetscCalloc1(n_interfaces + 1, &inter_V);
//			}
//			else if ((nline > 0) && (nline < nheader))
//			{
//				if (nunits - 1 <= n_interfaces) {ErrorInterfaces(rank, f_interfaces, 0); exit(1);}
//			}
//			else if (nunits != n_interfaces - 1) {ErrorInterfaces(rank, f_interfaces, 1); exit(1);}
//			nline += 1;
//		}
// 		End

		f_interfaces = fopen("interfaces.txt","r");
		int check_fscanf;

		fscanf(f_interfaces,"%s",str);
		if (strcmp (str,"C") == 0){
			for (PetscInt i=0;i<n_interfaces+1;i++){
				check_fscanf = fscanf(f_interfaces,"%lf",&inter_geoq[i]);
				if (check_fscanf==0) {ErrorInterfaces(); exit(1);}
			}
		}
		else { ErrorInterfaces(); exit(1);}

		fscanf(f_interfaces,"%s",str);
		if (strcmp (str,"rho") == 0)
			for (PetscInt i=0;i<n_interfaces+1;i++){
				fscanf(f_interfaces,"%lf",&inter_rho[i]);
				inter_rho[i]/=rho0_scaled;
			}
		else { ErrorInterfaces(); exit(1);}

		fscanf(f_interfaces,"%s",str);
		if (strcmp (str,"H") == 0)
			for (PetscInt i=0;i<n_interfaces+1;i++)
				fscanf(f_interfaces,"%lf",&inter_H[i]);
		else { ErrorInterfaces(); exit(1);}

		fscanf(f_interfaces,"%s",str);
		if (strcmp (str,"A") == 0)
			for (PetscInt i=0;i<n_interfaces+1;i++)
				fscanf(f_interfaces,"%lf",&inter_A[i]);
		else { ErrorInterfaces(); exit(1);}

		fscanf(f_interfaces,"%s",str);
		if (strcmp (str,"n") == 0)
			for (PetscInt i=0;i<n_interfaces+1;i++)
				fscanf(f_interfaces,"%lf",&inter_n[i]);
		else { ErrorInterfaces(); exit(1);}

		fscanf(f_interfaces,"%s",str);
		if (strcmp (str,"Q") == 0)
			for (PetscInt i=0;i<n_interfaces+1;i++)
				fscanf(f_interfaces,"%lf",&inter_Q[i]);
		else { ErrorInterfaces(); exit(1);}

		fscanf(f_interfaces,"%s",str);
		if (strcmp (str,"V") == 0)
			for (PetscInt i=0;i<n_interfaces+1;i++)
				fscanf(f_interfaces,"%lf",&inter_V[i]);
		else { ErrorInterfaces(); exit(1);}

		// Strain softening
		if (WITH_NON_LINEAR == 1 && PLASTICITY == 1)
		{
			fscanf(f_interfaces,"%s",str);
			if (strcmp (str,"weakening_seed") == 0)
				for (PetscInt i=0;i<n_interfaces+1;i++)
					fscanf(f_interfaces,"%lf",&weakening_seed[i]);
			else { ErrorInterfaces(); exit(1);}

			fscanf(f_interfaces,"%s",str);
			if (strcmp (str,"cohesion_min") == 0)
				for (PetscInt i=0;i<n_interfaces+1;i++)
					fscanf(f_interfaces,"%lf",&cohesion_min[i]);
			else { ErrorInterfaces(); exit(1);}

			fscanf(f_interfaces,"%s",str);
			if (strcmp (str,"cohesion_max") == 0)
				for (PetscInt i=0;i<n_interfaces+1;i++)
					fscanf(f_interfaces,"%lf",&cohesion_max[i]);
			else { ErrorInterfaces(); exit(1);}

			fscanf(f_interfaces,"%s",str);
			if (strcmp (str,"friction_angle_min") == 0)
				for (PetscInt i=0;i<n_interfaces+1;i++)
					fscanf(f_interfaces,"%lf",&friction_angle_min[i]);
			else { ErrorInterfaces(); exit(1);}

			fscanf(f_interfaces,"%s",str);
			if (strcmp (str,"friction_angle_max") == 0)
				for (PetscInt i=0;i<n_interfaces+1;i++)
					fscanf(f_interfaces,"%lf",&friction_angle_max[i]);
			else { ErrorInterfaces(); exit(1);}
		}

		if (dimensions == 2) {
			for (PetscInt i=0; i<Nx; i++){
				for (PetscInt j=0; j<n_interfaces; j++){
					fscanf(f_interfaces,"%lf",&interfaces[j*Nx+i]);
					interfaces[j*Nx+i] /= h0_scaled;
				}
			}
		} else {
			for (PetscInt i=0; i<Nx*Ny; i++){
				for (PetscInt j=0; j<n_interfaces; j++){
					fscanf(f_interfaces,"%lf",&interfaces[j*Nx*Ny+i]);
					interfaces[j*Nx*Ny+i] /= h0_scaled;
				}
			}
		}

		printf("Layers:\n  ");
		for (PetscInt i=0;i<n_interfaces+1;i++)
			printf("%9d ",i);

		printf("\nGeoq: ");
		for (PetscInt i=0;i<n_interfaces+1;i++)
			printf("%.3e ",inter_geoq[i]);
		printf("\n Rho: ");
		for (PetscInt i=0;i<n_interfaces+1;i++)
			printf("%.3e ",inter_rho[i]);
		printf("\n   H: ");
		for (PetscInt i=0;i<n_interfaces+1;i++)
			printf("%.3e ",inter_H[i]);
		printf("\n   A: ");
		for (PetscInt i=0;i<n_interfaces+1;i++)
			printf("%.3e ",inter_A[i]);
		printf("\n   n: ");
		for (PetscInt i=0;i<n_interfaces+1;i++)
			printf(" %lf ",inter_n[i]);
		printf("\n   Q: ");
		for (PetscInt i=0;i<n_interfaces+1;i++)
			printf("%.3e ",inter_Q[i]);
		printf("\n   V: ");
		for (PetscInt i=0;i<n_interfaces+1;i++)
			printf("%.3e ",inter_V[i]);
		if (WITH_NON_LINEAR==1 && PLASTICITY==1)
		{
			printf("\nseed: ");
			for (PetscInt i=0;i<n_interfaces+1;i++)
				printf("%.3f ", weakening_seed[i]);
			printf("\ncmin: ");
			for (PetscInt i=0;i<n_interfaces+1;i++)
				printf("%.3f ", cohesion_min[i]);
			printf("\ncmax: ");
			for (PetscInt i=0;i<n_interfaces+1;i++)
				printf("%.3f ", cohesion_max[i]);
			printf("\nfmin: ");
			for (PetscInt i=0;i<n_interfaces+1;i++)
				printf("%.3f ", friction_angle_min[i]);
			printf("\nfmax: ");
			for (PetscInt i=0;i<n_interfaces+1;i++)
				printf("%.3f ", friction_angle_max[i]);
		}
		printf("\n\n");

		fclose(f_interfaces); // Close file
	}

	if (n_interfaces>0 && interfaces_from_ascii==1) {
		if (dimensions == 2) {
			MPI_Bcast(interfaces,Nx*n_interfaces,MPIU_SCALAR,0,PETSC_COMM_WORLD);
		} else {
			MPI_Bcast(interfaces,Nx*Ny*n_interfaces,MPIU_SCALAR,0,PETSC_COMM_WORLD);
		}
	}
	MPI_Bcast(inter_geoq,n_interfaces+1,MPIU_SCALAR,0,PETSC_COMM_WORLD);
	MPI_Bcast(inter_rho,n_interfaces+1,MPIU_SCALAR,0,PETSC_COMM_WORLD);
	MPI_Bcast(inter_H,n_interfaces+1,MPIU_SCALAR,0,PETSC_COMM_WORLD);
	MPI_Bcast(inter_A,n_interfaces+1,MPIU_SCALAR,0,PETSC_COMM_WORLD);
	MPI_Bcast(inter_n,n_interfaces+1,MPIU_SCALAR,0,PETSC_COMM_WORLD);
	MPI_Bcast(inter_Q,n_interfaces+1,MPIU_SCALAR,0,PETSC_COMM_WORLD);
	MPI_Bcast(inter_V,n_interfaces+1,MPIU_SCALAR,0,PETSC_COMM_WORLD);
	MPI_Bcast(weakening_seed,n_interfaces+1,MPIU_SCALAR,0,PETSC_COMM_WORLD);
	MPI_Bcast(cohesion_min,n_interfaces+1,MPIU_SCALAR,0,PETSC_COMM_WORLD);
	MPI_Bcast(cohesion_max,n_interfaces+1,MPIU_SCALAR,0,PETSC_COMM_WORLD);
	MPI_Bcast(friction_angle_min,n_interfaces+1,MPIU_SCALAR,0,PETSC_COMM_WORLD);
	MPI_Bcast(friction_angle_max,n_interfaces+1,MPIU_SCALAR,0,PETSC_COMM_WORLD);

	// Broadcast, special cases
//	MPI_Bcast(&n_interfaces,1,MPI_INT,0,PETSC_COMM_WORLD); // Broadcast after interfaces.txt

	// Variable [B]oundary [C]onditions for the [V]elocity
	if (variable_bcv==1){
		FILE *f_bcv;
		f_bcv = fopen("scale_bcv.txt","r");
		if (f_bcv==NULL) {
			PetscPrintf(PETSC_COMM_WORLD,"\n\n\n\nscale_bcv.txt not found\n\n\n\n");
			exit(1);
		}
		if (rank==0){
			fscanf(f_bcv,"%d",&n_var_bcv);
		}
		MPI_Bcast(&n_var_bcv,1,MPI_INT,0,PETSC_COMM_WORLD);
		PetscCalloc1(n_var_bcv,&var_bcv_scale);
		PetscCalloc1(n_var_bcv,&var_bcv_time);
		if (rank==0){
			printf("Variable boundary condition for velocity\n");
			for (int i=0;i<n_var_bcv;i++){
				fscanf(f_bcv,"%lf%lf",&var_bcv_time[i],&var_bcv_scale[i]);
				printf("%lf Myr, scale factor = %lf\n",var_bcv_time[i],var_bcv_scale[i]);
			}
			printf("\n\n");
			fclose(f_bcv);
		}
		MPI_Bcast(var_bcv_time,n_var_bcv,MPIU_SCALAR,0,PETSC_COMM_WORLD);
		MPI_Bcast(var_bcv_scale,n_var_bcv,MPIU_SCALAR,0,PETSC_COMM_WORLD);
	}

	if (multi_velocity==1){
		FILE *f_mv;
		f_mv = fopen("multi_veloc.txt","r");
		if (f_mv==NULL){
			PetscPrintf(PETSC_COMM_WORLD,"\n\n\n\nmulti_veloc.txt not found\n\n\n\n");
			exit(1);
		}
		if (rank==0){
			fscanf(f_mv,"%d",&n_mv);
		}
		MPI_Bcast(&n_mv,1,MPI_INT,0,PETSC_COMM_WORLD);
		PetscCalloc1(n_mv,&mv_time);
		if (rank==0){
			printf("Multi velocity files\n");
			for (int i=0;i<n_mv;i++){
				fscanf(f_mv,"%lf",&mv_time[i]);
				printf("%lf Myr\n",mv_time[i]);
			}
			fclose(f_mv);
			printf("\n");
		}
		MPI_Bcast(mv_time,n_mv,MPIU_SCALAR,0,PETSC_COMM_WORLD);
	}

	PetscPrintf(PETSC_COMM_WORLD,"climate_change %d\n\n",climate_change);
	if (climate_change==1){
		FILE *f_climate;
		f_climate = fopen("climate.txt","r");
		if (f_climate==NULL) {
			PetscPrintf(PETSC_COMM_WORLD,"\n\n\n\nclimate.txt not found\n\n\n\n");
			exit(1);
		}
		if (rank==0){
			fscanf(f_climate,"%d",&n_var_climate);
		}
		MPI_Bcast(&n_var_climate,1,MPI_INT,0,PETSC_COMM_WORLD);
		PetscCalloc1(n_var_climate,&var_climate_scale);
		PetscCalloc1(n_var_climate,&var_climate_time);
		if (rank==0){
			printf("Variable climate\n");
			for (int i=0;i<n_var_climate;i++){
				fscanf(f_climate,"%lf%lf",&var_climate_time[i],&var_climate_scale[i]);
				printf("%lf %lf\n",var_climate_time[i],var_climate_scale[i]);
			}
			printf("\n");
			fclose(f_climate);
		}
		MPI_Bcast(var_climate_time,n_var_climate,MPIU_SCALAR,0,PETSC_COMM_WORLD);
		MPI_Bcast(var_climate_scale,n_var_climate,MPIU_SCALAR,0,PETSC_COMM_WORLD);
	}

	PetscFunctionReturn(0);

}

/* Filename: topo_var.txt */
/* File format: */
/*   #number_of_profiles*/
/*   t0 h h h ... h (#Nx h values) */
/*   t1 h h h ... h (#Nx h values) */
PetscErrorCode load_topo_var(int rank)
{
    PetscInt i;
    PetscInt j;
    FILE *f_topo_var;


    f_topo_var = fopen("topo_var.txt", "r");
    if (f_topo_var == NULL) {
        PetscPrintf(PETSC_COMM_WORLD, "\n\n\n\ntopo_var.txt not found\n\n\n\n");
        exit(-1);
    }

    sp_n_profiles = 0;

    if (rank == 0) {
        fscanf(f_topo_var, "%ld", &sp_n_profiles);
    }

    MPI_Barrier(PETSC_COMM_WORLD);
    MPI_Bcast(&sp_n_profiles, 1, MPI_LONG, 0, PETSC_COMM_WORLD);

    if (sp_n_profiles == 0) {
        PetscPrintf(PETSC_COMM_WORLD, "\n\n\n\nerror topo_var.txt\n\n\n\n");
        exit(-1);
    }

    PetscSynchronizedPrintf(PETSC_COMM_WORLD, "rank=%d Nx=%d sp_n_profiles=%d\n", rank, Nx, sp_n_profiles);
    PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
    MPI_Barrier(PETSC_COMM_WORLD);

    MPI_Barrier(PETSC_COMM_WORLD);
    PetscCalloc1(sp_n_profiles, &topo_var_time);
    PetscCalloc1(Nx*sp_n_profiles, &topo_var_rate);

    if (rank == 0) {
        for(i = 0; i < sp_n_profiles; i++) {
            fscanf(f_topo_var, "%lf", &topo_var_time[i]);

            for(j = 0; j < Nx; j++) {
                fscanf(f_topo_var, "%lf", &topo_var_rate[j + i*Nx]);
            }
        }
    }

    fclose(f_topo_var);

    MPI_Barrier(PETSC_COMM_WORLD);
    MPI_Bcast(topo_var_time, sp_n_profiles, MPIU_SCALAR, 0, PETSC_COMM_WORLD);
    MPI_Bcast(topo_var_rate, Nx*sp_n_profiles, MPIU_SCALAR, 0, PETSC_COMM_WORLD);

    PetscFunctionReturn(0);
}

void ErrorInterfaces(){
	printf("\n\n\n\n");
	printf("Problem to read interfaces_creep.txt:\n");
	printf("Check if the number of interfaces is\n");
	printf("compatible to the one indicated in the\n");
	printf("parameters file.\n\n\n\n");
}

// Interface file error
//void ErrorInterfaces(int rank, const char fname[], int flag)
//{
//	if (rank == 0)
//	{
//		if (flag == 0)
//		{
////			print_stars();
//			fprintf(stderr, "ERROR. The number of values in the header of the <%s> file \n", fname);
//			fprintf(stderr, "is INCONSISTENT to define a number of interfaces.\n");
////			print_stars();
//		}
//		else if (flag == 1)
//		{
////			print_stars();
//			fprintf(stderr, "ERROR. The number of columns for the data in <%s> file \n", fname);
//			fprintf(stderr, "DOES NOT match the requisites for the number of interfaces.\n");
////			print_stars();
//		}
//		else if (flag == 2)
//		{
////			print_stars();
//			fprintf(stderr, "ERROR. The number of lines for the data in <%s> file \n", fname);
//			fprintf(stderr, "DOES NOT match the requisites for the number of interfaces.\n");
////			print_stars();
//		}
//	}
//}

// Function that returns 1 if tkn_v equals str_a; 0 tkn_v equals str_b; and
// exit with code (1) if none.
int check_a_b(char tkn_w[], char tkn_v[], const char str_a[], const char str_b[])
{
	int value;
	if (strcmp(tkn_v, str_a) == 0) value = 1;
	else if (strcmp(tkn_v, str_b) == 0) value = 0;
	else
	{
		fprintf(stderr, "ERROR. Unrecognized value <%s> for <%s> in the parameter file.\n", tkn_v, tkn_w);
		exit(1);
	}
	return value;
}
PetscBool check_a_b_bool(char tkn_w[], char tkn_v[], const char str_a[], const char str_b[])
{
	PetscBool value;
	if (strcmp(tkn_v, str_a) == 0) value = PETSC_TRUE;
	else if (strcmp(tkn_v, str_b) == 0) value = PETSC_FALSE;
	else
	{
		fprintf(stderr, "ERROR. Unrecognized value <%s> for <%s> in the parameter file.\n", tkn_v, tkn_w);
		exit(1);
	}
	return value;
}
