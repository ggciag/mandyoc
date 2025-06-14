static char help[] = "\n\nMANDYOC: MANtle DYnamics simulatOr Code\n\n"\
"Flags:\n\n"\
"   -seed [int]:                  specify one (or more, comma separated) layer for weak plastic criterium (seed layer)\n"\
"                                 default value: no layer specified\n\n"\
"   -strain_seed [float]:         specify one (or more, comma separated) value for the seed layer strain\n"\
"                                 default value: 0.5\n\n"\
"   -strain_seed_constant [bool]: specify if the strain_seed value is used as a constant value (default)\n"\
"                                 or if it is summed over the random initial strain\n"\
"                                 default value: true\n\n"\
"";

/* MANDYOC: MANtle DYnamics simulatOr Code*/
/* Geophysics Department IAG/USP          */
/* Victor Sacek                           */
/* Jamison F. Assuncao                    */
/* Agustina Pesce                         */
/* Rafael M. da Silva                     */
/* [Contributed by Dave May]              */

#include <petscksp.h>
#include <petscdmda.h>
#include <petsctime.h>
#include "petscsys.h"
#include "header.h"

// Petsc prototypes
PetscErrorCode create_thermal(int dimensions, PetscInt mx, PetscInt my, PetscInt mz, PetscInt Px, PetscInt Py, PetscInt Pz);
PetscErrorCode build_thermal(int dimensions);
PetscErrorCode solve_thermal(int dimensions);
PetscErrorCode destroy_thermal_();
PetscErrorCode write_all_(int cont,Vec u, char *variable_name, PetscInt binary_out);
PetscErrorCode write_pressure(int cont, PetscInt binary_out);
PetscErrorCode write_geoq_(int cont, PetscInt binary_out);
PetscErrorCode create_veloc(int dimensions, PetscInt mx, PetscInt my, PetscInt mz, PetscInt Px, PetscInt Py, PetscInt Pz);
PetscErrorCode createSwarm_2d();
PetscErrorCode createSwarm_3d();
PetscErrorCode moveSwarm(int dimensions, PetscReal dt);
PetscErrorCode Swarm_add_remove_2d();
PetscErrorCode Swarm_add_remove_3d();
PetscErrorCode SwarmViewGP_2d(DM dms,const char prefix[]);
PetscErrorCode SwarmViewGP_3d(DM dms,const char prefix[]);
PetscErrorCode Init_Veloc(int dimensions);
PetscErrorCode reader(int rank, const char fName[]);
PetscErrorCode write_veloc(int cont, PetscInt binary_out);
PetscErrorCode write_veloc_cond(int cont, PetscInt binary_out);
PetscErrorCode destroy_veloc();
double Calc_dt_calor(int rank);
PetscErrorCode write_tempo(int cont);
PetscErrorCode calc_drho();
PetscErrorCode veloc_total(int dimensions);
PetscErrorCode rescaleVeloc(Vec Veloc_fut, double tempo);
PetscErrorCode multi_veloc_change(Vec Veloc_fut,double tempo);
PetscErrorCode parse_options(int rank);

PetscErrorCode sp_create_surface_swarm_2d();
PetscErrorCode sp_move_surface_swarm(PetscInt dimensions, PetscReal dt);
PetscErrorCode sp_surface_swarm_interpolation();
PetscErrorCode sp_evaluate_surface_processes(PetscInt dimensions, PetscReal dt);
PetscErrorCode sp_update_surface_swarm_particles_properties();
PetscErrorCode sp_destroy();
PetscErrorCode sp_view_2d(DM dm, const char prefix[]);

int main(int argc,char **args)
{
	PetscErrorCode ierr;
	char prefix[PETSC_MAX_PATH_LEN];

	PetscFunctionBeginUser;

	ierr = PetscInitialize(&argc,&args,(char*)0,help);CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD,"           __  __              _   _   _____   __     __   ____     _____ \n");
	PetscPrintf(PETSC_COMM_WORLD,"          |  \\/  |     /\\     | \\ | | |  __ \\  \\ \\   / /  / __ \\   / ____|\n");
	PetscPrintf(PETSC_COMM_WORLD,"          | \\  / |    /  \\    |  \\| | | |  | |  \\ \\_/ /  | |  | | | |     \n");
	PetscPrintf(PETSC_COMM_WORLD,"          | |\\/| |   / /\\ \\   | . ` | | |  | |   \\   /   | |  | | | |     \n");
	PetscPrintf(PETSC_COMM_WORLD,"          | |  | |  / ____ \\  | |\\  | | |__| |    | |    | |__| | | |____ \n");
	PetscPrintf(PETSC_COMM_WORLD,"          |_|  |_| /_/    \\_\\ |_| \\_| |_____/     |_|     \\____/   \\_____|\n");
	PetscPrintf(PETSC_COMM_WORLD,"                                                                          \n");

	PetscPrintf(PETSC_COMM_WORLD,"===================================================================================\n");
	PetscPrintf(PETSC_COMM_WORLD,"=   MANDYOC: MANtle DYnamics simulatOr Code.\n");
	PetscPrintf(PETSC_COMM_WORLD,"===================================================================================\n");

	#ifndef GIT_VERSION
	ierr = PetscPrintf(PETSC_COMM_WORLD, "*** Git version: %s ***\n\n", GIT_VERSION);CHKERRQ(ierr);
	#endif

	PetscBool flags;
	ierr = PetscOptionsHasName(NULL,NULL,"-flags",&flags);
	if (flags){
		PetscPrintf(PETSC_COMM_WORLD,"%s",help);
		exit(1);
	}

	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

	seed = rank;

	// Parse command line options
	ierr = parse_options(rank); CHKERRQ(ierr);

	// Read ASCII files
	ierr = reader(rank, "param.txt"); CHKERRQ(ierr);

	// Check if the number of interfaces in "param.txt" is higher than then number read from command line
	if (seed_layer_set && seed_layer_size > (n_interfaces + 1)) {
		PetscPrintf(PETSC_COMM_WORLD, "Error: The number of layers specified in command line \"-seed\" command is higher than the number in \"param.txt\".\n");
	}

	// Check surface processes configuration
	if (sp_surface_processes == PETSC_TRUE && (sp_mode != 1)) {
		PetscPrintf(PETSC_COMM_WORLD, "Error: Invalid \"sp_mode\" in \"param.txt\".\n");
		exit(1);
	}

	PetscPrintf(PETSC_COMM_WORLD, "Number of seed layers: %d\n", seed_layer_size);
	for (int k = 0; k < seed_layer_size; k++) {
		PetscPrintf(PETSC_COMM_WORLD, "seed layer: %d - strain: %lf\n", seed_layer[k], strain_seed_layer[k]);
	}
	PetscPrintf(PETSC_COMM_WORLD, "\n");

	// surface processes
	if (sp_surface_processes == PETSC_TRUE && sp_mode == 1) {
		if (set_sp_d_c == PETSC_FALSE) {
			ierr = PetscPrintf(PETSC_COMM_WORLD, "'sp_mode=1' (diffusion) using default value: sp_d_c %e\n", sp_d_c); CHKERRQ(ierr);
		} else {
			ierr = PetscPrintf(PETSC_COMM_WORLD, "'sp_mode=1' (diffusion) using custom value: sp_d_c %e\n", sp_d_c); CHKERRQ(ierr);
		}
	}

	// Update elements aux constants
	if (dimensions == 3) {
		T_NE = 8;
		T_GN = 1;
		DIMEN = 3;
		GaussQuad = 27;
		V_NE = 8;
		V_GN = 3;
		V_GT = V_NE*V_GN;
	}

	PetscLogDouble Tempo1,Tempo2;
	PetscTime(&Tempo1);
	char variable_name[100];

	dx_const = Lx/(Nx-1);
	dy_const = Ly/(Ny-1);
	dz_const = depth/(Nz-1);

	//PetscPrintf(PETSC_COMM_WORLD, "dx=%lf dz=%lf\n",dx_const,dz_const);

	ierr = create_thermal(dimensions, Nx-1, Ny-1, Nz-1, Px, Py, Pz);CHKERRQ(ierr);

	sprintf(variable_name,"temperature");
	ierr = write_all_(-1,Temper, variable_name, binary_output);

	ierr = create_veloc(dimensions, Nx-1, Ny-1, Nz-1, Px, Py, Pz);CHKERRQ(ierr);

	if (geoq_on){
		PetscPrintf(PETSC_COMM_WORLD,"\nSwarm (creating)\n");
		if (dimensions == 2) {
			ierr = createSwarm_2d(); CHKERRQ(ierr);
		} else {
			ierr = createSwarm_3d(); CHKERRQ(ierr);
		}
		PetscPrintf(PETSC_COMM_WORLD,"Swarm: done\n");
	}

	if (dimensions == 2 && (sp_surface_tracking)) {
		PetscPrintf(PETSC_COMM_WORLD,"\nSurface Processes Swarm (creating)\n");

		sp_create_surface_swarm_2d();

		PetscPrintf(PETSC_COMM_WORLD,"Surface Processes Swarm: done\n");
	}

//	PetscPrintf(PETSC_COMM_SELF,"********** <rank:%d> <particles_per_ele:%d>\n", rank, particles_per_ele); -> conversar com Victor sobre
//	PetscPrintf(PETSC_COMM_SELF,"********** <rank:%d> <layers:%d>\n", rank, layers); -> conversar com Victor sobre

	if (dimensions == 3) {
			calc_drho(); // CHECK: does need this call here?
	}

	// Gerya p. 215
	if (visc_MAX>visc_MIN && initial_dynamic_range>0){
		double visc_contrast = PetscLog10Real(visc_MAX/visc_MIN);

		double visc_mean = PetscPowReal(10.0,PetscLog10Real(visc_MIN)+visc_contrast/2);

		int n_visc=0;

		visc_MIN_comp = visc_mean;
		visc_MAX_comp = visc_mean;

		PetscPrintf(PETSC_COMM_WORLD,"\n\nViscosity range: %.3lg %.3lg\n\n",visc_MIN_comp,visc_MAX_comp);

		ierr = veloc_total(dimensions); CHKERRQ(ierr);

		while ((visc_MIN_comp!=visc_MIN) && (visc_MAX_comp!=visc_MAX)){

			visc_MIN_comp = visc_mean*PetscPowReal(10.0,-n_visc*1.0);
			visc_MAX_comp = visc_mean*PetscPowReal(10.0,n_visc*1.0);

			if (visc_MIN_comp<visc_MIN) visc_MIN_comp=visc_MIN;
			if (visc_MAX_comp>visc_MAX) visc_MAX_comp=visc_MAX;

			PetscPrintf(PETSC_COMM_WORLD,"\n\nViscosity range: %.3lg %.3lg\n\n",visc_MIN_comp,visc_MAX_comp);

			ierr = veloc_total(dimensions); CHKERRQ(ierr);

			n_visc++;
		}
	}
	else {
		if (visc_MAX==visc_MIN) visc_MAX = visc_MIN*1.0001;  //avoiding the problem to the f2 in the denominator (Gerya...)
		visc_MIN_comp = visc_MIN;
		visc_MAX_comp = visc_MAX;
		ierr = veloc_total(dimensions); CHKERRQ(ierr);
	}

	PetscPrintf(PETSC_COMM_WORLD,"Solution of the pressure and velocity fields: done\n");

	PetscPrintf(PETSC_COMM_WORLD,"\nWriting output files:\n");
	ierr = write_veloc(tcont,binary_output);
	ierr = write_veloc_cond(tcont,binary_output);

	sprintf(variable_name,"temperature");
	ierr = write_all_(tcont,Temper,variable_name, binary_output);
	ierr = write_pressure(tcont,binary_output);
	ierr = write_geoq_(tcont,binary_output);
	ierr = write_tempo(tcont);

	if (dimensions == 2 && sp_surface_tracking) {
		ierr = PetscSNPrintf(prefix, PETSC_MAX_PATH_LEN-1,"surface_%d", tcont); CHKERRQ(ierr);
		ierr = sp_view_2d(dms_s, prefix); CHKERRQ(ierr);
	}

	VecCopy(Veloc_fut,Veloc);
	dt_calor_sec = Calc_dt_calor(rank);

	if (initial_print_step > 0) {
		print_step_aux = print_step;
		print_step = initial_print_step;

		PetscPrintf(PETSC_COMM_WORLD, "\n\n*** Using custom initial print_step = %d until %.3g Myr\n\n", print_step, initial_print_max_time);
	}

	for (tempo = dt_calor,tcont=1;tempo<=timeMAX && tcont<=stepMAX;tempo+=dt_calor, tcont++){
		if ((tempo > initial_print_max_time || fabs(tempo-initial_print_max_time) < 0.0001) && initial_print_step > 0) {
			initial_print_step = 0;
			print_step = print_step_aux;
			PetscPrintf(PETSC_COMM_WORLD, "\n\n*** Restored print_step = %d\n\n", print_step);
		}
		if (non_dim==1){
			PetscPrintf(PETSC_COMM_WORLD,"\n\nstep = %d, time = %.3g, dt = %.3g\n",tcont,tempo,dt_calor);
		}
		else{
			PetscPrintf(PETSC_COMM_WORLD,"\n\nstep = %d, time = %.3g yr, dt = %.3g yr\n",tcont,tempo,dt_calor);
		}


		//PetscPrintf(PETSC_COMM_WORLD,"next sp %.3g Myr\n\n", sp_eval_time);

		if (dimensions == 2) {
			ierr = rescaleVeloc(Veloc_fut,tempo);
			ierr = multi_veloc_change(Veloc_fut,tempo);
		}

		ierr = build_thermal(dimensions);CHKERRQ(ierr);

		ierr = solve_thermal(dimensions);CHKERRQ(ierr);

		ierr = veloc_total(dimensions); CHKERRQ(ierr);

		if (geoq_on){
			if (RK4==1){
				VecCopy(Veloc_fut,Veloc_weight);
				ierr = moveSwarm(dimensions, dt_calor_sec);
			}
			else {
				for (PetscInt cont=0, max_cont=4;cont<max_cont; cont++){
					double fac = (1.0/max_cont)*(0.5+cont);
					//PetscPrintf(PETSC_COMM_WORLD,"%f %f\n",fac,(1.0-fac));
					//        x   ->   y
					VecCopy(Veloc,Veloc_weight);
					//			y			a		b		x
					VecAXPBY(Veloc_weight, fac, (1.0-fac),Veloc_fut); //y = a*x + b*y
					ierr = moveSwarm(dimensions, dt_calor_sec/max_cont);

					if (dimensions == 2 && sp_surface_tracking) {
						ierr = sp_move_surface_swarm(dimensions, dt_calor_sec/max_cont); CHKERRQ(ierr);
					}
				}
			}

			if (dimensions == 2) {
				Swarm_add_remove_2d();
			} else {
				Swarm_add_remove_3d();
			}
		}

		if (dimensions == 2 && sp_surface_tracking) {
			ierr = sp_surface_swarm_interpolation(); CHKERRQ(ierr);

			if (sp_surface_processes) {
				ierr = sp_evaluate_surface_processes(dimensions, dt_calor_sec); CHKERRQ(ierr);
				ierr = sp_update_surface_swarm_particles_properties(); CHKERRQ(ierr);
			}
		}

		if (tcont%print_step==0){
			PetscPrintf(PETSC_COMM_WORLD,"\nWriting output files:\n");
			sprintf(variable_name,"temperature");
			ierr = write_all_(tcont,Temper,variable_name,binary_output);
			ierr = write_geoq_(tcont,binary_output);
			ierr = write_veloc(tcont,binary_output);
			ierr = write_pressure(tcont,binary_output);
			ierr = write_tempo(tcont);
			PetscSNPrintf(prefix,PETSC_MAX_PATH_LEN-1,"step_%d",tcont);
			if (geoq_on){
				if (print_step_files==1){
					if (dimensions == 2) {
						ierr = SwarmViewGP_2d(dms,prefix);CHKERRQ(ierr);
					} else {
						ierr = SwarmViewGP_3d(dms,prefix);CHKERRQ(ierr);
					}
				}

				if (dimensions == 2 && sp_surface_tracking) {
					ierr = PetscSNPrintf(prefix, PETSC_MAX_PATH_LEN-1,"surface_%d", tcont); CHKERRQ(ierr);
					ierr = sp_view_2d(dms_s, prefix); CHKERRQ(ierr);
				}
			}
		}

		dt_calor_sec = Calc_dt_calor(rank);

	}
	PetscPrintf(PETSC_COMM_WORLD, "write\n");



	ierr = destroy_thermal_();CHKERRQ(ierr);

	destroy_veloc();

	if (dimensions == 2 && sp_surface_tracking) {
		sp_destroy();
	}

	PetscTime(&Tempo2);

	ierr = PetscFinalize();

	PetscFunctionReturn(0);
}


double Calc_dt_calor(int rank) {
	//////calc dt
	PetscInt ind_v_max,ind_v_min,ind_v_mod;
	PetscReal min_v,max_v,max_mod_v,dh_v_mod;

	VecMax(Veloc_fut,&ind_v_min,&max_v);
	VecMin(Veloc_fut,&ind_v_max,&min_v);
	max_mod_v = fabs(max_v);
	ind_v_mod = ind_v_max;
	if (max_mod_v<fabs(min_v)){
		max_mod_v = fabs(min_v);
		ind_v_mod = ind_v_min;
	}

	if (dimensions == 2) {
		dh_v_mod = dx_const; //check: only in 2d
		if (ind_v_mod%2==1) dh_v_mod = dz_const; //check: only in 2d
		dt_calor = 0.1*(dh_v_mod/max_mod_v)/(seg_per_ano*sub_division_time_step);
	} else {
		if (ind_v_mod%3==0) dh_v_mod = dx_const;
		if (ind_v_mod%3==1) dh_v_mod = dy_const;
		if (ind_v_mod%3==2) dh_v_mod = dz_const;
		//PetscPrintf(PETSC_COMM_WORLD, "dt = %g",(dh_v_mod/max_mod_v)/seg_per_ano);
		dt_calor = 0.1*(dh_v_mod/max_mod_v)/(seg_per_ano*sub_division_time_step);
	}

	if (dt_calor>dt_MAX) dt_calor=dt_MAX;
	dt_calor_sec = dt_calor*seg_per_ano;

	if (non_dim==1) dt_calor*=seg_per_ano;

	return dt_calor_sec;
}


PetscErrorCode write_tempo(int cont){

	int rank;

	PetscFunctionBeginUser;

	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

	PetscViewer viewer;

	char nome[100];

	sprintf(nome,"time_%d.txt",cont);

	PetscReal aa[1];

	aa[0]=tempo;

	PetscViewerASCIIOpen(PETSC_COMM_WORLD,nome,&viewer);
	PetscRealView(1,aa,viewer);
	PetscViewerDestroy(&viewer);




	PetscFunctionReturn(0);

}


PetscErrorCode rescaleVeloc(Vec Veloc_fut,double tempo)
{
	PetscErrorCode ierr;

	PetscFunctionBeginUser;

	if (cont_var_bcv<n_var_bcv){
		if (tempo>1.0E6*var_bcv_time[cont_var_bcv]){
			PetscScalar v_norm;
			VecNorm(Veloc_fut,NORM_1,&v_norm);
			PetscPrintf(PETSC_COMM_WORLD,"v_norm = %lf\n\n",v_norm);
			VecScale(Veloc_fut,var_bcv_scale[cont_var_bcv]);
			VecScale(Veloc_0,var_bcv_scale[cont_var_bcv]);
			cont_var_bcv++;
			VecNorm(Veloc_fut,NORM_1,&v_norm);
			PetscPrintf(PETSC_COMM_WORLD,"v_norm = %lf\n\n",v_norm);

			PetscPrintf(PETSC_COMM_WORLD,"velocity rescale\n\n");

			if (visc_MAX>visc_MIN && initial_dynamic_range>0){
				double visc_contrast = PetscLog10Real(visc_MAX/visc_MIN);

				double visc_mean = PetscPowReal(10.0,PetscLog10Real(visc_MIN)+visc_contrast/2);

				int n_visc=0;

				visc_MIN_comp = visc_mean;
				visc_MAX_comp = visc_mean;

				PetscPrintf(PETSC_COMM_WORLD,"\n\nViscosity range: %.3lg %.3lg\n\n",visc_MIN_comp,visc_MAX_comp);

				ierr = veloc_total(dimensions); CHKERRQ(ierr);

				while ((visc_MIN_comp!=visc_MIN) && (visc_MAX_comp!=visc_MAX)){

					visc_MIN_comp = visc_mean*PetscPowReal(10.0,-n_visc*1.0);
					visc_MAX_comp = visc_mean*PetscPowReal(10.0,n_visc*1.0);

					if (visc_MIN_comp<visc_MIN) visc_MIN_comp=visc_MIN;
					if (visc_MAX_comp>visc_MAX) visc_MAX_comp=visc_MAX;

					PetscPrintf(PETSC_COMM_WORLD,"\n\nViscosity range: %.3lg %.3lg\n\n",visc_MIN_comp,visc_MAX_comp);

					ierr = veloc_total(dimensions); CHKERRQ(ierr);

					n_visc++;
				}
			}
		}
	}

	PetscFunctionReturn(0);
}

PetscErrorCode multi_veloc_change(Vec Veloc_fut,double tempo){
	PetscErrorCode ierr;

	PetscFunctionBeginUser;

	if (cont_mv<n_mv){
		if (tempo>1.0E6*mv_time[cont_mv]){
			cont_mv++;
			Init_Veloc(dimensions);
			if (visc_MAX>visc_MIN && initial_dynamic_range>0){
				double visc_contrast = PetscLog10Real(visc_MAX/visc_MIN);

				double visc_mean = PetscPowReal(10.0,PetscLog10Real(visc_MIN)+visc_contrast/2);

				int n_visc=0;

				visc_MIN_comp = visc_mean;
				visc_MAX_comp = visc_mean;

				PetscPrintf(PETSC_COMM_WORLD,"\n\nViscosity range: %.3lg %.3lg\n\n",visc_MIN_comp,visc_MAX_comp);

				ierr = veloc_total(dimensions); CHKERRQ(ierr);

				while ((visc_MIN_comp!=visc_MIN) && (visc_MAX_comp!=visc_MAX)){

					visc_MIN_comp = visc_mean*PetscPowReal(10.0,-n_visc*1.0);
					visc_MAX_comp = visc_mean*PetscPowReal(10.0,n_visc*1.0);

					if (visc_MIN_comp<visc_MIN) visc_MIN_comp=visc_MIN;
					if (visc_MAX_comp>visc_MAX) visc_MAX_comp=visc_MAX;

					PetscPrintf(PETSC_COMM_WORLD,"\n\nViscosity range: %.3lg %.3lg\n\n",visc_MIN_comp,visc_MAX_comp);

					ierr = veloc_total(dimensions); CHKERRQ(ierr);

					n_visc++;
				}
			}
		}
	}
	PetscFunctionReturn(0);
}
