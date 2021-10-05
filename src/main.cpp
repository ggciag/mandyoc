static char help[] = "\n\nMANDYOC: MANtle DYnamics simulatOr Code\n\n"\
"Flags:\n\n"\
"   -seed [int]:          specify one (or more, comma separated) layer for weak plastic criterium (seed layer)\n"\
"                         default value: no layer specified\n\n"\
"   -strain_seed [float]: specify one (or more, comma separated) value for the seed layer strain\n"\
"                         default value: 2.0\n\n"\
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
PetscErrorCode create_thermal_2d(PetscInt mx,PetscInt mz,PetscInt Px,PetscInt Pz);
PetscErrorCode build_thermal_3d();
PetscErrorCode solve_thermal_3d();
PetscErrorCode destroy_thermal_();
PetscErrorCode write_all_(int cont,Vec u, char *variable_name, PetscInt binary_out);
PetscErrorCode write_pressure(int cont, PetscInt binary_out);
PetscErrorCode write_geoq_(int cont, PetscInt binary_out);
PetscErrorCode create_veloc_2d(PetscInt mx,PetscInt mz,PetscInt Px,PetscInt Pz);
PetscErrorCode createSwarm();
PetscErrorCode moveSwarm(PetscReal dt);
PetscErrorCode Swarm_add_remove();
PetscErrorCode SwarmViewGP(DM dms,const char prefix[]);
PetscErrorCode Init_Veloc();
PetscErrorCode reader(int rank, const char fName[]);
PetscErrorCode write_veloc_3d(int cont, PetscInt binary_out);
PetscErrorCode write_veloc_cond(int cont, PetscInt binary_out);
PetscErrorCode destroy_veloc_3d();
PetscErrorCode Calc_dt_calor();
PetscErrorCode write_tempo(int cont);
PetscErrorCode veloc_total();
PetscErrorCode rescaleVeloc(Vec Veloc_fut, double tempo);
PetscErrorCode multi_veloc_change(Vec Veloc_fut,double tempo);
PetscErrorCode sp_create_surface_vec();
PetscErrorCode sp_interpolate_surface_particles_to_vec();
PetscErrorCode evaluate_surface_processes();
PetscErrorCode sp_write_surface_vec(PetscInt i);
PetscErrorCode sp_destroy();
PetscErrorCode load_topo_var(int rank);
PetscErrorCode rescalePrecipitation(double tempo);

int main(int argc,char **args)
{
	PetscErrorCode ierr;
	char prefix[PETSC_MAX_PATH_LEN];
	PetscInt Px,Pz;

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
	ierr = PetscOptionsGetInt(NULL,NULL,"-dimen",&DIMEN,NULL);CHKERRQ(ierr);
	if (DIMEN!=2 && DIMEN!=3){
		PetscPrintf(PETSC_COMM_WORLD,"Incorrect dimension.\n");
		exit(-3);
	}

	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

	seed = rank;
	
	// Read ASCII files
	reader(rank, "param.txt");

	PetscLogDouble Tempo1,Tempo2;
	PetscTime(&Tempo1);
	char variable_name[100];

	Px = Pz = PETSC_DECIDE;
	ierr = PetscOptionsGetInt(NULL,NULL,"-Px",&Px,NULL);CHKERRQ(ierr);
	Pz = Px;
	ierr = PetscOptionsGetInt(NULL,NULL,"-Pz",&Pz,NULL);CHKERRQ(ierr);

	if (n_interfaces>0 && interfaces_from_ascii==1){
		ierr = PetscCalloc1(n_interfaces, &seed_layer); CHKERRQ(ierr);
		seed_layer_size = n_interfaces;
		ierr = PetscOptionsGetIntArray(NULL,NULL,"-seed",seed_layer,&seed_layer_size,&seed_layer_set); CHKERRQ(ierr);

		ierr = PetscCalloc1(n_interfaces, &strain_seed_layer); CHKERRQ(ierr);
		strain_seed_layer_size = n_interfaces;
		ierr = PetscOptionsGetRealArray(NULL,NULL,"-strain_seed",strain_seed_layer,&strain_seed_layer_size,&strain_seed_layer_set); CHKERRQ(ierr);
		if (strain_seed_layer_set == PETSC_TRUE && seed_layer_set == PETSC_FALSE) {
			PetscPrintf(PETSC_COMM_WORLD,"Specify the seed layer with the flag -seed (required by -strain_seed)\n");
			exit(1);
		}
		if (strain_seed_layer_set == PETSC_TRUE && seed_layer_set == PETSC_TRUE && seed_layer_size != strain_seed_layer_size) {
			PetscPrintf(PETSC_COMM_WORLD,"Specify the same number of values in the list for flags -seed and -strain_seed\n");
			exit(1);
		}
		if (strain_seed_layer_set == PETSC_FALSE && seed_layer_set == PETSC_TRUE) {
			PetscPrintf(PETSC_COMM_WORLD,"Using default value '2.0' for -strain_seed (for all seed layers)\n");
			for (int k = 0; k < seed_layer_size; k++) {
				strain_seed_layer[k] = 2.0;
			}
		}
		PetscPrintf(PETSC_COMM_WORLD,"Number of seed layers: %d\n", seed_layer_size);
		for (int k = 0; k < seed_layer_size; k++) {
			PetscPrintf(PETSC_COMM_WORLD,"seed layer: %d - strain: %lf\n", seed_layer[k], strain_seed_layer[k]);
		}
		PetscPrintf(PETSC_COMM_WORLD,"\n");
	}

	h_air=-1.0;
	ierr = PetscOptionsGetReal(NULL,NULL,"-h_air",&h_air,NULL);CHKERRQ(ierr);
	if (pressure_in_rheol==0 && h_air<0.0){
		PetscPrintf(PETSC_COMM_WORLD,"Specify the thickness of the air layer with the flag -h_air\n");
		PetscPrintf(PETSC_COMM_WORLD,"(you adopted depth dependent rheology: -pressure_in_rheol 0)\n");
		exit(1);
	}
	else h_air=0.0;


	if (sp_surface_processes && sp_surface_tracking && sp_mode == 1) load_topo_var(rank);

	
	if (sp_mode == 2 && PETSC_FALSE == set_sp_d_c) {
		ierr = PetscPrintf(PETSC_COMM_WORLD,"-sp_mode 2 (diffusion) using default value: sp_d_c %e\n", sp_d_c); CHKERRQ(ierr);
	} else if (sp_mode == 2) {
		ierr = PetscPrintf(PETSC_COMM_WORLD,"-sp_mode 2 (diffusion) using custom value: sp_d_c %e\n", sp_d_c); CHKERRQ(ierr);
	} else if (sp_mode == 3) {
		ierr = PetscPrintf(PETSC_COMM_WORLD,"-sp_mode 3 (fluvial erosion) using K_fluvial: %e and sea_level %e\n", K_fluvial,sea_level); CHKERRQ(ierr);
	}else if (sp_mode == 4) {
		ierr = PetscPrintf(PETSC_COMM_WORLD,"-sp_mode 4 (fluvial erosion mode 2) using K_fluvial: %e and sea_level %e\n", K_fluvial,sea_level); CHKERRQ(ierr);
	}
	
	dx_const = Lx/(Nx-1);
	dz_const = depth/(Nz-1);

	//if (rank==0) printf("dx=%lf dz=%lf\n",dx_const,dz_const);

	ierr = create_thermal_2d(Nx-1,Nz-1,Px,Pz);CHKERRQ(ierr);

	sprintf(variable_name,"temperature");
	ierr = write_all_(-1,Temper, variable_name, binary_output);

	ierr = create_veloc_2d(Nx-1,Nz-1,Px,Pz);CHKERRQ(ierr);

	if (geoq_on){
		PetscPrintf(PETSC_COMM_WORLD,"\nSwarm (creating)\n");
		ierr = createSwarm();CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD,"Swarm: done\n");
	}
	
//	PetscPrintf(PETSC_COMM_SELF,"********** <rank:%d> <particles_per_ele:%d>\n", rank, particles_per_ele); -> conversar com Victor sobre
//	PetscPrintf(PETSC_COMM_SELF,"********** <rank:%d> <layers:%d>\n", rank, layers); -> conversar com Victor sobre
	
	// Surface Processes Swarm
	if (geoq_on && sp_surface_tracking && n_interfaces>0 && interfaces_from_ascii==1) {
		PetscPrintf(PETSC_COMM_WORLD, "\nSP Swarm (creating)\n");
		ierr = sp_create_surface_vec(); CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD, "SP Swarm: done\n");
		ierr = sp_interpolate_surface_particles_to_vec(); CHKERRQ(ierr);
	}
	
	// Gerya p. 215
	if (visc_MAX>visc_MIN && initial_dynamic_range>0){
		double visc_contrast = PetscLog10Real(visc_MAX/visc_MIN);

		double visc_mean = PetscPowReal(10.0,PetscLog10Real(visc_MIN)+visc_contrast/2);

		int n_visc=0;

		visc_MIN_comp = visc_mean;
		visc_MAX_comp = visc_mean;

		PetscPrintf(PETSC_COMM_WORLD,"\n\nViscosity range: %.3lg %.3lg\n\n",visc_MIN_comp,visc_MAX_comp);

		ierr = veloc_total(); CHKERRQ(ierr);

		while ((visc_MIN_comp!=visc_MIN) && (visc_MAX_comp!=visc_MAX)){

			visc_MIN_comp = visc_mean*PetscPowReal(10.0,-n_visc*1.0);
			visc_MAX_comp = visc_mean*PetscPowReal(10.0,n_visc*1.0);

			if (visc_MIN_comp<visc_MIN) visc_MIN_comp=visc_MIN;
			if (visc_MAX_comp>visc_MAX) visc_MAX_comp=visc_MAX;

			PetscPrintf(PETSC_COMM_WORLD,"\n\nViscosity range: %.3lg %.3lg\n\n",visc_MIN_comp,visc_MAX_comp);

			ierr = veloc_total(); CHKERRQ(ierr);

			n_visc++;
		}
	}
	else {
		if (visc_MAX==visc_MIN) visc_MAX = visc_MIN*1.0001;  //avoiding the problem to the f2 in the denominator (Gerya...)
		visc_MIN_comp = visc_MIN;
		visc_MAX_comp = visc_MAX;
		ierr = veloc_total(); CHKERRQ(ierr);
	}

	PetscPrintf(PETSC_COMM_WORLD,"Solution of the pressure and velocity fields: done\n");

	ierr = write_veloc_3d(tcont,binary_output);
	ierr = write_veloc_cond(tcont,binary_output);
	
	sprintf(variable_name,"temperature");
	ierr = write_all_(tcont,Temper,variable_name, binary_output);
	ierr = write_pressure(tcont,binary_output);
	ierr = write_geoq_(tcont,binary_output);
	ierr = write_tempo(tcont);

	if (sp_surface_tracking && geoq_on && n_interfaces>0 && interfaces_from_ascii==1) {
		sp_write_surface_vec(tcont);
	}

	VecCopy(Veloc_fut,Veloc);
	ierr = Calc_dt_calor();

	if (initial_print_step > 0) {
		print_step_aux = print_step;
		print_step = initial_print_step;

		PetscPrintf(PETSC_COMM_WORLD, "\n\n*** Using custom initial print_step = %d until %.3g Myr\n\n", print_step, initial_print_max_time);
	}

	if (sp_surface_processes && PETSC_FALSE == set_sp_dt) {
		sp_dt = 10.0 * dt_calor;
		PetscPrintf(PETSC_COMM_WORLD, "\n\n*** Using default sp_dt\n\n");
	}

	sp_eval_time = sp_dt;
	sp_last_eval_time = 0.0;

	for (tempo = dt_calor,tcont=1;tempo<=timeMAX && tcont<=stepMAX;tempo+=dt_calor, tcont++){
		if ((tempo > initial_print_max_time || fabs(tempo-initial_print_max_time) < 0.0001) && initial_print_step > 0) {
			initial_print_step = 0;
			print_step = print_step_aux;
			PetscPrintf(PETSC_COMM_WORLD, "\n\n*** Restored print_step = %d\n\n", print_step);
		}

		PetscPrintf(PETSC_COMM_WORLD,"\n\nstep = %d, time = %.3g yr, dt = %.3g yr\n",tcont,tempo,dt_calor);

		//PetscPrintf(PETSC_COMM_WORLD,"next sp %.3g Myr\n\n", sp_eval_time);

		ierr = rescaleVeloc(Veloc_fut,tempo);
		ierr = multi_veloc_change(Veloc_fut,tempo);

		ierr = build_thermal_3d();CHKERRQ(ierr);

		ierr = solve_thermal_3d();CHKERRQ(ierr);

		ierr = veloc_total(); CHKERRQ(ierr);

		if (sp_surface_processes && geoq_on && n_interfaces>0 && interfaces_from_ascii==1 && (tempo > sp_eval_time || fabs(tempo-sp_eval_time) < 0.0001)) {
			PetscPrintf(PETSC_COMM_WORLD,"\nEvaluating sp...\n");

			ierr = rescalePrecipitation(tempo);

			evaluate_surface_processes();

			sp_eval_time += sp_dt;
			sp_last_eval_time = tempo;
		}

		if (geoq_on){
			if (RK4==1){
				VecCopy(Veloc_fut,Veloc_weight);
				ierr = moveSwarm(dt_calor_sec);
			}
			else {
				for (PetscInt cont=0, max_cont=4;cont<max_cont; cont++){
					double fac = (1.0/max_cont)*(0.5+cont);
					//PetscPrintf(PETSC_COMM_WORLD,"%f %f\n",fac,(1.0-fac));
					//        x   ->   y
					VecCopy(Veloc,Veloc_weight);
					//			y			a		b		x
					VecAXPBY(Veloc_weight, fac, (1.0-fac),Veloc_fut); //y = a*x + b*y
					ierr = moveSwarm(dt_calor_sec/max_cont);
				}
			}
			Swarm_add_remove();
		}

		if (sp_surface_tracking && geoq_on && n_interfaces>0 && interfaces_from_ascii==1) {
			ierr = sp_interpolate_surface_particles_to_vec(); CHKERRQ(ierr);
		}

		if (tcont%print_step==0){
			sprintf(variable_name,"temperature");
			ierr = write_all_(tcont,Temper,variable_name,binary_output);
			ierr = write_geoq_(tcont,binary_output);
			ierr = write_veloc_3d(tcont,binary_output);
			ierr = write_pressure(tcont,binary_output);
			ierr = write_tempo(tcont);
			PetscSNPrintf(prefix,PETSC_MAX_PATH_LEN-1,"step_%d",tcont);
			if (geoq_on){
				if (print_step_files==1){
					ierr = SwarmViewGP(dms,prefix);CHKERRQ(ierr);
				}

				if (sp_surface_tracking && n_interfaces>0 && interfaces_from_ascii==1) {
					ierr = sp_write_surface_vec(tcont); CHKERRQ(ierr);
				}
			}
		}

		ierr = Calc_dt_calor();

	}
	if (rank==0) printf("write\n");



	ierr = destroy_thermal_();CHKERRQ(ierr);

	destroy_veloc_3d();

	if (geoq_on && n_interfaces>0 && interfaces_from_ascii==1)	sp_destroy();

	PetscTime(&Tempo2);


	ierr = PetscFinalize();
	return 0;
}


PetscErrorCode Calc_dt_calor(){

	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

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
	dh_v_mod = dx_const; //check: only in 2d
	if (ind_v_mod%2==1) dh_v_mod = dz_const; //check: only in 2d
	//if (rank==0) printf("dt = %g",(dh_v_mod/max_mod_v)/seg_per_ano);
	dt_calor = 0.1*(dh_v_mod/max_mod_v)/(seg_per_ano*sub_division_time_step);
	if (dt_calor>dt_MAX) dt_calor=dt_MAX;

	dt_calor_sec = dt_calor*seg_per_ano;
	

	PetscFunctionReturn(0);

}


PetscErrorCode write_tempo(int cont){

	int rank;
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

				ierr = veloc_total(); CHKERRQ(ierr);

				while ((visc_MIN_comp!=visc_MIN) && (visc_MAX_comp!=visc_MAX)){

					visc_MIN_comp = visc_mean*PetscPowReal(10.0,-n_visc*1.0);
					visc_MAX_comp = visc_mean*PetscPowReal(10.0,n_visc*1.0);

					if (visc_MIN_comp<visc_MIN) visc_MIN_comp=visc_MIN;
					if (visc_MAX_comp>visc_MAX) visc_MAX_comp=visc_MAX;

					PetscPrintf(PETSC_COMM_WORLD,"\n\nViscosity range: %.3lg %.3lg\n\n",visc_MIN_comp,visc_MAX_comp);

					ierr = veloc_total(); CHKERRQ(ierr);

					n_visc++;
				}
			}
		}
	}

	PetscFunctionReturn(0);
}

PetscErrorCode multi_veloc_change(Vec Veloc_fut,double tempo){
	PetscErrorCode ierr;
	if (cont_mv<n_mv){
		if (tempo>1.0E6*mv_time[cont_mv]){
			cont_mv++;
			Init_Veloc();
			if (visc_MAX>visc_MIN && initial_dynamic_range>0){
				double visc_contrast = PetscLog10Real(visc_MAX/visc_MIN);

				double visc_mean = PetscPowReal(10.0,PetscLog10Real(visc_MIN)+visc_contrast/2);

				int n_visc=0;

				visc_MIN_comp = visc_mean;
				visc_MAX_comp = visc_mean;

				PetscPrintf(PETSC_COMM_WORLD,"\n\nViscosity range: %.3lg %.3lg\n\n",visc_MIN_comp,visc_MAX_comp);

				ierr = veloc_total(); CHKERRQ(ierr);

				while ((visc_MIN_comp!=visc_MIN) && (visc_MAX_comp!=visc_MAX)){

					visc_MIN_comp = visc_mean*PetscPowReal(10.0,-n_visc*1.0);
					visc_MAX_comp = visc_mean*PetscPowReal(10.0,n_visc*1.0);

					if (visc_MIN_comp<visc_MIN) visc_MIN_comp=visc_MIN;
					if (visc_MAX_comp>visc_MAX) visc_MAX_comp=visc_MAX;

					PetscPrintf(PETSC_COMM_WORLD,"\n\nViscosity range: %.3lg %.3lg\n\n",visc_MIN_comp,visc_MAX_comp);

					ierr = veloc_total(); CHKERRQ(ierr);

					n_visc++;
				}
			}
		}
	}
	PetscFunctionReturn(0);
}


PetscErrorCode rescalePrecipitation(double tempo)
{
	//PetscErrorCode ierr;
	if (cont_var_climate<n_var_climate){
		if (tempo>1.0E6*var_climate_time[cont_var_climate]){
			prec_factor=var_climate_scale[cont_var_climate];
			cont_var_climate++;
		}
	}

	PetscFunctionReturn(0);
}
