static char help[] = "\n\nMANDYOC: MANtle DYnamics simulatOr Code\n\n"\
"Flags:\n\n"\
"   -Px [int]:            specify the number of cores in the x-direction of the domain\n"\
"                         default value: PETSC_DECIDE\n\n"\
"   -Pz [int]:            specify the number of cores in the z-direction of the domain\n"\
"                         default value: PETSC_DECIDE\n\n"\
"                         (the product of Px*Pz must be equal to the total number of cores)\n\n"
"   -rtol [float]:        the absolute size of the residual norm (relevant only for iterative methods)\n"\
"                         default value: PETSC_DEFAULT = 1e-5\n\n"\
"   -te [0 or 1]:         initial temperature field from an external file (1) or internally calculated (0)\n"\
"                         default value: 0\n\n"\
"   -ve [0 or 1]:         initial velocity field from an external file (1) or null (0)\n"\
"                         default value: 0\n\n"\
"   -bcve [0 or 1]:       boundary condition for the velocity field from an external file (1) or defined from the parameter file (0)\n"\
"                         default value: 0\n\n"\
"   -denok [float]:       residual norm in the Uzawa scheme\n"\
"                         default value: 1e-4\n\n"\
"   -print_visc [0 or 1]: print (1) or not (0) the viscosity field\n"\
"                         default value: 0\n\n"\
"   -visc_const_per_element [0 or 1]:\n"\
"                         assume the viscosity constant (1) of linearly variable (0) in each finite element\n"\
"                         default value: 0\n\n"\
"   -visc_harmonic_mean [0 or 1]:\n"\
"                         assume the harmonic mean (1) or arithmetic mean (0) for the effective viscosity\n"\
"                         default value: 1\n\n"\
"   -particles_per_ele [int]:\n"\
"                         number of Lagrangian particles per element\n"\
"                         default value: 81\n\n"\
"   -particles_perturb_factor [float: 0.0-1.0]:\n"\
"                         amplitude of the perturbation of the initial position of the Lagrangian particles.\n"\
"                         0.0 means no perturbation (particles aligned in a grid).\n"\
"                         default value: 0.5\n\n"\
"   -free_surface_stab [0 or 1]:\n"\
"                         use (1) or not (0) of the Free Surface Stabilization Algorithm (Kaus et al., 2010)\n"\
"                         default value: 1\n\n"\
"   -theta_FSSA [float: 0.0-1.0]:\n"\
"                         weight factor for the Free Surface Stabilization Algorithm (Kaus et al., 2010)\n"\
"                         only relevant if -free_surface_stab 1\n"\
"                         default value: 0.5\n\n"\
"   -sub_division_time_step [float]:\n"\
"                         factor that rescale the time step: new_time_step = time_step/sub_division_time_step\n"\
"                         default value: 1.0\n\n"\
"   -print_step_files [0 or 1]:\n"\
"                         print (1) or not (0) the 'step files' containing the position of the Lagrangian particles\n"\
"                         default value: 1\n\n"\
"   -checkered [0 or 1]:  in the 'step files', save only one particle per element (0) or\n"\
"                         a sequence of particles in the horizontal direction and in the vertical one per element (checkered) (1)\n"\
"                         default value: 0\n\n"\
"   -direct_solver [0 or 1]:\n"\
"                         use direct (1) or iterative (0) solver for the momentum equation\n"\
"                         default value: 1\n\n"\
"   -RK4 [0 or 1]:        UNDER CONSTRUCTION (ONLY EULER IS WORKING)\n"\
"                         use Euler (0) or 4th order Runge-Kutta (1) scheme for particles advection\n"\
"                         default value: 0\n\n"\
"   -xi_min [float]:      convergence criterium for non-linear problems\n"\
"                         default value: 1e-14\n\n"\
"   -seed [int]:          specify one layer for weak plastic criterium (seed layer)\n"\
"                         default value: no layer specified\n\n"\
"   -pressure_in_rheol [0 or 1]:\n"\
"                         (0) if the effective viscosity is depth dependent\n"\
"                         (1) if the effective viscosity is pressure dependent\n"\
"                         this flag is mandatory even is the adopted rheology is not depth/pressure dependent!\n"\
"                         no default value\n\n"\
"   -h_air [float]:       if pressure_in_rheol 1, the user must specify the thickness of the sticky air layer.\n"\
"                         only relevant if pressure_in_rheol 1 and is mandatory in this case\n"\
"                         no default value\n\n"\
"   -initial_dynamic_range [0 or 1]:\n"\
"                         adopt (1) or not (0) a progressive increase in the viscosity range, solving the momentum equation\n"\
"                         in each range until the viscosity range is from visc_MIN to visc_MAX. Only used in the first step.\n"\
"                         Useful to make a smooth convergence of the solution in the case of large viscosity contrasts\n"\
"                         (for details, see 'Introduction to Numerical Geodynamic Modelling' 1st edition, Taras Gerya, p. 215)\n"\
"                         default value: 1\n\n"\
"   -periodic_boundary [0 or 1]:\n"\
"                         if (1) assume periodic boundary condition in the x-direction\n"\
"                         default value: 0\n\n"\
"   -variable_bcv [0 or 1]:\n"\
"                         if (1) use the file scale_bcv.txt to rescale the velocity field.\n"\
"                         The scale_bcv.txt file is composed of two columns, with the first line indicating the number of instants:\n"\
"                         [number of instants]\n"\
"                         [instant to change the velocity in Myr] [scale factor]\n"\
"                         default value: 0\n\n"\
"";




/* Contributed by Dave May */

#include <petscksp.h>
#include <petscdmda.h>

#include <petsctime.h>

#include "petscsys.h"

#include "header.h"


PetscErrorCode create_thermal_2d(PetscInt mx,PetscInt mz,PetscInt Px,PetscInt Pz);

PetscErrorCode build_thermal_3d();

PetscErrorCode solve_thermal_3d();

PetscErrorCode destroy_thermal_();

PetscErrorCode write_thermal_(int cont);

PetscErrorCode write_pressure(int cont);

PetscErrorCode write_geoq_(int cont);

PetscErrorCode create_veloc_3d(PetscInt mx,PetscInt mz,PetscInt Px,PetscInt Pz);

PetscErrorCode createSwarm();
PetscErrorCode moveSwarm(PetscReal dt);
PetscErrorCode Swarm_add_remove();
PetscErrorCode SwarmViewGP(DM dms,const char prefix[]);



PetscErrorCode reader(int rank);

PetscErrorCode write_veloc_3d(int cont);
PetscErrorCode write_veloc_cond(int cont);

PetscErrorCode destroy_veloc_3d();

PetscErrorCode Calc_dt_calor();

PetscErrorCode write_tempo(int cont);

PetscErrorCode veloc_total();

PetscErrorCode rescaleVeloc(Vec Veloc_fut, double tempo);




int main(int argc,char **args)
{
	PetscErrorCode ierr;
	char prefix[PETSC_MAX_PATH_LEN];
	PetscInt       Px,Pz;
	
	ierr = PetscInitialize(&argc,&args,(char*)0,help);CHKERRQ(ierr);


	PetscBool      flags;
	ierr = PetscOptionsHasName(NULL,NULL,"-flags",&flags);
	if (flags){
		PetscPrintf(PETSC_COMM_WORLD,"%s",help);
		exit(1);
	}

	variable_bcv=0;
	ierr = PetscOptionsGetInt(NULL,NULL,"-variable_bcv",&variable_bcv,NULL);CHKERRQ(ierr);
	
	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	
	seed = rank;
	
	reader(rank);
	
	PetscLogDouble Tempo1,Tempo2;
	
	PetscTime(&Tempo1);
	
	
	Px   = Pz = PETSC_DECIDE;
	ierr = PetscOptionsGetInt(NULL,NULL,"-Px",&Px,NULL);CHKERRQ(ierr);
	Pz = Px;
	ierr = PetscOptionsGetInt(NULL,NULL,"-Pz",&Pz,NULL);CHKERRQ(ierr);
	
	rtol = PETSC_DEFAULT;
	ierr = PetscOptionsGetReal(NULL,NULL,"-rtol",&rtol,NULL);CHKERRQ(ierr);
	
	temper_extern = 0;
	ierr = PetscOptionsGetInt(NULL,NULL,"-te",&temper_extern,NULL);CHKERRQ(ierr);
	
	veloc_extern = 0;
	ierr = PetscOptionsGetInt(NULL,NULL,"-ve",&veloc_extern,NULL);CHKERRQ(ierr);
	
	bcv_extern = 0;
	ierr = PetscOptionsGetInt(NULL,NULL,"-bcve",&bcv_extern,NULL);CHKERRQ(ierr);
	
	
	denok_min = 1.0E-4;
	ierr = PetscOptionsGetReal(NULL,NULL,"-denok",&denok_min,NULL);CHKERRQ(ierr);
	
	print_visc = 0;
	ierr = PetscOptionsGetInt(NULL,NULL,"-print_visc",&print_visc,NULL);CHKERRQ(ierr);
	
	visc_const_per_element=0;
	ierr = PetscOptionsGetInt(NULL,NULL,"-visc_const_per_element",&visc_const_per_element,NULL);CHKERRQ(ierr);
	
	visc_harmonic_mean=1;
	ierr = PetscOptionsGetInt(NULL,NULL,"-visc_harmonic_mean",&visc_harmonic_mean,NULL);CHKERRQ(ierr);

	particles_per_ele=81;
	ierr = PetscOptionsGetInt(NULL,NULL,"-particles_per_ele",&particles_per_ele,NULL);CHKERRQ(ierr);

	particles_perturb_factor=0.5;
	ierr = PetscOptionsGetReal(NULL,NULL,"-particles_perturb_factor",&particles_perturb_factor,NULL);CHKERRQ(ierr);
	if (particles_perturb_factor>1.0) particles_perturb_factor = 1.0;
	if (particles_perturb_factor<0.0) particles_perturb_factor = 0.0;

	free_surface_stab=1;
	ierr = PetscOptionsGetInt(NULL,NULL,"-free_surface_stab",&free_surface_stab,NULL);CHKERRQ(ierr);

	sub_division_time_step=1.0;
	ierr = PetscOptionsGetReal(NULL,NULL,"-sub_division_time_step",&sub_division_time_step,NULL);CHKERRQ(ierr);

	print_step_files=1;
	ierr = PetscOptionsGetInt(NULL,NULL,"-print_step_files",&print_step_files,NULL);CHKERRQ(ierr);

	theta_FSSA=0.5;
	ierr = PetscOptionsGetReal(NULL,NULL,"-theta_FSSA",&theta_FSSA,NULL);CHKERRQ(ierr);

	direct_solver=1;
	ierr = PetscOptionsGetInt(NULL,NULL,"-direct_solver",&direct_solver,NULL);CHKERRQ(ierr);

	RK4=0;
	ierr = PetscOptionsGetInt(NULL,NULL,"-RK4",&RK4,NULL);CHKERRQ(ierr);

	Xi_min=1.0E-14;
	ierr = PetscOptionsGetReal(NULL,NULL,"-xi_min",&Xi_min,NULL);CHKERRQ(ierr);

	seed_layer=-1;
	ierr = PetscOptionsGetInt(NULL,NULL,"-seed",&seed_layer,NULL);CHKERRQ(ierr);

	checkered=0;
	ierr = PetscOptionsGetInt(NULL,NULL,"-checkered",&checkered,NULL);CHKERRQ(ierr);

	pressure_in_rheol=2;
	ierr = PetscOptionsGetInt(NULL,NULL,"-pressure_in_rheol",&pressure_in_rheol,NULL);CHKERRQ(ierr);
	if (pressure_in_rheol==2) {
		PetscPrintf(PETSC_COMM_WORLD,"Specify pressure_in_rheol!\n\n");
		exit(1);
	}

	h_air=-1.0;
	ierr = PetscOptionsGetReal(NULL,NULL,"-h_air",&h_air,NULL);CHKERRQ(ierr);
	if (pressure_in_rheol==0 && h_air<0.0){
		PetscPrintf(PETSC_COMM_WORLD,"Specify the thickness of the air layer with the flag -h_air\n");
		PetscPrintf(PETSC_COMM_WORLD,"(you adopted depth dependent rheology: -pressure_in_rheol 0)\n");
		exit(1);
	}
	else h_air=0.0;

	initial_dynamic_range=0;
	ierr = PetscOptionsGetInt(NULL,NULL,"-initial_dynamic_range",&initial_dynamic_range,NULL);CHKERRQ(ierr);


	periodic_boundary=0;
	ierr = PetscOptionsGetInt(NULL,NULL,"-periodic_boundary",&periodic_boundary,NULL);CHKERRQ(ierr);
	
	dx_const = Lx/(Nx-1);
	dz_const = depth/(Nz-1);
	
	if (rank==0) printf("%lf %lf\n",dx_const,dz_const);
	
	
	ierr = create_thermal_2d(Nx-1,Nz-1,Px,Pz);CHKERRQ(ierr);

	ierr = write_thermal_(-1);
	
	ierr = create_veloc_3d(Nx-1,Nz-1,Px,Pz);CHKERRQ(ierr);

	if (geoq_on){
		PetscPrintf(PETSC_COMM_WORLD,"Swarm INICIO\n");
		ierr = createSwarm();CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD,"Swarm FIM\n");
	}
	

	// Gerya p. 215
	if (visc_MAX>visc_MIN && initial_dynamic_range>0){
		double visc_contrast = PetscLog10Real(visc_MAX/visc_MIN);

		double visc_mean = PetscPowReal(10.0,PetscLog10Real(visc_MIN)+visc_contrast/2);

		int n_visc=0;

		visc_MIN_comp = visc_mean;
		visc_MAX_comp = visc_mean;

		PetscPrintf(PETSC_COMM_WORLD,"\n\n%.3lg %.3lg\n\n",visc_MIN_comp,visc_MAX_comp);

		ierr = veloc_total(); CHKERRQ(ierr);

		while ((visc_MIN_comp!=visc_MIN) && (visc_MAX_comp!=visc_MAX)){

			visc_MIN_comp = visc_mean*PetscPowReal(10.0,-n_visc*1.0);
			visc_MAX_comp = visc_mean*PetscPowReal(10.0,n_visc*1.0);

			if (visc_MIN_comp<visc_MIN) visc_MIN_comp=visc_MIN;
			if (visc_MAX_comp>visc_MAX) visc_MAX_comp=visc_MAX;

			PetscPrintf(PETSC_COMM_WORLD,"\n\n%.3lg %.3lg\n\n",visc_MIN_comp,visc_MAX_comp);

			ierr = veloc_total(); CHKERRQ(ierr);

			n_visc++;
		}
	}
	else {
		visc_MIN_comp = visc_MIN;
		visc_MAX_comp = visc_MAX;
		ierr = veloc_total(); CHKERRQ(ierr);
	}
	
	PetscPrintf(PETSC_COMM_WORLD,"passou veloc_total\n");
	
	ierr = write_veloc_3d(tcont);
	ierr = write_veloc_cond(tcont);
	ierr = write_thermal_(tcont);
	ierr = write_pressure(tcont);
	ierr = write_geoq_(tcont);
	ierr = write_tempo(tcont);
	
	VecCopy(Veloc_fut,Veloc);
	
	PetscPrintf(PETSC_COMM_WORLD,"passou impressao\n");
	
	ierr = Calc_dt_calor();
	
	//float aux_le;
	
	
	
	for (tempo = dt_calor,tcont=1;tempo<=timeMAX && tcont<=stepMAX;tempo+=dt_calor, tcont++){
		PetscPrintf(PETSC_COMM_WORLD,"\n\nstep = %d, time = %.3g Myr, dt = %.3g Myr\n",tcont,tempo,dt_calor);

		ierr = rescaleVeloc(Veloc_fut,tempo);
		
		ierr = build_thermal_3d();CHKERRQ(ierr);

		ierr = solve_thermal_3d();CHKERRQ(ierr);
		
		ierr = veloc_total(); CHKERRQ(ierr);
		
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
			//exit(1);
		}
		
		if (tcont%print_step==0){
			ierr = write_thermal_(tcont);
			ierr = write_geoq_(tcont);
			ierr = write_veloc_3d(tcont);
			ierr = write_pressure(tcont);
			ierr = write_tempo(tcont);
			PetscSNPrintf(prefix,PETSC_MAX_PATH_LEN-1,"step_%d",tcont);
			if (geoq_on){
				if (print_step_files==1){
					ierr = SwarmViewGP(dms,prefix);CHKERRQ(ierr);
				}
			}
		}
		

		
		//if (rank==0) scanf("%f",&aux_le);
		
		//MPI_Barrier(PETSC_COMM_WORLD);
		
		
		ierr = Calc_dt_calor();
		
	}
	if (rank==0) printf("write\n");
	
	
	
	ierr = destroy_thermal_();CHKERRQ(ierr);
	
	destroy_veloc_3d();

	PetscTime(&Tempo2);
	
	//if (rank==0) printf("Tempo: %lf\n",Tempo2-Tempo1);
	
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
	//printf("max_v = %g\n",max_v);
	//printf("max_v = %g\n",min_v);
	max_mod_v = fabs(max_v);
	ind_v_mod = ind_v_max;
	if (max_mod_v<fabs(min_v)){
		max_mod_v = fabs(min_v);
		ind_v_mod = ind_v_min;
	}
	if (ind_v_mod%2==0) dh_v_mod = dx_const; //!!! 2d
	if (ind_v_mod%2==1) dh_v_mod = dz_const; //!!! 2d
	if (rank==0) printf("dt = %g",(dh_v_mod/max_mod_v)/seg_per_ano);
	dt_calor = 0.1*(dh_v_mod/max_mod_v)/(seg_per_ano*sub_division_time_step);
	if (dt_calor>dt_MAX) dt_calor=dt_MAX;
	
	//dt_calor=1000000.0; /// !!!! apenas teste
	dt_calor_sec = dt_calor*seg_per_ano;
	////////fim calc dt
	
	PetscFunctionReturn(0);
	
}


PetscErrorCode write_tempo(int cont){

	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	
	PetscViewer viewer;
	
	char nome[100];
	
	sprintf(nome,"Tempo_%d.txt",cont);
	
	PetscReal aa[1];
	
	aa[0]=tempo;
	
	PetscViewerASCIIOpen(PETSC_COMM_WORLD,nome,&viewer);
	PetscRealView(1,aa,viewer);
	PetscViewerDestroy(&viewer);
	
	if (rank==0) printf("Tempo: %lf\n",tempo);
	
	
	
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

				PetscPrintf(PETSC_COMM_WORLD,"\n\n%.3lg %.3lg\n\n",visc_MIN_comp,visc_MAX_comp);

				ierr = veloc_total(); CHKERRQ(ierr);

				while ((visc_MIN_comp!=visc_MIN) && (visc_MAX_comp!=visc_MAX)){

					visc_MIN_comp = visc_mean*PetscPowReal(10.0,-n_visc*1.0);
					visc_MAX_comp = visc_mean*PetscPowReal(10.0,n_visc*1.0);

					if (visc_MIN_comp<visc_MIN) visc_MIN_comp=visc_MIN;
					if (visc_MAX_comp>visc_MAX) visc_MAX_comp=visc_MAX;

					PetscPrintf(PETSC_COMM_WORLD,"\n\n%.3lg %.3lg\n\n",visc_MIN_comp,visc_MAX_comp);

					ierr = veloc_total(); CHKERRQ(ierr);

					n_visc++;
				}
			}
		}
	}
	
	PetscFunctionReturn(0);
}