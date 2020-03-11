static char help[] = "S1\n";

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





int main(int argc,char **args)
{
	PetscErrorCode ierr;
	char prefix[PETSC_MAX_PATH_LEN];
	PetscInt       Px,Pz;
	
	ierr = PetscInitialize(&argc,&args,(char*)0,help);CHKERRQ(ierr);
	
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

	free_surface_stab=1;
	ierr = PetscOptionsGetInt(NULL,NULL,"-free_surface_stab",&free_surface_stab,NULL);CHKERRQ(ierr);

	sub_division_time_step=1.0;
	ierr = PetscOptionsGetReal(NULL,NULL,"-sub_division_time_step",&sub_division_time_step,NULL);CHKERRQ(ierr);

	print_step_files=1;
	ierr = PetscOptionsGetInt(NULL,NULL,"-print_step_files",&print_step_files,NULL);CHKERRQ(ierr);

	theta_FSSA=0.5;
	ierr = PetscOptionsGetReal(NULL,NULL,"-theta_FSSA",&theta_FSSA,NULL);CHKERRQ(ierr);
	
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
	
	
	
	ierr = veloc_total(); CHKERRQ(ierr);
	
	
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
		
		
		
		ierr = build_thermal_3d();CHKERRQ(ierr);

		ierr = solve_thermal_3d();CHKERRQ(ierr);
		
		ierr = veloc_total(); CHKERRQ(ierr);
		
		if (geoq_on){
			for (PetscInt cont=0, max_cont=10;cont<max_cont; cont++){
				double fac = (1.0/max_cont)*(0.5+cont);
				//PetscPrintf(PETSC_COMM_WORLD,"%f %f\n",fac,(1.0-fac));
				//        x   ->   y
				VecCopy(Veloc,Veloc_weight);
				//			y			a		b		x
				VecAXPBY(Veloc_weight, fac, (1.0-fac),Veloc_fut); //y = a*x + b*y
				ierr = moveSwarm(dt_calor_sec/max_cont);
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
	//if (ind_v_mod%3==1) dh_v_mod = dy_const;
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

