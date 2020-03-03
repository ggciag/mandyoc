
#include <petscksp.h>

void ErrorInterfaces();


extern long Nx,Nz;
extern long layers;

extern double Lx, depth;

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

extern double H_lito;

extern int n_interfaces;
extern PetscScalar *interfaces;

extern PetscScalar *inter_rho;
extern PetscScalar *inter_geoq;
extern PetscScalar *inter_H;

extern PetscScalar *inter_A;
extern PetscScalar *inter_n;
extern PetscScalar *inter_Q;
extern PetscScalar *inter_V;

extern double H_per_mass;
extern double c_heat_capacity;

extern int T_initial_cond;
extern int rheol;

extern double beta_max;
extern double ramp_begin;
extern double ramp_end;

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

extern double h_air;

extern PetscInt WITH_NON_LINEAR;
extern PetscInt WITH_ADIABATIC_H;
extern PetscInt WITH_RADIOGENIC_H;

char str[100];

PetscErrorCode reader(int rank){
	if (rank==0){
		FILE *f_parametros;
		
		f_parametros = fopen("param_1.5.3_2D.txt","r");
		
		fscanf(f_parametros,"%ld %ld",&Nx,&Nz);
		fscanf(f_parametros,"%lg %lg",&Lx,&depth);
		
		layers=Nz;
		
		
		fscanf(f_parametros,"%s",str);
		if (strcmp (str,"mg") == 0) fscanf(f_parametros,"%d",&ContMult);
		else {printf("mg error\n"); exit(1);}
		
		
		fscanf(f_parametros,"%s",str);
		if (strcmp (str,"stepMAX") == 0) fscanf(f_parametros,"%ld",&stepMAX);
		else {printf("stepMAX error\n"); exit(1);}
		
		fscanf(f_parametros,"%s",str);
		if (strcmp (str,"timeMAX") == 0) fscanf(f_parametros,"%lf",&timeMAX);
		else {printf("timeMAX error\n"); exit(1);}
		
		fscanf(f_parametros,"%s",str);
		if (strcmp (str,"dt_MAX") == 0) fscanf(f_parametros,"%lf",&dt_MAX);
		else {printf("dt_MAX error\n"); exit(1);}
		
		fscanf(f_parametros,"%s",str);
		if (strcmp (str,"print_step") == 0) fscanf(f_parametros,"%ld",&print_step);
		else {printf("print_step error\n"); exit(1);}
		
		
		fscanf(f_parametros,"%s",str);
		if (strcmp (str,"visc") == 0) fscanf(f_parametros,"%lg",&visco_r);
		else {printf("visc error\n"); exit(1);}
		
		fscanf(f_parametros,"%s",str);
		if (strcmp (str,"visc_MAX") == 0) fscanf(f_parametros,"%lg",&visc_MAX);
		else {printf("visc_MAX error\n"); exit(1);}
		
		fscanf(f_parametros,"%s",str);
		if (strcmp (str,"visc_MIN") == 0) fscanf(f_parametros,"%lg",&visc_MIN);
		else {printf("visc_MIN error\n"); exit(1);}
		
		
		
		fscanf(f_parametros,"%s",str);
		if (strcmp (str,"n_interfaces") == 0) fscanf(f_parametros,"%d",&n_interfaces);
		else {printf("n_interfaces error\n"); exit(1);}
		
		
		
		
		fscanf(f_parametros,"%s",str);
		if (strcmp (str,"geoq_on") == 0) fscanf(f_parametros,"%d",&geoq_on);
		else {printf("geoq_on error\n"); exit(1);}
		
		fscanf(f_parametros,"%s",str);
		if (strcmp (str,"geoq_fac") == 0) fscanf(f_parametros,"%lg",&escala_viscosidade);
		else {printf("geoq_fac error\n"); exit(1);}
		
		fscanf(f_parametros,"%s",str);
		if (strcmp (str,"veloc") == 0) fscanf(f_parametros,"%lg",&veloc_superf);
		else {printf("veloc error\n"); exit(1);}
		
		
		
		fscanf(f_parametros,"%s",str);
		if (strcmp (str,"deltaT") == 0) fscanf(f_parametros,"%lg",&Delta_T);
		else {printf("deltaT error\n"); exit(1);}
		
		fscanf(f_parametros,"%s",str);
		if (strcmp (str,"alpha_exp_thermo") == 0) fscanf(f_parametros,"%lg",&alpha_exp_thermo);
		else {printf("alpha_exp_thermo error\n"); exit(1);}
		
		fscanf(f_parametros,"%s",str);
		if (strcmp (str,"kappa") == 0) fscanf(f_parametros,"%lg",&kappa);
		else {printf("kappa error\n"); exit(1);}
		
		fscanf(f_parametros,"%s",str);
		if (strcmp (str,"gravity") == 0) fscanf(f_parametros,"%lg",&gravity);
		else {printf("gravity error\n"); exit(1);}
		
		fscanf(f_parametros,"%s",str);
		if (strcmp (str,"rhom") == 0) fscanf(f_parametros,"%lg",&RHOM);
		else {printf("rhom error\n"); exit(1);}
		
		
		fscanf(f_parametros,"%s",str);
		if (strcmp (str,"H_per_mass") == 0) fscanf(f_parametros,"%lg",&H_per_mass);
		else {printf("H_per_mass error\n"); exit(1);}
		
		fscanf(f_parametros,"%s",str);
		if (strcmp (str,"c_heat_capacity") == 0) fscanf(f_parametros,"%lg",&c_heat_capacity);
		else {printf("c_heat_capacity error\n"); exit(1);}
		
		
		fscanf(f_parametros,"%s",str);
		if (strcmp (str,"non_linear") == 0) fscanf(f_parametros,"%d",&WITH_NON_LINEAR);
		else {printf("non_linear error\n"); exit(1);}
		
		fscanf(f_parametros,"%s",str);
		if (strcmp (str,"adiabatic_H") == 0) fscanf(f_parametros,"%d",&WITH_ADIABATIC_H);
		else {printf("adiabatic_H error\n"); exit(1);}
		
		fscanf(f_parametros,"%s",str);
		if (strcmp (str,"radiogenic_H") == 0) fscanf(f_parametros,"%d",&WITH_RADIOGENIC_H);
		else {printf("radiogenic_H error\n"); exit(1);}
		
		
		
		fscanf(f_parametros,"%s",str);
		if (strcmp (str,"bcv_top_normal") == 0) fscanf(f_parametros,"%d",&bcv_top_normal);
		else {printf("bcv_top_normal error\n"); exit(1);}
		
		fscanf(f_parametros,"%s",str);
		if (strcmp (str,"bcv_top_slip") == 0) fscanf(f_parametros,"%d",&bcv_top_slip);
		else {printf("bcv_top_slip error\n"); exit(1);}
		
		
		fscanf(f_parametros,"%s",str);
		if (strcmp (str,"bcv_bot_normal") == 0) fscanf(f_parametros,"%d",&bcv_bot_normal);
		else {printf("bcv_bot_normal error\n"); exit(1);}
		
		fscanf(f_parametros,"%s",str);
		if (strcmp (str,"bcv_bot_slip") == 0) fscanf(f_parametros,"%d",&bcv_bot_slip);
		else {printf("bcv_bot_slip error\n"); exit(1);}
		
		
		fscanf(f_parametros,"%s",str);
		if (strcmp (str,"bcv_left_normal") == 0) fscanf(f_parametros,"%d",&bcv_left_normal);
		else {printf("bcv_left_normal error\n"); exit(1);}
		
		fscanf(f_parametros,"%s",str);
		if (strcmp (str,"bcv_left_slip") == 0) fscanf(f_parametros,"%d",&bcv_left_slip);
		else {printf("bcv_left_slip error\n"); exit(1);}
		
		
		fscanf(f_parametros,"%s",str);
		if (strcmp (str,"bcv_right_normal") == 0) fscanf(f_parametros,"%d",&bcv_right_normal);
		else {printf("bcv_right_normal error\n"); exit(1);}
		
		fscanf(f_parametros,"%s",str);
		if (strcmp (str,"bcv_right_slip") == 0) fscanf(f_parametros,"%d",&bcv_right_slip);
		else {printf("bcv_right_slip error\n"); exit(1);}
		
		
		
		fscanf(f_parametros,"%s",str);
		if (strcmp (str,"bcT_top") == 0) fscanf(f_parametros,"%d",&bcT_top);
		else {printf("bcT_top error\n"); exit(1);}
		
		
		fscanf(f_parametros,"%s",str);
		if (strcmp (str,"bcT_bot") == 0) fscanf(f_parametros,"%d",&bcT_bot);
		else {printf("bcT_bot error\n"); exit(1);}
		
		
		fscanf(f_parametros,"%s",str);
		if (strcmp (str,"bcT_left") == 0) fscanf(f_parametros,"%d",&bcT_left);
		else {printf("bcT_left error\n"); exit(1);}
		
		
		fscanf(f_parametros,"%s",str);
		if (strcmp (str,"bcT_right") == 0) fscanf(f_parametros,"%d",&bcT_right);
		else {printf("bcT_right error\n"); exit(1);}
		
		
		
		
		
		fscanf(f_parametros,"%s",str);
		if (strcmp (str,"rheol") == 0) fscanf(f_parametros,"%d",&rheol);
		else {printf("rheol error\n"); exit(1);}
		
		fscanf(f_parametros,"%s",str);
		if (strcmp (str,"T_initial") == 0) fscanf(f_parametros,"%d",&T_initial_cond);
		else {printf("T_initial error\n"); exit(1);}
		
		
		fscanf(f_parametros,"%s",str);
		if (strcmp (str,"H_lito") == 0) fscanf(f_parametros,"%lf",&H_lito);
		else {printf("H_lito error\n"); exit(1);}
		
		fscanf(f_parametros,"%s",str);
		if (strcmp (str,"h_air") == 0) fscanf(f_parametros,"%lf",&h_air);
		else {printf("h_air error\n"); exit(1);}
		
		
		fscanf(f_parametros,"%s",str);
		if (strcmp (str,"beta_max") == 0) fscanf(f_parametros,"%lf",&beta_max);
		else {printf("beta_max error\n"); exit(1);}
		
		fscanf(f_parametros,"%s",str);
		if (strcmp (str,"ramp_begin") == 0) fscanf(f_parametros,"%lf",&ramp_begin);
		else {printf("ramp_begin error\n"); exit(1);}
		
		fscanf(f_parametros,"%s",str);
		if (strcmp (str,"ramp_end") == 0) fscanf(f_parametros,"%lf",&ramp_end);
		else {printf("ramp_end error\n"); exit(1);}
		
		fclose(f_parametros);
		
	}
	
	MPI_Bcast(&Nx,1,MPI_LONG,0,PETSC_COMM_WORLD);
	MPI_Bcast(&Nz,1,MPI_LONG,0,PETSC_COMM_WORLD);
	
	MPI_Bcast(&Lx,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
	MPI_Bcast(&depth,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
	
	MPI_Bcast(&ContMult,1,MPI_INT,0,PETSC_COMM_WORLD);
	
	MPI_Bcast(&stepMAX,1,MPI_LONG,0,PETSC_COMM_WORLD);
	
	MPI_Bcast(&timeMAX,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
	
	MPI_Bcast(&dt_MAX,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
	
	MPI_Bcast(&print_step,1,MPI_LONG,0,PETSC_COMM_WORLD);
	
	MPI_Bcast(&visco_r,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
	MPI_Bcast(&visc_MAX,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
	MPI_Bcast(&visc_MIN,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
	
	MPI_Bcast(&n_interfaces,1,MPI_INT,0,PETSC_COMM_WORLD);

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
	
	MPI_Bcast(&H_lito,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
	
	MPI_Bcast(&h_air,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
	
	MPI_Bcast(&beta_max,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
	
	MPI_Bcast(&ramp_begin,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
	MPI_Bcast(&ramp_end,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
	
	
	
	if (rank==0){
		printf("%ld %ld\n",Nx,Nz);
		printf("%lf %lf\n",Lx,depth);
		printf("%lf\n",beta_max);
		printf("%lf\n%lf\n",ramp_begin,ramp_end);
	}
	
	
	
	//!!!x if (n_interfaces>0){
		PetscCalloc1(Nx*n_interfaces,&interfaces);
		PetscCalloc1(n_interfaces+1,&inter_geoq);
		PetscCalloc1(n_interfaces+1,&inter_rho);
		PetscCalloc1(n_interfaces+1,&inter_H);
		
		PetscCalloc1(n_interfaces+1,&inter_A);
		PetscCalloc1(n_interfaces+1,&inter_n);
		PetscCalloc1(n_interfaces+1,&inter_Q);
		PetscCalloc1(n_interfaces+1,&inter_V);
		
		FILE *f_inter;
		f_inter = fopen("interfaces_creep.txt","r");
		if (f_inter==NULL) {
			printf("\n\n\n\ninterfaces_creep.txt not found\n\n\n\n");
			PetscFunctionReturn(-1);
		}
		if (rank==0){
			
			int check_fscanf;
			
			fscanf(f_inter,"%s",str);
			if (strcmp (str,"C") == 0){
				for (PetscInt i=0;i<n_interfaces+1;i++){
					check_fscanf = fscanf(f_inter,"%lf",&inter_geoq[i]);
					if (check_fscanf==0) {ErrorInterfaces(); exit(1);}
				}
			}
			else { ErrorInterfaces(); exit(1);}
			
			fscanf(f_inter,"%s",str);
			if (strcmp (str,"rho") == 0)
				for (PetscInt i=0;i<n_interfaces+1;i++)
					fscanf(f_inter,"%lf",&inter_rho[i]);
			else { ErrorInterfaces(); exit(1);}
			
			fscanf(f_inter,"%s",str);
			if (strcmp (str,"H") == 0)
				for (PetscInt i=0;i<n_interfaces+1;i++)
					fscanf(f_inter,"%lf",&inter_H[i]);
			else { ErrorInterfaces(); exit(1);}
			
			fscanf(f_inter,"%s",str);
			if (strcmp (str,"A") == 0)
				for (PetscInt i=0;i<n_interfaces+1;i++)
					fscanf(f_inter,"%lf",&inter_A[i]);
			else { ErrorInterfaces(); exit(1);}
			
			fscanf(f_inter,"%s",str);
			if (strcmp (str,"n") == 0)
				for (PetscInt i=0;i<n_interfaces+1;i++)
					fscanf(f_inter,"%lf",&inter_n[i]);
			else { ErrorInterfaces(); exit(1);}
			
			fscanf(f_inter,"%s",str);
			if (strcmp (str,"Q") == 0)
				for (PetscInt i=0;i<n_interfaces+1;i++)
					fscanf(f_inter,"%lf",&inter_Q[i]);
			else { ErrorInterfaces(); exit(1);}
			
			fscanf(f_inter,"%s",str);
			if (strcmp (str,"V") == 0)
				for (PetscInt i=0;i<n_interfaces+1;i++)
					fscanf(f_inter,"%lf",&inter_V[i]);
			else { ErrorInterfaces(); exit(1);}
			
			for (PetscInt i=0; i<Nx; i++){
				for (PetscInt j=0; j<n_interfaces; j++){
					fscanf(f_inter,"%lf",&interfaces[j*Nx+i]);
				}
			}
			
			printf("\nGeoq: ");
			for (PetscInt i=0;i<n_interfaces+1;i++)
				printf("%e ",inter_geoq[i]);
			printf("\nRho: ");
			for (PetscInt i=0;i<n_interfaces+1;i++)
				printf("%e ",inter_rho[i]);
			printf("\nH: ");
			for (PetscInt i=0;i<n_interfaces+1;i++)
				printf("%e ",inter_H[i]);
			printf("\nA: ");
			for (PetscInt i=0;i<n_interfaces+1;i++)
				printf("%e ",inter_A[i]);
			printf("\nn: ");
			for (PetscInt i=0;i<n_interfaces+1;i++)
				printf("%lf ",inter_n[i]);
			printf("\nQ: ");
			for (PetscInt i=0;i<n_interfaces+1;i++)
				printf("%e ",inter_Q[i]);
			printf("\nV: ");
			for (PetscInt i=0;i<n_interfaces+1;i++)
				printf("%e ",inter_V[i]);
			printf("\n\n");
		}
		
		MPI_Bcast(interfaces,Nx*n_interfaces,MPIU_SCALAR,0,PETSC_COMM_WORLD);
		MPI_Bcast(inter_geoq,n_interfaces+1,MPIU_SCALAR,0,PETSC_COMM_WORLD);
		MPI_Bcast(inter_rho,n_interfaces+1,MPIU_SCALAR,0,PETSC_COMM_WORLD);
		MPI_Bcast(inter_H,n_interfaces+1,MPIU_SCALAR,0,PETSC_COMM_WORLD);
		
		MPI_Bcast(inter_A,n_interfaces+1,MPIU_SCALAR,0,PETSC_COMM_WORLD);
		MPI_Bcast(inter_n,n_interfaces+1,MPIU_SCALAR,0,PETSC_COMM_WORLD);
		MPI_Bcast(inter_Q,n_interfaces+1,MPIU_SCALAR,0,PETSC_COMM_WORLD);
		MPI_Bcast(inter_V,n_interfaces+1,MPIU_SCALAR,0,PETSC_COMM_WORLD);
		
	//!!!x }
	
	PetscFunctionReturn(0);

}

void ErrorInterfaces(){
	printf("\n\n\n\n");
	printf("Problem to read interfaces_creep.txt:\n");
	printf("Check if the number of interfaces is\n");
	printf("compatible to the one indicated in the\n");
	printf("parameters file.\n\n\n\n");
}
