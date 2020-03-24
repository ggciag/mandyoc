#include <petscksp.h>
#include <petscdmda.h>

#include <petsctime.h>

PetscErrorCode montaKeThermal_general(PetscReal *Ke, PetscReal *Me, PetscReal *Fe);

PetscErrorCode DMDAGetElementCorners(DM da,PetscInt *sx,PetscInt *sz,PetscInt *mx,PetscInt *mz);

PetscErrorCode DMDAGetLocalElementSize(DM da,PetscInt *mxl,PetscInt *mzl);

PetscErrorCode DMDAGetElementCorners(DM da,PetscInt *sx,PetscInt *sz,PetscInt *mx,PetscInt *mz);

PetscErrorCode AssembleA_Thermal(Mat TA,DM thermal_da,PetscReal *TKe,PetscReal *TMe,PetscReal *TFe,
								 DM veloc_da, Vec Veloc_total);

PetscErrorCode AssembleF_Thermal(Vec F,DM thermal_da,PetscReal *TKe,PetscReal *TMe,PetscReal *TFe,
								 DM veloc_da, Vec Veloc_total);

PetscErrorCode Thermal_init(Vec F,DM thermal_da);

extern double Lx, Ly, depth;

extern PetscReal *NT;
extern PetscReal *NT_x;
extern PetscReal *NT_z;

extern long GaussQuad;

extern long T_NE;

extern PetscReal *TKe, *TCe, *TFe, *TCe_fut, *TMe, *Ttotal, *Ttotal_b;

extern PetscReal *T_vec_aux_ele;
extern PetscReal *T_vec_aux_ele_final;

extern PetscReal *v_vec_aux_ele;

extern long V_GT;
extern long V_GN;

extern Mat TA, TB;
extern Vec Tf, Temper;

extern Vec geoq,local_geoq;
extern Vec geoq_rho,local_geoq_rho;
extern Vec geoq_H,local_geoq_H;
extern Vec geoq_strain,local_geoq_strain;

extern Vec geoq_cont,local_geoq_cont;

extern Vec dRho;

extern Vec Pressure;

extern KSP T_ksp;

extern DM da_Thermal;

extern Vec Veloc, Veloc_fut;


extern DM da_Veloc;

extern Vec Temper_Cond;


extern int bcT_top;
extern int bcT_bot;
extern int bcT_left;
extern int bcT_right;

extern Vec local_FT;
extern Vec local_Temper;
extern Vec local_TC;

extern Vec local_dRho;

extern PetscReal rtol;


PetscErrorCode create_thermal_2d(PetscInt mx,PetscInt mz,PetscInt Px,PetscInt Pz)
{

	PetscInt       dof,stencil_width;
	PetscErrorCode ierr;
	
	
	
	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	
	PetscLogDouble Tempo1p,Tempo2p;
	
	PetscTime(&Tempo1p);
	
		

	
	//PetscFunctionBeginUser;
	
	
	dof           = 1;
	stencil_width = 1;
	ierr          = DMDACreate2d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,
								 mx+1,mz+1,Px,Pz,dof,stencil_width,NULL,NULL,&da_Thermal);CHKERRQ(ierr);
	ierr = DMSetFromOptions(da_Thermal);CHKERRQ(ierr);
	ierr = DMSetUp(da_Thermal);CHKERRQ(ierr);
	
	ierr = DMDASetFieldName(da_Thermal,0,"T");CHKERRQ(ierr);

	
	
	
	ierr = PetscCalloc1(GaussQuad*T_NE,&NT); CHKERRQ(ierr);
	ierr = PetscCalloc1(GaussQuad*T_NE,&NT_x); CHKERRQ(ierr);
	ierr = PetscCalloc1(GaussQuad*T_NE,&NT_z); CHKERRQ(ierr);
	
	ierr = PetscCalloc1(T_NE*T_NE,&TKe); CHKERRQ(ierr);
	ierr = PetscCalloc1(T_NE*T_NE,&TCe); CHKERRQ(ierr);
	ierr = PetscCalloc1(T_NE*T_NE,&TCe_fut); CHKERRQ(ierr);
	ierr = PetscCalloc1(T_NE*T_NE,&TMe); CHKERRQ(ierr);
	ierr = PetscCalloc1(T_NE*T_NE,&Ttotal); CHKERRQ(ierr);
	ierr = PetscCalloc1(T_NE*T_NE,&Ttotal_b); CHKERRQ(ierr);
	
	ierr = PetscCalloc1(T_NE,&TFe); CHKERRQ(ierr);
	
	ierr = PetscCalloc1(T_NE,&T_vec_aux_ele); CHKERRQ(ierr);
	ierr = PetscCalloc1(T_NE,&T_vec_aux_ele_final); CHKERRQ(ierr);
	
	ierr = PetscCalloc1(V_GT,&v_vec_aux_ele); CHKERRQ(ierr);
	
	montaKeThermal_general(TKe,TMe,TFe);
	
	
	ierr = DMDASetUniformCoordinates(da_Thermal,0.0,Lx,-depth,0.0,0.0,0.0);CHKERRQ(ierr);
	
		
	/* Generate a matrix with the correct non-zero pattern of type AIJ. This will work in parallel and serial */
	ierr = DMSetMatType(da_Thermal,MATAIJ);CHKERRQ(ierr);
	ierr = DMCreateMatrix(da_Thermal,&TA);CHKERRQ(ierr);
	ierr = DMCreateMatrix(da_Thermal,&TB);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da_Thermal,&Temper);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da_Thermal,&Tf);CHKERRQ(ierr);
	
	ierr = DMCreateGlobalVector(da_Thermal,&geoq);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da_Thermal,&geoq_rho);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da_Thermal,&geoq_H);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da_Thermal,&geoq_strain);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da_Thermal,&geoq_cont);CHKERRQ(ierr);
	
	ierr = DMCreateGlobalVector(da_Thermal,&dRho);CHKERRQ(ierr);
	
	ierr = Thermal_init(Temper,da_Thermal);
	
	//VecView(Temper,PETSC_VIEWER_STDOUT_WORLD);
	
	ierr = DMCreateLocalVector(da_Thermal,&local_FT);
	ierr = DMCreateLocalVector(da_Thermal,&local_Temper);
	ierr = DMCreateLocalVector(da_Thermal,&local_TC);

	ierr = DMCreateLocalVector(da_Thermal,&local_geoq);
	ierr = DMCreateLocalVector(da_Thermal,&local_geoq_rho);
	ierr = DMCreateLocalVector(da_Thermal,&local_geoq_H);
	ierr = DMCreateLocalVector(da_Thermal,&local_geoq_strain);
	ierr = DMCreateLocalVector(da_Thermal,&local_geoq_cont);
	
	ierr = DMCreateLocalVector(da_Thermal,&local_dRho);
	
	
	/*PetscViewer viewer;
	
	char nome[100];
	
	sprintf(nome,"Temper_0.txt");
	
	PetscViewerASCIIOpen(PETSC_COMM_WORLD,nome,&viewer);
	VecView(Temper,viewer);
	PetscViewerDestroy(&viewer);*/
	
	
	
	
	///////
	ierr = DMCreateGlobalVector(da_Thermal,&Temper_Cond);CHKERRQ(ierr);
	
	
	PetscScalar					**ff;
	PetscInt               M,P;
	
	ierr = DMDAGetInfo(da_Thermal,0,&M,&P,NULL,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	
	
	ierr = VecZeroEntries(local_FT);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da_Thermal,local_FT,&ff);CHKERRQ(ierr);
	
	PetscInt       sx,sz,mmx,mmz;
	PetscInt i,k;
	
	ierr = DMDAGetCorners(da_Thermal,&sx,&sz,NULL,&mmx,&mmz,NULL);CHKERRQ(ierr);
	
	for (k=sz; k<sz+mmz; k++) {
		for (i=sx; i<sx+mmx; i++) {
			ff[k][i] = 1.0;
			
			if (i==0   && bcT_left==1) ff[k][i] = 0.0;
			
			
			if (i==M-1 && bcT_right==1)ff[k][i] = 0.0;
			
			
			if (k==0   && bcT_bot==1) ff[k][i] = 0.0;
			
			
			if (k==P-1 && bcT_top==1) ff[k][i] = 0.0;
		}
	}
	
	ierr = DMDAVecRestoreArray(da_Thermal,local_FT,&ff);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(da_Thermal,local_FT,INSERT_VALUES,Temper_Cond);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(da_Thermal,local_FT,INSERT_VALUES,Temper_Cond);CHKERRQ(ierr);
	
	
	///////
	
	
	
	
	
	ierr = KSPCreate(PETSC_COMM_WORLD,&T_ksp);CHKERRQ(ierr);
	ierr = KSPSetOptionsPrefix(T_ksp,"thermal_"); /* stokes */ CHKERRQ(ierr);
	
	PetscTime(&Tempo2p);
	if (rank==0) printf("Thermal create : %lf\n",Tempo2p-Tempo1p);
	
	
	PetscFunctionReturn(0);
}


PetscErrorCode build_thermal_3d()
{

	PetscErrorCode ierr;
	
	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	
	PetscLogDouble Tempo1p,Tempo2p;
	
	PetscTime(&Tempo1p);
	if (rank==0) printf("TA,TB,Tf -> zero entries\n");
	ierr = MatZeroEntries(TA);CHKERRQ(ierr);
	if (rank==0) printf("passou TA\n");
	ierr = MatZeroEntries(TB);CHKERRQ(ierr);
	ierr = VecZeroEntries(Tf);CHKERRQ(ierr);
	
	if (rank==0) printf("build TA,Tf\n");
	ierr = AssembleA_Thermal(TA,da_Thermal,TKe,TMe,TFe,da_Veloc,Veloc_fut);CHKERRQ(ierr);
	//if (rank==0) printf("t\n");
	
	
	ierr = AssembleF_Thermal(Tf,da_Thermal,TKe,TMe,TFe,da_Veloc,Veloc);CHKERRQ(ierr);
	//if (rank==0) printf("t\n");

	
	PetscTime(&Tempo2p);
	if (rank==0) printf("Thermal build: %lf\n",Tempo2p-Tempo1p);
	
	PetscFunctionReturn(0);
		
}

PetscErrorCode solve_thermal_3d()
{
	PetscErrorCode ierr;
	PetscLogDouble Tempo1,Tempo2;
	
	int rank;
	
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	
	PetscTime(&Tempo1);
	
	/* SOLVE */
	
	//if (rank==0) printf("k\n");
	ierr = KSPSetOperators(T_ksp,TA,TA);CHKERRQ(ierr);
	//if (rank==0) printf("k\n");
	ierr = KSPSetFromOptions(T_ksp);CHKERRQ(ierr);
	//if (rank==0) printf("k\n");
	
	ierr = KSPSetTolerances(T_ksp,rtol,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	
	ierr = KSPSolve(T_ksp,Tf,Temper);CHKERRQ(ierr);

	PetscTime(&Tempo2);
	if (rank==0) printf("Thermal solve: %lf\n",Tempo2-Tempo1);
	
	PetscFunctionReturn(0);
	
}
/*
air_temperature()
{
	ierr = VecZeroEntries(local_Temper);CHKERRQ(ierr);
	
	ierr = DMGlobalToLocalBegin(thermal_da,Temper,INSERT_VALUES,local_Temper);
	ierr = DMGlobalToLocalEnd(thermal_da,Temper,INSERT_VALUES,local_Temper);
	
	ierr = DMDAVecGetArray(thermal_da,local_Temper,&tt);CHKERRQ(ierr);
	
	PetscInt       sx,sz,mmx,mmz;
	
	ierr = DMDAGetCorners(thermal_da,&sx,&sz,NULL,&mmx,&mmz,NULL);CHKERRQ(ierr);
	
	PetscReal xx,zz,t_inic;
	
	for (k=sz; k<sz+mmz; k++) {
		for (i=sx; i<sx+mmx; i++) {
			
		}
	}
	
	ierr = DMDAVecRestoreArray(thermal_da,local_Temper,&tt);CHKERRQ(ierr);
	
	ierr = DMLocalToGlobalBegin(thermal_da,local_Temper,INSERT_VALUES,Temper);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(thermal_da,local_Temper,INSERT_VALUES,Temper);CHKERRQ(ierr);
	
}
*/

PetscErrorCode destroy_thermal_()
{

	int rank;
	
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	PetscLogDouble Tempo1,Tempo2;
	PetscTime(&Tempo1);
	
	PetscErrorCode ierr;
	ierr = KSPDestroy(&T_ksp);CHKERRQ(ierr);
	ierr = VecDestroy(&Temper);CHKERRQ(ierr);
	ierr = VecDestroy(&Tf);CHKERRQ(ierr);
	ierr = MatDestroy(&TA);CHKERRQ(ierr);
	ierr = MatDestroy(&TB);CHKERRQ(ierr);
	ierr = DMDestroy(&da_Thermal);CHKERRQ(ierr);
	
	PetscTime(&Tempo2);
	if (rank==0) printf("Thermal destroy: %lf\n",Tempo2-Tempo1);
	
	
	PetscFunctionReturn(0);	
}

PetscErrorCode write_thermal_(int cont)
{
	int rank;
	
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	PetscLogDouble Tempo1,Tempo2;
	PetscTime(&Tempo1);
	
	PetscViewer viewer;
	
	char nome[100];
	
	sprintf(nome,"Temper_%d.txt",cont);
	
	PetscViewerASCIIOpen(PETSC_COMM_WORLD,nome,&viewer);
	VecView(Temper,viewer);
	PetscViewerDestroy(&viewer);
	
	PetscTime(&Tempo2);
	if (rank==0) printf("Thermal write: %lf\n",Tempo2-Tempo1);
	
	PetscFunctionReturn(0);	
}

PetscErrorCode write_pressure(int cont)
{
	int rank;
	
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	PetscLogDouble Tempo1,Tempo2;
	PetscTime(&Tempo1);
	
	PetscViewer viewer;
	
	char nome[100];
	
	sprintf(nome,"Pressure_%d.txt",cont);
	
	PetscViewerASCIIOpen(PETSC_COMM_WORLD,nome,&viewer);
	VecView(Pressure,viewer);
	PetscViewerDestroy(&viewer);
	
	PetscTime(&Tempo2);
	if (rank==0) printf("Pressure write: %lf\n",Tempo2-Tempo1);
	
	PetscFunctionReturn(0);
}

PetscErrorCode write_geoq_(int cont)
{
	int rank;
	
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	PetscLogDouble Tempo1,Tempo2;
	PetscTime(&Tempo1);
	
	PetscViewer viewer;
	
	char nome[100];
	
	sprintf(nome,"Geoq_%d.txt",cont);
	
	PetscViewerASCIIOpen(PETSC_COMM_WORLD,nome,&viewer);
	VecView(geoq,viewer);
	PetscViewerDestroy(&viewer);
	
	PetscTime(&Tempo2);
	if (rank==0) printf("Geoq write: %lf\n",Tempo2-Tempo1);
	
	PetscTime(&Tempo1);
	

	sprintf(nome,"Rho_%d.txt",cont);
	
	PetscViewerASCIIOpen(PETSC_COMM_WORLD,nome,&viewer);
	VecView(geoq_rho,viewer);
	PetscViewerDestroy(&viewer);
	
	PetscTime(&Tempo2);
	if (rank==0) printf("Rho write: %lf\n",Tempo2-Tempo1);
	
	
	sprintf(nome,"H_%d.txt",cont);
	
	PetscViewerASCIIOpen(PETSC_COMM_WORLD,nome,&viewer);
	VecView(geoq_H,viewer);
	PetscViewerDestroy(&viewer);
	
	PetscTime(&Tempo2);
	if (rank==0) printf("H write: %lf\n",Tempo2-Tempo1);
	
	sprintf(nome,"strain_%d.txt",cont);
	
	PetscViewerASCIIOpen(PETSC_COMM_WORLD,nome,&viewer);
	VecView(geoq_strain,viewer);
	PetscViewerDestroy(&viewer);
	
	PetscTime(&Tempo2);
	if (rank==0) printf("H write: %lf\n",Tempo2-Tempo1);
	
	PetscFunctionReturn(0);
}

