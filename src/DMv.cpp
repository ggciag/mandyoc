#include <petscksp.h>
#include <petscdmda.h>
#include <petscmath.h>
#include <petsctime.h>


extern int bcv_top_normal;
extern int bcv_top_slip;

extern int bcv_bot_normal;
extern int bcv_bot_slip;

extern int bcv_left_normal;
extern int bcv_left_slip;

extern int bcv_right_normal;
extern int bcv_right_slip;

extern PetscInt bcv_extern;

extern PetscInt visc_const_per_element;


typedef struct {
	PetscScalar u;
	PetscScalar w;
	//PetscScalar p;
} Stokes;

PetscErrorCode ascii2bin(char *s1, char *s2);

PetscErrorCode AssembleA_Veloc(Mat A,Mat AG,DM veloc_da, DM temper_da);

PetscErrorCode AssembleF_Veloc(Vec F,DM veloc_da,DM drho_da, Vec FP);

PetscErrorCode montaKeVeloc_general(PetscReal *KeG, double dx_const, double dz_const);

PetscReal montaKeVeloc_simplif(PetscReal *Ke,PetscReal *KeG, PetscReal *geoq_ele);

PetscErrorCode montaCeVeloc(PetscReal *Ce);

PetscErrorCode montafeVeloc(PetscReal *fMe);

PetscErrorCode calc_drho();

PetscErrorCode calc_pressure();

PetscErrorCode write_veloc_3d(int cont);

PetscErrorCode write_veloc_cond(int cont);

PetscErrorCode write_pressure(int cont);

PetscErrorCode Init_Veloc();

PetscErrorCode shift_pressure();


PetscErrorCode moveSwarm(PetscReal dt);
PetscErrorCode Swarm2Mesh();

extern double r06;
extern double r8p9;
extern double r5p9;

extern long V_NE, V_GN, V_GT;

extern long GaussQuad;

extern Mat VA, VB, VG;
extern Vec Vf,Vf_P, Veloc, Veloc_fut, Veloc_weight,Veloc_0;

extern Vec Veloc_step1;
extern Vec Veloc_step2;


extern Vec rk_vec2;

extern Vec rk_vec;
extern Vec sk_vec;
extern Vec gs_vec;
extern Vec uk_vec;

extern Vec zk_vec;
extern Vec zk_vec2;

extern Vec Veloc_Cond;

extern Vec Pressure;

extern int PRESSURE_INIT;

extern DM da_Veloc;
extern DM da_Thermal;

extern KSP V_ksp;

extern double Lx, Ly, depth;

extern PetscReal *Ke_veloc;
extern PetscReal *Ke_veloc_final;
extern PetscReal *Ke_veloc_general;

extern PetscReal *VCe;
extern PetscReal *VfMe;

extern PetscReal *Vfe;


extern double dx_const;
extern double dz_const;


extern Vec local_V;
extern Vec local_VC;
extern Vec local_FV;
extern Vec local_FP;
extern Vec local_P;

extern PetscReal rtol;

extern PetscReal denok_min;



extern Vec Precon;
extern Vec local_Precon;

extern double visc_aux_MAX;
extern double visc_aux_MIN;

extern double dz_const;

extern Vec Adiag;

extern PetscInt direct_solver;

extern PetscInt Verif_first_veloc;


PetscErrorCode create_veloc_3d(PetscInt mx,PetscInt mz,PetscInt Px,PetscInt Pz)
{

	PetscInt       dof,stencil_width;

	PetscErrorCode ierr;
	
	
	
	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	
	PetscLogDouble Tempo1p,Tempo2p;
	
	PetscTime(&Tempo1p);
	
	
	PetscFunctionBeginUser;
	
	
	dof           = 2; //modif
	stencil_width = 1;
	ierr          = DMDACreate2d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,
								 mx+1,mz+1,Px,Pz,dof,stencil_width,NULL,NULL,&da_Veloc);CHKERRQ(ierr);
	ierr = DMSetFromOptions(da_Veloc);CHKERRQ(ierr);
	ierr = DMSetUp(da_Veloc);CHKERRQ(ierr);
	
	ierr = DMDASetFieldName(da_Veloc,0,"V_x");CHKERRQ(ierr);
	ierr = DMDASetFieldName(da_Veloc,1,"V_z");CHKERRQ(ierr);
	//ierr = DMDASetFieldName(da_Veloc,3,"P");CHKERRQ(ierr);
	
	//printf("a\n");
	
	
	
	
	ierr = PetscCalloc1(V_GT*V_GT,&Ke_veloc); CHKERRQ(ierr);
	ierr = PetscCalloc1(V_GT*V_GT,&Ke_veloc_final); CHKERRQ(ierr);

	ierr = PetscCalloc1(V_GT*V_GT*GaussQuad,&Ke_veloc_general); CHKERRQ(ierr);
	//ierr = PetscCalloc1(V_GT,&Indices_Ke_veloc_I); CHKERRQ(ierr);
	//ierr = PetscCalloc1(V_GT,&Indices_Ke_veloc_J); CHKERRQ(ierr);
	
	ierr = PetscCalloc1(V_GT,&Vfe); CHKERRQ(ierr);
	//ierr = PetscCalloc1(V_GT,&Indices_Vfe); CHKERRQ(ierr);
	ierr = PetscCalloc1(V_GT*V_NE,&VfMe); CHKERRQ(ierr); // 24 x 8 no elemento hexa
	ierr = PetscCalloc1(V_GT,&VCe); CHKERRQ(ierr);
	
	
	ierr = DMDASetUniformCoordinates(da_Veloc,0.0,Lx,-depth,0.0,-1,-1);CHKERRQ(ierr); //!!!! 2d : the last two are ignored
	

	
	//printf("a\n");
	
	

	ierr = DMSetMatType(da_Veloc,MATAIJ);CHKERRQ(ierr);
	ierr = DMCreateMatrix(da_Veloc,&VA);CHKERRQ(ierr);
	ierr = DMCreateMatrix(da_Veloc,&VB);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da_Veloc,&Veloc);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da_Veloc,&Veloc_fut);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da_Veloc,&Veloc_weight);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da_Veloc,&Veloc_0);CHKERRQ(ierr);
	
	ierr = DMCreateGlobalVector(da_Veloc,&Veloc_step1);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da_Veloc,&Veloc_step2);CHKERRQ(ierr);
	
	ierr = DMCreateGlobalVector(da_Veloc,&Vf);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da_Veloc,&Vf_P);CHKERRQ(ierr);
	
	ierr = DMCreateGlobalVector(da_Veloc,&Pressure);CHKERRQ(ierr);
	
	ierr = DMCreateMatrix(da_Veloc,&VG);CHKERRQ(ierr);
	
	ierr = DMCreateGlobalVector(da_Veloc,&rk_vec2);CHKERRQ(ierr);
	
	ierr = DMCreateGlobalVector(da_Veloc,&rk_vec);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da_Veloc,&sk_vec);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da_Veloc,&gs_vec);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da_Veloc,&uk_vec);CHKERRQ(ierr);
	
	ierr = DMCreateGlobalVector(da_Veloc,&Adiag);CHKERRQ(ierr);
	
	ierr = DMCreateGlobalVector(da_Veloc,&zk_vec);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da_Veloc,&zk_vec2);CHKERRQ(ierr);
	
	ierr = DMCreateGlobalVector(da_Veloc,&Precon);CHKERRQ(ierr);

	ierr = DMCreateGlobalVector(da_Veloc,&Veloc_Cond);CHKERRQ(ierr);
	
	
	ierr = DMCreateLocalVector(da_Veloc,&local_V);CHKERRQ(ierr);
	ierr = DMCreateLocalVector(da_Veloc,&local_VC);CHKERRQ(ierr);
	ierr = DMCreateLocalVector(da_Veloc,&local_FV);CHKERRQ(ierr);
	ierr = DMCreateLocalVector(da_Veloc,&local_FP);CHKERRQ(ierr);
	ierr = DMCreateLocalVector(da_Veloc,&local_P);CHKERRQ(ierr);
	
	ierr = DMCreateLocalVector(da_Veloc,&local_Precon);CHKERRQ(ierr);
	
	
	Stokes					**ff;
	PetscInt               M,P;
	
	ierr = DMDAGetInfo(da_Veloc,0,&M,&P,NULL,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	
	
	ierr = VecZeroEntries(local_FV);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da_Veloc,local_FV,&ff);CHKERRQ(ierr);
	
	PetscInt       sx,sz,mmx,mmz;
	PetscInt i,k,t;
	
	ierr = DMDAGetCorners(da_Veloc,&sx,&sz,NULL,&mmx,&mmz,NULL);CHKERRQ(ierr);
	
	PetscInt ix[1];
	PetscScalar y[1];
	
	PetscInt low,high;
	
	if (bcv_extern==1){
		char s1[100],s2[100];
		
		sprintf(s1,"bcv_0_3D.txt");
		sprintf(s2,"bcv_init.bin");
		
		if (rank==0){
			ierr = ascii2bin(s1,s2); CHKERRQ(ierr);
		}
		MPI_Barrier(PETSC_COMM_WORLD);
		
		
		PetscInt size0;
		PetscViewer    viewer;
		
		VecGetSize(Veloc_Cond,&size0);
		
		
		Vec Fprov;
		
		//PetscPrintf(PETSC_COMM_WORLD,"size = %d\n",size);
		
		//PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Temper_init.bin",FILE_MODE_READ,&viewer);
		PetscViewerBinaryOpen(PETSC_COMM_WORLD,s2,FILE_MODE_READ,&viewer);
		VecCreate(PETSC_COMM_WORLD,&Fprov);
		VecLoad(Fprov,viewer);
		PetscViewerDestroy(&viewer);
		
		
		VecGetOwnershipRange(Fprov,&low,&high);
		
		printf("%d %d\n",low,high);
		
		
		Vec FN;
		
		DMDACreateNaturalVector(da_Veloc,&FN);
		
		
		for (t=low;t<high;t++){
			ix[0] = t;
			VecGetValues(Fprov,1,ix,y);
			VecSetValue(FN,t,y[0], INSERT_VALUES);
		}
		
		VecAssemblyBegin(FN);
		VecAssemblyEnd(FN);
		
		DMDANaturalToGlobalBegin(da_Veloc,FN,INSERT_VALUES,Veloc_Cond);
		DMDANaturalToGlobalEnd(da_Veloc,FN,INSERT_VALUES,Veloc_Cond);
		
		
		//VecView(F,PETSC_VIEWER_STDOUT_WORLD);
		
		PetscBarrier(NULL);
		
		VecAssemblyBegin(Veloc_Cond);
		VecAssemblyEnd(Veloc_Cond);

		char nome[100];
		PetscViewer viewer2;
		sprintf(nome,"Init_bcve.txt");
		PetscViewerASCIIOpen(PETSC_COMM_WORLD,nome,&viewer2);
		VecView(Veloc_Cond,viewer2);
		PetscViewerDestroy(&viewer2);
		//VecView(F,PETSC_VIEWER_STDOUT_WORLD);
		
	}
	else {
		for (k=sz; k<sz+mmz; k++) {
			for (i=sx; i<sx+mmx; i++) {
				ff[k][i].u = 1.0;
				ff[k][i].w = 1.0;
				
				if (i==0   && bcv_left_normal==1) ff[k][i].u = 0.0;
				if (i==0   && bcv_left_slip==1) {
					ff[k][i].w = 0.0;
				}
				
				if (i==M-1 && bcv_right_normal==1)ff[k][i].u = 0.0;
				if (i==M-1   && bcv_right_slip==1) {
					ff[k][i].w = 0.0;
				}
				
				if (k==0   && bcv_bot_normal==1) ff[k][i].w = 0.0;
				if (k==0   && bcv_bot_slip==1){
					ff[k][i].u = 0.0;
				}
				
				if (k==P-1 && bcv_top_normal==1) ff[k][i].w = 0.0;
				if (k==P-1 && bcv_top_slip==1){
					ff[k][i].u = 0.0;
				}
				
			}
		}
		ierr = DMDAVecRestoreArray(da_Veloc,local_FV,&ff);CHKERRQ(ierr);
		ierr = DMLocalToGlobalBegin(da_Veloc,local_FV,INSERT_VALUES,Veloc_Cond);CHKERRQ(ierr);
		ierr = DMLocalToGlobalEnd(da_Veloc,local_FV,INSERT_VALUES,Veloc_Cond);CHKERRQ(ierr);
	}
	
	
	
	if (rank==0) printf("Init_veloc\n");
	Init_Veloc();
	if (rank==0) printf("Init_veloc: fim\n");
	
	char nome[100];
	PetscViewer viewer;
	sprintf(nome,"Init_veloc.txt");
	
	PetscViewerASCIIOpen(PETSC_COMM_WORLD,nome,&viewer);
	VecView(Veloc,viewer);
	PetscViewerDestroy(&viewer);
	
	int ind;
	PetscReal r;
	
	VecMax(Veloc_Cond,&ind,&r);
	
	if (rank==0) printf("VecMax Veloc_Cond: %f\n",r);
	
	VecSum(Veloc_Cond,&r);
	
	if (rank==0) printf("VecSum Veloc_Cond: %f\n",r);
	
	
	
	ierr = KSPCreate(PETSC_COMM_WORLD,&V_ksp);CHKERRQ(ierr);
  ierr = KSPSetDM(V_ksp,da_Veloc);CHKERRQ(ierr);
  ierr = KSPSetDMActive(V_ksp,PETSC_FALSE);CHKERRQ(ierr);
	ierr = KSPSetOptionsPrefix(V_ksp,"veloc_"); CHKERRQ(ierr);
	
	PetscTime(&Tempo2p);
	if (rank==0) printf("Velocity create : %lf\n",Tempo2p-Tempo1p);
	
	montaKeVeloc_general(Ke_veloc_general,dx_const,dz_const);
	
	montaCeVeloc(VCe);
	montafeVeloc(VfMe);
	
	
	PetscFunctionReturn(0);
	 
}

PetscErrorCode build_veloc_3d()
{
	
	PetscErrorCode ierr;
	
	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	
	PetscLogDouble Tempo1p,Tempo2p;
	
	PetscTime(&Tempo1p);
	if (rank==0) printf("VA,VB,Vf -> zero entries\n");
	ierr = MatZeroEntries(VA);CHKERRQ(ierr);
	if (rank==0) printf("passou VA\n");
	ierr = MatZeroEntries(VB);CHKERRQ(ierr);
	ierr = VecZeroEntries(Vf);CHKERRQ(ierr);
	ierr = VecZeroEntries(Vf_P);CHKERRQ(ierr);
	
	ierr = VecZeroEntries(Precon);CHKERRQ(ierr);
	if (Verif_first_veloc==1)	VecCopy(Veloc_fut,Veloc_weight);
	Verif_first_veloc=1;

	if (PRESSURE_INIT==0){
		PRESSURE_INIT=1;
		ierr = calc_pressure();
		ierr = shift_pressure();
		write_pressure(-1);
	}

	ierr = moveSwarm(0.0);
	ierr = Swarm2Mesh();
	
	if (rank==0) printf("build VA,Vf\n");
	ierr = AssembleA_Veloc(VA,VG,da_Veloc,da_Thermal);CHKERRQ(ierr);
	//if (rank==0) printf("t\n");
	
	ierr = VecReciprocal(Precon);
	
	ierr = calc_drho();CHKERRQ(ierr);
	
	
	ierr = AssembleF_Veloc(Vf,da_Veloc,da_Thermal,Vf_P);CHKERRQ(ierr);
	//if (rank==0) printf("t\n");
	
	
	PetscTime(&Tempo2p);
	if (rank==0) printf("Velocity build : %lf\n",Tempo2p-Tempo1p);
	
	PetscFunctionReturn(0);
	
}




PetscErrorCode solve_veloc_3d()
{
	PetscErrorCode ierr;
	PetscLogDouble Tempo1,Tempo2;

	PC V_pc;
	
	int rank;
	
	
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	
	PetscTime(&Tempo1);
	
	
	//VecCopy(Veloc_fut,Veloc);tem que vir antes?? está no veloc_total agora!!!ok
	
	/* SOLVE */
	
	//if (rank==0) printf("k\n");
	ierr = KSPSetOperators(V_ksp,VA,VA);CHKERRQ(ierr);
	if (direct_solver==1){
		ierr = KSPGetPC(V_ksp,&V_pc);CHKERRQ(ierr);
		ierr = PCSetType(V_pc,PCLU);CHKERRQ(ierr);
		PCFactorSetMatSolverType(V_pc,MATSOLVERMUMPS);
	}

	//if (rank==0) printf("k\n");
	ierr = KSPSetFromOptions(V_ksp);CHKERRQ(ierr);
	//if (rank==0) printf("k\n");
	ierr = KSPSetInitialGuessNonzero(V_ksp,PETSC_TRUE);
	
	ierr = KSPSetTolerances(V_ksp,rtol,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	
	
	////////
	
	PetscReal denok=0,betak,alphak;
	
	PetscInt maxk = 400,k;
	
	ierr = KSPSolve(V_ksp,Vf_P,Veloc_fut);CHKERRQ(ierr); /// Solve KV0 = F
	
	PetscInt its;
	
	ierr = KSPGetIterationNumber(V_ksp,&its);CHKERRQ(ierr);
	
	VecPointwiseMult(Veloc_fut,Veloc_fut,Veloc_Cond); ///!!!ok ... apenas zerando nas b.c.
	VecAXPY(Veloc_fut,1.0,Veloc_0); ///!!!ok apenas colocando os valores de Veloc_0 em Veloc nas b.c.
	
	//write_veloc_3d(101);
	
	
	ierr = MatMultTranspose(VG,Veloc_fut,rk_vec2);CHKERRQ(ierr); /// r0 = G^T V0
	
	
	ierr = VecPointwiseMult(zk_vec2,Precon,rk_vec2); /// z0 = Precon*r0 = M^-1 *r0 <----
	
	ierr = VecDot(rk_vec2,rk_vec2,&denok);CHKERRQ(ierr); //denok = r0^2 ?
	
	if (rank==0) printf("denok = %lg, its = %d\n",denok, its);
	
	
	
	for (k=1;k<maxk && denok>denok_min;k++){ /// while denok>denok_min ?
		if (k==1) VecCopy(zk_vec2,sk_vec); ///if k=1: s1 = zk <----
		else {
			VecCopy(zk_vec2,zk_vec);
			
			ierr = VecPointwiseMult(zk_vec2,Precon,rk_vec2); /// z_(k-1) = Precon*r_(k-1) <----
			
			VecDot(zk_vec2,rk_vec2,&betak);
			VecDot(zk_vec,rk_vec,&denok);
			betak=betak/denok; /// bk = (z_(k-1)*r_(k-1))/(z_(k-2)*r_(k-2)) <----
			VecAYPX(sk_vec,betak,zk_vec2); /// sk = z_(k-1) + b*s_(k-1) <----
		}
		ierr = MatMult(VG,sk_vec,gs_vec);
		
		VecPointwiseMult(gs_vec,gs_vec,Veloc_Cond);
		
		VecDot(gs_vec,gs_vec,&denok);
		
		
		
		KSPSolve(V_ksp,gs_vec,uk_vec); /// K uk = G sk
		VecPointwiseMult(uk_vec,uk_vec,Veloc_Cond);//!!!ok adicionado agora para zerar as condicoes de contorno
		
		ierr = KSPGetIterationNumber(V_ksp,&its);CHKERRQ(ierr);
		
		VecDot(zk_vec2,rk_vec2,&alphak);
		VecDot(gs_vec,uk_vec,&denok);
		
		alphak=alphak/denok; /// ak = (r_(k-1) z_(k-1))/((G sk)uk) <----
		
		VecAXPY(Pressure,alphak,sk_vec); /// Pk = P_(k-1) + ak*sk
		
		VecAXPY(Veloc_fut,-alphak,uk_vec); /// Vk = V_(k-1) - ak*uk
		
		VecCopy(rk_vec2,rk_vec);
		
		ierr = MatMultTranspose(VG,uk_vec,rk_vec2);CHKERRQ(ierr);
		
		VecAYPX(rk_vec2,-alphak,rk_vec); ///rk = r_(k-1) - ak*(G uk)
		
		VecDot(rk_vec2,rk_vec2,&denok);
		
		if (rank==0) printf("denok = %lg, k=%d, its = %d\n",denok,k,its);
		
	}

	shift_pressure();
	
	////////
	
	
	PetscTime(&Tempo2);
	if (rank==0) printf("Velocity solve: %lf\n",Tempo2-Tempo1);
	
	PetscFunctionReturn(0);
	
}


PetscErrorCode destroy_veloc_3d()
{
	
	int rank;
	
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	PetscLogDouble Tempo1,Tempo2;
	PetscTime(&Tempo1);
	
	PetscErrorCode ierr;
	ierr = KSPDestroy(&V_ksp);CHKERRQ(ierr);
	ierr = VecDestroy(&Veloc_weight);CHKERRQ(ierr);
	ierr = VecDestroy(&Veloc_fut);CHKERRQ(ierr);
	ierr = VecDestroy(&Veloc);CHKERRQ(ierr);
	ierr = VecDestroy(&Vf);CHKERRQ(ierr);
	ierr = VecDestroy(&Vf_P);CHKERRQ(ierr);
	ierr = MatDestroy(&VA);CHKERRQ(ierr);
	ierr = MatDestroy(&VB);CHKERRQ(ierr);
	ierr = DMDestroy(&da_Veloc);CHKERRQ(ierr);
	
	PetscTime(&Tempo2);
	if (rank==0) printf("Velocity destroy: %lf\n",Tempo2-Tempo1);
	
	
	PetscFunctionReturn(0);
}

PetscErrorCode write_veloc_3d(int cont)
{
	int rank;
	
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	PetscLogDouble Tempo1,Tempo2;
	PetscTime(&Tempo1);
	
	PetscViewer viewer;
	
	char nome[100];
	
	sprintf(nome,"Veloc_fut_%d.txt",cont);
	
	PetscViewerASCIIOpen(PETSC_COMM_WORLD,nome,&viewer);
	VecView(Veloc_fut,viewer);
	PetscViewerDestroy(&viewer);
	
	PetscTime(&Tempo2);
	if (rank==0) printf("Velocity write: %lf\n",Tempo2-Tempo1);
	
	PetscFunctionReturn(0);
}

PetscErrorCode write_veloc_cond(int cont)
{
	int rank;
	
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	PetscLogDouble Tempo1,Tempo2;
	PetscTime(&Tempo1);
	
	PetscViewer viewer;
	
	char nome[100];
	
	sprintf(nome,"Veloc_Cond_%d.txt",cont);
	
	PetscViewerASCIIOpen(PETSC_COMM_WORLD,nome,&viewer);
	VecView(Veloc_Cond,viewer);
	PetscViewerDestroy(&viewer);
	
	PetscTime(&Tempo2);
	if (rank==0) printf("Velocity Cond write: %lf\n",Tempo2-Tempo1);
	
	PetscFunctionReturn(0);
}










PetscErrorCode montaKeVeloc_general(PetscReal *KeG, double dx_const, double dz_const){
	
	long i,j,ii;
	long aux;
	
	
	
	double kx,kz;
	
	double ex,ez;
	
	long cont;
	
	double N_x[V_NE];
	double N_z[V_NE];
	
	double SN[3][V_GT];
	
	for (i=0;i<3;i++){
		for (j=0;j<V_GT;j++){
			SN[i][j]=0;
		}
	}
	
	
	long point;
	
	for (point=0;point<GaussQuad;point++){
		for (i=0;i<V_GT;i++){
			for (j=0;j<V_GT;j++){
				KeG[(i*V_GT+j)+point*V_GT*V_GT]=0.0;
			}
		}
	}
	
	double Hx,Hz,prodH;
	
	
	point=0;
	for (kz=-r06; kz<=r06; kz+=r06){
		if (kz==0) Hz=r8p9;
		else Hz=r5p9;
		for (kx=-r06; kx<=r06; kx+=r06){
			if (kx==0) Hx=r8p9;
			else Hx=r5p9;
			prodH = Hx*Hz;
			cont=0;
			for (ez=-1.;ez<=1.;ez+=2.){
				for (ex=-1.;ex<=1.;ex+=2.){
					//N[cont]=(1+ex*kx)*(1+ez*kz)/4.0;
					N_x[cont]=ex*(1+ez*kz)/2.0/dx_const; //!!! 2d
					N_z[cont]=(1+ex*kx)*ez/2.0/dz_const; //!!! 2d
					
					cont++;
				}
			}
			
			
			for (j=0;j<V_NE;j++){
				aux = j*V_GN;
				SN[0][aux  ]=N_x[j];
				SN[1][aux+1]=N_z[j];
				
				SN[2][aux  ]=N_z[j];SN[2][aux+1]=N_x[j];
				
				//Nd[aux  ]=N_x[j];
				//Nd[aux+1]=N_y[j];
				//Nd[aux+2]=N_z[j];
			}
			
			
			for (i=0;i<V_GT;i++){
				for (j=0;j<V_GT;j++){
					for (ii=0;ii<2;ii++){
						KeG[(i*V_GT+j)+point*V_GT*V_GT]+=prodH*2*SN[ii][i]*SN[ii][j];
					}
					for (;ii<3;ii++){
						KeG[(i*V_GT+j)+point*V_GT*V_GT]+=prodH*SN[ii][i]*SN[ii][j];
					}
					
				}
			}
			
			point++;
			
		}
	
	}
	
	/*for (i=0;i<V_NE*GaussQuad;i++){
		PetscPrintf(PETSC_COMM_WORLD,"%lg %lg %lg\n",N_x_Gauss[i],N_y_Gauss[i],N_z_Gauss[i]);
	}
	exit(1);*/
	
	
	PetscFunctionReturn(0);
	
}



PetscReal montaKeVeloc_simplif(PetscReal *Ke,PetscReal *KeG, PetscReal *geoq_ele){
	
	long i,j;
	
	double Visc_local,Geoq_local;
	double visc_meio;
	
	//PetscErrorCode ierr=0;
	
	double kx,kz;
	
	double ex,ez;
	
	long cont;
	
	
	
	//PetscReal Visc_ele[V_NE];
	
	//for (i=0;i<V_NE;i++) Visc_ele[i]=visco_const;
	
	//VecGetValues(visc_vec,V_NE,&Hexa_thermal[t*V_NE],Visc_ele);
	
	//for (i=0;i<V_NE;i++) printf("%g ",Visc_ele[i]);
	//printf("\n");
	
	double Hx,Hz,prodH;
	
	long point=0;//,cont_p=0;
	
	//PetscReal strain[6];
	
	for (i=0;i<V_GT*V_GT;i++) Ke[i]=0.0;
	
	if (visc_const_per_element==0){
		
		for (kz=-r06; kz<=r06; kz+=r06){
			if (kz==0) Hz=r8p9;
			else Hz=r5p9;
			
			for (kx=-r06; kx<=r06; kx+=r06){
				if (kx==0) Hx=r8p9;
				else Hx=r5p9;
				
				Geoq_local = 0.0;
				
				prodH = Hx*Hz;
				cont=0;
				for (ez=-1.;ez<=1.;ez+=2.){
					for (ex=-1.;ex<=1.;ex+=2.){
						Geoq_local+=geoq_ele[cont]*(1+ex*kx)*(1+ez*kz)/4.0;
						cont++;
					}
				}
				
				Visc_local = Geoq_local;
				
				if (kx==0 && kz==0) visc_meio = Visc_local;
				
				
				if (Visc_local<visc_aux_MIN) visc_aux_MIN=Visc_local;
				if (Visc_local>visc_aux_MAX) visc_aux_MAX=Visc_local;
				
				
				for (i=0;i<V_GT;i++){
					for (j=0;j<V_GT;j++){
						Ke[i*V_GT+j]+=KeG[(i*V_GT+j)+point]*Visc_local;
					}
				}
				
				point+=V_GT*V_GT;
			}
		
		}
	}
	else {
		double Visc_mean;

		Geoq_local = geoq_ele[0];
				
		Visc_local = Geoq_local;
		
		Visc_mean = Visc_local;
		
		
		if (Visc_local<visc_aux_MIN) visc_aux_MIN=Visc_local;
		if (Visc_local>visc_aux_MAX) visc_aux_MAX=Visc_local;
		
		for (kz=-r06; kz<=r06; kz+=r06){
			if (kz==0) Hz=r8p9;
			else Hz=r5p9;
			
			for (kx=-r06; kx<=r06; kx+=r06){
				if (kx==0) Hx=r8p9;
				else Hx=r5p9;
				
				
				for (i=0;i<V_GT;i++){
					for (j=0;j<V_GT;j++){
						Ke[i*V_GT+j]+=KeG[(i*V_GT+j)+point];
					}
				}
				
				point+=V_GT*V_GT;
			}
		
		}
		visc_meio = Visc_mean;
		
		for (i=0;i<V_GT;i++){
			for (j=0;j<V_GT;j++){
				Ke[i*V_GT+j]*=Visc_mean;
			}
		}
	}
	
	
	
	//printf("Visc_min = %lg; Visc_max = %lg\n",Visc_min,Visc_max);
	
	return(visc_meio);

}


PetscErrorCode montafeVeloc(PetscReal *fMe)
{
	long i,j;
	
	double kx,kz;
	
	double ex,ez;
	
	long cont;
	
	double N[V_NE];
	
	for (i=0;i<V_GT;i++){
		for (j=0;j<V_NE;j++){
			fMe[i*V_NE+j]=0.0;
		}
	}
	
	double Hx,Hz,prodH;
	for (kz=-r06; kz<=r06; kz+=r06){
		if (kz==0) Hz=r8p9;
		else Hz=r5p9;
		
		for (kx=-r06; kx<=r06; kx+=r06){
			if (kx==0) Hx=r8p9;
			else Hx=r5p9;
			
			
			prodH = Hx*Hz;
			cont=0;
			for (ez=-1.;ez<=1.;ez+=2.){
				for (ex=-1.;ex<=1.;ex+=2.){
					N[cont]=(1+ex*kx)*(1+ez*kz)/4.0;
					cont++;
				}
			}
			
			
			for (i=0;i<V_NE;i++){
				for (j=0;j<V_NE;j++){
					fMe[(i*V_GN+1)*V_NE+j]+=prodH*N[i]*N[j]; ///!!!! em 2d é +1
				}
			}
			
		}
	
	}
	
	for (i=0;i<V_GT*V_NE;i++){
		fMe[i]*=dx_const*dz_const;
	}

	
	PetscFunctionReturn(0);
}


PetscErrorCode montaCeVeloc(PetscReal *Ce){
	
	long i;
	
	/*double dx,dy,dz;
	 
	 //
	 
	 dx = xyz_thermal[Hexa_thermal[t][2]][0]-xyz_thermal[Hexa_thermal[t][0]][0];
	 dy = xyz_thermal[Hexa_thermal[t][1]][1]-xyz_thermal[Hexa_thermal[t][0]][1];
	 dz = xyz_thermal[Hexa_thermal[t][4]][2]-xyz_thermal[Hexa_thermal[t][0]][2];
	 
	 //*/
	
	double kx,kz;
	
	double ex,ez;
	
	long cont;
	
	double N_x[V_NE];
	double N_z[V_NE];
	
	for (i=0;i<V_GT;i++){
		Ce[i]=0.0;
	}
	
	double Hx,Hz,prodH;
	
	for (kz=-r06; kz<=r06; kz+=r06){
		if (kz==0) Hz=r8p9;
		else Hz=r5p9;
		
		for (kx=-r06; kx<=r06; kx+=r06){
			if (kx==0) Hx=r8p9;
			else Hx=r5p9;
			
			
			prodH = Hx*Hz;
			cont=0;
			for (ez=-1.;ez<=1.;ez+=2.){
				for (ex=-1.;ex<=1.;ex+=2.){
					//N[cont]=(1+ex*kx)*(1+ey*ky)*(1+ez*kz)/4.0;
					N_x[cont]=ex*(1+ez*kz)/2.0/dx_const;
					N_z[cont]=(1+ex*kx)*ez/2.0/dz_const;
					cont++;
				}
			}
			
			for (i=0;i<V_NE;i++){
				Ce[i*V_GN]+=prodH*N_x[i]; // 1/8 funcao de forma da pressao (constante no elemento) isso no 3D
				Ce[i*V_GN+1]+=prodH*N_z[i];
			}
			
		}
		
	}
	
	///adicionado
	for (i=0;i<V_GT;i++){
		Ce[i]*=-dx_const*dz_const;
	}
	
	PetscFunctionReturn(0);

	
}
