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
} Stokes2d;

typedef struct {
	PetscScalar u;
	PetscScalar v;
	PetscScalar w;
} Stokes3d;

PetscErrorCode ascii2bin(char *s1, char *s2);

PetscErrorCode AssembleA_Veloc_2d(Mat A,Mat AG,DM veloc_da, DM temper_da);
PetscErrorCode AssembleF_Veloc_2d(Vec F,DM veloc_da,DM drho_da, Vec FP);

PetscErrorCode AssembleA_Veloc_3d(Mat A,Mat AG,DM veloc_da, DM temper_da);
PetscErrorCode AssembleF_Veloc_3d(Vec F,DM veloc_da,DM drho_da,Vec FP);

PetscErrorCode montaKeVeloc_general_2d(PetscReal *KeG, double dx_const, double dz_const);
PetscErrorCode montaKeVeloc_general_3d(PetscReal *KeG, double dx_const, double dy_const, double dz_const);

PetscReal montaKeVeloc_simplif_2d(PetscReal *Ke,PetscReal *KeG, PetscReal *geoq_ele);

PetscErrorCode montaCeVeloc2d(PetscReal *Ce);
PetscErrorCode montafeVeloc2d(PetscReal *fMe);

PetscErrorCode montaCeVeloc3d(PetscReal *Ce);
PetscErrorCode montafeVeloc3d(PetscReal *fMe);

PetscErrorCode calc_drho();

PetscErrorCode calc_pressure_2d();
PetscErrorCode shift_pressure_2d();

PetscErrorCode calc_pressure_3d();
PetscErrorCode shift_pressure_3d();

PetscErrorCode write_veloc(int cont, PetscInt binary_out);

PetscErrorCode write_veloc_cond(int cont, PetscInt binary_out);

PetscErrorCode write_pressure(int cont, PetscInt binary_out);

PetscErrorCode Init_Veloc(int dimensions);

PetscErrorCode write_all_(int cont,Vec u, char *variable_name, PetscInt binary_out);

PetscErrorCode moveSwarm(int dimensions, PetscReal dt);
PetscErrorCode Swarm2Mesh_2d();
PetscErrorCode Swarm2Mesh_3d();

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
extern double dy_const;
extern double dz_const;

extern PetscReal *N_x_Gauss;
extern PetscReal *N_y_Gauss;
extern PetscReal *N_z_Gauss;

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

extern PetscInt periodic_boundary;

extern int n_interfaces;

extern PetscInt binary_output;

PetscErrorCode create_veloc(int dimensions, PetscInt mx, PetscInt my, PetscInt mz, PetscInt Px, PetscInt Py, PetscInt Pz)
{
	PetscInt       dof,stencil_width;
	DMBoundaryType boundary_type;
	PetscErrorCode ierr;

	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

	PetscLogDouble Tempo1p,Tempo2p;

	PetscTime(&Tempo1p);

	PetscFunctionBeginUser;

	if (periodic_boundary == 0) {
		stencil_width = 1;
		boundary_type = DM_BOUNDARY_NONE;
	} else {
		stencil_width = 2;
		boundary_type = DM_BOUNDARY_PERIODIC;
	}

	if (dimensions == 2) {
		dof = 2;
		ierr = DMDACreate2d(PETSC_COMM_WORLD, boundary_type, DM_BOUNDARY_NONE, DMDA_STENCIL_BOX,
							mx+1, mz+1, Px, Pz, dof, stencil_width, NULL, NULL, &da_Veloc); CHKERRQ(ierr);
	} else {
		if (boundary_type == DM_BOUNDARY_PERIODIC) {
			//
			// NOTE: periodic_boundary is currently not used in 3D
			//
			ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\nWARNING: periodic_boundary is currently not supported in 3D mode\n\n");
		}

		dof = 3;
		ierr = DMDACreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_BOX,
							mx+1, my+1, mz+1, Px, Py, Pz, dof, stencil_width, NULL, NULL, NULL, &da_Veloc); CHKERRQ(ierr);
	}

	ierr = DMSetFromOptions(da_Veloc);CHKERRQ(ierr);

	ierr = DMSetUp(da_Veloc);CHKERRQ(ierr);

	ierr = DMDASetFieldName(da_Veloc, 0, "V_x"); CHKERRQ(ierr);
	if (dimensions == 2) {
		ierr = DMDASetFieldName(da_Veloc, 1, "V_z"); CHKERRQ(ierr);
	} else {
		ierr = DMDASetFieldName(da_Veloc, 1, "V_y"); CHKERRQ(ierr);
		ierr = DMDASetFieldName(da_Veloc, 2, "V_z"); CHKERRQ(ierr);
	}

	ierr = PetscCalloc1(V_GT*V_GT,&Ke_veloc); CHKERRQ(ierr);
	ierr = PetscCalloc1(V_GT*V_GT,&Ke_veloc_final); CHKERRQ(ierr);

	ierr = PetscCalloc1(V_GT*V_GT*GaussQuad,&Ke_veloc_general); CHKERRQ(ierr);

	if (dimensions == 3) {
		ierr = PetscCalloc1(GaussQuad*V_NE, &N_x_Gauss); CHKERRQ(ierr);
		ierr = PetscCalloc1(GaussQuad*V_NE, &N_y_Gauss); CHKERRQ(ierr);
		ierr = PetscCalloc1(GaussQuad*V_NE, &N_z_Gauss); CHKERRQ(ierr);
	}

	ierr = PetscCalloc1(V_GT,&Vfe); CHKERRQ(ierr);

	ierr = PetscCalloc1(V_GT*V_NE,&VfMe); CHKERRQ(ierr);
	ierr = PetscCalloc1(V_GT,&VCe); CHKERRQ(ierr);

	if (dimensions == 2) {
		ierr = DMDASetUniformCoordinates(da_Veloc, 0.0, Lx, -depth, 0.0, 0.0, 0.0);CHKERRQ(ierr);
	} else {
		ierr = DMDASetUniformCoordinates(da_Veloc, 0.0, Lx, 0.0, Ly, -depth, 0.0);CHKERRQ(ierr);
	}

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


	Stokes2d **ff2d;
	Stokes3d ***ff3d;
	PetscInt M, N, P;
	PetscInt sx, sy, sz, mmx, mmy, mmz;
	PetscInt i, j, k, t;


	ierr = VecZeroEntries(local_FV); CHKERRQ(ierr);

	if (dimensions == 2) {
		ierr = DMDAGetInfo(da_Veloc, 0, &M, &P, NULL, 0, 0, 0,  0, 0, 0, 0, 0, 0); CHKERRQ(ierr);
		ierr = DMDAGetCorners(da_Veloc, &sx, &sz, NULL, &mmx, &mmz, NULL); CHKERRQ(ierr);
		ierr = DMDAVecGetArray(da_Veloc, local_FV, &ff2d); CHKERRQ(ierr);
	} else {
		ierr = DMDAGetInfo(da_Veloc, 0, &M, &N, &P, 0, 0, 0,  0, 0, 0, 0, 0, 0); CHKERRQ(ierr);
		ierr = DMDAGetCorners(da_Veloc, &sx, &sy, &sz, &mmx, &mmy, &mmz); CHKERRQ(ierr);
		ierr = DMDAVecGetArray(da_Veloc, local_FV, &ff3d); CHKERRQ(ierr);
	}


	PetscInt ix[1];
	PetscScalar y[1];

	PetscInt low,high;

	if (bcv_extern == 1) {
		char s1[100],s2[100];

		sprintf(s1,"bcv_0.txt");
		sprintf(s2,"bcv_init.bin");

		if (rank==0){
			ierr = ascii2bin(s1,s2); CHKERRQ(ierr);
		}
		MPI_Barrier(PETSC_COMM_WORLD);


		PetscInt size0;
		PetscViewer    viewer;

		VecGetSize(Veloc_Cond,&size0);


		Vec Fprov;


		PetscViewerBinaryOpen(PETSC_COMM_WORLD,s2,FILE_MODE_READ,&viewer);
		VecCreate(PETSC_COMM_WORLD,&Fprov);
		VecLoad(Fprov,viewer);
		PetscViewerDestroy(&viewer);


		VecGetOwnershipRange(Fprov,&low,&high);

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

	}
	else {
		if (dimensions == 2) {
			for (k=sz; k<sz+mmz; k++) {
				for (i=sx; i<sx+mmx; i++) {
					ff2d[k][i].u = 1.0;
					ff2d[k][i].w = 1.0;

					if (periodic_boundary==0){
						if (i==0   && bcv_left_normal==1) ff2d[k][i].u = 0.0;
						if (i==0   && bcv_left_slip==1) {
							ff2d[k][i].w = 0.0;
						}

						if (i==M-1 && bcv_right_normal==1)ff2d[k][i].u = 0.0;
						if (i==M-1   && bcv_right_slip==1) {
							ff2d[k][i].w = 0.0;
						}
					}

					if (k==0   && bcv_bot_normal==1) ff2d[k][i].w = 0.0;
					if (k==0   && bcv_bot_slip==1){
						ff2d[k][i].u = 0.0;
					}

					if (k==P-1 && bcv_top_normal==1) ff2d[k][i].w = 0.0;
					if (k==P-1 && bcv_top_slip==1){
						ff2d[k][i].u = 0.0;
					}

				}
			}
			ierr = DMDAVecRestoreArray(da_Veloc,local_FV,&ff2d);CHKERRQ(ierr);
		} else {
			for (k=sz; k<sz+mmz; k++) {
				for (j=sy; j<sy+mmy; j++) {
					for (i=sx; i<sx+mmx; i++) {
						ff3d[k][j][i].u = 1.0;
						ff3d[k][j][i].v = 1.0;
						ff3d[k][j][i].w = 1.0;

						if (i==0   && bcv_left_normal==1) ff3d[k][j][i].u = 0.0;
						if (i==0   && bcv_left_slip==1) {
							ff3d[k][j][i].v = 0.0;
							ff3d[k][j][i].w = 0.0;
						}

						if (i==M-1 && bcv_right_normal==1)ff3d[k][j][i].u = 0.0;
						if (i==M-1   && bcv_right_slip==1) {
							ff3d[k][j][i].v = 0.0;
							ff3d[k][j][i].w = 0.0;
						}

						if (j==0 || j==N-1) ff3d[k][j][i].v = 0.0;

						if (k==0   && bcv_bot_normal==1) ff3d[k][j][i].w = 0.0;
						if (k==0   && bcv_bot_slip==1){
							ff3d[k][j][i].u = 0.0;
							ff3d[k][j][i].v = 0.0;
						}

						if (k==P-1 && bcv_top_normal==1) ff3d[k][j][i].w = 0.0;
						if (k==P-1 && bcv_top_slip==1){
							ff3d[k][j][i].u = 0.0;
							ff3d[k][j][i].v = 0.0;
						}

					}
				}
			}
			ierr = DMDAVecRestoreArray(da_Veloc, local_FV, &ff3d);CHKERRQ(ierr);
		}

		ierr = DMLocalToGlobalBegin(da_Veloc, local_FV, INSERT_VALUES, Veloc_Cond); CHKERRQ(ierr);
		ierr = DMLocalToGlobalEnd(da_Veloc, local_FV, INSERT_VALUES, Veloc_Cond); CHKERRQ(ierr);
	}

	Init_Veloc(dimensions);

	int ind;
	PetscReal r;

	VecMax(Veloc_Cond,&ind,&r);

	VecSum(Veloc_Cond,&r);

	ierr = KSPCreate(PETSC_COMM_WORLD,&V_ksp);CHKERRQ(ierr);
	ierr = KSPSetDM(V_ksp,da_Veloc);CHKERRQ(ierr); // NOTE: check for 3d
	ierr = KSPSetDMActive(V_ksp,PETSC_FALSE);CHKERRQ(ierr); // NOTE: check for 3d
	ierr = KSPSetOptionsPrefix(V_ksp,"veloc_"); CHKERRQ(ierr);

	PetscTime(&Tempo2p);
	PetscPrintf(PETSC_COMM_WORLD, "Velocity field (creation): %lf s\n",Tempo2p-Tempo1p);

	if (dimensions == 2) {
		montaKeVeloc_general_2d(Ke_veloc_general, dx_const, dz_const);
		montaCeVeloc2d(VCe);
		montafeVeloc2d(VfMe);
	} else {
		montaKeVeloc_general_3d(Ke_veloc_general, dx_const, dy_const, dz_const);
		montaCeVeloc3d(VCe);
		montafeVeloc3d(VfMe);
	}

	PetscFunctionReturn(0);
}

PetscErrorCode build_veloc(int dimensions)
{
	PetscErrorCode ierr;
	PetscLogDouble Tempo1p, Tempo2p;

	PetscTime(&Tempo1p);

	ierr = MatZeroEntries(VA);CHKERRQ(ierr);

	ierr = MatZeroEntries(VB);CHKERRQ(ierr);
	ierr = VecZeroEntries(Vf);CHKERRQ(ierr);
	ierr = VecZeroEntries(Vf_P);CHKERRQ(ierr);

	ierr = VecZeroEntries(Precon);CHKERRQ(ierr);
	if (Verif_first_veloc==1)	VecCopy(Veloc_fut,Veloc_weight);
	Verif_first_veloc=1;

	if (PRESSURE_INIT == 0 && n_interfaces > 0) {
		PRESSURE_INIT = 1;

		if (dimensions == 2) {
			ierr = calc_pressure_2d();
			ierr = shift_pressure_2d();
		} else {
			ierr = calc_pressure_3d();
			ierr = shift_pressure_3d();
		}

		write_pressure(-1,binary_output);
	}

	ierr = moveSwarm(dimensions, 0.0);

	if (dimensions == 2) {
		ierr = Swarm2Mesh_2d();

		ierr = AssembleA_Veloc_2d(VA,VG,da_Veloc,da_Thermal); CHKERRQ(ierr);
	} else {
		ierr = Swarm2Mesh_3d();

		ierr = AssembleA_Veloc_3d(VA,VG,da_Veloc,da_Thermal); CHKERRQ(ierr);
	}

	ierr = VecReciprocal(Precon);

	ierr = calc_drho();CHKERRQ(ierr);

	if (dimensions == 2) {
		ierr = AssembleF_Veloc_2d(Vf,da_Veloc,da_Thermal,Vf_P); CHKERRQ(ierr);
	} else {
		ierr = AssembleF_Veloc_3d(Vf,da_Veloc,da_Thermal,Vf_P); CHKERRQ(ierr);
	}

	PetscTime(&Tempo2p);
	PetscPrintf(PETSC_COMM_WORLD, "  Velocity field (building): %lf s\n",Tempo2p-Tempo1p);

	PetscFunctionReturn(0);
}

PetscErrorCode solve_veloc(int dimensions)
{
	PetscErrorCode ierr;
	PetscLogDouble Tempo1,Tempo2;

	PC V_pc;

	PetscTime(&Tempo1);


	/* SOLVE */

	//PetscPrintf(PETSC_COMM_WORLD, "k\n");
	ierr = KSPSetOperators(V_ksp,VA,VA);CHKERRQ(ierr);
	if (direct_solver==1){
		ierr = KSPGetPC(V_ksp,&V_pc);CHKERRQ(ierr);
		ierr = PCSetType(V_pc,PCLU);CHKERRQ(ierr);
		PCFactorSetMatSolverType(V_pc,MATSOLVERMUMPS);
	}

	//PetscPrintf(PETSC_COMM_WORLD, "k\n");
	ierr = KSPSetFromOptions(V_ksp);CHKERRQ(ierr);
	//PetscPrintf(PETSC_COMM_WORLD, "k\n");
	ierr = KSPSetInitialGuessNonzero(V_ksp,PETSC_TRUE);

	ierr = KSPSetTolerances(V_ksp,rtol,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);


	////////

	PetscReal denok=0,betak,alphak;

	PetscInt maxk = 400,k;

	ierr = KSPSolve(V_ksp,Vf_P,Veloc_fut);CHKERRQ(ierr); /// Solve KV0 = F

	PetscInt its;

	ierr = KSPGetIterationNumber(V_ksp,&its);CHKERRQ(ierr);

	VecPointwiseMult(Veloc_fut,Veloc_fut,Veloc_Cond); ///check: zero values at the b.c.
	VecAXPY(Veloc_fut,1.0,Veloc_0); ///check: applying Veloc_0 in Veloc at the b.c.

	//write_veloc(101);


	ierr = MatMultTranspose(VG,Veloc_fut,rk_vec2);CHKERRQ(ierr); /// r0 = G^T V0


	ierr = VecPointwiseMult(zk_vec2,Precon,rk_vec2); /// z0 = Precon*r0 = M^-1 *r0 <----

	ierr = VecDot(rk_vec2,rk_vec2,&denok);CHKERRQ(ierr); //denok = r0^2 ?

	PetscPrintf(PETSC_COMM_WORLD, "    Uzawa iteration: 0, denok = %lg, its = %d\n",denok, its);



	for (k=1;k<maxk && denok>denok_min;k++){ /// while denok>denok_min
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
		VecPointwiseMult(uk_vec,uk_vec,Veloc_Cond);//zero values at the b.c.

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

		PetscPrintf(PETSC_COMM_WORLD, "    Uzawa iteration: %d, denok = %lg, its = %d\n",k,denok,its);
	}

	if (dimensions == 2) {
		shift_pressure_2d();
	} else {
		shift_pressure_3d();
	}

	PetscTime(&Tempo2);
	PetscPrintf(PETSC_COMM_WORLD, "    Velocity field (solution): %lf s\n",Tempo2-Tempo1);

	PetscFunctionReturn(0);
}


PetscErrorCode destroy_veloc()
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
	PetscPrintf(PETSC_COMM_WORLD, "Velocity field (destroying): %lf s\n",Tempo2-Tempo1);

	PetscFunctionReturn(0);
}

PetscErrorCode write_veloc(int cont, PetscInt binary_out)
{
	char variable_name[100];

	sprintf(variable_name,"velocity");
	write_all_(cont,Veloc_fut,variable_name,binary_out);
	PetscFunctionReturn(0);
}

PetscErrorCode write_veloc_cond(int cont, PetscInt binary_out)
{
	char variable_name[100];

	sprintf(variable_name,"bc_velocity");
	write_all_(cont,Veloc_Cond,variable_name,binary_out);

	PetscFunctionReturn(0);
}

PetscErrorCode montaKeVeloc_general_2d(PetscReal *KeG, double dx_const, double dz_const){
	long i, j, ii;
	long aux;
	double kx, kz;
	double ex, ez;
	long cont;
	double N_x[V_NE];
	double N_z[V_NE];
	double SN[3][V_GT];
	long point;
	double Hx, Hz, prodH;

	for (i=0; i<3; i++){
		for (j=0; j<V_GT; j++){
			SN[i][j] = 0;
		}
	}

	for (point=0; point<GaussQuad; point++){
		for (i=0; i<V_GT; i++){
			for (j=0; j<V_GT; j++){
				KeG[(i*V_GT+j) + point*V_GT*V_GT] = 0.0;
			}
		}
	}

	point = 0;
	for (kz=-r06; kz<=r06; kz+=r06){
		if (kz==0) Hz = r8p9;
		else Hz = r5p9;
		for (kx=-r06; kx<=r06; kx+=r06){
			if (kx==0) Hx = r8p9;
			else Hx = r5p9;
			prodH = Hx*Hz;
			cont=0;
			for (ez=-1.;ez<=1.;ez+=2.){
				for (ex=-1.;ex<=1.;ex+=2.){
					N_x[cont]=ex*(1+ez*kz)/2.0/dx_const;
					N_z[cont]=(1+ex*kx)*ez/2.0/dz_const;

					cont++;
				}
			}

			for (j=0;j<V_NE;j++){
				aux = j*V_GN;
				SN[0][aux  ] = N_x[j];
				SN[1][aux+1] = N_z[j];
				SN[2][aux  ] = N_z[j]; SN[2][aux+1] = N_x[j];
			}

			for (i=0; i<V_GT; i++){
				for (j=0; j<V_GT; j++){
					for (ii=0; ii<2; ii++){
						KeG[(i*V_GT+j) + point*V_GT*V_GT] += prodH*2*SN[ii][i]*SN[ii][j];
					}
					for (; ii<3; ii++){
						KeG[(i*V_GT+j) + point*V_GT*V_GT] += prodH*SN[ii][i]*SN[ii][j];
					}
				}
			}

			point++;
		}
	}

	PetscFunctionReturn(0);
}

PetscErrorCode montaKeVeloc_general_3d(PetscReal *KeG, double dx_const, double dy_const, double dz_const){
	long i, j, ii;
	long aux;
	double kx, ky, kz;
	double ex, ey, ez;
	long cont;
	double N_x[V_NE];
	double N_y[V_NE];
	double N_z[V_NE];
	double SN[6][V_GT];
	long point;
	double Hx, Hy, Hz, prodH;

	for (i=0; i<6; i++){
		for (j=0; j<V_GT; j++){
			SN[i][j]=0;
		}
	}

	for (point=0;point<GaussQuad;point++){
		for (i=0;i<V_GT;i++){
			for (j=0;j<V_GT;j++){
				KeG[(i*V_GT+j)+point*V_GT*V_GT]=0.0;
			}
		}
	}

	point = 0;
	for (kz=-r06; kz<=r06; kz+=r06){
		if (kz==0) Hz=r8p9;
		else Hz = r5p9;
		for (ky=-r06; ky<=r06; ky+=r06){
			if (ky==0) Hy = r8p9;
			else Hy = r5p9;
			for (kx=-r06; kx<=r06; kx+=r06){
				if (kx==0) Hx = r8p9;
				else Hx = r5p9;

				prodH = Hx*Hy*Hz;
				cont = 0;
				for (ez=-1.;ez<=1.;ez+=2.){
					for (ey=-1.;ey<=1.;ey+=2.){
						for (ex=-1.;ex<=1.;ex+=2.){
							N_x[cont]=ex*(1+ey*ky)*(1+ez*kz)/4.0/dx_const;
							N_y[cont]=(1+ex*kx)*ey*(1+ez*kz)/4.0/dy_const;
							N_z[cont]=(1+ex*kx)*(1+ey*ky)*ez/4.0/dz_const;

							N_x_Gauss[cont+point*8] = ex*(1+ey*ky)*(1+ez*kz)/4.0/dx_const;
							N_y_Gauss[cont+point*8] = (1+ex*kx)*ey*(1+ez*kz)/4.0/dy_const;
							N_z_Gauss[cont+point*8] = (1+ex*kx)*(1+ey*ky)*ez/4.0/dz_const;

							cont++;
						}
					}
				}

				for (j=0;j<V_NE;j++){
					aux = j*V_GN;
					SN[0][aux  ] = N_x[j];
					SN[1][aux+1] = N_y[j];
					SN[2][aux+2] = N_z[j];

					SN[3][aux  ] = N_y[j]; SN[3][aux+1] = N_x[j];
					SN[4][aux+1] = N_z[j]; SN[4][aux+2] = N_y[j];
					SN[5][aux  ] = N_z[j]; SN[5][aux+2] = N_x[j];
				}


				for (i=0; i<V_GT; i++){
					for (j=0; j<V_GT; j++){
						for (ii=0; ii<3; ii++){
							KeG[(i*V_GT+j)+point*V_GT*V_GT] += prodH*2*SN[ii][i]*SN[ii][j];
						}
						for (; ii<6; ii++){
							KeG[(i*V_GT+j) + point*V_GT*V_GT] += prodH*SN[ii][i]*SN[ii][j];
						}
					}
				}

				point++;
			}
		}
	}

	PetscFunctionReturn(0);
}


PetscReal montaKeVeloc_simplif_2d(PetscReal *Ke,PetscReal *KeG, PetscReal *geoq_ele){

	long i,j;

	double Visc_local,Geoq_local;
	double visc_meio=-1.0;

	double kx,kz;

	double ex,ez;

	long cont;



	long point=0;

	for (i=0;i<V_GT*V_GT;i++) Ke[i]=0.0;

	if (visc_const_per_element==0){

		for (kz=-r06; kz<=r06; kz+=r06){

			for (kx=-r06; kx<=r06; kx+=r06){

				Geoq_local = 0.0;

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

			for (kx=-r06; kx<=r06; kx+=r06){

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


	return(visc_meio);

}

PetscReal montaKeVeloc_simplif_3d(PetscReal *Ke,PetscReal *KeG, PetscReal *geoq_ele){

	long i,j;

	double Visc_local,Geoq_local;
	double visc_meio;

	//PetscErrorCode ierr=0;

	double kx,ky,kz;

	double ex,ey,ez;

	long cont;



	long point=0;


	for (i=0;i<V_GT*V_GT;i++) Ke[i]=0.0;

	if (visc_const_per_element==0){

		for (kz=-r06; kz<=r06; kz+=r06){

			for (ky=-r06; ky<=r06; ky+=r06){

				for (kx=-r06; kx<=r06; kx+=r06){

					Geoq_local = 0.0;

					cont=0;
					for (ez=-1.;ez<=1.;ez+=2.){
						for (ey=-1.;ey<=1.;ey+=2.){
							for (ex=-1.;ex<=1.;ex+=2.){
								Geoq_local+=geoq_ele[cont]*(1+ex*kx)*(1+ey*ky)*(1+ez*kz)/8.0;
								cont++;
							}
						}
					}


					Visc_local = Geoq_local;

					if (kx==0 && ky==0 && kz==0) visc_meio = Visc_local;


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
	}
	else {
		double Visc_mean;

		Geoq_local = geoq_ele[0];

		Visc_local = Geoq_local;

		Visc_mean = Visc_local;


		if (Visc_local<visc_aux_MIN) visc_aux_MIN=Visc_local;
		if (Visc_local>visc_aux_MAX) visc_aux_MAX=Visc_local;

		for (kz=-r06; kz<=r06; kz+=r06){
			for (ky=-r06; ky<=r06; ky+=r06){
				for (kx=-r06; kx<=r06; kx+=r06){


					for (i=0;i<V_GT;i++){
						for (j=0;j<V_GT;j++){
							Ke[i*V_GT+j]+=KeG[(i*V_GT+j)+point];
						}
					}

					point+=V_GT*V_GT;
				}
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

PetscErrorCode montafeVeloc2d(PetscReal *fMe)
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
					fMe[(i*V_GN+1)*V_NE+j]+=prodH*N[i]*N[j]; ///check: in 2d it is +1
				}
			}

		}

	}

	for (i=0;i<V_GT*V_NE;i++){
		fMe[i]*=dx_const*dz_const;
	}


	PetscFunctionReturn(0);
}

PetscErrorCode montafeVeloc3d(PetscReal *fMe)
{
	long i,j;

	double kx,ky,kz;

	double ex,ey,ez;

	long cont;

	double N[V_NE];

	for (i=0;i<V_GT;i++){
		for (j=0;j<V_NE;j++){
			fMe[i*V_NE+j]=0.0;
		}
	}

	double Hx,Hy,Hz,prodH;
	for (kz=-r06; kz<=r06; kz+=r06){
		if (kz==0) Hz=r8p9;
		else Hz=r5p9;
		for (ky=-r06; ky<=r06; ky+=r06){
			if (ky==0) Hy=r8p9;
			else Hy=r5p9;
			for (kx=-r06; kx<=r06; kx+=r06){
				if (kx==0) Hx=r8p9;
				else Hx=r5p9;


				prodH = Hx*Hy*Hz;
				cont=0;
				for (ez=-1.;ez<=1.;ez+=2.){
					for (ey=-1.;ey<=1.;ey+=2.){
						for (ex=-1.;ex<=1.;ex+=2.){
							N[cont]=(1+ex*kx)*(1+ey*ky)*(1+ez*kz)/8.0;
							cont++;
						}
					}
				}


				for (i=0;i<V_NE;i++){
					for (j=0;j<V_NE;j++){
						fMe[(i*V_GN+2)*V_NE+j]+=prodH*N[i]*N[j];
					}
				}



			}
		}
	}

	for (i=0;i<V_GT*V_NE;i++){
		fMe[i]*=dx_const*dy_const*dz_const;
	}


	PetscFunctionReturn(0);
}


PetscErrorCode montaCeVeloc2d(PetscReal *Ce) {

	long i;


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
				Ce[i*V_GN]+=prodH*N_x[i]; // 1/8 shape function for pressure (constant in the element) in 3D
				Ce[i*V_GN+1]+=prodH*N_z[i];
			}

		}

	}

	for (i=0;i<V_GT;i++){
		Ce[i]*=-dx_const*dz_const;
	}

	PetscFunctionReturn(0);


}

PetscErrorCode montaCeVeloc3d(PetscReal *Ce) {

	long i;

	/*double dx,dy,dz;

	 //

	 dx = xyz_thermal[Hexa_thermal[t][2]][0]-xyz_thermal[Hexa_thermal[t][0]][0];
	 dy = xyz_thermal[Hexa_thermal[t][1]][1]-xyz_thermal[Hexa_thermal[t][0]][1];
	 dz = xyz_thermal[Hexa_thermal[t][4]][2]-xyz_thermal[Hexa_thermal[t][0]][2];

	 //*/

	double kx,ky,kz;

	double ex,ey,ez;

	long cont;

	double N_x[V_NE];
	double N_y[V_NE];
	double N_z[V_NE];

	for (i=0;i<V_GT;i++){
		Ce[i]=0.0;
	}

	double Hx,Hy,Hz,prodH;

	for (kz=-r06; kz<=r06; kz+=r06){
		if (kz==0) Hz=r8p9;
		else Hz=r5p9;
		for (ky=-r06; ky<=r06; ky+=r06){
			if (ky==0) Hy=r8p9;
			else Hy=r5p9;
			for (kx=-r06; kx<=r06; kx+=r06){
				if (kx==0) Hx=r8p9;
				else Hx=r5p9;


				prodH = Hx*Hy*Hz;
				cont=0;
				for (ez=-1.;ez<=1.;ez+=2.){
					for (ey=-1.;ey<=1.;ey+=2.){
						for (ex=-1.;ex<=1.;ex+=2.){
							//N[cont]=(1+ex*kx)*(1+ey*ky)*(1+ez*kz)/8.0;
							N_x[cont]=ex*(1+ey*ky)*(1+ez*kz)/4.0/dx_const;
							N_y[cont]=(1+ex*kx)*ey*(1+ez*kz)/4.0/dy_const;
							N_z[cont]=(1+ex*kx)*(1+ey*ky)*ez/4.0/dz_const;
							cont++;
						}
					}
				}

				for (i=0;i<V_NE;i++){
					Ce[i*V_GN]+=prodH*N_x[i]; // 1/8 funcao de forma da pressao (constante no elemento)
					Ce[i*V_GN+1]+=prodH*N_y[i];
					Ce[i*V_GN+2]+=prodH*N_z[i];
				}

			}
		}
	}

	///adicionado
	for (i=0;i<V_GT;i++){
		Ce[i]*=-dx_const*dy_const*dz_const;
	}

	PetscFunctionReturn(0);


}
