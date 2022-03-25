#include <petscksp.h>
#include <petscdmda.h>

#include <petsctime.h>

PetscErrorCode montaKeThermal_general_2d(PetscReal *Ke, PetscReal *Me, PetscReal *Fe);
PetscErrorCode montaKeThermal_general_3d(PetscReal *Ke, PetscReal *Me, PetscReal *Fe);

PetscErrorCode DMDAGetLocalElementSize_2d(DM da,PetscInt *mxl,PetscInt *mzl);

PetscErrorCode AssembleA_Thermal_2d(Mat TA,DM thermal_da,PetscReal *TKe,PetscReal *TMe,PetscReal *TFe,
								 DM veloc_da, Vec Veloc_total);

PetscErrorCode AssembleF_Thermal_2d(Vec F,DM thermal_da,PetscReal *TKe,PetscReal *TMe,PetscReal *TFe,
								 DM veloc_da, Vec Veloc_total);

PetscErrorCode AssembleA_Thermal_3d(Mat A,DM thermal_da,PetscReal *TKe,PetscReal *TMe,PetscReal *TFe,
								 DM veloc_da, Vec Veloc_total);

PetscErrorCode AssembleF_Thermal_3d(Vec F,DM thermal_da,PetscReal *TKe,PetscReal *TMe,PetscReal *TFe,
								 DM veloc_da, Vec Veloc_total);

PetscErrorCode Thermal_init_2d(Vec F,DM thermal_da);

PetscErrorCode Heat_flow_at_the_base();

extern double Lx, Ly, depth;

extern double dz_const;

extern PetscReal *NT;
extern PetscReal *NT_x;
extern PetscReal *NT_y;
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
extern Vec Tf, Temper, Temper_0;

extern Vec Pressure_aux;

extern Vec geoq,local_geoq;
extern Vec geoq_rho,local_geoq_rho;
extern Vec geoq_H,local_geoq_H;
extern Vec geoq_strain,local_geoq_strain;
extern Vec geoq_strain_rate,local_geoq_strain_rate;

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

extern Vec local_P_aux;

extern Vec local_dRho;

extern PetscReal rtol;

extern PetscInt periodic_boundary;

extern PetscScalar basal_heat;
extern double kappa;
extern double RHOM;
extern double c_heat_capacity;

PetscErrorCode create_thermal(int dimensions, PetscInt mx, PetscInt my, PetscInt mz, PetscInt Px, PetscInt Py, PetscInt Pz)
{

	PetscInt       dof,stencil_width;
	DMBoundaryType boundary_type;
	PetscErrorCode ierr;



	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

	PetscLogDouble Tempo1p,Tempo2p;

	PetscTime(&Tempo1p);

	dof = 1;

	if (periodic_boundary == 0) {
		stencil_width = 1;
		boundary_type = DM_BOUNDARY_NONE;
	} else {
		stencil_width = 2;
		boundary_type = DM_BOUNDARY_PERIODIC;
	}

	if (dimensions == 2) {
		ierr = DMDACreate2d(PETSC_COMM_WORLD, boundary_type, DM_BOUNDARY_NONE, DMDA_STENCIL_BOX,
							mx+1, mz+1, Px, Pz, dof, stencil_width, NULL, NULL, &da_Thermal); CHKERRQ(ierr);
	} else {
		//
		// NOTE: periodic_boundary is currently not used in 3D
		//
		ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\nWARNING: periodic_boundary is currently not supported in 3D mode\n\n");

		ierr = DMDACreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_BOX,
							mx+1, my+1, mz+1, Px, Py, Pz, dof, stencil_width, NULL, NULL, NULL, &da_Thermal); CHKERRQ(ierr);
	}

	ierr = DMSetFromOptions(da_Thermal);CHKERRQ(ierr);
	ierr = DMSetUp(da_Thermal);CHKERRQ(ierr);

	ierr = DMDASetFieldName(da_Thermal,0,"T");CHKERRQ(ierr);




	ierr = PetscCalloc1(GaussQuad*T_NE,&NT); CHKERRQ(ierr);
	ierr = PetscCalloc1(GaussQuad*T_NE,&NT_x); CHKERRQ(ierr);
	if (dimensions == 3) {
		ierr = PetscCalloc1(GaussQuad*T_NE,&NT_y); CHKERRQ(ierr);
	}
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

	if (dimensions == 2) {
		montaKeThermal_general_2d(TKe,TMe,TFe);
	} else {
		montaKeThermal_general_3d(TKe,TMe,TFe);
	}

	if (dimensions == 2) {
		ierr = DMDASetUniformCoordinates(da_Thermal, 0.0, Lx, -depth, 0.0, 0.0, 0.0); CHKERRQ(ierr);
	} else {
		ierr = DMDASetUniformCoordinates(da_Thermal, 0.0, Lx, 0.0, Ly, -depth, 0.0); CHKERRQ(ierr);
	}

	/* Generate a matrix with the correct non-zero pattern of type AIJ. This will work in parallel and serial */
	ierr = DMSetMatType(da_Thermal,MATAIJ);CHKERRQ(ierr);
	ierr = DMCreateMatrix(da_Thermal,&TA);CHKERRQ(ierr);
	ierr = DMCreateMatrix(da_Thermal,&TB);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da_Thermal,&Temper);CHKERRQ(ierr);
	if (dimensions == 2) {
		ierr = DMCreateGlobalVector(da_Thermal,&Temper_0);CHKERRQ(ierr);
	}
	ierr = DMCreateGlobalVector(da_Thermal,&Tf);CHKERRQ(ierr);

	ierr = DMCreateGlobalVector(da_Thermal,&Pressure_aux);CHKERRQ(ierr);



	ierr = DMCreateGlobalVector(da_Thermal,&geoq);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da_Thermal,&geoq_rho);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da_Thermal,&geoq_H);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da_Thermal,&geoq_strain);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da_Thermal,&geoq_strain_rate);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da_Thermal,&geoq_cont);CHKERRQ(ierr);

	ierr = DMCreateGlobalVector(da_Thermal,&dRho);CHKERRQ(ierr);

	ierr = Thermal_init_2d(Temper,da_Thermal);

	ierr = DMCreateLocalVector(da_Thermal,&local_FT);
	ierr = DMCreateLocalVector(da_Thermal,&local_Temper);
	ierr = DMCreateLocalVector(da_Thermal,&local_TC);

	ierr = DMCreateLocalVector(da_Thermal,&local_P_aux);

	ierr = DMCreateLocalVector(da_Thermal,&local_geoq);
	ierr = DMCreateLocalVector(da_Thermal,&local_geoq_rho);
	ierr = DMCreateLocalVector(da_Thermal,&local_geoq_H);
	ierr = DMCreateLocalVector(da_Thermal,&local_geoq_strain);
	ierr = DMCreateLocalVector(da_Thermal,&local_geoq_strain_rate);
	ierr = DMCreateLocalVector(da_Thermal,&local_geoq_cont);

	ierr = DMCreateLocalVector(da_Thermal,&local_dRho);


	ierr = DMCreateGlobalVector(da_Thermal,&Temper_Cond);CHKERRQ(ierr);


	PetscScalar **ff2d;
	PetscScalar ***ff3d;
	PetscInt M, N, P;
	PetscInt sx, sy, sz, mmx, mmy, mmz;
	PetscInt i, j, k;

	ierr = VecZeroEntries(local_FT); CHKERRQ(ierr);

	if (dimensions == 2) {
		ierr = DMDAGetInfo(da_Thermal, 0, &M, &P, NULL, 0, 0, 0, 0, 0, 0, 0, 0, 0); CHKERRQ(ierr);
		ierr = DMDAVecGetArray(da_Thermal, local_FT, &ff2d); CHKERRQ(ierr);
		ierr = DMDAGetCorners(da_Thermal, &sx, &sz, NULL, &mmx, &mmz, NULL); CHKERRQ(ierr);

		for (k=sz; k<sz+mmz; k++) {
			for (i=sx; i<sx+mmx; i++) {
				ff2d[k][i] = 1.0;

				if (periodic_boundary==0){
					if (i==0   && bcT_left==1) ff2d[k][i] = 0.0;

					if (i==M-1 && bcT_right==1)ff2d[k][i] = 0.0;
				}

				if (k==0   && bcT_bot==1) ff2d[k][i] = 0.0;

				if (k==P-1 && bcT_top==1) ff2d[k][i] = 0.0;
			}
		}

		ierr = DMDAVecRestoreArray(da_Thermal, local_FT, &ff2d); CHKERRQ(ierr);
	} else {
		ierr = DMDAGetInfo(da_Thermal, 0, &M, &N, &P, 0, 0, 0,  0, 0, 0, 0, 0, 0);CHKERRQ(ierr);
		ierr = DMDAVecGetArray(da_Thermal, local_FT, &ff3d);CHKERRQ(ierr);
		ierr = DMDAGetCorners(da_Thermal, &sx, &sy, &sz, &mmx, &mmy, &mmz);CHKERRQ(ierr);

		for (k=sz; k<sz+mmz; k++) {
			for (j=sy; j<sy+mmy; j++) {
				for (i=sx; i<sx+mmx; i++) {
					ff3d[k][j][i] = 1.0;

					if (i==0   && bcT_left==1) ff3d[k][j][i] = 0.0;

					if (i==M-1 && bcT_right==1)ff3d[k][j][i] = 0.0;

					if (k==0   && bcT_bot==1) ff3d[k][j][i] = 0.0;

					if (k==P-1 && bcT_top==1) ff3d[k][j][i] = 0.0;
				}
			}
		}

		ierr = DMDAVecRestoreArray(da_Thermal, local_FT, &ff3d);CHKERRQ(ierr);
	}

	ierr = DMLocalToGlobalBegin(da_Thermal, local_FT, INSERT_VALUES, Temper_Cond); CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(da_Thermal, local_FT, INSERT_VALUES, Temper_Cond); CHKERRQ(ierr);


	///////

	if (dimensions == 2) {
		VecCopy(Temper,Temper_0);

		VecScale(Temper_Cond, -1.0);
		VecShift(Temper_Cond,  1.0);
		VecPointwiseMult(Temper_0,Temper_0,Temper_Cond);
		VecScale(Temper_Cond, -1.0);
		VecShift(Temper_Cond,  1.0);
	}



	ierr = KSPCreate(PETSC_COMM_WORLD,&T_ksp);CHKERRQ(ierr);
	ierr = KSPSetOptionsPrefix(T_ksp,"thermal_"); /* stokes */ CHKERRQ(ierr);

	PetscTime(&Tempo2p);
	PetscPrintf(PETSC_COMM_WORLD, "temperature (creation): %lf s\n",Tempo2p-Tempo1p);

	PetscFunctionReturn(0);
}


PetscErrorCode build_thermal(int dimensions)
{

	PetscErrorCode ierr;

	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

	PetscLogDouble Tempo1p,Tempo2p;

	PetscTime(&Tempo1p);

	ierr = MatZeroEntries(TA);CHKERRQ(ierr);

	ierr = MatZeroEntries(TB);CHKERRQ(ierr);
	ierr = VecZeroEntries(Tf);CHKERRQ(ierr);

	if (dimensions == 2) {
		ierr = AssembleA_Thermal_2d(TA,da_Thermal,TKe,TMe,TFe,da_Veloc,Veloc_fut);CHKERRQ(ierr);
		ierr = AssembleF_Thermal_2d(Tf,da_Thermal,TKe,TMe,TFe,da_Veloc,Veloc);CHKERRQ(ierr);
	} else {
		ierr = AssembleA_Thermal_3d(TA,da_Thermal,TKe,TMe,TFe,da_Veloc,Veloc_fut);CHKERRQ(ierr);
		ierr = AssembleF_Thermal_3d(Tf,da_Thermal,TKe,TMe,TFe,da_Veloc,Veloc);CHKERRQ(ierr);
	}

	PetscTime(&Tempo2p);
	PetscPrintf(PETSC_COMM_WORLD, "Thermal (building): %lf s\n",Tempo2p-Tempo1p);

	PetscFunctionReturn(0);

}

PetscErrorCode solve_thermal(int dimensions)
{
	PetscErrorCode ierr;
	PetscLogDouble Tempo1,Tempo2;

	int rank;

	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

	PetscTime(&Tempo1);

	/* SOLVE */

	ierr = KSPSetOperators(T_ksp,TA,TA);CHKERRQ(ierr);

	ierr = KSPSetFromOptions(T_ksp);CHKERRQ(ierr);


	ierr = KSPSetTolerances(T_ksp,rtol,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);

	ierr = KSPSolve(T_ksp,Tf,Temper);CHKERRQ(ierr);

	if (dimensions == 2) {
		VecPointwiseMult(Temper,Temper,Temper_Cond); ///zero at the b.c.
		VecAXPY(Temper,1.0,Temper_0); ///applying Teloc_0 in Teloc at the b.c.

		if (basal_heat>0) Heat_flow_at_the_base();
	}

	PetscTime(&Tempo2);
	PetscPrintf(PETSC_COMM_WORLD, "Thermal (solution): %lf s\n",Tempo2-Tempo1);

	PetscFunctionReturn(0);
}


PetscErrorCode Heat_flow_at_the_base(){
	PetscErrorCode ierr;
	PetscScalar  **TT;

	ierr = DMGlobalToLocalBegin(da_Thermal,Temper,INSERT_VALUES,local_Temper);
	ierr = DMGlobalToLocalEnd(  da_Thermal,Temper,INSERT_VALUES,local_Temper);

	ierr = DMDAVecGetArray(da_Thermal,local_Temper,&TT);CHKERRQ(ierr);

	PetscInt       sx,sz,mmx,mmz;

	ierr = DMDAGetCorners(da_Thermal,&sx,&sz,NULL,&mmx,&mmz,NULL);CHKERRQ(ierr);

	int k,i;

	PetscScalar condutivity = kappa*RHOM*c_heat_capacity;
	PetscScalar delta_T_basal = basal_heat*dz_const/condutivity;

	for (k=sz; k<sz+mmz; k++) {
		for (i=sx; i<sx+mmx; i++) {
			if (k==0){
				TT[k][i]=TT[k+1][i]+delta_T_basal;
			}
		}
	}

	ierr = DMDAVecRestoreArray(da_Thermal,local_Temper,&TT);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(da_Thermal,local_Temper,INSERT_VALUES,Temper);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(da_Thermal,local_Temper,INSERT_VALUES,Temper);CHKERRQ(ierr);

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
	PetscPrintf(PETSC_COMM_WORLD, "Thermal (destroying): %lf\n",Tempo2-Tempo1);


	PetscFunctionReturn(0);
}

PetscErrorCode write_all_(int cont,Vec u, char *variable_name, PetscInt binary_out)
{
	int rank;

	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	PetscLogDouble Tempo1,Tempo2;
	PetscTime(&Tempo1);

	PetscViewer viewer;

	char nome[100];

	if (binary_out==0){
		sprintf(nome,"%s_%d.txt",variable_name,cont);
		PetscViewerASCIIOpen(PETSC_COMM_WORLD,nome,&viewer);
	}
	else {
		sprintf(nome,"%s_%d.bin",variable_name,cont);
		PetscViewerBinaryOpen(PETSC_COMM_WORLD,nome,FILE_MODE_WRITE,&viewer);
	}

	VecView(u,viewer);
	PetscViewerDestroy(&viewer);

	PetscTime(&Tempo2);
	if (rank==0 && cont>=0) printf("%s (writing): %lf s\n",variable_name,Tempo2-Tempo1);

	PetscFunctionReturn(0);
}

PetscErrorCode write_pressure(int cont, PetscInt binary_out)
{
	char variable_name[100];

	sprintf(variable_name,"pressure");
	write_all_(cont,Pressure_aux,variable_name,binary_out);

	PetscFunctionReturn(0);
}

PetscErrorCode write_geoq_(int cont, PetscInt binary_out)
{

	char variable_name[100];

	sprintf(variable_name,"viscosity");
	write_all_(cont,geoq,variable_name,binary_out);

	sprintf(variable_name,"density");
	write_all_(cont,geoq_rho,variable_name,binary_out);

	sprintf(variable_name,"heat");
	write_all_(cont,geoq_H,variable_name,binary_out);

	sprintf(variable_name,"strain");
	write_all_(cont,geoq_strain,variable_name,binary_out);

	sprintf(variable_name,"strain_rate");
	write_all_(cont,geoq_strain_rate,variable_name,binary_out);

	PetscFunctionReturn(0);
}
