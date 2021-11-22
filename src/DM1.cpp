#include <petscksp.h>
#include <petscdmda.h>

#include <petsctime.h>

PetscErrorCode montaKeThermal_simplif(double *Ke_local,double *Ke,double  kappa_eff);
PetscErrorCode ascii2bin(char *s1, char *s2);
PetscErrorCode mean_value_periodic_boundary(DM da,Vec F,Vec local_F, PetscScalar **ff,int esc);

extern double alpha_exp_thermo;
extern double gravity;

extern int DIMEN;

extern PetscReal *TCe;
extern PetscReal *TCe_fut;

extern long Nx,Nz;

extern long T_NE;
extern PetscReal *Ttotal;
extern PetscReal *Ttotal_b;

extern double alpha_thermal;
extern double comp_alpha_thermal;
extern double dt_calor_sec;

extern PetscReal *T_vec_aux_ele;
extern PetscReal *T_vec_aux_ele_final;

extern PetscReal *v_vec_aux_ele;

extern Vec Temper;

extern Vec Temper_Cond;

extern double Delta_T;

extern double Lx, Ly, depth;

extern double h_air;

extern int T_initial_cond;

extern double seg_per_ano;

extern double kappa;


extern double beta_max;
extern double ramp_begin;
extern double ramp_end;

extern double H_lito;

extern Vec local_FT;
extern Vec local_Temper;

extern Vec geoq_H,local_geoq_H;

extern Vec local_TC;

extern Vec local_V;

extern PetscInt temper_extern;
extern PetscInt veloc_extern;
extern PetscInt bcv_extern;

double Thermal_profile(double t, double zz);

PetscReal Temper3(double xx,double zz);
PetscReal Temper4(double xx,double zz);
PetscReal Temper4_3D(double xx,double yy,double zz);


extern PetscInt WITH_NON_LINEAR;
extern PetscInt WITH_ADIABATIC_H;
extern PetscInt WITH_RADIOGENIC_H;

extern unsigned int seed;

extern PetscInt periodic_boundary;

extern PetscInt high_kappa_in_asthenosphere;


typedef struct {
	PetscScalar u;
	PetscScalar w;
	//PetscScalar p; check for the 3d version
} Stokes;


PetscErrorCode DMDAGetLocalElementSize(DM da,PetscInt *mxl,PetscInt *mzl)
{
	PetscInt       m,p,M,P;
	PetscInt       sx,sz;
	PetscErrorCode ierr;
	
	PetscFunctionBeginUser;
	ierr = DMDAGetInfo(da,0,&M,&P,NULL,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = DMDAGetCorners(da,&sx,&sz,NULL,&m,&p,NULL);CHKERRQ(ierr);
	
	if (mxl) {
		*mxl = m;
		if ((sx+m) == M) *mxl = m-1;  /* last proc */
	}
	if (mzl) {
		*mzl = p;
		if ((sz+p) == P) *mzl = p-1;  /* last proc */
	}
	PetscFunctionReturn(0);
}

PetscErrorCode DMDAGetLocalElementSize_3d(DM da,PetscInt *mxl,PetscInt *myl,PetscInt *mzl)
{
	PetscInt       m,n,p,M,N,P;
	PetscInt       sx,sy,sz;
	PetscErrorCode ierr;
	
	PetscFunctionBeginUser;
	ierr = DMDAGetInfo(da,0,&M,&N,&P,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = DMDAGetCorners(da,&sx,&sy,&sz,&m,&n,&p);CHKERRQ(ierr);
	
	if (mxl) {
		*mxl = m;
		if ((sx+m) == M) *mxl = m-1;  /* last proc */
	}
	if (myl) {
		*myl = n;
		if ((sy+n) == N) *myl = n-1;  /* last proc */
	}
	if (mzl) {
		*mzl = p;
		if ((sz+p) == P) *mzl = p-1;  /* last proc */
	}
	PetscFunctionReturn(0);
}


PetscErrorCode DMDAGetElementCorners(DM da,PetscInt *sx,PetscInt *sz,PetscInt *mx,PetscInt *mz)
{
	PetscInt       si,sk;
	PetscErrorCode ierr;
	
	PetscFunctionBeginUser;
	ierr = DMDAGetGhostCorners(da,&si,&sk,NULL,0,0,0);CHKERRQ(ierr);
	
	if (sx) {
		*sx = si;
		if (si != 0) *sx = si+1+periodic_boundary;
	}
	if (sz) {
		*sz = sk;
		if (sk != 0) *sz = sk+1+periodic_boundary;
	}
	ierr = DMDAGetLocalElementSize(da,mx,mz);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

PetscErrorCode DMDAGetElementCorners3d(DM da,PetscInt *sx,PetscInt *sy,PetscInt *sz,PetscInt *mx,PetscInt *my,PetscInt *mz)
{
	PetscInt       si,sj,sk;
	PetscErrorCode ierr;
	
	PetscFunctionBeginUser;
	ierr = DMDAGetGhostCorners(da,&si,&sj,&sk,0,0,0);CHKERRQ(ierr);
	
	if (sx) {
		*sx = si;
		if (si != 0) *sx = si+1;
	}
	if (sy) {
		*sy = sj;
		if (sj != 0) *sy = sj+1;
	}
	if (sz) {
		*sz = sk;
		if (sk != 0) *sz = sk+1;
	}
	ierr = DMDAGetLocalElementSize_3d(da,mx,my,mz);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}


PetscErrorCode AssembleA_Thermal(Mat A,DM thermal_da,PetscReal *TKe,PetscReal *TMe,PetscReal *TFe,
								 DM veloc_da, Vec Veloc_total)
{
	
	PetscErrorCode ierr;
	
	PetscInt               sex,sez,mx,mz;
	PetscInt               ei,ek;
	
	
	Stokes					**VV;
	
	ierr = VecZeroEntries(local_V);CHKERRQ(ierr);
	
	ierr = DMGlobalToLocalBegin(veloc_da,Veloc_total,INSERT_VALUES,local_V);
	ierr = DMGlobalToLocalEnd(  veloc_da,Veloc_total,INSERT_VALUES,local_V);
	
	ierr = DMDAVecGetArray(veloc_da,local_V,&VV);CHKERRQ(ierr);
	
	
	
	////////
	
	PetscScalar					**TTC;
	
	ierr = VecZeroEntries(local_TC);CHKERRQ(ierr);
	
	ierr = DMGlobalToLocalBegin(thermal_da,Temper_Cond,INSERT_VALUES,local_TC);
	ierr = DMGlobalToLocalEnd(  thermal_da,Temper_Cond,INSERT_VALUES,local_TC);
	
	ierr = DMDAVecGetArray(thermal_da,local_TC,&TTC);CHKERRQ(ierr);
	
	////////

	PetscScalar **tt;

	ierr = VecZeroEntries(local_Temper);CHKERRQ(ierr);
	
	ierr = DMGlobalToLocalBegin(thermal_da,Temper,INSERT_VALUES,local_Temper);
	ierr = DMGlobalToLocalEnd(thermal_da,Temper,INSERT_VALUES,local_Temper);
	
	ierr = DMDAVecGetArray(thermal_da,local_Temper,&tt);CHKERRQ(ierr);
	
	
	
	PetscInt               M,P;

	
	MatStencil ind1[1],ind2[1],ind[4];
	
	PetscScalar u[4*4],val_cond[1];
	
	PetscFunctionBeginUser;
	ierr = DMDAGetInfo(thermal_da,0,&M,&P,NULL,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	
	
	
	
	
	PetscInt i,j,k,c;
	
	PetscInt       sx,sz,mmx,mmz;
	
	ierr = DMDAGetCorners(thermal_da,&sx,&sz,NULL,&mmx,&mmz,NULL);CHKERRQ(ierr);
	
	for (k=sz; k<sz+mmz; k++) {
		for (i=sx; i<sx+mmx; i++) {
			if (periodic_boundary==1 && i==Nx-1){
				ind1[0].i = i;
				ind1[0].j = k;
				val_cond[0]=1.0;
				ierr = MatSetValuesStencil(A,1,ind1,1,ind1,val_cond,ADD_VALUES);

				val_cond[0]=-1.0;
				ind2[0].i = Nx;
				ind2[0].j = k;
				ierr = MatSetValuesStencil(A,1,ind1,1,ind2,val_cond,ADD_VALUES);
			}
			else {
				ind1[0].i = i;
				ind1[0].j = k;
				
				val_cond[0] = 1.0-TTC[k][i];
				ierr = MatSetValuesStencil(A,1,ind1,1,ind1,val_cond,ADD_VALUES);
			}
		}
	}
	
	ierr = DMDAGetElementCorners(thermal_da,&sex,&sez,&mx,&mz);CHKERRQ(ierr);

	PetscReal tt_ele[4],tt_mean,tt_ref=1330.0;
	double kappa_eff;
	
	for (ek = sez; ek < sez+mz; ek++) {
		for (ei = sex; ei < sex+mx; ei++) {
			
			ind[0].i=ei  ; ind[0].j=ek  ;
			ind[1].i=ei+1; ind[1].j=ek  ;
			ind[2].i=ei  ; ind[2].j=ek+1;
			ind[3].i=ei+1; ind[3].j=ek+1;
			
			for (i=0;i<4;i++){
				v_vec_aux_ele[i*2+0] = VV[ind[i].j][ind[i].i].u;
				v_vec_aux_ele[i*2+1] = VV[ind[i].j][ind[i].i].w;
			}

			kappa_eff = kappa;
			if (high_kappa_in_asthenosphere==1){
				tt_ele[0] = tt[ek][ei];
				tt_ele[1] = tt[ek][ei+1];
				tt_ele[2] = tt[ek+1][ei];
				tt_ele[3] = tt[ek+1][ei+1];

				tt_mean = (tt_ele[0]+tt_ele[1]+tt_ele[2]+tt_ele[3])/4;

				if (tt_mean > tt_ref){
					kappa_eff=kappa*(1.0 + (tt_mean-tt_ref)*5);
				}
				if (kappa_eff>kappa*50) kappa_eff = kappa*50;
			}
			
			montaKeThermal_simplif(TCe_fut, TKe, kappa_eff);//check: modify the simplified version
			
			for (c=0;c<T_NE*T_NE;c++) Ttotal[c] = TMe[c] + alpha_thermal*dt_calor_sec*TCe_fut[c];
			
			if (periodic_boundary==1){
				if (ei+1==Nx-1){
					ind[1].i=Nx;
					ind[3].i=Nx;
				}
			}
			

			for (i=0;i<4;i++){
				if (TTC[ind[i].j][ind[i].i]==0){
					for (j=0;j<4;j++) u[i*4+j]=0.0;
				}
				else {
					for (j=0;j<4;j++) u[i*4+j]=Ttotal[i*4+j];
				}
			}
			
			ierr = MatSetValuesStencil(A,4,ind,4,ind,u,ADD_VALUES);
		}
	
	}
	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
	
	ierr = DMDAVecRestoreArray(veloc_da,local_V,&VV);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(thermal_da,local_TC,&TTC);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(thermal_da,local_Temper,&tt);CHKERRQ(ierr);
	
	
	PetscFunctionReturn(0);
}

PetscErrorCode AssembleF_Thermal(Vec F,DM thermal_da,PetscReal *TKe,PetscReal *TMe,PetscReal *TFe,
								 DM veloc_da, Vec Veloc_total)
{

	PetscScalar              **ff,**tt,**HH;
	PetscInt               M,P;
	PetscErrorCode         ierr;
	
	PetscInt               sex,sez,mx,mz;
	PetscInt               ei,ek;
	
	
	PetscFunctionBeginUser;
	ierr = DMDAGetInfo(thermal_da,0,&M,&P,NULL,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	

	
	
	Stokes					**VV;

	ierr = VecZeroEntries(local_V);CHKERRQ(ierr);
	
	ierr = DMGlobalToLocalBegin(veloc_da,Veloc_total,INSERT_VALUES,local_V);
	ierr = DMGlobalToLocalEnd(  veloc_da,Veloc_total,INSERT_VALUES,local_V);
	
	ierr = DMDAVecGetArray(veloc_da,local_V,&VV);CHKERRQ(ierr);
	
	
	////////
	
	PetscScalar					**TTC;
	
	
	ierr = VecZeroEntries(local_TC);CHKERRQ(ierr);
	
	ierr = DMGlobalToLocalBegin(thermal_da,Temper_Cond,INSERT_VALUES,local_TC);
	ierr = DMGlobalToLocalEnd(  thermal_da,Temper_Cond,INSERT_VALUES,local_TC);
	
	ierr = DMDAVecGetArray(thermal_da,local_TC,&TTC);CHKERRQ(ierr);
	
	////////
	
	
	
	/* get acces to the vector */

	ierr = VecZeroEntries(local_FT);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(thermal_da,local_FT,&ff);CHKERRQ(ierr);
	
	

	
	ierr = VecZeroEntries(local_Temper);CHKERRQ(ierr);
	
	ierr = DMGlobalToLocalBegin(thermal_da,Temper,INSERT_VALUES,local_Temper);
	ierr = DMGlobalToLocalEnd(thermal_da,Temper,INSERT_VALUES,local_Temper);
	
	ierr = DMDAVecGetArray(thermal_da,local_Temper,&tt);CHKERRQ(ierr);
	
	
	
	ierr = VecZeroEntries(local_geoq_H);CHKERRQ(ierr);
	
	ierr = DMGlobalToLocalBegin(thermal_da,geoq_H,INSERT_VALUES,local_geoq_H);
	ierr = DMGlobalToLocalEnd(thermal_da,geoq_H,INSERT_VALUES,local_geoq_H);
	
	ierr = DMDAVecGetArray(thermal_da,local_geoq_H,&HH);CHKERRQ(ierr);
	
	
	PetscInt i,j,k,c;
	
	MatStencil ind[4];
	
	PetscReal H_efetivo;
	
	
	PetscReal tt_ele[4],tt_mean,tt_ref=1330.0;
	double kappa_eff;

	ierr = DMDAGetElementCorners(thermal_da,&sex,&sez,&mx,&mz);CHKERRQ(ierr);
	for (ek = sez; ek < sez+mz; ek++) {
		for (ei = sex; ei < sex+mx; ei++) {
			
			
			ind[0].i=ei  ; ind[0].j=ek  ;
			ind[1].i=ei+1; ind[1].j=ek  ;
			ind[2].i=ei  ; ind[2].j=ek+1;
			ind[3].i=ei+1; ind[3].j=ek+1;

			
			
			for (i=0;i<4;i++){
				v_vec_aux_ele[i*2+0] = VV[ind[i].j][ind[i].i].u;
				v_vec_aux_ele[i*2+1] = VV[ind[i].j][ind[i].i].w;
			}

			kappa_eff = kappa;
			if (high_kappa_in_asthenosphere==1){
				tt_ele[0] = tt[ek][ei];
				tt_ele[1] = tt[ek][ei+1];
				tt_ele[2] = tt[ek+1][ei];
				tt_ele[3] = tt[ek+1][ei+1];

				tt_mean = (tt_ele[0]+tt_ele[1]+tt_ele[2]+tt_ele[3])/4;

				if (tt_mean > tt_ref){
					kappa_eff=kappa*(1.0 + (tt_mean-tt_ref)*5);
				}
				if (kappa_eff>kappa*50) kappa_eff = kappa*50;
			}

			montaKeThermal_simplif(TCe, TKe, kappa_eff);//check: modify the simplified version
			
			for (c=0;c<T_NE*T_NE;c++) Ttotal_b[c] = TMe[c] - comp_alpha_thermal*dt_calor_sec*TCe[c];
			
			
			
			
			H_efetivo = 0.0;
			
			if (WITH_RADIOGENIC_H==1){
				double H_mean=0;
				for (j=0;j<T_NE;j++) H_mean+=HH[ind[j].j][ind[j].i];
				H_mean/=T_NE;// Variable radiogenic heat flow here
			
				H_efetivo = H_mean;
			}
			
			if (WITH_ADIABATIC_H==1){
				double T_mean=0;
				for (j=0;j<T_NE;j++) T_mean+=tt[ind[j].j][ind[j].i];
				T_mean/=T_NE;
				T_mean+=273.0; //Converting to Kelvin
				
				double Vz_mean=0;
				for (j=0;j<T_NE;j++) Vz_mean+=VV[ind[j].j][ind[j].i].w;
				Vz_mean/=T_NE;
				
				H_efetivo += -T_mean*alpha_exp_thermo*gravity*Vz_mean;
			}
			
			
			
			
			
			
			for (c=0;c<T_NE;c++){
				T_vec_aux_ele_final[c]=dt_calor_sec*TFe[c]*H_efetivo;
				for (j=0;j<T_NE;j++){
					T_vec_aux_ele_final[c]+=tt[ind[j].j][ind[j].i]*Ttotal_b[c*T_NE+j];
				}
			}
			
			for (i=0;i<4;i++){
				if (TTC[ind[i].j][ind[i].i]==0){
					T_vec_aux_ele_final[i]=0.0;
				}
			}
			
			for (c=0;c<T_NE;c++){
				ff[ind[c].j][ind[c].i] += T_vec_aux_ele_final[c];
			}
			
			
		}
	}
	
	
	
	PetscInt       sx,sz,mmx,mmz;
	
	ierr = DMDAGetCorners(thermal_da,&sx,&sz,NULL,&mmx,&mmz,NULL);CHKERRQ(ierr);
	
	for (k=sz; k<sz+mmz; k++) {
		
		for (i=sx; i<sx+mmx; i++) {
			if (periodic_boundary==1 && i==Nx-1) ;
			else {
			if (TTC[k][i]==0){
				ff[k][i]=tt[k][i];
			}
			}
		}
	
	}
	
	ierr = DMDAVecRestoreArray(thermal_da,local_FT,&ff);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(thermal_da,local_FT,ADD_VALUES,F);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(thermal_da,local_FT,ADD_VALUES,F);CHKERRQ(ierr);
	
	ierr = DMDAVecRestoreArray(veloc_da,local_V,&VV);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(thermal_da,local_TC,&TTC);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(thermal_da,local_Temper,&tt);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(thermal_da,local_geoq_H,&HH);CHKERRQ(ierr);

	if (periodic_boundary==1){

		ierr = DMGlobalToLocalBegin(thermal_da,F,INSERT_VALUES,local_FT);
		ierr = DMGlobalToLocalEnd(  thermal_da,F,INSERT_VALUES,local_FT);
		ierr = DMDAVecGetArray(thermal_da,local_FT,&ff);CHKERRQ(ierr);

		for (k=sz; k<sz+mmz; k++) {
			for (i=sx; i<sx+mmx; i++) {
				if (i==0){
					ff[k][i] += ff[k][i-1];
				}
			}
		}
		ierr = DMDAVecRestoreArray(thermal_da,local_FT,&ff);CHKERRQ(ierr);
		ierr = DMLocalToGlobalBegin(thermal_da,local_FT,INSERT_VALUES,F);CHKERRQ(ierr);
		ierr = DMLocalToGlobalEnd(thermal_da,local_FT,INSERT_VALUES,F);CHKERRQ(ierr);

		ierr = DMGlobalToLocalBegin(thermal_da,F,INSERT_VALUES,local_FT);
		ierr = DMGlobalToLocalEnd(  thermal_da,F,INSERT_VALUES,local_FT);
		ierr = DMDAVecGetArray(thermal_da,local_FT,&ff);CHKERRQ(ierr);

		for (k=sz; k<sz+mmz; k++) {
			for (i=sx; i<sx+mmx; i++) {
				if (i==Nx-1){
					ff[k][i] = 0;
				}
			}
		}

		ierr = DMDAVecRestoreArray(thermal_da,local_FT,&ff);CHKERRQ(ierr);
		ierr = DMLocalToGlobalBegin(thermal_da,local_FT,INSERT_VALUES,F);CHKERRQ(ierr);
		ierr = DMLocalToGlobalEnd(thermal_da,local_FT,INSERT_VALUES,F);CHKERRQ(ierr);

	}
	

	PetscFunctionReturn(0);
}

#include <math.h>


PetscErrorCode Thermal_init(Vec F,DM thermal_da)
{
	if (DIMEN==2){
		Vec                    local_F;
		PetscScalar              **ff;
		PetscInt               M,P;
		PetscErrorCode         ierr;
		
		PetscInt i,k,t;
		
		
		PetscFunctionBeginUser;
		ierr = DMDAGetInfo(thermal_da,0,&M,&P,NULL,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
		
		PetscMPIInt rank;
		ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
		
		
		PetscReal temper_aux=0.0,t1_aux;

		PetscInt ix[1];
		PetscScalar y[1];
		
		PetscInt low,high;
		
		if (temper_extern==1){
			
			char s1[100],s2[100];
			
			sprintf(s1,"input_temperature_0.txt");
			sprintf(s2,"Temper_init.bin");
			
			if (rank==0){
				ierr = ascii2bin(s1,s2); CHKERRQ(ierr);
			}
			MPI_Barrier(PETSC_COMM_WORLD);
			
			
			PetscInt size0;
			PetscViewer    viewer;
			
			VecGetSize(F,&size0);
			
			
			Vec Fprov;
			
			PetscViewerBinaryOpen(PETSC_COMM_WORLD,s2,FILE_MODE_READ,&viewer);
			VecCreate(PETSC_COMM_WORLD,&Fprov);
			VecLoad(Fprov,viewer);
			PetscViewerDestroy(&viewer);
			
			
			VecGetOwnershipRange(Fprov,&low,&high);
			
			
			Vec FN;
			
			DMDACreateNaturalVector(thermal_da,&FN);
			
			
			for (t=low;t<high;t++){
				ix[0] = t;
				VecGetValues(Fprov,1,ix,y);
				VecSetValue(FN,t,y[0], INSERT_VALUES);
			}
			
			VecAssemblyBegin(FN);
			VecAssemblyEnd(FN);
			
			DMDANaturalToGlobalBegin(thermal_da,FN,INSERT_VALUES,F);
			DMDANaturalToGlobalEnd(thermal_da,FN,INSERT_VALUES,F);
			
			
			PetscBarrier(NULL);
			
			VecAssemblyBegin(F);
			VecAssemblyEnd(F);
			
			
			
			
		}
		else{
			/* get access to the vector */
			ierr = DMGetLocalVector(thermal_da,&local_F);CHKERRQ(ierr);
			ierr = VecZeroEntries(local_F);CHKERRQ(ierr);
			ierr = DMDAVecGetArray(thermal_da,local_F,&ff);CHKERRQ(ierr);
			
			PetscInt       sx,sz,mmx,mmz;
			
			ierr = DMDAGetCorners(thermal_da,&sx,&sz,NULL,&mmx,&mmz,NULL);CHKERRQ(ierr);
			
			PetscReal xx,zz,t_inic;
			
			for (k=sz; k<sz+mmz; k++) {
				for (i=sx; i<sx+mmx; i++) {
					
					xx = i*Lx/(M-1);
					zz = -(P-1-k)*depth/(P-1);
					
					zz += h_air;
					
					if (T_initial_cond==0){
						temper_aux=(Delta_T*(P-1-k))/(P-1) + 100*cos(i*3.14159/(M-1));
					}
					
					if (T_initial_cond==1){
						temper_aux=(Delta_T*(P-1-k))/(P-1) + 100*cos(i*3.14159/(M-1));
					}
					if (T_initial_cond==2){
						if (xx<1500.E3) t_inic = 1001.0E6*seg_per_ano;
						else {
							if (xx<2000.E3){
								t_inic = 1000.0E6*seg_per_ano;
							}
							else t_inic = 1001.0E6*seg_per_ano;
						}
						
						temper_aux = Thermal_profile(t_inic, zz);
					}
					if (T_initial_cond==3){
						temper_aux = Temper3(xx,zz);
						
						temper_aux-=(float)rand_r(&seed)/RAND_MAX;
					}
					
					if (T_initial_cond==8){
						temper_aux = Temper4(xx,zz);
					}
					
					
					if (T_initial_cond==9){
						temper_aux=(depth/H_lito)*(Delta_T*(P-1-k))/(P-1);
						
						if (temper_aux>Delta_T) temper_aux=Delta_T;
						
						temper_aux-=(float)rand_r(&seed)/RAND_MAX;
					}
					
					
					if (T_initial_cond==74){
						t1_aux = 0.5*tan((P-1-k)*3./(P-1)-1.5)/tan(1.5)+0.5;
						temper_aux=Delta_T*t1_aux + Delta_T*0.2*tan(i*3./(M-1)-1.5);
					}
					
					if (T_initial_cond==141){
						double z_aux = -depth*(1.0-k*1.0/(Nz-1));
						double x_aux = Lx*(i*1.0/(Nx-1));
						
						if (z_aux>-depth*0.1) temper_aux=0;
						else {
							if (z_aux>-depth*0.2) temper_aux=Delta_T*(1.0-(z_aux+depth*0.2)/(depth*0.1));
							else temper_aux=Delta_T;
						}
						
						if (z_aux>-depth*(1-0.5) && z_aux<=-depth*(1-0.6) && x_aux<=Lx*0.15){
							temper_aux=Delta_T*0.2;
						}
						
						
						
						if (temper_aux>Delta_T) temper_aux=Delta_T;
						
						temper_aux-=(float)rand_r(&seed)/RAND_MAX;
					}
					
					
					if (k==P-1) temper_aux=0.0;
					if (k==0) temper_aux=Delta_T;
					
					if (temper_aux>Delta_T) temper_aux=Delta_T;
					if (temper_aux<0.0) temper_aux=0.0;
					
					ff[k][i]=temper_aux;
				}
			
			}
			

			ierr = DMDAVecRestoreArray(thermal_da,local_F,&ff);CHKERRQ(ierr);
			ierr = DMLocalToGlobalBegin(thermal_da,local_F,INSERT_VALUES,F);CHKERRQ(ierr);
			ierr = DMLocalToGlobalEnd(thermal_da,local_F,INSERT_VALUES,F);CHKERRQ(ierr);
			ierr = DMRestoreLocalVector(thermal_da,&local_F);CHKERRQ(ierr);

		}
			
		if (periodic_boundary==1){
			ierr = mean_value_periodic_boundary(thermal_da,F,local_F,ff,2);
		}
	}
	else {

		Vec                    local_F;
		PetscScalar              ***ff;
		PetscInt               M,N,P;
		PetscErrorCode         ierr;
		
		PetscInt i,j,k,t;
		
		
		PetscFunctionBeginUser;
		ierr = DMDAGetInfo(thermal_da,0,&M,&N,&P,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
		
		PetscMPIInt rank;
		ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
		
		
		PetscReal temper_aux,t1_aux;

		PetscInt ix[1];
		PetscScalar y[1];
		
		PetscInt low,high;
		
		if (temper_extern==1){
			
			char s1[100],s2[100];
			
			sprintf(s1,"input_temperature_0.txt");
			sprintf(s2,"Temper_init.bin");
			
			if (rank==0){
				ierr = ascii2bin(s1,s2); CHKERRQ(ierr);
			}
			MPI_Barrier(PETSC_COMM_WORLD);
			
			
			PetscInt size0;
			PetscViewer    viewer;
			
			VecGetSize(F,&size0);
			
			
			Vec Fprov;
			
			//PetscPrintf(PETSC_COMM_WORLD,"size = %d\n",size);
			
			//PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Temper_init.bin",FILE_MODE_READ,&viewer);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD,s2,FILE_MODE_READ,&viewer);
			VecCreate(PETSC_COMM_WORLD,&Fprov);
			VecLoad(Fprov,viewer);
			PetscViewerDestroy(&viewer);
			
			
			VecGetOwnershipRange(Fprov,&low,&high);
			
			
			Vec FN;
			
			DMDACreateNaturalVector(thermal_da,&FN);
			
			
			for (t=low;t<high;t++){
				ix[0] = t;
				VecGetValues(Fprov,1,ix,y);
				VecSetValue(FN,t,y[0], INSERT_VALUES);
			}
			
			VecAssemblyBegin(FN);
			VecAssemblyEnd(FN);
			
			DMDANaturalToGlobalBegin(thermal_da,FN,INSERT_VALUES,F);
			DMDANaturalToGlobalEnd(thermal_da,FN,INSERT_VALUES,F);
			
			
			//VecView(F,PETSC_VIEWER_STDOUT_WORLD);
			
			PetscBarrier(NULL);
			
			VecAssemblyBegin(F);
			VecAssemblyEnd(F);
			
			
			
			
		}
		else{
			/* get acces to the vector */
			ierr = DMGetLocalVector(thermal_da,&local_F);CHKERRQ(ierr);
			ierr = VecZeroEntries(local_F);CHKERRQ(ierr);
			ierr = DMDAVecGetArray(thermal_da,local_F,&ff);CHKERRQ(ierr);
			
			PetscInt       sx,sy,sz,mmx,mmy,mmz;
			
			ierr = DMDAGetCorners(thermal_da,&sx,&sy,&sz,&mmx,&mmy,&mmz);CHKERRQ(ierr);
			
			PetscReal xx,yy,zz,t_inic;
			
			for (k=sz; k<sz+mmz; k++) {
				for (j=sy; j<sy+mmy; j++) {
					for (i=sx; i<sx+mmx; i++) {
						
						xx = i*Lx/(M-1);
						yy = j*Ly/(N-1);
						zz = -(P-1-k)*depth/(P-1);
						
						zz += h_air;
						
						if (T_initial_cond==0){
							temper_aux=(Delta_T*(P-1-k))/(P-1) + 100*cos(i*3.14159/(M-1));
						}
						
						if (T_initial_cond==1){
							temper_aux=(Delta_T*(P-1-k))/(P-1) + 100*cos(j*3.14159/(N-1))*cos(i*3.14159/(M-1));
						}
						if (T_initial_cond==2){
							if (xx<1500.E3) t_inic = 1001.0E6*seg_per_ano;
							else {
								if (xx<2000.E3){
									t_inic = 1000.0E6*seg_per_ano;
								}
								else t_inic = 1001.0E6*seg_per_ano;
							}
							
							temper_aux = Thermal_profile(t_inic, zz);
						}
						if (T_initial_cond==3){
							temper_aux = Temper3(xx,zz);
							
							temper_aux-=(float)rand_r(&seed)/RAND_MAX;
						}
						
						if (T_initial_cond==8){
							temper_aux = Temper4_3D(xx,yy,zz);
						}
						
						
						if (T_initial_cond==9){
							temper_aux=(depth/H_lito)*(Delta_T*(P-1-k))/(P-1);
							
							if (temper_aux>Delta_T) temper_aux=Delta_T;
							
							temper_aux-=(float)rand_r(&seed)/RAND_MAX;
						}
						
						
						if (T_initial_cond==74){
							t1_aux = 0.5*tan((P-1-k)*3./(P-1)-1.5)/tan(1.5)+0.5;
							temper_aux=Delta_T*t1_aux + Delta_T*0.2*tan(i*3./(M-1)-1.5);
						}
						
						if (T_initial_cond==141){
							double z_aux = -depth*(1.0-k*1.0/(Nz-1));
							double x_aux = Lx*(i*1.0/(Nx-1));
							//if (array[p*3+2]>-depth*(1-0.52) && array[p*3+2]<=-depth*(1-0.55)
							/*if ((i==0 || i==Nx-1) && (z_aux>=-depth*(1.0-0.5+0.1)) && (z_aux<-depth*(1-0.6-0.1))){
								VV[k][j][i].u=300.0/seg_per_ano;
							}*/
							
							if (z_aux>-depth*0.1) temper_aux=0;
							else {
								if (z_aux>-depth*0.2) temper_aux=Delta_T*(1.0-(z_aux+depth*0.2)/(depth*0.1));
								else temper_aux=Delta_T;
							}
							
							if (z_aux>-depth*(1-0.5) && z_aux<=-depth*(1-0.6) && x_aux<=Lx*0.15){
								temper_aux=Delta_T*0.2;
							}
							
							
							
							if (temper_aux>Delta_T) temper_aux=Delta_T;
							
							temper_aux-=(float)rand_r(&seed)/RAND_MAX;
						}
						
						
						if (k==P-1) temper_aux=0.0;
						if (k==0) temper_aux=Delta_T;
						
						if (temper_aux>Delta_T) temper_aux=Delta_T;
						if (temper_aux<0.0) temper_aux=0.0;
						
						ff[k][j][i]=temper_aux;
					}
				}
			}
			
			
			ierr = DMDAVecRestoreArray(thermal_da,local_F,&ff);CHKERRQ(ierr);
			ierr = DMLocalToGlobalBegin(thermal_da,local_F,INSERT_VALUES,F);CHKERRQ(ierr);
			ierr = DMLocalToGlobalEnd(thermal_da,local_F,INSERT_VALUES,F);CHKERRQ(ierr);
			ierr = DMRestoreLocalVector(thermal_da,&local_F);CHKERRQ(ierr);
			
			
		}
	}
		
	
	PetscFunctionReturn(0);
	
}

PetscReal Temper3(double xx,double zz){
	PetscReal Temper;
	
	double t_inic = 1001.0E6*seg_per_ano;
	Temper = Thermal_profile(t_inic, zz);
	
	if (xx<ramp_begin) Temper*=1.0;
	else {
		if (xx<ramp_end){
			Temper*=(1.0+(beta_max-1.0)*(xx-ramp_begin)/(ramp_end-ramp_begin));
		}
		else Temper*=beta_max;
	}
	
	if (Temper>Delta_T) Temper = Delta_T;
	if (Temper<0.0) Temper = 0;
	
	
	return(Temper);
}

PetscReal Temper4(double xx,double zz){
	PetscReal Temper;
	
	double t_inic = 1001.0E6*seg_per_ano;
	Temper = Thermal_profile(t_inic, zz);
	
	Temper*=(1.0+(beta_max-1.0)*(cos(10*xx/Lx)+1.0)/4.0);
	
	if (Temper>Delta_T) Temper = Delta_T;
	if (Temper<0.0) Temper = 0;
	
	
	return(Temper);
}

PetscReal Temper4_3D(double xx,double yy,double zz){
	PetscReal Temper;
	
	double t_inic = 1001.0E6*seg_per_ano;
	Temper = Thermal_profile(t_inic, zz);
	
	Temper*=(1.0+(beta_max-1.0)*(cos(10*xx/Lx)+1.0)*(cos(10*yy/Ly)+1.0)/4.0);
	
	if (Temper>Delta_T) Temper = Delta_T;
	if (Temper<0.0) Temper = 0;
	
	
	return(Temper);
}



double Thermal_profile(double t, double zz){
	
	double T_q = Delta_T;
	double a_q = H_lito;
	
	double PI = 3.14159;
	
	double z = -zz;
	
	double T = z/a_q;
	
	if (z>0){
		for (long n=1;n<100;n++){
			T+=(2.0/PI)*sin(n*PI*z/a_q)*exp(-n*n*PI*PI*kappa*t/(a_q*a_q))/n;
		}
	}
	else T=0.0;
	
	
	T*=T_q;
	
	if (T>T_q) T=T_q;
	
	
	return(T);
}

PetscErrorCode ascii2bin(char *ss1, char *ss2){
	FILE *entra;
	
	entra = fopen(ss1,"r");
	char c,s[100],s1[100];
	fscanf(entra,"%c",&c);
	while(c!='\n') fscanf(entra,"%c",&c);
	fscanf(entra,"%c",&c);
	while(c!='\n') fscanf(entra,"%c",&c);
	fscanf(entra,"%c",&c);
	while(c!='\n') fscanf(entra,"%c",&c);
	fscanf(entra,"%c",&c);
	while(c!='\n') fscanf(entra,"%c",&c);
	
	PetscInt m=0;
	while (!feof(entra)){
		fscanf(entra,"%s",s);
		if (s[0]=='P') {
			fscanf(entra,"%s",s1);
		}
		else {
			m++;
		}
	}
	fclose(entra);
	
	m--;
	
	
	Vec u;
	PetscScalar    v;
	PetscInt	n;
	
	n=m;
	
	VecCreateSeq(PETSC_COMM_SELF,n,&u);
	VecSetFromOptions(u);
	
	entra = fopen(ss1,"r");
	
	fscanf(entra,"%c",&c);
	while(c!='\n') fscanf(entra,"%c",&c);
	fscanf(entra,"%c",&c);
	while(c!='\n') fscanf(entra,"%c",&c);
	fscanf(entra,"%c",&c);
	while(c!='\n') fscanf(entra,"%c",&c);
	fscanf(entra,"%c",&c);
	while(c!='\n') fscanf(entra,"%c",&c);
	
	m=0;
	while (m<n){
		fscanf(entra,"%s",s);
		if (s[0]=='P') {
			fscanf(entra,"%s",s1);
		}
		else {
			v = atof(s);
			VecSetValues(u,1,&m,&v,INSERT_VALUES);
			m++;
		}
	}
	fclose(entra);
	
	VecAssemblyBegin(u);
	VecAssemblyEnd(u);
	
	PetscViewer    viewer;
	
	//PetscPrintf(PETSC_COMM_SELF,"writing vector in binary to vector.dat ...\n");
	PetscViewerBinaryOpen(PETSC_COMM_SELF,ss2,FILE_MODE_WRITE,&viewer);
	VecView(u,viewer);
	PetscViewerDestroy(&viewer);
	VecDestroy(&u);
	
	PetscFunctionReturn(0);

}

PetscErrorCode mean_value_periodic_boundary(DM da,Vec F,Vec local_F, PetscScalar **ff,int esc){

	PetscErrorCode ierr;

	ierr = DMGetLocalVector(da,&local_F);CHKERRQ(ierr);

	ierr = DMGlobalToLocalBegin(da,F,INSERT_VALUES,local_F);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(  da,F,INSERT_VALUES,local_F);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da,local_F,&ff);CHKERRQ(ierr);

	PetscInt       sx,sz,mmx,mmz;

	ierr = DMDAGetCorners(da,&sx,&sz,NULL,&mmx,&mmz,NULL);CHKERRQ(ierr);

	int i,k;

	for (k=sz; k<sz+mmz; k++) {
		for (i=sx; i<sx+mmx; i++) {
			if (i==Nx-1){
				ff[k][Nx-1]= (ff[k][Nx-1] + ff[k][Nx])/esc;
			}
			if (i==0){
				ff[k][0]= (ff[k][-1] + ff[k][0])/esc;
			}
		}
	}

	ierr = DMDAVecRestoreArray(da,local_F,&ff);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(da,local_F,INSERT_VALUES,F);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(da,local_F,INSERT_VALUES,F);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(da,&local_F);CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD,"periodic boundary: done\n");



	PetscFunctionReturn(0);

}