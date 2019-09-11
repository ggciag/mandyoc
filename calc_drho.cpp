#include <petscksp.h>
#include <petscdmda.h>

extern Vec dRho;
extern Vec Temper;
extern Vec geoq_rho;

extern double alpha_exp_thermo;
extern double gravity;
extern double RHOM;

PetscErrorCode calc_drho()
{
	PetscErrorCode ierr=0;
	
	VecCopy(Temper, dRho); CHKERRQ(ierr);
	
	PetscScalar a;
	
	a=alpha_exp_thermo*gravity*RHOM; CHKERRQ(ierr);//!!! checar
	
	//VecScale(dRho,a); CHKERRQ(ierr);
	VecAXPBY(dRho,-gravity,a,geoq_rho); CHKERRQ(ierr);
	
	return ierr;
	
}

PetscErrorCode write_pressure(int cont);

extern Vec local_P;
extern Vec Pressure;

extern Vec local_geoq_rho;

extern DM da_Veloc;

extern DM da_Thermal;

extern long Nx,Nz;

extern long V_GT;
extern long T_NE;

extern double Lx, depth;


extern int n_interfaces;
extern PetscScalar *interfaces;

extern PetscScalar *inter_rho;

typedef struct {
	PetscScalar u;
	PetscScalar w;
} Stokes;

PetscErrorCode DMDAGetElementCorners(DM da,PetscInt *sx,PetscInt *sz,PetscInt *mx,PetscInt *mz);

PetscErrorCode calc_pressure()
{
	Stokes					**pp;
	
	PetscErrorCode         ierr;
	
	ierr = VecZeroEntries(local_P);CHKERRQ(ierr);
	
	ierr = DMGlobalToLocalBegin(da_Veloc,Pressure,INSERT_VALUES,local_P);
	ierr = DMGlobalToLocalEnd(  da_Veloc,Pressure,INSERT_VALUES,local_P);
	
	ierr = DMDAVecGetArray(da_Veloc,local_P,&pp);CHKERRQ(ierr);
	
	
	
	//PetscScalar             **qq_rho;
	
	//ierr = DMGlobalToLocalBegin(da_Thermal,geoq_rho,INSERT_VALUES,local_geoq_rho);
	//ierr = DMGlobalToLocalEnd(  da_Thermal,geoq_rho,INSERT_VALUES,local_geoq_rho);
	
	//ierr = DMDAVecGetArray(da_Thermal,local_geoq_rho,&qq_rho);CHKERRQ(ierr);
	
	PetscInt i,j,k,c,n,g;
	
	MatStencil indr[T_NE],ind[V_GT];
	
	int M,P;
	double xx,zz;
	
	PetscFunctionBeginUser;
	ierr = DMDAGetInfo(da_Thermal,0,&M,&P,NULL,0,0,0, 0,0,0,0,0,0);CHKERRQ(ierr);
	
	printf("%d %d P\n",M,P);
	
	//
	
	PetscInt               sex1,sez1,mx1,mz1;
	
	PetscInt               sex,sez,mx,mz;
	PetscInt               ei,ek;
	
	ierr = DMDAGetElementCorners(da_Thermal,&sex,&sez,&mx,&mz);CHKERRQ(ierr);
	
	ierr = DMDAGetElementCorners(da_Veloc,&sex1,&sez1,&mx1,&mz1);CHKERRQ(ierr);
	
	if ((sex-sex1!=0)||(sez-sez1!=0)||(mx-mx1!=0)||(mz-mz1!=0)){
		printf("%d %d %d %d\n",sex-sex1,sez-sez1,mx-mx1,mz-mz1);
		SETERRQ1(PETSC_COMM_WORLD,1,"Wrong partition (temper,velocity)\n",1);
	}
	
	printf("%d %d %d %d Pressure\n",sez,sez+mz,sex,sex+mx);
	//printf("%d %d %d %d Pressure\n",sez,sez+mz,sex,sex+mx);
	
	PetscReal interp_interfaces[n_interfaces];
	
	for (ek = sez; ek < sez+mz; ek++) {
		for (ei = sex; ei < sex+mx; ei++) {
			indr[0].i=ei  ; indr[0].j=ek  ;
			indr[1].i=ei+1; indr[1].j=ek  ;
			indr[2].i=ei  ; indr[2].j=ek+1;
			indr[3].i=ei+1; indr[3].j=ek+1;
			
			xx = ei*Lx/(M-1);
			zz = -(P-1-ek)*depth/(P-1);
			
			
			for (int in=0;in<n_interfaces;in++){
				float rfac = (0.5);
				interp_interfaces[in] = interfaces[ei + Nx*in] * rfac;
				
				rfac = (0.5);
				interp_interfaces[in] += interfaces[(ei+1) + Nx*in] * rfac;
				
			}
			
			double pressure_cumulat = 0;
			
			int verif=0;
			/*for (int in=0;in<n_interfaces && verif==0;in++){
				if (zz<interp_interfaces[in]){
					verif=1;
					pressure_cumulat += inter_rho[in]*gravity*
											(interp_interfaces[in]-zz);
				}
				else{
					pressure_cumulat += inter_rho[in]*gravity*
											(interp_interfaces[in+1]-interp_interfaces[in]);
				}
			}
			if (verif==0){
				pressure_cumulat += inter_rho[n_interfaces]*gravity*(-zz);
			}*/
			/*int in;
			for (in=n_interfaces; in>0 && verif==0;in--){
				if (zz>=interp_interfaces[in]){
					verif=1;
					pressure_cumulat += inter_rho[in]*gravity*
					(interp_interfaces[in]-zz);
				}
				else{
					pressure_cumulat += inter_rho[in]*gravity*
					(interp_interfaces[in+1]-interp_interfaces[in]);
				}
			}
			if (verif==0){
				pressure_cumulat += inter_rho[0]*gravity*(-zz);
			}*/
			
			
			if (zz>interp_interfaces[n_interfaces-1]){
				pressure_cumulat += inter_rho[n_interfaces]*gravity*(-zz);
			}
			else{
				pressure_cumulat += inter_rho[n_interfaces]*gravity*
									(-interp_interfaces[n_interfaces-1]);
				for (int in=n_interfaces-1; in>0 && verif==0;in--){
					if (zz>interp_interfaces[in-1]){
						pressure_cumulat += inter_rho[in]*gravity*
											(interp_interfaces[in]-zz);
						verif=1;
					}
					else{
						pressure_cumulat += inter_rho[in]*gravity*
											(interp_interfaces[in]-interp_interfaces[in-1]);
					}
				}
				if (verif==0){
					pressure_cumulat += inter_rho[0]*gravity*(interp_interfaces[0]-zz);
				}
			}
			
			pp[ek][ei].u = pressure_cumulat;
			
			//printf("%d ",in);
			//if (ek==1) printf("pp = %lf\n",pressure_cumulat);
			
		}
	}
	
	
	//ierr = DMDAVecRestoreArray(da_Thermal,local_geoq_rho,&qq_rho);CHKERRQ(ierr);
	
	ierr = DMDAVecRestoreArray(da_Veloc,local_P,&pp);CHKERRQ(ierr);
	
	ierr = DMLocalToGlobalBegin(da_Veloc,local_P,INSERT_VALUES,Pressure);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(da_Veloc,local_P,INSERT_VALUES,Pressure);CHKERRQ(ierr);
	
	write_pressure(-1);
	
	return ierr;
	
}
