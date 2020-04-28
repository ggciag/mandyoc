#include <petscdm.h>
#include <petscdmda.h>
#include <petscdmswarm.h>

PetscErrorCode mean_value_periodic_boundary(DM da,Vec F,Vec local_F, PetscScalar **ff,int esc);

extern DM dms;

extern DM da_Thermal;

extern Vec geoq,local_geoq;
extern Vec geoq_rho,local_geoq_rho;
extern Vec geoq_cont,local_geoq_cont;
extern Vec geoq_H,local_geoq_H;
extern Vec geoq_strain, local_geoq_strain;

extern Vec Temper,local_Temper;

extern long Nx,Nz;

extern double dx_const;
extern double dz_const;

extern double Lx, depth;

extern PetscInt visc_harmonic_mean;
extern PetscInt visc_const_per_element;

extern PetscScalar *inter_rho;
extern PetscScalar *inter_H;

extern PetscInt periodic_boundary;

PetscErrorCode Swarm2Mesh(){

	PetscErrorCode ierr;
	PetscScalar             **qq,**qq_cont,**qq_rho,**TT,**qq_H,**qq_strain;
	
	ierr = VecSet(geoq,0.0);CHKERRQ(ierr);
	ierr = VecSet(geoq_rho,0.0);CHKERRQ(ierr);
	ierr = VecSet(geoq_cont,0.0);CHKERRQ(ierr);
	ierr = VecSet(geoq_H,0.0);CHKERRQ(ierr);
	ierr = VecSet(geoq_strain,0.0);CHKERRQ(ierr);
	
	ierr = VecZeroEntries(local_geoq);CHKERRQ(ierr);
	ierr = VecZeroEntries(local_geoq_rho);CHKERRQ(ierr);
	ierr = VecZeroEntries(local_geoq_H);CHKERRQ(ierr);
	ierr = VecZeroEntries(local_geoq_strain);CHKERRQ(ierr);
	
	ierr = DMGlobalToLocalBegin(da_Thermal,geoq,INSERT_VALUES,local_geoq);
	ierr = DMGlobalToLocalEnd(  da_Thermal,geoq,INSERT_VALUES,local_geoq);
	
	ierr = DMGlobalToLocalBegin(da_Thermal,geoq_rho,INSERT_VALUES,local_geoq_rho);
	ierr = DMGlobalToLocalEnd(  da_Thermal,geoq_rho,INSERT_VALUES,local_geoq_rho);
	
	ierr = DMGlobalToLocalBegin(da_Thermal,geoq_H,INSERT_VALUES,local_geoq_H);
	ierr = DMGlobalToLocalEnd(  da_Thermal,geoq_H,INSERT_VALUES,local_geoq_H);
	
	ierr = DMGlobalToLocalBegin(da_Thermal,geoq_strain,INSERT_VALUES,local_geoq_strain);
	ierr = DMGlobalToLocalEnd(  da_Thermal,geoq_strain,INSERT_VALUES,local_geoq_strain);
	
	
	ierr = DMDAVecGetArray(da_Thermal,local_geoq,&qq);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da_Thermal,local_geoq_rho,&qq_rho);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da_Thermal,local_geoq_H,&qq_H);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da_Thermal,local_geoq_strain,&qq_strain);CHKERRQ(ierr);
	
	
	ierr = DMGlobalToLocalBegin(da_Thermal,geoq_cont,INSERT_VALUES,local_geoq_cont);
	ierr = DMGlobalToLocalEnd(  da_Thermal,geoq_cont,INSERT_VALUES,local_geoq_cont);
	
	ierr = DMDAVecGetArray(da_Thermal,local_geoq_cont,&qq_cont);CHKERRQ(ierr);
	
	

	PetscInt nlocal,bs,p;
	
	PetscReal *array;
	PetscReal *geoq_fac;
	PetscReal *strain_fac;
	PetscInt *layer_array;
	
	ierr = DMSwarmGetLocalSize(dms,&nlocal);CHKERRQ(ierr);
	
	ierr = DMSwarmGetField(dms,DMSwarmPICField_coor,&bs,NULL,(void**)&array);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"geoq_fac",NULL,NULL,(void**)&geoq_fac);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"layer",&bs,NULL,(void**)&layer_array);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"strain_fac",NULL,NULL,(void**)&strain_fac);CHKERRQ(ierr);

	for (p=0; p<nlocal; p++) {
		PetscReal cx,cz;
		PetscReal rx,rz,rfac;
		PetscInt i,k;
		
		cx = array[2*p];
		cz = array[2*p+1];
		
		i = (int)(cx/dx_const);
		k = (int)((cz+depth)/dz_const);
		
		
		
		if (i<0 || i>=Nx-1) {printf("estranho i=%d\n",i); exit(1);}
		if (k<0 || k>=Nz-1) {printf("estranho k=%d\n",k); exit(1);}
		
		if (i==Nx-1) i=Nx-2;
		if (k==Nz-1) k=Nz-2;
		
		rx = (cx-i*dx_const)/dx_const;
		rz = (cz-(-depth+k*dz_const))/dz_const;
		
		if (rx<0 || rx>1) {printf("estranho rx=%f\n",rx); exit(1);}
		if (rz<0 || rz>1) {printf("estranho rz=%f\n",rz); exit(1);}
		
		
		
		
		rfac = (1.0-rx)*(1.0-rz);
		qq_rho	[k][i] += rfac*inter_rho[layer_array[p]];
		qq_H	[k][i] += rfac*inter_H[layer_array[p]];
		qq_strain[k][i] += rfac*strain_fac[p];
		qq_cont	[k][i] += rfac;
		
		rfac = (rx)*(1.0-rz);
		qq_rho	[k][i+1] += rfac*inter_rho[layer_array[p]];
		qq_H	[k][i+1] += rfac*inter_H[layer_array[p]];
		qq_strain[k][i+1] += rfac*strain_fac[p];
		qq_cont	[k][i+1] += rfac;
		
		rfac = (1.0-rx)*(rz);
		qq_rho	[k+1][i] += rfac*inter_rho[layer_array[p]];
		qq_H	[k+1][i] += rfac*inter_H[layer_array[p]];
		qq_strain[k+1][i] += rfac*strain_fac[p];
		qq_cont	[k+1][i] += rfac;
		
		rfac = (rx)*(rz);
		qq_rho	[k+1][i+1] += rfac*inter_rho[layer_array[p]];
		qq_H	[k+1][i+1] += rfac*inter_H[layer_array[p]];
		qq_strain[k+1][i+1] += rfac*strain_fac[p];
		qq_cont	[k+1][i+1] += rfac;
	}

	if (visc_const_per_element==0){
		
		if (visc_harmonic_mean==1){
			for (p=0; p<nlocal; p++) {
				PetscReal cx,cz;
				PetscReal rx,rz,rfac;
				PetscInt i,k;
				
				cx = array[2*p];
				cz = array[2*p+1];
				
				i = (int)(cx/dx_const);
				k = (int)((cz+depth)/dz_const);
				
				
				
				if (i<0 || i>=Nx-1) {printf("estranho i=%d\n",i); exit(1);}
				if (k<0 || k>=Nz-1) {printf("estranho k=%d\n",k); exit(1);}
				
				if (i==Nx-1) i=Nx-2;
				if (k==Nz-1) k=Nz-2;
				
				rx = (cx-i*dx_const)/dx_const;
				rz = (cz-(-depth+k*dz_const))/dz_const;
				
				if (rx<0 || rx>1) {printf("estranho rx=%f\n",rx); exit(1);}
				if (rz<0 || rz>1) {printf("estranho rz=%f\n",rz); exit(1);}
				
				rfac = (1.0-rx)*(1.0-rz);
				qq		[k][i] += rfac/geoq_fac[p]; //<--- harmonic
				
				rfac = (rx)*(1.0-rz);
				qq		[k][i+1] += rfac/geoq_fac[p]; //<--- harmonic
				
				rfac = (1.0-rx)*(rz);
				qq		[k+1][i] += rfac/geoq_fac[p]; //<--- harmonic
				
				rfac = (rx)*(rz);
				qq		[k+1][i+1] += rfac/geoq_fac[p]; //<--- harmonic
			}
		}
		else{
			for (p=0; p<nlocal; p++) {
				PetscReal cx,cz;
				PetscReal rx,rz,rfac;
				PetscInt i,k;
				
				cx = array[2*p];
				cz = array[2*p+1];
				
				i = (int)(cx/dx_const);
				k = (int)((cz+depth)/dz_const);
				
				if (i<0 || i>=Nx-1) {printf("estranho i=%d\n",i); exit(1);}
				if (k<0 || k>=Nz-1) {printf("estranho k=%d\n",k); exit(1);}
				
				if (i==Nx-1) i=Nx-2;
				if (k==Nz-1) k=Nz-2;
				
				rx = (cx-i*dx_const)/dx_const;
				rz = (cz-(-depth+k*dz_const))/dz_const;
				
				if (rx<0 || rx>1) {printf("estranho rx=%f\n",rx); exit(1);}
				if (rz<0 || rz>1) {printf("estranho rz=%f\n",rz); exit(1);}
				
				
				rfac = (1.0-rx)*(1.0-rz);
				qq		[k][i] += rfac*geoq_fac[p];
				
				rfac = (rx)*(1.0-rz);
				qq		[k][i+1] += rfac*geoq_fac[p];
				
				rfac = (1.0-rx)*(rz);
				qq		[k+1][i] += rfac*geoq_fac[p];
				
				rfac = (rx)*(rz);
				qq		[k+1][i+1] += rfac*geoq_fac[p];
				
			}
		}
	}
	
	ierr = DMDAVecRestoreArray(da_Thermal,local_geoq,&qq);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(da_Thermal,local_geoq,ADD_VALUES,geoq);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(da_Thermal,local_geoq,ADD_VALUES,geoq);CHKERRQ(ierr);
	
	ierr = DMDAVecRestoreArray(da_Thermal,local_geoq_rho,&qq_rho);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(da_Thermal,local_geoq_rho,ADD_VALUES,geoq_rho);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(da_Thermal,local_geoq_rho,ADD_VALUES,geoq_rho);CHKERRQ(ierr);
	
	ierr = DMDAVecRestoreArray(da_Thermal,local_geoq_H,&qq_H);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(da_Thermal,local_geoq_H,ADD_VALUES,geoq_H);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(da_Thermal,local_geoq_H,ADD_VALUES,geoq_H);CHKERRQ(ierr);
	
	ierr = DMDAVecRestoreArray(da_Thermal,local_geoq_strain,&qq_strain);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(da_Thermal,local_geoq_strain,ADD_VALUES,geoq_strain);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(da_Thermal,local_geoq_strain,ADD_VALUES,geoq_strain);CHKERRQ(ierr);
	
	ierr = DMDAVecRestoreArray(da_Thermal,local_geoq_cont,&qq_cont);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(da_Thermal,local_geoq_cont,ADD_VALUES,geoq_cont);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(da_Thermal,local_geoq_cont,ADD_VALUES,geoq_cont);CHKERRQ(ierr);

	if (periodic_boundary==1){
		ierr = mean_value_periodic_boundary(da_Thermal,geoq_rho,local_geoq_rho,qq_rho,1);
		ierr = mean_value_periodic_boundary(da_Thermal,geoq_H,local_geoq_H,qq_H,1);
		ierr = mean_value_periodic_boundary(da_Thermal,geoq,local_geoq,qq,1);
		ierr = mean_value_periodic_boundary(da_Thermal,geoq_strain,local_geoq_strain,qq_strain,1);
		ierr = mean_value_periodic_boundary(da_Thermal,geoq_cont,local_geoq_cont,qq_cont,1);
	}
	
	//VecPointwiseMax(geoq_cont,geoq_cont,geoqOnes);
	VecPointwiseDivide(geoq,geoq,geoq_cont);
	if (visc_harmonic_mean==1) VecReciprocal(geoq); //<--- harmonic
		
	VecPointwiseDivide(geoq_rho,geoq_rho,geoq_cont);
	VecPointwiseDivide(geoq_H,geoq_H,geoq_cont);
	VecPointwiseDivide(geoq_strain,geoq_strain,geoq_cont);
	//VecPointwiseMax(geoq,geoq,geoqOnes);

	if (visc_const_per_element==1){
		ierr = VecSet(geoq,0.0);CHKERRQ(ierr);
		ierr = VecSet(geoq_cont,0.0);CHKERRQ(ierr);

		ierr = DMGlobalToLocalBegin(da_Thermal,geoq,INSERT_VALUES,local_geoq);
		ierr = DMGlobalToLocalEnd(  da_Thermal,geoq,INSERT_VALUES,local_geoq);
		ierr = DMDAVecGetArray(da_Thermal,local_geoq,&qq);CHKERRQ(ierr);

		
		ierr = DMGlobalToLocalBegin(da_Thermal,geoq_cont,INSERT_VALUES,local_geoq_cont);
		ierr = DMGlobalToLocalEnd(  da_Thermal,geoq_cont,INSERT_VALUES,local_geoq_cont);
		ierr = DMDAVecGetArray(da_Thermal,local_geoq_cont,&qq_cont);CHKERRQ(ierr);
		if (visc_harmonic_mean==1){
			for (p=0; p<nlocal; p++) {
				PetscReal cx,cz;
				PetscInt i,k;
				
				cx = array[2*p];
				cz = array[2*p+1];
				
				i = (int)(cx/dx_const);
				k = (int)((cz+depth)/dz_const);
				
				qq[k][i] += 1.0/geoq_fac[p];
				qq_cont[k][i] += 1.0;
			}
		}
		else {
			for (p=0; p<nlocal; p++) {
				PetscReal cx,cz;
				PetscInt i,k;
				
				cx = array[2*p];
				cz = array[2*p+1];
				
				i = (int)(cx/dx_const);
				k = (int)((cz+depth)/dz_const);
				
				qq[k][i] += geoq_fac[p];
				qq_cont[k][i] += 1.0;
			}
		}


		ierr = DMDAVecRestoreArray(da_Thermal,local_geoq,&qq);CHKERRQ(ierr);
		ierr = DMLocalToGlobalBegin(da_Thermal,local_geoq,ADD_VALUES,geoq);CHKERRQ(ierr);
		ierr = DMLocalToGlobalEnd(da_Thermal,local_geoq,ADD_VALUES,geoq);CHKERRQ(ierr);
		
		ierr = DMDAVecRestoreArray(da_Thermal,local_geoq_cont,&qq_cont);CHKERRQ(ierr);
		ierr = DMLocalToGlobalBegin(da_Thermal,local_geoq_cont,ADD_VALUES,geoq_cont);CHKERRQ(ierr);
		ierr = DMLocalToGlobalEnd(da_Thermal,local_geoq_cont,ADD_VALUES,geoq_cont);CHKERRQ(ierr);

		VecPointwiseDivide(geoq,geoq,geoq_cont);
		if (visc_harmonic_mean==1) VecReciprocal(geoq); //<--- harmonic

	}


	
	ierr = DMSwarmRestoreField(dms,DMSwarmPICField_coor,&bs,NULL,(void**)&array);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms,"geoq_fac",NULL,NULL,(void**)&geoq_fac);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms,"layer",&bs,NULL,(void**)&layer_array);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms,"strain_fac",NULL,NULL,(void**)&strain_fac);CHKERRQ(ierr);
	
	
	ierr = DMGlobalToLocalBegin(da_Thermal,Temper,INSERT_VALUES,local_Temper);
	ierr = DMGlobalToLocalEnd(  da_Thermal,Temper,INSERT_VALUES,local_Temper);
	
	ierr = DMDAVecGetArray(da_Thermal,local_Temper,&TT);CHKERRQ(ierr);
	
	ierr = DMGlobalToLocalBegin(da_Thermal,geoq_rho,INSERT_VALUES,local_geoq_rho);
	ierr = DMGlobalToLocalEnd(  da_Thermal,geoq_rho,INSERT_VALUES,local_geoq_rho);
	
	ierr = DMDAVecGetArray(da_Thermal,local_geoq_rho,&qq_rho);CHKERRQ(ierr);

	
	PetscInt       sx,sz,mmx,mmz;
	
	ierr = DMDAGetCorners(da_Thermal,&sx,&sz,NULL,&mmx,&mmz,NULL);CHKERRQ(ierr);
	
	int k,i;
	
	for (k=sz; k<sz+mmz; k++) {
		for (i=sx; i<sx+mmx; i++) {
			if (qq_rho[k][i]<100.0){ ///air!!! ar!!!
				TT[k][i]=0.0;
			}
		}
	}
	
	ierr = DMDAVecRestoreArray(da_Thermal,local_geoq_rho,&qq_rho);CHKERRQ(ierr);
	
	ierr = DMDAVecRestoreArray(da_Thermal,local_Temper,&TT);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(da_Thermal,local_Temper,INSERT_VALUES,Temper);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(da_Thermal,local_Temper,INSERT_VALUES,Temper);CHKERRQ(ierr);
	
	
	PetscFunctionReturn(0);
	
}
