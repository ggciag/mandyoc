#include <petscdm.h>
#include <petscdmda.h>
#include <petscdmswarm.h>

PetscErrorCode mean_value_periodic_boundary_2d(DM da,Vec F,Vec local_F, PetscScalar **ff,int esc);

extern DM dms;

extern DM da_Thermal;

extern Vec geoq,local_geoq;
extern Vec geoq_rho,local_geoq_rho;
extern Vec geoq_cont,local_geoq_cont;
extern Vec geoq_H,local_geoq_H;
extern Vec geoq_strain, local_geoq_strain;
extern Vec geoq_strain_rate, local_geoq_strain_rate;
extern Vec geoq_kappa, local_geoq_kappa;

extern Vec Temper,local_Temper;

extern long Nx,Ny,Nz;

extern double dx_const;
extern double dy_const;
extern double dz_const;

extern double Lx, Ly, depth;

extern PetscInt visc_harmonic_mean;
extern PetscInt visc_const_per_element;

extern PetscInt sticky_blanket_air;

extern PetscScalar *inter_rho;
extern PetscScalar *inter_H;
extern PetscScalar *conductivity;

extern PetscInt periodic_boundary;

extern PetscInt WITH_SHEAR_H;
extern PetscReal Xi_shear;

extern double visc_MIN;

extern PetscReal air_threshold_density;

extern PetscReal rho0_scaled;

extern PetscReal epsilon_x;

//Shear heating only in 2D!!! Must be implemented for 3D!!!

PetscErrorCode Swarm2Mesh_2d(){

	PetscErrorCode ierr;
	PetscScalar             **qq,**qq_cont,**qq_rho,**TT,**qq_H,**qq_strain,**qq_kappa;
	PetscScalar				**qq_strain_rate;

	ierr = VecSet(geoq,0.0);CHKERRQ(ierr);
	ierr = VecSet(geoq_rho,0.0);CHKERRQ(ierr);
	ierr = VecSet(geoq_cont,0.0);CHKERRQ(ierr);
	ierr = VecSet(geoq_H,0.0);CHKERRQ(ierr);
	ierr = VecSet(geoq_strain,0.0);CHKERRQ(ierr);
	ierr = VecSet(geoq_strain_rate,0.0);CHKERRQ(ierr);
	

	ierr = VecZeroEntries(local_geoq);CHKERRQ(ierr);
	ierr = VecZeroEntries(local_geoq_rho);CHKERRQ(ierr);
	ierr = VecZeroEntries(local_geoq_H);CHKERRQ(ierr);
	ierr = VecZeroEntries(local_geoq_strain);CHKERRQ(ierr);
	ierr = VecZeroEntries(local_geoq_strain_rate);CHKERRQ(ierr);
	

	ierr = DMGlobalToLocalBegin(da_Thermal,geoq,INSERT_VALUES,local_geoq);
	ierr = DMGlobalToLocalEnd(  da_Thermal,geoq,INSERT_VALUES,local_geoq);

	ierr = DMGlobalToLocalBegin(da_Thermal,geoq_rho,INSERT_VALUES,local_geoq_rho);
	ierr = DMGlobalToLocalEnd(  da_Thermal,geoq_rho,INSERT_VALUES,local_geoq_rho);

	ierr = DMGlobalToLocalBegin(da_Thermal,geoq_H,INSERT_VALUES,local_geoq_H);
	ierr = DMGlobalToLocalEnd(  da_Thermal,geoq_H,INSERT_VALUES,local_geoq_H);

	ierr = DMGlobalToLocalBegin(da_Thermal,geoq_strain,INSERT_VALUES,local_geoq_strain);
	ierr = DMGlobalToLocalEnd(  da_Thermal,geoq_strain,INSERT_VALUES,local_geoq_strain);

	ierr = DMGlobalToLocalBegin(da_Thermal,geoq_strain_rate,INSERT_VALUES,local_geoq_strain_rate);
	ierr = DMGlobalToLocalEnd(  da_Thermal,geoq_strain_rate,INSERT_VALUES,local_geoq_strain_rate);

	ierr = DMDAVecGetArray(da_Thermal,local_geoq,&qq);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da_Thermal,local_geoq_rho,&qq_rho);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da_Thermal,local_geoq_H,&qq_H);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da_Thermal,local_geoq_strain,&qq_strain);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da_Thermal,local_geoq_strain_rate,&qq_strain_rate);CHKERRQ(ierr);


	ierr = DMGlobalToLocalBegin(da_Thermal,geoq_cont,INSERT_VALUES,local_geoq_cont);
	ierr = DMGlobalToLocalEnd(  da_Thermal,geoq_cont,INSERT_VALUES,local_geoq_cont);

	ierr = DMDAVecGetArray(da_Thermal,local_geoq_cont,&qq_cont);CHKERRQ(ierr);



	PetscInt nlocal,bs,p;

	PetscReal *array;
	PetscReal *geoq_fac;
	PetscReal *strain_fac;
	PetscReal *strain_rate_fac;
	PetscInt *layer_array;

	PetscReal inter_H_aux;

	ierr = DMSwarmGetLocalSize(dms,&nlocal);CHKERRQ(ierr);

	ierr = DMSwarmGetField(dms,DMSwarmPICField_coor,&bs,NULL,(void**)&array);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"geoq_fac",NULL,NULL,(void**)&geoq_fac);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"layer",&bs,NULL,(void**)&layer_array);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"strain_fac",NULL,NULL,(void**)&strain_fac);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"strain_rate_fac",NULL,NULL,(void**)&strain_rate_fac);CHKERRQ(ierr);

	for (p=0; p<nlocal; p++) {
		PetscReal cx,cz;
		PetscReal rx,rz,rfac;
		PetscInt i,k;

		cx = array[2*p];
		if (cx>=Lx) {
			printf("Swarm2Mesh - outside: cx=%lf>=%lf\n",cx,Lx);
			cx=Lx-epsilon_x;
		}
		if (cx<=0.0) {
			printf("Swarm2Mesh - outside: cx=%lf<=0.0\n",cx);
			cx=epsilon_x;
		}
		cz = array[2*p+1];
		if (cz>=0){
			printf("Swarm2Mesh - outside: cz=%lf>=0.0\n",cz);
			cz=-epsilon_x;
		}
		if (cz<=-depth){
			printf("Swarm2Mesh - outside: cz=%lf<=-%lf\n",cz,depth);
			cz=-depth+epsilon_x;
		}

		i = (int)(cx/dx_const);
		k = (int)((cz+depth)/dz_const);



		if (i<0 || i>=Nx-1) {printf("weird i=%d cx=%f , Swarm2Mesh\n",i,cx); exit(1);}
		if (k<0 || k>=Nz-1) {printf("weird k=%d cz=%f , Swarm2Mesh\n",k,cz); exit(1);}

		if (i==Nx-1) i=Nx-2;
		if (k==Nz-1) k=Nz-2;

		rx = (cx-i*dx_const)/dx_const;
		rz = (cz-(-depth+k*dz_const))/dz_const;

		if (rx<0 || rx>1) {printf("weird rx=%f , Swarm2Mesh\n",rx); exit(1);}
		if (rz<0 || rz>1) {printf("weird rz=%f , Swarm2Mesh\n",rz); exit(1);}

		if (WITH_SHEAR_H == 1){
			inter_H_aux = inter_H[layer_array[p]] + Xi_shear*4*geoq_fac[p]*strain_rate_fac[p]*strain_rate_fac[p]/inter_rho[layer_array[p]];
		}
		else {
			inter_H_aux = inter_H[layer_array[p]];
		}


		rfac = (1.0-rx)*(1.0-rz);
		qq_rho	[k][i] += rfac*inter_rho[layer_array[p]];
		qq_H	[k][i] += rfac*inter_H_aux;
		qq_strain[k][i] += rfac*strain_fac[p];
		qq_strain_rate[k][i] += rfac*strain_rate_fac[p];
		qq_cont	[k][i] += rfac;

		rfac = (rx)*(1.0-rz);
		qq_rho	[k][i+1] += rfac*inter_rho[layer_array[p]];
		qq_H	[k][i+1] += rfac*inter_H_aux;
		qq_strain[k][i+1] += rfac*strain_fac[p];
		qq_strain_rate[k][i+1] += rfac*strain_rate_fac[p];
		qq_cont	[k][i+1] += rfac;

		rfac = (1.0-rx)*(rz);
		qq_rho	[k+1][i] += rfac*inter_rho[layer_array[p]];
		qq_H	[k+1][i] += rfac*inter_H_aux;
		qq_strain[k+1][i] += rfac*strain_fac[p];
		qq_strain_rate[k+1][i] += rfac*strain_rate_fac[p];
		qq_cont	[k+1][i] += rfac;

		rfac = (rx)*(rz);
		qq_rho	[k+1][i+1] += rfac*inter_rho[layer_array[p]];
		qq_H	[k+1][i+1] += rfac*inter_H_aux;
		qq_strain[k+1][i+1] += rfac*strain_fac[p];
		qq_strain_rate[k+1][i+1] += rfac*strain_rate_fac[p];
		qq_cont	[k+1][i+1] += rfac;
	}

	if (visc_const_per_element==0){

		PetscPrintf(PETSC_COMM_WORLD,"\n\nThe option \nviscosity_per_element = variable \nis no longer available.\n");
		exit(-2);

		if (visc_harmonic_mean==1){
			for (p=0; p<nlocal; p++) {
				PetscReal cx,cz;
				PetscReal rx,rz,rfac;
				PetscInt i,k;

				cx = array[2*p];
				if (cx>=Lx) {
					printf("Swarm2Mesh visc_harmonic_mean 1 - outside: cx=%lf>=%lf\n",cx,Lx);
					cx=Lx-epsilon_x;
				}
				if (cx<=0.0) {
					printf("Swarm2Mesh visc_harmonic_mean 1 - outside: cx=%lf<=0.0\n",cx);
					cx=epsilon_x;
				}
				cz = array[2*p+1];
				if (cz>=0){
					printf("Swarm2Mesh visc_harmonic_mean 1 - outside: cz=%lf>=0.0\n",cz);
					cz=-epsilon_x;
				}
				if (cz<=-depth){
					printf("Swarm2Mesh visc_harmonic_mean 1 - outside: cz=%lf<=-%lf\n",cz,depth);
					cz=-depth+epsilon_x;
				}

				i = (int)(cx/dx_const);
				k = (int)((cz+depth)/dz_const);



				if (i<0 || i>=Nx-1) {printf("weird i=%d cx=%f\n",i,cx); exit(1);}
				if (k<0 || k>=Nz-1) {printf("weird k=%d cz=%f\n",k,cz); exit(1);}

				if (i==Nx-1) i=Nx-2;
				if (k==Nz-1) k=Nz-2;

				rx = (cx-i*dx_const)/dx_const;
				rz = (cz-(-depth+k*dz_const))/dz_const;

				if (rx<0 || rx>1) {printf("weird rx=%f\n",rx); exit(1);}
				if (rz<0 || rz>1) {printf("weird rz=%f\n",rz); exit(1);}

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
				if (cx>=Lx) {
					printf("Swarm2Mesh visc_harmonic_mean 0 - outside: cx=%lf>=%lf\n",cx,Lx);
					cx=Lx-epsilon_x;
				}
				if (cx<=0.0) {
					printf("Swarm2Mesh visc_harmonic_mean 0 - outside: cx=%lf<=0.0\n",cx);
					cx=epsilon_x;
				}
				cz = array[2*p+1];
				if (cz>=0){
					printf("Swarm2Mesh visc_harmonic_mean 0 - outside: cz=%lf>=0.0\n",cz);
					cz=-epsilon_x;
				}
				if (cz<=-depth){
					printf("Swarm2Mesh visc_harmonic_mean 0 - outside: cz=%lf<=-%lf\n",cz,depth);
					cz=-depth+epsilon_x;
				}

				i = (int)(cx/dx_const);
				k = (int)((cz+depth)/dz_const);

				if (i<0 || i>=Nx-1) {printf("weird i=%d\n",i); exit(1);}
				if (k<0 || k>=Nz-1) {printf("weird k=%d\n",k); exit(1);}

				if (i==Nx-1) i=Nx-2;
				if (k==Nz-1) k=Nz-2;

				rx = (cx-i*dx_const)/dx_const;
				rz = (cz-(-depth+k*dz_const))/dz_const;

				if (rx<0 || rx>1) {printf("weird rx=%f\n",rx); exit(1);}
				if (rz<0 || rz>1) {printf("weird rz=%f\n",rz); exit(1);}


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

	ierr = DMDAVecRestoreArray(da_Thermal,local_geoq_strain_rate,&qq_strain_rate);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(da_Thermal,local_geoq_strain_rate,ADD_VALUES,geoq_strain_rate);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(da_Thermal,local_geoq_strain_rate,ADD_VALUES,geoq_strain_rate);CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(da_Thermal,local_geoq_cont,&qq_cont);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(da_Thermal,local_geoq_cont,ADD_VALUES,geoq_cont);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(da_Thermal,local_geoq_cont,ADD_VALUES,geoq_cont);CHKERRQ(ierr);


	if (periodic_boundary==1){
		ierr = mean_value_periodic_boundary_2d(da_Thermal,geoq_rho,local_geoq_rho,qq_rho,1);
		ierr = mean_value_periodic_boundary_2d(da_Thermal,geoq_H,local_geoq_H,qq_H,1);
		ierr = mean_value_periodic_boundary_2d(da_Thermal,geoq,local_geoq,qq,1);
		ierr = mean_value_periodic_boundary_2d(da_Thermal,geoq_strain,local_geoq_strain,qq_strain,1);
		ierr = mean_value_periodic_boundary_2d(da_Thermal,geoq_strain_rate,local_geoq_strain_rate,qq_strain_rate,1);
		ierr = mean_value_periodic_boundary_2d(da_Thermal,geoq_cont,local_geoq_cont,qq_cont,1);
	}


	VecPointwiseDivide(geoq,geoq,geoq_cont);
	if (visc_harmonic_mean==1) VecReciprocal(geoq); //<--- harmonic

	VecPointwiseDivide(geoq_rho,geoq_rho,geoq_cont);
	VecPointwiseDivide(geoq_H,geoq_H,geoq_cont);
	VecPointwiseDivide(geoq_strain,geoq_strain,geoq_cont);
	VecPointwiseDivide(geoq_strain_rate,geoq_strain_rate,geoq_cont);

	if (visc_const_per_element==1){
		ierr = VecSet(geoq,0.0);CHKERRQ(ierr);
		ierr = VecSet(geoq_cont,0.0);CHKERRQ(ierr);
		ierr = VecSet(geoq_kappa,0.0);CHKERRQ(ierr);

		ierr = DMGlobalToLocalBegin(da_Thermal,geoq,INSERT_VALUES,local_geoq);
		ierr = DMGlobalToLocalEnd(  da_Thermal,geoq,INSERT_VALUES,local_geoq);
		ierr = DMDAVecGetArray(da_Thermal,local_geoq,&qq);CHKERRQ(ierr);

		ierr = DMGlobalToLocalBegin(da_Thermal,geoq_kappa,INSERT_VALUES,local_geoq_kappa);
		ierr = DMGlobalToLocalEnd(  da_Thermal,geoq_kappa,INSERT_VALUES,local_geoq_kappa);
		ierr = DMDAVecGetArray(da_Thermal,local_geoq_kappa,&qq_kappa);CHKERRQ(ierr);

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
				qq_kappa[k][i] += 1.0/conductivity[layer_array[p]];
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
				qq_kappa[k][i] += conductivity[layer_array[p]];
				qq_cont[k][i] += 1.0;
			}
		}


		ierr = DMDAVecRestoreArray(da_Thermal,local_geoq,&qq);CHKERRQ(ierr);
		ierr = DMLocalToGlobalBegin(da_Thermal,local_geoq,ADD_VALUES,geoq);CHKERRQ(ierr);
		ierr = DMLocalToGlobalEnd(da_Thermal,local_geoq,ADD_VALUES,geoq);CHKERRQ(ierr);

		ierr = DMDAVecRestoreArray(da_Thermal,local_geoq_cont,&qq_cont);CHKERRQ(ierr);
		ierr = DMLocalToGlobalBegin(da_Thermal,local_geoq_cont,ADD_VALUES,geoq_cont);CHKERRQ(ierr);
		ierr = DMLocalToGlobalEnd(da_Thermal,local_geoq_cont,ADD_VALUES,geoq_cont);CHKERRQ(ierr);

		ierr = DMDAVecRestoreArray(da_Thermal,local_geoq_kappa,&qq_kappa);CHKERRQ(ierr);
		ierr = DMLocalToGlobalBegin(da_Thermal,local_geoq_kappa,ADD_VALUES,geoq_kappa);CHKERRQ(ierr);
		ierr = DMLocalToGlobalEnd(da_Thermal,local_geoq_kappa,ADD_VALUES,geoq_kappa);CHKERRQ(ierr);

		VecPointwiseDivide(geoq,geoq,geoq_cont);
		VecPointwiseDivide(geoq_kappa,geoq_kappa,geoq_cont);
		if (visc_harmonic_mean==1) VecReciprocal(geoq); //<--- harmonic
		if (visc_harmonic_mean==1) VecReciprocal(geoq_kappa); //<--- harmonic

	}



	ierr = DMSwarmRestoreField(dms,DMSwarmPICField_coor,&bs,NULL,(void**)&array);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms,"geoq_fac",NULL,NULL,(void**)&geoq_fac);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms,"layer",&bs,NULL,(void**)&layer_array);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms,"strain_fac",NULL,NULL,(void**)&strain_fac);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms,"strain_rate_fac",NULL,NULL,(void**)&strain_rate_fac);CHKERRQ(ierr);


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
			if (qq_rho[k][i]<air_threshold_density){ //check: force temperature zero for low density material
				TT[k][i]=0.0;
			}
		}
	}

	ierr = DMDAVecRestoreArray(da_Thermal,local_geoq_rho,&qq_rho);CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(da_Thermal,local_Temper,&TT);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(da_Thermal,local_Temper,INSERT_VALUES,Temper);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(da_Thermal,local_Temper,INSERT_VALUES,Temper);CHKERRQ(ierr);

	if (sticky_blanket_air==1){
		if (visc_const_per_element==1){

			ierr = DMGlobalToLocalBegin(da_Thermal,geoq,INSERT_VALUES,local_geoq);
			ierr = DMGlobalToLocalEnd(  da_Thermal,geoq,INSERT_VALUES,local_geoq);
			ierr = DMDAVecGetArray(da_Thermal,local_geoq,&qq);CHKERRQ(ierr);

			ierr = DMGlobalToLocalBegin(da_Thermal,geoq_rho,INSERT_VALUES,local_geoq_rho);
			ierr = DMGlobalToLocalEnd(  da_Thermal,geoq_rho,INSERT_VALUES,local_geoq_rho);
			ierr = DMDAVecGetArray(da_Thermal,local_geoq_rho,&qq_rho);CHKERRQ(ierr);
			PetscReal rho_min=5.0/rho0_scaled,rho_max=2700.0*0.8/rho0_scaled, rho_mean; //check: in the future, make this free to the user modify
			for (k=sz; k<sz+mmz-1; k++) {
				for (i=sx; i<sx+mmx-1; i++) {
					rho_mean = (qq_rho[k][i] + qq_rho[k][i+1] + qq_rho[k+1][i] + qq_rho[k+1][i+1])/4.0;
					if (rho_mean>rho_min && rho_mean<rho_max){
						if (qq[k][i]<visc_MIN*50) qq[k][i]=visc_MIN*50;
						if (qq[k][i+1]<visc_MIN*50) qq[k][i+1]=visc_MIN*50;
					}
				}
			}

			ierr = DMDAVecRestoreArray(da_Thermal,local_geoq,&qq);CHKERRQ(ierr);
			ierr = DMLocalToGlobalBegin(da_Thermal,local_geoq,INSERT_VALUES,geoq);CHKERRQ(ierr);
			ierr = DMLocalToGlobalEnd(da_Thermal,local_geoq,INSERT_VALUES,geoq);CHKERRQ(ierr);


			ierr = DMDAVecRestoreArray(da_Thermal,local_geoq_rho,&qq_rho);CHKERRQ(ierr);

		}
	}


	PetscFunctionReturn(0);

}

PetscErrorCode Swarm2Mesh_3d(){ ///!!!  geoq_kappa for 3D models have to be implemented!!!

	PetscErrorCode ierr;
	PetscScalar             ***qq,***qq_cont,***qq_rho,***TT,***qq_H,***qq_strain;
	PetscScalar				***qq_strain_rate;

	ierr = VecSet(geoq,0.0);CHKERRQ(ierr);
	ierr = VecSet(geoq_rho,0.0);CHKERRQ(ierr);
	ierr = VecSet(geoq_cont,0.0);CHKERRQ(ierr);
	ierr = VecSet(geoq_H,0.0);CHKERRQ(ierr);
	ierr = VecSet(geoq_strain,0.0);CHKERRQ(ierr);
	ierr = VecSet(geoq_strain_rate,0.0);CHKERRQ(ierr);

	ierr = VecZeroEntries(local_geoq);CHKERRQ(ierr);
	ierr = VecZeroEntries(local_geoq_rho);CHKERRQ(ierr);
	ierr = VecZeroEntries(local_geoq_H);CHKERRQ(ierr);
	ierr = VecZeroEntries(local_geoq_strain);CHKERRQ(ierr);
	ierr = VecZeroEntries(local_geoq_strain_rate);CHKERRQ(ierr);

	ierr = DMGlobalToLocalBegin(da_Thermal,geoq,INSERT_VALUES,local_geoq);
	ierr = DMGlobalToLocalEnd(  da_Thermal,geoq,INSERT_VALUES,local_geoq);

	ierr = DMGlobalToLocalBegin(da_Thermal,geoq_rho,INSERT_VALUES,local_geoq_rho);
	ierr = DMGlobalToLocalEnd(  da_Thermal,geoq_rho,INSERT_VALUES,local_geoq_rho);

	ierr = DMGlobalToLocalBegin(da_Thermal,geoq_H,INSERT_VALUES,local_geoq_H);
	ierr = DMGlobalToLocalEnd(  da_Thermal,geoq_H,INSERT_VALUES,local_geoq_H);

	ierr = DMGlobalToLocalBegin(da_Thermal,geoq_strain,INSERT_VALUES,local_geoq_strain);
	ierr = DMGlobalToLocalEnd(  da_Thermal,geoq_strain,INSERT_VALUES,local_geoq_strain);

	ierr = DMGlobalToLocalBegin(da_Thermal,geoq_strain_rate,INSERT_VALUES,local_geoq_strain_rate);
	ierr = DMGlobalToLocalEnd(  da_Thermal,geoq_strain_rate,INSERT_VALUES,local_geoq_strain_rate);


	ierr = DMDAVecGetArray(da_Thermal,local_geoq,&qq);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da_Thermal,local_geoq_rho,&qq_rho);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da_Thermal,local_geoq_H,&qq_H);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da_Thermal,local_geoq_strain,&qq_strain);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da_Thermal,local_geoq_strain_rate,&qq_strain_rate);CHKERRQ(ierr);


	ierr = DMGlobalToLocalBegin(da_Thermal,geoq_cont,INSERT_VALUES,local_geoq_cont);
	ierr = DMGlobalToLocalEnd(  da_Thermal,geoq_cont,INSERT_VALUES,local_geoq_cont);

	ierr = DMDAVecGetArray(da_Thermal,local_geoq_cont,&qq_cont);CHKERRQ(ierr);



	PetscInt nlocal,bs,p;

	PetscReal *array;
	PetscReal *geoq_fac;
	PetscReal *strain_fac;
	PetscReal *strain_rate_fac;
	PetscInt *layer_array;

	ierr = DMSwarmGetLocalSize(dms,&nlocal);CHKERRQ(ierr);

	ierr = DMSwarmGetField(dms,DMSwarmPICField_coor,&bs,NULL,(void**)&array);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"geoq_fac",NULL,NULL,(void**)&geoq_fac);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"layer",&bs,NULL,(void**)&layer_array);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"strain_fac",NULL,NULL,(void**)&strain_fac);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"strain_rate_fac",NULL,NULL,(void**)&strain_rate_fac);CHKERRQ(ierr);

	if (visc_harmonic_mean==1){
		for (p=0; p<nlocal; p++) {
			PetscReal cx,cy,cz;
			PetscReal rx,ry,rz,rfac;
			PetscInt i,j,k;

			cx = array[3*p];
			cy = array[3*p+1];
			cz = array[3*p+2];

			if (cx>=Lx) {
				printf("moveSwarm in 3D - outside: cx=%lf>=%lf\n",cx,Lx);
				cx=Lx-epsilon_x;
			}
			if (cx<=0.0) {
				printf("moveSwarm in 3D - outside: cx=%lf<=0.0\n",cx);
				cx=epsilon_x;
			}
			if (cy>=Ly) {
				printf("moveSwarm in 3D - outside: cy=%lf>=%lf\n",cy,Ly);
				cy=Ly-epsilon_x;
			}
			if (cy<=0.0) {
				printf("moveSwarm in 3D - outside: cy=%lf<=0.0\n",cy);
				cy=epsilon_x;
			}
			if (cz>=0){
				printf("moveSwarm in 3D - outside: cz=%lf>=0.0\n",cz);
				cz=-epsilon_x;
			}
			if (cz<=-depth){
				printf("moveSwarm in 3D - outside: cz=%lf<=-%lf\n",cz,depth);
				cz=-depth+epsilon_x;
			}

			i = (int)(cx/dx_const);
			j = (int)(cy/dy_const);
			k = (int)((cz+depth)/dz_const);



			if (i<0 || i>=Nx-1) {printf("estranho i=%d\n",i); exit(1);}
			if (j<0 || j>=Ny-1) {printf("estranho j=%d\n",j); exit(1);}
			if (k<0 || k>=Nz-1) {printf("estranho k=%d\n",k); exit(1);}

			if (i==Nx-1) i=Nx-2;
			if (j==Ny-1) j=Ny-2;
			if (k==Nz-1) k=Nz-2;

			rx = (cx-i*dx_const)/dx_const;
			ry = (cy-j*dy_const)/dy_const;
			rz = (cz-(-depth+k*dz_const))/dz_const;

			if (rx<0 || rx>1) {printf("estranho rx=%f\n",rx); exit(1);}
			if (ry<0 || ry>1) {printf("estranho ry=%f\n",ry); exit(1);}
			if (rz<0 || rz>1) {printf("estranho rz=%f\n",rz); exit(1);}


			rfac = (1.0-rx)*(1.0-ry)*(1.0-rz);
			qq		[k][j][i] += rfac/geoq_fac[p]; //<--- harmonic
			qq_rho	[k][j][i] += rfac*inter_rho[layer_array[p]];
			qq_H	[k][j][i] += rfac*inter_H[layer_array[p]];
			qq_strain[k][j][i] += rfac*strain_fac[p];
			qq_strain_rate[k][j][i] += rfac*strain_rate_fac[p];
			qq_cont	[k][j][i] += rfac;

			rfac = (rx)*(1.0-ry)*(1.0-rz);
			qq		[k][j][i+1] += rfac/geoq_fac[p]; //<--- harmonic
			qq_rho	[k][j][i+1] += rfac*inter_rho[layer_array[p]];
			qq_H	[k][j][i+1] += rfac*inter_H[layer_array[p]];
			qq_strain[k][j][i+1] += rfac*strain_fac[p];
			qq_strain_rate[k][j][i+1] += rfac*strain_rate_fac[p];
			qq_cont	[k][j][i+1] += rfac;

			rfac = (1.0-rx)*(ry)*(1.0-rz);
			qq		[k][j+1][i] += rfac/geoq_fac[p]; //<--- harmonic
			qq_rho	[k][j+1][i] += rfac*inter_rho[layer_array[p]];
			qq_H	[k][j+1][i] += rfac*inter_H[layer_array[p]];
			qq_strain[k][j+1][i] += rfac*strain_fac[p];
			qq_strain_rate[k][j+1][i] += rfac*strain_rate_fac[p];
			qq_cont	[k][j+1][i] += rfac;

			rfac = (rx)*(ry)*(1.0-rz);
			qq		[k][j+1][i+1] += rfac/geoq_fac[p]; //<--- harmonic
			qq_rho	[k][j+1][i+1] += rfac*inter_rho[layer_array[p]];
			qq_H	[k][j+1][i+1] += rfac*inter_H[layer_array[p]];
			qq_strain[k][j+1][i+1] += rfac*strain_fac[p];
			qq_strain_rate[k][j+1][i+1] += rfac*strain_rate_fac[p];
			qq_cont	[k][j+1][i+1] += rfac;



			rfac = (1.0-rx)*(1.0-ry)*(rz);
			qq		[k+1][j][i] += rfac/geoq_fac[p]; //<--- harmonic
			qq_rho	[k+1][j][i] += rfac*inter_rho[layer_array[p]];
			qq_H	[k+1][j][i] += rfac*inter_H[layer_array[p]];
			qq_strain[k+1][j][i] += rfac*strain_fac[p];
			qq_strain_rate[k+1][j][i] += rfac*strain_rate_fac[p];
			qq_cont	[k+1][j][i] += rfac;

			rfac = (rx)*(1.0-ry)*(rz);
			qq		[k+1][j][i+1] += rfac/geoq_fac[p]; //<--- harmonic
			qq_rho	[k+1][j][i+1] += rfac*inter_rho[layer_array[p]];
			qq_H	[k+1][j][i+1] += rfac*inter_H[layer_array[p]];
			qq_strain[k+1][j][i+1] += rfac*strain_fac[p];
			qq_strain_rate[k+1][j][i+1] += rfac*strain_rate_fac[p];
			qq_cont	[k+1][j][i+1] += rfac;

			rfac = (1.0-rx)*(ry)*(rz);
			qq		[k+1][j+1][i] += rfac/geoq_fac[p]; //<--- harmonic
			qq_rho	[k+1][j+1][i] += rfac*inter_rho[layer_array[p]];
			qq_H	[k+1][j+1][i] += rfac*inter_H[layer_array[p]];
			qq_strain[k+1][j+1][i] += rfac*strain_fac[p];
			qq_strain_rate[k+1][j+1][i] += rfac*strain_rate_fac[p];
			qq_cont	[k+1][j+1][i] += rfac;

			rfac = (rx)*(ry)*(rz);
			qq		[k+1][j+1][i+1] += rfac/geoq_fac[p]; //<--- harmonic
			qq_rho	[k+1][j+1][i+1] += rfac*inter_rho[layer_array[p]];
			qq_H	[k+1][j+1][i+1] += rfac*inter_H[layer_array[p]];
			qq_strain[k+1][j+1][i+1] += rfac*strain_fac[p];
			qq_strain_rate[k+1][j+1][i+1] += rfac*strain_rate_fac[p];
			qq_cont	[k+1][j+1][i+1] += rfac;
		}
	}
	else{
		for (p=0; p<nlocal; p++) {
			PetscReal cx,cy,cz;
			PetscReal rx,ry,rz,rfac;
			PetscInt i,j,k;

			cx = array[3*p];
			cy = array[3*p+1];
			cz = array[3*p+2];

			if (cx>=Lx) {
				printf("moveSwarm in 3D - outside: cx=%lf>=%lf\n",cx,Lx);
				cx=Lx-epsilon_x;
			}
			if (cx<=0.0) {
				printf("moveSwarm in 3D - outside: cx=%lf<=0.0\n",cx);
				cx=epsilon_x;
			}
			if (cy>=Ly) {
				printf("moveSwarm in 3D - outside: cy=%lf>=%lf\n",cy,Ly);
				cy=Ly-epsilon_x;
			}
			if (cy<=0.0) {
				printf("moveSwarm in 3D - outside: cy=%lf<=0.0\n",cy);
				cy=epsilon_x;
			}
			if (cz>=0){
				printf("moveSwarm in 3D - outside: cz=%lf>=0.0\n",cz);
				cz=-epsilon_x;
			}
			if (cz<=-depth){
				printf("moveSwarm in 3D - outside: cz=%lf<=-%lf\n",cz,depth);
				cz=-depth+epsilon_x;
			}

			i = (int)(cx/dx_const);
			j = (int)(cy/dy_const);
			k = (int)((cz+depth)/dz_const);



			if (i<0 || i>=Nx-1) {printf("estranho i=%d\n",i); exit(1);}
			if (j<0 || j>=Ny-1) {printf("estranho j=%d\n",j); exit(1);}
			if (k<0 || k>=Nz-1) {printf("estranho k=%d\n",k); exit(1);}

			if (i==Nx-1) i=Nx-2;
			if (j==Ny-1) j=Ny-2;
			if (k==Nz-1) k=Nz-2;

			rx = (cx-i*dx_const)/dx_const;
			ry = (cy-j*dy_const)/dy_const;
			rz = (cz-(-depth+k*dz_const))/dz_const;

			if (rx<0 || rx>1) {printf("estranho rx=%f\n",rx); exit(1);}
			if (ry<0 || ry>1) {printf("estranho ry=%f\n",ry); exit(1);}
			if (rz<0 || rz>1) {printf("estranho rz=%f\n",rz); exit(1);}


			rfac = (1.0-rx)*(1.0-ry)*(1.0-rz);
			qq		[k][j][i] += rfac*geoq_fac[p];
			qq_rho	[k][j][i] += rfac*inter_rho[layer_array[p]];
			qq_H	[k][j][i] += rfac*inter_H[layer_array[p]];
			qq_strain[k][j][i] += rfac*strain_fac[p];
			qq_strain_rate[k][j][i] += rfac*strain_rate_fac[p];
			qq_cont	[k][j][i] += rfac;

			rfac = (rx)*(1.0-ry)*(1.0-rz);
			qq		[k][j][i+1] += rfac*geoq_fac[p];
			qq_rho	[k][j][i+1] += rfac*inter_rho[layer_array[p]];
			qq_H	[k][j][i+1] += rfac*inter_H[layer_array[p]];
			qq_strain[k][j][i+1] += rfac*strain_fac[p];
			qq_strain_rate[k][j][i+1] += rfac*strain_rate_fac[p];
			qq_cont	[k][j][i+1] += rfac;

			rfac = (1.0-rx)*(ry)*(1.0-rz);
			qq		[k][j+1][i] += rfac*geoq_fac[p];
			qq_rho	[k][j+1][i] += rfac*inter_rho[layer_array[p]];
			qq_H	[k][j+1][i] += rfac*inter_H[layer_array[p]];
			qq_strain[k][j+1][i] += rfac*strain_fac[p];
			qq_strain_rate[k][j+1][i] += rfac*strain_rate_fac[p];
			qq_cont	[k][j+1][i] += rfac;

			rfac = (rx)*(ry)*(1.0-rz);
			qq		[k][j+1][i+1] += rfac*geoq_fac[p];
			qq_rho	[k][j+1][i+1] += rfac*inter_rho[layer_array[p]];
			qq_H	[k][j+1][i+1] += rfac*inter_H[layer_array[p]];
			qq_strain[k][j+1][i+1] += rfac*strain_fac[p];
			qq_strain_rate[k][j+1][i+1] += rfac*strain_rate_fac[p];
			qq_cont	[k][j+1][i+1] += rfac;



			rfac = (1.0-rx)*(1.0-ry)*(rz);
			qq		[k+1][j][i] += rfac*geoq_fac[p];
			qq_rho	[k+1][j][i] += rfac*inter_rho[layer_array[p]];
			qq_H	[k+1][j][i] += rfac*inter_H[layer_array[p]];
			qq_strain[k+1][j][i] += rfac*strain_fac[p];
			qq_strain_rate[k+1][j][i] += rfac*strain_rate_fac[p];
			qq_cont	[k+1][j][i] += rfac;

			rfac = (rx)*(1.0-ry)*(rz);
			qq		[k+1][j][i+1] += rfac*geoq_fac[p];
			qq_rho	[k+1][j][i+1] += rfac*inter_rho[layer_array[p]];
			qq_H	[k+1][j][i+1] += rfac*inter_H[layer_array[p]];
			qq_strain[k+1][j][i+1] += rfac*strain_fac[p];
			qq_strain_rate[k+1][j][i+1] += rfac*strain_rate_fac[p];
			qq_cont	[k+1][j][i+1] += rfac;

			rfac = (1.0-rx)*(ry)*(rz);
			qq		[k+1][j+1][i] += rfac*geoq_fac[p];
			qq_rho	[k+1][j+1][i] += rfac*inter_rho[layer_array[p]];
			qq_H	[k+1][j+1][i] += rfac*inter_H[layer_array[p]];
			qq_strain[k+1][j+1][i] += rfac*strain_fac[p];
			qq_strain_rate[k+1][j+1][i] += rfac*strain_rate_fac[p];
			qq_cont	[k+1][j+1][i] += rfac;

			rfac = (rx)*(ry)*(rz);
			qq		[k+1][j+1][i+1] += rfac*geoq_fac[p];
			qq_rho	[k+1][j+1][i+1] += rfac*inter_rho[layer_array[p]];
			qq_H	[k+1][j+1][i+1] += rfac*inter_H[layer_array[p]];
			qq_strain[k+1][j+1][i+1] += rfac*strain_fac[p];
			qq_strain_rate[k+1][j+1][i+1] += rfac*strain_rate_fac[p];
			qq_cont	[k+1][j+1][i+1] += rfac;
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

	ierr = DMDAVecRestoreArray(da_Thermal,local_geoq_strain_rate,&qq_strain_rate);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(da_Thermal,local_geoq_strain_rate,ADD_VALUES,geoq_strain_rate);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(da_Thermal,local_geoq_strain_rate,ADD_VALUES,geoq_strain_rate);CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(da_Thermal,local_geoq_cont,&qq_cont);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(da_Thermal,local_geoq_cont,ADD_VALUES,geoq_cont);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(da_Thermal,local_geoq_cont,ADD_VALUES,geoq_cont);CHKERRQ(ierr);

	//VecPointwiseMax(geoq_cont,geoq_cont,geoqOnes);
	VecPointwiseDivide(geoq,geoq,geoq_cont);
	if (visc_harmonic_mean==1) VecReciprocal(geoq); //<--- harmonic

	VecPointwiseDivide(geoq_rho,geoq_rho,geoq_cont);
	VecPointwiseDivide(geoq_H,geoq_H,geoq_cont);
	VecPointwiseDivide(geoq_strain,geoq_strain,geoq_cont);
	VecPointwiseDivide(geoq_strain_rate,geoq_strain_rate,geoq_cont);
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
				PetscReal cx,cy,cz;
				PetscInt i,j,k;

				cx = array[3*p];
				cy = array[3*p+1];
				cz = array[3*p+2];

				if (cx>=Lx) {
					printf("moveSwarm in 3D - outside: cx=%lf>=%lf\n",cx,Lx);
					cx=Lx-epsilon_x;
				}
				if (cx<=0.0) {
					printf("moveSwarm in 3D - outside: cx=%lf<=0.0\n",cx);
					cx=epsilon_x;
				}
				if (cy>=Ly) {
					printf("moveSwarm in 3D - outside: cy=%lf>=%lf\n",cy,Ly);
					cy=Ly-epsilon_x;
				}
				if (cy<=0.0) {
					printf("moveSwarm in 3D - outside: cy=%lf<=0.0\n",cy);
					cy=epsilon_x;
				}
				if (cz>=0){
					printf("moveSwarm in 3D - outside: cz=%lf>=0.0\n",cz);
					cz=-epsilon_x;
				}
				if (cz<=-depth){
					printf("moveSwarm in 3D - outside: cz=%lf<=-%lf\n",cz,depth);
					cz=-depth+epsilon_x;
				}

				i = (int)(cx/dx_const);
				j = (int)(cy/dy_const);
				k = (int)((cz+depth)/dz_const);

				qq[k][j][i] += 1.0/geoq_fac[p];
				qq_cont[k][j][i] += 1.0;
			}
		}
		else {
			for (p=0; p<nlocal; p++) {
				PetscReal cx,cy,cz;
				PetscInt i,j,k;

				cx = array[3*p];
				cy = array[3*p+1];
				cz = array[3*p+2];

				if (cx>=Lx) {
					printf("moveSwarm in 3D - outside: cx=%lf>=%lf\n",cx,Lx);
					cx=Lx-epsilon_x;
				}
				if (cx<=0.0) {
					printf("moveSwarm in 3D - outside: cx=%lf<=0.0\n",cx);
					cx=epsilon_x;
				}
				if (cy>=Ly) {
					printf("moveSwarm in 3D - outside: cy=%lf>=%lf\n",cy,Ly);
					cy=Ly-epsilon_x;
				}
				if (cy<=0.0) {
					printf("moveSwarm in 3D - outside: cy=%lf<=0.0\n",cy);
					cy=epsilon_x;
				}
				if (cz>=0){
					printf("moveSwarm in 3D - outside: cz=%lf>=0.0\n",cz);
					cz=-epsilon_x;
				}
				if (cz<=-depth){
					printf("moveSwarm in 3D - outside: cz=%lf<=-%lf\n",cz,depth);
					cz=-depth+epsilon_x;
				}

				i = (int)(cx/dx_const);
				j = (int)(cy/dy_const);
				k = (int)((cz+depth)/dz_const);

				qq[k][j][i] += geoq_fac[p];
				qq_cont[k][j][i] += 1.0;
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
	//ierr = DMSwarmRestoreField(dms,"rho_fac",NULL,NULL,(void**)&rho_fac);CHKERRQ(ierr);
	//ierr = DMSwarmRestoreField(dms,"H_fac",NULL,NULL,(void**)&H_fac);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms,"layer",&bs,NULL,(void**)&layer_array);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms,"strain_fac",NULL,NULL,(void**)&strain_fac);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms,"strain_rate_fac",NULL,NULL,(void**)&strain_rate_fac);CHKERRQ(ierr);

	ierr = DMGlobalToLocalBegin(da_Thermal,Temper,INSERT_VALUES,local_Temper);
	ierr = DMGlobalToLocalEnd(  da_Thermal,Temper,INSERT_VALUES,local_Temper);

	ierr = DMDAVecGetArray(da_Thermal,local_Temper,&TT);CHKERRQ(ierr);

	ierr = DMGlobalToLocalBegin(da_Thermal,geoq_rho,INSERT_VALUES,local_geoq_rho);
	ierr = DMGlobalToLocalEnd(  da_Thermal,geoq_rho,INSERT_VALUES,local_geoq_rho);

	ierr = DMDAVecGetArray(da_Thermal,local_geoq_rho,&qq_rho);CHKERRQ(ierr);


	PetscInt       sx,sy,sz,mmx,mmy,mmz;

	ierr = DMDAGetCorners(da_Thermal,&sx,&sy,&sz,&mmx,&mmy,&mmz);CHKERRQ(ierr);

	int k,j,i;

	for (k=sz; k<sz+mmz; k++) {
		for (j=sy; j<sy+mmy; j++) {
			for (i=sx; i<sx+mmx; i++) {
				if (qq_rho[k][j][i]<air_threshold_density){
					TT[k][j][i]=0.0;
				}
			}
		}
	}

	ierr = DMDAVecRestoreArray(da_Thermal,local_geoq_rho,&qq_rho);CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(da_Thermal,local_Temper,&TT);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(da_Thermal,local_Temper,INSERT_VALUES,Temper);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(da_Thermal,local_Temper,INSERT_VALUES,Temper);CHKERRQ(ierr);

	if (sticky_blanket_air==1){
		if (visc_const_per_element==1){

			ierr = DMGlobalToLocalBegin(da_Thermal,geoq,INSERT_VALUES,local_geoq);
			ierr = DMGlobalToLocalEnd(  da_Thermal,geoq,INSERT_VALUES,local_geoq);
			ierr = DMDAVecGetArray(da_Thermal,local_geoq,&qq);CHKERRQ(ierr);

			ierr = DMGlobalToLocalBegin(da_Thermal,geoq_rho,INSERT_VALUES,local_geoq_rho);
			ierr = DMGlobalToLocalEnd(  da_Thermal,geoq_rho,INSERT_VALUES,local_geoq_rho);
			ierr = DMDAVecGetArray(da_Thermal,local_geoq_rho,&qq_rho);CHKERRQ(ierr);
			PetscReal rho_min=5.0/rho0_scaled,rho_max=2700.0*0.8/rho0_scaled, rho_mean; //check: in the future, make this free to the user modify
			for (k=sz; k<sz+mmz-1; k++) {
				for (j=sy; j<sy+mmy-1; j++) {
					for (i=sx; i<sx+mmx-1; i++) {
						rho_mean = (qq_rho[k][j][i] + qq_rho[k][j][i+1] + qq_rho[k+1][j][i] + qq_rho[k+1][j][i+1])/8.0;
						rho_mean+= (qq_rho[k][j+1][i] + qq_rho[k][j+1][i+1] + qq_rho[k+1][j+1][i] + qq_rho[k+1][j+1][i+1])/8.0;
						if (rho_mean>rho_min && rho_mean<rho_max){
							if (qq[k][j][i]<visc_MIN*50) qq[k][j][i]=visc_MIN*50;
							if (qq[k][j][i+1]<visc_MIN*50) qq[k][j][i+1]=visc_MIN*50;
							if (qq[k][j+1][i]<visc_MIN*50) qq[k][j+1][i]=visc_MIN*50;
							if (qq[k][j+1][i+1]<visc_MIN*50) qq[k][j+1][i+1]=visc_MIN*50;
						}
					}
				}
			}

			ierr = DMDAVecRestoreArray(da_Thermal,local_geoq,&qq);CHKERRQ(ierr);
			ierr = DMLocalToGlobalBegin(da_Thermal,local_geoq,INSERT_VALUES,geoq);CHKERRQ(ierr);
			ierr = DMLocalToGlobalEnd(da_Thermal,local_geoq,INSERT_VALUES,geoq);CHKERRQ(ierr);


			ierr = DMDAVecRestoreArray(da_Thermal,local_geoq_rho,&qq_rho);CHKERRQ(ierr);

		}
	}


	PetscFunctionReturn(0);

}
