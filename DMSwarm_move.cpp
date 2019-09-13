//#include <petscksp.h>
//#include <petscksp.h>
//#include <petscmath.h>
///#include <petscdmda.h>
///#include <petscdmswarm.h>
#include <petscsf.h>
//#include <petscdm.h>
#include <petscksp.h>
#include <petscdmda.h>
//#include <petscdmshell.h>
#include <petscdmswarm.h>
#include <petsc/private/dmimpl.h>
#include <petscmath.h>
#include "petscsys.h"

typedef struct {
	PetscScalar u;
	PetscScalar w;
	//PetscScalar p;
} Stokes;

//PetscErrorCode SwarmViewGP(DM dms,const char prefix[]);

double calc_visco_ponto(double T,double x, double z,double geoq_ponto,double e2_inva,double strain_cumulate,
						double A, double n_exp, double QE, double VE);

extern DM dms;

extern DM da_Veloc;

extern DM da_Thermal;

extern long Nx,Nz;

extern double dx_const;
extern double dz_const;

extern double Lx, depth;

extern Vec local_V,Veloc_weight;

extern Vec local_Temper,Temper;

extern Vec geoq_cont,local_geoq_cont;

extern PetscInt particles_per_ele;
extern PetscInt cont_particles;

extern long V_NE;

extern PetscInt particles_add_remove;

extern PetscInt *ppp;
extern PetscInt *p_remove;
extern PetscInt *p_i;

extern PetscReal *p_add_coor;
extern PetscReal *p_add_r;
extern PetscReal *p_add_r_rho;
extern PetscReal *p_add_r_H;
extern PetscInt *p_add_i;
extern PetscInt *p_add_layer;
extern PetscReal *p_add_r_strain;

extern unsigned int seed;

extern PetscScalar *inter_geoq;
extern PetscScalar *inter_A;
extern PetscScalar *inter_n;
extern PetscScalar *inter_Q;
extern PetscScalar *inter_V;

PetscErrorCode moveSwarm(PetscReal dt)
{
	PetscErrorCode ierr=0;
	//Velocity
	Stokes					**VV;
	ierr = VecZeroEntries(local_V);CHKERRQ(ierr);
	
	ierr = DMGlobalToLocalBegin(da_Veloc,Veloc_weight,INSERT_VALUES,local_V);
	ierr = DMGlobalToLocalEnd(  da_Veloc,Veloc_weight,INSERT_VALUES,local_V);
	
	ierr = DMDAVecGetArray(da_Veloc,local_V,&VV);CHKERRQ(ierr);
	
	
	//Temperature
	PetscScalar             **tt;
	
	ierr = VecZeroEntries(local_Temper);CHKERRQ(ierr);
	
	ierr = DMGlobalToLocalBegin(da_Thermal,Temper,INSERT_VALUES,local_Temper);
	ierr = DMGlobalToLocalEnd(  da_Thermal,Temper,INSERT_VALUES,local_Temper);
	
	ierr = DMDAVecGetArray(da_Thermal,local_Temper,&tt);CHKERRQ(ierr);
	
	
	
	PetscInt nlocal,bs,p;
	
	PetscReal *array;
	
	ierr = DMSwarmGetLocalSize(dms,&nlocal);CHKERRQ(ierr);
	
	ierr = DMSwarmGetField(dms,DMSwarmPICField_coor,&bs,NULL,(void**)&array);CHKERRQ(ierr);
	
	PetscReal N_x0[V_NE], N_z0[V_NE],strain[6],E2_invariant;
	
	
	PetscReal *strain_fac;
	PetscReal *rarray;
	PetscInt *layer_array;
	ierr = DMSwarmGetField(dms,"strain_fac",&bs,NULL,(void**)&strain_fac);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"geoq_fac",&bs,NULL,(void**)&rarray);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"layer",&bs,NULL,(void**)&layer_array);CHKERRQ(ierr);
	
	for (p=0; p<nlocal; p++) {
		PetscReal cx,cz,vx,vz,tp;
		PetscReal rx,rz,rfac;
		PetscInt i,k;
		PetscInt ii,kk;
		
		PetscReal kx,kz,ex,ez;
		
		cx = array[2*p];
		cz = array[2*p+1];
		
		i = (int)(cx/dx_const);
		k = (int)((cz+depth)/dz_const);
		
		
		
		if (i<0 || i>=Nx-1) {printf("estranho i=%d\n",i); exit(1);}
		if (k<0 || k>=Nz-1) {printf("estranho k=%d\n",k); exit(1);}
		
		if (i==Nx-1) i=Nx-2;
		if (k==Nz-1) k=Nz-2;
		
		
		
		//VV[k][j][i].u + ;
		
		rx = (cx-i*dx_const)/dx_const;
		rz = (cz-(-depth+k*dz_const))/dz_const;
		
		if (rx<0 || rx>1) {printf("estranho rx=%f\n",rx); exit(1);}
		if (rz<0 || rz>1) {printf("estranho rz=%f\n",rz); exit(1);}
		
		rfac = (1.0-rx)*(1.0-rz);
		vx = VV[k][i].u * rfac;
		vz = VV[k][i].w * rfac;
		tp = tt[k][i] * rfac;
		
		rfac = (rx)*(1.0-rz);
		vx += VV[k][i+1].u * rfac;
		vz += VV[k][i+1].w * rfac;
		tp += tt[k][i+1] * rfac;
		
		rfac = (1.0-rx)*(rz);
		vx += VV[k+1][i].u * rfac;
		vz += VV[k+1][i].w * rfac;
		tp += tt[k+1][i] * rfac;
		
		rfac = (rx)*(rz);
		vx += VV[k+1][i+1].u * rfac;
		vz += VV[k+1][i+1].w * rfac;
		tp += tt[k+1][i+1] * rfac;
		
		///////// strain
		
		
		kx = 2*rx-1;
		kz = 2*rz-1;
		
		
		PetscInt cont = 0;
		for (ez=-1.;ez<=1.;ez+=2.){
			for (ex=-1.;ex<=1.;ex+=2.){
				//N[cont]=(1+ex*kx)*(1+ey*ky)*(1+ez*kz)/4.0;
				N_x0[cont]=ex*(1+ez*kz)/2.0/dx_const;//!!!2d checar
				N_z0[cont]=(1+ex*kx)*ez/2.0/dz_const;//!!!2d checar
				cont++;
			}
		}
		
		cont=0;
		
		for (ii=0;ii<6;ii++) strain[ii]=0.0;
		E2_invariant=0.0;
		
		for (kk=k;kk<=k+1;kk++){
			for (ii=i;ii<=i+1;ii++){
				strain[0] += N_x0[cont]*VV[kk][ii].u;
				strain[1] += 0;//N_y0[cont]*VV[kk][jj][ii].v;
				strain[2] += N_z0[cont]*VV[kk][ii].w;
				
				//strain[3] += 0.0;//!!!2d checar
				//strain[4] += 0.0;//!!!2d checar
				strain[5] += N_x0[cont]*VV[kk][ii].w + N_z0[cont]*VV[kk][ii].u;
				
				
				cont++;
			}
		}
		
		
		
		
		
		E2_invariant = (strain[0]-strain[1])*(strain[0]-strain[1]);
		E2_invariant+= (strain[1]-strain[2])*(strain[1]-strain[2]);
		E2_invariant+= (strain[2]-strain[0])*(strain[2]-strain[0]);
		E2_invariant/=6.0;
		E2_invariant+= strain[3]*strain[3];
		E2_invariant+= strain[4]*strain[4];
		E2_invariant+= strain[5]*strain[5];
		
		
		
		
		strain_fac[p]+= dt*PetscSqrtReal(E2_invariant);//original!!!!
		//strain_fac[p]= PetscSqrtReal(E2_invariant);//!!!! não é o cumulativo! apenas o instantaneo.
		
		
		rarray[p] = calc_visco_ponto(tp,cx,cz,inter_geoq[layer_array[p]],PetscSqrtReal(E2_invariant),strain_fac[p]/*!!!!checar*/,
									 inter_A[layer_array[p]], inter_n[layer_array[p]], inter_Q[layer_array[p]], inter_V[layer_array[p]]);
		
		
		//dt=0.01;
		//vx = 0;
		//vy = -cz;
		//vz = cy;
		
		array[2*p  ] += dt * vx;
		array[2*p+1] += dt * vz;
		
		
		
		
	}
	ierr = DMSwarmRestoreField(dms,"geoq_fac",&bs,NULL,(void**)&rarray);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms,"strain_fac",&bs,NULL,(void**)&strain_fac);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms,"layer",&bs,NULL,(void**)&layer_array);CHKERRQ(ierr);
	
	 
	ierr = DMSwarmRestoreField(dms,DMSwarmPICField_coor,&bs,NULL,(void**)&array);CHKERRQ(ierr);
	
	ierr = DMSwarmMigrate(dms,PETSC_TRUE);CHKERRQ(ierr);
	
	//PetscSNPrintf(prefix,PETSC_MAX_PATH_LEN-1,"step%d",tk);
	//ierr = SwarmViewGP(dms,prefix);CHKERRQ(ierr);
	
	//ierr = SwarmViewGP(dms,"step1");CHKERRQ(ierr);
	
	ierr = DMDAVecRestoreArray(da_Veloc,local_V,&VV);CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(da_Thermal,local_Temper,&tt);CHKERRQ(ierr);
	
	//exit(1);
	
	PetscFunctionReturn(0);
	
}



PetscErrorCode Swarm_add_remove()
{
	PetscErrorCode ierr=0;
	
	PetscInt *carray;
	
	
	PetscScalar             **qq_cont;
	
	ierr = VecSet(geoq_cont,0.0);CHKERRQ(ierr);
	
	ierr = DMGlobalToLocalBegin(da_Thermal,geoq_cont,INSERT_VALUES,local_geoq_cont);
	ierr = DMGlobalToLocalEnd(  da_Thermal,geoq_cont,INSERT_VALUES,local_geoq_cont);
	
	ierr = DMDAVecGetArray(da_Thermal,local_geoq_cont,&qq_cont);CHKERRQ(ierr);
	
	
	PetscInt nlocal,bs,p;
	
	PetscReal *array;
	PetscInt *iarray;
	PetscInt *layer_array;
	PetscReal *rarray;
	PetscReal *rarray_rho;
	PetscReal *rarray_H;
	PetscReal *strain_fac;
	
	ierr = DMSwarmGetLocalSize(dms,&nlocal);CHKERRQ(ierr);
	
	ierr = DMSwarmGetField(dms,DMSwarmPICField_coor,&bs,NULL,(void**)&array);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"cont",&bs,NULL,(void**)&carray);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"itag",&bs,NULL,(void**)&iarray);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"layer",&bs,NULL,(void**)&layer_array);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"geoq_fac",&bs,NULL,(void**)&rarray);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"rho_fac",&bs,NULL,(void**)&rarray_rho);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"H_fac",&bs,NULL,(void**)&rarray_H);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"strain_fac",&bs,NULL,(void**)&strain_fac);CHKERRQ(ierr);
	
	PetscInt Mx=0,mx=10000,Mz=0,mz=10000;
	PetscInt       sx,sz,mmx,mmz;
	
	ierr = DMDAGetCorners(da_Thermal,&sx,&sz,NULL,&mmx,&mmz,NULL);CHKERRQ(ierr);

	PetscReal cx,cz,dx,dz;
	PetscInt i,k;
	PetscReal cx_v[10],cz_v[10];
	
	for (p=0; p<nlocal; p++) {

		
		cx = array[2*p];
		cz = array[2*p+1];
		
		i = (int)(cx/dx_const);
		k = (int)((cz+depth)/dz_const);
		
		if (i<0 || i>=Nx-1) {printf("estranho i=%d\n",i); exit(1);}
		if (k<0 || k>=Nz-1) {printf("estranho k=%d\n",k); exit(1);}
		
		if (i==Nx-1) i=Nx-2;
		if (k==Nz-1) k=Nz-2;
		
		qq_cont[k][i] += 1.0;
		
		carray[p] = k*Nx + i;
		
		if (Mx<i) Mx=i;
		if (mx>i) mx=i;
		
		if (Mz<k) Mz=k;
		if (mz>k) mz=k;
	
	}
	
	PetscInt max_particles_per_ele=particles_per_ele+particles_per_ele/10+2;
	PetscInt min_particles_per_ele=particles_per_ele-particles_per_ele/10-2;
	
	PetscInt kk,pp;
	PetscInt cont_p;
	
	PetscInt cont_p_remove=0;
	
	
	
	PetscInt cont_p_add=0;
	
	
	
	
	
	PetscReal dist,dist_p;
	PetscInt chosen;
	
	//PetscRandom rand;
	
	PetscReal rx,rz,xx,zz;
	
	//ierr = PetscRandomCreate(PETSC_COMM_SELF,&rand);CHKERRQ(ierr);
	//ierr = PetscRandomSetType(rand,PETSCRAND48);CHKERRQ(ierr);
	//ierr = PetscRandomSetInterval(rand,-1.0,1.0);CHKERRQ(ierr);
	
	
	for (k=mz; k<=Mz; k++){
		for (i=mx; i<=Mx; i++){
			if (qq_cont[k][i]>max_particles_per_ele){
				qq_cont[k][i] -= 1.0;
				
				
				cont_p=0;
				kk = k*Nx + i;
				for (p=0; p<nlocal; p++){
					if (carray[p]==kk){
						ppp[cont_p]=p;
						cont_p++;
					}
				}
				
				dist = 1.0E40;
				for (p=0; p<cont_p; p++){
					xx = array[ppp[p]*2];
					zz = array[ppp[p]*2+1];
					dist_p=1.0E40;
					
					for (pp=0; pp<cont_p; pp++){
						dx = xx - array[ppp[pp]*2];
						dz = zz - array[ppp[pp]*2+1];
						if (dist_p>dx*dx+dz*dz && p!=pp)
							dist_p=dx*dx+dz*dz;
					}
					if (dist_p<dist){
						dist = dist_p;
						chosen = ppp[p];
					}
				}
				
				p_remove[cont_p_remove]=chosen;
				cont_p_remove++;
				
				
				if (cont_p_remove>particles_add_remove){
					printf("MUITO1\n");
					exit(1);
				}
				
				
				//printf("REMOVEU %d %d %d: %ld!\n",k,j,i,(long)qq_cont[k][j][i]);
			}
			
			if (qq_cont[k][i]<min_particles_per_ele){
				qq_cont[k][i] += 1.0;
				cont_p=0;
				
				kk = k*Nx + i;
				for (p=0; p<nlocal; p++){
					if (carray[p]==kk){
						ppp[cont_p]=p;
						cont_p++;
					}
				}
				for (pp=0;pp<10;pp++){
					//ierr = PetscRandomGetValueReal(rand,&rx);CHKERRQ(ierr);
					//ierr = PetscRandomGetValueReal(rand,&ry);CHKERRQ(ierr);
					//ierr = PetscRandomGetValueReal(rand,&rz);CHKERRQ(ierr);
					rx = 2.0*(float)rand_r(&seed)/RAND_MAX-1.0;
					rz = 2.0*(float)rand_r(&seed)/RAND_MAX-1.0;
					
					cx_v[pp] = i*dx_const + (0.5*rx+0.5)*dx_const;
					cz_v[pp] = k*dz_const - depth + (0.5*rz+0.5)*dz_const;
					
				}
				
				dist = 0;
				int p_prox,p_prox_total;
				for (pp=0;pp<10;pp++){
					cx = cx_v[pp];
					cz = cz_v[pp];
					dist_p = 1.0E30;
					for (p=0;p<cont_p;p++){
						dx = cx - array[ppp[p]*2];
						dz = cz - array[ppp[p]*2+1];
						
						if (dx*dx+dz*dz<dist_p){
							p_prox = ppp[p];
							dist_p = dx*dx+dz*dz;
						}
					}
					if (dist<dist_p){
						p_prox_total = p_prox;
						dist=dist_p;
						chosen = pp;
					}
				}
				p_add_coor[cont_p_add*2] = cx_v[chosen];
				p_add_coor[cont_p_add*2+1] = cz_v[chosen];
				
				p_add_i[cont_p_add] = cont_particles%particles_per_ele;
				cont_particles++;
				
				p_add_layer[cont_p_add] = layer_array[p_prox_total];
				
				p_add_r[cont_p_add] = rarray[p_prox_total];
				p_add_r_rho[cont_p_add] = rarray_rho[p_prox_total];
				p_add_r_H[cont_p_add] = rarray_H[p_prox_total];
				p_add_r_strain[cont_p_add] = strain_fac[p_prox_total];
				
				//printf("ADDED %d %d %d: !\n",k,j,i);
				//printf("ADDED %lf %lf %lf: !\n",cx_v[chosen],cy_v[chosen],cz_v[chosen]);
				
				
				cont_p_add++;
				if (cont_p_add>particles_add_remove){
					printf("MUITO2\n");
					exit(1);
				}
			}
			
		}
	
	}
	
	
	//ierr = PetscRandomDestroy(&rand);CHKERRQ(ierr);
	
	//printf("%d %d   %d %d   %d %d\n",mx,Mx,my,My,mz,Mz);
	
	//printf("b: %d %d   %d %d   %d %d\n",sx,sx+mmx-1,sy,sy+mmy-1,sz,sz+mmz-1);
	
	
	
	ierr = DMSwarmRestoreField(dms,"cont",&bs,NULL,(void**)&carray);CHKERRQ(ierr);
	
	ierr = DMSwarmRestoreField(dms,DMSwarmPICField_coor,&bs,NULL,(void**)&array);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms,"geoq_fac",&bs,NULL,(void**)&rarray);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms,"rho_fac",&bs,NULL,(void**)&rarray_rho);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms,"H_fac",&bs,NULL,(void**)&rarray_H);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms,"itag",&bs,NULL,(void**)&iarray);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms,"layer",&bs,NULL,(void**)&layer_array);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms,"strain_fac",&bs,NULL,(void**)&strain_fac);CHKERRQ(ierr);
	
	ierr = DMDAVecRestoreArray(da_Thermal,local_geoq_cont,&qq_cont);CHKERRQ(ierr);
	
	ierr = DMSwarmGetLocalSize(dms,&nlocal);CHKERRQ(ierr);
	printf("nlocal_%d_antes\n",nlocal);
	
	
	for (pp=0; pp<cont_p_remove; pp++){
		p_i[pp]=pp;
	}
	
	if (cont_p_remove>0)
		PetscSortIntWithPermutation(cont_p_remove,p_remove,p_i);
	
	//for (pp=0; pp<cont_p_remove; pp++){
	for (pp=cont_p_remove-1; pp>=0; pp--){
		DMSwarmRemovePointAtIndex(dms,p_remove[p_i[pp]]);
	}
	
	ierr = DMSwarmGetLocalSize(dms,&nlocal);CHKERRQ(ierr);
	printf("nlocal_%d %d %d_depois\n",nlocal,cont_p_remove,particles_add_remove);
	
	if (cont_p_add>0){
		ierr = DMSwarmAddNPoints(dms,cont_p_add);
		
		ierr = DMSwarmGetField(dms,DMSwarmPICField_coor,&bs,NULL,(void**)&array);CHKERRQ(ierr);
		ierr = DMSwarmGetField(dms,"itag",&bs,NULL,(void**)&iarray);CHKERRQ(ierr);
		ierr = DMSwarmGetField(dms,"layer",&bs,NULL,(void**)&layer_array);CHKERRQ(ierr);
		ierr = DMSwarmGetField(dms,"geoq_fac",&bs,NULL,(void**)&rarray);CHKERRQ(ierr);
		ierr = DMSwarmGetField(dms,"rho_fac",&bs,NULL,(void**)&rarray_rho);CHKERRQ(ierr);
		ierr = DMSwarmGetField(dms,"H_fac",&bs,NULL,(void**)&rarray_H);CHKERRQ(ierr);
		ierr = DMSwarmGetField(dms,"strain_fac",&bs,NULL,(void**)&strain_fac);CHKERRQ(ierr);
		
		for (pp=0; pp<cont_p_add; pp++){
			array[(nlocal+pp)*2] = p_add_coor[pp*2];
			array[(nlocal+pp)*2+1] = p_add_coor[pp*2+1];
			
			rarray[nlocal+pp] = p_add_r[pp];
			
			rarray_rho[nlocal+pp] = p_add_r_rho[pp];
			
			rarray_H[nlocal+pp] = p_add_r_H[pp];
			
			strain_fac[nlocal+pp] = p_add_r_strain[pp];
			
			iarray[nlocal+pp] = p_add_i[pp];
			layer_array[nlocal+pp] = p_add_layer[pp];
		}
		
		ierr = DMSwarmRestoreField(dms,DMSwarmPICField_coor,&bs,NULL,(void**)&array);CHKERRQ(ierr);
		ierr = DMSwarmRestoreField(dms,"geoq_fac",&bs,NULL,(void**)&rarray);CHKERRQ(ierr);
		ierr = DMSwarmRestoreField(dms,"rho_fac",&bs,NULL,(void**)&rarray_rho);CHKERRQ(ierr);
		ierr = DMSwarmRestoreField(dms,"H_fac",&bs,NULL,(void**)&rarray_H);CHKERRQ(ierr);
		ierr = DMSwarmRestoreField(dms,"itag",&bs,NULL,(void**)&iarray);CHKERRQ(ierr);
		ierr = DMSwarmRestoreField(dms,"layer",&bs,NULL,(void**)&layer_array);CHKERRQ(ierr);
		ierr = DMSwarmRestoreField(dms,"strain_fac",&bs,NULL,(void**)&strain_fac);CHKERRQ(ierr);
		
	}
	
	ierr = DMSwarmGetLocalSize(dms,&nlocal);CHKERRQ(ierr);
	printf("nlocal_%d %d %d_depois2\n",nlocal,cont_p_add,particles_add_remove);
	
	
	//exit(1);
	PetscFunctionReturn(0);
}

