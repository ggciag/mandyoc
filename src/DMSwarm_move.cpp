#include <petscsf.h>
#include <petscksp.h>
#include <petscdmda.h>
#include <petscdmswarm.h>
#include <petsc/private/dmimpl.h>
#include <petscmath.h>
#include "petscsys.h"

typedef struct {
	PetscScalar u;
	PetscScalar w;
	//PetscScalar p; //check for 3D future version
} Stokes;

double calc_visco_ponto(double T, double P, double x, double z,double geoq_ponto,double e2_inva,double strain_cumulate,
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

extern Vec local_P_aux;
extern Vec Pressure_aux;

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
extern PetscInt *p_add_i;
extern PetscInt *p_add_layer;
extern PetscReal *p_add_r_strain;
extern PetscReal *p_add_r_strain_rate;

extern unsigned int seed;

extern PetscScalar *inter_geoq;
extern PetscScalar *inter_A;
extern PetscScalar *inter_n;
extern PetscScalar *inter_Q;
extern PetscScalar *inter_V;

extern double e2_aux_MAX;
extern double e2_aux_MIN;

extern PetscInt RK4;

extern PetscInt periodic_boundary;


PetscReal linear_interpolation(PetscReal rx, PetscReal rz,PetscScalar V0, PetscScalar V1, PetscScalar V2, PetscScalar V3){
	PetscReal rfac,vx;
	rfac = (1.0-rx)*(1.0-rz);
	vx = V0 * rfac;
	
	rfac = (rx)*(1.0-rz);
	vx += V1 * rfac;
	
	rfac = (1.0-rx)*(rz);
	vx += V2 * rfac;
	
	rfac = (rx)*(rz);
	vx += V3 * rfac;

	return (vx);
}

PetscInt get_i(PetscReal cx){
	PetscInt i = (int)(cx/dx_const);
	if (i<0 || i>=Nx-1) {printf("weird i=%d get_i cx = %lf\n",i,cx); exit(1);}
	if (i==Nx-1) i=Nx-2;

	return i;
}

PetscInt get_k(PetscReal cz){
	PetscInt k = (int)((cz+depth)/dz_const);
	if (k<0 || k>=Nz-1) {printf("weird k=%d get_k cz=%lf\n",k,cz); exit(1);}
	if (k==Nz-1) k=Nz-2;

	return k;
}

PetscReal get_rx(PetscReal cx, PetscInt i){
	PetscReal rx = (cx-i*dx_const)/dx_const;
	if (rx<0 || rx>1) {printf("weird rx=%f\n",rx); exit(1);}

	return rx;
}

PetscReal get_rz(PetscReal cz, PetscInt k){
	PetscReal rz = (cz-(-depth+k*dz_const))/dz_const;
	if (rz<0 || rz>1) {printf("weird rz=%f\n",rz); exit(1);}

	return rz;
}

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

	//Pressure
	PetscScalar					**pp_aux;
	
	ierr = DMGlobalToLocalBegin(da_Thermal,Pressure_aux,INSERT_VALUES,local_P_aux);
	ierr = DMGlobalToLocalEnd(  da_Thermal,Pressure_aux,INSERT_VALUES,local_P_aux);
	ierr = DMDAVecGetArray(da_Thermal,local_P_aux,&pp_aux);CHKERRQ(ierr);
	
	
	PetscInt nlocal,bs,p;
	
	PetscReal *array;
	
	ierr = DMSwarmGetLocalSize(dms,&nlocal);CHKERRQ(ierr);
	
	ierr = DMSwarmGetField(dms,DMSwarmPICField_coor,&bs,NULL,(void**)&array);CHKERRQ(ierr);
	
	PetscReal N_x0[V_NE], N_z0[V_NE],strain[6],E2_invariant;
	
	
	PetscReal *strain_fac;
	PetscReal *strain_rate_fac;
	PetscReal *rarray;
	PetscInt *layer_array;
	ierr = DMSwarmGetField(dms,"strain_fac",&bs,NULL,(void**)&strain_fac);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"strain_rate_fac",&bs,NULL,(void**)&strain_rate_fac);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"geoq_fac",&bs,NULL,(void**)&rarray);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"layer",&bs,NULL,(void**)&layer_array);CHKERRQ(ierr);

	e2_aux_MAX = 0.0;
	e2_aux_MIN = 1.0E50;
	
	for (p=0; p<nlocal; p++) {
		PetscReal cx,cz,vx,vz,tp,Pp;

		PetscReal vxA,vxB,vxC,vxD;
		PetscReal vzA,vzB,vzC,vzD;

		PetscReal xA,xB,xC,xD;
		PetscReal zA,zB,zC,zD;

		PetscReal rx,rz;
		PetscInt i,k;
		PetscInt ii,kk;
		
		PetscReal kx,kz,ex,ez;
		PetscReal epsilon_x = 1.0E-7;

		if (RK4==1){
			//fourth-order Runge-Kutta scheme
			cx = array[2*p];
			cz = array[2*p+1];
			i = get_i(cx);
			k = get_k(cz);
			rx = get_rx(cx,i);
			rz = get_rz(cz,k);

			tp = linear_interpolation(rx,rz,tt[k][i],tt[k][i+1],tt[k+1][i],tt[k+1][i+1]);
			Pp = linear_interpolation(rx,rz,pp_aux[k][i],pp_aux[k][i+1],pp_aux[k+1][i],pp_aux[k+1][i+1]);

			xA = cx;
			zA = cz;
			vxA = linear_interpolation(rx,rz,VV[k][i].u,VV[k][i+1].u,VV[k+1][i].u,VV[k+1][i+1].u);
			vzA = linear_interpolation(rx,rz,VV[k][i].w,VV[k][i+1].w,VV[k+1][i].w,VV[k+1][i+1].w);

			xB = xA + vxA*dt/2.0;
			zB = zA + vzA*dt/2.0;
			if (xB>Lx || xB<0 || zB>0 || zB<-depth){
				array[2*p  ] = xB;
				array[2*p+1] = zB;
				break;
			}
			i = get_i(xB);
			k = get_k(zB);
			rx = get_rx(xB,i);
			rz = get_rz(zB,k);
			vxB = linear_interpolation(rx,rz,VV[k][i].u,VV[k][i+1].u,VV[k+1][i].u,VV[k+1][i+1].u);
			vzB = linear_interpolation(rx,rz,VV[k][i].w,VV[k][i+1].w,VV[k+1][i].w,VV[k+1][i+1].w);

			xC = xA + vxB*dt/2.0;
			zC = zA + vzB*dt/2.0;
			if (xC>Lx || xC<0 || zC>0 || zC<-depth){
				array[2*p  ] = xC;
				array[2*p+1] = zC;
				break;
			}
			i = get_i(xC);
			k = get_k(zC);
			rx = get_rx(xC,i);
			rz = get_rz(zC,k);
			vxC = linear_interpolation(rx,rz,VV[k][i].u,VV[k][i+1].u,VV[k+1][i].u,VV[k+1][i+1].u);
			vzC = linear_interpolation(rx,rz,VV[k][i].w,VV[k][i+1].w,VV[k+1][i].w,VV[k+1][i+1].w);

			xD = xA + vxC*dt;
			zD = zA + vzC*dt;
			if (xD>Lx || xD<0 || zD>0 || zD<-depth){
				array[2*p  ] = xD;
				array[2*p+1] = zD;
				break;
			}
			i = get_i(xD);
			k = get_k(zD);
			rx = get_rx(xD,i);
			rz = get_rz(zD,k);
			vxD = linear_interpolation(rx,rz,VV[k][i].u,VV[k][i+1].u,VV[k+1][i].u,VV[k+1][i+1].u);
			vzD = linear_interpolation(rx,rz,VV[k][i].w,VV[k][i+1].w,VV[k+1][i].w,VV[k+1][i+1].w);

			vx = (vxA + 2*vxB + 2*vxC + vxD)/6.0;
			vz = (vzA + 2*vzB + 2*vzC + vzD)/6.0;

			array[2*p  ] += dt * vx;
			array[2*p+1] += dt * vz;
		}
		else {
			cx = array[2*p];
			if (cx>=Lx) {
				printf("moveSwarm - outside: cx=%lf>=%lf\n",cx,Lx);
				cx=Lx-epsilon_x;
			}
			if (cx<=0.0) {
				printf("moveSwarm - outside: cx=%lf<=0.0\n",cx);
				cx=epsilon_x;
			}

			cz = array[2*p+1];
			if (cz>=0){
				printf("moveSwarm - outside: cz=%lf>=0.0\n",cz);
				cz=-epsilon_x;
			}
			if (cz<=-depth){
				printf("moveSwarm - outside: cz=%lf<=-%lf\n",cz,depth);
				cz=-depth+epsilon_x;
			}

			i = get_i(cx);
			k = get_k(cz);
			rx = get_rx(cx,i);
			rz = get_rz(cz,k);

			tp = linear_interpolation(rx,rz,tt[k][i],tt[k][i+1],tt[k+1][i],tt[k+1][i+1]);
			Pp = linear_interpolation(rx,rz,pp_aux[k][i],pp_aux[k][i+1],pp_aux[k+1][i],pp_aux[k+1][i+1]);

			vx = linear_interpolation(rx,rz,VV[k][i].u,VV[k][i+1].u,VV[k+1][i].u,VV[k+1][i+1].u);
			vz = linear_interpolation(rx,rz,VV[k][i].w,VV[k][i+1].w,VV[k+1][i].w,VV[k+1][i+1].w);

			array[2*p  ] += dt * vx;
			array[2*p+1] += dt * vz;
		}

		if (periodic_boundary==1){
			if (array[2*p]>Lx) array[2*p]-=Lx;
			if (array[2*p]<0) array[2*p]+=Lx;
		}

		
		
		kx = 2*rx-1;
		kz = 2*rz-1;
		
		
		PetscInt cont = 0;
		for (ez=-1.;ez<=1.;ez+=2.){
			for (ex=-1.;ex<=1.;ex+=2.){
				//N[cont]=(1+ex*kx)*(1+ey*ky)*(1+ez*kz)/4.0;
				N_x0[cont]=ex*(1+ez*kz)/2.0/dx_const;//check: only for 2D model
				N_z0[cont]=(1+ex*kx)*ez/2.0/dz_const;//check: only for 2D model
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
				
				//strain[3] += 0.0;//check: only for 2D model
				//strain[4] += 0.0;//check: only for 2D model
				strain[5] += N_x0[cont]*VV[kk][ii].w + N_z0[cont]*VV[kk][ii].u;
				
				
				cont++;
			}
		}
		
		

		PetscReal strain_mean = (strain[0] + strain[1] + strain[2])/3.0;
		strain[0]-=strain_mean;
		strain[1]-=strain_mean;
		strain[2]-=strain_mean;
		E2_invariant = 0;
		for (ii=0;ii<3;ii++) E2_invariant+=strain[ii]*strain[ii];
		for (ii=3;ii<6;ii++) E2_invariant+=strain[ii]*strain[ii]/2.0;
		
		E2_invariant = PetscSqrtReal(E2_invariant/2.0);

		if (E2_invariant<e2_aux_MIN) e2_aux_MIN=E2_invariant;
		if (E2_invariant>e2_aux_MAX) e2_aux_MAX=E2_invariant;
		
		
		strain_fac[p]+= dt*E2_invariant;//cumulative strain
		strain_rate_fac[p] = E2_invariant;
		
	
		rarray[p] = calc_visco_ponto(tp,Pp,cx,cz,inter_geoq[layer_array[p]],E2_invariant,strain_fac[p],
									 inter_A[layer_array[p]], inter_n[layer_array[p]], inter_Q[layer_array[p]], inter_V[layer_array[p]]);
		
		
		
		
	}


	ierr = DMSwarmRestoreField(dms,"geoq_fac",&bs,NULL,(void**)&rarray);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms,"strain_fac",&bs,NULL,(void**)&strain_fac);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms,"strain_rate_fac",&bs,NULL,(void**)&strain_rate_fac);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms,"layer",&bs,NULL,(void**)&layer_array);CHKERRQ(ierr);
	
	 
	ierr = DMSwarmRestoreField(dms,DMSwarmPICField_coor,&bs,NULL,(void**)&array);CHKERRQ(ierr);
	
	ierr = DMSwarmMigrate(dms,PETSC_TRUE);CHKERRQ(ierr);
	
	
	ierr = DMDAVecRestoreArray(da_Veloc,local_V,&VV);CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(da_Thermal,local_Temper,&tt);CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(da_Thermal,local_P_aux,&pp_aux);CHKERRQ(ierr);
	
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
	PetscReal *strain_fac;
	PetscReal *strain_rate_fac;
	
	ierr = DMSwarmGetLocalSize(dms,&nlocal);CHKERRQ(ierr);
	
	ierr = DMSwarmGetField(dms,DMSwarmPICField_coor,&bs,NULL,(void**)&array);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"cont",&bs,NULL,(void**)&carray);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"itag",&bs,NULL,(void**)&iarray);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"layer",&bs,NULL,(void**)&layer_array);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"geoq_fac",&bs,NULL,(void**)&rarray);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"strain_fac",&bs,NULL,(void**)&strain_fac);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"strain_rate_fac",&bs,NULL,(void**)&strain_rate_fac);CHKERRQ(ierr);
	
	PetscInt Mx=0,mx=10000,Mz=0,mz=10000;
	PetscInt       sx,sz,mmx,mmz;
	
	ierr = DMDAGetCorners(da_Thermal,&sx,&sz,NULL,&mmx,&mmz,NULL);CHKERRQ(ierr);

	PetscReal cx,cz,dx,dz;
	PetscInt i,k;
	PetscReal cx_v[10],cz_v[10];
	PetscReal epsilon_x = 1.0E-7;
	
	for (p=0; p<nlocal; p++) {

		
		cx = array[2*p];
		if (cx>=Lx) {
			printf("Swarm_add_remove - outside: cx=%lf>=%lf\n",cx,Lx);
			cx=Lx-epsilon_x;
		}
		if (cx<=0.0) {
			printf("Swarm_add_remove - outside: cx=%lf<=0.0\n",cx);
			cx=epsilon_x;
		}

		cz = array[2*p+1];
		if (cz>=0){
			printf("Swarm_add_remove - outside: cz=%lf>=0.0\n",cz);
			cz=-epsilon_x;
		}
		if (cz<=-depth){
			printf("Swarm_add_remove - outside: cz=%lf<=-%lf\n",cz,depth);
			cz=-depth+epsilon_x;
		}

		i = (int)(cx/dx_const);
		k = (int)((cz+depth)/dz_const);
		
		if (i<0 || i>=Nx-1) {printf("weird i=%d  add remove cx = %lf\n",i,cx); exit(1);}
		if (k<0 || k>=Nz-1) {printf("weird k=%d  add remove cz = %lf\n",k,cz); exit(1);}
		
		if (i==Nx-1) i=Nx-2;
		if (k==Nz-1) k=Nz-2;
		
		qq_cont[k][i] += 1.0;
		
		carray[p] = k*Nx + i;
		
		if (Mx<i) Mx=i;
		if (mx>i) mx=i;
		
		if (Mz<k) Mz=k;
		if (mz>k) mz=k;
	
	}

	//printf("Swarm move: %d %d %d %d\n",mx,Mx,mz,Mz);
	
	PetscInt max_particles_per_ele=particles_per_ele+particles_per_ele/10+2;
	PetscInt min_particles_per_ele=particles_per_ele-particles_per_ele/10-2;
	
	PetscInt kk,pp;
	PetscInt cont_p;
	
	PetscInt cont_p_remove=0;
	
	
	
	PetscInt cont_p_add=0;
	
	
	
	
	
	PetscReal dist,dist_p;
	PetscInt chosen=0;
	
	PetscReal rx,rz,xx,zz;
	
	
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
					rx = 2.0*(float)rand_r(&seed)/RAND_MAX-1.0;
					rz = 2.0*(float)rand_r(&seed)/RAND_MAX-1.0;
					
					cx_v[pp] = i*dx_const + (0.5*rx+0.5)*dx_const;
					cz_v[pp] = k*dz_const - depth + (0.5*rz+0.5)*dz_const;
					
				}
				
				dist = 0;
				int p_prox,p_prox_total=0;
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
				p_add_r_strain[cont_p_add] = strain_fac[p_prox_total];
				p_add_r_strain_rate[cont_p_add] = strain_rate_fac[p_prox_total];
				
				
				cont_p_add++;
				if (cont_p_add>particles_add_remove){
					printf("MUITO2\n");
					exit(1);
				}
			}
			
		}
	
	}

	//PetscPrintf(PETSC_COMM_WORLD,"Swarm move: 1\n");
	
	
	ierr = DMSwarmRestoreField(dms,"cont",&bs,NULL,(void**)&carray);CHKERRQ(ierr);
	
	ierr = DMSwarmRestoreField(dms,DMSwarmPICField_coor,&bs,NULL,(void**)&array);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms,"geoq_fac",&bs,NULL,(void**)&rarray);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms,"itag",&bs,NULL,(void**)&iarray);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms,"layer",&bs,NULL,(void**)&layer_array);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms,"strain_fac",&bs,NULL,(void**)&strain_fac);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms,"strain_rate_fac",&bs,NULL,(void**)&strain_rate_fac);CHKERRQ(ierr);
	
	ierr = DMDAVecRestoreArray(da_Thermal,local_geoq_cont,&qq_cont);CHKERRQ(ierr);
	
	ierr = DMSwarmGetLocalSize(dms,&nlocal);CHKERRQ(ierr);
	
	
	
	for (pp=0; pp<cont_p_remove; pp++){
		p_i[pp]=pp;
	}
	
	if (cont_p_remove>0)
		PetscSortIntWithPermutation(cont_p_remove,p_remove,p_i);
	
	for (pp=cont_p_remove-1; pp>=0; pp--){
		DMSwarmRemovePointAtIndex(dms,p_remove[p_i[pp]]);
	}
	
	ierr = DMSwarmGetLocalSize(dms,&nlocal);CHKERRQ(ierr);
	
	if (cont_p_add>0){
		ierr = DMSwarmAddNPoints(dms,cont_p_add);
		
		ierr = DMSwarmGetField(dms,DMSwarmPICField_coor,&bs,NULL,(void**)&array);CHKERRQ(ierr);
		ierr = DMSwarmGetField(dms,"itag",&bs,NULL,(void**)&iarray);CHKERRQ(ierr);
		ierr = DMSwarmGetField(dms,"layer",&bs,NULL,(void**)&layer_array);CHKERRQ(ierr);
		ierr = DMSwarmGetField(dms,"geoq_fac",&bs,NULL,(void**)&rarray);CHKERRQ(ierr);
		ierr = DMSwarmGetField(dms,"strain_fac",&bs,NULL,(void**)&strain_fac);CHKERRQ(ierr);
		ierr = DMSwarmGetField(dms,"strain_rate_fac",&bs,NULL,(void**)&strain_rate_fac);CHKERRQ(ierr);
		
		for (pp=0; pp<cont_p_add; pp++){
			array[(nlocal+pp)*2] = p_add_coor[pp*2];
			array[(nlocal+pp)*2+1] = p_add_coor[pp*2+1];
			
			rarray[nlocal+pp] = p_add_r[pp];
			
			strain_fac[nlocal+pp] = p_add_r_strain[pp];
			strain_rate_fac[nlocal+pp] = p_add_r_strain_rate[pp];
			
			iarray[nlocal+pp] = p_add_i[pp];
			layer_array[nlocal+pp] = p_add_layer[pp];
		}
		
		ierr = DMSwarmRestoreField(dms,DMSwarmPICField_coor,&bs,NULL,(void**)&array);CHKERRQ(ierr);
		ierr = DMSwarmRestoreField(dms,"geoq_fac",&bs,NULL,(void**)&rarray);CHKERRQ(ierr);
		ierr = DMSwarmRestoreField(dms,"itag",&bs,NULL,(void**)&iarray);CHKERRQ(ierr);
		ierr = DMSwarmRestoreField(dms,"layer",&bs,NULL,(void**)&layer_array);CHKERRQ(ierr);
		ierr = DMSwarmRestoreField(dms,"strain_fac",&bs,NULL,(void**)&strain_fac);CHKERRQ(ierr);
		ierr = DMSwarmRestoreField(dms,"strain_rate_fac",&bs,NULL,(void**)&strain_rate_fac);CHKERRQ(ierr);
		
	}
	
	ierr = DMSwarmGetLocalSize(dms,&nlocal);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"Swarm advection: done\n\n");
	
	PetscFunctionReturn(0);
}

