#include <petscsf.h>
#include <petscdm.h>
#include <petscksp.h>
#include <petscdmda.h>
#include <petscdmshell.h>
#include <petscdmswarm.h>
#include <petsc/private/dmimpl.h>
#include <petscmath.h>

extern DM dmcell;

extern DM dms;

extern DM da_Veloc;

extern long Nx,Ny;

extern double dx_const;
extern double dy_const;
extern double dz_const;

extern double Lx, Ly, depth;

extern double H_lito;
extern double escala_viscosidade;

extern PetscInt particles_per_ele;
extern PetscInt nx_ppe;
extern PetscInt ny_ppe;
extern PetscInt nz_ppe;

extern double H_per_mass;

extern double h_air;

extern double RHOM;

extern int n_interfaces;
extern PetscScalar *interfaces;

extern PetscScalar *inter_rho;
extern PetscScalar *inter_geoq;
extern PetscScalar *inter_H;

extern PetscInt particles_add_remove;

extern PetscReal particles_perturb_factor;

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
extern PetscReal *p_add_r_strain_rate;

extern unsigned int seed;

extern PetscInt print_step_files;

extern PetscReal random_initial_strain;

// From param.txt and interfaces.txt
extern PetscInt WITH_NON_LINEAR;
extern PetscInt PLASTICITY;
extern PetscScalar *weakening_seed;

extern PetscInt *seed_layer;
extern PetscInt seed_layer_size;
extern PetscBool seed_layer_set;


extern PetscReal *strain_seed_layer;
extern PetscInt strain_seed_layer_size;
extern PetscBool strain_seed_layer_set;
extern PetscBool strain_seed_constant;

extern PetscReal epsilon_x;


PetscErrorCode _DMLocatePoints_DMDARegular_IS_3d(DM dm,Vec pos,IS *iscell)
{

	PetscInt p,n,bs,npoints,si,sj,sk,milocal,mjlocal,mklocal,mx,my,mz;
	DM dmregular;
	PetscInt *cellidx;
	PetscScalar *coor;
	PetscReal dx,dy,dz;
	PetscErrorCode ierr;
	PetscMPIInt rank;

	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
	ierr = VecGetLocalSize(pos,&n);CHKERRQ(ierr);
	ierr = VecGetBlockSize(pos,&bs);CHKERRQ(ierr);
	npoints = n/bs;

	//printf("npoints: %d\n",npoints);

	PetscMalloc1(npoints,&cellidx);

	ierr = DMGetApplicationContext(dm,(void**)&dmregular);CHKERRQ(ierr);
	ierr = DMDAGetCorners(dmregular,&si,&sj,&sk,&milocal,&mjlocal,&mklocal);CHKERRQ(ierr);
	ierr = DMDAGetInfo(dmregular,NULL,&mx,&my,&mz,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);



	dx = dx_const;
	dy = dy_const;
	dz = dz_const;

	ierr = VecGetArray(pos,&coor);CHKERRQ(ierr);
	for (p=0; p<npoints; p++) {
		PetscReal coorx,coory,coorz;
		PetscInt mi,mj,mk;

		coorx = coor[3*p];
		coory = coor[3*p+1];
		coorz = coor[3*p+2];

		mi = (PetscInt)( (coorx)/dx );
		mj = (PetscInt)( (coory)/dy );
		mk = (PetscInt)( (coorz+depth)/dz );


		cellidx[p] = DMLOCATEPOINT_POINT_NOT_FOUND;

		if ((mk >= sk) && (mk < sk + mklocal)){
			if ((mj >= sj) && (mj < sj + mjlocal)) {
				if ((mi >= si) && (mi < si + milocal)) {
					cellidx[p] = (mi-si) + (mj-sj) * milocal + (mk-sk) * milocal * mjlocal;
				}
			}
		}
		if (coorx < 0.0) cellidx[p] = DMLOCATEPOINT_POINT_NOT_FOUND;
		if (coorx >  Lx) cellidx[p] = DMLOCATEPOINT_POINT_NOT_FOUND;
		if (coory < 0.0) cellidx[p] = DMLOCATEPOINT_POINT_NOT_FOUND;
		if (coory >  Ly) cellidx[p] = DMLOCATEPOINT_POINT_NOT_FOUND;
		if (coorz < -depth) cellidx[p] = DMLOCATEPOINT_POINT_NOT_FOUND;
		if (coorz > 0.0) cellidx[p] = DMLOCATEPOINT_POINT_NOT_FOUND;
	}
	ierr = VecRestoreArray(pos,&coor);CHKERRQ(ierr);

	ierr = ISCreateGeneral(PETSC_COMM_SELF,npoints,cellidx,PETSC_OWN_POINTER,iscell);CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

PetscErrorCode DMLocatePoints_DMDARegular_3d(DM dm,Vec pos,DMPointLocationType ltype, PetscSF cellSF)
{
	IS iscell;
	PetscSFNode *cells;
	PetscInt p,bs,npoints,nfound;
	const PetscInt *boxCells;
	PetscErrorCode ierr;

	ierr = _DMLocatePoints_DMDARegular_IS_3d(dm,pos,&iscell);CHKERRQ(ierr);

	//PetscPrintf(PETSC_COMM_WORLD,"teste swarm\n");

	ierr = VecGetLocalSize(pos,&npoints);CHKERRQ(ierr);
	ierr = VecGetBlockSize(pos,&bs);CHKERRQ(ierr);
	npoints = npoints / bs;

	ierr = PetscMalloc1(npoints, &cells);CHKERRQ(ierr);
	ierr = ISGetIndices(iscell, &boxCells);CHKERRQ(ierr);

	for (p=0; p<npoints; p++) {
		cells[p].rank  = 0;
		cells[p].index = DMLOCATEPOINT_POINT_NOT_FOUND;

		cells[p].index = boxCells[p];
	}
	ierr = ISRestoreIndices(iscell, &boxCells);CHKERRQ(ierr);
	ierr = ISDestroy(&iscell);CHKERRQ(ierr);
	nfound = npoints;
	ierr = PetscSFSetGraph(cellSF, npoints, nfound, NULL, PETSC_OWN_POINTER, cells, PETSC_OWN_POINTER);CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

PetscErrorCode DMGetNeighbors_DMDARegular_3d(DM dm,PetscInt *nneighbors,const PetscMPIInt **neighbors)
{
	DM dmregular;
	PetscErrorCode ierr;

	ierr = DMGetApplicationContext(dm,(void**)&dmregular);CHKERRQ(ierr);
	ierr = DMGetNeighbors(dmregular,nneighbors,neighbors);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

PetscErrorCode SwarmViewGP_3d(DM dms,const char prefix[])
{
	PetscReal *array;
	PetscInt *iarray;
	PetscInt *layer_array;

	PetscReal *strain_fac;
	PetscInt npoints,p,bs;
	FILE *fp;
	char name[PETSC_MAX_PATH_LEN];
	PetscMPIInt rank;
	PetscErrorCode ierr;

	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
	PetscSNPrintf(name,PETSC_MAX_PATH_LEN-1,"%s-rank_new%d.txt",prefix,rank);
	fp = fopen(name,"w");
	if (!fp) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Cannot open file %s",name);
	ierr = DMSwarmGetLocalSize(dms,&npoints);CHKERRQ(ierr);
	printf("npoints = %d\n",npoints);
	ierr = DMSwarmGetField(dms,DMSwarmPICField_coor,&bs,NULL,(void**)&array);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"itag",NULL,NULL,(void**)&iarray);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"layer",NULL,NULL,(void**)&layer_array);CHKERRQ(ierr);

	ierr = DMSwarmGetField(dms,"strain_fac",NULL,NULL,(void**)&strain_fac);CHKERRQ(ierr);
	for (p=0; p<npoints; p++) {
		if (iarray[p]>9999)
			fprintf(fp,"%+1.4e %+1.4e %+1.4e %d %d %1.4e\n",
					array[3*p],array[3*p+1],array[3*p+2],
					iarray[p],layer_array[p],(double)strain_fac[p]);
	}
	ierr = DMSwarmRestoreField(dms,"itag",NULL,NULL,(void**)&iarray);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms,"layer",NULL,NULL,(void**)&layer_array);CHKERRQ(ierr);

	ierr = DMSwarmRestoreField(dms,"strain_fac",NULL,NULL,(void**)&strain_fac);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms,DMSwarmPICField_coor,&bs,NULL,(void**)&array);CHKERRQ(ierr);
	fclose(fp);
	PetscFunctionReturn(0);
}




PetscErrorCode createSwarm_3d()
{
	PetscErrorCode ierr=0;

	PetscMPIInt rank;
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);

	PetscInt bs,nlocal,p;

	PetscReal *array;
	PetscInt *iarray;
	PetscInt *layer_array;
	PetscReal *rarray;

	PetscReal *strain_array;

	PetscInt nx_part = (int)PetscCbrtReal(particles_per_ele*dx_const*dx_const/(dy_const*dz_const));
	PetscInt ny_part = (int)PetscCbrtReal(particles_per_ele*dy_const*dy_const/(dz_const*dx_const));
	PetscInt nz_part = (int)PetscCbrtReal(particles_per_ele*dz_const*dz_const/(dx_const*dy_const));

	if (nx_ppe > 0) {nx_part = nx_ppe;}
	if (ny_ppe > 0) {ny_part = ny_ppe;}
	if (nz_ppe > 0) {nz_part = nz_ppe;}

	particles_per_ele = nx_part*ny_part*nz_part;

	PetscPrintf(PETSC_COMM_WORLD,"particles per element in x:  %d\n",nx_part);
	PetscPrintf(PETSC_COMM_WORLD,"particles per element in y:  %d\n",ny_part);
	PetscPrintf(PETSC_COMM_WORLD,"particles per element in z:  %d\n",nz_part);
	PetscPrintf(PETSC_COMM_WORLD,"total particles per element: %d\n\n",particles_per_ele);

	//PetscRandom rand;


	ierr = DMShellCreate(PETSC_COMM_WORLD,&dmcell);CHKERRQ(ierr);
	ierr = DMSetApplicationContext(dmcell,(void*)da_Veloc);CHKERRQ(ierr);
	dmcell->ops->locatepoints = DMLocatePoints_DMDARegular_3d;
	dmcell->ops->getneighbors = DMGetNeighbors_DMDARegular_3d;
	PetscPrintf(PETSC_COMM_WORLD,"teste swarm externo\n");

	/* Create the swarm */
	ierr = DMCreate(PETSC_COMM_WORLD,&dms);CHKERRQ(ierr);
	ierr = DMSetType(dms,DMSWARM);CHKERRQ(ierr);
	ierr = DMSetDimension(dms,3);CHKERRQ(ierr);

	ierr = DMSwarmSetType(dms,DMSWARM_PIC);CHKERRQ(ierr);
	ierr = DMSwarmSetCellDM(dms,dmcell);CHKERRQ(ierr);

	/* init fields */
	ierr = DMSwarmRegisterPetscDatatypeField(dms,"itag",1,PETSC_INT);CHKERRQ(ierr);
	ierr = DMSwarmRegisterPetscDatatypeField(dms,"layer",1,PETSC_INT);CHKERRQ(ierr);
	ierr = DMSwarmRegisterPetscDatatypeField(dms,"geoq_fac",1,PETSC_REAL);CHKERRQ(ierr);

	ierr = DMSwarmRegisterPetscDatatypeField(dms,"strain_fac",1,PETSC_REAL);CHKERRQ(ierr);
	ierr = DMSwarmRegisterPetscDatatypeField(dms,"strain_rate_fac",1,PETSC_REAL);CHKERRQ(ierr);
	ierr = DMSwarmRegisterPetscDatatypeField(dms,"cont",1,PETSC_INT);CHKERRQ(ierr);
	ierr = DMSwarmFinalizeFieldRegister(dms);CHKERRQ(ierr);

	{
		PetscInt si,sj,sk,milocal,mjlocal,mklocal;
		PetscReal *LA_coors;
		Vec coors;
		PetscInt cnt;

		ierr = DMDAGetCorners(da_Veloc,&si,&sj,&sk,&milocal,&mjlocal,&mklocal);CHKERRQ(ierr);

		printf("%d: %d %d %d %d %d %d\n",rank,milocal,mjlocal,mklocal,si,sj,sk);

		ierr = DMGetCoordinates(da_Veloc,&coors);CHKERRQ(ierr);
		/*VecView(coors,PETSC_VIEWER_STDOUT_WORLD);*/
		ierr = VecGetArray(coors,&LA_coors);CHKERRQ(ierr);

		ierr = DMSwarmSetLocalSizes(dms,milocal*mjlocal*mklocal*(particles_per_ele),4);CHKERRQ(ierr);
		ierr = DMSwarmGetLocalSize(dms,&nlocal);CHKERRQ(ierr);

		//particles_add_remove = milocal*mjlocal*mklocal;
		if (particles_per_ele > 50) {
			particles_add_remove = milocal*mjlocal*mklocal * (PetscInt)(particles_per_ele/50);
		} else {
			particles_add_remove = milocal*mjlocal*mklocal;
		}

		printf("%d: %d\nparticles_add_remove = %d",rank,nlocal,particles_add_remove);

		ierr = DMSwarmGetField(dms,DMSwarmPICField_coor,&bs,NULL,(void**)&array);CHKERRQ(ierr);

		printf("bs = %d\n",bs);

		cnt = 0;

		for (p=0; p<nlocal/particles_per_ele; p++) {
			PetscReal px,py,pz,rx,ry,rz;

			for (int contx=0;contx<nx_part;contx++){
				for (int conty=0;conty<ny_part;conty++){
					for (int contz=0;contz<nz_part;contz++){


						rx = 1.0*(float)rand_r(&seed)/RAND_MAX-0.5;
						ry = 1.0*(float)rand_r(&seed)/RAND_MAX-0.5;
						rz = 1.0*(float)rand_r(&seed)/RAND_MAX-0.5;

						px = LA_coors[3*p+0] + (0.5+contx+rx*particles_perturb_factor)*(dx_const/nx_part);
						py = LA_coors[3*p+1] + (0.5+conty+ry*particles_perturb_factor)*(dy_const/ny_part);
						pz = LA_coors[3*p+2] + (0.5+contz+rz*particles_perturb_factor)*(dz_const/nz_part);


						if ((px>=0) && (px<=Lx) && (py>=0) && (py<=Ly) && (pz>=-depth) && (pz<=0)) {
							array[bs*cnt+0] = px;
							array[bs*cnt+1] = py;
							array[bs*cnt+2] = pz;
							cnt++;
						}
					}
				}
			}
		}

		//ierr = PetscRandomDestroy(&rand);CHKERRQ(ierr);
		ierr = DMSwarmRestoreField(dms,DMSwarmPICField_coor,&bs,NULL,(void**)&array);CHKERRQ(ierr);
		ierr = VecRestoreArray(coors,&LA_coors);CHKERRQ(ierr);
		ierr = DMSwarmSetLocalSizes(dms,cnt,4);CHKERRQ(ierr);

		ierr = DMSwarmGetLocalSize(dms,&nlocal);CHKERRQ(ierr);
		ierr = DMSwarmGetField(dms,"itag",&bs,NULL,(void**)&iarray);CHKERRQ(ierr);
		ierr = DMSwarmGetField(dms,"layer",&bs,NULL,(void**)&layer_array);CHKERRQ(ierr);
		ierr = DMSwarmGetField(dms,"strain_fac",&bs,NULL,(void**)&strain_array);CHKERRQ(ierr);

		ierr = DMSwarmGetField(dms,DMSwarmPICField_coor,&bs,NULL,(void**)&array);CHKERRQ(ierr);
		for (p=0; p<nlocal; p++) {
			iarray[p] = p%particles_per_ele;
			if (p%particles_per_ele==0){
				iarray[p] = 10000 + p/particles_per_ele + 1000000*rank;
			}
		}

		ierr = DMSwarmRestoreField(dms,"itag",&bs,NULL,(void**)&iarray);CHKERRQ(ierr);

		ierr = DMSwarmGetField(dms,"geoq_fac",&bs,NULL,(void**)&rarray);CHKERRQ(ierr);


		if (n_interfaces==0){
			for (p=0; p<nlocal; p++){
				rarray[p] = escala_viscosidade;
				//rarray_rho[p] = RHOM;
				//rarray_H[p] = H_per_mass;
				layer_array[p] = 0;
				/*rarray[p] = 1.0;
				rarray_rho[p] = 2700.0;
				rarray_H[p] = 0.0;//!!!!*/

			}
		}
		else{
			PetscReal cx,cy,cz;
			PetscReal rx,ry,rfac;
			PetscReal interp_interfaces[n_interfaces];
			PetscInt in,verif,i,j,k;

			unsigned int seed_strain;

			for (p=0; p<nlocal; p++){
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
				seed_strain=(k*Nx+i)*(k*Nx+i);

				if (i<0 || i>=Nx-1) {printf("estranho i=%d\n",i); exit(1);}
				if (j<0 || j>=Ny-1) {printf("estranho j=%d\n",j); exit(1);}


				if (i==Nx-1) i=Nx-2;
				if (j==Ny-1) j=Ny-2;


				rx = (cx-i*dx_const)/dx_const;
				ry = (cy-j*dy_const)/dy_const;


				if (rx<0 || rx>1) {printf("estranho rx=%f\n",rx); exit(1);}
				if (ry<0 || ry>1) {printf("estranho ry=%f\n",ry); exit(1);}


				for (in=0;in<n_interfaces;in++){
					rfac = (1.0-rx)*(1.0-ry);
					interp_interfaces[in] = interfaces[j*Nx+i + Nx*Ny*in] * rfac;

					rfac = (rx)*(1.0-ry);
					interp_interfaces[in] += interfaces[j*Nx+(i+1) + Nx*Ny*in] * rfac;

					rfac = (1.0-rx)*(ry);
					interp_interfaces[in] += interfaces[(j+1)*Nx+i + Nx*Ny*in] * rfac;

					rfac = (rx)*(ry);
					interp_interfaces[in] += interfaces[(j+1)*Nx+(i+1) + Nx*Ny*in] * rfac;
				}

				verif=0;
				for (in=0;in<n_interfaces && verif==0;in++){
					if (cz<interp_interfaces[in]){
						verif=1;
						rarray[p] = inter_geoq[in];
						//rarray_rho[p] = inter_rho[in];
						//rarray_H[p] = inter_H[in];
						layer_array[p] = in;
					}
				}
				if (verif==0){
					rarray[p] = inter_geoq[n_interfaces];
					//rarray_rho[p] = inter_rho[n_interfaces];
					//rarray_H[p] = inter_H[n_interfaces];
					layer_array[p] = n_interfaces;
					//printf("entrei!\n");
				}

				// Load initial strain_array with weakening_seed to or random value
				rand_r(&seed_strain);
				if (weakening_seed[layer_array[p]] >= 0) {
					strain_array[p] = weakening_seed[layer_array[p]];
				} 
				else
				{
					strain_array[p] = random_initial_strain*(float)rand_r(&seed_strain)/RAND_MAX;
				}
<<<<<<< HEAD

=======
				else 
				{
					if (WITH_NON_LINEAR==1 && PLASTICITY==1)
					{
						if (weakening_seed[layer_array[p]] >= 0) {
							strain_array[p] = weakening_seed[layer_array[p]];
						} 
					}
				}
>>>>>>> 1ffdf59e4044d78aac86fb293026d1f3d4c3d62a
				/////!!!!
				/*if (rank==0){
					printf("\n");
					for (in=0;in<n_interfaces;in++)
						printf("Interface %d %f\n",in,interp_interfaces[in]);

					printf("%lf %lf %lf\n",cx,cy,cz);


					exit(-1);
				}*/
				/////
			}
		}

		ierr = DMSwarmRestoreField(dms,"layer",&bs,NULL,(void**)&layer_array);CHKERRQ(ierr);

		ierr = DMSwarmRestoreField(dms,"geoq_fac",&bs,NULL,(void**)&rarray);CHKERRQ(ierr);
		//ierr = DMSwarmRestoreField(dms,"rho_fac",&bs,NULL,(void**)&rarray_rho);CHKERRQ(ierr);
		//ierr = DMSwarmRestoreField(dms,"H_fac",&bs,NULL,(void**)&rarray_H);CHKERRQ(ierr);

		ierr = DMSwarmRestoreField(dms,"strain_fac",&bs,NULL,(void**)&strain_array);CHKERRQ(ierr);

		ierr = DMSwarmRestoreField(dms,DMSwarmPICField_coor,&bs,NULL,(void**)&array);CHKERRQ(ierr);

	}

	ierr = DMView(dms,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

	if (print_step_files==1){
		ierr = SwarmViewGP_3d(dms,"step_0");CHKERRQ(ierr);
	}


	MPI_Barrier(PETSC_COMM_WORLD);


	ierr = PetscCalloc1(particles_add_remove ,&ppp);
	ierr = PetscCalloc1(particles_add_remove ,&p_remove);
	ierr = PetscCalloc1(particles_add_remove ,&p_i);

	ierr = PetscCalloc1(particles_add_remove*3,&p_add_coor);
	ierr = PetscCalloc1(particles_add_remove ,&p_add_r);
	//ierr = PetscCalloc1(particles_add_remove ,&p_add_r_rho);
	//ierr = PetscCalloc1(particles_add_remove ,&p_add_r_H);
	ierr = PetscCalloc1(particles_add_remove ,&p_add_i);
	ierr = PetscCalloc1(particles_add_remove ,&p_add_layer);
	ierr = PetscCalloc1(particles_add_remove ,&p_add_r_strain);
	ierr = PetscCalloc1(particles_add_remove ,&p_add_r_strain_rate);



	PetscFunctionReturn(0);

}
