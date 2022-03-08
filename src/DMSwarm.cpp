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

extern long Nx;

extern double dx_const;
extern double dz_const;

extern double Lx, depth;

extern double H_lito;
extern double escala_viscosidade;

extern PetscInt particles_per_ele;
extern PetscInt nx_ppe;
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
extern PetscInt *p_add_i;
extern PetscInt *p_add_layer;
extern PetscReal *p_add_r_strain;
extern PetscReal *p_add_r_strain_rate;

extern unsigned int seed;

extern PetscInt print_step_files;

extern PetscInt *seed_layer;
extern PetscInt seed_layer_size;
extern PetscBool seed_layer_set;


extern PetscReal *strain_seed_layer;
extern PetscInt strain_seed_layer_size;
extern PetscBool strain_seed_layer_set;

extern PetscInt checkered;

extern PetscBool plot_sediment;

extern PetscReal random_initial_strain;

extern PetscInt binary_output;



PetscErrorCode _DMLocatePoints_DMDARegular_IS(DM dm,Vec pos,IS *iscell)
{

	PetscInt p,n,bs,npoints,si,sk,milocal,mklocal,mx,mz;
	DM dmregular;
	PetscInt *cellidx;
	PetscScalar *coor;
	PetscReal dx,dz;
	PetscErrorCode ierr;
	PetscMPIInt rank;

	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
	ierr = VecGetLocalSize(pos,&n);CHKERRQ(ierr);
	ierr = VecGetBlockSize(pos,&bs);CHKERRQ(ierr);
	npoints = n/bs;

	PetscMalloc1(npoints,&cellidx);

	ierr = DMGetApplicationContext(dm,(void**)&dmregular);CHKERRQ(ierr);
	ierr = DMDAGetCorners(dmregular,&si,&sk,NULL,&milocal,&mklocal,NULL);CHKERRQ(ierr);
	ierr = DMDAGetInfo(dmregular,NULL,&mx,&mz,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);



	dx = dx_const;
	dz = dz_const;

	ierr = VecGetArray(pos,&coor);CHKERRQ(ierr);
	for (p=0; p<npoints; p++) {
		PetscReal coorx,coorz;
		PetscInt mi,mk;

		coorx = coor[2*p];
		coorz = coor[2*p+1];

		mi = (PetscInt)( (coorx)/dx );
		mk = (PetscInt)( (coorz+depth)/dz );


		cellidx[p] = DMLOCATEPOINT_POINT_NOT_FOUND;

		if ((mk >= sk) && (mk < sk + mklocal)){
			if ((mi >= si) && (mi < si + milocal)) {
				cellidx[p] = (mi-si) + (mk-sk) * milocal;
			}
		}
		if (coorx < 0.0) cellidx[p] = DMLOCATEPOINT_POINT_NOT_FOUND;
		if (coorx >  Lx) cellidx[p] = DMLOCATEPOINT_POINT_NOT_FOUND;
		if (coorz < -depth) cellidx[p] = DMLOCATEPOINT_POINT_NOT_FOUND;
		if (coorz > 0.0) cellidx[p] = DMLOCATEPOINT_POINT_NOT_FOUND;
	}
	ierr = VecRestoreArray(pos,&coor);CHKERRQ(ierr);

	ierr = ISCreateGeneral(PETSC_COMM_SELF,npoints,cellidx,PETSC_OWN_POINTER,iscell);CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

PetscErrorCode DMLocatePoints_DMDARegular(DM dm,Vec pos,DMPointLocationType ltype, PetscSF cellSF)
{
	IS iscell;
	PetscSFNode *cells;
	PetscInt p,bs,npoints,nfound;
	const PetscInt *boxCells;
	PetscErrorCode ierr;

	ierr = _DMLocatePoints_DMDARegular_IS(dm,pos,&iscell);CHKERRQ(ierr);

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

PetscErrorCode DMGetNeighbors_DMDARegular(DM dm,PetscInt *nneighbors,const PetscMPIInt **neighbors)
{
	DM dmregular;
	PetscErrorCode ierr;

	ierr = DMGetApplicationContext(dm,(void**)&dmregular);CHKERRQ(ierr);
	ierr = DMGetNeighbors(dmregular,nneighbors,neighbors);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

PetscErrorCode SwarmViewGP(DM dms,const char prefix[])
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
	if (binary_output==0){
		PetscSNPrintf(name,PETSC_MAX_PATH_LEN-1,"%s_%d.txt",prefix,rank);
		fp = fopen(name,"w");
	}
	else{
		PetscSNPrintf(name,PETSC_MAX_PATH_LEN-1,"%s_%d.bin",prefix,rank);
		fp = fopen(name,"wb");
	}


	if (!fp) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Cannot open file %s",name);
	ierr = DMSwarmGetLocalSize(dms,&npoints);CHKERRQ(ierr);
	//printf("npoints = %d\n",npoints);
	ierr = DMSwarmGetField(dms,DMSwarmPICField_coor,&bs,NULL,(void**)&array);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"itag",NULL,NULL,(void**)&iarray);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"layer",NULL,NULL,(void**)&layer_array);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"strain_fac",NULL,NULL,(void**)&strain_fac);CHKERRQ(ierr);
	if (binary_output==0){
		for (p=0; p<npoints; p++) {
			if (iarray[p]>9999 || (PETSC_TRUE == plot_sediment && layer_array[p] == n_interfaces - 1))
				fprintf(fp,"%+1.5e %+1.5e %d %d %1.4e\n",
						array[2*p],array[2*p+1],
						iarray[p],layer_array[p],(double)strain_fac[p]);
		}
	}
	else{
		int cont_print=0;
		for (p=0;p<npoints; p++) {
			if (iarray[p]>9999 || (PETSC_TRUE == plot_sediment && layer_array[p] == n_interfaces - 1)) cont_print++;
		}
		fwrite(&cont_print,sizeof(cont_print),1,fp);
		for (p=0;p<npoints; p++) {
			if (iarray[p]>9999 || (PETSC_TRUE == plot_sediment && layer_array[p] == n_interfaces - 1)){
				fwrite(&array[2*p],sizeof(array[0]),2,fp);
			}
		}
		for (p=0;p<npoints; p++) {
			if (iarray[p]>9999 || (PETSC_TRUE == plot_sediment && layer_array[p] == n_interfaces - 1)){
				fwrite(&iarray[p],sizeof(iarray[0]),1,fp);
			}
		}
		for (p=0;p<npoints; p++) {
			if (iarray[p]>9999 || (PETSC_TRUE == plot_sediment && layer_array[p] == n_interfaces - 1)){
				fwrite(&layer_array[p],sizeof(layer_array[0]),1,fp);
			}
		}
		for (p=0;p<npoints; p++) {
			if (iarray[p]>9999 || (PETSC_TRUE == plot_sediment && layer_array[p] == n_interfaces - 1)){
				fwrite(&strain_fac[p],sizeof(strain_fac[0]),1,fp);
			}
		}
	}
	ierr = DMSwarmRestoreField(dms,"itag",NULL,NULL,(void**)&iarray);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms,"layer",NULL,NULL,(void**)&layer_array);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms,"strain_fac",NULL,NULL,(void**)&strain_fac);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms,DMSwarmPICField_coor,&bs,NULL,(void**)&array);CHKERRQ(ierr);
	fclose(fp);
	PetscFunctionReturn(0);
}




PetscErrorCode createSwarm()
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


	PetscInt nz_part = (int)PetscSqrtReal(particles_per_ele*dz_const/dx_const);
	PetscInt nx_part = (int)(particles_per_ele/nz_part);

	if (nx_ppe > 0) {
		nx_part = nx_ppe;
	}

	if (nz_ppe > 0) {
		nz_part = nz_ppe;
	}

	particles_per_ele = nx_part*nz_part;


	PetscPrintf(PETSC_COMM_WORLD,"particles per element in x:  %d\n",nx_part);
	PetscPrintf(PETSC_COMM_WORLD,"particles per element in z:  %d\n",nz_part);
	PetscPrintf(PETSC_COMM_WORLD,"total particles per element: %d\n\n",particles_per_ele);

	ierr = DMShellCreate(PETSC_COMM_WORLD,&dmcell);CHKERRQ(ierr);
	ierr = DMSetApplicationContext(dmcell,(void*)da_Veloc);CHKERRQ(ierr);
	dmcell->ops->locatepoints = DMLocatePoints_DMDARegular;
	dmcell->ops->getneighbors = DMGetNeighbors_DMDARegular;
	//PetscPrintf(PETSC_COMM_WORLD,"teste swarm externo\n");

	/* Create the swarm */
	ierr = DMCreate(PETSC_COMM_WORLD,&dms);CHKERRQ(ierr);
	ierr = DMSetType(dms,DMSWARM);CHKERRQ(ierr);
	ierr = DMSetDimension(dms,2);CHKERRQ(ierr);

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
		PetscInt si,sk,milocal,mklocal;
		PetscReal *LA_coors;
		Vec coors;
		PetscInt cnt;

		ierr = DMDAGetCorners(da_Veloc,&si,&sk,NULL,&milocal,&mklocal,NULL);CHKERRQ(ierr);


		ierr = DMGetCoordinates(da_Veloc,&coors);CHKERRQ(ierr);

		ierr = VecGetArray(coors,&LA_coors);CHKERRQ(ierr);

		ierr = DMSwarmSetLocalSizes(dms,milocal*mklocal*(particles_per_ele),4);CHKERRQ(ierr);
		ierr = DMSwarmGetLocalSize(dms,&nlocal);CHKERRQ(ierr);

		if (particles_per_ele > 50) {
			particles_add_remove = milocal*mklocal * (PetscInt)(particles_per_ele/50);
		} else {
			particles_add_remove = milocal*mklocal;
		}

		ierr = DMSwarmGetField(dms,DMSwarmPICField_coor,&bs,NULL,(void**)&array);CHKERRQ(ierr);

		cnt = 0;

		for (p=0; p<nlocal/particles_per_ele; p++) {
			PetscReal px,pz,rx,rz;

			for (int contx=0;contx<nx_part;contx++){
				for (int contz=0;contz<nz_part;contz++){
					rx = 1.0*(float)rand_r(&seed)/RAND_MAX-0.5;
					rz = 1.0*(float)rand_r(&seed)/RAND_MAX-0.5;

					px = LA_coors[2*p+0] + (0.5+contx+rx*particles_perturb_factor)*(dx_const/nx_part);
					pz = LA_coors[2*p+1] + (0.5+contz+rz*particles_perturb_factor)*(dz_const/nz_part);
					if ((px>=0) && (px<=Lx) && (pz>=-depth) && (pz<=0)) {
						array[bs*cnt+0] = px;
						array[bs*cnt+1] = pz;
						cnt++;
					}
				}
			}

		}

		ierr = DMSwarmRestoreField(dms,DMSwarmPICField_coor,&bs,NULL,(void**)&array);CHKERRQ(ierr);
		ierr = VecRestoreArray(coors,&LA_coors);CHKERRQ(ierr);
		ierr = DMSwarmSetLocalSizes(dms,cnt,4);CHKERRQ(ierr);

		ierr = DMSwarmGetLocalSize(dms,&nlocal);CHKERRQ(ierr);
		ierr = DMSwarmGetField(dms,"itag",&bs,NULL,(void**)&iarray);CHKERRQ(ierr);
		ierr = DMSwarmGetField(dms,"layer",&bs,NULL,(void**)&layer_array);CHKERRQ(ierr);
		ierr = DMSwarmGetField(dms,"strain_fac",&bs,NULL,(void**)&strain_array);CHKERRQ(ierr);

		ierr = DMSwarmGetField(dms,DMSwarmPICField_coor,&bs,NULL,(void**)&array);CHKERRQ(ierr);
		if (checkered==0){
			for (p=0; p<nlocal; p++) {
				iarray[p] = p%particles_per_ele;
				if (p%particles_per_ele==0){
					iarray[p] = 10000 + p/particles_per_ele + 1000000*rank;
				}
			}
		}
		else{
			for (p=0; p<nlocal/particles_per_ele; p++) {
				int cont_part_ele=0;
				for (int contx=0;contx<nx_part;contx++){
					for (int contz=0;contz<nz_part;contz++){
						iarray[p*particles_per_ele+cont_part_ele] = 1;
						if (contx==0 || contz==0){
							iarray[p*particles_per_ele+cont_part_ele] = 10000 + p + 1000000*rank;
						}
						cont_part_ele++;
					}
				}
			}
		}

		ierr = DMSwarmRestoreField(dms,"itag",&bs,NULL,(void**)&iarray);CHKERRQ(ierr);

		ierr = DMSwarmGetField(dms,"geoq_fac",&bs,NULL,(void**)&rarray);CHKERRQ(ierr);

		if (n_interfaces==0){
			for (p=0; p<nlocal; p++){
				rarray[p] = escala_viscosidade;
				layer_array[p] = 0;
			}
		}
		else{
			PetscReal cx,cz;
			PetscReal rx,rfac;
			PetscReal interp_interfaces[n_interfaces];
			PetscInt in,verif,i,k;

			unsigned int seed_strain;

			for (p=0; p<nlocal; p++){
				cx = array[2*p];
				cz = array[2*p+1];

				i = (int)(cx/dx_const);
				k = (int)((cz+depth)/dz_const);
				seed_strain=(k*Nx+i)*(k*Nx+i);

				if (i<0 || i>=Nx-1) {printf("weird: i=%d create cx = %lf\n",i,cx); exit(1);}

				if (i==Nx-1) i=Nx-2;

				rx = (cx-i*dx_const)/dx_const;

				if (rx<0 || rx>1) {printf("weird rx=%f\n",rx); exit(1);}

				for (in=0;in<n_interfaces;in++){
					rfac = (1.0-rx);
					interp_interfaces[in] = interfaces[i + Nx*in] * rfac;

					rfac = (rx);
					interp_interfaces[in] += interfaces[(i+1) + Nx*in] * rfac;

				}

				verif=0;
				for (in=0;in<n_interfaces && verif==0;in++){
					if (cz<interp_interfaces[in]){
						verif=1;
						rarray[p] = inter_geoq[in];
						layer_array[p] = in;
					}
				}
				if (verif==0){
					rarray[p] = inter_geoq[n_interfaces];
					layer_array[p] = n_interfaces;
				}

				rand_r(&seed_strain);
				strain_array[p]=random_initial_strain*(float)rand_r(&seed_strain)/RAND_MAX;

				if (seed_layer_set == PETSC_TRUE) {
					if (seed_layer_size == 1 && layer_array[p] == seed_layer[0]) {
						strain_array[p] = strain_seed_layer[0];
					}
					else {
						for (int k = 0; k < seed_layer_size; k++) {
							if (layer_array[p] == seed_layer[k]) {
								strain_array[p] = strain_seed_layer[k];
							}
						}
					}
				}


			}
		}

		ierr = DMSwarmRestoreField(dms,"layer",&bs,NULL,(void**)&layer_array);CHKERRQ(ierr);

		ierr = DMSwarmRestoreField(dms,"geoq_fac",&bs,NULL,(void**)&rarray);CHKERRQ(ierr);

		ierr = DMSwarmRestoreField(dms,"strain_fac",&bs,NULL,(void**)&strain_array);CHKERRQ(ierr);

		ierr = DMSwarmRestoreField(dms,DMSwarmPICField_coor,&bs,NULL,(void**)&array);CHKERRQ(ierr);

	}

	ierr = DMView(dms,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

	if (print_step_files==1){
		ierr = SwarmViewGP(dms,"step_0");CHKERRQ(ierr);
	}


	MPI_Barrier(PETSC_COMM_WORLD);


	ierr = PetscCalloc1(particles_add_remove ,&ppp);
	ierr = PetscCalloc1(particles_add_remove ,&p_remove);
	ierr = PetscCalloc1(particles_add_remove ,&p_i);

	ierr = PetscCalloc1(particles_add_remove*2,&p_add_coor);
	ierr = PetscCalloc1(particles_add_remove ,&p_add_r);
	ierr = PetscCalloc1(particles_add_remove ,&p_add_i);
	ierr = PetscCalloc1(particles_add_remove ,&p_add_layer);
	ierr = PetscCalloc1(particles_add_remove ,&p_add_r_strain);
	ierr = PetscCalloc1(particles_add_remove ,&p_add_r_strain_rate);



	PetscFunctionReturn(0);

}
