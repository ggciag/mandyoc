#include <petscsf.h>
#include <petscdm.h>
#include <petscksp.h>
#include <petscdmda.h>
#include <petscdmshell.h>
#include <petscdmswarm.h>
#include <petsc/private/dmimpl.h>
#include <petscmath.h>
//#include <petscsys.h>

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

extern double H_per_mass;

extern double h_air;

extern double RHOM;

extern int n_interfaces;
extern PetscScalar *interfaces;

extern PetscScalar *inter_rho;
extern PetscScalar *inter_geoq;
extern PetscScalar *inter_H;

extern PetscInt particles_add_remove;

extern PetscInt *ppp;
extern PetscInt *p_remove;
extern PetscInt *p_i;

extern PetscReal *p_add_coor;
extern PetscReal *p_add_r;
extern PetscInt *p_add_i;
extern PetscInt *p_add_layer;
extern PetscReal *p_add_r_strain;

extern unsigned int seed;

extern PetscInt print_step_files;


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
	
	//printf("npoints: %d\n",npoints);
	
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
	FILE *fp;//,*fp2;
	char name[PETSC_MAX_PATH_LEN];
	PetscMPIInt rank;
	PetscErrorCode ierr;
	
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
	PetscSNPrintf(name,PETSC_MAX_PATH_LEN-1,"%s-rank_new%d.txt",prefix,rank);
	fp = fopen(name,"w");
	//PetscSNPrintf(name,PETSC_MAX_PATH_LEN-1,"surf_%s-rank_new%d.txt",prefix,rank);
	//fp2 = fopen(name,"w");
	if (!fp) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Cannot open file %s",name);
	ierr = DMSwarmGetLocalSize(dms,&npoints);CHKERRQ(ierr);
	printf("npoints = %d\n",npoints);
	ierr = DMSwarmGetField(dms,DMSwarmPICField_coor,&bs,NULL,(void**)&array);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"itag",NULL,NULL,(void**)&iarray);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"layer",NULL,NULL,(void**)&layer_array);CHKERRQ(ierr);
	ierr = DMSwarmGetField(dms,"strain_fac",NULL,NULL,(void**)&strain_fac);CHKERRQ(ierr);
	for (p=0; p<npoints; p++) {
		if (iarray[p]>999)
			fprintf(fp,"%+1.4e %+1.4e %d %d %1.4e\n",
					array[2*p],array[2*p+1],
					iarray[p],layer_array[p],(double)strain_fac[p]);
		//if ((array[2*p+1]>-105.0E3)&&(array[2*p+1]<-95.0E3))
		//	fprintf(fp2,"%+1.4e %+1.4e %d %d %1.4e\n",
		//			array[2*p],array[2*p+1],
		//			iarray[p],layer_array[p],(double)strain_fac[p]);
	}
	ierr = DMSwarmRestoreField(dms,"itag",NULL,NULL,(void**)&iarray);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms,"layer",NULL,NULL,(void**)&layer_array);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms,"strain_fac",NULL,NULL,(void**)&strain_fac);CHKERRQ(ierr);
	ierr = DMSwarmRestoreField(dms,DMSwarmPICField_coor,&bs,NULL,(void**)&array);CHKERRQ(ierr);
	fclose(fp);
	//fclose(fp2);
	PetscFunctionReturn(0);
}




PetscErrorCode createSwarm()
{
	PetscErrorCode ierr=0;

	PetscMPIInt rank;
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
	
	PetscInt bs,nlocal,p,cont;
	
	PetscReal *array;
	PetscInt *iarray;
	PetscInt *layer_array;
	PetscReal *rarray;
	
	//PetscRandom rand;
	
	
	ierr = DMShellCreate(PETSC_COMM_WORLD,&dmcell);CHKERRQ(ierr);
	ierr = DMSetApplicationContext(dmcell,(void*)da_Veloc);CHKERRQ(ierr);
	dmcell->ops->locatepoints = DMLocatePoints_DMDARegular;
	dmcell->ops->getneighbors = DMGetNeighbors_DMDARegular;
	PetscPrintf(PETSC_COMM_WORLD,"teste swarm externo\n");
	
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
	ierr = DMSwarmRegisterPetscDatatypeField(dms,"cont",1,PETSC_INT);CHKERRQ(ierr);
	ierr = DMSwarmFinalizeFieldRegister(dms);CHKERRQ(ierr);
	
	{
		PetscInt si,sk,milocal,mklocal;
		PetscReal *LA_coors;
		Vec coors;
		PetscInt cnt;
		
		ierr = DMDAGetCorners(da_Veloc,&si,&sk,NULL,&milocal,&mklocal,NULL);CHKERRQ(ierr);
		
		printf("%d: %d %d %d %d\n",rank,milocal,mklocal,si,sk);
		
		ierr = DMGetCoordinates(da_Veloc,&coors);CHKERRQ(ierr);
		/*VecView(coors,PETSC_VIEWER_STDOUT_WORLD);*/
		ierr = VecGetArray(coors,&LA_coors);CHKERRQ(ierr);
		
		ierr = DMSwarmSetLocalSizes(dms,milocal*mklocal*(particles_per_ele),4);CHKERRQ(ierr);
		ierr = DMSwarmGetLocalSize(dms,&nlocal);CHKERRQ(ierr);
		
		particles_add_remove = milocal*mklocal;
		
		printf("%d: %d\nparticles_add_remove = %d",rank,nlocal,particles_add_remove);
		
		ierr = DMSwarmGetField(dms,DMSwarmPICField_coor,&bs,NULL,(void**)&array);CHKERRQ(ierr);
		
		printf("bs = %d\n",bs);
		
		
		/*if (rank==0){
			for (int cont=0;cont<10;cont++)
				printf("array %f\n",array[cont]);
		}*/
		
		cnt = 0;
		//ierr = PetscRandomCreate(PETSC_COMM_SELF,&rand);CHKERRQ(ierr);
		//ierr = PetscRandomSetType(rand,PETSCRAND48);CHKERRQ(ierr);
		//ierr = PetscRandomSetInterval(rand,-1.0,1.0);CHKERRQ(ierr);
		
		for (p=0; p<nlocal/particles_per_ele; p++) {
			PetscReal px,pz,rx,rz;
			
			for (cont=0;cont<particles_per_ele;cont++){
				
				//ierr = PetscRandomGetValueReal(rand,&rx);CHKERRQ(ierr);
				//ierr = PetscRandomGetValueReal(rand,&ry);CHKERRQ(ierr);
				//ierr = PetscRandomGetValueReal(rand,&rz);CHKERRQ(ierr);
				rx = 2.0*(float)rand_r(&seed)/RAND_MAX-1.0;
				rz = 2.0*(float)rand_r(&seed)/RAND_MAX-1.0;
				
				px = LA_coors[2*p+0] + (0.5*rx+0.5)*dx_const;
				pz = LA_coors[2*p+1] + (0.5*rz+0.5)*dz_const;
				
				if ((px>=0) && (px<=Lx) && (pz>=-depth) && (pz<=0)) {
					array[bs*cnt+0] = px;
					array[bs*cnt+1] = pz;
					cnt++;
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
		
		ierr = DMSwarmGetField(dms,DMSwarmPICField_coor,&bs,NULL,(void**)&array);CHKERRQ(ierr);
		for (p=0; p<nlocal; p++) {
			//iarray[p] = (PetscInt)rank;
			iarray[p] = p%particles_per_ele;
			if (p%particles_per_ele==0){
				iarray[p] = 1000 + p + 1000000*rank;
			}
		}
		
		ierr = DMSwarmRestoreField(dms,"itag",&bs,NULL,(void**)&iarray);CHKERRQ(ierr);
		
		ierr = DMSwarmGetField(dms,"geoq_fac",&bs,NULL,(void**)&rarray);CHKERRQ(ierr);
		
		if (n_interfaces==0){
			for (p=0; p<nlocal; p++){
				rarray[p] = escala_viscosidade;
				layer_array[p] = 0;
				/*rarray[p] = 1.0;!!!!*/
				
			}
		}
		else{
			PetscReal cx,cz;
			PetscReal rx,rfac;
			PetscReal interp_interfaces[n_interfaces];
			PetscInt in,verif,i;
			
			for (p=0; p<nlocal; p++){
				cx = array[2*p];
				cz = array[2*p+1];

				i = (int)(cx/dx_const);
				
				if (i<0 || i>=Nx-1) {printf("estranho i=%d\n",i); exit(1);}
				
				if (i==Nx-1) i=Nx-2;
				
				rx = (cx-i*dx_const)/dx_const;
				
				if (rx<0 || rx>1) {printf("estranho rx=%f\n",rx); exit(1);}
				
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
					//printf("entrei!\n");
				}
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
	
	
	
	PetscFunctionReturn(0);
	
}
