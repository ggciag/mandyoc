#include <petscksp.h>

extern long Nx;
extern long Ny;
extern long Nz;

extern long V_NE;
extern long V_GT;
extern long V_GN;

extern double r06;
extern double r8p9;
extern double r5p9;


extern double dx_const;
extern double dy_const;
extern double dz_const;

extern Vec v_vec;
extern Vec v_vec_fut;

extern PetscReal *v_vec_aux_ele;

extern double dt_calor_sec;

extern long GaussQuad;

extern double H_per_mass;
extern double c_heat_capacity;


extern PetscReal *NT;
extern PetscReal *NT_x;
extern PetscReal *NT_y;
extern PetscReal *NT_z;


PetscErrorCode montaKeThermal_simplif(double *Ke_local,double *Ke,double kappa_eff){
	long i,j,c;
	
	PetscErrorCode ierr=0;
	

	
	//check for 3D
	/*for (c=0;c<V_NE;c++){
		indice_aux_vec_ele[c*V_GN+0]=Hexa_thermal[t*V_NE+c]*3;
		indice_aux_vec_ele[c*V_GN+1]=Hexa_thermal[t*V_NE+c]*3+1;
		indice_aux_vec_ele[c*V_GN+2]=Hexa_thermal[t*V_NE+c]*3+2;
	}
	if (cond_fut==0) VecGetValues(v_vec,V_GT,indice_aux_vec_ele,v_vec_aux_ele);
	else VecGetValues(v_vec_fut,V_GT,indice_aux_vec_ele,v_vec_aux_ele);*/
	//for (c=0;c<V_GT;c++) v_vec_aux_ele[c]=0.0;
	///////////
	
	double volume;
	
	double vv[2],vvx[V_NE],vvz[V_NE];
	
	vv[0]=0.0;
	vv[1]=0.0;


	for (c=0;c<V_NE;c++){
		vvx[c]=v_vec_aux_ele[c*V_GN  ];
		vvz[c]=v_vec_aux_ele[c*V_GN+1];
		
		vv[0]+=vvx[c];
		vv[1]+=vvz[c];
	}
	vv[0]/=V_NE;
	vv[1]/=V_NE;
	
	volume=dx_const*dz_const;
	
	double kx,kz;
	
	for (i=0;i<V_NE;i++){
		for (j=0;j<V_NE;j++){
			Ke_local[i*V_NE+j]=0.0;
		}
	}
	

	double mod_vv = sqrt(vv[0]*vv[0]+vv[1]*vv[1]);
	double unit_vv[2];
	unit_vv[0]=vv[0]/mod_vv;
	unit_vv[1]=vv[1]/mod_vv;
	
	
	double eps[2],a_aux[2];

	a_aux[0]= vv[0]*dx_const/(kappa_eff*2.0);
	a_aux[1]= vv[1]*dz_const/(kappa_eff*2.0);
	
	
	for (i=0;i<2;i++){
		if (a_aux[i]<-1.0)		eps[i]=-1.-1./a_aux[i];
		else{
			if (a_aux[i]<=1.0)	eps[i]=0;
			else				eps[i]=1.-1./a_aux[i];
		}
	}
	
	
	double kappa_arti=0.0;
	
	a_aux[0]= vv[0]*dx_const;
	a_aux[1]= vv[1]*dz_const;
	
	for (i=0;i<2;i++)	kappa_arti += eps[i]*a_aux[i]/2.0;
	
	double Hx,Hz,prodH;
	
	long point=0;
	
	

	
	double NT_V[V_NE],NT_u[V_NE];
	
	for (kz=-r06; kz<=r06; kz+=r06){
		if (kz==0) Hz=r8p9;
		else Hz=r5p9;
	
		for (kx=-r06; kx<=r06; kx+=r06){
			if (kx==0) Hx=r8p9;
			else Hx=r5p9;
			
			
			prodH = Hx*Hz;
			
			for (i=0;i<V_NE;i++) NT_u[i] = (kappa_arti*(NT_x[i+point]*unit_vv[0]+NT_z[i+point]*unit_vv[1])/mod_vv);
			for (j=0;j<V_NE;j++) NT_V[j] = prodH*volume*(NT_x[j+point]*vv[0]+NT_z[j+point]*vv[1]);
			
			
			for (i=0;i<V_NE;i++){
				for (j=0;j<V_NE;j++){
					
					
					Ke_local[i*V_NE+j]+=  NT[i+point]*NT_V[j];
					
					if (mod_vv>0) Ke_local[i*V_NE+j]+= NT_u[i]*NT_V[j];
					
				}
			}
			
			
			point +=V_NE;
			
		}
	
	}
	

	
	//exit(0);
	
	for (i=0;i<V_NE*V_NE;i++) Ke_local[i]+=Ke[i]*kappa_eff;
	
	return (ierr);
	
	
	
}


PetscErrorCode montaKeThermal_simplif3d(double *Ke_local,double *Ke,double kappa_eff){
	long i,j,c;
	
	PetscErrorCode ierr=0;
	

	
	//modificar!!!!
	/*for (c=0;c<V_NE;c++){
		indice_aux_vec_ele[c*V_GN+0]=Hexa_thermal[t*V_NE+c]*3;
		indice_aux_vec_ele[c*V_GN+1]=Hexa_thermal[t*V_NE+c]*3+1;
		indice_aux_vec_ele[c*V_GN+2]=Hexa_thermal[t*V_NE+c]*3+2;
	}
	if (cond_fut==0) VecGetValues(v_vec,V_GT,indice_aux_vec_ele,v_vec_aux_ele);
	else VecGetValues(v_vec_fut,V_GT,indice_aux_vec_ele,v_vec_aux_ele);*/
	//for (c=0;c<V_GT;c++) v_vec_aux_ele[c]=0.0;
	///////////
	
	double volume;
	
	double vv[3],vvx[V_NE],vvy[V_NE],vvz[V_NE];
	
	vv[0]=0.0;
	vv[1]=0.0;
	vv[2]=0.0;

	for (c=0;c<V_NE;c++){
		vvx[c]=v_vec_aux_ele[c*V_GN  ];
		vvy[c]=v_vec_aux_ele[c*V_GN+1];
		vvz[c]=v_vec_aux_ele[c*V_GN+2];
		
		vv[0]+=vvx[c];
		vv[1]+=vvy[c];
		vv[2]+=vvz[c];
	}
	vv[0]/=8.0;
	vv[1]/=8.0;
	vv[2]/=8.0;
	
	volume=dx_const*dy_const*dz_const;
	
	double kx,ky,kz;
	
	for (i=0;i<V_NE;i++){
		for (j=0;j<V_NE;j++){
			Ke_local[i*V_NE+j]=0.0;
		}
	}
	

	double mod_vv = sqrt(vv[0]*vv[0]+vv[1]*vv[1]+vv[2]*vv[2]);
	double unit_vv[3];
	unit_vv[0]=vv[0]/mod_vv;
	unit_vv[1]=vv[1]/mod_vv;
	unit_vv[2]=vv[2]/mod_vv;
	
	
	double eps[3],a_aux[3];
	
	a_aux[0]= vv[0]*dx_const/(kappa_eff*2.0);
	a_aux[1]= vv[1]*dy_const/(kappa_eff*2.0);
	a_aux[2]= vv[2]*dz_const/(kappa_eff*2.0);
	
	
	for (i=0;i<3;i++){
		if (a_aux[i]<-1.0)		eps[i]=-1.-1./a_aux[i];
		else{
			if (a_aux[i]<=1.0)	eps[i]=0;
			else				eps[i]=1.-1./a_aux[i];
		}
	}
	
	
	double kappa_arti=0.0;
	
	a_aux[0]= vv[0]*dx_const;
	a_aux[1]= vv[1]*dy_const;
	a_aux[2]= vv[2]*dz_const;
	
	for (i=0;i<3;i++)	kappa_arti += eps[i]*a_aux[i]/2.0;
	
	double Hx,Hy,Hz,prodH;
	
	long point=0;
	
	

	
	double NT_V[V_NE],NT_u[V_NE];
	
	for (kz=-r06; kz<=r06; kz+=r06){
		if (kz==0) Hz=r8p9;
		else Hz=r5p9;
		for (ky=-r06; ky<=r06; ky+=r06){
			if (ky==0) Hy=r8p9;
			else Hy=r5p9;
			for (kx=-r06; kx<=r06; kx+=r06){
				if (kx==0) Hx=r8p9;
				else Hx=r5p9;
				
				
				prodH = Hx*Hy*Hz;
				
				for (i=0;i<V_NE;i++) NT_u[i] = (kappa_arti*(NT_x[i+point]*unit_vv[0]+NT_y[i+point]*unit_vv[1]+NT_z[i+point]*unit_vv[2])/mod_vv);
				for (j=0;j<V_NE;j++) NT_V[j] = prodH*volume*(NT_x[j+point]*vv[0]+NT_y[j+point]*vv[1]+NT_z[j+point]*vv[2]);
				
				
				for (i=0;i<V_NE;i++){
					for (j=0;j<V_NE;j++){
						
						
						Ke_local[i*V_NE+j]+=  NT[i+point]*NT_V[j];
						
						if (mod_vv>0) Ke_local[i*V_NE+j]+= NT_u[i]*NT_V[j];
						
					}
				}
				
				
				point +=V_NE;
				
			}
		}
	}
	

	
	//exit(0);
	
	for (i=0;i<V_NE*V_NE;i++) Ke_local[i]+=Ke[i]*kappa_eff;
	
	return (ierr);

}







PetscErrorCode montaKeThermal_general(PetscReal *Ke, PetscReal *Me, PetscReal *Fe)
{
	long i,j;
	
	PetscErrorCode ierr=0;
	
	
	
	double volume;
	
	
	volume=dx_const*dz_const;
	
	double kx,kz;
	
	double ex,ez;
	
	long cont;
	
	double N[V_NE];
	double N_x[V_NE];
	double N_z[V_NE];
	
	for (i=0;i<V_NE;i++){
		Fe[i]=0.0;
		for (j=0;j<V_NE;j++){
			Me[i*V_NE+j]=0.0;
			Ke[i*V_NE+j]=0.0;
		}
	}
	
	
	double Hx,Hz,prodH;
	
	long point=0;
	
	for (kz=-r06; kz<=r06; kz+=r06){
		if (kz==0) Hz=r8p9;
		else Hz=r5p9;
		
		for (kx=-r06; kx<=r06; kx+=r06){
			if (kx==0) Hx=r8p9;
			else Hx=r5p9;
			
			
			prodH = Hx*Hz;
			cont=0;
			for (ez=-1.;ez<=1.;ez+=2.){
				for (ex=-1.;ex<=1.;ex+=2.){
					N[cont]=(1+ex*kx)*(1+ez*kz)/4.0;
					N_x[cont]=ex*(1+ez*kz)/2.0/dx_const;
					N_z[cont]=(1+ex*kx)*ez/2.0/dz_const;
					
					cont++;
				}
			}
			
			
			for (i=0;i<V_NE;i++){
				
				NT[i+point*V_NE] = N[i];
				NT_x[i+point*V_NE] = N_x[i];
				NT_z[i+point*V_NE] = N_z[i];
				
				if (c_heat_capacity>0.0) Fe[i] += prodH*volume*N[i]/c_heat_capacity;///H_per_mass removed from here: ok
				
				for (j=0;j<V_NE;j++){
					Me[i*V_NE+j]+=prodH*volume*N[i]*N[j];
					
					
					Ke[i*V_NE+j]+= prodH*volume*(N_x[i]*N_x[j]+N_z[i]*N_z[j]);
				}
			}
			
			point++;
			
		}
		
	}
	
	return (ierr);
	
}


PetscErrorCode montaKeThermal_general3d(PetscReal *Ke, PetscReal *Me, PetscReal *Fe)
{
	long i,j;
	
	PetscErrorCode ierr=0;
	
	
	
	double volume;
	
	
	volume=dx_const*dy_const*dz_const;
	
	double kx,ky,kz;
	
	double ex,ey,ez;
	
	long cont;
	
	double N[V_NE];
	double N_x[V_NE];
	double N_y[V_NE];
	double N_z[V_NE];
	
	for (i=0;i<V_NE;i++){
		Fe[i]=0.0;
		for (j=0;j<V_NE;j++){
			Me[i*V_NE+j]=0.0;
			Ke[i*V_NE+j]=0.0;
		}
	}
	
	

	
	double Hx,Hy,Hz,prodH;
	
	long point=0;
	
	for (kz=-r06; kz<=r06; kz+=r06){
		if (kz==0) Hz=r8p9;
		else Hz=r5p9;
		for (ky=-r06; ky<=r06; ky+=r06){
			if (ky==0) Hy=r8p9;
			else Hy=r5p9;
			for (kx=-r06; kx<=r06; kx+=r06){
				if (kx==0) Hx=r8p9;
				else Hx=r5p9;
				
				
				prodH = Hx*Hy*Hz;
				cont=0;
				for (ez=-1.;ez<=1.;ez+=2.){
					for (ey=-1.;ey<=1.;ey+=2.){
						for (ex=-1.;ex<=1.;ex+=2.){
							N[cont]=(1+ex*kx)*(1+ey*ky)*(1+ez*kz)/8.0;
							N_x[cont]=ex*(1+ey*ky)*(1+ez*kz)/4.0/dx_const;
							N_y[cont]=(1+ex*kx)*ey*(1+ez*kz)/4.0/dy_const;
							N_z[cont]=(1+ex*kx)*(1+ey*ky)*ez/4.0/dz_const;
							
							cont++;
						}
					}
				}
				
				
				for (i=0;i<V_NE;i++){
					
					NT[i+point*V_NE] = N[i];
					NT_x[i+point*V_NE] = N_x[i];
					NT_y[i+point*V_NE] = N_y[i];
					NT_z[i+point*V_NE] = N_z[i];
					
					if (c_heat_capacity>0.0) Fe[i] += prodH*volume*N[i]/c_heat_capacity;///H_per_mass removido daqui!!!ok
					
					for (j=0;j<V_NE;j++){
						Me[i*V_NE+j]+=prodH*volume*N[i]*N[j];
						
						
						Ke[i*V_NE+j]+= prodH*volume*(N_x[i]*N_x[j]+N_y[i]*N_y[j]+N_z[i]*N_z[j]);
					}
				}
				
				point++;
				
			}
		}
	}
	
	return (ierr);
}
