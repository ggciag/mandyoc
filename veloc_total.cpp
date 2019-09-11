#include <petscksp.h>
#include <petscmath.h>

PetscErrorCode build_veloc_3d();

PetscErrorCode solve_veloc_3d();

extern Vec Veloc;
extern Vec Veloc_fut;
extern Vec Veloc_step1;
extern Vec Veloc_step2;

extern PetscInt WITH_NON_LINEAR;

extern int tcont;
extern double visc_MAX;
extern double visc_MIN;


PetscErrorCode veloc_total()
{
	PetscErrorCode ierr;
	
	PetscInt n;
	
	VecCopy(Veloc_fut,Veloc);
	
	if (tcont==0){
		double visc_MAX_aux = visc_MAX;
		double visc_MIN_aux = visc_MIN;
		
		visc_MAX = 1.0E23;
		visc_MIN = 1.0E22;
		
		while (visc_MAX<visc_MAX_aux && visc_MIN>visc_MIN_aux){
			ierr = build_veloc_3d();CHKERRQ(ierr);
			
			ierr = solve_veloc_3d();CHKERRQ(ierr);
			
			visc_MAX*=5;
			visc_MIN/=5;
			
			if (visc_MAX>visc_MAX_aux) visc_MAX = visc_MAX_aux;
			if (visc_MIN<visc_MIN_aux) visc_MIN = visc_MIN_aux;
		}
	}
	
	
	
	if (WITH_NON_LINEAR==0){
		ierr = build_veloc_3d();CHKERRQ(ierr);
		
		ierr = solve_veloc_3d();CHKERRQ(ierr);
	}
	else {
	
		VecCopy(Veloc_fut,Veloc_step1);
		
		VecGetSize(Veloc_fut,&n);
		
		
		PetscReal VM1,VM2,sig1,sig2,vivi;
		
		PetscReal Xi=50000.0;
		
		for (int step=0; step<700 && Xi>1.0E-14; step++){
			
			ierr = build_veloc_3d();CHKERRQ(ierr);
			
			ierr = solve_veloc_3d();CHKERRQ(ierr);
			
			VecCopy(Veloc_fut,Veloc_step2);
			
			
			VecSum(Veloc_step1,&VM1);
			VM1/=n;
			VecShift(Veloc_step1,-VM1);

			
			VecSum(Veloc_step2,&VM2);
			VM2/=n;
			VecShift(Veloc_step2,-VM2);
			
			
			
			VecDot(Veloc_step1,Veloc_step1,&sig1);
			VecDot(Veloc_step2,Veloc_step2,&sig2);
			
			VecDot(Veloc_step1,Veloc_step2,&vivi);
			//vivi/=n;
			
			Xi = 1.0 - vivi/PetscSqrtReal(sig1*sig2);
			
			PetscPrintf(PETSC_COMM_WORLD,"------------------------------------------Xi = %lg %d\n",Xi,step);
			
			VecCopy(Veloc_fut,Veloc_step1);
		
		}
	}
	
	PetscFunctionReturn(0);
	
}
