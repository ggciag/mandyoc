#include <petscksp.h>

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
	
	a=alpha_exp_thermo*gravity*RHOM; CHKERRQ(ierr);
	
	//VecScale(dRho,a); CHKERRQ(ierr);
	VecAXPBY(dRho,-gravity,a,geoq_rho); CHKERRQ(ierr);
	
	return ierr;
	
}