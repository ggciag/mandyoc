#include <petscksp.h>
#include <stdio.h>
#include <petscmath.h>

extern int rheol;

extern double visc_MAX;
extern double visc_MIN;

extern int geoq_on;

extern double visco_r;

extern double Delta_T;

extern PetscInt WITH_NON_LINEAR;


double strain_softening(double strain, double f1, double f2){
	double fac;
	
	double st1=0.05,st2=1.05;

	if (strain<st1) fac=f1;
	else{
		if (strain>st2) fac=f2;
		else{
			fac = f1 - (f1-f2)*(strain-st1)/(st2-st1);
		}
	}
	
	return(fac);
}



double calc_visco_ponto(double T,double z,double geoq_ponto,double e2_inva,double strain_cumulate,
						double A, double n_exp, double QE, double VE){
	
	double visco_real;
	
	if (e2_inva<1.0E-18) e2_inva=1.0E-18; ///!!!! e2_inva min
	
	
	if (rheol==0)	visco_real = visco_r;
	
	if (rheol==1){
		double r=20.0;
		double Q = 225.0/log(r)-0.25*log(r);
		double G = 15./log(r)-0.5;
		
		return(geoq_ponto*visco_r*exp(  Q/(T/Delta_T+G) - Q/(0.5+G)    ));
	}
	
	if (rheol==2){
		double R = 8.31;     // J/mol/K
		double E = 120000.0; // J/mol
		
		
		double Tb = Delta_T+273.0;
		
		double aux = E*(1.0/(T+273.0)-1.0/Tb)/R;
		visco_real = visco_r*exp(aux);
		
	}
	
	
	if (rheol==3){
		double R = 8.31;     // J/mol/K
		double E = 120000.0; // J/mol
		
		
		double Tb = 1300.0+273.0;
		
		double aux = E*(1.0/(T+273.0)-1.0/Tb)/R;
		visco_real = visco_r*exp(aux);
		
	}
	
	if (rheol==4){
		double R = 8.31;     // J/mol/K
		double E = 120000.0; // J/mol
		
		
		double Tb = Delta_T+273.0;
		
		double aux = E*(1.0/(T+273.0)-1.0/Tb)/R;
		visco_real = visco_r*exp(aux);
		
	}
	
	if (rheol==5){
		double R = 8.31;     // J/mol/K
		double E = 240000.0; // J/mol
		
		double b = 1.0E7;
		
		
		double Tb = Delta_T+273.0;
		
		double aux = -(T+273)*E/(R*Tb*Tb);
		visco_real = visco_r*b*exp(aux);
		
	}
	
	if (rheol==6){
		double R = 8.3144;     // J/mol/K
		double E = 240000.0; // J/mol
		
		double b = 1.0E7;
		
		
		double Tb = Delta_T+273.0;
		
		double aux = -(T+273)*E/(R*Tb*Tb);
		visco_real = visco_r*b*exp(aux);
		
		
	}
	
	if (rheol==7){
		double R = 8.3144;
		double TK = T+273.0;
		
		//double aux = -(T+273);
		visco_real = visco_r*A*exp(-(QE+VE*10.0*3300.*(-z))/(R*TK));
	}
	
	if (rheol==8){
		
		double TK = T+273.0;
		
		visco_real = visco_r*exp(-QE*TK + VE*(-z));
	}
	
	if (WITH_NON_LINEAR==1){
		if (rheol==9){
			double R = 8.3144;
			
			double TK = T+273.;
			
			
			
			visco_real = pow(A,-1./n_exp)*pow(e2_inva,(1.-n_exp)/(2*n_exp))*exp((QE+VE*10.0*3300.*(-z))/(n_exp*R*TK));
			//printf("%e %e %.1f %e %e %e\n",A,e2_inva,n_exp,QE,VE,visco_real);
		}
	}
	
	
	if (rheol>9){
		printf("rheol error: larger than maximum available option\n");
		exit(1);
	}
	
	if (geoq_on)
		visco_real *= geoq_ponto;
	
	
	///!!!!
	if (WITH_NON_LINEAR==1){
		//double c0 = 1.0;//!!!!
		//double mu = 0.01;//!!!!
		//double c0 = 22.0E6;//!!!! Petersen et al. (2010)
		//double mu = 0.58778;//!!!!
		
		double c0 = strain_softening(strain_cumulate,20.0E6,4.0E6);
		double mu = strain_softening(strain_cumulate,0.261799,0.034906);
		
		double tau_yield = c0*cos(mu) + sin(mu)*10.0*3300.*(-z);//!!!!
		
		double visco_yield = visc_MAX;
		
		if (e2_inva>0) visco_yield = tau_yield/e2_inva;
		
		if (visco_real>visco_yield) visco_real = visco_yield;
		
		/*visco_real = visco_r;
		if (2*visco_real*e2_inva-1.0>0){
			visco_real = 1.0/(2*(e2_inva));
		}*/
		//if (e2_inva>0)	visco_real = 1.0/(2*(e2_inva));//!!!! rigid plastic
	}
	
	
	if (visco_real>visc_MAX) visco_real=visc_MAX;
	if (visco_real<visc_MIN) visco_real=visc_MIN;
	
	return(visco_real);
	
	return(0);
	
	
	
}
