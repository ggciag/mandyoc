#include <petscksp.h>
#include <stdio.h>
#include <petscmath.h>

extern int rheol;
extern double visc_MAX;
extern double visc_MIN;
extern double visc_MAX_comp;
extern double visc_MIN_comp;
extern int geoq_on;
extern double visco_r;
extern double Delta_T;
extern PetscInt WITH_NON_LINEAR;
extern PetscInt pressure_in_rheol;
extern PetscReal pressure_const;
extern double h_air;
extern int tcont;

extern PetscReal visc0_scaled;
extern PetscReal h0_scaled;
extern PetscReal strain_rate0_scaled;
extern PetscReal pressure0_scaled;

double strain_softening(double strain, double f1, double f2)
{
	double fac;
	double st1 = 0.05;
	double st2 = 1.05;

	if (strain<st1) {fac = f1;}
	else
	{
		if (strain>st2) {fac = f2;}
		else
		{
			fac = f1 - (f1-f2) * (strain-st1) / (st2-st1);
		}
	}
	
	return(fac);
}

double calc_visco_ponto(double T,double P, double x, double z,double geoq_ponto,double e2_inva,double strain_cumulate,
						double A, double n_exp, double QE, double VE){

	e2_inva *= strain_rate0_scaled;
	P *= pressure0_scaled;
	z *= h0_scaled;

	
	double visco_real = visc_MIN;
	double depth = 0.0;
	if (pressure_in_rheol==0){
		depth = -(z + h_air);
		if (depth<0.0) depth=0.0;
	}
	if (P<0.0) {P = 0.0;}

	if (pressure_const>=0.0) P = pressure_const;

	
	if (e2_inva<1.0E-36) e2_inva=1.0E-36;
	
	
	if (rheol==0)	visco_real = visco_r;
	
	if (rheol==1){
		double r = 20.0;
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
		
		if (pressure_in_rheol==0) {
			visco_real = visco_r*A*exp(-(QE+VE*10.0*3300.*(depth))/(R*TK));
		}
		else {
			visco_real = visco_r*A*exp(-(QE+VE*P)/(R*TK));
		}
	}
	
	if (rheol==8){
		
		double TK = T+273.0;
		
		visco_real = visco_r*exp(-QE*TK + VE*(-z));
	}
	
	if (WITH_NON_LINEAR==1){
		if (rheol==9){
			double R = 8.3144;
			
			double TK = T+273.;
			
			if (pressure_in_rheol==0) {
				visco_real = pow(A,-1./n_exp)*pow(e2_inva,(1.-n_exp)/(n_exp))*exp((QE+VE*10.0*3300.*(depth))/(n_exp*R*TK));
			}
			else {
				visco_real = pow(A,-1./n_exp)*pow(e2_inva,(1.-n_exp)/(n_exp))*exp((QE+VE*P)/(n_exp*R*TK));
			}
		}
		if (rheol==19){
			double R = 8.3144;
			
			double TK = T+273.;

			visco_real = pow(A,-1./n_exp)*pow(e2_inva,(1.-n_exp)/(n_exp))*exp((QE+VE*P)/(n_exp*R*TK));

			if (n_exp>2.99 && n_exp<3.01){
				//diffusion creep
				//A = Adiff/mu = 5.3E15/8.0E10 = 66250.0;
				double A1 = 66250.0;
				//double m1 = -2.5;
				//hb = (h/b)^(-m)
				double hb_power_m = 5656854249492380.0;
				//double b1 = 5.0E-10;
				//double mu1 = 8.0E10;
				//double h1 = 0.001;

				double QE1 = 240000.0;
				double VE1 = 5.0E-6;

				double visc_diffu = (hb_power_m*exp((QE1+VE1*P)/(R*TK)))/(A1);

				if (visco_real>visc_diffu) visco_real = visc_diffu;

			}

		}
	}
	

	if (rheol==10){
		double beta = 6.907755279;
		double DT = 1000.0;

		visco_real = visco_r * exp(-(beta*T/DT));
	}
	
	if (geoq_on)
		visco_real *= geoq_ponto;
	
	
	
	if (WITH_NON_LINEAR==1){
		//double c0 = 1.0;// Petersen et al. (2010) plastic criterium
		//double mu = 0.01;//
		//double c0 = 22.0E6;// 
		//double mu = 0.58778;//
		
		double c0 = strain_softening(strain_cumulate,20.0E6,4.0E6);
		double mu = strain_softening(strain_cumulate,0.261799,0.034906);
		double tau_yield;
		if (pressure_in_rheol==0){
			tau_yield = c0*cos(mu) + sin(mu)*10.0*3300.*(depth);
		}
		else {
			tau_yield = c0*cos(mu) + sin(mu)*P;
		}
		
		double visco_yield = visc_MAX;
		
		if (e2_inva>0) visco_yield = tau_yield/(2*e2_inva);
		
		if (visco_real>visco_yield) visco_real = visco_yield;


		if (rheol == 70){
			visco_real = visco_r;
			if (2*visco_real*e2_inva-1.0>0){
				visco_real = 1.0/(2*(e2_inva));
			}
		}
		
	}


	//printf("%lf %lg %lg %lg\n",z,P,e2_inva,visco_real);
	visco_real /=visc0_scaled;
	
	
	if (visco_real>visc_MAX) visco_real=visc_MAX;
	if (visco_real<visc_MIN) visco_real=visc_MIN;

	
	double f1 = PetscLogReal(visc_MAX_comp/visc_MIN_comp);
	double f2 = PetscLogReal(visc_MAX/visc_MIN);
	double f3 = PetscLogReal(visco_real/visc_MIN);
	visco_real = visc_MIN_comp*PetscExpReal(f1*f3/f2);

	
	return(visco_real);
	
	return(0);
}
