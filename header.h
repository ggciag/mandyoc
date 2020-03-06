

long Nx,Nz;
long layers;

long T_NE = 4;
long T_GN = 1;

long DIMEN = 2;

long GaussQuad = 9;


long V_NE = 4;
long V_GN = 2;
long V_GT = V_NE*V_GN;

double Lx, depth;

PetscScalar *interfaces;
int n_interfaces;

PetscScalar *inter_rho;
PetscScalar *inter_geoq;
PetscScalar *inter_H;

PetscScalar *inter_A;
PetscScalar *inter_n;
PetscScalar *inter_Q;
PetscScalar *inter_V;

PetscInt visc_harmonic_mean;


/////////

int tcont=0;


/////////

double seg_per_ano = 365.0*24.0*3600.0;

double dt_calor = 40000.0;
double dt_calor_sec=dt_calor*seg_per_ano;

double tempo=0;


double alpha_thermal=0.5;
double comp_alpha_thermal = 1.0 - alpha_thermal;



PetscReal rtol;

PetscInt temper_extern;
PetscInt veloc_extern;
PetscInt bcv_extern;

PetscInt visc_const_per_element;

/////////

double dx_const;
double dz_const;

int ContMult;

long stepMAX;
double timeMAX;
double dt_MAX;

long print_step;

double visco_r;

double visc_MAX;
double visc_MIN;


double visc_aux_MAX;
double visc_aux_MIN;

double e2_aux_MAX;
double e2_aux_MIN;

double escala_viscosidade;

double veloc_superf;

double RHOM;
double alpha_exp_thermo;
double kappa;


double gravity;

double Delta_T;

double H_lito;

double h_air;

double H_per_mass;
double c_heat_capacity;

int T_initial_cond;
int rheol;

double beta_max;
double ramp_begin;
double ramp_end;

int bcv_top_normal;
int bcv_top_slip;

int bcv_bot_normal;
int bcv_bot_slip;

int bcv_left_normal;
int bcv_left_slip;

int bcv_right_normal;
int bcv_right_slip;

int bcT_top;

int bcT_bot;

int bcT_left;

int bcT_right;

///////

PetscReal *TKe, *TCe, *TFe, *TCe_fut, *TMe, *Ttotal, *Ttotal_b;


PetscReal *T_vec_aux_ele;
PetscReal *T_vec_aux_ele_final;



double r06 = 0.7745966692414834; //sqrt(0.6)
double r8p9 = 8.0/9.0;
double r5p9 = 5.0/9.0;


PetscReal *NT;
PetscReal *NT_x;
PetscReal *NT_z;


////////

Vec v_vec;
Vec v_vec_fut;

PetscInt *indice_aux_vec_ele;

PetscReal *v_vec_aux_ele;

/////////

Mat TA, TB;
Vec Tf, Temper;

Vec dRho;

DM da_Thermal;

KSP T_ksp;

Vec Temper_Cond;

Vec local_FT;
Vec local_Temper;
Vec local_TC;

int geoq_on;

Vec geoq;
Vec local_geoq;

Vec geoq_rho;
Vec local_geoq_rho;

Vec geoq_H;
Vec local_geoq_H;

Vec geoq_cont;
Vec local_geoq_cont;

Vec geoq_strain;
Vec local_geoq_strain;


/////////

PetscReal *N_x_Gauss;
PetscReal *N_z_Gauss;

/////////

PetscReal denok_min;

PetscInt print_visc;

Mat VA, VB, VG;
Vec Vf, Veloc, Veloc_fut,Veloc_weight,Veloc_0;

Vec Adiag;

Vec Veloc_step1, Veloc_step2;

Vec Vf_P;

Vec Pressure;

int PRESSURE_INIT=0;

DM da_Veloc;

KSP V_ksp;

PetscReal *Vfe;

PetscReal *Ke_veloc;
PetscReal *Ke_veloc_final;

PetscReal *Ke_veloc_general;

PetscReal *VCe;
PetscReal *VfMe;

PetscInt Verif_VG=0;


Vec Precon;

Vec local_Precon;

Vec rk_vec2;

Vec rk_vec;
Vec sk_vec;
Vec gs_vec;
Vec uk_vec;

Vec Veloc_Cond;


Vec zk_vec;
Vec zk_vec2;

Vec local_V;
Vec local_VC;
Vec local_FV;
Vec local_FP;
Vec local_P;

Vec local_dRho;

///////

DM dmcell;

DM dms;

PetscInt particles_per_ele;

PetscInt particles_add_remove;

PetscInt *ppp;
PetscInt *p_remove;
PetscInt *p_i;

PetscReal *p_add_coor;
PetscReal *p_add_r;
PetscInt *p_add_i;
PetscInt *p_add_layer;
PetscReal *p_add_r_strain;



PetscInt cont_particles=0;

PetscInt free_surface_stab;

PetscReal sub_division_time_step;

PetscInt print_step_files;

unsigned int seed;


////// Flags
PetscInt WITH_NON_LINEAR = 0; // Controla o uso da reologia plástica e/ou viscosidade não linear
PetscInt WITH_ADIABATIC_H = 0;     // Controla a adição do calor adiabático
PetscInt WITH_RADIOGENIC_H = 0;    // Controla a adição do calor radiogênico
