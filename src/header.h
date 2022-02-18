#ifndef GIT_VERSION
#define GIT_VERSION "git-version-unavailable"
#endif

// Parameter file numerical variables
PetscReal rtol = PETSC_DEFAULT;
PetscReal denok_min = 1.0E-4;
PetscInt particles_per_ele = 81;
PetscReal theta_FSSA = 0.5;
PetscReal sub_division_time_step = 1.0;
PetscReal particles_perturb_factor = 0.5;
PetscInt sp_mode = 1;
PetscReal Xi_min = 1.0E-14;
PetscReal random_initial_strain = 0;
PetscReal pressure_const = -1.0;
PetscInt nx_ppe = 0;
PetscInt nz_ppe = 0;
PetscInt initial_print_step = 0;
PetscReal initial_print_max_time = 1.0E6;
PetscScalar K_fluvial = 2.0E-7;
PetscScalar m_fluvial = 1.0;
PetscScalar sea_level = 0.0;
PetscScalar basal_heat = -1.0;
PetscReal sp_dt = 0.0;
PetscScalar sp_d_c = 0.0; 
// Parameter file boolean variables
PetscInt WITH_NON_LINEAR = 0; // 1=True, 0=False
PetscInt WITH_ADIABATIC_H = 0; // 1=True, 0=False
PetscInt WITH_RADIOGENIC_H = 0; // 1=True, 0=False
PetscInt direct_solver = 1; // 1=direct, 0=iterative
PetscInt visc_const_per_element = 0; // 1=constant, 0=variable
PetscInt visc_harmonic_mean = 1; // 1=harmoninc, 0=arithmetic
PetscInt pressure_in_rheol = 0; // 1=pressure, 0=depth
PetscInt variable_bcv = 0; // 0=False, 1=True
PetscInt interfaces_from_ascii = 0; // 1=True, 0=False
PetscInt temper_extern = 0; // 1=True, 0=False
PetscInt veloc_extern = 0; // 1=True, 0=False
PetscInt bcv_extern = 0; // 1=True, 0=False
PetscInt binary_output = 0; // 1=True, 0=False
PetscInt sticky_blanket_air = 0; // 1=True, 0=False
PetscInt multi_velocity = 0; // 1=True, 0=False
PetscInt precipitation_profile = 0; // 1=True, 0=False
PetscInt climate_change = 0; // 1=True, 0=False
PetscInt free_surface_stab = 1; // 1=True, 0=False
PetscInt print_step_files = 1; // 1=True, 0=False
PetscInt RK4 = 0; // 0=Euler, 1=Runge-Kutta (not working yet!)
PetscInt checkered = 0; // 1=True, 0=False
PetscInt initial_dynamic_range = 0; // 1=True, 0=False
PetscInt periodic_boundary = 0; // 1=True, 0=False
PetscInt high_kappa_in_asthenosphere = 0; // 1=True, 0=False
// Will be added to param.txt
PetscBool sp_surface_tracking = PETSC_FALSE; // PETSC_TRUE/PETSC_FALSE
PetscBool sp_surface_processes = PETSC_FALSE; // PETSC_TRUE/PETSC_FALSE
PetscBool set_sp_dt = PETSC_FALSE; // PETSC_TRUE/PETSC_FALSE
PetscBool set_sp_d_c = PETSC_FALSE; // PETSC_TRUE/PETSC_FALSE
PetscBool plot_sediment = PETSC_FALSE; // PETSC_TRUE/PETSC_FALSE
PetscBool a2l = PETSC_TRUE; // PETSC_TRUE/PETSC_FALSE
// Parameter file native C variables
long Nx = -1;
long Nz = -1;
long layers;
double Lx, depth;
int n_interfaces = 0;
int ContMult;
long stepMAX;
double timeMAX;
double dt_MAX;
long print_step;
double visco_r;
double visc_MAX;
double visc_MIN;
int geoq_on = 1; // 1=on, 2=off
double escala_viscosidade;
double veloc_superf;
double RHOM;
double alpha_exp_thermo;
double kappa;
double gravity;
double Delta_T;
double H_per_mass;
double c_heat_capacity;
int T_initial_cond;
int rheol;
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
// End of parameter file variables

double visc_MAX_comp;
double visc_MIN_comp;
double visc_aux_MAX;
double visc_aux_MIN;

long T_NE = 4;
long T_GN = 1;

long DIMEN = 2;

long GaussQuad = 9;

long V_NE = 4;
long V_GN = 2;
long V_GT = V_NE*V_GN;

// Interfaces file variables
PetscScalar *interfaces;
PetscScalar *inter_rho;
PetscScalar *inter_geoq;
PetscScalar *inter_H;
PetscScalar *inter_A;
PetscScalar *inter_n;
PetscScalar *inter_Q;
PetscScalar *inter_V;

int tcont=0;

double seg_per_ano = 365.0*24.0*3600.0;

double dt_calor = 0.0;
double dt_calor_sec=dt_calor*seg_per_ano;

double tempo=0;


double alpha_thermal=0.5;
double comp_alpha_thermal = 1.0 - alpha_thermal;

PetscInt i_veloc=0;

double dx_const;
double dz_const;


double e2_aux_MAX;
double e2_aux_MIN;

double H_lito;

double h_air;

double beta_max;
double ramp_begin;
double ramp_end;

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
Vec Tf, Temper, Temper_0;

Vec dRho;

DM da_Thermal;

KSP T_ksp;

Vec Temper_Cond;

Vec local_FT;
Vec local_Temper;
Vec local_TC;

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


Vec geoq_strain_rate;
Vec local_geoq_strain_rate;

Mat VA, VB, VG;
Vec Vf, Veloc, Veloc_fut,Veloc_weight,Veloc_0;

Vec Adiag;

Vec Veloc_step1, Veloc_step2;

Vec Vf_P;

Vec Pressure;
Vec Pressure_aux;

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

PetscInt Verif_first_veloc=0;


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
Vec local_P_aux;

Vec local_dRho;

///////

DM dmcell;

DM dms;

PetscInt particles_add_remove;

PetscInt *ppp;
PetscInt *p_remove;
PetscInt *p_i;

PetscReal *p_add_coor;
PetscReal *p_add_r;
PetscInt *p_add_i;
PetscInt *p_add_layer;
PetscReal *p_add_r_strain;
PetscReal *p_add_r_strain_rate;

PetscInt cont_particles=0;

unsigned int seed;

PetscInt *seed_layer;
PetscInt seed_layer_size;
PetscBool seed_layer_set = PETSC_FALSE;

PetscReal *strain_seed_layer;
PetscInt strain_seed_layer_size;
PetscBool strain_seed_layer_set = PETSC_FALSE;

PetscScalar *var_bcv_time;
PetscScalar *var_bcv_scale;
PetscInt n_var_bcv=0;
PetscInt cont_var_bcv=0;

PetscScalar *mv_time;
PetscInt n_mv=0;
PetscInt cont_mv=0;

////// Flags

PetscInt print_step_aux;

Vec sp_surface_global;
Vec sp_surface_global_n;
Vec sp_surface_coords_global;
Vec sp_top_surface_global;
Vec sp_bot_surface_global;
Vec sp_mean_surface_global;
Vec sp_top_surface_global_n;
Vec sp_bot_surface_global_n;
PetscReal sp_eval_time;
PetscReal sp_last_eval_time;

long sp_n_profiles;
PetscScalar *topo_var_time;
PetscScalar *topo_var_rate;

Vec sp_surface_global_aux;

PetscScalar *global_surface_array_helper;
PetscScalar *global_surface_array_helper_aux;

PetscScalar prec_factor=1.0;

PetscScalar *var_climate_time;
PetscScalar *var_climate_scale;
PetscInt n_var_climate;
PetscInt cont_var_climate=0;
