"""
Create the input file to simulate the an Keken et al. (1997) model
"""
import numpy as np
import matplotlib.pyplot as plt

Nx = 81
Nz = 81

Lx = 92.42
Lz = 100.0

xn = np.linspace(0, Lx, Nx)

liton = -0.8 * Lz + 0.02 * Lz * np.cos(np.pi * xn / Lx)

liton = liton * 1000


# Create the interface file

layer_properties = [
    "C    1.0     1.0",
    "rho -1000.   0.",
    "H    0.0e-12 0.0e-12",
    "A    0.0     0.0",
    "n    0.0     0.0",
    "Q    0.0     0.0",
    "V    0.0     0.0",
]

with open("interfaces.txt", "w") as f:
    for line in layer_properties:
        line = line.strip()
        if len(line):
            f.write(" ".join(line.split()) + "\n")

    for line in np.arange(Nx):
        f.write("%lf\n" % (liton[line]))


##############################################################################
# Parameters file
##############################################################################
params = f"""
# Geometry
nx                                  = {Nx}          # n. of nodes in the longitudinal direction
nz                                  = {Nz}          # n. of nodes in the vertical direction
lx                                  = 91420.0       # extent in the longitudinal direction
lz                                  = 100000.       # extent in the vertical direction
particles_per_element_x             = 0             # default is 0
particles_per_element_z             = 0             # default is 0

# Simulation options
multigrid                           = 1             # ok -> soon to be on the command line only
solver                              = direct        # default is direct [direct/iterative]
denok                               = 1.0e-15       # default is 1.0e-4
particles_per_element               = 1100          # default is 81
particles_perturb_factor            = 0.0           # default is 0.5 [values are between 0 and 1]
rtol                                = 1.0e-5        # the absolute size of the residual norm (relevant only for iterative methods), default is 1.0E-5
RK4                                 = Euler         # default is Euler [Euler/Runge-Kutta]
Xi_min                              = 1.0e-14       # default is 1.0e-14
random_initial_strain               = 0.0           # default is 0.0
pressure_const                      = -1.0          # default is -1.0 (not used)
initial_dynamic_range               = False         # default is False [True/False]
periodic_boundary                   = False         # default is False [True/False]
high_kappa_in_asthenosphere         = False         # default is False [True/False]
K_fluvial                           = 2.0e-7        # default is 2.0e-7
m_fluvial                           = 1.0           # default is 1.0
sea_level                           = 0.0           # default is 0.0
basal_heat                          = -1.0          # default is -1.0

# Surface processes
sp_surface_tracking                 = False         # default is False [True/False]
sp_surface_processes                = False         # default is False [True/False]
sp_dt                               = 0             # default is 0
sp_d_c                              = 0             # default is 0
plot_sediment                       = False         # default is False [True/False]
a2l                                 = True          # default is True [True/False]

free_surface_stab                   = True          # default is True [True/False]
theta_FSSA                          = 0.5           # default is 0.5 (only relevant when free_surface_stab = True)

# Time constrains
step_max                            = 4000          # maximum time-step of the simulation
time_max                            = 140000.0e6    # maximum time of the simulation [years]
dt_max                              = 1000.0e6      # maximum time between steps of the simulation [years]
step_print                          = 20            # make file every <step_print>
sub_division_time_step              = 1.0           # default is 1.0
initial_print_step                  = 0             # default is 0
initial_print_max_time              = 1.0e6         # default is 1.0e6 [years]

# Viscosity
viscosity_reference                 = 1.0e21        # reference viscosity [Pa.s]
viscosity_max                       = 1.0e25        # maximum viscosity [Pa.s]
viscosity_min                       = 1.0e17        # minimum viscosity [Pa.s]
viscosity_per_element               = constant      # default is variable [constant/variable]
viscosity_mean_method               = harmonic      # default is harmonic [harmonic/arithmetic]
viscosity_dependence                = pressure      # default is depth [pressure/depth]

# External ASCII inputs/outputs
interfaces_from_ascii               = True          # default is False [True/False]
n_interfaces                        = 1             # n. of interfaces in the interfaces.txt file
variable_bcv                        = False         # default is False [True/False]
temperature_from_ascii              = False         # default is False [True/False]
velocity_from_ascii                 = False         # default is False [True/False]
binary_output                       = False         # default is False [True/False]
sticky_blanket_air                  = False         # default is False [True/False]
precipitation_profile_from_ascii    = False         # default is False [True/False]
climate_change_from_ascii           = False         # default is False [True/False]


print_step_files                    = True          # default is True [True/False]
checkered                           = False         # print one element in the print_step_filesdefault is False [True/False]

sp_mode                             = 1             # default is 1 [0/1/2]

geoq                                = on            # ok
geoq_fac                            = 1.0           # ok

# Physical parameters
temperature_difference              = 0.            # ok
thermal_expansion_coefficient       = 3.28e-5       # ok
thermal_diffusivity_coefficient     = 1.0e-6        # ok
gravity_acceleration                = 10.0          # ok
density_mantle                      = 3300.         # ok
external_heat                       = 0.0e-12       # ok
heat_capacity                       = 1250.         # ok

non_linear_method                   = off           # ok
adiabatic_component                 = off           # ok
radiogenic_component                = off           # ok

# Velocity boundary conditions
top_normal_velocity                 = fixed         # ok
top_tangential_velocity             = fixed         # ok
bot_normal_velocity                 = fixed         # ok
bot_tangential_velocity             = fixed         # ok
left_normal_velocity                = fixed         # ok
left_tangential_velocity            = free          # ok
right_normal_velocity               = fixed         # ok
right_tangential_velocity           = free          # ok

surface_velocity                    = 0.0e-2        # ok
multi_velocity                      = False         # default is False [True/False]

# Temperature boundary conditions
top_temperature                     = fixed         # ok
bot_temperature                     = fixed         # ok
left_temperature                    = free          # ok
right_temperature                   = free          # ok

rheology_model                      = 0             # ok
T_initial                           = 0             # ok
"""

# Create the parameter file
with open("param.txt", "w") as f:
    for line in params.split("\n"):
        f.write(line + "\n")
