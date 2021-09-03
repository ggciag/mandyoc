import numpy as np
import matplotlib.pyplot as plt


# total model horizontal extent (m)
Lx = 2800 * 1.0e3

# total model vertical extent (m)
Lz = 850 * 1.0e3

# number of points in horizontal direction
Nx = 281

# number of points in vertical direction
Nz = 169

# thickness of sticky air layer (m)
H_sa = 150 * 1.0e3

# thickness of lid (m)
H_lid = 100 * 1.0e3

# plume diameter (m)
D_plume = 100 * 1.0e3

# plume center (above the bottom) (m)
H_plume = 300 * 1.0e3

x = np.linspace(0, Lx, Nx)
z = np.linspace(Lz, 0, Nz)
X, Z = np.meshgrid(x, z)

##############################################################################
# Interfaces (bottom first)
##############################################################################

x_plume_left = int(((Lx - D_plume) / 2) / (Lx / (Nx - 1)))
x_plume_right = int(((Lx + D_plume) / 2) / (Lx / (Nx - 1)))

plume_bottom = -1 * np.ones(Nx) * (Lz - H_plume)
plume_top = -1 * np.ones(Nx) * (Lz - H_plume)

plume_bottom[x_plume_left:x_plume_right] = -1 * (Lz - H_plume) - np.sqrt(
    np.power(D_plume / 2, 2) - np.power(x[x_plume_left:x_plume_right] - Lx / 2, 2)
)
plume_top[x_plume_left:x_plume_right] = -1 * (Lz - H_plume) + np.sqrt(
    np.power(D_plume / 2, 2) - np.power(x[x_plume_left:x_plume_right] - Lx / 2, 2)
)

interfaces = {
    "plume_bottom": plume_bottom,
    "plume_top": plume_top,
    "mantle_top": -1 * (np.ones(Nx) * (H_sa + H_lid)),
    "lid_top": -1 * (np.ones(Nx) * (H_sa)),
}

with open("interfaces.txt", "w") as f:
    layer_properties = f"""
        C   1.0    0.1    1.0    100.0  0.01
        rho 3300.0 3200.0 3300.0 3300.0 0.0
        H   0.0    0.0    0.0    0.0    0.0
        A   0.0    0.0    0.0    0.0    0.0
        n   0.0    0.0    0.0    0.0    0.0
        Q   0.0    0.0    0.0    0.0    0.0
        V   0.0    0.0    0.0    0.0    0.0
    """

    for line in layer_properties.split("\n"):
        line = line.strip()
        if len(line):
            f.write(" ".join(line.split()) + "\n")

    # layer interfaces
    data = np.array(tuple(interfaces.values())).T
    np.savetxt(f, data, fmt="%f")

##############################################################################
# Plot interfaces
##############################################################################
fig, ax = plt.subplots(figsize=(16, 8))

for label, layer in interfaces.items():
    ax.plot(x / 1.0e3, layer / 1.0e3, label=f"{label}")

ax.set_xticks(np.arange(0, Lx / 1.0e3 + 1, 100))
ax.set_yticks(np.arange(-Lz / 1.0e3, 0 + 1, 50))

ax.set_xlim([0, Lx / 1.0e3])
ax.set_ylim([-Lz / 1.0e3, 0])

plt.legend()

plt.savefig("interfaces.png")
plt.close()

##############################################################################
# Parameters file
##############################################################################
params = f"""# Geometry
nx                                  = {Nx}           # Number of elements in the longitudinal direction
nz                                  = {Nz}           # Number of elements in the vertical direction
lx                                  = {Lx:.1e}       # Extent in the longitudinal direction
lz                                  = {Lz:.1e}       # Extent in the vertical direction

# Simulation options
multigrid                           = 1
solver                              = direct        # default is direct [direct/iterative]
denok                               = 1.0E-15       # default is 1.0E-4
particles_per_element               = 600           # default is 81
particles_per_element_x             = 6             # Number of particles per element in longitudinal (default is 0)
particles_per_element_z             = 100           # Number of particles per element in vertical (default is 0)
particles_perturb_factor            = 0.7           # default is 0.5 [values are between 0 and 1]
rtol                                = 1.0E-7        # the absolute size of the residual norm (relevant only for iterative methods), default is 1.0E-5
RK4                                 = Euler         # default is Euler [Euler/Runge-Kutta]
Xi_min                              = 1.0E-14       # default is 1.0E-14
random_initial_strain               = 0.0           # default is 0.0
pressure_const                      = -1.0          # default is -1.0 (not used) - useful only in horizontal 2D models
initial_dynamic_range               = False         # default is False [True/False]
periodic_boundary                   = False         # default is False [True/False]
high_kappa_in_asthenosphere         = False         # default is False [True/False]
K_fluvial                           = 2.0E-7        # default is 2.0E-7
m_fluvial                           = 1.0           # default is 1.0
sea_level                           = 0.0           # default is 0.0
basal_heat                          = -1.0          # default is -1.0

# Surface processes
sp_surface_tracking                 = True          # default is False [True/False]
sp_surface_processes                = False         # default is False [True/False]
sp_mode                             = 0             # default is 1 [0/1/2]
sp_dt                               = 0             # default is 0.0
sp_d_c                              = 0.0           # default is 0.0
plot_sediment                       = False         # default is False [True/False]
a2l                                 = True          # default is True [True/False]

free_surface_stab                   = True          # default is True [True/False]
theta_FSSA                          = 0.5           # default is 0.5 (only relevant when free_surface_stab = True)

# Time constrains
step_max                            = 7000          # Maximum time-step of the simulation
time_max                            = 21.0E6        # Maximum time of the simulation [years]
dt_max                              = 10.0E3        # Maximum time between steps of the simulation [years]
step_print                          = 10            # Make file every <step_print>
sub_division_time_step              = 1.0           # default is 1.0
initial_print_step                  = 2             # default is 0
initial_print_max_time              = 1.0E6         # default is 1.0E6 [years]

# Viscosity
viscosity_reference                 = 1.0E21        # Reference viscosity [Pa.s]
viscosity_max                       = 1.0E23        # Maximum viscosity [Pa.s]
viscosity_min                       = 1.0E19        # Minimum viscosity [Pa.s]
viscosity_per_element               = constant      # default is variable [constant/variable]
viscosity_mean_method               = arithmetic    # default is harmonic [harmonic/arithmetic]
viscosity_dependence                = pressure      # default is depth [pressure/depth]

# External ASCII inputs/outputs
interfaces_from_ascii               = True          # default is False [True/False]
n_interfaces                        = {len(interfaces.keys())}             # Number of interfaces in the 'interfaces.txt' file
variable_bcv                        = False         # default is False [True/False]
temperature_from_ascii              = False         # default is False [True/False]
velocity_from_ascii                 = False         # default is False [True/False]
binary_output                       = False         # default is False [True/False]
sticky_blanket_air                  = False         # default is False [True/False]
precipitation_profile_from_ascii    = False         # default is False [True/False]
climate_change_from_ascii           = False         # default is False [True/False]

print_step_files                    = True          # default is True [True/False]
checkered                           = False         # Print one element in the print_step_files (default is False [True/False])

geoq                                = on
geoq_fac                            = 100.0

# Physical parameters
temperature_difference              = 0.0
thermal_expansion_coefficient       = 3.28E-5
thermal_diffusivity_coefficient     = 1.0E-6
gravity_acceleration                = 10.0
density_mantle                      = 3300.0
external_heat                       = 0.0E-12
heat_capacity                       = 1250.0

non_linear_method                   = off
adiabatic_component                 = off
radiogenic_component                = on

# Velocity boundary conditions
top_normal_velocity                 = fixed
top_tangential_velocity             = free
bot_normal_velocity                 = fixed
bot_tangential_velocity             = fixed
left_normal_velocity                = fixed
left_tangential_velocity            = free
right_normal_velocity               = fixed
right_tangential_velocity           = free

surface_velocity                    = 0.0
multi_velocity                      = False         # default is False [True/False]

# Temperature boundary conditions
top_temperature                     = fixed
bot_temperature                     = fixed
left_temperature                    = fixed
right_temperature                   = fixed

rheology_model                      = 0
T_initial                           = 0"""

with open("param.txt", "w") as f:
    for line in params.split("\n"):
        f.write(line + "\n")
