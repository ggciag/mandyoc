"""
This example simulates the evolution of divergent margins, taking into account the plastic rheology and the sin-rift geodynamics
"""

import numpy as np
import matplotlib.pyplot as plt

# total model horizontal extent (m)
Lx = 1600 * 1.0e3
# total model vertical extent (m)
Lz = 300 * 1.0e3
# number of points in horizontal direction
Nx = 1601
# number of points in vertical direction
Nz = 301
# thickness of sticky air layer (m)
H_sa = 40 * 1.0e3
# thickness of lower crust (m)
H_lower_crust = 20 * 1.0e3
# thickness of upper crust (m)
H_upper_crust = 20 * 1.0e3
# total thickness of lithosphere (m)
H_litho = 130 * 1.0e3
# seed depth bellow base of lower crust (m)
seed_depth = 13 * 1.0e3

x = np.linspace(0, Lx, Nx)
z = np.linspace(Lz, 0, Nz)
X, Z = np.meshgrid(x, z)


##############################################################################
# Interfaces (bottom first)
##############################################################################
interfaces = {
    "litho": np.ones(Nx) * (H_litho + H_sa),
    "seed_base": np.ones(Nx) * (seed_depth + H_lower_crust + H_upper_crust + H_sa),
    "seed_top": np.ones(Nx) * (seed_depth + H_lower_crust + H_upper_crust + H_sa),
    "lower_crust": np.ones(Nx) * (H_lower_crust + H_upper_crust + H_sa),
    "upper_crust": np.ones(Nx) * (H_upper_crust + H_sa),
    "air": np.ones(Nx) * (H_sa),
}

# seed thickness (m)
H_seed = 6 * 1.0e3
# seed horizontal position (m)
x_seed = 750 * 1.0e3
# seed: number of points of horizontal extent
n_seed = 2

interfaces["seed_base"][
    int(Nx * x_seed // Lx - n_seed // 2) : int(Nx * x_seed // Lx + n_seed // 2)
] = (
    interfaces["seed_base"][
        int(Nx * x_seed // Lx - n_seed // 2) : int(Nx * x_seed // Lx + n_seed // 2)
    ]
    + H_seed // 2
)
interfaces["seed_top"][
    int(Nx * x_seed // Lx - n_seed // 2) : int(Nx * x_seed // Lx + n_seed // 2)
] = (
    interfaces["seed_top"][
        int(Nx * x_seed // Lx - n_seed // 2) : int(Nx * x_seed // Lx + n_seed // 2)
    ]
    - H_seed // 2
)

Huc = 2.5e-6 / 2700.0  # old 9.2e-10
Hlc = 0.8e-6 / 2800.0  # old 2.9e-10

# Create the interface file
with open("interfaces.txt", "w") as f:
    layer_properties = f"""
        C   1.0       1.0        0.1        1.0        1.0         1.0         1.0
        rho 3378.0    3354.0     3354.0     3354.0     2800.0      2700.0      1.0
        H   0.0       9.0e-12    9.0e-12    9.0e-12    {Hlc}       {Huc}       0.0
        A   1.393e-14 2.4168E-15 2.4168E-15 2.4168E-15 8.574e-28   8.574e-28   1.0e-18
        n   3.0       3.5        3.5        3.5        4.0         4.0         1.0
        Q   429.0e3   540.0E3    540.0E3    540.0E3    222.0e3     222.0e3     0.0
        V   15.0e-6   25.0e-6    25.0e-6    25.0e-6    0.0         0.0         0.0
    """

    for line in layer_properties.split("\n"):
        line = line.strip()
        if len(line):
            f.write(" ".join(line.split()) + "\n")

    # layer interfaces
    data = -1 * np.array(tuple(interfaces.values())).T
    np.savetxt(f, data, fmt="%.1f")

# Plot interfaces
##############################################################################
fig, ax = plt.subplots(figsize=(16, 8))

for label, layer in interfaces.items():
    ax.plot(x, layer, label=f"{label}")

ax.set_xticks(np.arange(0, Lx + 1, 100 * 1.0e3))
ax.set_yticks(np.arange(-Lz, 0 + 1, 10 * 1.0e3))

ax.set_xlim([0, Lx])
ax.set_ylim([-Lz, 0])

plt.legend()

# plt.show()
plt.close()


plt.figure()
for label, layer in interfaces.items():
    print(label, ":", np.size(layer))
    plt.plot(x, layer, label=f"{label}")
ax.set_xlim([0, Lx])
ax.set_ylim([Lz, 0])
plt.legend()
plt.savefig("interfaces_teste.png")
plt.close()


##############################################################################
# Parameters file
##############################################################################
params = f"""
nx = {Nx}
nz = {Nz}
lx = {Lx}
lz = {Lz}


# Simulation options
multigrid                           = 1             # ok -> soon to be on the command line only
solver                              = direct        # default is direct [direct/iterative]
denok                               = 1.0E-15       # default is 1.0E-4
particles_per_element               = 40          # default is 81
particles_perturb_factor            = 0.7           # default is 0.5 [values are between 0 and 1]
rtol                                = 1.0E-7        # the absolute size of the residual norm (relevant only for iterative methods), default is 1.0E-5
RK4                                 = Euler         # default is Euler [Euler/Runge-Kutta]
Xi_min                              = 1.0E-7       # default is 1.0E-14
random_initial_strain               = 0.3           # default is 0.0
pressure_const                      = -1.0          # default is -1.0 (not used) - useful only in horizontal 2D models
initial_dynamic_range               = True         # default is False [True/False]
periodic_boundary                   = False         # default is False [True/False]
high_kappa_in_asthenosphere         = False         # default is False [True/False]
K_fluvial                           = 2.0E-7        # default is 2.0E-7
m_fluvial                           = 1.0           # default is 1.0
sea_level                           = 0.0           # default is 0.0
basal_heat                          = 0.0          # default is -1.0

# Surface processes
sp_surface_tracking                 = False         # default is False [True/False]
sp_surface_processes                = False         # default is False [True/False]
sp_dt                               = 1.0E5        # default is 0.0
sp_d_c                              = 1.0          # default is 0.0
plot_sediment                       = False         # default is False [True/False]
a2l                                 = True          # default is True [True/False]

free_surface_stab                   = True          # default is True [True/False]
theta_FSSA                          = 0.5           # default is 0.5 (only relevant when free_surface_stab = True)

# Time constrains
step_max                            = 5000          # Maximum time-step of the simulation
time_max                            = 100.0E6     # Maximum time of the simulation [years]
dt_max                              = 10.0E3      # Maximum time between steps of the simulation [years]
step_print                          = 10            # Make file every <step_print>
sub_division_time_step              = 0.5           # default is 1.0
initial_print_step                  = 0             # default is 0
initial_print_max_time              = 1.0E6         # default is 1.0E6 [years]

# Viscosity
viscosity_reference                 = 1.0E26        # Reference viscosity [Pa.s]
viscosity_max                       = 1.0E25        # Maximum viscosity [Pa.s]
viscosity_min                       = 1.0E18        # Minimum viscosity [Pa.s]
viscosity_per_element               = constant      # default is variable [constant/variable]
viscosity_mean_method               = arithmetic      # default is harmonic [harmonic/arithmetic]
viscosity_dependence                = pressure      # default is depth [pressure/depth]

# External ASCII inputs/outputs
interfaces_from_ascii               = True          # default is False [True/False]
n_interfaces                        = {len(interfaces.keys())}           # Number of interfaces int the interfaces.txt file
variable_bcv                        = False         # default is False [True/False]
temperature_from_ascii              = True         # default is False [True/False]
velocity_from_ascii                 = True         # default is False [True/False]
binary_output                       = False         # default is False [True/False]
sticky_blanket_air                  = True         # default is False [True/False]
precipitation_profile_from_ascii    = False         # default is False [True/False]
climate_change_from_ascii           = False         # default is False [True/False]


print_step_files                    = True          # default is True [True/False]
checkered                           = False         # Print one element in the print_step_files (default is False [True/False])

sp_mode                             = 5             # default is 1 [0/1/2]

geoq                                = on            # ok
geoq_fac                            = 100.0           # ok

# Physical parameters
temperature_difference              = 1500.         # ok
thermal_expansion_coefficient       = 3.28E-5       # ok
thermal_diffusivity_coefficient     = 1.0E-6        # ok
gravity_acceleration                = 10.0          # ok
density_mantle                      = 3300.         # ok
external_heat                       = 0.0E-12       # ok
heat_capacity                       = 1250.         # ok

non_linear_method                   = on            # ok
adiabatic_component                 = on            # ok
radiogenic_component                = on            # ok

# Velocity boundary conditions
top_normal_velocity                 = fixed         # ok
top_tangential_velocity             = free          # ok
bot_normal_velocity                 = fixed         # ok
bot_tangential_velocity             = free          # ok
left_normal_velocity                = fixed         # ok
left_tangential_velocity            = fixed         # ok
right_normal_velocity               = fixed         # ok
right_tangential_velocity           = fixed         # ok

surface_velocity                    = 0.0E-2        # ok
multi_velocity                      = False         # default is False [True/False]

# Temperature boundary conditions
top_temperature                     = fixed         # ok
bot_temperature                     = fixed         # ok
left_temperature                    = fixed          # ok
right_temperature                   = fixed          # ok

rheology_model                      = 9             # ok
T_initial                           = 3             # ok

"""
# Create the parameter file
with open("param.txt", "w") as f:
    for line in params.split("\n"):
        line = line.strip()
        if len(line):
            f.write(" ".join(line.split()) + "\n")

##############################################################################
# Initial temperature field
##############################################################################

T = 1300 * (z - H_sa) / (H_litho)  # Temperatura

Ta = 1262 / np.exp(-10 * 3.28e-5 * (z - H_sa) / 1250)

T[T < 0.0] = 0.0
T[T > Ta] = Ta[T > Ta]

kappa = 1.0e-6

ccapacity = 1250

H = np.zeros_like(T)

cond = (z >= H_sa) & (z < H_upper_crust + H_sa)  # upper crust
H[cond] = Huc

cond = (z >= H_upper_crust + H_sa) & (
    z < H_lower_crust + H_upper_crust + H_sa
)  # lower crust
H[cond] = Hlc

Taux = np.copy(T)
t = 0
dt = 10000
dt_sec = dt * 365 * 24 * 3600
cond = (z > H_sa + H_litho) | (T == 0)  # (T > 1300) | (T == 0)
dz = Lz / (Nz - 1)

while t < 500.0e6:
    T[1:-1] += (
        kappa * dt_sec * ((T[2:] + T[:-2] - 2 * T[1:-1]) / dz ** 2)
        + H[1:-1] * dt_sec / ccapacity
    )
    T[cond] = Taux[cond]
    t = t + dt

T = np.ones_like(X) * T[:, None]

print(np.shape(T))

# Save the initial temperature file
np.savetxt("input_temperature_0.txt", np.reshape(T, (Nx * Nz)), header="T1\nT2\nT3\nT4")

# Plot temperature field and thermal profile
##############################################################################

fig, (ax0, ax1) = plt.subplots(nrows=1, ncols=2, figsize=(16, 8))

ax0.contour(X / 1.0e3, (Z - H_sa) / 1.0e3, T, levels=np.arange(0, 1610, 100))
ax0.set_ylim((Lz - H_sa) / 1.0e3, -H_sa / 1000)
ax0.set_xlabel("km")
ax0.set_ylabel("km")
ax1.set_xlabel("$^\circ$C")

ax1.plot(T[:, 0], (z - H_sa) / 1.0e3, "-k")

code = 0
for label in list(interfaces.keys()):
    code += 1
    color = "C" + str(code)
    ax1.hlines(
        (interfaces[label][0] - H_sa) / 1.0e3,
        np.min(T[:, 0]),
        np.max(T[:, 0]),
        label=f"{label}",
        color=color,
    )

ax1.set_ylim((Lz - H_sa) / 1.0e3, -H_sa / 1000)

plt.legend()
plt.savefig("temperature_field.png")
# plt.show()
plt.close()


##############################################################################
# Boundary condition - velocity
##############################################################################

fac_air = 10.0e3

# 1 cm/ano
vL = 0.005 / (365 * 24 * 3600)  # m/s

h_v_const = 150.0e3  # espessura com velocidade constante
ha = Lz - H_sa - h_v_const  # diferenca

vR = 2 * vL * (h_v_const + ha) / ha  # garante que integral seja zero

VX = np.zeros_like(X)
cond = (Z > h_v_const + H_sa) & (X == 0)
VX[cond] = vR * (Z[cond] - h_v_const - H_sa) / ha

cond = (Z > h_v_const + H_sa) & (X == Lx)
VX[cond] = -vR * (Z[cond] - h_v_const - H_sa) / ha

cond = X == Lx
VX[cond] += +2 * vL

cond = Z <= H_sa - fac_air
VX[cond] = 0

# print(np.sum(VX))

v0 = VX[(X == 0)]
vf = VX[(X == Lx)]
sv0 = np.sum(v0[1:-1]) + (v0[0] + v0[-1]) / 2.0
svf = np.sum(vf[1:-1]) + (vf[0] + vf[-1]) / 2.0
# print(sv0, svf, svf - sv0)

diff = (svf - sv0) * dz

vv = -diff / Lx
# print(vv, diff, svf, sv0, dz, Lx)

VZ = np.zeros_like(X)

cond = Z == 0
VZ[cond] = vv

# print(np.sum(v0))
# print(np.sum(vf))
# print(np.sum(vf)-np.sum(v0))

VVX = np.copy(np.reshape(VX, Nx * Nz))
VVZ = np.copy(np.reshape(VZ, Nx * Nz))

v = np.zeros((2, Nx * Nz))

v[0, :] = VVX
v[1, :] = VVZ

v = np.reshape(v.T, (np.size(v)))

# Create the initial velocity file
np.savetxt("input_velocity_0.txt", v, header="v1\nv2\nv3\nv4")

# Plot veolocity
##############################################################################
fig, (ax0, ax1) = plt.subplots(nrows=1, ncols=2, figsize=(9, 9))

ax0.plot(VX[:, 0], (z - H_sa) / 1000, "k-", label="left side")
ax1.plot(VZ[:, 0], (z - H_sa) / 1000, "k-", label="left side")

ax0.plot(VX[:, -1], (z - H_sa) / 1000, "r-", label="right side")
ax1.plot(VZ[:, -1], (z - H_sa) / 1000, "r-", label="right side")

ax0.legend()
ax1.legend()

ax0_xlim = ax0.get_xlim()
ax1_xlim = ax1.get_xlim()

ax0.set_yticks(np.arange(0, Lz / 1000, 20))
ax1.set_yticks(np.arange(0, Lz / 1000, 20))

ax0.set_ylim([Lz / 1000 - H_sa / 1000, -H_sa / 1000])
ax1.set_ylim([Lz / 1000 - H_sa / 1000, -H_sa / 1000])

ax0.set_xlim([-8e-10, 8e-10])
ax1.set_xlim([-8e-10, 8e-10])

ax0.set_xlabel("(m/s)")
ax0.set_ylabel("Depth (km)")

ax0.set_title("horizontal component of velocity")

ax1.set_title("vertical component of velocity")

plt.savefig("velocity.png")
plt.close()
