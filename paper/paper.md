---
title: 'Mandyoc: A finite element code to simulate thermo-chemical convection in parallel'
tags:
  - PETSc
  - mantle convection
  - lithosphere geodynamics
  - finite element

authors:
  - name: Victor Sacek
    orcid: 0000-0003-0872-7098
    affiliation: 1
  - name: Jamison F. Assunção
    affiliation: 1
  - name: Agustina Pesce
    affiliation: 2
  - name: Rafael Monteiro da Silva
    affiliation: 1
affiliations:
 - name: Instituto de Astronomia, Geofísica e Ciências Atmosféricas, Universidade de São Paulo, Brazil
   index: 1
 - name: Instituto Geofísico Sismológico Ing. Volponi, Universidad Nacional de San Juan, CONICET, Argentina
   index: 2
date: 10 June 2021
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
#aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
#aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

`Mandyoc` is a finite element code written in C dedicated to simulate thermo-chemical convection in the interior of terrestrial planets. Different linear and non-linear rheologies can be adopted, appropriately simulating the strain and stress pattern in the Earth's crust and mantle, both in extensional or collisional tectonics. 
Additionally, the code allows variations of boundary condition for the velocity field in space and time, simulating different pulses of tectonism in the same numerical scenario. 


# Statement of need

Although `Mandyoc` is the acronym for MANtle DYnamics simulatOr Code, the code is not limited to mantle convection, but it is designed to simulate general thermo-chemical convection in Stokes flow, taking different compositional layers into account, also appropriate to simulate Earth's lithospheric dynamics in the geological timescale.

Previous versions of the code was applied to study the evolution of continental margins, showing the interaction between continental lithosphere with the asthenospheric mantle [@sacek2017post;@salazar2021lateral].


# Mathematics


`Mandyoc` solves the equations for conservation of mass, momentum and energy using the Finite Element Method assuming the extended Boussinesq approximation, respectively:

$$u_{i,i}=0$$

$$ \sigma_{ij,j} + g_i \rho = 0 $$

$$ \frac{\partial T}{\partial t} + u_i T_{,i} = \kappa T_{,ii} + \frac{H}{c_p\rho} + u_i g_i \alpha T / c_p$$

where

$$\sigma_{ij} = -P \delta_{ij} + \eta \left( u_{i,j} + u_{j,i}\right)$$

$$\rho = \rho_0 \left( 1 - \alpha (T-T_0)\right)$$

$u_i$ is the component $i$ of the velocity field,
$T$ is temperature, $t$ is time, $\kappa$ is the thermal diffusivity, $H$ is the volumetric heat production, $c_p$ is the specific heat capacity, $g$ is gravity, $\rho$ is the effective rock density dependent on temperature and composition, $\rho_0$ is the reference rock density at temperature $T_0$, $\alpha$ is the coefficient of thermal expansion, $P$ is the dynamic pressure, $\eta$ is the effective viscosity, and $\delta_{ij}$ is the Kronecker delta.

The code is fully parallelized using the Portable, Extensible Toolkit for Scientific Computation (PETSc) [@petsc-efficient;@petsc-user-ref;@petsc-web-page].
The present version of the code can simulate thermochemical convection
using different rheological formulations: Newtonian flow, non-linear viscous flow or visco-plastic deformation.
For example, the lithosphere can be simulated as a combination of different visco-plastic layers in which the
effective viscosity depends on a nonlinear power law viscous rheology and a plastic yield criterium, like the
Drucker-Prager criterium. Additionally, strain softening is implemented to facilitate the localization of strain in
the plastic regime during, for example, lithospheric stretching.

The composition and strain history is tracked by particles present in the interior of the finite element. The exchange of particles among the subdomains of the model is efficiently parallelized in PETSc using DMSwarm [@may2017dmswarm].

The free surface of the Earth can be simulated and is numerically stabilized using the Free Surface Stabilization
Algorithm [@kaus2010stabilization]. Surface processes of erosion and sedimentation can also be incorporated in the
thermo-mechanical model. Complex boundary conditions for the velocity field, variable both in space and time,
can be adopted by the user to simulate different episodes of tectonism. Different benchmarks are available in the
repository and can be reproduced by the user (e.g. thermochemical convection – @van1997comparison; plume-lithosphere interaction – @crameri2012comparison).

As an example of application of `Mandyoc`, \autoref{fig:rift} presents snapshots of one numerical scenario of lithospheric stretching imposing a divergent flow direction, resulting in rifting and break-up. In this example, the upper crust, lower crust, lithospheric mantle and asthenosphere present different rheology and density, resulting in faulting mainly in the upper crust and part of the lithospheric mantle. Additioanally, deformations in the lower crust and at the base of the lithospheric mantle are accomodated by ductil creep flow. This example can be reproduced from the repository.


# Figures

![`Mandyoc` example of application of the thermo-mechanical model to simulate the stretching of the lithosphere, assuming different rheologies. Details can be found in the repository.\label{fig:rift}](JOSS_figure.png)

# Acknowledgements

We acknowledge contributions from Dave May for the correct implementation of the multigrid algorithm. This project was sponsored by FAPESP (Processes 2017/24870-5 and 2019/23246-1) and Petrobras (Process 2017/00461-9).

# References
