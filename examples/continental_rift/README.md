# Example: Continental rift

This example simulates the evolution of divergent margins, taking into account the plastic rheology and the sin-rift geodynamics.

The domain of the model comprises 1600 x 300 km<sup>2</sup>, composed of a regular mesh with square elements of 1 x 1 km<sup>2</sup>.
The boundary conditions for the velocity field simulate the lithospheric stretching
assuming a reference frame fixed on the lithospheric plate on the left side of the model,
and the plate on the right side moves rightward with a velocity of 1 cm/year.
The velocity field in the left and right boundaries of the model is chosen to ensure conservation of mass
and is symmetrical if the adopted reference frame movies to the right with a velocity of 0.5 cm/year relative to the left plate.
Additionally, free slip condition was assumed on the top and bottom of the numerical domain.
To simulate the free surface, the "sticky air" approach (e.g. Gerya and Yuen, 2003b) is adopted,
taking into account a 40-km thick layer with a relatively low viscosity material but with a compatible density with the atmospheric air.
The initial temperature structure is only depth dependent and is 0 ºC at the surface and 1300 ºC at the base of the lithosphere at 130 km.
Please, see the `generate_input_files.py` to check out how the initial temperature structure in the interior of the lithosphere is obtained.

To avoid artifacts created by a homogeneous rheology, a random perturbation of the initial strain in each finite element of the model (e.g. Brune et al., 2014) is applied.
This random perturbation follows a normal distribution in which the mean initial strain is 0.25 with a standard deviation of 0.08.
Additionally, to ensure the nucleation of rifting at the center of the numerical domain,
a weak seed (e.g. Huismans and Beaumont, 2003) is present in the lithospheric mantle with a constant initial strain of 0.3.


## How to run this example

### Dependencies

To run this example you need to install:

* `numpy`
* `matplotlib`
* `pandas`
* `sys`
* `glob`

### Generate the input files

To generate the interface, velocity, temperature and parameter file, you need to run `generate_input_file.py` as:
```
python generate_input_files.py
```

### Run the model

In this example, `mandyoc` use the following flags:

* `-seed 0,2`

* `-strain_seed 0.0,1.0`

You can run the model as:
```
~/petsc/arch-label-optimized/bin/mpirun -n NUMBER_OF_CORES ../../mandyoc -seed 0,2 -strain_seed 0.0,1.0
```

Or you can use the script called `run.sh` to run it:
```
sh run.sh
```
__You have to change `NUMBER_OF_CORES`.__ 

Remember that `PETSc` is installed by default in `~/`. 
If you have changed it you must adjust the path to _mpirun_ accordingly.

### Plot results

To plot the result, run `view_litho_strain_temper.py` as:
```
python view_litho_strain_temper.py
```


## References
Brune S., Heine C., Pérez-Gussinyé M., Sobolev S. V., Rift migration explains continental margin asymmetry and crustal hyper-extension,
Nature communications, 2014, vol. 5, p. 1.

Gerya T. V., Yuen D. A., Characteristics-based marker-in-cell method with conservative finite-differences schemes for modeling geological flows with strongly variable transport properties, Physics of the Earth and Planetary Interiors, 2003a, vol. 140, p. 293

Huismans R. S., Beaumont C., Symmetric and asymmetric lithospheric extension: Relative effects of frictional-plastic and viscous strain softening, Journal of Geophysical Research: Solid Earth, 2003, vol. 108
