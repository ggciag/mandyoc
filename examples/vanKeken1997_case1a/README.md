# Example: van Keken

The set of simulations proposed by van Keken et al. (1997) compares
several methods of studying two dimensional thermochemical convection,
where the Boussinesq approximation and infinite Prandtl number are used.

The simulation consists of two layers, where a buoyant thin layer is
under a denser thicker package.

## How to run this example

### Dependencies

To run this example you need to install:

* `numpy`

### Generate the input files

To generate the interface, velocity, temperature and parameter file, you need to run `generate_interface.py` as:
```
python generate_interfaces.py
```

### Run the model

Now, you can run the model as:
```
~/petsc/arch-label-optimized/bin/mpirun -n NUMBER_OF_CORES ../../mandyoc
```
Or you can use the script called `run.sh` to run it:
```
sh run.sh
```
__You have to change `NUMBER_OF_CORES`.__ 

Remember that `PETSc` is installed by default in `~/`. 
If you have changed it you must adjust the path to _mpirun_ accordingly.


### Plot results

To plot the result, run `view_density.py` as:
```
python view_density.py
```

## References

PE Van Keken, SD King, H Schmeling, UR Christensen, D Neumeister, and
M-P Doin. A comparison of methods for the modeling of thermochemical
convection. Journal of Geophysical Research: Solid Earth,
102(B10):22477–22495, 1997.
