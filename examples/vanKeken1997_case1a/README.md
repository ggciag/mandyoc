# Example: van Keken

The set of simulations proposed by van Keken et al. (1997) compares
several methods of studying two dimensional thermochemical convection,
where the Boussinesq approximation and infinite Prandtl number are used.

The simulation consists of two layers, where a buoyant thin layer is
under a denser thicker package.

## How to run this example

### Generate the input files

To generate the interface, velocity, temperature and parameter file, you
need to run `generate_interface.py` as:

    python generate_interfaces.py

### Run the model

Now, you can run the model as:

    ${PETSC_DIR}/${PETSC_ARCH}/bin/mpirun -n number_of_cores ../../mandyoc

**You have to change `number_of_cores`.** If `PETSC_DIR` and `PETSC_ARCH` environment variables are not defined, you must adjust the path to _mpirun_ accordingly.

Or you can use the script called `run.sh` to run it. __Remember to update the `NUMBER_OF_CORES`__ variable in the script.

    sh run.sh

### Plot results

To plot the result, run `view_density.py` as:

    python view_density.py

## References

PE Van Keken, SD King, H Schmeling, UR Christensen, D Neumeister, and
M-P Doin. A comparison of methods for the modeling of thermochemical
convection. Journal of Geophysical Research: Solid Earth,
102(B10):22477–22495, 1997.
