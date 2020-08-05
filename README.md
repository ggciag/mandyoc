# MANDYOC #

MANtle DYnamics simulatOr Code

### Steps to compile the code ###
* [Install PETSc](https://www.mcs.anl.gov/petsc/) (copy the config_LU.sh script in
the PETSc folder and run this script)

* Clone the MANDYOC repository.

* In the mandyoc/src/ directory, modify the first two lines of the makefile to
the appropriate directory for PETSc.

* Compile the code in the mandyoc/src/ directory:

    ```bash
    make all
    ```


### Steps to run the van Keken et al. (1997) example ###
* From the src/ folder, copy the executable to the examples/vanKeken1997/  folder:

    ```bash
    cp mandyoc ../examples/vanKeken1997/

    cd ../examples/vanKeken1997/
    ```

* Modify the path of the mpirun in run.sh.

* Run the run.sh script:

    ```bash
    sh run.sh
    ```

* To visualize the evolution of the density structure, run:

    ```bash
    ipython rho_imshow.py
    ```


## How to visualize the results using Xarray
Using the [tapIOca](https://github.com/aguspesce/tapioca) repository, you can
transform all output files from MAMDYOC code to an only
[xarray.Dataset](https://xarray.pydata.org/en/stable/).

To know how to use this tool, go to [tapIOca](https://github.com/aguspesce/tapioca)
web site.
There you will find the installation data and example of how to use it.
