# MANDYOC #

MANtle DYnamics simulatOr Code

### Steps to compile the code ###
* [Install PETSc](https://www.mcs.anl.gov/petsc/) (copy the `config_LU.sh` script in
the PETSc folder and run this script)

* Clone the MANDYOC repository.

* In the `mandyoc/src/` directory, modify the first two lines of the `Makefile` to
the appropriate directory for PETSc.

* Compile the code in the `mandyoc/src/` directory:

    ```
    make all
    ```


### Steps to run the van Keken et al. (1997) example ###
* From the `src/` folder, copy the executable to the `examples/vanKeken1997/` folder:

    ```
    cp mandyoc ../examples/vanKeken1997/
    ```

    ```
    cd ../examples/vanKeken1997/
    ```

* Modify the path of the `mpirun` in `run.sh`.

* Run the `run.sh` script:

    ```
    sh run.sh
    ```

* To visualize the evolution of the density structure, run:

    ```
    ipython rho_imshow.py
    ```


## How to visualize the results using Xarray
You can easily transform all output files produced by MANDYOC to a single
[xarray.Dataset](https://xarray.pydata.org/en/stable/)
through [tapIOca](https://github.com/aguspesce/tapioca).
This makes managing and plotting data much easier.

You can find installation instructions and examples at the
[tapIOca website](https://github.com/aguspesce/tapioca).