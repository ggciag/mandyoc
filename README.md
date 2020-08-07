# MANDYOC #

MANtle DYnamics simulatOr Code

### Steps to compile the code ###
* [Install PETSc](https://www.mcs.anl.gov/petsc/) (copy the `config_LU.sh`
script in the PETSc folder and run this script)

* Clone the MANDYOC repository.

* In the `mandyoc/src/` directory, modify the first two lines of the
`Makefile` to the appropriate directory for PETSc.

* Compile the code in the `mandyoc/src/` directory:

    ```
    make all
    ```


### Steps to run the van Keken et al. (1997) example ###
* From the `src/` folder, copy the executable to the `examples/vanKeken1997/`
 folder:

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
You can easily load all output files produced by MANDYOC in Python through
[tapIOca](https://github.com/aguspesce/tapioca). It organize all output data
into a single [xarray.Dataset](https://xarray.pydata.org/en/stable/), which
makes much easier to manage, inspect and visualize the results.

You can find installation instructions and examples at the
[tapIOca website](https://github.com/aguspesce/tapioca).

## License
This is free software, you can redistribute it and/or modify it under the terms
 of the **BSD 3-clause License**.
 A copy of this license is provided in
 [LICENSE](/https://bitbucket.org/victorsacek/mandyoc/src/master/LICENSE).
