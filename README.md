# MANDYOC

MANtle DYnamics simulatOr Code

## Installation

Before installing MANDYOC, we need to _install its dependencies_:

- [PETSc](https://www.mcs.anl.gov/petsc/)
- `gcc`
- `make`
- `git`

Optional dependencies:

- `gfortran`

### Building and installing PETSc

First, _download_ the latest release of PETSc [from their
[website](https://www.mcs.anl.gov/petsc/download/index.html) or clone their
repository:

```
git clone -b release https://gitlab.com/petsc/petsc.git $HOME/petsc
```

By default, we will download/clone in `~/petsc`.

Second, _configure the PETSc build_ and set up the installation directory.
By default, we will install PETSc in `~/petsc`.

```
cd $HOME/petsc
./configure \
  PETSC_DIR=$HOME/petsc \
  PETSC_ARCH=arch-label-optimized \
  --with-debugging=0 \
  --with-cc=gcc \
  --with-cxx=g++ \
  --download-fblaslapack \
  --download-mpich \
  --download-hdf5 \
  --download-superlu_dist \
  --download-metis \
  --download-parmetis \
  --download-mumps \
  --download-scalapack \
  --download-cmake \
  COPTFLAGS='-O3 -march=native -mtune=native' \
  CXXOPTFLAGS='-O3 -march=native -mtune=native'
```

> Note: If using `gfortran` optional dependency add the options
`--with-fc=gfortran` and `FOPTFLAGS='-O3 -march=native -mtune=native'`
to the PETSc build configuration above.

> Note: If you are build a development version of MANDYOC you can build
a debug version of PETSc by setting `--with-debugging=1` and removing the
COPTFLAGS, CXXOPTFLAGS (and FOPTFLAGS) flags.
In this case, you may set `PETSC_ARCH=arch-label-debug`.

_Check_ the installation with:

```
make all check
```

### Build and install MANDYOC

Build and install MANDYOC by running:

```
make all
```

> Note: To print MANDYOC runtime options, run mandyoc with `-flags` command line argument.

## Examples

### Steps to run the van Keken et al. (1997) example

- From the `src/` folder, copy the executable to the `examples/vanKeken1997/`
  folder:

  ```
  cp mandyoc ../examples/vanKeken1997/
  ```

  ```
  cd ../examples/vanKeken1997/
  ```

- Modify the path of the `mpirun` in `run.sh`.

- Run the `run.sh` script:

  ```
  sh run.sh
  ```

- To visualize the evolution of the density structure, run:

  ```
  ipython rho_imshow.py
  ```

## Contributing

### Code of conduct

Please note that this project is released with a
[Contributor Code of Conduct](https://bitbucket.org/victorsacek/mandyoc/src/master/CODE_OF_CONDUCT.md).
By participating in this project you agree to abide by its terms.

### Contributing Guidelines

**Contributions to MANDYOC are welcome**.
If you have an issue, a bug, a code contribution or a documentation
contribution, _thanks for helping to improve MANDYOC!_
Check out our
[Contributing guide](https://bitbucket.org/victorsacek/mandyoc/src/master/CONTRIBUTING.md).

## How to visualize the results using Xarray

You can easily load all output files produced by MANDYOC in Python through
[tapIOca](https://github.com/aguspesce/tapioca).
It organize all output data into a single
[xarray.Dataset](https://xarray.pydata.org/en/stable/), which makes much easier
to manage, inspect and visualize the results.

You can find installation instructions and examples at the
[tapIOca website](https://github.com/aguspesce/tapioca).

## License

This is free software, you can redistribute it and/or modify it under the terms
of the **BSD 3-clause License**.
A copy of this license is provided in
[LICENSE](/https://bitbucket.org/victorsacek/mandyoc/src/master/LICENSE).
