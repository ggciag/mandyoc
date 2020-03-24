# MANDYOC #

MANtle DYnamics simulatOr Code

### Steps to compile the code ###

* [Install PETSc](https://www.mcs.anl.gov/petsc/)
* Clone the repository
* In the mandyoc/src/ directory, modify the first two lines of the makefile to the appropriate directory for PETSc.
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

