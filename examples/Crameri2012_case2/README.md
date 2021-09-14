# Example: Crameri et al. (2012) - Case 2

The Case 2 experiment presented by Crameri et al. (2012) evaluates the sticky air method to obtain a numerical surface topography in geodynamic modelling.

The experiment analyses the change in topography due to the rising of a mantle plume.
The model setup consists of a 2800 km by 850 km box with a 150 km sticky air layer on the top of the model.
The mantle thickness is 600 km with a 100 km thick lithosphere.
The lithosphere density is 3300 kg/m<sup>3</sup> with viscosity 10<sup>23</sup> Pa.s,
the mantle density is 3300 kg/m<sup>3</sup> with viscosity 10<sup>21</sup> Pa.s
and the mantle plume density is 3200 kg/m<sup>3</sup> with viscosity 10<sup>20</sup> Pa.s.

Initially, the center of the plume is horizontally centered and 300 km above the base of the model.
At the top, the sticky air layer has density 0 kg/m<sup>3</sup> with viscosity 10<sup>19</sup> Pa.s.
A free slip boundary condition is applied to the upper boundary of the sticky air layer and the vertical sides of the model and the base are kept fixed.
There is no temperature difference, and the geodynamic evolution is guided solely by compositional density differences.

## How to run this example

### Generate the input files

To generate the interface, velocity, temperature and parameter file, you need to run `generate_input_file.py` as:
```
python generate_input_files.py
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


### Plot viscosity results

To plot the result, run `plot_viscosity.py` as:
```
python plot_viscosity.py
```

## Compare the result with other programs

You can compare the result obtained using the Mandyoc code with the results of other codes, such as UNDERWORLD, I2VES, and STAGYY, for the same model, run:
```
python plot_max_surface.py
```

## References

F Crameri, H Schmeling, GJ Golabek, T Duretz, R Orendt, SJH Buiter, DA May, BJP Kaus, TV Gerya, and PJ Tackley. A comparison of numerical
surface topography calculations in geodynamic modelling: an evaluation of the 'sticky air' method.
Geophysical Journal International, 189(1):38–54, 2012.
