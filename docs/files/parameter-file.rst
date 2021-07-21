.. _parameterfile:

Parameter file 
==============

The parameter file ``param.txt`` contains the information that is necessary for the simulation to run. 

#. Geometry

    * Number of elements in each direction 
    * The extension of the model 

#. Time constrains

    * The time-step of the simulation
    * Maximum time-step of the simulation
    * Maximum time of the simulation (in seconds)
    * Maximum time between steps of the simulation (in seconds)

#. Viscosity values
    
    * The reference viscosity
    * The minimum viscosity
    * The maximum viscosity

#. Physical parameters    

    * Gravity acceleration
    * Heat capacity
    * Coefficient of thermal conductivity 
    * Coefficient of thermal expansion

#. Boundary conditions

    * Temperature boundary condition at all borders
    * Tangential velocity boundary conditions at all borders
    * Normal velocity boundary conditions at all borders

#. Rheology and temperature model numbers

#. External initial temperature from ASCII file

#. External initial interfaces from ASCII file

#. External initial velocity files from ASCII file

#. Type of solver

    * Direct
    * Iteractive 

#. Simulation options

    * denok
    * rtol

Below there is an example of a parameter file ``param.txt`` used to make a simulation with MANDYOC.

.. literalinclude:: src/param.txt
   :language: text
   :linenos: