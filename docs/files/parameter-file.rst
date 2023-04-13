.. _parameterfile:

Parameter file 
==============

.. role:: raw-html(raw)
    :format: html

The parameter file ``param.txt`` contains the information that is necessary for the simulation to run. 

#. Geometry

    * nx
        default: :raw-html:`<br />` 
        type: integer :raw-html:`<br />` 
        unit: :raw-html:`<br />` 
        definition: number of nodes in the horizontal direction
    * nz
        default: :raw-html:`<br />` 
        type: integer :raw-html:`<br />` 
        unit: :raw-html:`<br />` 
        definition: number of nodes in the vertical direction
    * lx                                  
        default: :raw-html:`<br />` 
        type: real :raw-html:`<br />` 
        unit: m :raw-html:`<br />` 
        definition: extent in the horizontal direction
    * lz                                  
        default: :raw-html:`<br />`
        type: real :raw-html:`<br />`
        unit: m :raw-html:`<br />`
        definition: extent in the vertical direction

#. Simulation options

    * solver                              
        default: direct :raw-html:`<br />`
        type: direct/iterative :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set the solver to be direct or iterative
    * denok                               
        default: 1.0e-4 :raw-html:`<br />`
        type: real number :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: tolerance criterion for the Uzawa's scheme
    * rtol                                
        default: 1.0e-5 :raw-html:`<br />`
        type: real number :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: the absolute size of the residual norm (relevant only for iterative methods)
    * RK4                                 
        default: Euler  :raw-html:`<br />`
        type: Euler/Runge-Kutta :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: particles advection method
    * Xi_min                              
        default: 1.0e-14        :raw-html:`<br />`
        type: real number :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: tolerance criterion for the convergence of the non-linear flow
    * random_initial_strain               
        default: 0.0            :raw-html:`<br />`
        type: real number :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: non-dimensional value for the initial strain perturbation for the entire domain
    * pressure_const                      
        default: -1.0 (i.e. not used) :raw-html:`<br />`
        type: real number :raw-html:`<br />`
        unit: Pa :raw-html:`<br />`
        definition: set constant pressure value for the domain (relevant when 2-D is plain view)
    * initial_dynamic_range               
        default: False  :raw-html:`<br />`
        type: True/False :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: method to smoothen convergence of the velocity field in scenarios with wide viscosity range, see Gerya (2019) :cite:`gerya2019`
    * periodic_boundary                   
        default: False  :raw-html:`<br />`
        type: True/False :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: allows simulation with periodic boundary in the horizontal direction
    * high_kappa_in_asthenosphere         
        default: False  :raw-html:`<br />`
        type: True/False :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: mimics high heat transport in the asthenosphere increasing its thermal diffusivity coefficient
    * basal_heat                          
        default: -1.0 (i.e. not used)  :raw-html:`<br />`
        type: real number :raw-html:`<br />`
        unit: W/m^2 :raw-html:`<br />`
        definition: set basal heat flux value

#. Particles options

    * particles_per_element               
        default: 81           :raw-html:`<br />`
        type: integer :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: number of Lagrangian particles in each element
    * particles_per_element_x             
        default: 0 (automatic calculation) :raw-html:`<br />`
        type: integer :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: number of Lagrangian particles in the horizontal direction
    * particles_per_element_z             
        default: 0 (automatic calculation) :raw-html:`<br />`
        type: integer :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: number of Lagrangian particles in the vertical direction
    * particles_perturb_factor            
        default: 0.5  :raw-html:`<br />`
        type: real number between 0 and 1 :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: indicates the amount of perturbation of the initial location of the particles relative to a regular grid distribution. 

#. Surface processes

    * sp_surface_tracking                 
        default: False  :raw-html:`<br />`
        type: True/False :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: allows free surface tracking across time and outputs it
    * sea_level                           
        default: 0.0            :raw-html:`<br />`
        type: real number :raw-html:`<br />`
        unit: m :raw-html:`<br />`
        definition: sea level used to limit the influence of the surface process
    * sp_surface_processes                
        default: False  :raw-html:`<br />`
        type: True/False :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: allows erosion and sedimentation simulation
    * sp_dt                               
        default: 0              :raw-html:`<br />`
        type: real number :raw-html:`<br />`
        unit: years :raw-html:`<br />`
        definition: time step for surface processes simulation
    * a2l                                 
        default: True  :raw-html:`<br />`
        type:  True/False :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: allows the conversion of air particles to land particles during sedimentation
    * sp_mode                             
        default: 1              :raw-html:`<br />`
        type: integer :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: specify the surface processes method
    * free_surface_stab                   
        default: True  :raw-html:`<br />`
        type: True/False :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set if the free surface stabilization algorithm (FSSA) is used, see Kaus et al. (2010) :cite:`kaus2010`
    * theta_FSSA                          
        default: 0.5            :raw-html:`<br />`
        type: real number :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: weight of the influence of the FSSA method (only relevant when <free_surface_stab> is True)
    * sticky_blanket_air
        default: False  :raw-html:`<br />`
        type: True/False :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: allows the increase of viscosity for the first air layer of particles
    * precipitation_profile_from_ascii    
        default: False  :raw-html:`<br />`
        type: True/False :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set if precipitation profile along the horizontal axis is read from an ASCII file
    * climate_change_from_ascii           
        default: False  :raw-html:`<br />`
        type: True/False :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: if True, re-scales through time the precipitation profile using an ASCII file

#. Time constrains
    
    * step_max                            
        default: :raw-html:`<br />`
        type: integer :raw-html:`<br />`
        unit: steps :raw-html:`<br />`
        definition: maximum time-step of the simulation
    * time_max                            
        default:  :raw-html:`<br />`
        type: real number :raw-html:`<br />`
        unit: years :raw-html:`<br />`
        definition: maximum time of the simulation
    * dt_max                              
        default: :raw-html:`<br />`
        type: real number :raw-html:`<br />`
        unit: years :raw-html:`<br />`
        definition: maximum time between steps of the simulation 
    * step_print                          
        default:              :raw-html:`<br />`
        type: integer :raw-html:`<br />`
        unit: steps :raw-html:`<br />`
        definition: make output files every <step_print>
    * sub_division_time_step              
        default: 1.0            :raw-html:`<br />`
        type: real number :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: re-scale value for the calculated time-step
    * initial_print_step                  
        default: 0 (i.e. not used) :raw-html:`<br />`
        type: integer :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: <step_print> used until <initial_print_max_time>
    * initial_print_max_time              
        default: 1.0e6 :raw-html:`<br />`
        type: real number :raw-html:`<br />`
        unit: years :raw-html:`<br />`
        definition: maximum time to make output files every <initial_print_step>

#. Viscosity

    * viscosity_reference                 
        default:  :raw-html:`<br />`
        type: real number :raw-html:`<br />`
        unit: Pa.s :raw-html:`<br />`
        definition: reference mantle viscosity 
    * viscosity_max                       
        default:  :raw-html:`<br />`
        type: real number :raw-html:`<br />`
        unit: Pa.s :raw-html:`<br />`
        definition: maximum viscosity during simulation 
    * viscosity_min                       
        default:  :raw-html:`<br />`
        type: real number :raw-html:`<br />`
        unit: Pa.s :raw-html:`<br />`
        definition: minimum viscosity during simulation 
    * viscosity_per_element               
        default: constant  :raw-html:`<br />`
        type: constant/variable :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: sets if viscosity is constant or linearly variable for every element
    * viscosity_mean_method               
        default: harmonic  :raw-html:`<br />`
        type: harmonic/arithmetic :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: defines method do calculate the viscosity for each element
    * viscosity_dependence                
        default: depth  :raw-html:`<br />`
        type: pressure/depth :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: defines if viscosity depends on pressure or depth

#. External ASCII inputs/outputs

    * interfaces_from_ascii               
        default: False  :raw-html:`<br />`
        type: True/False :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set if interfaces between lithologies are read from an ASCII file (interfaces.txt)
    * n_interfaces                        
        default:        :raw-html:`<br />`
        type: :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set the number of interfaces to be read from the interfaces ASCII file (interfaces.txt) 
    * temperature_from_ascii              
        default: False  :raw-html:`<br />`
        type: True/False :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set if initial temperature is read from an ASCII file (input_temperature_0.txt)
    * velocity_from_ascii                 
        default: False  :raw-html:`<br />`
        type: True/False :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set if initial velocity field is read from an ASCII file (input_velocity_0.txt)
    * variable_bcv                        
        default: False  :raw-html:`<br />`
        type: True/False :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: allows velocity field re-scaling through time according to an ASCII file (scale_bcv.txt)
    * multi_velocity                      
        default: False  :raw-html:`<br />`
        type: True/False :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set if boundary velocities can change with time from ASCII file(s) (multi_veloc.txt and additional input_velocity_[X].txt files)
    * binary_output                       
        default: False  :raw-html:`<br />`
        type: True/False :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set if output is in binary format
    * print_step_files                    
        default: True  :raw-html:`<br />`
        type: True/False :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set if the particles position are printed to an output file

#. Physical parameters

    * temperature_difference              
        default:  :raw-html:`<br />`
        type: real number :raw-html:`<br />`
        unit: K :raw-html:`<br />`
        definition: temperature difference between the top and bottom of the model (relevant if <temperature_from_ascii> is False) 
    * thermal_expansion_coefficient       
        default:  :raw-html:`<br />`
        type: real number :raw-html:`<br />`
        unit: 1/K :raw-html:`<br />`
        definition: value for the coefficient of thermal expansion
    * thermal_diffusivity_coefficient     
        default:  :raw-html:`<br />`
        type: real number :raw-html:`<br />`
        unit: m^2/s :raw-html:`<br />`
        definition: value for the coefficient of thermal diffusivity 
    * gravity_acceleration                
        default:  :raw-html:`<br />`
        type: real number :raw-html:`<br />`
        unit: m/s^2 :raw-html:`<br />`
        definition: value for the gravity acceleration 
    * density_mantle                      
        default:  :raw-html:`<br />`
        type: real number :raw-html:`<br />`
        unit: kg/m^3 :raw-html:`<br />`
        definition: value for the mantle reference density 
    * heat_capacity                       
        default:  :raw-html:`<br />`
        type: real number :raw-html:`<br />`
        unit: J/K :raw-html:`<br />`
        definition: value for the heat capacity 
    * non_linear_method                   
        default:  :raw-html:`<br />`
        type: on/off :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set if non linear method is used for the momentum equation 
    * adiabatic_component                 
        default:  :raw-html:`<br />`
        type: on/off :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set if adiabatic heating/cooling is active 
    * radiogenic_component                
        default:  :raw-html:`<br />`
        type: on/off :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set if radiogenic heating is active 

#. Strain softening parameters

    * weakeing_min                 
        default: 0.05 :raw-html:`<br />`
        type: real number :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: value of the accumulated strain where strain softening starts :cite:`salazarmora2018`
    * weakening_max             
        default: 1.05 :raw-html:`<br />`
        type: real number :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: value of the accumulated strain where strain softening stops :cite:`salazarmora2018`

#. Velocity boundary conditions

    * top_normal_velocity                 
        default:            :raw-html:`<br />`
        type: fixed/free :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set the normal velocity on the top side of the model to be fixed or free 
    * top_tangential_velocity             
        default:            :raw-html:`<br />`
        type: fixed/free :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set the tangential velocity on the top side of the model to be fixed or free
    * bot_normal_velocity                 
        default:            :raw-html:`<br />`
        type: fixed/free :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set the normal velocity on the bottom side of the model to be fixed or free
    * bot_tangential_velocity             
        default:            :raw-html:`<br />`
        type: fixed/free :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set the tangential velocity on the bot side of the model to be fixed or free
    * left_normal_velocity                
        default:            :raw-html:`<br />`
        type: fixed/free :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set the normal velocity on the left side of the model to be fixed or free
    * left_tangential_velocity            
        default:             :raw-html:`<br />`
        type: fixed/free :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set the tangential velocity on the left side of the model to be fixed or free
    * right_normal_velocity               
        default:            :raw-html:`<br />`
        type: fixed/free :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set the normal velocity on the right side of the model to be fixed or free
    * right_tangential_velocity           
        default:             :raw-html:`<br />`
        type: fixed/free :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set the tangential velocity on the right side of the model to be fixed or free

#. Temperature boundary conditions

    * top_temperature                     
        default:            :raw-html:`<br />`
        type: fixed/free :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set temperature on the top side of the model to be fixed or free
    * bot_temperature                     
        default:            :raw-html:`<br />`
        type: fixed/free :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set temperature on the bottom side of the model to be fixed or free
    * left_temperature                    
        default:             :raw-html:`<br />`
        type: fixed/free :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set temperature on the left side of the model to be fixed or free
    * right_temperature                   
        default:             :raw-html:`<br />`
        type: fixed/free :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set temperature on the right side of the model to be fixed or free
    * rheology_model                      
        default:  :raw-html:`<br />`
        type: integer :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: flag number of a pre-defined rheology model to use during simulation
    * T_initial                           
        default:  :raw-html:`<br />`
        type: integer :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: flag number of a pre-defined temperature model to use during simulation (relevant when <temperature_from_ascii> is False)

Below there is an example of a parameter file ``param.txt`` used to make a simulation with *Mandyoc*.

.. literalinclude:: src/param.txt
   :language: text
   :linenos:



