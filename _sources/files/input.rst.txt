.. _inputfiles:

Input files
===========

*Mandyoc* accepts ASCII files to use as initial conditions for the simulation, such as an initial temperature file ``input_temperature_0.txt``, a file containing the initial interfaces geometry ``interfaces.txt`` and a file containing the initial velocity field at the faces of the model ``input_velocity_0.txt``. The next subsections will help you creating each one of these files.


ASCII temperature file
----------------------

..
   For both 2-D and 3-D grids, the initial temperature configuration can be provided as an ASCII file called ``intial-temperature.txt``. Considering :math:`x` to be the longitudinal direction, :math:`y` the vertical direction, and :math:`z` the latitudinal direction, the next subsections will guide you through the understanding of the initial temperature file structure.

For a 2-D grid, the initial temperature configuration can be provided as an ASCII file called ``input_temperature_0.txt``. Considering :math:`x` to be the longitudinal direction and :math:`y` the vertical direction, the next subsections will guide you through the understanding of the initial temperature file structure.

2-D temperature grid
********************

Consider a 2-D grid where the number of nodes in the :math:`x` and :math:`y` directions are :math:`nx\geq 2` and :math:`ny\geq 2`, respectively, with :math:`(nx-1)` elements in the :math:`x` direction and :math:`(ny-1)` elements in the :math:`y` direction. Each node in the grid can be identified with a coordinates pair :math:`(x_n, y_m)`, where :math:`n` and :math:`m` are natural numbers, and :math:`x_n` and :math:`y_m` are the :math:`x` and :math:`y` positions (in meters) of any grid node. :numref:`coordinates2d` shows a scheme for a 2-D grid, where the red dots represent the grid nodes and every node is labeled with a :math:`(x_n, y_m)` pair. The dotted lines represent any number of intermediate nodes.

.. Consider a 2-D grid with :math:`(nx-1)>2` elements in the :math:`x` direction and :math:`(ny-1)>2` elements in the :math:`y` direction. The number of nodes in the :math:`x` and :math:`y` directions are :math:`nx` and :math:`ny`, respectively, and each node can be identified with a pair :math:`(x_n, y_m)`, where :math:`n \leq (nx-1) \in \mathbb{N}` and :math:`m \leq (ny-1) \in \mathbb{N}`, where :math:`x_n` and :math:`y_m` are the :math:`x` and :math:`y` positions (in meters) of any grid node. :numref:`coordinates2d` shows a scheme for a 2-D grid, where the red dots represent the grid nodes, the the dotted lines represent any number of element nodes, and every node is labeled with a :math:`(x_n, y_m)` pair.

.. _coordinates2d:

.. figure:: figs/coordinates-2d.png
   :width: 90%
   :align: center
   :alt: Coordinates

   2-D grid scheme. The red dots represent the grid nodes where the values for the temperature must be defined, the dotted lines represent any number of grid nodes, for simplification, and the labels indicate the :math:`x` and :math:`y` position of each node.

..
   .. note::
      Because of the way *Mandyoc* deals with the coordinates in 2-D, the node at :math:`(x_0,y_0)` is always at the origin :math:`(0,0)`. This is also true for the 3-D grid, where the :math:`(x_0,y_0,z_0)` is at :math:`(0,0,0)`.

.. note::
   Because of the way *Mandyoc* deals with the coordinates in 2-D, the node at :math:`(x_0,y_0)` is always at the origin :math:`(0,0)`.

The example below shows how the ``input_temperature_0.txt`` file must be written for the 2-D grid in :numref:`coordinates2d`: after writing four lines of comments, the file must include the temperature :math:`T(x_n, y_m)` at each node starting at the bottom left of the model and going up in the :math:`y` direction. Intuitively, the temperature is given in horizontal layers for every :math:`x_n`.

.. literalinclude:: src/initial-temperature-2d.txt
   :language: text
   :linenos:

.. tip::
   It is always helpful to use the four comment lines to write down any useful information for later. Because these lines are simply skipped, there is no rule to write them, the '#' symbol is used only by convention.

..
   3-D temperature grid
   ********************

   Constructing the ``input_temperature.txt`` file for a three dimensional grid (:numref:`coordinates3d`) is very similar to the 2-D case, except it is necessary to provide the temperature of each horizontal layer of nodes :math:`xz` for each :math:`y=y_0,y_1,...,y_{ny-1}`. For the 3-D grid in :numref:`coordinates3d`, the number of nodes in the :math:`z` direction is called :math:`nz` and the temperature at each node is :math:`T(x_n, y_m, z_p)`, where :math:`p \leq (nz-1) \in \mathbb{N}`.

   .. _coordinates3d:

   .. figure:: figs/coordinates-3d.png
      :width: 90%
      :align: center
      :alt: Coordinates

      3-D grid scheme. The red dots represent the grid nodes where the values for the temperature must be defined. The dotted lines represent any number of grid nodes, for simplification.

   The example below shows the ``input_temperature.txt`` file for the generic 3-D grid in :numref:`coordinates3d`. Following the four lines of comments, the temperature value at each node must be provided in each line.

   .. literalinclude:: src/initial-temperature-3d.txt
      :language: text
      :linenos:

ASCII interfaces file
---------------------

..
   The ``input_interfaces.txt`` file, for both 2-D and 3-D grids, starts with seven lines of variables that are used by the rheology models. The variables that are in the interfaces are those that might be used in :eq:`power-law` or :eq:`frank-kamenetskii`: :math:`C`, :math:`\rho_r`, :math:`H`, :math:`A`, :math:`n`, :math:`Q`, and :math:`V`. Each one of these variables possesses a number of values (columns) that are assigned to each lithological unit in the model. Every interface is a boundary between two lithological units, therefore the number of columns is 1 plus the number of interfaces set in the ``param.txt`` file (see the :doc:`parameter file section<parameter-file>`), such that if there are :math:`i` interfaces, :math:`i+1` values for each variable **must** be provided in the ``input_interfaces.txt`` file. More will be discussed about the order of these values shortly.

The ``interfaces.txt`` file, for a 2-D grid, starts with seven lines of variables that are used by the rheology models. The variables that are in the interfaces are those that might be used in :eq:`power-law` or :eq:`frank-kamenetskii`: :math:`C`, :math:`\rho_r`, :math:`H`, :math:`A`, :math:`n`, :math:`Q`, and :math:`V`. Each one of these variables possesses a number of values (columns) that are assigned to each lithological unit in the model. Every interface is a boundary between two lithological units, therefore the number of columns is 1 plus the number of interfaces set in the ``param.txt`` file (see the :doc:`parameter file section<parameter-file>`), such that if there are :math:`i` interfaces, :math:`i+1` values for each variable **must** be provided in the ``interfaces.txt`` file. More will be discussed about the order of these values shortly.


2-D Initial interfaces
**********************

Below, the example corresponds to a 2-D grid with two interfaces and, therefore, three lithological units. The first column contains the vertical positions **in meters** of every grid node :math:`y_m` that corresponds to the **deepest** interface boundary, starting at :math:`x=x_0` on line 8, and ending at :math:`x=x_{nx-1}` on line :math:`nx+7`. The second column contains the vertical position of every :math:`y_m` that corresponds to the second interface boundary. When defining the interfaces, it is rather common for them to "touch". Because of that, all the interfaces must be provided in a "tetris" manner, where interfaces that are collinear in parts fit the interface below.

.. note::
   Because the interfaces are defined linearly between the nodes, it is important to define them properly, so every point inside the grid can be attributed to a lithological unit.

.. literalinclude:: src/initial-interfaces-2d.txt
   :language: text
   :linenos:

The values of the variables in the first seven lines are the values of the lithological units bound by the interfaces. For the 2-D grid with two interfaces, the first values of each variable refers to the lithological unit below the first interface, the second value of each variable refers to the lithological unit between the first and second interfaces, and the third value of each variable refers to the lithological unit above the second interface. If more interfaces are added, there will be more units bounded by an upper and a lower interface.

To consider different values of internal cohesion :math:`c_0` and internal angle of friction :math:`\varphi` for each layer, the interface file can be written as the file below, where *weakening_seed*, *cohesion_min*, *cohesion_max*, *friction_angle_min*, and *friction_angle_max* are given.

.. literalinclude:: src/new-initial-interfaces-2d.txt
   :language: text
   :linenos:

The *weakening_seed* value for each layer represents its initial accumulated strain value (seed), and when it is negative, a random seed will be used instead. When *cohesion_min* and *cohesion_max* are different, that corresponding layer will be linearly softened between those values. Notice that a constant cohesion can be set if the values are simply equal. The same logic aplies to *friction_angle_min* and *friction_angle_max*.

The default interval of accumualted strain where strain softening occurs is 0.05 and 1.05, cf. :cite:`salazarmora2018`. To use different values, the user must provide them in the ``param.txt`` withe the keywords *weakeing_min* and *weakening_max*.

..
   3-D Initial interfaces
   **********************

   The 3-D grid ``input_interfaces.txt``  file is very similar to the 2-D grid one. The difference lies that after completing the interfaces for :math:`z=z_0`, another set of values is added subsequently to :math:`z=z_1` and so on, until :math:`z=z_{nz-1}`. The ``input_interfaces.txt`` file below shows and example of such case.

   .. literalinclude:: src/initial-interfaces-3d.txt
      :language: text
      :linenos:


ASCII velocities files
----------------------

..
   For both 2-D and 3-D grids, the ``input_velocity.txt`` file possesses four lines of headers followed by the velocity data for each node of the model.

An initial velocity file allows an initial velocity field to be established so the simulated fluid presents an initial momentum. Respecting the conservation of mass equation when designing the initial velocity field is fundamental so the simulation can run properly.

Additionally, the velocity boundary conditions defined in the ``param.txt`` (see the :doc:`parameter file section<parameter-file>`) will be used to define the normal and tangential velocities on the faces of the model during the simulation.

2-D initial velocity
********************

When ``velocity_from_ascii = True`` is set in the ``param.txt`` file, an initial velocity field must be provided in an ASCII file called ``input_velocity_0.txt``.

For a 2-D grid, the ``input_velocity_0.txt`` file possesses four lines of headers followed by the velocity data for each node of the model. The velocity must be defined for both :math:`x` and :math:`y` directions. Taking the example in the file below, the ``input_velocity_0.txt`` file must be written in such a way that the values for the velocity in the :math:`x` direction :math:`v_x` and the velocity in the :math:`y` direction :math:`v_y` for the node :math:`(x_0,y_0)` corresponds to lines 5 and 6, respectively, the velocity values for the :math:`(x_1,y_0)` corresponds to lines 7 and 8, and so on until :math:`(x_{nx-1},y_0)`. Once all the nodes in the :math:`y_0` are written, the values for all :math:`x_n` are inserted for :math:`y_1`, and so on until :math:`y_{ny-1}`.

.. literalinclude:: src/initial-velocity-2d.txt
   :language: text
   :linenos:

2-D multiple initial velocities
*******************************

The user can provide multiple velocity fields to be set in different instants of the simulation. By setting ``multi_velocity = True`` in the ``param.txt`` file, an ASCII file called ``multi_veloc.txt`` must be provided following the structure shown in the example file below. The first line of the file contains the amount of instants where a different velocity field will be used and the other lines contain the time, in millions of years, that each velocity field will be adopted. In the example below, a different velocity field will be used at three instants, the files must follow the structure presented in the last section and should be named in sequence: ``input_velocity_1.txt``, ``input_velocity_2.txt`` and ``input_velocity_3.txt``. The names of the files are labeled after the three different instants set in the first line of the ``multi_veloc.txt`` file.

.. literalinclude:: src/multi_veloc.txt
   :language: text
   :linenos:

.. note::
   The ``input_velocity_0.txt`` is used for the beginning of the simulation. The multiple initial velocities that can be set during the simulation start at ``input_velocity_1.txt`` and end at ``input_velocity_X.txt``, where :math:`X` is the number of different instants that the velocity field will be changed. 

2-D re-scalable boundary condition
**********************************

When ``variable_bcv = True`` is set in the ``param.txt``file, an ASCII file called ``scale_bcv.txt`` is used to define a sequence of re-scaling steps for the velocity boundary conditions. The example file below is an example of how such file can be written. In the example file, the first line corresponds to the number of instants during the simulation that the velocity field will be re-scaled. The first column contains the time, in millions of years, that the velocity field will be changed and the second column contains the scale value. For the example below, in the instant :math:`10.0` Myr the velocity field will be multiplied by :math:`0.5`, effectively halving it; at :math:`25.0` Myr the current velocity field will have its direction changed, hence multiplied by :math:`-1.0`; and finally at :math:`50` My the velocity field will have its direction changed again and doubled by multiplying it by :math:`-2.0`.

.. literalinclude:: src/scale_bcv.txt
   :language: text
   :linenos:


..
   3-D initial velocity
   ********************

   Similarly to the 2-D grid, the ``input_velocity_0.txt`` file for a 3-D grid must contain the velocity values for each :math:`x`, :math:`y` and :math:`z` components. Therefore, the resulting files differs from the 2-D case as shown in the example below, which specifies the velocity for each node in the :math:`xz` plane for each depth :math:`y_m`.

   .. literalinclude:: src/initial-velocity-3d.txt
      :language: text
      :linenos:


