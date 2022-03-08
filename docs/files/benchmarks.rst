.. _benchmarks:

Benchmarks
==========

In order to test the accuracy of the *Mandyoc* code, its results can be compared to benchmark studies. The following subsections will present the procedures and results for some well established modeling problems. Such cases should provide information about the applicability and performance of the code as well as some of its limitations.

The scripts to build and run these numerical experiments are located in the `examples <https://github.com/ggciag/mandyoc/tree/main/examples>`_ folder of the Mandyoc repository.
Inside each example folder, you find a ``README.md`` file with detailed explanation and instrutions on how to run the experiment.

van Keken et al. (1997) :cite:`vankeken1997`
--------------------------------------------

.. note::

  Mandyoc version used for this benchmark: v0.1.4.

The set of simulations proposed by van Keken et al. (1997) :cite:`vankeken1997` compares several methods of studying two dimensional thermochemical convection, where the Boussinesq approximation and infinite Prandtl number are used.

For this benchmark, the first case is presented with the three distinct variations, as proposed by the article. The simulation consists of two layers, where a buoyant thin layer is under a denser thicker package. The problem can be interpreted as a salt layer under a sediment package, and the interface between the layers is defined by the :eq:`interfacerayleigh` below.

.. math::
    :label: interfacerayleigh

    y=-0.8 \lambda_y + 0.02 \cos{\frac{\pi x }{\lambda_x}}

where :math:`\lambda_x` and :math:`\lambda_y` are the horizontal and vertical lengths of the simulated 2-D box, respectively.

The simulations are carried out in a Cartesian box where the fluid is isothermal and Rayleigh-Taylor instability is expected for the proposed setup. The table below lists the parameters used to run this scenario together with the distinction between the three proposed simulation (*case 1a*, *1b* and *1c*).

.. list-table:: Parameters used for the Rayleigh-Taylor instability simulation.
    :header-rows: 1
    :widths: 30 20 20
    :align: center

    * - Parameter
      - .. centered:: Symbol
      - Value
    * - Horizontal length
      - .. centered:: :math:`\lambda_x`
      - 1.0000
    * - Vertical length
      - .. centered:: :math:`\lambda_y`
      - 0.9142
    * - Thermal diffusion coefficient
      - .. centered:: :math:`\kappa`
      - :math:`1.0\times 10^{-6}`
    * - Gravity acceleration
      - .. centered:: :math:`g`
      - :math:`10`
    * - Reference viscosity
      - .. centered:: :math:`\eta_r`
      - :math:`1.0\times 10^{21}`
    * - Buoyant layer viscosity
      - .. centered:: :math:`\eta_0`
      - | :math:`1.00\times\eta_r` (case 1a)
        | :math:`0.10\times\eta_r` (case 1b)
        | :math:`0.01\times\eta_r` (case 1c)

Results for *case 1a*
*********************

For the *case 1a* where :math:`\eta_0/\eta_r=1.00`, :numref:`vankekenCase1aEvolution` below compares the evolution of the isoviscous Rayleigh-Taylor instability between the van Keken et al. (1997) :cite:`vankeken1997` and the *Mandyoc*. The time steps shown for the *Mandyoc* code are the closest the simulation could provide, considering the chosen simulation parameters.

.. _vankekenCase1aEvolution:

.. figure:: figs/vankeken-snaps-1a.png
  :align: center
  :width: 80%
  :alt: Results

Evolution of the isoviscous Rayleigh-Taylor instability for :math:`\eta_0/\eta_r=1.00`. The best result presented by van Keken et al. (1997) :cite:`vankeken1997` are on the left and the *Mandyoc* results are on the right.

.. note::
  Because of the different methods used by van Keken et al. (1997) :cite:`vankeken1997` and *Mandyoc*, the *Mandyoc* results for the evolution of the isoviscous Rayleigh-Taylor instability presents its data colored instead of contoured.

:numref:`vankekenCase1aGraph` below compares the change of the :math:`v_{rms}` with time, showing the results from van Keken et al. (1997) :cite:`vankeken1997` in gray and *Mandyoc* in black.

.. _vankekenCase1aGraph:

.. figure:: figs/vrms-1a.png
  :align: center
  :width: 100%
  :alt: Results

  Evolution of the :math:`v_{rms}` for :math:`\eta_0/\eta_r=1.00`. The van Keken et al. (1997) :cite:`vankeken1997` result is shown in black and the *Mandyoc* code result is shown in gray.

Results for *case 1b*
*********************

For the *case 1b* where :math:`\eta_0/\eta_r=0.10`, :numref:`vankekenCase1bEvolution` compares the evolution of the isoviscous Rayleigh-Taylor instability between van Keken et al. (1997) :cite:`vankeken1997` and *Mandyoc*. The time steps shown for the *Mandyoc* code are the closest the simulation could provide, considering the chosen simulation parameters.

.. _vankekenCase1bEvolution:

.. figure:: figs/vankeken-snaps-1b.png
  :align: center
  :width: 80%
  :alt: Results

  Evolution of the isoviscous Rayleigh-Taylor instability for :math:`\eta_0/\eta_r=0.10`. The best result presented by van Keken et al. (1997) :cite:`vankeken1997` are on the left and the *Mandyoc* results are on the right.

:numref:`vankekenCase1bGraph` below compares the change of the :math:`v_{rms}` with time, showing the results from van Keken et al. (1997) :cite:`vankeken1997` in gray and *Mandyoc* in black.

.. _vankekenCase1bGraph:

.. figure:: figs/vrms-1b.png
  :align: center
  :width: 100%
  :alt: Results

  Evolution of the :math:`v_{rms}` for :math:`\eta_0/\eta_r=0.10`. The van Keken et al. (1997) :cite:`vankeken1997` result is shown in black and the *Mandyoc* code result is shown in gray.

Results for *case 1c*
*********************

For the *case 1c* where :math:`\eta_0/\eta_r=0.01`, :numref:`vankekenCase1cEvolution` compares the evolution of the isoviscous Rayleigh-Taylor instability between van Keken et al. (1997) :cite:`vankeken1997` and *Mandyoc*. The time steps shown for the *Mandyoc* code are the closest the simulation could provide, considering the chosen simulation parameters.

.. _vankekenCase1cEvolution:

.. figure:: figs/vankeken-snaps-1c.png
  :align: center
  :width: 80%
  :alt: Results

  Evolution of the isoviscous Rayleigh-Taylor instability for :math:`\eta_0/\eta_r=0.01`. The best result presented by van Keken et al. (1997) :cite:`vankeken1997` are on the left and the *Mandyoc* results are on the right.

:numref:`vankekenCase1cGraph` below compares the change of the :math:`v_{rms}` with time, showing the results from van Keken et al. (1997) :cite:`vankeken1997` in gray and *Mandyoc* in black.

.. _vankekenCase1cGraph:

.. figure:: figs/vrms-1c.png
  :align: center
  :width: 100%
  :alt: Results

  Evolution of the :math:`v_{rms}` for :math:`\eta_0/\eta_r=0.01`. The van Keken et al. (1997) :cite:`vankeken1997` result is shown in black and the *Mandyoc* code result is shown in gray.

Crameri et al. (2012) :cite:`crameri2012`
-----------------------------------------

.. note::

  Mandyoc version used for this benchmark: v0.1.4.

The *Case 2* experiment presented by Crameri et al. (2012) :cite:`crameri2012` evaluates the *sticky air* method to obtain a numerical surface topography in geodynamic modelling.

The experiment analyses the change in topography due to the rising of a mantle plume.
The model setup (:numref:`crameri_setup`) consists of a :math:`2800 \, \mathrm{km}` by :math:`850 \, \mathrm{km}` box with a :math:`150 \, \mathrm{km}` sticky air layer on the top of the model.
The mantle thickness is :math:`600 \, \mathrm{km}` with a :math:`100 \, \mathrm{km}` thick lithosphere.
The lithosphere density is :math:`3300 \, \mathrm{kg/m}^3` with viscosity :math:`10^{23} \, \mathrm{Pa\,s}`,
the mantle density is :math:`3300 \, \mathrm{kg/m}^3` with viscosity :math:`10^{21} \, \mathrm{Pa\,s}`
and the mantle plume density is :math:`3200 \, \mathrm{kg/m}^3` with viscosity :math:`10^{20} \, \mathrm{Pa\,s}`.
Initially, the center of the plume is horizontally centered and :math:`300 \, \mathrm{km}` above the base of the model.
At the top, the sticky air layer has density :math:`0 \, \mathrm{kg/m}^3` with viscosity :math:`10^{19} \, \mathrm{Pa\,s}`.
A free slip boundary condition is applied to the upper boundary of the sticky air layer and the vertical sides of the model and the base is kept fixed.
There is no temperature difference, and the geodynamic evolution is guided solely by compositional density differences.

.. _crameri_setup:

.. figure:: figs/crameri-et-al-2012-case-2-setup.png
	:width: 90%
	:align: center
	:alt: Crameri case 2 model setup

	*Case 2* model setup to evaluate the sticky air method. Extracted from Crameri et al. (2012) :cite:`crameri2012`.

From the results of this experiment reproduced in MANDYOC we obtain the maximum topography with time, similar to Fig. 6a of Crameri et al. (2012) :cite:`crameri2012`, presented in :numref:`maximum_topography`.
The models used for comparison are: UNDERWORLD :cite:`moresi2003`, STAGYY :cite:`tackley1993` and I2VIS :cite:`gerya2003`.

.. _maximum_topography:

.. figure:: figs/crameri-et-al-2012-case-2-comparison.png
   :width: 100%
   :align: center
   :alt: Comparison of MANDYOC results

   Comparison of the maximum topography with time for the *Case 2* (:numref:`crameri_setup`) model setup from Crameri et al. (2012) :cite:`crameri2012`.
