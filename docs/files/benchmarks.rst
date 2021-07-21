.. _benchmarks:

Benchmarks
==========

In order to test the accuracy of the MANDYOC code, its results can be compared to benchmark studies. The following subsections will present the procedures and results for some well established modeling problems. Such cases should provide information about the applicability and performance of the code as well as some of its limitations.

van Keken et al. (1997) :cite:`vankeken1997`
--------------------------------------------

The set of simulations proposed by van Keken et al. (1997) :cite:`vankeken1997` compares several methods of studying two dimensional thermochemical convection, where the Boussinesq approximation and infinite Prandtl number are used. 

The first simulation consists of two layers, where a buoyant thin layer is under a denser thicker package. The problem can be interpreted as a salt layer under a sediment package, and the interface between the layers is defined by the :eq:`interfacerayleigh` below.

.. math::
    :label: interfacerayleigh

    y=-0.8 \lambda_y + 0.02 \cos{\frac{\pi x }{\lambda_x}}

where :math:`\lambda_x` and :math:`\lambda_y` are the horizontal and vertical lengths of the simulated 2-D box, respectively. 

The simulation is carried out in a Cartesian box where the fluid is isothermal and Rayleigh-Taylor instability is expected for the proposed setup. The table below lists the parameters used to run this scenario.

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
      - :math:`\eta_0`
      - | (case 1a)
        | (case 1b)
        | (case 1c)

Lorem ipsum dolor sit amet.

.. figure:: figs/vankeken-snaps-1a.pdf
   :width: 80%
   :align: center
   :alt: Results

   Evolution of the isoviscous Rayleigh-Taylor instability for :math:`\eta_0/\eta_r=1.00`. The best result presented by van Keken et al. (1997) :cite:`vankeken1997` are on the left and the MANDYOC results are on the right. 

.. figure:: figs/vrms-1a.pdf
   :width: 100%
   :align: center
   :alt: Results

   Evolution of the :math:`v_{rms}` for :math:`\eta_0/\eta_r=1.00`. The van Keken et al. (1997) :cite:`vankeken1997` result is shown in black and the MANDYOC code result is shown in gray.

.. figure:: figs/vankeken-snaps-1b.pdf
   :width: 80%
   :align: center
   :alt: Results

   Evolution of the isoviscous Rayleigh-Taylor instability for :math:`\eta_0/\eta_r=0.10`. The best result presented by van Keken et al. (1997) :cite:`vankeken1997` are on the left and the MANDYOC results are on the right. 

.. figure:: figs/vrms-1b.pdf
   :width: 100%
   :align: center
   :alt: Results

   Evolution of the :math:`v_{rms}` for :math:`\eta_0/\eta_r=0.10`. The van Keken et al. (1997) :cite:`vankeken1997` result is shown in black and the MANDYOC code result is shown in gray.

.. figure:: figs/vankeken-snaps-1c.pdf
   :width: 80%
   :align: center
   :alt: Results

   Evolution of the isoviscous Rayleigh-Taylor instability for :math:`\eta_0/\eta_r=0.01`. The best result presented by van Keken et al. (1997) :cite:`vankeken1997` are on the left and the MANDYOC results are on the right. 

.. figure:: figs/vrms-1c.pdf
   :width: 100%
   :align: center
   :alt: Results

   Evolution of the :math:`v_{rms}` for :math:`\eta_0/\eta_r=0.01`. The van Keken et al. (1997) :cite:`vankeken1997` result is shown in black and the MANDYOC code result is shown in gray.
