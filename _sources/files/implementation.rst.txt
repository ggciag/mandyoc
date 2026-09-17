.. _numericalimplementation:

Numerical implementation
========================

The following sections show how the numerical methods were implemented to solve the equations of conservation of mass, momentum and energy that were presented in the :ref:`basictheory` section.

.. _massmomentumimplementation:

Mass and momentum equations
---------------------------

For a solution domain :math:`\Omega`, finding the flow velocity :math:`u_i=g_i+v_i` and pressure :math:`P` states the Galerkin weak formulation for Stokes' flow, where :math:`g_i` is a specified boundary velocity, :math:`v_i` belongs to a set of functions :math:`\mathcal{V}` in which every function is equal to zero on the boundary where the :math:`ith` components of the velocity is :math:`g_i`, and :math:`P` belongs to a set of functions :math:`\mathcal{P}` such that the equations of conservation of mass and momentum can be written as follows :cite:`zhong2007`:

.. math::
    :label: weak-formulation-1

    \int_{\Omega}{w_{i,j}\sigma_{ij}d\Omega} -
    \int_{\Omega}{qu_{i,i}d\Omega} =
    \int_{\Omega}{w_i f_id\Omega} +
    \sum_{i=1}^{n_{sd}}{\int_{\Gamma_{h_i}}{w_i h_i d\Gamma}}

where :math:`w_i \in \mathcal{V}` and :math:`q \in \mathcal{P}` are weighting functions and the boundary conditions are defined:

.. math::
    :label: boundary-conditions

    u_i = g_i \text{ on } \Gamma_{g_i}, \sigma_{ij}n_j = h_i \text{ on } \Gamma_{h_i}

where :math:`\Gamma_{h_i}` is the boundary where the :math:`ith` components of the forces are set to be :math:`h_i`, :math:`n_j` is the normal vector at the boundary :math:`\Gamma_{h_i}`, and :math:`f_i=g\rho_0 \alpha \delta_{i3}` :cite:`hughes2000`. The weak formulation can be re-written as below:

.. math::
    :label: weak-formulation-2

    \int_{\Omega}{w_{i,j}c_{ijkl}v_{k,l}d\Omega} -
    \int_{\Omega}{qv_{i,i}d\Omega} -
    \int_{\Omega}{w_{i,i}Pd\Omega} = \\
    \int_{\Omega}{w_{i}f_{i}d\Omega} +
    \sum_{i=1}^{n_{sd}}{\int_{\Gamma_{h_i}}{w_{i}h_{i}d\Gamma}} - 
    \int_{\Omega}{w_{i,j}c_{ijkl}g_{k,l}d\Omega}
    
where :math:`c_{ijkl}=\eta(\delta_{ik}\delta_{jl}+\delta_{il}\delta_{jk})` is obtained from the stress tensor equation (see :eq:`stress-tensor` in the :doc:`basic theory section<theory>`).

The velocity field, the pressure field and the weighting functions shape functions that interpolate the grid points at every node are given:

.. math::
    :label: v-shape-function

    \mathbf{v} = 
    v_i \mathbf{e}_i = 
    \sum_{A\in \Omega^{v}-\Gamma^{v}_{g_i}}{N_A v_{iA}\mathbf{e}_i}

.. math::
    :label: w-shape-function

    \mathbf{w} =
    w_i \mathbf{e}_i =
    \sum_{A\in \Omega^{v}-\Gamma^{v}_{g_i}}{N_A w_{iA}\mathbf{e}_i}

.. math::
    :label: g-shape-function

    \mathbf{g} =
    \sum_{A\in \Omega^{v}-\Gamma^{v}_{g_i}}{N_A g_{iA}\mathbf{e}_i}

.. math::
    :label: p-shape-functionn

    P = \sum_{B\in \Omega^{p}}{M_B P_B}

.. math::
    :label: q-shape-function

    q = \sum_{B\in \Omega^{p}}{M_B q_B}

where :math:`N_A` is the shape function for the velocity at node A, :math:`M_B` is the shape function for the pressure at node B, :math:`\Omega^{v}` is the velocity nodes set, :math:`\Omega^{p}` is the pressure nodes set and :math:`\Gamma^{g}_{g_i}` is the velocity nodes set along the boundary :math:`\Gamma_{g_i}`. *Mandyoc* defines :math:`\Omega^{p}` at the center of each element, while :math:`\Omega^{v}` is defined at every element vertex. This avoids spurious flow solutions and numerical instabilities :cite:`zhong2007` and it keeps the velocity shape functions one order higher than the pressure shape functions, a common strategy used in finite element modeling of incompressible media :cite:`hughes2000`.

From the shape functions above and the Galerkin weak formulation (:eq:`weak-formulation-2`), the following expression can be obtained:

.. math::
    :label: implication-1

    \sum_{B\in\Omega^{v}-\Gamma^{v}_{g_j}}{\Big( \mathbf{e}^{T}_{i} \int_{\Omega}{B^T_A D B_B d\Omega \mathbf{e}_j v_{jB}} \Big)} - 
    \sum_{B\in \Omega^p}{\Big( \mathbf{e}_i \int_{\Omega}{N_{A,i} M_B d\Omega P_{B}} \Big)} = \\
    \int_{\Omega}{N_A \mathbf{e}_i f_i d\Omega} +
    \sum_{i=1}^{n_{sd}}{\int_{\Gamma_{h_i}}{N_A \mathbf{e}_i h_i d\Gamma}} -
    \sum_{B\in \Gamma^{v}_{g_j}}{\Big( \mathbf{e}^{T}_{i} \int_{\Omega}{B^{T}_{A} D B^{T}_{B} d\Omega \mathbf{e}_j g_{jB}} \Big)}
    
and also:

.. math::
    :label: implication-2

    \sum_{B\in \Omega^{v} - \Gamma^{v}_{g_j}}{\int_{\Omega}{M_A N_{Bj} d\Omega \mathbf{e}_j v_{jB}}} = 0

The matrix representation of these two equations above can be presented:

.. math::
    :label: matrix-representation

    \begin{bmatrix}
        K & G \\
        G^T & 0
    \end{bmatrix}
    \left\{
        \begin{array}{c}
            V \\
            P
        \end{array}
    \right\} = 
    \left\{
        \begin{array}{c}
            F \\
            0
        \end{array}
    \right\} 

where :math:`V` is the vector of velocity values at :math:`\Omega^v`, :math:`P` is the vector of pressure values at :math:`\Omega^p`, :math:`F` is the resulting vector of the rigth-hand side of equations :eq:`implication-1` or :eq:`implication-2`, :math:`K` is the stiffness matrix, :math:`G` is the discrete gradient operator, and :math:`G^T` is the discrete divergence operator :cite:`zhong2007`. :math:`K`, :math:`G` and :math:`G^T` are derived from the first and second terms of :eq:`implication-1` and :eq:`implication-2`.

The matrix operator :math:`B` from :eq:`implication-1` contains the spatial derivatives of the shapes function :math:`N` and, for a 2-D plane, can be written as:

.. math::
    :label: B_A

    B_A = 
    \begin{bmatrix}
        N_{A,1} & 0 \\
        0 & N_{A,2} \\
        N_{A,2} & N_{A,1}
    \end{bmatrix}

Again from :eq:`implication-1` and for 2-D plane strain problems, the effective viscosity matrix :math:`D` can be written:

.. math::
    :label: D

    D = 
    \begin{bmatrix}
        2\eta & 0 & 0 \\
        0 & 2\eta & 0 \\
        0 & 0 & \eta
    \end{bmatrix}

The stiffness matrix :math:`K` and the gradient operator :math:`G` can be written:

.. math::
    :label: K

    K_{lm} = \mathbf{e}^T_{i} \int_{\Omega}{B^T_A D B_B d\Omega \mathbf{e}_j}

.. math::
    :label: G

    G_{lm} = \mathbf{e}_i \int_{\Omega}{N_A M_B d\Omega \mathbf{e}_j}

where the subscripts :math:`A` and :math:`B` are the global velocity node numbers, :math:`i` and :math:`j` are the degree of freedom per grid node, ranging from :math:`1` to :math:`n_{sd}`, :math:`l` and :math:`m` are global equation numbers for the velocity ranging from :math:`1` to :math:`n_v n_{sd}`, where :math:`n_{v}` is the number of velocity nodes in the grid.

.. _energyimplementation:

Energy equation
---------------

It is possible to represent the energy conservation equation (:eq:`energy-conservation`) as a finite element problem for a solution domain :math:`\Omega_V` as follow:

.. math::
    :label: fe-energy-conservation

    \mathbf{M} \mathbf{\dot{a}}_T + (\mathbf{K_a} + \mathbf{K_c})\mathbf{a}_T = \mathbf{F}

where :math:`\mathbf{M}`, :math:`\mathbf{K}_a`, :math:`\mathbf{K}_c` and :math:`\mathbf{F}` are written below:

.. math::
    :label: M

    \mathbf{M} = \int_{\Omega_V}{\mathbf{N}^T_V \rho_0 c_p \mathbf{N}_V d\Omega_V}

.. math::
    :label: Ka

    \mathbf{K}_a = \int_{\Omega_V}{\mathbf{N}^T_V  \rho_0 c_p \mathbf{v} \cdot \mathbf{B}_V d\Omega_v}

.. math::
    :label: Kc

    \mathbf{K}_c = \int_{\Omega_V}{\mathbf{B}^T_V \rho_0 c_p \mathbf{v} \cdot \mathbf{B}_v d\Omega_v}

.. math::
    :label: F

    \mathbf{F} = \int_{\Omega_V}{\mathbf{N}^T_V \Big( \frac{H}{c_p} - \frac{\alpha T g u_3}{c_p} \Big) d\Omega_v}

where :math:`\mathbf{N}_V` is a row vector of shape functions, :math:`\mathbf{a}_T` is a column vector of the unknown temperature parameters, :math:`\mathbf{\dot{a}}_T` is its time derivative, and :math:`\mathbf{B}_V \equiv \nabla \mathbf{N}_V`. 

.. note::
    The superscript :math:`T` represents the transpose of the matrix, while the non-superscript :math:`T` represents the temperature.

:math:`\mathbf{M}` and :math:`\mathbf{K}_c` are symmetric, but :math:`\mathbf{K}_a` is not. This asymmetry decreases the accuracy of the solution when advection is more dominant than conduction (:cite:`zienkiewicz2000`, chapter 2). To increase numerical accuracy and stability, *Mandyoc* uses the streamline upwind Petrov-Galerkin process to modify :math:`\mathbf{K}_a` to :math:`\mathbf{K}_a^*` :cite:`zienkiewicz2000,hughes1979,hughes1982`:

.. math::
    :label: Ka-star

    \mathbf{K}_a^* = \int_{\Omega_V}{\mathbf{N}^{*T}_V \rho_0 c_p \mathbf{v} \cdot \mathbf{B}_V d\Omega_V}

where :math:`N^{*}_{Vi}` is:

.. math::
    :label: N-star

    N^{*}_{Vi} = N_{Vi} + \frac{\alpha_{opt} h^e \mathbf{v} \cdot \nabla N_{Vi}}{2|\mathbf{v}|} 

where :math:`h^e` is the characteristic element size in the advection velocity direction (:math:`\mathbf{v}`), and :math:`\alpha_{opt}` is:

.. math::
    :label: alpha-opt

    \alpha_{opt} = \coth{Pe} - \frac{1}{Pe} \text{, where } Pe = \frac{|\mathbf{v}|h^e}{2 \kappa \rho_0 c_p}

The time discretization is made by an implicit scheme :cite:`braun2003`:

.. math::
    :label: time-discretization

    \frac{\mathbf{a}_T(t+\Delta t) - \mathbf{a}_T(t)}{\Delta t} = 
    \theta \mathbf{\dot{a}}_T(t+\Delta t)+(1-\theta)\mathbf{a}_T(t)

where :math:`\theta=0.5` is a weighting parameter. Multiplying both sides of :eq:`time-discretization` by :math:`\mathbf{M}(t+\Delta t)` and assuming :math:`\mathbf{M}(t+\Delta t) \approx \mathbf{M}(t)`:

.. math::
    :label: time-discretization-2

    \mathbf{M}(t+\Delta t) \frac{\mathbf{a}_T(t+\Delta t) - \mathbf{a}_T(t)}{\Delta t} = 
    \theta [\mathbf{F}(t+\Delta t) - \mathbf{K}_T(t+\Delta t) \mathbf{a}_T(t+\Delta t)] + \\
    (1-\theta)[\mathbf{F}(t)-\mathbf{K}_T(t)\mathbf{a}_T(t)]

where :math:`\mathbf{K}_T=\mathbf{K}^*_a+\mathbf{K}_c`. 

Rearranging :eq:`time-discretization-2` allows to rewrite it in a numerical form to solve the energy equation.

.. math::
    :label: numerical-form-energy-equation

    [\mathbf{M}(t+\Delta t) + \Delta t \theta \mathbf{K}_T(t+\Delta t)]\mathbf{a}_T(t)(t+\Delta t) = \\
    [\mathbf{M}(t+\Delta t)- \Delta t(1-\theta)\mathbf{K}_T(t)]\mathbf{a}_T(t) + \\
    \Delta t[\theta \mathbf{F}(t+\Delta t) + (1-\theta)\mathbf{F}(t)]
    
Free surface
------------

*Mandyoc* uses the *Free Surface Stabilization Algorithm* :cite:`kaus2010` to modify the Stokes equation and avoid numerical instabilities that can occur on the surface of the model. The up-and-down oscillations around the steady state ("sloshing instability" or "drunken sailer effect") would require small time steps, which would increase running time by unfeasible amounts.

The modification is done on the stiffness matrix :math:`K_e`, such that :math:`\tilde{K}_e=K_e+L_e`, where the correction :math:`L_e` is evaluated at the boundary :math:`\Gamma_e` of each finite element. The correction is given by the :eq:`Le` below:

.. math::
    :label: Le

    L_e = \int_{\Gamma_e}{\mathbf{N} \Theta \Delta\rho \Delta t \mathbf{g} \mathbf{n} d\Gamma}

where :math:`\mathbf{N}` is the element shape function, :math:`0 \leq \Theta \leq 1` is a wight factor of the correction term, :math:`\Delta \rho` is the density contrast between the two mediums, :math:`\Delta t` is the numerical integration time step, :math:`\mathbf{g}` is the gravity acceleration vector, and :math:`\mathbf{n}` is the normal vector to the element.

.. _rheologysection:

Rheology
--------

Considering a visco-plastic model, the effective viscosity :math:`\eta` follows the formulation described by Moresi and Solomatov (1998) :cite:`moresi1998`, which combines plastic deformation and viscous deformation:

.. math::
    :label: effective-eta
    
    \eta 
    = \min{(\eta_{plas},\eta_{visc})}
    = \min{\bigg(\frac{\tau_{yield}}{2\dot{\varepsilon}_{II}},\eta_{visc}\bigg)}

where :math:`\tau_{yield}` is the rupture tension and :math:`\dot{\varepsilon}_{II}=(\dot{\varepsilon}_{ij}'\dot{\varepsilon}_{ij}'/2)^{1/2}` is the second invariant of the deviatoric strain rate tensor, and :math:`\eta_{plas}` and :math:`\eta_{visc}` are the plastic and ductile viscosities.

Plastic deformation
*******************

The plastic deformation can be calculated using the Byerlee Law :cite:`byerlee1968` to compute :math:`\tau_{yield}` and :math:`\eta_{plas}`. :eq:`byerlee-law` shows the relationship implemented in *Mandyoc*.

.. math::
    :label: byerlee-law

    \tau_{yield} = c_{0}+\mu \rho g z

where :math:`c_{0}` in the internal cohesion of the rock, :math:`\mu` is the friction coefficient, :math:`\rho` is the density and :math:`z` is the depth.

Alternatively, the user can choose to use the the Druker-Prager criterion :cite:`drucker-prager1952`, which is presented in :eq:`drucker-prager`.

.. math::
    :label: drucker-prager

    \tau_{yield} = c_0 \cos{\varphi} + P \sin{\varphi}

where :math:`\varphi` is the internal angle of friction.

Viscous deformation
*******************

*Mandyoc* contains several rheology models that the user can choose for viscous deformation. Among them, two will be discussed here. 

The ductile rheology can be simulated using the Frank-Kamenetskii approximation, following the formulation described by Solomatov and Moresi (2000) :cite:`solomatov2000`, where the viscosity is a function of the temperature :math:`T` as in the equation below. Such formulation was also used by Sacek (2017) :cite:`sacek2017` in an earlier *Mandyoc* version.

.. math::
    :label: frank-kamenetskii

    \eta_{visc}(T) = C \eta_r b^* \exp{(-\gamma T)}

where :math:`\eta_r` is the reference viscosity, :math:`C` is a compositional factor to scale the effective viscosity, and :math:`b^*` and :math:`\gamma = E_a / RT^2_b` are constants, which in turn, :math:`E_a` is the activation energy, :math:`R` is the gas constant, and :math:`T_b` is the basal temperature.

Additionally, the rheology can also be considered to follow a power law, as a function of the temperature :math:`T`, compositional factor :math:`C`, pressure :math:`P` and strain rate :math:`\varepsilon` as follows:

.. math::
    :label: power-law

    \eta_{visc} = C A^{\frac{-1}{n}} \dot{\varepsilon}^{\frac{1-n}{n}} \exp{\frac{Q+V P}{nRT}}

where :math:`A` is a pre-exponential scale factor, :math:`n` is the power law exponent, :math:`\dot{\varepsilon}` is the square root of the second invariant of the strain rate tensor, :math:`Q` is the activation energy, and :math:`V` is the activation volume. The values of :math:`A`, :math:`n`, :math:`Q`, and :math:`V` are measured under laboratory conditions :cite:`karato1993,gleason1995`.


Non-linear iterations
*********************

When the non-linear option is chosen by the user, the effective viscosity is dependent on the velocity field, which is iteratively updated following the algorithm described by Thieulot (2014) :cite:`thieulot2014`. In this algorithm, the velocity and effective viscosity field are iteratively updated until the following convergence criterion is satisfied:

.. math::
    :label: non-linear

    \chi_f = 1 - \frac{\langle (f^i-\langle f^i\rangle) \cdot (f^{i+1}-\langle f^{i+1}\rangle) \rangle}{|f^i-\langle f^i\rangle| \ |f^{i+1}-\langle f^{i+1}\rangle|} \le tol 

where :math:`f` represents an array with all the nodal values of the velocity components, :math:`tol` is a tolerance factor, and :math:`\langle f\rangle` is the mean value of :math:`f`. The superscript :math:`i` and :math:`i+1` indicate two consecutive iterations in the same time step.  In each iteration, the momentum and mass equations are calculated with an updated effective viscosity field.





