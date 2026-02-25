.. raw:: html

   <!--
   type: slide

   slideOptions:
     transition: concave
     theme: moon 
   -->

Lecture: General review
=======================

.. raw:: html

   <!--
   - Overview
       - overview of physics (hydro, MHD, RHD, gravity, …) in ramses, and organisation of the code in directories 
       - EXERCISE: Students explore the directory structure and figure out together what physical processes are modelled. 
   - General hydro + gravity
       - explain assuming uniform grid: start in a given state, compute evolution and update (→ introduce the order of things and explain why). Hydro, then gravity, then combined. 
       - EXERCISE: introduce basic data structure (unew, uold). Students rewrite in pseudo-code the basic solver steps discussed above using these variables. 
       - AMR = multiple levels, with refinement and de-refinement. Explain strategies, how and when it is done. 
       - `amr_step`: recursiveness (that comes from timestepping … ) to implement the level by level strategy recursively. Explain upload.Try to be convincing. 
       - EXERCISE: students write pseudo-code again with the recursive loop? Then compare it to the real `amr_step` / or extract pseudo code from messy `amr_step.f90`.

   Starting from the full overview of the physics found in the discussion, go through amr_step and explain why things are done in this order.


   [TOC]
   -->

1. Architecture of the code
===========================

.. container:: success

   **Topics:** basic structure, different directories, introduction of
   AMR step

1.1 The program RAMSES
----------------------

The root of the program RAMSES is found in *amr/ramses.f90*:

.. code:: fortran=

   program ramses
    call read_params      ! Read run parameters
    call adaptive_loop    ! Start time integration
   end program ramses

First, the routine ``read_params`` will load the parameters from the
namelist that was given as input by the user. Then, the routine
``adaptive_loop`` is called.

It is found in *amr/adaptive_loop.f90* and structured as follows:

.. code:: fortran

   subroutine adaptive_loop
     ! Initialize the simulation
     call init_amr
     call init_time
     if(hydro)call init_hydro
     ...
     ! Main time loop
     do
       ! Make new refinements level 1 to levelmin
       ...
       ! Call base level
       call amr_step(levelmin,1)
       ! Do some other stuff
       ...
       ! Print some info
       ...
     end do
   end subroutine adaptive_loop

First, the simulation is initialized: arrays are allocated and set to
appropriate initial values, initial conditions are calculated or read
from file, … Then the main time loop is started, which will evolve the
simulation in time.

The core of RAMSES is the recursive routine ``amr_step`` found in the
file ``amr/amr_step.f90``. In this routine, all individual physics
components are called in a specific order.

.. container:: info

   **Exercise:** Look into the file ``amr/amr_step.f90``. Can you make a
   list of which physical processes are modelled in ramses? Take a few
   minutes to explore the different directories of the code. Can you
   find out where the code of each physical process is?

1.2 An overview of amr_step
---------------------------

A simplified schematic version of the core routine ``amr_step`` shows
the structure (see ``amr/amr_step.f90``):

.. code:: fortran=

   recursive subroutine amr_step(ilevel,icount)

      call refine
      call load_balance
       
      ... ! Some sink and particle stuff
      
      if(time to output) call dump_all
      
      if (conditions are met)
         call kinetic_feedback             ! feedback from stars
         OR
         call make_stellar_from_sinks
         call make_sn_stellar              ! feedback from sinks
      end if
      
      if(poisson) call rho_fine   ! calc density field for Poisson source term
      
      ... ! Some particle stuff
      
      ! Gravity update: compute grav potential and acceleration
      if(poisson)then
         ...
         call phi_fine_cg(ilevel,icount) OR multigrid_fine(ilevel,icount)
         call force_fine(ilevel,icount)
         ...
     end if
      
     if(rt .and. rt_star/sink) call update_star/sink_RT_feedback(ilevel)
      
     call calc_turb_forcing(ilevel)   ! turbulence forcing
      
     call newdt_fine(ilevel)          ! Compute new time step
    
     if(hydro)call set_unew(ilevel)   ! set unew = uold
     if(rt)call rt_set_unew(ilevel) 
    
     ! --- Recursive call to amr_step ---
     ...
     !-----------------------------------

     if(conditions met) call thermal_feedback(ilevel)  ! feedback from stars

     if(sink.and.hydro) call grow_sink(ilevel,.false.)  ! sink accretion

     ! Hydro step: solve hydro and add source terms
     if((hydro).and.(.not.static_gas))then
        call godunov_fine(ilevel)
        ...
     endif
     
    ! Do RT/Chemistry step -> works on uold 
     if(rt .and. rt_advect) then
        call rt_step(ilevel)
     else
        call cooling_fine(ilevel)
     endif
     
     if(pic) call move_fine(ilevel)  ! Move particles

     if(conditions met)call star_formation(ilevel)

     ...  ! Update physical and virtual boundaries

     if(MHD) call diffusion  ! Magnetic diffusion step

     if(conditions met) call flag_fine  ! Compute refinement map

     ... ! particle stuff

     if(conditions met)call create_sink  ! Sink production
    
   end subroutine amr_step

Things are done is a specific order. The reasons for this will become
more clear in the course of these lectures.

2. Hydrodynamics
================

2.1 The Euler equations of hydrodynamics
----------------------------------------

RAMSES uses conservative variables that are updated with the
conservative equations (mass, momentum, total energy)
:math:`\frac{\partial \rho}{\partial t}  + \nabla\cdot \left[\rho\textbf{u} \right]  =  0 \\
\frac{\partial \rho\textbf{u}}{\partial t}  + \nabla \cdot\left[\rho \textbf{u}\otimes \textbf{u} + P \mathbb{I} \right] = 0\\
\frac{\partial E_\mathrm{T}}{\partial t} + \nabla\cdot \left[\textbf{u}\left( E_\mathrm{T} + P \right)\right] = 0,`
with :math:`\rho` the density, :math:`\textbf{u}` the velocity,
:math:`E_\mathrm{T}=e + 1/2\rho
\textbf{u}^2` the total energy, :math:`e` the gas thermal energy, and
:math:`P=(\gamma-1)e` the gas pressure. :math:`\gamma` is the adiabatic
index

The system can be written in the canonical form
:math:`\frac{\partial \mathbb{U}}{\partial t} + \nabla\cdot\mathbb{F}(\mathbb{U}) = \mathbb{S},`
where :math:`\mathbb{U}` is the vector of conservative variables
:math:`\mathbb{U}=\left[ \begin{array}{c}
\rho \\
\rho\textbf{u} \\
E_\mathrm{T} \\
\end{array} \right],` :math:`\mathbb{F}` is the flux, and
:math:`\mathbb{S}` the source terms. :math:`\mathbb{U}` are the
fundamentals variables in RAMSES which are stored in ``uold`` and
``unew`` (see below). For practical reasons, RAMSES often switches to
primitive variables :math:`\mathbb{Q}=\left[ \begin{array}{c}
\rho \\
\textbf{u} \\
P \\
\end{array} \right].` Indeed, it is easier to use the dual space of the
primitive variables in the integration of the hyperbolic solver. In
addition, when for instance feedback is activates, mass, momentum or
thermal energy are added, which requires to use primitive variables.

Additionnal variables can also be handled in RAMSES. These are
non-thermal energies :math:`E_\mathrm{NT}` (cosmic rays, radiative
energy) and passive scalars :math:`X` (metals, chemical species,
tracers, etc…). Passive scalars are variables that are passively
advected with the flow. These variables are integrated with the
following evolution equations

:math:`\frac{\partial \rho X}{\partial t}  + \nabla\cdot \left[\rho X\textbf{u} \right]  =  0`
and
:math:`\frac{\partial E_\mathrm{NT}}{\partial t}  + \nabla\cdot \left[E_\mathrm{NT} \textbf{u} \right]  =  -P_\mathrm{NT}\nabla.\textbf{u}`
with :math:`P_\mathrm{NT}=(\gamma_\mathrm{rad}-1)E_\mathrm{NT}`.

Magnetic fields :math:`\textbf{B}` are stored at each cell faces. They
are updated using the induction equation in the ideal MHD limit.

:math:`\frac{\partial \textbf{B}}{\partial t}  - \nabla\times \left[\textbf{u} \times \textbf{B}\right]  = 0.`

2.2 Hydro variables in RAMSES
-----------------------------

The arrays ``uold`` and ``unew``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In RAMSES, all hydro variables are stored together in the
two-dimensional arrays ``uold`` and ``unew``, which contain the value of
each variable in each cell of the AMR grid. They are defined in
``hydro_commons``

.. code:: fortran

   real(dp),allocatable,dimension(:,:)::uold,unew

and allocated in ``init_hydro``

.. code:: fortran

   allocate(uold(1:ncell,1:nvar))
   allocate(unew(1:ncell,1:nvar))

The arrays ``uold`` and ``unew`` are used for different purposes.
``uold`` stores the state that is used as input to integrate forward in
time.\ ``unew`` is a temporary array that is used to sum up flux
contributions in all directions, contributions for source terms, or
feedback. Note that ``unew`` can be desynchronized in time when AMR and
adaptive timesteps are used. In that case, ``unew`` stores the
contribution from the fine level to the coarse level. Importantly, this
means that ``unew`` should never be used to access the state variables,
or to communicate state variables of level :math:`\ell+1` or
:math:`\ell-1` between MPI domains.

Number of hydro variables ``nvar``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The total amount of independent variables stored on the AMR grid is set
at compile time in the Makefile by the parameter ``nvar``. It includes
\* the Euler variables: density, velocity and pressure \* the magnetic
field vector (on the left cell face), when compiled with MHD \* an
optional number of non-thermal energies (``NENER``) \* an optional
number of passive scalars (``NMETALS`` and ``NPSCAL``)

Remark that ``NMETALS`` and ``NPSCAL`` are not passed to the source
code. For the hydro-solver there is no distinction between these types
of variables and both are treated as passive scalars.

.. warning::

   Remark that ``uold`` and ``unew`` actually contain the
   **conservative** Euler variables :math:`\mathbb{U}`: density,
   momentum and total energy. For **primitive** variables (density,
   velocity and pressure) a conversion needs to be made first. See
   subroutine ``ctoprim``.

Accessing variables in ``uold`` and ``unew``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Several parameters are used to keep track of the number of different
types of variables and their indices in the arrays ``uold`` and
``unew``.

The number of **Euler variables** is indicated by ``neul`` in the code.
We always have density and pressure, but the amount of velocities to
keep track of depends on the number of dimensions of the simulation and
whether the code is compile with the HYDRO or MHD solver. For HYDRO, we
have ``neul = ndim+2`` with ``ndim`` the amount of spatial dimensions of
the simulation. For MHD simulations, we always need to keep track of the
three velocities and so ``neul = 5``. They are defined in
``hydro_parameters.f90``.

When using the MHD solver, we need additional variables to keep track of
the magnetic field in the three dimensions. Contrary to the hydro
variables which are defined in the center of each cell, the magnetic
field is defined in the cell faces. The magnetic field (on the left cell
face) for the three spatial directions is added after the Euler
variables. The number of Euler variables with addition of the magnetic
field is indicated by ``nhydro`` in the code. This makes it easy to loop
over them:

::

   do i=1,nhydro

We have - for HYDRO: ``nhydro=neul`` - for MHD: ``nhydro=neul+3=8``

Additionally, the magnetic field on the right cell face is added at the
very end of the variable array. This means that when following the
evolution of magnetic fields, we store 6 additional variables. For
convenience, the code defines the parameter ``nvar_all`` which, in
addition to the independent ``nvar`` variables, includes the right
magentic field. So we have - for HYDRO: ``nvar_all = nvar`` - for MHD:
``nvar_all = nvar+3``

To compute the cell-center magnetic field, for example for outputting,
we do

.. code:: fortran=

   Bx(i) = 0.5*(uold(i,neul+1) + uold(i,nvar+1))
   By(i) = 0.5*(uold(i,neul+2) + uold(i,nvar+2))
   Bz(i) = 0.5*(uold(i,neul+3) + uold(i,nvar+3))

RAMSES allow to include additional variables in the simulation: passive
scalars and non-thermal energies. Passive scalars can be used to keep
track of metals and star formation recipe variables. These specific
scalars are accessed through the indices ``imetal``, ``idelay``,
``ivirial1``, ``ivirial2``, ``ixion``, ``ichem``. These are used in
different star formation recipes. In the hydro-solver, these are evolved
as regular passive scalars.

.. warning::

   **Summary** To access the Euler variables in ``uold`` and ``unew``,
   use the indices:

   -  density: ``1``
   -  total energy or pressure: ``neul``
   -  momentum or velocities (HYDRO case): ``2`` up to ``1+ndim``
   -  momentum or velocities (MHD case): ``2, 3, 4``

   To access the magnetic field: - on the left side of the cells:
   ``neul + 1``, ``neul+2``, ``neul+3`` - on the right side of the
   cells: ``nvar+1``, ``nvar+2``, ``nvar+3``

   Additional variables are stored after the left magnetic field in the
   following order: - NENER: ``inener=nhydro+1`` up to ``nhydro+nener``
   - passive scalars: ``nhydro+nener+1`` up to ``nvar``

.. container:: info

   **Exercise**: How to add a field variable to ``uold`` and ``unew``?
   List all the things that need changing (allocation?, initialisation?)
   :::spoiler **Solution** \* nvar in makefile \* how to output metadata
   info, input (don’t forget conversion to primitive vars) \* a recipe
   to convert conservative to primitive variable in the routine
   ``ctoprim``

2.3 The hydro step in ``amr_step``
----------------------------------

.. container:: info

   **Exercise**: Where are ‘uold’ and ‘unew’ altered? (ignore MPI
   communications for now) Because RAMSES makes use of common arrays
   which are globally defined and accessible by all parts of the code,
   it can be tricky to see which routines use these arrays as input or
   alter their state. Rewrite ``amr_step`` in pseudo-code, indicating
   where the arrays ``uold`` and ``unew`` are updated. ::: spoiler
   **Solution**

   .. code:: fortran=

        ...
        ! Compute new time step
        call newdt_fine(ilevel)

        ! Set unew = uold
        if(hydro)call set_unew(ilevel)

        !---------------------------
        ! Recursive call to amr_step
        !---------------------------
        ...
        
        !-----------
        ! Hydro step
        !-----------
        if(hydro)then
           ! Hyperbolic solver - add flux to unew
           call godunov_fine(ilevel)

           ...
           ! Add gravity source terms to unew
           if(poisson) call add_gravity_source_terms(ilevel)

           ! Add non conservative pdV terms to unew
           ! for thermal and/or non-thermal energies
           if(pressure_fix.OR.nener>0) call add_pdv_source_terms(ilevel)

           ! Set uold = unew
           call set_uold(ilevel)

           ...

        endif
      ...

2.4 Solving hydro with the Finite Volume Method (FVM)
-----------------------------------------------------

Several numerical methods exist to solve the hydrodynamic equations on a
grid. RAMSES uses a finite volume method, more specifically an explicit
second-order predictor-corrector Godunov scheme (see further), to
integrate this conservative system of equations. This subset of methods
is particularly well-suited for solving hyperbolic systems of
conservation laws like the Euler equations. The method ensures
conservation of mass, momentum, and energy across cell boundaries by
construction.

In finite volume methods, the computational domain is discretized into
control volumes :math:`V_i` called cells. Remark that when using AMR,
these volumes can have different sizes. Writing the system equation in
integral form results in

:math:`\frac{d}{dt} \int_{V_i} \mathbb{U} \, dV + \int_{\partial V_i} \mathbb{F} \cdot \hat{n} dS, dA = \int_{V_i} \mathbb{S} \, dV`

When defining the cell-averaged conserved quantities as

:math:``
  \mathbb{U}_i(t) = \frac{1}{|V_i|} \int_{V_i} \mathbb{U}(\mathbf{x}, t) \, dV ,
  ``

and the numerical flux through face :math:`f` as :math:`\mathbb{F}_f`,
the semi-discrete update of the conservative variables becomes:

:math:``
\frac{d \mathbb{U}_i}{dt} = -\frac{1}{|V_i|} \sum_{f \in  \text{faces}} \mathbb{F}_f A_f + \mathbb{S}_i
``

with :math:`A_f` the area of face :math:`f` and :math:`\mathbb{S}_i` the
cell-averaged source term.

Numerical fluxes :math:`\mathbb{F}_f` are obtained by solving Riemann
problems at each interface between cells (see further), using the states
on the left :math:`\mathbb{U}_L` and right :math:`\mathbb{U}_R` of the
interface. For this, the interface states themselves first have to be
interpolated from the cell-centered values. To compute the fluxes RAMSES
used the MUSCL-Hancock scheme, which is a second-order Godunov method
(more on this further).

Discretizing further, the conservative variable update becomes (in 1D):

:math:``
\mathbb{U}_i^{n+1} = \mathbb{U}_i^{n} - \frac{\Delta t}{\Delta x} (\mathbb{F}_{i+1/2} -  \mathbb{F}_{i-1/2}) + \Delta t \, \mathbb{S}_i
``

where :math:`\mathbb{F}_{i\pm1/2}` are the fluxes going through opposing
faces of the cell.

In the code, the entry point for hydro solver is ``godunov_fine``. This
routine uses the standard RAMSES structure to gather ``nvector`` cells
to be send of to the core calculation routine ``godfine1``.

Inside ``godfine1``, the first step is to gather a stencil of 6 x 6 x 6
cells (in 3D) neighboring cells (this is needed to have second order
accurary in the slope calculation, see later). If there is AMR, the
neighboring cells may not be on the same refinement level and an
interpolation is performed. In that case, the coarse level is virtually
refined at the fine level, and the hydro variables are interpolated
(done by the routine ``interpol_fine``).

|image1|

Then, the subroutine ``unsplit`` is called, which calculates the fluxes.
The name comes from the implementation being ‘unsplit’, meaning that the
fluxes are computed simultaniously for all spatial directions (see
further). The fluxes are always computed at the left side of the cell.

Finally, at the end of ``godfine1``, ``unew`` is updated using the
fluxes from all directions (unsplit method).

How source terms are treated will be discussed further.

2.5 Unsplit MUSCL-Hancock scheme for computing the fluxes
---------------------------------------------------------

By default, RAMSES uses the MUSCL-Hancock scheme for computing the
numerical fluxes across cell faces. This is a predictor-corrector
extension of Godunov’s method that allows for second-order accuracy by
using \* piecewise linear reconstruction of the cell states, in contrast
to Godunov’s original piecewise constant \* a half step prediction for
time evolution.

The scheme is extended to multiple dimension in an unsplit fashion,
which is why the corresponding entry subroutine is named ``unsplit``.
This means that the fluxes in all directions are computed at the same
time, that is in one timestep. Transverse corrections are included. This
approach avoids splitting errors, better preserves symmetry and is
overall more accurate for multi-dimensional simulations.

The steps for obtaining the flux are as follows. \* Convert conservative
cell-centered variables to primitive variables (``ctoprim``) \*
Calculate the limited slopes (TVD) for the primitive variables that will
be used to reconstruct the state at the cell edges and evaluate space
derivative for the time evolution below (MUSCL part, ``uslope``) \*
Evolve (or ``trace``) the cell centered states forward in time for half
a time step and then project on cell faces (Hancock part). \* Solve the
Riemann problem at each interface using predicted left/right states to
obtain the fluxes :math:`F_{i+1/2}` ``cmpflxm`` (“compute flux minus”,
because only the left flux is calculated)

These fluxes are then used to update the conserved variables. Below, we
go into eqch step in more detail. You can visit the following page for a
more complete description of the MUSCL-Hancock scheme:
https://ammar-hakim.org/sj/hancock-muscl.html

Switch to primitive variables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, the conservative variables are converted into their primitive
form. This is done by the subroutine ``ctoprim``, which has as input the
cell-centered conservative variables :math:`\mathbb{U}` and as output
the corresponding primitive variables :math:`\mathbb{Q}`.

While not strictly required, it has become standard practice to use
primitive variables for higher-order Godunov schemes like the
MUSCL-Hancock scheme. This choice improves both accuracy and stability
of the method. Primitive variables tend to vary more smoothly,
especially across contact discontinuities and shocks, making them better
suited for linear reconstruction and slope limiting. Reconstructing in
conserved variables can introduce spurious oscillations or unphysical
states such as negative densities or pressures. By switching to
primitive variables during reconstruction and time prediction, the
scheme maintains physical realism while achieving second-order accuracy.

A subtlety arises when gravity is included in the simulation. In that
case, half a gravity predictor step is applied to the velocity in this
routine. See more further.

Compute the slopes
~~~~~~~~~~~~~~~~~~

In the MUSCL-Hancock scheme, the goal of the reconstruction step is to
approximate the solution inside each cell using a linear profile,
instead of a constant value. Given the cell-averaged primitive variables
:math:`Q_i`\ ​, the slope :math:`\Delta Q_i` is computed​ using the
cell-centered values of the neighboring cells. A slope limiter is
applied to control oscillations near discontinuities:

:math:`\Delta Q_i = \text{Limiter}(Q_{i+1}-Q_i, Q_i - Q_{i-1})`

|image2|

For example, the average slope without limiting would be simply

:math:`\Delta Q_i = \frac{(Q_{i+1} - Q_i) + (Q_i - Q_{i-1})}{2}`, which
can introduce new extrema in the reconstructed field. This configuration
is prone to oscillations. Total variation diminishing (TVD) limiters
prevent the apparition of these spurious oscillations.

In RAMSES, calculating the limited slope is done by the subroutine
``uslope``. Several slope limiters are implemented (minmod, monotonized
centered, etc…).

Predictor step: evolution of cell boundary values for half a time step
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The predictor step advances these states forward by half a timestep to
account for their local evolution (Source term part) and reconstruct the
state at the left and right interface of the cell. This increases
temporal accuracy to second order.

First, the cell centered primitive variables are advanced for hal a
timestep. The evolution is governed by the Euler equations, linearized
around the reconstructed primitive states. In practice, this is done by
applying the method of lines to estimate the time derivative and then
updating the interface states. In 1D

:math:`\rho^{n+1/2}_i=\rho^n_i - (u^n_i\Delta_x \rho^n_i+\rho^n_i\Delta_x u^n_i)\Delta t/\Delta x`
:math:`u^{n+1/2}_i=u^n_i - (u^n_i\Delta_x u^n_i+\Delta_x P^n_i/\rho^n_i)\Delta t/\Delta x`

In the case of multiple dimensions, transverse corrections are applied,
accounting for the coupling between spatial directions. These are
cross-derivative terms that arise when wave propagation in one direction
is affected by gradients in perpendicular directions, which is crucial
for preserving accuracy and symmetry in 2D and 3D flows.

In 2D
:math:`\rho^{n+1/2}_i=\rho^n_i - (u^n_i\Delta_x \rho^n_i+\rho^n_i\Delta_x u^n_i)\Delta t/\Delta x - (u^n_i\Delta_y \rho^n_i+\rho^n_i\Delta_y u^n_i)\Delta t/\Delta y`

:math:`u^{n+1/2}_i=u^n_i - (u^n_i\Delta_x u^n_i+\Delta_x P^n_i/\rho^n_i)\Delta t/\Delta x - (v^n_i\Delta_y u^n_i)\Delta t/\Delta y`

:math:`v^{n+1/2}_i=v^n_i - (u^n_i\Delta_x v^n_i)\Delta t/\Delta x - (v^n_i\Delta_y v^n_i +\Delta_y P^n_i/\rho^n_i)\Delta t/\Delta y`

Then, we define left and right interface states at the boundaries of
each cell:

:math:`Q_i^m = Q_i^{n+1/2} - \frac{1}{2} \Delta Q_i\\
Q_i^p = Q_i^{n+1/2} + \frac{1}{2} \Delta Q_i` These reconstructed values
represent the linearly extrapolated states at each interface, from
within cell ii. This reconstruction is applied in each spatial direction
independently.

The predictor step is perfomed in the subroutine ``trace``

Solving the Riemann problem to obtain the fluxes through cell interfaces
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The last step is to compute the flux through the cell faces given the
states to the left and right of the interface. This is a Riemann
problem, i.e. a discontinuity with a value on the left and another value
on the right, which can be solved numerically using Riemann solvers.
RAMSES users have the choice between several different Riemann solvers.
For more info on Riemann solvers, see elsewhere.

In the code, the subroutine ``cmpflxm`` is called for each spatial
direction. Inside ``cmpflxm``, the requested Riemann solver is called.
The output array ``flux`` is then updated with the 1D flux for the
requested spatial direction.

If AMR is active, the coarse level is updated during the fine level
update at fine-to-coarse boundary.

.. container:: info

   **Exercise**: 1/ Find in the code where the coarse level is virtually
   refined. Does it imply some interpolation? 2/ Find in the code where
   the coarse level update is done? What does it imply for conservative
   variable evolution? What about coarse-to-fine boundary? :::spoiler
   **Solution** 2/ in ‘umuscl.f90’. Quantities are conserved. The
   coarse-to-fine boundary is never considered, the flux being set to
   zero.

3. Gravity
==========

3.1 The Poisson equation for gravity
------------------------------------

Poisson’s equation for gravity is written as follows

:math:`\nabla^2 \phi = -4 \pi G \rho`

where :math:`\phi` is the gravitational potention and :math:`\rho` the
total density field, taking into account the gas and particles. The
gravitational force is given by

:math:`\bf{f}= -\nabla \phi`

The steps for solving gravity are as followed: \* determine the Poisson
source term \* compute the gravitational potential by solving the
Poisson equation \* calculate the gravitational force (or acceleration)
by taking the gradient of the potential \* apply the force to the gas
and particles.

We will address each of these steps in more detail below.

The routines to compute the gravitational potential and force can be
found in the subdirectory *poisson/*. Since source term and force
application deal with either gas or particles, the routines related to
these steps are either found in the directory *pm/* or *hydro/*.

3.2 Gravity variables in RAMSES
-------------------------------

The variables for the gravity module are defined in ``poisson_commons``:

::

    real(dp),allocatable,dimension(:)  ::phi,phi_old       ! Potential
    real(dp),allocatable,dimension(:)  ::rho               ! Density
    real(dp),allocatable,dimension(:,:)::f                 ! 3-force

and allocated in ``init_poisson``:

::

    allocate(rho (1:ncell))
    allocate(phi (1:ncell))
    allocate(phi_old (1:ncell))
    allocate(f   (1:ncell,1:3))
    rho=0; phi=0; f=0

They consist of \* the gravitational force ``f``, a vector with ``ndim``
dimensions, \* the gravitational potential ``phi`` and a copy of the old
state ``phi_old``, \* the total density distribution ``rho``, including
gas and particles.

3.3 Computing the Poisson source term
-------------------------------------

The first step is to compute the Poisson source term on the grid. We
need to take into account both the contribution from the gas and the
particles. This implies converting the collection of free moving
particle masses to a density field on the grid. There are several
numerical schemes to do this. In RAMSES, we use the cloud-in-cell (CIC)
scheme, which will be detailed in the chapter on Particles, section XX.
If the simulation has both gas and paticles, the contributions of both
will be added together, resulting in a total density field (stored in
``rho``) which will be used as the Poisson source term.

The routine responsible for this step is ``rho_fine``, which is called
in ``amr_step``:

.. code:: fortran=

   ! Compute poisson source term (i.e. the density field)
   call rho_fine(ilevel,icount)

3.4 Solving for the gravitational potential and force
-----------------------------------------------------

Once we have the density field, we can feed it to our Poisson solver of
choice, which can be either Multigrid (``multigrid_fine``) or conjugent
gradient (CG - ``phi_fine_cg``). These routines calculate the
gravitational potential and store it in ``phi``. The inner workings of
these solver is quite complex, so we will not go into the details here.

The gravitational force can then be determined by computing the gradient
of ``phi``. This is done by the routine ``force_fine``. The resulting
3-force is stored in the array ``f``.

We can find these two step in ``amr_step``:

.. code:: fortran=

   ! Compute gravitational potential
   if(ilevel>levelmin)then
      if(ilevel .ge. cg_levelmin) then
         call phi_fine_cg(ilevel,icount)
      else
         call multigrid_fine(ilevel,icount)
      end if
   else
      call multigrid_fine(levelmin,icount)
   end if
   !when there is no old potential...
   if (nstep==0)call save_phi_old(ilevel)

   ! Compute gravitational acceleration
   call force_fine(ilevel,icount)

3.5 Applying the gravitational force on the gas
-----------------------------------------------

When there is gravity, a source term :math:`\mathbb{S}` is added to the
Euler equations to account for the gravitational acceleration:

:math:`\frac{\partial \mathbb{U}}{\partial t} + \nabla\cdot\mathbb{F}(\mathbb{U}) = \mathbb{S}`.

In ramses, the gravitational source term is calculated as

:math:`\mathbb{S}^{n+\frac{1}{2}}_{i}=
\left[
\begin{array}{l}
0 \\[6pt]
\dfrac{\rho_i^{n}\nabla \phi_i^{n} + \rho_i^{n+1}\nabla \phi_i^{\,n+1}}{2} \\[10pt]
\dfrac{(\rho \textbf{u})_i^{n}\nabla \phi_i^{n} + (\rho \textbf{u})_i^{n+1}\nabla \phi_i^{\,n+1}}{2}
\end{array}
\right]
\;\cdot`

In the code, :math:`\mathbb{S}` is referred to as the gravity source
term (not to be confused with the Poisson source term). For the gas,
gravitional acceleration is taken into account using a Velvet scheme
(time centered).

.. container:: info

   **Exercise**: Where in ``amr_step`` is the gravitationnal
   acceleration source term integrated when you have hydro? Write the
   corresponding pseudo-code (assuming there are no particles).
   :::spoiler **Solution** The acceleration is added in four places, but
   with a subtile change of sign in one of the calls. Equation 13 in
   Teyssier (2002) is done with ``add_gravity_source_terms(ilevel)``
   (index n) and ``line 19`` for n+1.

   The parts in ``amr_step`` relevant for the gravity calculation in the
   case of hydro (ignoring particles) can be summarized as follows:

   .. code:: fortran=

      ! GRAVITY UPDATE
      if(poisson)then
         ! Save old potential for time-extrapolation at level boundaries
         call save_phi_old(ilevel)

         ! Compute poisson source term (i.e. the density field)
         call rho_fine(ilevel,icount)

         ! Remove gravity source term with half time step and old force (u+0.5*f*dt)
         if(hydro) call synchro_hydro_fine(ilevel,-0.5*dtnew(ilevel),1)

         ! Compute gravitational potential using multigrid of CG method
         call multigrid_fine(ilevel,icount) OR phi_fine_cg(ilevel,icount)

         ! Compute gravitational acceleration
         call force_fine(ilevel,icount)

         ! Add gravity source term with half time step and new force
         if(hydro) call synchro_hydro_fine(ilevel,+0.5*dtnew(ilevel),1)
      end if
      ...

      ! Compute new time step
      call newdt_fine(ilevel)

      ! Set unew equal to uold
      if(hydro)call set_unew(ilevel)

      ! RECURSIVE STEP TO AMR_STEP
      ...

      ! HYDRO STEP
      if(hydro)then
         ! Solve hydro
         ...
         ! Add gravity source terms with half a time step to unew
         if(poisson)call add_gravity_source_terms(ilevel)
         
         ! Set uold equal to unew
         call set_uold(ilevel)
         
         ! Add gravity source term with half time step and old force
         ! in order to complete the time step
         if(poisson)call synchro_hydro_fine(ilevel,+0.5*dtnew(ilevel),1)
         ...
      end if

   Remark that the routine ``synchro_hydro_fine()`` alters ``uold``,
   while ``add_gravity_source_terms()`` alters ``unew``.

   Half a timestep is added at the end of the global time step to
   syncrhonize all levels and to make outputs at the beginning of the
   next timestep. This contribution is then removed ater the dump.

3.5 Applying the gravitational force on the particles
-----------------------------------------------------

For the particles, the gravitational acceleration is taken into account
using a leap-frog integrator (see Section 2.2.5 in Teyssier (2002). The
particle positions and velocities are first updated by a predictor step
(“kick-drift”):

:math:`v_{n+1/2}=v_n + \frac{1}{2} f_n \Delta t`

:math:`x_{n+1}=x_n + v_{n+1/2} \Delta t`

where :math:`f` is the gravitational acceleration. This is then followed
by a corrector step (“kick”):

:math:`v_{n+1} = v_{n+1/2} + \frac{1}{2} f_{n+1} \Delta t`

Finding where these updates are performed can be a little tricky and
counter-intuitive. The predictor step is done by ``move_fine``, while
the corrector step is done by ``synchro_fine`` (it *synchronises* the
velocities to the current time). For the corrector step, the
gravitational acceleration at time :math:`t^{n+1}` is needed. For this
reason, it is postponed until the *next* time step, right after the new
gravitational force has been calculated using the updated particle
positions. In ``amr_step``, we thus find ``synchro_fine`` *before*
``move_fine``:

.. code:: fortran=

   ! Gravity update
   if(poisson)then
      ! Compute gravitational potential using multigrid of CG method
      call multigrid_fine(ilevel,icount) OR phi_fine_cg(ilevel,icount)

      ! Compute gravitational acceleration
      call force_fine(ilevel,icount)

      ! Synchronize remaining particles for gravity
      if(pic) call synchro_fine(ilevel)

   end if
   ...
   ! Compute new time step
   call newdt_fine(ilevel)

   ! RECURSIVE STEP TO AMR_STEP
   ...

   ! Move particles
   if(pic) call move_fine(ilevel) ! Only remaining particles  

Because the gravitational force is known on the grid, not at the
particle positions, we need to apply the inverse CIC scheme (see chapter
on Particles) when updating particle velocities. In summary, the force
acting on the particle will be interpolated from the cells with which
the particle “overlaps”. Once this is done, the particle velocities can
be updated using the time-step of the level on which the particle lives.

4. Detailed dive into ``amr_step``
==================================

Because RAMSES is an AMR code with the possiblitly for adaptive
timestepping, there are certain complexities in ``amr_step``.

recursivity
-----------

.. code:: fortran=

   recursive subroutine amr_step(ilevel,icount)
     ...
     ! do things in the beginning
     ...
     !---------------------------
     ! Recursive call to amr_step
     !---------------------------
     if(ilevel<nlevelmax)then
        if(numbtot(1,ilevel+1)>0)then  !if there is stuff in next level
           if(nsubcycle(ilevel)==2)then
              call amr_step(ilevel+1,1)
              call amr_step(ilevel+1,2)
           else
              call amr_step(ilevel+1,1)
           endif
        else
           ! Otherwise, update time and finer level time-step
           dtold(ilevel+1)=dtnew(ilevel)/dble(nsubcycle(ilevel))
           dtnew(ilevel+1)=dtnew(ilevel)/dble(nsubcycle(ilevel))
           call update_time(ilevel)
        end if
     else
        call update_time(ilevel)
     end if
     ...  
     ! do things at the end
     ...
   end subroutine amr_step

.. container:: info

   **Exercise** Write down the calls to ``amr_step`` assuming there are
   3 refinement levels. A) First assume there is no subcycling
   (nsubcycle(ilevel)==1 for all levels). B) Now assume you have
   subcycling for all levels. :::spoiler **Solution** If we split the
   computations done in amr_step into two parts: stuff done before the
   recursive call to amr_step and stuff done after

   .. code:: fortran=

      recursive subroutine amr_step(ilevel,icount)

         ! calc phi(ilevel), set unew(ilevel)=uold(ilevel), ..., calc dt(ilevel)
         stuff_before(ilevel)
         
         ! recursive call
         if(ilevel<nlevelmax)then
            if(nsubcycle(ilevel)==2)then
               call amr_step(ilevel+1,1)
               call amr_step(ilevel+1,2)
            else
               call amr_step(ilevel+1,1)
            endif
         else
           call update_time(ilevel)
         end if
         
         ! solve hydro, set uold(ilevel)=unew(ilevel), ...
         stuff_after(ilevel)

      end subroutine amr_step

   A) Without subcycling

   ::

      call amr_step(l-1)
              stuff_before(l-1)
          call amr_step(l)
              stuff_before(l)
              call amr_step(l+1)
                  stuff_before(l+1)
                  t = t+dt(l+1)
                  stuff_after(l+1)
              stuff_after(l)
          stuff_after(l-1)

   At the end we have advanced by dt(l+1)

   B) With subcycling

   ::

      call amr_step(l-1,1)
              stuff_before(l-1)
          call amr_step(l,1)
              stuff_before(l)
              call amr_step(l+1,1)
                  stuff_before(l+1)
                  t = t+dt(l+1)
                  stuff_after(l+1)
              call amr_step(l+1,2)
                  stuff_before(l+1)
                  t = t+dt(l+1)
                  stuff_after(l+1)
              stuff_after(l)
          call amr_step(l,2)
              stuff_before(l)
              call amr_step(l+1,1)
                  stuff_before(l+1)
                  t = t+dt(l+1)
                  stuff_after(l+1)
              call amr_step(l+1,2)
                  stuff_before(l+1)
                  t = t+dt(l+1)
                  stuff_after(l+1)
              stuff_after(l)
          stuff_after(l-1)

   At the end we have advanced by 4 level l+1 timesteps dt(l+1)

Time stepping
-------------

RAMSES enables adaptive time stepping where each AMR level evolves with
individual timesteps. Though, the following rule always applies:

:math:`\Delta t^{\ell}=\Delta t^{\ell+1}_1+\Delta t^{\ell+1}_2`

An example of time stepping with two levels in the figure below
(Credits: Romain Teyssier) |image3|

Level 2 is updated first with first with a time step of size
:math:`\Delta t^{\ell+1}_1` and second with
:math:`\Delta t^{\ell+1}_2`.The coarse level :math:`\ell=1` is frozen
during fine level solves (one order of accuracy down !). The fine flux
are averaged in time at coarse fine boundaries. Then level :math:`\ell`
is updated.

:math:`\bf{F}^{n+1/2,\ell}_{i+1/2,j}=\frac{1}{\Delta t_1^{\ell+1}+\Delta t_2^{\ell+1}}\left( \Delta t_1^{\ell+1}\frac{\bf{F}^{n+1/4,\ell+1}_{i+1/2,j-1/4}+\bf{F}^{n+1/4,\ell+1}_{i+1/2,j+1/4}}{2} + \Delta t_2^{\ell+1}\frac{\bf{F}^{n+3/4,\ell+1}_{i+1/2,j-1/4}+\bf{F}^{n+3/4,\ell+1}_{i+1/2,j+1/4}}{2}    \right)`

The timestep is computed in ``pm/newdt_fine.f90``. For more information,
see Section 2.4 in the RAMSES paper (Teyssier 2002).

.. |image1| image:: https://codimd.math.cnrs.fr/uploads/upload_67cf05e8e166e3e215c15c36b1f5c711.png
.. |image2| image:: https://codimd.math.cnrs.fr/uploads/upload_7f9476551d953186c06bafd5044a02c0.png
.. |image3| image:: https://codimd.math.cnrs.fr/uploads/upload_6308ef0625fe29ca63f27341366485d5.png
