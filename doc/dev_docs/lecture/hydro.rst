Hydrodynamics
================

.. contents::


1. The Euler equations of hydrodynamics
----------------------------------------

RAMSES uses conservative variables that are updated with the conservative equations (mass, momentum, total energy)

:math:`\frac{\partial \rho}{\partial t}  + \nabla\cdot \left[\rho\textbf{u} \right]  =  0 \\
\frac{\partial \rho\textbf{u}}{\partial t}  + \nabla \cdot\left[\rho \textbf{u}\otimes \textbf{u} + P \mathbb{I} \right] = 0\\
\frac{\partial E_\mathrm{T}}{\partial t} + \nabla\cdot \left[\textbf{u}\left( E_\mathrm{T} + P \right)\right] = 0,`

with
:math:`\rho` the density,
:math:`\textbf{u}` the velocity,
:math:`E_\mathrm{T}=e + 1/2\rho \textbf{u}^2` the total energy,
:math:`e` the gas thermal energy, and
:math:`P=(\gamma-1)e` the gas pressure.
:math:`\gamma` is the adiabatic index.

The system can be written in the canonical form

:math:`\frac{\partial \mathbb{U}}{\partial t} + \nabla\cdot\mathbb{F}(\mathbb{U}) = \mathbb{S},`

where :math:`\mathbb{U}` is the vector of conservative variables

:math:`\mathbb{U}=\left[ \begin{array}{c}
\rho \\
\rho\textbf{u} \\
E_\mathrm{T} \\
\end{array} \right],`

:math:`\mathbb{F}` is the flux, and
:math:`\mathbb{S}` the source terms.
:math:`\mathbb{U}` are the fundamentals variables in RAMSES which are stored in ``uold`` and ``unew`` (see below).
For practical reasons, RAMSES often switches to primitive variables

:math:`\mathbb{Q}=\left[ \begin{array}{c}
\rho \\
\textbf{u} \\
P \\
\end{array} \right].`

Indeed, it is easier to use the dual space of the
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

:math:`\frac{\partial E_\mathrm{NT}}{\partial t}  + \nabla\cdot \left[E_\mathrm{NT} \textbf{u} \right]  =  -P_\mathrm{NT}\nabla.\textbf{u}`

with :math:`P_\mathrm{NT}=(\gamma_\mathrm{rad}-1)E_\mathrm{NT}`.

Magnetic fields :math:`\textbf{B}` are stored at each cell faces. They
are updated using the induction equation in the ideal MHD limit.

:math:`\frac{\partial \textbf{B}}{\partial t}  - \nabla\times \left[\textbf{u} \times \textbf{B}\right]  = 0.`

2. Hydro variables in RAMSES
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

* the Euler variables: density, velocity and pressure
* the magnetic field vector (on the left cell face), when compiled with MHD
* an optional number of non-thermal energies (``NENER``)
* an optional number of passive scalars (``NMETALS`` and ``NPSCAL``)

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

.. code:: fortran

   do i=1,nhydro

We have

* for HYDRO: ``nhydro=neul``
* for MHD: ``nhydro=neul+3=8``

Additionally, the magnetic field on the right cell face is added at the
very end of the variable array. This means that when following the
evolution of magnetic fields, we store 6 additional variables. For
convenience, the code defines the parameter ``nvar_all`` which, in
addition to the independent ``nvar`` variables, includes the right
magentic field. So we have

* for HYDRO: ``nvar_all = nvar``
* for MHD: ``nvar_all = nvar+3``

To compute the cell-center magnetic field, for example for outputting,
we do

.. code:: fortran

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

.. Admonition:: **Summary**

   To access the Euler variables in ``uold`` and ``unew``,
   use the indices:

   -  density: ``1``
   -  total energy or pressure: ``neul``
   -  momentum or velocities (HYDRO case): ``2`` up to ``1+ndim``
   -  momentum or velocities (MHD case): ``2, 3, 4``

   To access the magnetic field:

   - on the left side of the cells: ``neul + 1``, ``neul+2``, ``neul+3``
   - on the right side of the cells: ``nvar+1``, ``nvar+2``, ``nvar+3``

   Additional variables are stored after the left magnetic field in the
   following order:

   - NENER: ``inener=nhydro+1`` up to ``nhydro+nener``
   - passive scalars: ``nhydro+nener+1`` up to ``nvar``

.. admonition:: Exercise

   How to add a field variable to ``uold`` and ``unew``?
   List all the things that need changing (allocation?, initialisation?)

   .. admonition:: **Solution**
      :class: dropdown

      * nvar in makefile
      * how to output metadata info, input (don’t forget conversion to primitive vars)
      * a recipe to convert conservative to primitive variable in the routine ``ctoprim``

3. The hydro step in ``amr_step``
----------------------------------

.. admonition:: Exercise

   Where are ‘uold’ and ‘unew’ altered? Ignore MPI communications for now.
   Because RAMSES makes use of common arrays
   which are globally defined and accessible by all parts of the code,
   it can be tricky to see which routines use these arrays as input or
   alter their state. Rewrite ``amr_step`` in pseudo-code, indicating
   where the arrays ``uold`` and ``unew`` are updated.

   .. admonition:: **Solution**
      :class: dropdown

      .. code:: fortran

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
