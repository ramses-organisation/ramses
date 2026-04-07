3. Gravity
==========

.. contents::


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

