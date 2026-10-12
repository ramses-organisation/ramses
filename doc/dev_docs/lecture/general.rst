
Structure of the code
==============================

.. contents::

1. The program RAMSES
----------------------

The root of the program RAMSES is found in *amr/ramses.f90*:

.. code:: fortran

   program ramses
    call read_params      ! Read run parameters
    call adaptive_loop    ! Start time integration
   end program ramses

First, the routine ``read_params`` will load the parameters from the
namelist that was given as input by the user. Then, the routine
``adaptive_loop`` is called. It is found in *amr/adaptive_loop.f90* and structured as follows:

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

The simulation starts with the initialization: arrays are allocated and set to
appropriate initial values, initial conditions are calculated or read
from file,… After that, the main time loop is started, which will evolve the
simulation in time.

The core of RAMSES is the recursive routine ``amr_step`` found in the
file ``amr/amr_step.f90``. In this routine, all individual physics
components are called in a specific order.

.. admonition:: Exercise

   Look into the file ``amr/amr_step.f90``. Can you make a
   list of which physical processes are modelled in ramses? Take a few
   minutes to explore the different directories of the code. Can you
   find out where the code of each physical process is?

2. An overview of ``amr_step``
---------------------------

A simplified schematic version of the core routine ``amr_step`` shows
the structure (see ``amr/amr_step.f90``):

.. code:: fortran

   recursive subroutine amr_step(ilevel,icount)

      call refine
      call load_balance

      ... ! Some sink and particle stuff

      if(time_to_output) call dump_all

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

Things are done in a specific order. The reasons for this will become
more clear over the course of these lectures.


3. Recursivity of ``amr_step``
------------------------------

.. code:: fortran

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

.. admonition:: **Exercise**

   Write down the calls to ``amr_step`` assuming there are
   3 refinement levels.
      A) First assume there is no subcycling (nsubcycle(ilevel)==1 for all levels).
      B) Now assume you have subcycling for all levels.

   .. admonition:: **Solution**
      :class: dropdown

      If we split the
      computations done in amr_step into two parts: stuff done before the
      recursive call to amr_step and stuff done after

      .. code:: fortran

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

4. Time stepping
----------------

RAMSES enables adaptive time stepping where each AMR level evolves with
individual timesteps. Though, the following rule always applies:

:math:`\Delta t^{\ell}=\Delta t^{\ell+1}_1+\Delta t^{\ell+1}_2`

An example of time stepping with two levels in the figure below
(Credits: Romain Teyssier)

|image3|

Level 2 is updated first with first with a time step of size
:math:`\Delta t^{\ell+1}_1` and second with
:math:`\Delta t^{\ell+1}_2`.The coarse level :math:`\ell=1` is frozen
during fine level solves (one order of accuracy down !). The fine flux
are averaged in time at coarse fine boundaries. Then level :math:`\ell`
is updated.

:math:`\bf{F}^{n+1/2,\ell}_{i+1/2,j}=\frac{1}{\Delta t_1^{\ell+1}+\Delta t_2^{\ell+1}}\left( \Delta t_1^{\ell+1}\frac{\bf{F}^{n+1/4,\ell+1}_{i+1/2,j-1/4}+\bf{F}^{n+1/4,\ell+1}_{i+1/2,j+1/4}}{2} + \Delta t_2^{\ell+1}\frac{\bf{F}^{n+3/4,\ell+1}_{i+1/2,j-1/4}+\bf{F}^{n+3/4,\ell+1}_{i+1/2,j+1/4}}{2}    \right)`

The timestep is computed in ``pm/newdt_fine.f90``. For more information,
see Section 2.4 in the RAMSES paper (Teyssier 2002).

.. |image3| image:: https://codimd.math.cnrs.fr/uploads/upload_6308ef0625fe29ca63f27341366485d5.png
