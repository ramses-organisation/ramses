
1. Architecture of the code
===========================

.. Admonition:: **Topics**

   basic code structure, introduction of AMR step

.. contents::

1.1 The program RAMSES
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

1.2 An overview of amr_step
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
