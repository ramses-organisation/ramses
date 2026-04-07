4. Detailed dive into ``amr_step``
==================================

.. contents::

Because RAMSES is an AMR code with the possiblitly for adaptive
timestepping, there are certain complexities in ``amr_step``.

4.1 recursivity
---------------

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

4.2 Time stepping
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
