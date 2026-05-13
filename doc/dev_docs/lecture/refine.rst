2 Adaptive Mesh Refinement [WIP]
================================

   [name=Blaizot] It is somewhat unclear to me how the grid is
   initialised. Is this done from level 1 to levelmin(-1 ?) at the start
   ? If so, this is probably simple enough to provide the pseudo code as
   an example ? TC: I didn’t want to go into that here, but in section
   2.

In amr step
-----------

very beginning of amr step

.. code:: fortran

      ! Make new refinements and update boundaries
      if(ilevel==levelmin.or.icount>1)then
         do i=ilevel,nlevelmax
            ...
            call refine_fine(i)
         end do
       end if

2.1 Initialisation of the AMR grid
----------------------------------

-  ``make_grid_coarse``: create 1 new grid at level 1
-  ``kill_grid``: destroys grids at a specified level

Adding/removing grids
~~~~~~~~~~~~~~~~~~~~~

To creating and destroying of grids is handled by the routine
``refine``:

.. code:: fortran

   subroutine refine
     use amr_commons
     implicit none
     integer::ilevel

     call refine_coarse
     call build_comm(1)
     call make_virtual_fine_int(cpu_map(1),1)
     do ilevel=1,nlevelmax-1
        call refine_fine(ilevel)
        call build_comm(ilevel+1)
        call make_virtual_fine_int(cpu_map(1),ilevel+1)
     end do

   end subroutine refine

All the relevant routines can be found in the file
``amr/refine_utils.f90``: \* ``make_grid_coarse``: create 1 new grid at
level 1 \* ``kill_grid``: destroy a grid

Destroying grids
~~~~~~~~~~~~~~~~

Handled by the routine ``kill_grid``, which takes as input a list of
cells on a specified level L for which their child grid (containing
2\ :math:`^{ndim}` cells at level L+1) is to be destroyed.

This routine \* disconnects the grid from the linked list, \* adjusts
the grid variables to keep the tree up-to-date, \* resets *all* cell
variables for the destroyed cells at level L+1 (uold, unew, phi, …).

After this, the disconnected grids are added to the free memory linked
list.
