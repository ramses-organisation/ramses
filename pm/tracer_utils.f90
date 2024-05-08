! Hooks for RAMSES
module tracer_utils
   use amr_parameters      ! nx, ny, nz, dp, ngridmax, nvector, …
   use amr_commons         ! ncoarse, father, xg, son, myid, cpu_map, cpu_map2, ncpu, nbor
   use pm_commons          ! xp, tp, idp, levelp, headp, mp, localseed, numbp, nextp, move_flag
   use random, only        : ranf
   use hydro_commons, only : uold, if1, if2, jf1, jf2, kf1, kf2, nvar
#ifndef WITHOUTMPI
   use mpi_mod
#endif
   implicit none

   real(kind=dp), dimension(1:3) :: skip_loc
   real(kind=dp) :: scale
   integer :: nx_loc

   logical :: ddebug = .false.
   ! logical :: ddebug = .true.
   integer :: dbpart = -1 !105

   real(dp), allocatable, dimension(:) :: proba_yield  ! Working array

   contains

   subroutine initialize_skip_loc
      logical, save :: firstCall = .true.

      if (firstCall) then
         skip_loc=(/0.0d0, 0.0d0, 0.0d0/)
         if(ndim>0) skip_loc(1) = dble(icoarse_min)
         if(ndim>1) skip_loc(2) = dble(jcoarse_min)
         if(ndim>2) skip_loc(3) = dble(kcoarse_min)

         nx_loc=(icoarse_max-icoarse_min+1)
         scale = boxlen/dble(nx_loc)

         firstCall = .false.
      end if

   end subroutine initialize_skip_loc

   !---------------------------------------------------------------------------------
   ! Hooks before and after grid creation
   !---------------------------------------------------------------------------------
   subroutine pre_kill_grid_hook(ind_cell, ilevel, nn, ibound, boundary_region)

      integer, intent(in) :: nn, ilevel, ibound
      logical, intent(in) :: boundary_region
      integer, intent(in), dimension(1:nvector) :: ind_cell

      !######################
      ! Customize here
      !######################

      integer :: ipart, igrid, i, j, dim
      real(dp) :: dx, prevxp(1:ndim)

      real(dp), dimension(1:ndim) :: tmp_xp

      call initialize_skip_loc

      ! print*, 'in pre_kill_grid_hook'
      dx = 0.5D0**ilevel

      if (MC_tracer) then
         ! For all particles, recenter them in the center of the grid that's being deleted
         ! (becoming a cell)
         do j = 1, nn
            igrid = son(ind_cell(j))
            ipart = headp(igrid)

            do i = 1, numbp(igrid)
               if (is_gas_tracer(typep(ipart))) then
                  prevxp(:) = xp(ipart, 1:ndim)
                  do dim = 1, ndim
                     xp(ipart, dim) = (xg(igrid, dim) - skip_loc(dim)) * scale
                  end do
                  ! Attach to the new cell
                  partp(ipart) = ind_cell(j)
                  tmp_xp(:) = xp(ipart, :)
                  call safe_move(tmp_xp, prevxp(1:ndim))
                  xp(ipart, :) = tmp_xp(:)
               end if
               ipart = nextp(ipart)
            end do
         end do
      end if
   end subroutine pre_kill_grid_hook

   subroutine post_kill_grid_hook(ind_cell, ilevel, nn, ibound, boundary_region)
      use amr_commons

      integer, intent(in) :: nn, ilevel, ibound
      logical, intent(in) :: boundary_region
      integer, intent(in), dimension(:) :: ind_cell

      !######################
      ! Customize here
      !######################

      call initialize_skip_loc

   end subroutine post_kill_grid_hook

   !---------------------------------------------------------------------------------
   ! Hooks before and after grid creation
   !---------------------------------------------------------------------------------
   subroutine pre_make_grid_fine_hook(ind_grid, ind_cell, ind, &
      ilevel, nn, ibound, boundary_region)
      use amr_commons
      integer, intent(in) :: nn, ind, ilevel, ibound
      logical, intent(in) :: boundary_region
      integer, dimension(1:nvector), intent(in) :: ind_grid, ind_cell

      !######################
      ! Customize here
      !######################
      call initialize_skip_loc

   end subroutine pre_make_grid_fine_hook

   subroutine post_make_grid_fine_hook(ind_grid, ind_cell, ind, &
      ilevel, nn, ibound, boundary_region)
      use amr_commons
      use hydro_commons
      use pm_commons
      use random
      integer, intent(in) :: nn, ind, ilevel, ibound
      logical, intent(in) :: boundary_region
      integer, dimension(1:nvector), intent(in) :: ind_grid, ind_cell

      !######################
      ! Customize here
      !######################
      real(dp) :: dx, dxcoarse

      real(dp) :: x(1:ndim), rand, mass(0:twotondim), prevxp(1:ndim)
      integer :: j, i, fgrid, igrid, ipart, ison, iskip, icell
      integer :: loc(1:3)
      logical :: ok

      real(dp), dimension(1:ndim) :: tmp_xp

      call initialize_skip_loc

      dx = 0.5D0**ilevel           ! dx of the new level
      dxcoarse = 0.5D0**(ilevel-1) ! dx of the previous level

      ! print*, 'in post_make_grid_fine_hook'
      if (MC_tracer) then
         ! Compute the expected location of particles relative to xg in dx units
         loc(3) = (ind-1) / 4
         loc(2) = (ind-1-loc(3)*4) / 2
         loc(1) = (ind-1-loc(3)*4-loc(2)*2)

         do j = 1, nn
            fgrid = ind_grid(j)
            igrid = son(ind_cell(j))
            ipart = headp(fgrid)

            ! Load masses
            do ison = 1, twotondim
               iskip = ncoarse + (ison-1)*ngridmax
               icell = iskip + igrid
               mass(ison) = uold(icell, 1)
            end do

            mass(0) = sum(mass(1:twotondim))

            do i = 1, numbp(fgrid)
               if (is_gas_tracer(typep(ipart))) then

                  ! Check whether the particle was in the refined cell
                  x(1:ndim) = cellCenter(ind, fgrid, dxcoarse)

                  ok = all(xp(ipart, 1:ndim) == x(1:ndim))
                  ! If the particle is in refined cell, spread it accordingly
                  if (ok) then

                     ! Pick a random direction
                     call ranf(tracer_seed, rand)

                     do ison = 1, twotondim
                        if (rand < mass(ison) / mass(0)) then
                           ! Move particle to center of new cells
                           prevxp(:) = xp(ipart, :)
                           xp(ipart, :) = cellCenter(ison, igrid, dx)

                           tmp_xp(:) = xp(ipart, :)
                           call safe_move(tmp_xp, prevxp(1:ndim))
                           xp(ipart, :) = tmp_xp(:)
                           partp(ipart) = igrid + iskip + ison
                           exit
                        else
                           rand = rand - mass(ison) / mass(0)
                        end if
                     end do
                  end if
               end if
               ipart = nextp(ipart)
            end do
         end do
      end if
   end subroutine post_make_grid_fine_hook

   ! Local function that returns the *cell* center given a grid and the
   ! index of the cell within the grid.
   function cellCenter(ison, ind_grid, dx) result (x)
      implicit none

      integer, intent(in)  :: ison, ind_grid
      real(dp), intent(in) :: dx
      real(dp), dimension(1:ndim) :: x

      real(dp), dimension(1:3) :: xc

      integer :: ixn, iyn, izn, dim
      ! Get the location of neighbor cell in its grid
      izn = (ison-1)/4
      iyn = (ison-1-4*izn)/2
      ixn = (ison-1-4*izn-2*iyn)

      ! Compute the expected location of the particles
      xc(1) = (dble(ixn)-0.5d0)*dx
      xc(2) = (dble(iyn)-0.5d0)*dx
      xc(3) = (dble(izn)-0.5d0)*dx

      do dim = 1, ndim
         x(dim) = (xg(ind_grid, dim) + xc(dim) - skip_loc(dim)) * scale
      end do
   end function cellCenter

   ! Move particles taking into account periodic boundary conditions
   subroutine safe_move(new_xp, prev_xp)
      real(dp), dimension(1:ndim), intent(inout) :: new_xp
      real(dp), dimension(1:ndim), intent(in)    :: prev_xp

      integer :: dim

      do dim = 1, ndim
         if (new_xp(dim) - prev_xp(dim) > 0.5d0*scale) then
            new_xp(dim) = new_xp(dim) - scale
         else if (new_xp(dim) - prev_xp(dim) < -0.5d0*scale) then
            new_xp(dim) = new_xp(dim) + scale
         end if
      end do

   end subroutine safe_move

   ! Helper function to determine the relative level of neighbouring cells
   subroutine relative_level(ind_ngrid, ind_ncell, direction, level)
      integer, intent(in) :: direction
      integer, intent(in), dimension(0:twondim) :: ind_ncell ! the neighboring cells
      integer, intent(in), dimension(0:twondim) :: ind_ngrid ! the neighboring parent grids
      integer, intent(out) :: level

      integer :: pos, ind_ngrid_ncell

      ! If the neighbor cell is a grid
      if (son(ind_ncell(direction)) > 0 ) then
         level = +1
      else
         ! get the grid containing the neighbor cell
         pos = (ind_ncell(direction)-ncoarse-1)/ngridmax+1
         ind_ngrid_ncell = ind_ncell(direction)-ncoarse-(pos-1)*ngridmax

         ! if the father grid and the neighbor cell's grid are neighors / the same
         if (ind_ngrid_ncell == ind_ngrid(direction) .or. ind_ngrid_ncell == ind_ngrid(0)) then
            level = 0
         else
            level = -1
         end if

      end if

   end subroutine relative_level

   ! Compute all four cells (in 3D) on the cell pointing in
   ! the given direction
   subroutine get_cells_on_face(direction, locs)
      implicit none

      integer, intent(in) :: direction
      integer, dimension(1:twotondim/2), intent(out) :: locs

      integer, save, dimension(1:6, 1:4) :: mapping
      logical, save :: firstCall = .true.

      integer, save :: twotondimo2 = twotondim/2

      if (firstCall) then
         mapping(1, 1:4) = (/1, 3, 5, 7/) ! left cells
         mapping(2, 1:4) = (/2, 4, 6, 8/) ! right cells
         mapping(3, 1:4) = (/1, 2, 5, 6/) ! top cells
         mapping(4, 1:4) = (/3, 4, 7, 8/) ! bottom cells
         mapping(5, 1:4) = (/1, 2, 3, 4/) ! front cells
         mapping(6, 1:4) = (/5, 6, 7, 8/) ! back cells

         firstCall = .false.
      end if

      locs(1:twotondimo2) = mapping(direction, 1:twotondimo2)

   end subroutine get_cells_on_face

   ! This routines decrease by one the move_flag of the MC tracer at
   ! level ilevel
   subroutine reset_tracer_move_flag(ilevel)
      use pm_commons
      use amr_commons

      implicit none

      integer, intent(in) :: ilevel

      integer :: ipart, jpart, next_part, jgrid, npart1, igrid

      ! Loop over grids
      igrid = headl(myid, ilevel)
      do jgrid = 1, numbl(myid, ilevel)
         npart1 = numbp(igrid)  ! Number of particles in the grid
         if (npart1 > 0) then
            ipart = headp(igrid)
            ! Loop over particles
            do jpart = 1, npart1
               ! Save next particle  <---- Very important !!!
               next_part = nextp(ipart)

               if (is_tracer(typep(ipart))) then
                  move_flag(ipart) = max(move_flag(ipart) - 1, 0)
               end if
               ipart = next_part  ! Go to next particle
            end do
         end if
         igrid = next(igrid)   ! Go to next grid
      end do

    end subroutine reset_tracer_move_flag

end module tracer_utils
