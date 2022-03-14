subroutine move_tracer_fine(ilevel)
   use amr_commons, only: headl, numbl, next, nvector, myid
   use pm_commons, only: typep, headp, numbp, nextp, move_flag
   use pm_commons, only: is_gas_tracer, is_star_tracer, is_cloud_tracer, part_t
   implicit none

   integer, intent(in) :: ilevel

   integer :: ig, ip, igrid, jgrid, ipart, jpart, next_part, npart1, local_counter
   integer, dimension(nvector) :: ind_grid, ind_part, ind_grid_part
   type(part_t) :: part_type
   ig=0
   ip=0
   ind_grid=0
   ind_part=0
   ind_grid_part=0

   igrid=headl(myid,ilevel)
   do jgrid=1,numbl(myid,ilevel)
      npart1=numbp(igrid)  ! Number of particles in the grid
      if(npart1>0)then
         ig=ig+1

         ind_grid(ig)=igrid
         ipart=headp(igrid)
         local_counter=0
         ! Loop over particles
         do jpart=1,npart1
            ! Save next particle  <---- Very important !!!
            next_part=nextp(ipart)
            if(ig==0)then
               ig=1
               ind_grid(ig)=igrid
            end if

            ! call debug_part(ipart, '@move_fine')
            part_type = typep(ipart)
            if (is_gas_tracer(part_type) .and. move_flag(ipart) == 0) then
               local_counter=local_counter+1
               ip=ip+1
               ind_part(ip)=ipart
               ind_grid_part(ip)=ig
               if(ip==nvector)then
                  call move_gas_tracer(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
                  local_counter=0
                  ip=0
                  ig=0
               end if
            end if

            ipart=next_part  ! Go to next particle
         end do
         ! End loop over particles

         ! If there was no particle in the grid, remove the grid from the buffer
         if (local_counter == 0 .and. ig>0) then
            ig=ig-1
         end if
      end if
      igrid=next(igrid)   ! Go to next grid
   end do
   ! End loop over grids
   if(ip>0) call move_gas_tracer(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel) ! MC Tracer

end subroutine move_tracer_fine

subroutine move_gas_tracer(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use poisson_commons
  use hydro_commons, only: fluxes
  use tracer_utils, only: safe_move, relative_level, get_cells_on_face
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !------------------------------------------------------------
  ! This routine moves the tracer following the fluxes
  ! This routine is called by move_fine.
  !------------------------------------------------------------
  integer::i,j,idim,nx_loc
  real(dp)::dx,scale
  ! Grid-based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
  ! Particle-based arrays
  real(dp),dimension(1:nvector,1:ndim),save::new_xp,old_xp
  real(dp),dimension(1:nvector,1:3),save::x
  real(dp),dimension(1:3)::skip_loc

  ! MC tracer
  real(dp),dimension(1:nvector,1:twondim,1:twotondim),save::flux
  real(dp),dimension(1:nvector,1:twondim),save::flux_part
  real(dp),dimension(1:nvector),save::rand1, rand2, rand3
  real(dp),dimension(1:nvector,1:ndim),save::proba_correction
  real(dp),dimension(1:twotondim/2)::neighborflux
  integer,dimension(1:twotondim/2)::ncell
  real(dp)::proba1, proba2, proba3
  integer::ison,ipart,iskip,dir,ndir,itmp
  integer,dimension(1:nvector, 0:twondim), save :: ind_ngrid
  integer,dimension(1:nvector, 1:twotondim), save :: ind_cell

  integer,dimension(0:twondim), save :: tmp_ind_ngrid2
  integer,dimension(0:twondim), save :: tmp_ind_ncell2
  real(dp), dimension(1:ndim) :: tmp_old_xp, tmp_new_xp

  integer,dimension(1:nvector) :: ind_parent_part, ison_part, ind_ggrid_part
  integer,dimension(1:nvector, 1:twotondim, 0:twondim) :: ind_ncell
  integer,dimension(1:nvector, 0:twondim) :: tmp_ncell
  integer, dimension(1:nvector) :: tmp_ind_cell
  integer, dimension(1:twotondim/2) :: tmp_ind_ncell
  integer, dimension(1:nvector, 1:twotondim, 1:twondim) :: rel_lvl
  integer, dimension(1:nvector) :: new_partp
  real(dp) :: rand, rtmp ! temporary real
  real(dp), dimension(1:nvector) :: outflux, factor
  logical, dimension(1:nvector) :: move

  logical :: ok

  integer :: ix, iy, iz

  ! Mesh spacing in that level
  dx = 0.5D0**ilevel
  nx_loc = (icoarse_max - icoarse_min + 1)
  skip_loc = (/0.0d0, 0.0d0, 0.0d0/)
  if (ndim > 0) skip_loc(1) = dble(icoarse_min)
  if (ndim > 1) skip_loc(2) = dble(jcoarse_min)
  if (ndim > 2) skip_loc(3) = dble(kcoarse_min)
  scale = boxlen / dble(nx_loc)

  !=======================================!
  ! Grid specific code                    !
  !=======================================!
  ! Get neighbor grids
  call getnborgrids(ind_grid, ind_ngrid, ng)

  ! Compute the index of the cells of each grid
  do ison = 1, twotondim
     iskip = ncoarse + (ison - 1)*ngridmax
     do j = 1, ng
        ind_cell(j, ison) = ind_grid(j) + iskip
     end do
  end do

  ! Get each cell's neighbour in all 6 directions
  do ison = 1, twotondim
     do j = 1, ng
        tmp_ind_cell(j) = ind_cell(j, ison)
     end do
     call getnborfather(tmp_ind_cell, tmp_ncell, ng, ilevel)
     do dir = 0, twondim
        do j = 1, ng
           ind_ncell(j, ison, dir) = tmp_ncell(j, dir)
        end do
     end do
  end do

  ! Loop over dimension
  do ison = 1, twotondim
     iskip = ncoarse + (ison-1)*ngridmax

     ! Compute level of cell 'ison' in direction 'dir'
     do dir = 1, twondim
        do j = 1, ng
           tmp_ind_ngrid2(0:twondim) = ind_ngrid(j, 0:twondim)
           tmp_ind_ncell2(0:twondim) = ind_ncell(j, ison, :)
           call relative_level( &
                tmp_ind_ngrid2, tmp_ind_ncell2, dir, &
                rel_lvl(j, ison, dir))
        end do
     end do

     ! Get the flux for the cells in the grid
     do dir = 1, twondim
        do j = 1, ng
           flux(j, dir, ison) = fluxes(ind_cell(j, ison), dir)
        end do
     end do
  end do

  !=======================================!
  ! Particle specific code                !
  !=======================================!
  ! Create global index of grid
  do ipart = 1, np
     ind_ggrid_part(ipart) = ind_grid(ind_grid_part(ipart))
  end do

  ! Generate random numbers for each particle
  do ipart = 1, np
     call ranf(tracer_seed, rand1(ipart))
     call ranf(tracer_seed, rand2(ipart))
     call ranf(tracer_seed, rand3(ipart))
  end do

  ! Movable particles have a flag == 0
  do ipart = 1, np
     move(ipart) = .true.
  end do

  !=======================================!
  ! Force particles to be attached        !
  !=======================================!
  do idim = 1, ndim
     do j = 1, ng
        x0(j, idim) = xg(ind_grid(j), idim)
     end do
  end do

  ! Compute the location of the particle relative to its grid
  do idim = 1, ndim
     do ipart = 1, np
        ! Get the location in the grid in dx unit from center
        x(ipart, idim) = xp(ind_part(ipart), idim) / scale + skip_loc(idim)
        x(ipart, idim) = x(ipart, idim) - x0(ind_grid_part(ipart), idim)
        x(ipart, idim) = x(ipart, idim) / dx
     end do
  end do

  ! Reset ison_part to 1
  do ipart = 1, np
     ison_part(ipart) = 1
  end do

  ! Move particles to cell centers
  do idim = 1, ndim
     do ipart = 1, np
        ! Particle in the center of the grid
        if (x(ipart, idim) == 0.0D0) then
           call ranf(tracer_seed, rand)

           ! Project the particle either to the right or the left
           if (rand < 0.5) then
              ! do nothing
              x(ipart, idim) = -0.5D0
           else
              x(ipart, idim) = +0.5D0
              ison_part(ipart) = ison_part(ipart) + 2**(idim-1)
           end if

        else if (x(ipart, idim) < 0.0D0) then
           x(ipart, idim) = -0.5D0
        else if (x(ipart, idim) > 0.0D0) then
           x(ipart, idim) = +0.5D0
           ison_part(ipart) = ison_part(ipart) + 2**(idim-1)
        end if
     end do
  end do

  ! Recompute the index of the parent cell
  do ipart = 1, np
     iskip = ncoarse + (ison_part(ipart) - 1)*ngridmax
     new_partp(ipart) = ind_ggrid_part(ipart) + iskip
     ind_parent_part(ipart) = new_partp(ipart)
  end do

  ! Get the flux for each gas particle
  do ipart = 1, np
     do dir = 1, twondim
        flux_part(ipart, dir) = flux(ind_grid_part(ipart), dir, ison_part(ipart))
     end do
  end do

  do ipart = 1, np
     factor(ipart) = 1.0D0
  end do

  !=======================================!
  ! Move tracers attached to grid         !
  !=======================================!
  ! Compute the outgoing fluxes for each particle in a cell + mass of cell
  outflux(1:np) = 0
  do dir = 1, twondim
     do ipart = 1, np
        rtmp = flux_part(ipart, dir)
        if (rtmp < 0) outflux(ipart) = outflux(ipart) + rtmp
     end do
  end do

  ! Compute correction factor for probability of moving
  proba_correction = 1.0d0

  ! for each direction ...
  do dir = 1, twondim
     ! 'Reverse' the direction
     if (mod(dir, 2) == 1) then ! 1<->2, 3<->4, 5<->6
        ndir = dir + 1
     else
        ndir = dir - 1
     end if
     ! Get the index in grid of cells on the common face of neighbor grid
     call get_cells_on_face(ndir, tmp_ind_ncell(1:twotondim/2))

     ! ... for each particle ...
     do ipart = 1, np

        proba1 = -outflux(ipart)

        ! Store the relative level of the neighbor cell
        itmp = rel_lvl(ind_grid_part(ipart), ison_part(ipart), dir)

        ! ... decide whether it'll move ...
        if (itmp == -1) then
           ! Correct bug with subcycling
           ok = move(ipart) .and. rand1(ipart) < proba1 * nsubcycle(ilevel)
        else
           ok = move(ipart) .and. rand1(ipart) < proba1
        end if

        if (ok) then
           proba2 = flux_part(ipart, dir)/outflux(ipart) * proba_correction(ipart, 1+(dir-1)/2)

           ! ... pick a direction ...
           if (rand2(ipart) < proba2 .and. move(ipart)) then
              ! if (idp(ind_part(ipart)) == 24114) then
              !    print*, 'moving in direction', dir
              ! end if
              ! Tag the particle as moved
              move(ipart) = .false.

              ! === Move to coarser === !
              if (itmp == -1) then
                 factor(ipart) = 2.0D0
                 new_partp(ipart) = ind_ncell(ind_grid_part(ipart), ison_part(ipart), dir)

                 ! === Move to same level === !
              else if (itmp == 0) then
                 factor(ipart) = 1.0D0
                 new_partp(ipart) = ind_ncell(ind_grid_part(ipart), ison_part(ipart), dir)


                 ! === Move to finer === !
              else
                 factor(ipart) = 0.5D0
                 ! Get the fine-to-coarse flux (and invert the sign)
                 do i = 1, twotondim/2
                    ! Compute the neighbor cell index
                    iskip = ncoarse + (tmp_ind_ncell(i)-1)*ngridmax
                    ncell(i) = son(ind_ncell(ind_grid_part(ipart), ison_part(ipart), dir)) &
                         + iskip
                 end do

                 ! Compute the flux from neighbor cell
                 do i = 1, twotondim/2
                    neighborflux(i) = -fluxes(ncell(i), ndir) / twotondim
                 end do

                 ! Recompute the flux in the direction
                 flux_part(ipart, dir) = sum(neighborflux)

                 ! TODO: fix this!
                 if (flux_part(ipart, dir) == 0) then
                    flux_part(ipart, dir) = 1
                    do i = 1, twotondim/2
                       neighborflux(i) = 2d0/twotondim
                    end do
                 end if

                 ! Chose randomly the target cell
                 do i = 1, twotondim/2
                    proba3 = neighborflux(i) / flux_part(ipart, dir)

                    if (rand3(ipart) == 0) cycle
                    if (rand3(ipart) < proba3) then
                       new_partp(ipart) = ncell(i)

                       ! Prevent other directions
                       rand3(ipart) = 0
                    else
                       rand3(ipart) = rand3(ipart) - max(0d0, proba3)
                    end if
                 end do

              end if
           else
              rand2(ipart) = rand2(ipart) - max(0d0, proba2)
           end if
        end if
     end do
  end do

  ! Actually move the particles to their new location
  do ipart = 1, np
     ! Compute new location in grid + father grid position
     ison_part(ipart) = (new_partp(ipart)-ncoarse-1)/ngridmax + 1
     ind_ggrid_part(ipart) = new_partp(ipart) - ncoarse - (ison_part(ipart)-1)*ngridmax

     iz = (ison_part(ipart)-1)/4
     iy = (ison_part(ipart)-1-4*iz)/2
     ix = (ison_part(ipart)-1-4*iz-2*iy)

     x(ipart, 1) = (dble(ix)-0.5d0)
     x(ipart, 2) = (dble(iy)-0.5d0)
     x(ipart, 3) = (dble(iz)-0.5d0)

     do idim = 1, ndim
        new_xp(ipart, idim) = (xg(ind_ggrid_part(ipart), idim) - skip_loc(idim) &
             + x(ipart, idim)*dx*factor(ipart)) * scale
     end do
  end do

  ! Save old positions
  do idim = 1, ndim
     do ipart = 1, np
        old_xp(ipart, idim) = xp(ind_part(ipart), idim)
     end do
  end do

  ! Safely move particles -- taking care of boundaries
  do ipart = 1, np
     tmp_new_xp(:) = new_xp(ipart, :)
     tmp_old_xp(:) = old_xp(ipart, :)
     call safe_move(tmp_new_xp, tmp_old_xp)
     new_xp(ipart, :) = tmp_new_xp(:)
  end do

  ! Velocity field -- mean of previous speed plus current speed
  do idim = 1, ndim
     do ipart = 1, np
        ! Update speed
        vp(ind_part(ipart), idim) = &
             (new_xp(ipart, idim) - old_xp(ipart, idim)) / dtnew(ilevel)
     end do
  end do

  !--------------------------------------------------------------------
  ! Actually move particles
  !--------------------------------------------------------------------
  do idim = 1, ndim
     do ipart = 1, np
        xp(ind_part(ipart), idim) = new_xp(ipart, idim)
     end do
  end do

  ! Store the new parent (here a cell) of the particle
  do ipart = 1, np
     partp(ind_part(ipart)) = new_partp(ipart)
  end do

end subroutine move_gas_tracer
