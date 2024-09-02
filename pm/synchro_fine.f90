subroutine synchro_fine(ilevel)
    use pm_commons
    use amr_commons
    use mpi_mod
    implicit none
#ifndef WITHOUTMPI
    integer :: info
#endif
    integer :: ilevel
    ! --------------------------------------------------------------------
    ! This routine synchronizes particle velocity with particle
    ! position for ilevel particle only. If particle sits entirely
    ! in level ilevel, then use inverse CIC at fine level to compute
    ! the force. Otherwise, use coarse level force and coarse level CIC.
    ! --------------------------------------------------------------------
    integer :: igrid, jgrid, ipart, jpart
    integer :: ig, ip, npart1, isink, local_counter
    integer, dimension(1:nvector), save :: ind_grid, ind_part, ind_grid_part

    if (numbtot(1, ilevel) == 0) return
    if (verbose) write(*, 111) ilevel

    if (sink) then
        fsink_new = 0
    end if

    ! Synchronize velocity using CIC
    ig = 0
    ip = 0
    ! Loop over grids
    igrid = headl(myid, ilevel)
    do jgrid = 1, numbl(myid, ilevel)
        npart1 = numbp(igrid)  ! Number of particles in the grid
        if (npart1 > 0) then
            ig = ig + 1
            ind_grid(ig) = igrid
            ipart = headp(igrid)
            local_counter = 0
            ! Loop over particles
            do jpart = 1, npart1
                if (ig == 0) then
                    ig = 1
                    ind_grid(ig) = igrid
                end if
                ! MC TRACER PATCH
                if (is_not_tracer(typep(ipart))) then
                    local_counter = local_counter + 1
                    ip = ip + 1
                    ind_part(ip) = ipart
                    ind_grid_part(ip) = ig
                    if (ip == nvector) then
                        call sync(ind_grid, ind_part, ind_grid_part, ig, ip, ilevel)
                        local_counter = 0
                        ip = 0
                        ig = 0
                    end if
                else
                    ! Update level for tracer particles
                    levelp(ipart) = ilevel
                end if
                ipart = nextp(ipart)  ! Go to next particle
            end do
            ! End loop over particles
            ! If there was no particle in the grid, remove the grid from the buffer
            if (local_counter == 0 .and. ig > 0) then
                ig = ig - 1
            end if
        end if
        igrid = next(igrid)   ! Go to next grid
    end do
    ! End loop over grids
    if (ip > 0) call sync(ind_grid, ind_part, ind_grid_part, ig, ip, ilevel)

    ! sink cloud particles are used to average the grav. acceleration
    if (sink) then
        if (nsink > 0) then
#ifndef WITHOUTMPI
            call MPI_ALLREDUCE(fsink_new, fsink_all, nsinkmax * ndim, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
#else
            fsink_all = fsink_new
#endif
        end if
        do isink = 1, nsink
            if (.not. direct_force_sink(isink)) then
                fsink_partial(isink, 1:ndim, ilevel) = fsink_all(isink, 1:ndim)
            end if
        end do
    end if

111 format('   Entering synchro_fine for level ', I2)

end subroutine synchro_fine
!####################################################################
!####################################################################
!####################################################################
!####################################################################
subroutine synchro_fine_static(ilevel)
    use pm_commons
    use amr_commons
    use mpi_mod
    implicit none
#ifndef WITHOUTMPI
    integer :: info
#endif
    integer :: ilevel
    ! --------------------------------------------------------------------
    ! This routine synchronizes particle velocity with particle
    ! position for ilevel particle only. If particle sits entirely
    ! in level ilevel, then use inverse CIC at fine level to compute
    ! the force. Otherwise, use coarse level force and coarse level CIC.
    ! --------------------------------------------------------------------
    integer :: igrid, jgrid, ipart, jpart
    integer :: ig, ip, next_part, npart1, npart2, isink
    integer, dimension(1:nvector), save :: ind_grid, ind_part, ind_grid_part

    if (numbtot(1, ilevel) == 0) return
    if (verbose) write(*, 111) ilevel

    if (sink) then
        fsink_new = 0
        fsink_all = 0
    end if

    ! Synchronize velocity using CIC
    ig = 0
    ip = 0
    ! Loop over grids
    igrid = headl(myid, ilevel)
    do jgrid = 1, numbl(myid, ilevel)
        npart1 = numbp(igrid)  ! Number of particles in the grid
        npart2 = 0

        ! Count particles
        if (npart1 > 0) then
            ipart = headp(igrid)
            ! Loop over particles
            do jpart = 1, npart1
                ! Save next particle   <--- Very important !!!
                next_part = nextp(ipart)
                if (star) then
                    if ( (.not. static_DM .and. is_DM(typep(ipart))) .or. &
                        & (.not. static_stars .and. (is_star(typep(ipart)) .or. is_debris(typep(ipart))) )  ) then
                    ! FIXME: there should be a static_sink as well
                    npart2 = npart2 + 1
                end if
            else
                if (.not. static_dm) then
                    npart2 = npart2 + 1
                end if
            end if
            ipart = next_part  ! Go to next particle
        end do
    end if

    ! Gather star particles
    if (npart2 > 0) then
        ig = ig + 1
        ind_grid(ig) = igrid
        ipart = headp(igrid)
        ! Loop over particles
        do jpart = 1, npart1
            ! Save next particle   <--- Very important !!!
            next_part = nextp(ipart)
            ! Select particles
            if (star) then
                if ( (.not. static_DM .and. is_DM(typep(ipart))) .or. &
                    & (.not. static_stars .and. (is_star(typep(ipart)) .or. is_debris(typep(ipart))) )  ) then
                ! FIXME: what about sinks?
                if (ig == 0) then
                    ig = 1
                    ind_grid(ig) = igrid
                end if
                ip = ip + 1
                ind_part(ip) = ipart
                ind_grid_part(ip) = ig
            end if
        else
            if (.not. static_dm) then
                if (ig == 0) then
                    ig = 1
                    ind_grid(ig) = igrid
                end if
                ip = ip + 1
                ind_part(ip) = ipart
                ind_grid_part(ip) = ig
            end if
        end if
        if (ip == nvector) then
            call sync(ind_grid, ind_part, ind_grid_part, ig, ip, ilevel)
            ip = 0
            ig = 0
        end if
        ipart = next_part  ! Go to next particle
    end do
    ! End loop over particles
end if
igrid = next(igrid)   ! Go to next grid
end do
! End loop over grids
if (ip > 0) call sync(ind_grid, ind_part, ind_grid_part, ig, ip, ilevel)

! sink cloud particles are used to average the grav. acceleration
if (sink) then
    if (nsink > 0) then
#ifndef WITHOUTMPI
        call MPI_ALLREDUCE(fsink_new, fsink_all, nsinkmax * ndim, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
#else
        fsink_all = fsink_new
#endif
    end if
    do isink = 1, nsink
        if (.not. direct_force_sink(isink)) then
            fsink_partial(isink, 1:ndim, ilevel) = fsink_all(isink, 1:ndim)
        end if
    end do
end if

111 format('   Entering synchro_fine for level ', I2)

end subroutine synchro_fine_static
!####################################################################
!####################################################################
!####################################################################
!####################################################################
subroutine sync(ind_grid, ind_part, ind_grid_part, ng, np, ilevel)
    use amr_commons
    use pm_commons
    use poisson_commons
    implicit none
    integer :: ng, np, ilevel
    integer, dimension(1:nvector) :: ind_grid
    integer, dimension(1:nvector) :: ind_grid_part, ind_part
    !
    !
    !
    logical :: error
    integer :: i, j, ind, idim, nx_loc, isink
    real(dp) :: dx, scale
    ! Grid-based arrays
    real(dp), dimension(1:nvector, 1:ndim), save :: x0
    integer , dimension(1:nvector), save :: ind_cell
    integer , dimension(1:nvector, 1:threetondim), save :: nbors_father_cells
    integer , dimension(1:nvector, 1:twotondim), save :: nbors_father_grids
    ! Particle-based arrays
    logical , dimension(1:nvector), save :: ok
    real(dp), dimension(1:nvector), save :: dteff
    real(dp), dimension(1:nvector, 1:ndim), save :: x, ff, new_vp, dd, dg
    integer , dimension(1:nvector, 1:ndim), save :: ig, id, igg, igd, icg, icd
    real(dp), dimension(1:nvector, 1:twotondim), save :: vol
    integer , dimension(1:nvector, 1:twotondim), save :: igrid, icell, indp, kg
    real(dp), dimension(1:3) :: skip_loc

    ! Mesh spacing in that level
    dx = 0.5D0 ** ilevel
    nx_loc = (icoarse_max - icoarse_min + 1)
    skip_loc = (/ 0.0d0, 0.0d0, 0.0d0 /)
    if (ndim > 0) skip_loc(1) = dble(icoarse_min)
    if (ndim > 1) skip_loc(2) = dble(jcoarse_min)
    if (ndim > 2) skip_loc(3) = dble(kcoarse_min)
    scale = boxlen / dble(nx_loc)

    ! Lower left corner of 3x3x3 grid-cube
    do idim = 1, ndim
        do i = 1, ng
            x0(i, idim) = xg(ind_grid(i), idim) - 3.0D0 * dx
        end do
    end do

    ! Gather 27 neighboring father cells (should be present anytime !)
    do i = 1, ng
        ind_cell(i) = father(ind_grid(i))
    end do
    call get3cubefather(ind_cell, nbors_father_cells, nbors_father_grids, ng, ilevel)

    ! Rescale position at level ilevel
    do idim = 1, ndim
        do j = 1, np
            x(j, idim) = xp(ind_part(j), idim) / scale + skip_loc(idim)
        end do
    end do
    do idim = 1, ndim
        do j = 1, np
            x(j, idim) = x(j, idim) - x0(ind_grid_part(j), idim)
        end do
    end do
    do idim = 1, ndim
        do j = 1, np
            x(j, idim) = x(j, idim) / dx
        end do
    end do

    ! Check for illegal moves
    error = .false.
    do idim = 1, ndim
        do j = 1, np
            if (x(j, idim) < 0.5D0 .or. x(j, idim) > 5.5D0) error = .true.
        end do
    end do
    if (error) then
        write(*, *) 'problem in sync'
        do idim = 1, ndim
            do j = 1, np
                if (x(j, idim) < 0.5D0 .or. x(j, idim) > 5.5D0) then
                    write(*, *) x(j, 1:ndim)
                end if
            end do
        end do
        stop
    end if

    ! CIC at level ilevel (dd: right cloud boundary; dg: left cloud boundary)
    do idim = 1, ndim
        do j = 1, np
            dd(j, idim) = x(j, idim) + 0.5D0
            id(j, idim) = int(dd(j, idim))
            dd(j, idim) = dd(j, idim) - id(j, idim)
            dg(j, idim) = 1.0D0 - dd(j, idim)
            ig(j, idim) = id(j, idim) - 1
        end do
    end do

    ! Compute parent grids
    do idim = 1, ndim
        do j = 1, np
            igg(j, idim) = ig(j, idim) / 2
            igd(j, idim) = id(j, idim) / 2
        end do
    end do
#if NDIM == 1
    do j = 1, np
        kg(j, 1) = 1 + igg(j, 1)
        kg(j, 2) = 1 + igd(j, 1)
    end do
#endif
#if NDIM == 2
    do j = 1, np
        kg(j, 1) = 1 + igg(j, 1) + 3 * igg(j, 2)
        kg(j, 2) = 1 + igd(j, 1) + 3 * igg(j, 2)
        kg(j, 3) = 1 + igg(j, 1) + 3 * igd(j, 2)
        kg(j, 4) = 1 + igd(j, 1) + 3 * igd(j, 2)
    end do
#endif
#if NDIM == 3
    do j = 1, np
        kg(j, 1) = 1 + igg(j, 1) + 3 * igg(j, 2) + 9 * igg(j, 3)
        kg(j, 2) = 1 + igd(j, 1) + 3 * igg(j, 2) + 9 * igg(j, 3)
        kg(j, 3) = 1 + igg(j, 1) + 3 * igd(j, 2) + 9 * igg(j, 3)
        kg(j, 4) = 1 + igd(j, 1) + 3 * igd(j, 2) + 9 * igg(j, 3)
        kg(j, 5) = 1 + igg(j, 1) + 3 * igg(j, 2) + 9 * igd(j, 3)
        kg(j, 6) = 1 + igd(j, 1) + 3 * igg(j, 2) + 9 * igd(j, 3)
        kg(j, 7) = 1 + igg(j, 1) + 3 * igd(j, 2) + 9 * igd(j, 3)
        kg(j, 8) = 1 + igd(j, 1) + 3 * igd(j, 2) + 9 * igd(j, 3)
    end do
#endif
    do ind = 1, twotondim
        do j = 1, np
            igrid(j, ind) = son(nbors_father_cells(ind_grid_part(j), kg(j, ind)))
        end do
    end do

    ! Check if particles are entirely in level ilevel
    ok(1:np) = .true.
    do ind = 1, twotondim
        do j = 1, np
            ok(j) = ok(j) .and. igrid(j, ind) > 0
        end do
    end do

    ! If not, rescale position at level ilevel-1
    do idim = 1, ndim
        do j = 1, np
            if (.not. ok(j)) then
                x(j, idim) = x(j, idim) / 2.0D0
            end if
        end do
    end do
    ! If not, redo CIC at level ilevel-1
    do idim = 1, ndim
        do j = 1, np
            if (.not. ok(j)) then
                dd(j, idim) = x(j, idim) + 0.5D0
                id(j, idim) = int(dd(j, idim))
                dd(j, idim) = dd(j, idim) - id(j, idim)
                dg(j, idim) = 1.0D0 - dd(j, idim)
                ig(j, idim) = id(j, idim) - 1
            end if
        end do
    end do

    ! Compute parent cell position
    do idim = 1, ndim
        do j = 1, np
            if (ok(j)) then
                icg(j, idim) = ig(j, idim) - 2 * igg(j, idim)
                icd(j, idim) = id(j, idim) - 2 * igd(j, idim)
            else
                icg(j, idim) = ig(j, idim)
                icd(j, idim) = id(j, idim)
            end if
        end do
    end do
#if NDIM == 1
    do j = 1, np
        icell(j, 1) = 1 + icg(j, 1)
        icell(j, 2) = 1 + icd(j, 1)
    end do
#endif
#if NDIM == 2
    do j = 1, np
        if (ok(j)) then
            icell(j, 1) = 1 + icg(j, 1) + 2 * icg(j, 2)
            icell(j, 2) = 1 + icd(j, 1) + 2 * icg(j, 2)
            icell(j, 3) = 1 + icg(j, 1) + 2 * icd(j, 2)
            icell(j, 4) = 1 + icd(j, 1) + 2 * icd(j, 2)
        else
            icell(j, 1) = 1 + icg(j, 1) + 3 * icg(j, 2)
            icell(j, 2) = 1 + icd(j, 1) + 3 * icg(j, 2)
            icell(j, 3) = 1 + icg(j, 1) + 3 * icd(j, 2)
            icell(j, 4) = 1 + icd(j, 1) + 3 * icd(j, 2)
        end if
    end do
#endif
#if NDIM == 3
    do j = 1, np
        if (ok(j)) then
            icell(j, 1) = 1 + icg(j, 1) + 2 * icg(j, 2) + 4 * icg(j, 3)
            icell(j, 2) = 1 + icd(j, 1) + 2 * icg(j, 2) + 4 * icg(j, 3)
            icell(j, 3) = 1 + icg(j, 1) + 2 * icd(j, 2) + 4 * icg(j, 3)
            icell(j, 4) = 1 + icd(j, 1) + 2 * icd(j, 2) + 4 * icg(j, 3)
            icell(j, 5) = 1 + icg(j, 1) + 2 * icg(j, 2) + 4 * icd(j, 3)
            icell(j, 6) = 1 + icd(j, 1) + 2 * icg(j, 2) + 4 * icd(j, 3)
            icell(j, 7) = 1 + icg(j, 1) + 2 * icd(j, 2) + 4 * icd(j, 3)
            icell(j, 8) = 1 + icd(j, 1) + 2 * icd(j, 2) + 4 * icd(j, 3)
        else
            icell(j, 1) = 1 + icg(j, 1) + 3 * icg(j, 2) + 9 * icg(j, 3)
            icell(j, 2) = 1 + icd(j, 1) + 3 * icg(j, 2) + 9 * icg(j, 3)
            icell(j, 3) = 1 + icg(j, 1) + 3 * icd(j, 2) + 9 * icg(j, 3)
            icell(j, 4) = 1 + icd(j, 1) + 3 * icd(j, 2) + 9 * icg(j, 3)
            icell(j, 5) = 1 + icg(j, 1) + 3 * icg(j, 2) + 9 * icd(j, 3)
            icell(j, 6) = 1 + icd(j, 1) + 3 * icg(j, 2) + 9 * icd(j, 3)
            icell(j, 7) = 1 + icg(j, 1) + 3 * icd(j, 2) + 9 * icd(j, 3)
            icell(j, 8) = 1 + icd(j, 1) + 3 * icd(j, 2) + 9 * icd(j, 3)
        end if
    end do
#endif

    ! Compute parent cell adresses
    do ind = 1, twotondim
        do j = 1, np
            if (ok(j)) then
                indp(j, ind) = ncoarse + (icell(j, ind) - 1) * ngridmax + igrid(j, ind)
            else
                indp(j, ind) = nbors_father_cells(ind_grid_part(j), icell(j, ind))
            end if
        end do
    end do

    ! Compute cloud volumes
#if NDIM == 1
    do j = 1, np
        vol(j, 1) = dg(j, 1)
        vol(j, 2) = dd(j, 1)
    end do
#endif
#if NDIM == 2
    do j = 1, np
        vol(j, 1) = dg(j, 1) * dg(j, 2)
        vol(j, 2) = dd(j, 1) * dg(j, 2)
        vol(j, 3) = dg(j, 1) * dd(j, 2)
        vol(j, 4) = dd(j, 1) * dd(j, 2)
    end do
#endif
#if NDIM == 3
    do j = 1, np
        vol(j, 1) = dg(j, 1) * dg(j, 2) * dg(j, 3)
        vol(j, 2) = dd(j, 1) * dg(j, 2) * dg(j, 3)
        vol(j, 3) = dg(j, 1) * dd(j, 2) * dg(j, 3)
        vol(j, 4) = dd(j, 1) * dd(j, 2) * dg(j, 3)
        vol(j, 5) = dg(j, 1) * dg(j, 2) * dd(j, 3)
        vol(j, 6) = dd(j, 1) * dg(j, 2) * dd(j, 3)
        vol(j, 7) = dg(j, 1) * dd(j, 2) * dd(j, 3)
        vol(j, 8) = dd(j, 1) * dd(j, 2) * dd(j, 3)
    end do
#endif

    ! Gather 3-force
    ff(1:np, 1:ndim) = 0.0D0
    do ind = 1, twotondim
        do idim = 1, ndim
            do j = 1, np
                ff(j, idim) = ff(j, idim) + f(indp(j, ind), idim) * vol(j, ind)
            end do
        end do
    end do

    ! For sink particle only, store contribution to the sink force
    if (sink) then
        do idim = 1, ndim
            do j = 1, np
                if ( is_cloud(typep(ind_part(j))) ) then
                    isink = - idp(ind_part(j))
                    if (.not. direct_force_sink(isink)) then
                        fsink_new(isink, idim) = fsink_new(isink, idim) + ff(j, idim)
                    end if
                end if
            end do
        end do
    end if

    ! Compute individual time steps
    do j = 1, np
        if (levelp(ind_part(j)) >= ilevel) then
            dteff(j) = dtnew(levelp(ind_part(j)))
        else
            dteff(j) = dtold(levelp(ind_part(j)))
        end if
    end do

    ! Update particles level
    do j = 1, np
        levelp(ind_part(j)) = ilevel
    end do

    ! Update 3-velocity
    do idim = 1, ndim
        if (static) then
            do j = 1, np
                new_vp(j, idim) = ff(j, idim)
            end do
        else
            do j = 1, np
                new_vp(j, idim) = vp(ind_part(j), idim) + ff(j, idim) * 0.5D0 * dteff(j)
            end do
        end if
    end do
    do idim = 1, ndim
        do j = 1, np
            vp(ind_part(j), idim) = new_vp(j, idim)
        end do
    end do

    ! For sink particle only, overwrite cloud particle velocity with sink velocity
    if (sink) then
        do idim = 1, ndim
            do j = 1, np
                if ( is_cloud(typep(ind_part(j))) ) then
                    isink = - idp(ind_part(j))
                    ! Remember that vsink is half time step older than other particles
                    vp(ind_part(j), idim) = vsink(isink, idim)
                end if
            end do
        end do
    end if

end subroutine sync
