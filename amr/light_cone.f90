! =======================================================================
! =======================================================================
! =======================================================================
! =======================================================================
subroutine output_cone()
    use amr_commons
    use pm_commons
    use mpi_mod
    use file_module, ONLY: mkdir
    implicit none

#ifndef WITHOUTMPI
    integer :: info, info2, dummy_io
#endif

    integer, parameter :: tag=1118

    character(len=5) :: istep_str
    character(len=100) :: conedir, conecmd, conefile

    integer :: ilun, ipout, npout, npart_out
    character(LEN=80) :: fileloc
    character(LEN=5) :: nchar
    real(kind=8), dimension(1:3, 1:nvector), save :: pos, vel, var
    real(kind=8), dimension(:,:), allocatable :: posout, velout, varout
    real(kind=8), dimension(:), allocatable :: zout
    real(kind=8), dimension(:,:), allocatable :: tmparr
    real(sp), dimension(:,:), allocatable :: xp_out, vp_out, mp_out
    real(sp), dimension(:), allocatable :: zp_out
    real(kind=8) :: z1, z2, om0in, omLin, hubin, Lbox
    real(kind=8) :: observer(3), thetay, thetaz, theta, phi
    integer :: igrid, jgrid, ipart, jpart, idim, icpu, ilevel
    integer :: i, ip, npart1
    integer :: nalloc1, nalloc2
    integer, parameter :: mode = int(O'755')
    integer :: ierr

    integer, dimension(1:nvector), save :: ind_part
    logical :: opened
    opened = .false.

    if (nstep_coarse < 2) return

    z2 = 1 / aexp_old - 1d0
    z1 = 1 / aexp - 1d0

    if (z2 > zmax_cone) return
    if (abs(z2 - z1) < 1d-6) return

    theta = 25.
    phi = 17.
    thetay = thetay_cone
    thetaz = thetaz_cone
    om0in = omega_m
    omLin = omega_l
    hubin = h0 / 100.
    Lbox = boxlen_ini / hubin
    observer = (/ Lbox / 2.0, Lbox / 2.0, Lbox / 2.0 /)

    ilun = 3 * ncpu + myid + 103

    ! Determine the filename, dir, etc
    if (myid == 1) write(*, *) 'Computing and dumping lightcone'

    call title(nstep_coarse, istep_str)
    conedir = "cone_" // trim(istep_str) // "/"
    conecmd = "mkdir -p " // trim(conedir)
    if (.not. withoutmkdir) then
        ! if (myid==1) call system(conecmd)
        if (myid == 1) call mkdir(trim(conedir), mode, ierr)
    end if

#ifndef WITHOUTMPI
    call MPI_BARRIER(MPI_COMM_WORLD, info)
#endif

    conefile = trim(conedir) //'cone_'// trim(istep_str) //'.out'
    call title(myid, nchar)
    fileloc = TRIM(conefile) // TRIM(nchar)

    npart_out = 0
    ipout = 0
    npout = 0

    ! Pre-allocate arrays for particle selection -----
    nalloc1 = nvector
    allocate(posout(1:3, 1:nalloc1))
    allocate(velout(1:3, 1:nalloc1))
    allocate(varout(1:3, 1:nalloc1))
    allocate(zout(1:nalloc1))

    nalloc2 = nvector + nstride
    allocate(xp_out(1:nalloc2, 1:3))
    allocate(vp_out(1:nalloc2, 1:3))
    allocate(mp_out(1:nalloc2, 1:3))
    allocate(zp_out(1:nalloc2))

    allocate(tmparr(1:3, 1:nalloc2))
    ! ------------------------------------------------

    ! Wait for the token
#ifndef WITHOUTMPI
    if (IOGROUPSIZECONE > 0) then
        if (mod(myid - 1, IOGROUPSIZECONE) /= 0) then
            call MPI_RECV(dummy_io, 1, MPI_INTEGER, myid - 1 - 1, tag,&
                & MPI_COMM_WORLD, MPI_STATUS_IGNORE, info2)
        end if
    end if
#endif

    ilevel = levelmin
    ! Loop over cpus
    do icpu = 1, ncpu
        ! Loop over grids
        igrid = headl(icpu, ilevel)
        ip = 0
        do jgrid = 1, numbl(icpu, ilevel)
            npart1 = numbp(igrid)  ! Number of particles in the grid
            if (npart1 > 0) then
                ipart = headp(igrid)

                ! Loop over particles
                do jpart = 1, npart1
                    ip = ip + 1
                    ind_part(ip) = ipart
                    if (ip == nvector) then
                        ! Lower left corner of 3x3x3 grid-cube
                        do idim = 1, ndim
                            do i = 1, ip
                                pos(idim, i) = xp(ind_part(i), idim) * Lbox
                                vel(idim, i) = vp(ind_part(i), idim)
                                var(idim, i) = 0 ! No additional property
                            end do
                        end do
                        ! ===========================================================================
                        ! Count selection particles
                        call perform_my_selection(.true., z1, z2, &
                            &                           om0in, omLin, hubin, Lbox, &
                            &                           observer, thetay, thetaz, theta, phi, &
                            &                           pos, vel, var, ip, &
                            &                           posout, velout, varout, zout, npout,.false.)

                        call extend_arrays_if_needed()

                        ! Perform actual selection
                        call perform_my_selection(.false., z1, z2, &
                            &                           om0in, omLin, hubin, Lbox, &
                            &                           observer, thetay, thetaz, theta, phi, &
                            &                           pos, vel, var, ip, &
                            &                           posout, velout, varout, zout, npout,.false.)
                        ! ===========================================================================
                        if (npout > 0) then
                            do idim = 1, ndim
                                do i = 1, npout
                                    xp_out(ipout + i, idim) = real(posout(idim, i) / Lbox, kind = sp)
                                    vp_out(ipout + i, idim) = real(velout(idim, i), kind = sp)
                                    mp_out(ipout + i, idim) = real(varout(idim, i), kind = sp)
                                end do
                            end do
                            do i = 1, npout
                                zp_out(ipout + i) = real(zout(i), kind = sp)
                            end do
                            ipout = ipout + npout
                            npart_out = npart_out + npout
                        end if

                        ip = 0
                    end if
                    if (ipout >= nstride) then
                        if (.not. opened) then
                            open(ilun, file=TRIM(fileloc), form='unformatted')
                            rewind(ilun)
                            write(ilun) ncpu
                            write(ilun) nstride
                            write(ilun) npart
                            opened = .true.
                        end if
                        do idim = 1, ndim
                            write(ilun) xp_out(1:nstride, idim)
                            write(ilun) vp_out(1:nstride, idim)
                        end do
                        write(ilun) zp_out(1:nstride)
                        do idim = 1, ndim
                            do i = 1, ipout - nstride
                                xp_out(i, idim) = xp_out(i + nstride, idim)
                                vp_out(i, idim) = vp_out(i + nstride, idim)
                                mp_out(i, idim) = mp_out(i + nstride, idim)
                            end do
                        end do
                        do i = 1, ipout - nstride
                            zp_out(i) = zp_out(i + nstride)
                        end do
                        ipout = ipout - nstride
                    end if
                    ipart = nextp(ipart)  ! Go to next particle
                end do
                ! End loop over particles
            end if
            igrid = next(igrid)   ! Go to next grid
        end do
        ! End loop over grids

        if (ip > 0) then
            ! Lower left corner of 3x3x3 grid-cube
            do idim = 1, ndim
                do i = 1, ip
                    pos(idim, i) = xp(ind_part(i), idim) * Lbox
                    vel(idim, i) = vp(ind_part(i), idim)
                    var(idim, i) = 0 ! No additional property
                end do
            end do
            ! ===========================================================================
            ! Count selection particles
            call perform_my_selection(.true., z1, z2, &
                &                           om0in, omLin, hubin, Lbox, &
                &                           observer, thetay, thetaz, theta, phi, &
                &                           pos, vel, var, ip, &
                &                           posout, velout, varout, zout, npout,.false.)

            call extend_arrays_if_needed()

            ! Perform actual selection
            call perform_my_selection(.false., z1, z2, &
                &                           om0in, omLin, hubin, Lbox, &
                &                           observer, thetay, thetaz, theta, phi, &
                &                           pos, vel, var, ip, &
                &                           posout, velout, varout, zout, npout,.false.)
            ! ===========================================================================
            if (npout > 0) then
                do idim = 1, ndim
                    do i = 1, npout
                        xp_out(ipout + i, idim) = real(posout(idim, i) / Lbox, kind = sp)
                        vp_out(ipout + i, idim) = real(velout(idim, i), kind = sp)
                        mp_out(ipout + i, idim) = real(varout(idim, i), kind = sp)
                    end do
                end do
                do i = 1, npout
                    zp_out(ipout + i) = real(zout(i), kind = sp)
                end do
                ipout = ipout + npout
                npart_out = npart_out + npout
            end if
        end if
        if (ipout >= nstride) then
            if (.not. opened) then
                open(ilun, file=TRIM(fileloc), form='unformatted')
                rewind(ilun)
                write(ilun) ncpu
                write(ilun) nstride
                write(ilun) npart
                opened = .true.
            end if
            do idim = 1, ndim
                write(ilun) xp_out(1:nstride, idim)
                write(ilun) vp_out(1:nstride, idim)
            end do
            write(ilun) zp_out(1:nstride)
            do idim = 1, ndim
                do i = 1, ipout - nstride
                    xp_out(i, idim) = xp_out(i + nstride, idim)
                    vp_out(i, idim) = vp_out(i + nstride, idim)
                    mp_out(i, idim) = mp_out(i + nstride, idim)
                end do
            end do
            do i = 1, ipout - nstride
                zp_out(i) = zp_out(i + nstride)
            end do
            ipout = ipout - nstride
        end if
    end do
    ! End loop over cpus

    if (ipout > 0) then
        if (.not. opened) then
            open(ilun, file=TRIM(fileloc), form='unformatted')
            rewind(ilun)
            write(ilun) ncpu
            write(ilun) nstride
            write(ilun) npart
            opened = .true.
        end if
        do idim = 1, ndim
            write(ilun) xp_out(1:ipout, idim)
            write(ilun) vp_out(1:ipout, idim)
        end do
        write(ilun) zp_out(1:ipout)
    end if

    if (opened) close(ilun)

    if (verbose) write(*, *) 'cone output=', myid, npart_out

    if (npart_out > 0) then
        open(ilun, file=TRIM(fileloc) //".txt", form='formatted')
        rewind(ilun)
        write(ilun, *) ncpu
        write(ilun, *) nstride
        write(ilun, *) npart_out
        write(ilun, *) aexp_old
        write(ilun, *) aexp
        close(ilun)
    end if

    ! Send the token
#ifndef WITHOUTMPI
    if (IOGROUPSIZECONE > 0) then
        if (mod(myid, IOGROUPSIZECONE) /= 0 .and. (myid < ncpu)) then
            dummy_io = 1
            call MPI_SEND(dummy_io, 1, MPI_INTEGER, myid - 1 + 1, tag, &
                & MPI_COMM_WORLD, info2)
        end if
    end if
#endif

    if ((opened .and. (npart_out == 0)) .or. ((.not. opened) .and. (npart_out > 0))) then
        write(*, *) 'Error in output_cone'
        write(*, *) 'npart_out=', npart_out,'opened=', opened
        stop
    end if

contains

    ! Extends (deallocates and reallocates) the arrays
    ! posout, velout, varout, zout, xp_out, vp_out, mp_out and zp_out
    ! after npout has been updated, so they can hold enough particles
    !
    ! Reallocation is done in chunks of size alloc_chunk_size, to avoid
    ! reallocating too frequently.

    subroutine extend_arrays_if_needed()

        ! Allocation chunk size
        integer, parameter :: alloc_chunk_size = 100
        integer :: new_nalloc1, new_nalloc2
        integer :: nchunks1, nchunks2

        if (nalloc1 >= npout .and. nalloc2 >= npout + nstride) return


        ! Compute new array sizes
        nchunks1 = npout / alloc_chunk_size
        if (mod(npout, alloc_chunk_size) > 0) nchunks1 = nchunks1 + 1

        nchunks2 = (npout + nstride) / alloc_chunk_size
        if (mod(npout + nstride, alloc_chunk_size) > 0) nchunks2 = nchunks2 + 1

        new_nalloc1 = nchunks1 * alloc_chunk_size
        new_nalloc2 = nchunks2 * alloc_chunk_size

        ! Resize temp array
        deallocate(tmparr)
        allocate(tmparr(1:3, 1:max(new_nalloc1, new_nalloc2)))


        ! Resize xp_out, vp_out, mp_out, zp_out
        do idim = 1, ndim
            tmparr(idim, 1:nalloc2) = xp_out(1:nalloc2, idim)
        end do
        deallocate(xp_out); allocate(xp_out(1:new_nalloc2, 1:3))
        do idim = 1, ndim
            xp_out(1:nalloc2, idim) = real(tmparr(idim, 1:nalloc2), kind = sp)
        end do

        do idim = 1, ndim
            tmparr(idim, 1:nalloc2) = vp_out(1:nalloc2, idim)
        end do
        deallocate(vp_out); allocate(vp_out(1:new_nalloc2, 1:3))
        do idim = 1, ndim
            vp_out(1:nalloc2, idim) = real(tmparr(idim, 1:nalloc2), kind = sp)
        end do

        do idim = 1, ndim
            tmparr(idim, 1:nalloc2) = mp_out(1:nalloc2, idim)
        end do
        deallocate(mp_out); allocate(mp_out(1:new_nalloc2, 1:3))
        do idim = 1, ndim
            mp_out(1:nalloc2, idim) = real(tmparr(idim, 1:nalloc2), kind = sp)
        end do

        tmparr(1, 1:nalloc2) = zp_out(1:nalloc2)
        deallocate(zp_out); allocate(zp_out(1:new_nalloc2))
        zp_out(1:nalloc2) = real(tmparr(1, 1:nalloc2), kind = sp)

        nalloc2 = new_nalloc2


        ! Resize posout, velout, varout, zout
        do idim = 1, ndim
            tmparr(idim, 1:nalloc1) = posout(idim, 1:nalloc1)
            deallocate(posout); allocate(posout(1:3, 1:new_nalloc1))
        end do
        do idim = 1, ndim
            posout(idim, 1:nalloc1) = tmparr(idim, 1:nalloc1)
        end do

        do idim = 1, ndim
            tmparr(idim, 1:nalloc1) = velout(idim, 1:nalloc1)
        end do
        deallocate(velout); allocate(velout(1:3, 1:new_nalloc1))
        do idim = 1, ndim
            velout(idim, 1:nalloc1) = tmparr(idim, 1:nalloc1)
        end do

        do idim = 1, ndim
            tmparr(idim, 1:nalloc1) = varout(idim, 1:nalloc1)
        end do
        deallocate(varout); allocate(varout(1:3, 1:new_nalloc1))
        do idim = 1, ndim
            varout(idim, 1:nalloc1) = tmparr(idim, 1:nalloc1)
        end do

        tmparr(1, 1:nalloc1) = zout(1:nalloc1)
        deallocate(zout); allocate(zout(1:new_nalloc1))
        zout(1:nalloc1) = tmparr(1, 1:nalloc1)

        nalloc1 = new_nalloc1

    end subroutine extend_arrays_if_needed
end subroutine output_cone




subroutine perform_my_selection(justcount, z1, z2, &
        &                          om0in, omLin, hubin, Lbox, &
        &                          observer, thetay, thetaz, theta, phi, &
        &                          pos, vel, var, npart, &
        &                          posout, velout, varout, zout, npartout, verbose)
    ! ===========================================================================
    ! All the quantities below are real*8 except
    ! juscount : logical
    ! npart,npartout : integer*4
    !
    ! juscount : .true. to just count the particles to be selected. Needed
    ! to allocate appropriately the arrays posout,velout and zout
    ! The parameter npartout becomes then an output given the number
    ! of particles selected.
    ! .false. to perform the actual selection: npartout is an input
    ! then  posout, velout, zout are appropriate outputs.
    !
    ! z1,z2    : the lightcone part of interest is in between z1 and z2,
    ! with z1 < z2. If we consider a redshift z(t) where all the
    ! particles are synchrone, and if coarse timestep is
    ! a fixed dt, it is most appropriate to choose z1 and z2 such
    ! that z1=z(t+dt/2) and z2=z(t-dt/2) to have best accuracy.
    !
    ! om0in    : the value of the cosmological parameter omega0 (typically 0.3)
    !
    ! omLin    : the value of the cosmological constant lambda (typically 0.7)
    !
    ! hubin    : the value of H0/100, where H0 is the present Hubble constant
    ! in km/s/Mpc (typically 0.7)
    !
    ! Lbox     : the comoving size of the simulation box in Mpc (NOT in Mpc/h)
    !
    ! observer(3) : the observer position in the box in Mpc, assuming that
    ! coordinates are in [0,Lbox[
    !
    ! thetay   : half the opening angle in degrees of the lightcone along y direction
    ! (it should be obviously smaller than 90 degrees to avoid catastrophic
    ! behavior). The lightcone is assume to be aligned with x axis (after
    ! appropriates rotations given by angles theta and phi)
    !
    ! thetaz   : half the opening angle in degrees of the lightcone along z direction
    ! Given thetay and thetaz, the area of the survey is thus 4.thetay.thetaz
    !
    ! theta, phi : 2 angles in degrees defining a rotation to avoid alignement of
    ! the lightcone with the major axes of the simulation box.
    ! Example : theta=21, phi=17.
    !
    ! pos(3,npart) : comoving positions of the input particles in Mpc, assumed to be
    ! in [0,Lbox[.
    !
    ! vel(3,npart) : velocities of the input particles (in any unit, it does not
    ! matter)
    !
    ! var(3,npart) : additional properties of the input particles (in any unit, it does not
    ! matter)
    !
    ! npart    : number of input particles to be treated
    !
    ! posout(3,npartout) : output comoving positions of selected particles in Mpc.
    !
    ! velout(3,npartout) : output velocities of selected particles
    !
    ! varout(3,npartout) : output properties of selected particles
    !
    ! zout(npartout) : output redshift of selected particles
    !
    ! npartout : number of selected particles. To be computed appropriately,
    ! this routine must be called with juscount=.true., which will give
    ! npartout as an output. Then this routine must be called with
    ! juscount=.false. with the correct value of npartout, after having
    ! allocated correctly arrays posout,velout,zout.
    ! ===========================================================================
    use amr_parameters, ONLY: nvector
    use constants, ONLY: pi
    implicit none
    logical :: justcount, verbose
    integer :: npart, npartout
    integer :: myint
    real(kind=8) :: z1, z2, om0in, omLin, hubin, Lbox
    real(kind=8) :: Omega0, OmegaL, OmegaR, coverH0
    real(kind=8) :: observer(3), thetay, thetaz, theta, phi
    real(kind=8) :: pos(1:3, 1:nvector), vel(1:3, 1:nvector), var(1:3, 1:nvector)
    real(kind=8) :: posout(3, npartout), velout(3, npartout), varout(3, npartout), zout(npartout)
    real(kind=8) :: coord_distance
    real(kind=8) :: thetarad, phirad, thetayrad, thetazrad, tanybound, tanzbound
    real(kind=8) :: rot(3, 3), rotm1(3, 3), dist1, dist2
    real(kind=8) :: xcoordfr, ycoordfr, zcoordfr, xcoord, ycoord, zcoord
    real(kind=8) :: tany, tanz, dist, vxfr, vyfr, vzfr, dxtest1, dxtest2, facnorm
    real(kind=8) :: small=1d-5

    integer :: nrepxm, nrepxp, nrepym, nrepyp, nrepzm, nrepzp
    integer :: i, j, k, np, npartcount
    logical :: keep_part, fullsky

    if (verbose) write(*, *) 'Entering perform_my_selection'

    ! Initialize cosmological parameters
    call init_cosmo_cone(om0in, omLin, hubin, Omega0, OmegaL, OmegaR, coverH0)

    ! Convert angles in radians
    thetarad = theta * pi / 180.0d0
    phirad = phi * pi / 180.0d0

    ! Compute the rotation matrix and its inverse to be in the appropriate frame
    call compute_rotation_matrix(thetarad, phirad, rot, rotm1)

    ! Compute comoving distance of the photon planes from the observer
    ! dist1,dist2=integral of c.dt/a between zero and z1,z2
    dist1 = coord_distance(z1, Omega0, OmegaL, OmegaR, coverH0)
    dist2 = coord_distance(z2, Omega0, OmegaL, OmegaR, coverH0)

    ! Convert angles in radians
    thetayrad = thetay * pi / 180.0d0
    thetazrad = thetaz * pi / 180.0d0

    ! Compute the set of replica to be considered
    if (thetay < 44d0 .AND. thetaz < 44d0) then
        call compute_replica(thetayrad, thetazrad, dist1, dist2, observer, Lbox, rot, &
            &                       nrepxm, nrepxp, nrepym, nrepyp, nrepzm, nrepzp)
        fullsky = .false.
    else
        ! Compute how many replica are needed
        nrepxm = myint((observer(1) - dist2) / Lbox)
        nrepxp = myint((observer(1) + dist2) / Lbox)
        nrepym = myint((observer(2) - dist2) / Lbox)
        nrepyp = myint((observer(2) + dist2) / Lbox)
        nrepzm = myint((observer(3) - dist2) / Lbox)
        nrepzp = myint((observer(3) + dist2) / Lbox)
        fullsky = .true.
    end if

    facnorm = 1.0d0 / (dist2 - dist1)
    tanybound = tan(thetayrad)
    tanzbound = tan(thetazrad)

    npartcount = 0
    ! loop on all the replica of potential interest
    do k = nrepzm, nrepzp, 1
        do j = nrepym, nrepyp, 1
            do i = nrepxm, nrepxp, 1
                do np = 1, npart
                    xcoordfr = pos(1, np) + Lbox * dble(i) - observer(1)
                    ycoordfr = pos(2, np) + Lbox * dble(j) - observer(2)
                    zcoordfr = pos(3, np) + Lbox * dble(k) - observer(3)

                    ! Rotation to get in the framework of the photon plane
                    xcoord = xcoordfr * rotm1(1, 1) + &
                        & ycoordfr * rotm1(2, 1) + &
                        & zcoordfr * rotm1(3, 1)
                    ycoord = xcoordfr * rotm1(1, 2) + &
                        & ycoordfr * rotm1(2, 2) + &
                        & zcoordfr * rotm1(3, 2)
                    zcoord = xcoordfr * rotm1(1, 3) + &
                        & ycoordfr * rotm1(2, 3) + &
                        & zcoordfr * rotm1(3, 3)

                    dist = sqrt(xcoord ** 2 + ycoord ** 2 + zcoord ** 2)

                    keep_part = .false.

                    if (fullsky) then ! Full sky light cone
                        if (dist > dist1 .and. dist <= dist2) then
                            keep_part = .true.
                        end if
                    else ! Narrow light cone
                        if (xcoord > small) then ! To avoid divergences near the origin
                            tany = abs(ycoord / xcoord)
                            tanz = abs(zcoord / xcoord)
                            if (tany <= tanybound .and. tanz <= tanzbound &
                                &  .and. dist > dist1 .and. dist <= dist2) then
                            keep_part = .true.
                        end if
                    end if
                end if

                if (keep_part) then
                    ! This particle is good, we can add it to the list
                    npartcount = npartcount + 1

                    if (.not. justcount) then
                        posout(1, npartcount) = xcoord
                        posout(2, npartcount) = ycoord
                        posout(3, npartcount) = zcoord

                        ! Velocities are rotated
                        vxfr = vel(1, np)
                        vyfr = vel(2, np)
                        vzfr = vel(3, np)
                        velout(1, npartcount) = vxfr * rotm1(1, 1) + &
                            &               vyfr * rotm1(2, 1) + &
                            &               vzfr * rotm1(3, 1)
                        velout(2, npartcount) = vxfr * rotm1(1, 2) + &
                            &               vyfr * rotm1(2, 2) + &
                            &               vzfr * rotm1(3, 2)
                        velout(3, npartcount) = vxfr * rotm1(1, 3) + &
                            &               vyfr * rotm1(2, 3) + &
                            &               vzfr * rotm1(3, 3)

                        ! Additional properties are just passed
                        varout(1, npartcount) = var(1, np)
                        varout(2, npartcount) = var(2, np)
                        varout(3, npartcount) = var(3, np)

                        ! Compute the redshift of the particle using linear
                        ! interpolation
                        dxtest1 = dist - dist1
                        dxtest2 = dist2 - dist
                        zout(npartcount) = (dxtest1 * z2 + dxtest2 * z1) * facnorm
                    end if
                end if

            end do
        end do
    end do
end do
npartout = npartcount
if (verbose) write(*, *) 'End of perform_my_selection'
end subroutine perform_my_selection

! ===========================================================================
subroutine compute_rotation_matrix(thetashiftrad, phishiftrad, rot, rotm1)
    ! ===========================================================================
    ! Rotations matrixes used to perform the calculations.
    ! theta and phi are expressed in radians
    ! ===========================================================================
    implicit none
    real(kind=8) :: thetashiftrad, phishiftrad
    real(kind=8) :: rot(3, 3), rotm1(3, 3)

    integer :: i, j

    rot(1, 1) = cos(thetashiftrad) * cos(phishiftrad)
    rot(1, 2) = cos(thetashiftrad) * sin(phishiftrad)
    rot(1, 3) = - sin(thetashiftrad)
    rot(2, 1) = - sin(phishiftrad)
    rot(2, 2) = cos(phishiftrad)
    rot(2, 3) = 0.0d0
    rot(3, 1) = cos(phishiftrad) * sin(thetashiftrad)
    rot(3, 2) = sin(phishiftrad) * sin(thetashiftrad)
    rot(3, 3) = cos(thetashiftrad)
    do j = 1, 3
        do i = 1, 3
            rotm1(i, j) = rot(j, i)
        end do
    end do
end subroutine compute_rotation_matrix

! ===========================================================================
subroutine compute_minimum_polygon(x1, x2, thetayrad, thetazrad, sl)
    ! ===========================================================================
    ! A slice of photons between redshifts z1 and z2 corresponding to coordinates
    ! x1 and x2 at its center and of opening angles thetay and thetaz is considered.
    ! We compute the coordinates of the eights points of the mimimum (simple)
    ! polygon containing it.
    ! ===========================================================================
    implicit none
    real(kind=8) :: x1, x2, thetayrad, thetazrad, sl(3, 8)

    ! Part of the polygon close to the observer
    sl(1, 1:4) = x1 / sqrt(1.0d0 + tan(thetayrad) ** 2 + tan(thetazrad) ** 2)
    sl(2, 1) = - sl(1, 1) * tan(thetayrad)
    sl(3, 1) = - sl(1, 1) * tan(thetazrad)
    sl(2, 2) = sl(1, 1) * tan(thetayrad)
    sl(3, 2) = - sl(1, 1) * tan(thetazrad)
    sl(2, 3) = - sl(1, 1) * tan(thetayrad)
    sl(3, 3) = sl(1, 1) * tan(thetazrad)
    sl(2, 4) = sl(1, 1) * tan(thetayrad)
    sl(3, 4) = sl(1, 1) * tan(thetazrad)


    ! Part of the polygon far away from the observer
    sl(1, 5:8) = x2
    sl(2, 5) = - x2 * tan(thetayrad)
    sl(3, 5) = - x2 * tan(thetazrad)
    sl(2, 6) = x2 * tan(thetayrad)
    sl(3, 6) = - x2 * tan(thetazrad)
    sl(2, 7) = - x2 * tan(thetayrad)
    sl(3, 7) = x2 * tan(thetazrad)
    sl(2, 8) = x2 * tan(thetayrad)
    sl(3, 8) = x2 * tan(thetazrad)
end subroutine compute_minimum_polygon

! ===========================================================================
subroutine compute_replica(thetayrad, thetazrad, dist1, dist2, observer, Lbox, rot, &
        &                           nrepxm, nrepxp, nrepym, nrepyp, nrepzm, nrepzp)
    ! ===========================================================================
    ! 2*theta1 and 2*theta2 are the opening angles of the lightcone in degrees.
    ! The observer position is expressed in comoving Mpc, as well as the simulation
    ! box size Lbox. Furthermore, the positions of particles inside the simulation
    ! are supposed to be in [0,Lbox[.
    ! z1 and z2 are the redshifts of the successive photon planes, z1 < z2
    ! ===========================================================================
    implicit none
    real(kind=8) :: thetayrad, thetazrad, observer(3), Lbox, rot(3, 3), dist1, dist2
    integer :: nrepxm, nrepxp, nrepym, nrepyp, nrepzm, nrepzp
    integer :: myint
    real(kind=8) :: sl(3, 8), slfr(3)
    real(kind=8) :: xplmin=0, xplmax=0, yplmin=0, yplmax=0, zplmin=0, zplmax=0
    integer :: i, j

    ! Compute the minimum polygon containing the 2 plans of photons (which
    ! are slightly curved)
    call compute_minimum_polygon(dist1, dist2, thetayrad, thetazrad, sl)

    ! Rotate the minimum polygon in the reference frame of the simulation
    do j = 1, 8
        do i = 1, 3
            slfr(i) = sl(1, j) * rot(1, i) &
                & + sl(2, j) * rot(2, i) &
                & + sl(3, j) * rot(3, i)
        end do
        if (j == 1) then
            xplmin = slfr(1)
            xplmax = xplmin
            yplmin = slfr(2)
            yplmax = yplmin
            zplmin = slfr(3)
            zplmax = zplmin
        else
            xplmin = min(xplmin, slfr(1))
            xplmax = max(xplmax, slfr(1))
            yplmin = min(yplmin, slfr(2))
            yplmax = max(yplmax, slfr(2))
            zplmin = min(zplmin, slfr(3))
            zplmax = max(zplmax, slfr(3))
        end if
    end do

    ! Uses the fact that a cube will contain the minimum polygon if and only
    ! if all its edges are contained in the cube to compute the relevant
    ! replica
    nrepxm = myint((xplmin + observer(1)) / Lbox)
    nrepxp = myint((xplmax + observer(1)) / Lbox)
    nrepym = myint((yplmin + observer(2)) / Lbox)
    nrepyp = myint((yplmax + observer(2)) / Lbox)
    nrepzm = myint((zplmin + observer(3)) / Lbox)
    nrepzp = myint((zplmax + observer(3)) / Lbox)
end subroutine compute_replica


! ===================
! cone cosmo routines
! ===================


! ===========================================================================
subroutine init_cosmo_cone(om0in, omLin, hubin, Omega0, OmegaL, OmegaR, coverH0)
    ! ===========================================================================
    ! om0in : the value of omega0
    ! omLin : the value of Lambda
    ! We MUST have omega0+Lambda=1.0d0
    ! hubin : the value of H0/100 where H0 is the present Hubble constant
    ! in km/s/Mpc
    ! ===========================================================================
    implicit none
    real(kind=8) :: om0in, omLin, hubin
    real(kind=8) :: Omega0, OmegaL, OmegaR, coverH0
    real(kind=8) :: verysmall=1d-6
    omega0 = om0in
    omegaL = omLin
    omegaR = 1.0d0 - omega0 - omegaL
    if (abs(omegaR) > verysmall) then
        write(*, *) 'ERROR in propagate_photons, init_cosmo.'
        write(*, *) 'This routine works only for flat universes, omega0+Lambda=1.'
        STOP
    end if
    coverH0 = 299792.5d0 / (100.0d0 * hubin)
end subroutine init_cosmo_cone


! ===========================================================================
function coord_distance(zz, Omega0, OmegaL, OmegaR, coverH0)
    ! ===========================================================================
    implicit none
    real(kind=8) :: z, res, coord_distance, zz
    real(kind=8) :: Omega0, OmegaL, OmegaR, coverH0
    z = abs(zz)
    call qromb(0d0, z, res, omega0, omegaL, OmegaR)
    coord_distance = coverH0 * res
    if (zz < 0) coord_distance = - coord_distance
end function coord_distance

! ===========================================================================
function funcE(z, Omega0, OmegaL, OmegaR)
    ! ===========================================================================
    implicit none
    real(kind=8) :: funcE, z, HsurH0
    real(kind=8) :: omega0, omegaL, OmegaR

    funcE = 1d0 / HsurH0(z, Omega0, OmegaL, OmegaR)
end function funcE

! ===========================================================================
function HsurH0(z, omega0, omegaL, OmegaR)
    ! ===========================================================================
    implicit none
    real(kind=8) :: z, omega0, omegaL, OmegaR, HsurH0
    HsurH0 = sqrt(Omega0 * (1d0 + z) ** 3 + OmegaR * (1d0 + z) ** 2 + OmegaL)
end function HsurH0


! ===========================================================================
SUBROUTINE qromb(a, b, ss, omega0, omegaL, OmegaR)
    ! ===========================================================================
    implicit none
    integer :: JMAX, JMAXP, K, KM
    real(kind=8) :: a, b, ss, EPS, omega0, omegaL, OmegaR
    PARAMETER (EPS = 1d-6, JMAX = 20, JMAXP = JMAX + 1, K = 5, KM = K - 1)
    ! USES polint,trapzd
    integer :: j
    real(kind=8) :: dss, h(JMAXP), s(JMAXP)
    h(1) = 1.
    do j = 1, JMAX
        call trapzd(a, b, s(j), j, omega0, omegaL, OmegaR)
        if (j >= K) then
            call polint(h(j - KM), s(j - KM), K, 0d0, ss, dss)
            if (abs(dss) <= EPS * abs(ss)) return
        end if
        s(j + 1) = s(j)
        h(j + 1) = 0.25 * h(j)
    end do

    print *, 'too many steps in qromb'
END SUBROUTINE qromb

! ===========================================================================
SUBROUTINE polint(xa, ya, n, x, y, dy)
    ! ===========================================================================
    implicit none
    integer :: n, NMAX
    real(kind=8) :: dy, x, y, xa(n), ya(n)
    PARAMETER (NMAX = 10)
    integer :: i, m, ns
    real(kind=8) :: den, dif, dift, ho, hp, w, c(NMAX), d(NMAX)
    ns = 1
    dif = abs(x - xa(1))
    do i = 1, n
        dift = abs(x - xa(i))
        if (dift < dif) then
            ns = i
            dif = dift
        end if
        c(i) = ya(i)
        d(i) = ya(i)
    end do
    y = ya(ns)
    ns = ns - 1
    do m = 1, n - 1
        do i = 1, n - m
            ho = xa(i) - x
            hp = xa(i + m) - x
            w = c(i + 1) - d(i)
            den = ho - hp
            if (den == 0.) print *, 'failure in polint'
            den = w / den
            d(i) = hp * den
            c(i) = ho * den
        end do
        if (2 * ns < n - m) then
            dy = c(ns + 1)
        else
            dy = d(ns)
            ns = ns - 1
        end if
        y = y + dy
    end do
    return
END SUBROUTINE polint

! ===========================================================================
SUBROUTINE trapzd(a, b, s, n, omega0, omegaL, OmegaR)
    ! ===========================================================================
    implicit none
    integer :: n
    real(kind=8) :: a, b, s, funcE, omega0, omegaL, OmegaR
    integer :: it, j
    real(kind=8) :: del, sum, tnm, x
    if (n == 1) then
        s = 0.5 * (b - a) * (funcE(a, omega0, omegaL, OmegaR) + funcE(b, omega0, omegaL, OmegaR))
    else
        it = 2 ** (n - 2)
        tnm = it
        del = (b - a) / tnm
        x = a + 0.5 * del
        sum = 0.
        do j = 1, it
            sum = sum + funcE(x, omega0, omegaL, OmegaR)
            x = x + del
        end do
        s = 0.5 * (s + (b - a) * sum / tnm)
    end if
    return
END SUBROUTINE trapzd
! =======================================================================
function myint(x)
    ! =======================================================================
    ! The REAL int function
    ! =======================================================================
    real(kind=8) :: x
    integer :: myint

    if (x >= 0.0d0) then
        myint = int(x)
    else
        myint = int(x) - 1
    end if
end function myint
