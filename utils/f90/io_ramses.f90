module io_ramses
    use random
    use utils

    integer, dimension(:), allocatable :: idout
    real(KIND=8), dimension(:,:), allocatable :: xout, vout
    real(KIND=8), dimension(:), allocatable :: mout, ageout, metout

    real(KIND=8), dimension(:,:), allocatable :: xp
    real(KIND=8), dimension(:,:), allocatable :: varp

contains

    subroutine getcell(xcell, ycell, zcell, varcell, levcell, ncell, nvarin, repository, &
            & levelmax, verbose)
        implicit none
        integer :: ncell, nvarin
        character(LEN=128) :: repository
        real(kind=8), dimension(1:ncell) :: xcell, ycell, zcell
        real(kind=8), dimension(1:ncell, 1:nvarin) :: varcell
        integer, dimension(1:ncell) :: levcell
        integer, optional :: levelmax
        logical, optional :: verbose
        ! --------------------------------------------------------------------------
        ! This routine compute the values of the hydro variables at the position
        ! of a set of input mesh points.
        ! INPUT ARGUMENTS:
        ! xcell,ycell,zcell: one dimensional arrays containing the cell coordinates.
        ! ncell:             an integer specifying the number of cell in the arrays.
        ! nvarin:            an integer specifying the number of required hydro variables.
        ! repository:        a character string containing the full path of the ramses output
        ! directory
        ! OUTPUT ARGUMENTS:
        ! varcell:           a 2-dimensional double array containing the cell variables.
        ! levcell:           a one-dimensional integer array containing the cell level of
        ! refinement.
        ! OPTIONAL ARGUMENTS:
        ! verbose:           logical variable activating verbosity
        ! levelmax:          integer variable specifying the required depth in the AMR tree.
        !
        ! CALLING SEQUENCE:
        ! use io_ramses
        ! call getcell(x,y,z,v,l,1000,6,'output_00200',levelmax=14,verbose=.false.)
        !
        ! R. Teyssier 24/04/09
        ! --------------------------------------------------------------------------
        integer, dimension(:), allocatable :: cpu_cell, icell, jcell, kcell, ind_sort_cell
        real(kind=8), dimension(:), allocatable :: order_cell, dcpu

        integer :: ndim, i, j, k, twotondim, indcell
        integer :: ivar, ncpu, lmax=0, levelmin, ifirst, ilast, istart
        integer :: nx, ny, nz, ilevel, igrid
        integer :: nlevelmax
        integer :: ind, ipos
        integer :: ngridmax, icpu, ncpu_read
        integer :: itop
        real :: boxlen
        real :: t
        integer :: ngridtot
        integer :: nvarh
        integer :: nboundary, ngrid_current
        integer :: ix, iy, iz, impi, bit_length, maxdom, ii, jj, kk
        real(KIND=8) :: dx, xx, yy, zz
        real(KIND=8), dimension(:,:), allocatable :: xg
        real(KIND=8), dimension(:,:,:), allocatable :: var
        integer, dimension(:)  , allocatable :: idp, sontop
        integer, dimension(:,:), allocatable :: son, ngridfile, ngridlevel, ngridbound
        real(KIND=8), dimension(1:8, 1:3) :: xc
        real(KIND=8), dimension(1:3) :: xbound=(/ 0d0, 0d0, 0d0 /)
        character(LEN=5) :: nchar, ncharcpu
        character(LEN=80) :: ordering
        character(LEN=128) :: nomfich
        logical :: ok
        real(kind=8), dimension(:), allocatable :: bound_key, xdp
        logical, dimension(:), allocatable :: cpu_read
        integer, dimension(:), allocatable :: cpu_list
        logical :: verbosity=.false.

        if (present(verbose)) verbosity = verbose

        ! -----------------------------------------------
        ! Allocate cell based arrays
        ! -----------------------------------------------
        allocate(cpu_cell(1:ncell), icell(1:ncell), jcell(1:ncell), kcell(1:ncell))
        allocate(ind_sort_cell(1:ncell))
        allocate(order_cell(1:ncell), dcpu(1:ncell))

        ! -----------------------------------------------
        ! Lecture du fichier hydro au format RAMSES
        ! -----------------------------------------------
        ipos = INDEX(repository,'output_')
        nchar = repository(ipos + 7:ipos + 13)
        nomfich = TRIM(repository) //'/hydro_'// TRIM(nchar) //'.out00001'
        inquire(file = nomfich, exist = ok) ! verify input file
        if ( .not. ok ) then
            print *, TRIM(nomfich) //' not found.'
            stop
        end if
        nomfich = TRIM(repository) //'/amr_'// TRIM(nchar) //'.out00001'
        inquire(file = nomfich, exist = ok) ! verify input file
        if ( .not. ok ) then
            print *, TRIM(nomfich) //' not found.'
            stop
        end if

        nomfich = TRIM(repository) //'/amr_'// TRIM(nchar) //'.out00001'
        open(unit=10, file=nomfich, status='old', form='unformatted')
        read(10) ncpu
        read(10) ndim
        read(10) nx, ny, nz
        read(10) nlevelmax
        read(10) ngridmax
        read(10) nboundary
        read(10) ngrid_current
        read(10) boxlen
        close(10)
        twotondim = 2 ** ndim
        xbound = (/ dble(nx / 2), dble(ny / 2), dble(nz / 2) /)
        allocate(sontop(1:nx * ny * nz))
        allocate(ngridfile(1:ncpu + nboundary, 1:nlevelmax))
        allocate(ngridlevel(1:ncpu, 1:nlevelmax))
        if (nboundary > 0) then
            allocate(ngridbound(1:nboundary, 1:nlevelmax))
        end if

        nomfich = TRIM(repository) //'/info_'// TRIM(nchar) //'.txt'
        inquire(file = nomfich, exist = ok) ! verify input file
        if ( .not. ok ) then
            print *, TRIM(nomfich) //' not found.'
            stop
        end if
        open(unit=10, file=nomfich, form='formatted', status='old')
        read(10,*)
        read(10,*)
        read(10,'("levelmin    =",I11)') levelmin
        read(10,*)
        read(10,*)
        read(10,*)
        read(10,*)

        read(10,*)
        read(10,'("time        =",E23.15)') t
        read(10,*)
        read(10,*)
        read(10,*)
        read(10,*)
        read(10,*)
        read(10,*)
        read(10,*)
        read(10,*)
        read(10,*)
        read(10,*)

        read(10,'("ordering type=",A80)') ordering
        if (verbosity) write(*, '(" ordering type=",A20)') TRIM(ordering)
        read(10,*)
        allocate(cpu_list(1:ncpu))
        if (TRIM(ordering) == 'hilbert') then
            allocate(bound_key(0:ncpu))
            allocate(cpu_read(1:ncpu))
            cpu_read = .false.
            do impi = 1, ncpu
                read(10,'(I8,1X,E23.15,1X,E23.15)') i, bound_key(impi - 1), bound_key(impi)
            end do
        end if
        close(10)

        lmax = nlevelmax
        if (present(levelmax)) lmax = MAX(MIN(levelmax, nlevelmax), 1)

        if (verbosity) write(*, *) 'time=', t
        if (TRIM(ordering) == 'hilbert') then

            maxdom = 2 ** (nlevelmax + 1)
            do bit_length = 1, 32
                maxdom = maxdom / 2
                if (maxdom <= 1) exit
            end do
            maxdom = 2 ** (nlevelmax + 1)

            do i = 1, ncell
                icell(i) = int(xcell(i) * dble(maxdom))
                jcell(i) = int(ycell(i) * dble(maxdom))
                kcell(i) = int(zcell(i) * dble(maxdom))
            end do

            call hilbert3d(icell, jcell, kcell, order_cell, bit_length, ncell)

            do i = 1, ncell
                cpu_cell(i) = ncpu
                do impi = 1, ncpu
                    if (   order_cell(i) >= bound_key(impi - 1) .and.&
                        & order_cell(i) < bound_key(impi  )) then
                    cpu_cell(i) = impi
                    cpu_read(impi) = .true.
                end if
            end do
        end do

        ncpu_read = 0
        do j = 1, ncpu
            if (cpu_read(j)) then
                ncpu_read = ncpu_read + 1
                cpu_list(ncpu_read) = j
            end if
        end do
    else
        ncpu_read = ncpu
        do j = 1, ncpu
            cpu_list(j) = j
        end do
    end  if

    do i = 1, ncell
        dcpu(i) = cpu_cell(i)
    end do
    call quick_sort(dcpu, ind_sort_cell, ncell)

    if (verbosity) then
        write(*, *) 'Processor list'
        do i = 1, ncpu_read
            write(*, *) cpu_list(i)
        end do
    end if

    ! -----------------------------------------------
    ! Read up variable from the AMR grid
    ! ----------------------------------------------

    ! Loop over processor files
    ifirst = 1
    do k = 1, ncpu_read
        icpu = cpu_list(k)
        ilast = ifirst
        do while (ilast < ncell .and. cpu_cell(ind_sort_cell(ilast)) == icpu)
            ilast = ilast + 1
        end do
        if (ilast < ncell) ilast = ilast - 1

        call title(icpu, ncharcpu)

        ! Open AMR file and skip header
        nomfich = TRIM(repository) //'/amr_'// TRIM(nchar) //'.out'// TRIM(ncharcpu)
        open(unit=10, file=nomfich, status='old', form='unformatted')
        if (verbosity) write(*, *) 'Processing file '// TRIM(nomfich)
        do i = 1, 21
            read(10)
        end do
        ! Read grid numbers
        read(10) ngridlevel
        ngridfile(1:ncpu, 1:nlevelmax) = ngridlevel
        read(10)
        if (nboundary > 0) then
            do i = 1, 2
                read(10)
            end do
            read(10) ngridbound
            ngridfile(ncpu + 1:ncpu + nboundary, 1:nlevelmax) = ngridbound
        end if
        read(10)
        ! ROM: comment the single follwing line for old stuff
        read(10)
        if (TRIM(ordering) == 'bisection') then
            do i = 1, 5
                read(10)
            end do
        else
            read(10)
        end if
        read(10) sontop
        read(10)
        read(10)

        ! Open HYDRO file and skip header
        nomfich = TRIM(repository) //'/hydro_'// TRIM(nchar) //'.out'// TRIM(ncharcpu)
        open(unit=11, file=nomfich, status='old', form='unformatted')
        read(11)
        read(11) nvarh
        read(11)
        read(11)
        read(11)
        read(11)

        ! Compute total number of grids
        ngridtot = 0
        do j = 1, ncpu + nboundary
            do ilevel = 1, nlevelmax
                ngridtot = ngridtot + ngridfile(j, ilevel)
            end do
        end do
        if (verbosity) write(*, *) 'Found ', ngridtot,' grid'

        ! Allocate grid-based arrays
        allocate(xg (1:ngridtot, 1:ndim))
        allocate(son(1:ngridtot, 1:twotondim))
        allocate(var(1:ngridtot, 1:twotondim, 1:nvarh))

        ! Loop over levels
        istart = 0
        do ilevel = 1, lmax

            ! Geometry
            dx = 0.5 ** ilevel
            do ind = 1, twotondim
                iz = (ind - 1) / 4
                iy = (ind - 1 - 4 * iz) / 2
                ix = (ind - 1 - 2 * iy - 4 * iz)
                xc(ind, 1) = (dble(ix) - 0.5D0) * dx
                xc(ind, 2) = (dble(iy) - 0.5D0) * dx
                xc(ind, 3) = (dble(iz) - 0.5D0) * dx
            end do

            ! Loop over domains
            do j = 1, nboundary + ncpu

                ! Read AMR data
                if (ngridfile(j, ilevel) > 0) then
                    allocate(xdp(1:ngridfile(j, ilevel)))
                    allocate(idp(1:ngridfile(j, ilevel)))
                    read(10) ! Skip grid index
                    read(10) ! Skip next index
                    read(10) ! Skip prev index
                    ! Read grid center
                    do ind = 1, ndim
                        read(10) xdp
                        do i = 1, ngridfile(j, ilevel)
                            xg(istart + i, ind) = xdp(i)
                        end do
                    end do
                    read(10) ! Skip father index
                    do ind = 1, 2 * ndim
                        read(10) ! Skip nbor index
                    end do
                    ! Read son index
                    do ind = 1, twotondim
                        read(10) idp
                        do i = 1, ngridfile(j, ilevel)
                            son(istart + i, ind) = idp(i)
                        end do
                    end do
                    ! Skip cpu map
                    do ind = 1, twotondim
                        read(10)
                    end do
                    ! Skip refinement map
                    do ind = 1, twotondim
                        read(10)
                    end do
                end if

                ! Read HYDRO data
                read(11)
                read(11)
                if (ngridfile(j, ilevel) > 0) then
                    ! Read hydro variables
                    do ind = 1, twotondim
                        do ivar = 1, nvarh
                            read(11) xdp
                            do i = 1, ngridfile(j, ilevel)
                                var(istart + i, ind, ivar) = xdp(i)
                            end do
                        end do
                    end do
                    istart = istart + ngridfile(j, ilevel)
                    deallocate(xdp)
                    deallocate(idp)
                end if
            end do
        end do
        ! End loop over levels

        close(10)
        close(11)

        if (verbosity) write(*, *) 'End of file'

        if (verbosity) write(*, *) 'Will process cell ', ifirst,' to ', ilast

        do i = ifirst, ilast
            indcell = ind_sort_cell(i)
            xx = xcell(indcell) + xbound(1)
            yy = ycell(indcell) + xbound(2)
            zz = zcell(indcell) + xbound(3)
            itop = 1 + (nx / 2) + (ny / 2) * nx + (nz / 2) * nx * ny
            igrid = sontop(itop)
            do ilevel = 1, nlevelmax
                ii = 1; jj = 1; kk = 1
                if (xx < xg(igrid, 1)) ii = 0
                if (yy < xg(igrid, 2)) jj = 0
                if (zz < xg(igrid, 3)) kk = 0
                ind = 1 + ii + 2 * jj + 4 * kk
                if (son(igrid, ind) == 0 .or. ilevel == lmax) exit
                igrid = son(igrid, ind)
            end do
            do ivar = 1, min(nvarh, nvarin)
                varcell(indcell, ivar) = var(igrid, ind, ivar)
                levcell(indcell) = ilevel
            end do
        end do

        ifirst = ilast + 1
        deallocate(xg)
        deallocate(var)
        deallocate(son)
    end do
    ! End loop over cpu

end subroutine getcell

! =======================================================================
! =======================================================================
! =======================================================================

subroutine gaspart3(ncpu, ncpu_read, cpu_list, repository, ordering, ndummypart, facdens,&
        & lmin, lmax, xmin, xmax, ymin, ymax, zmin, zmax, mdm, partmass, averdens,&
        & denspartcount)
    implicit none

    integer :: ndummypart, ncpu, ndim, lmax, lmin, nlevelmax, nx, ny, nz, ngrid_current, nvarh, ix, iy, iz, nmm
    integer :: ipos, i, j, k, ngridmax, nboundary, twotondim, ngridtot, ilevel, istart, ind, ivar
    integer :: l, npartlocal, partcount, respart, denspartcount, nmin, nmax, pc, respc, indexcelltmp
    integer :: ncpu_read, icpu, ngridactual, ii, jj, kk, icell, ilowdtot, ncpufull, indexcelltot
    real(kind=8) :: xmin, xmax, ymin, ymax, zmin, zmax, boxlen, gamma, dx, dl
    real(kind=8), dimension(1:3) :: xc
    real(kind=8) :: volume, facdens, xx, mdm
    real(KIND=8) :: partmass, averdens, rnpartlocal, massleft, masslefttot
    integer, dimension(1:ncpu) :: cpu_list

    real(KIND=8), dimension(1:3) :: xbound=(/ 0d0, 0d0, 0d0 /)
    character(LEN=128) :: repository, nomfich
    character(LEN=80) :: ordering
    character(LEN=5) :: ncharcpu, nchar
    real(kind=8), dimension(:), allocatable :: xdp, mleft, partm
    real(KIND=8), dimension(:,:), allocatable :: xxdp
    real(KIND=8), dimension(:,:,:), allocatable :: vvdp
    integer, dimension(:), allocatable :: idp, ilowd, nfake, locind, indexcell, levmaxlev, rindex, flagcell
    integer, dimension(:,:), allocatable :: ngridfile, ngridlevel, ngridbound
    integer, dimension(:,:), allocatable :: sdp

    integer , dimension(1:1, 1:IRandNumSize) :: allseed
    integer , dimension(1:IRandNumSize) :: localseed
    integer :: iseed=0, poisson

    ! Initialize random number generator
    call rans(1, iseed, allseed)
    localseed = allseed(1, 1:IRandNumSize)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Compute gas particle mass
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    partmass = 0d0
    averdens = 0d0
    volume = 0d0
    indexcelltot = 0
    ndummypart = 0

    ngridactual = 0
    do k = 1, ncpu
        call title(k, ncharcpu)

        ! Open AMR file and skip header
        ipos = INDEX(repository,'output_')
        nchar = repository(ipos + 7:ipos + 13)
        nomfich = TRIM(repository) //'/amr_'// TRIM(nchar) //'.out'// TRIM(ncharcpu)
        open(unit=10, file=nomfich, status='old', form='unformatted')
        read(10)
        read(10) ndim
        read(10) nx, ny, nz
        read(10) nlevelmax
        read(10) ngridmax
        read(10) nboundary
        read(10) ngrid_current
        read(10) boxlen
        ngridactual = ngridactual + ngrid_current
        close(10)
    end do

    allocate(indexcell(1:ncpu_read))
    allocate(locind(1:ncpu_read))
    allocate(nfake(1:ncpu_read))
    allocate(partm(1:ncpu_read))
    allocate(levmaxlev(1:ncpu_read))

    do k = 1, ncpu_read
        icpu = cpu_list(k)
        call title(icpu, ncharcpu)

        nfake(k) = 0
        partm(k) = 0d0
        indexcell(k) = 0
        levmaxlev(k) = - 100

        ! Open AMR file and skip header
        ipos = INDEX(repository,'output_')
        nchar = repository(ipos + 7:ipos + 13)
        nomfich = TRIM(repository) //'/amr_'// TRIM(nchar) //'.out'// TRIM(ncharcpu)
        open(unit=10, file=nomfich, status='old', form='unformatted')
        read(10)
        read(10) ndim
        read(10) nx, ny, nz
        read(10) nlevelmax
        read(10) ngridmax
        read(10) nboundary
        read(10) ngrid_current
        read(10) boxlen
        do i = 1, 13
            read(10)
        end do
        twotondim = 2 ** ndim
        xbound = (/ dble(nx / 2), dble(ny / 2), dble(nz / 2) /)

        ! Open HYDRO file and skip header
        nomfich = TRIM(repository) //'/hydro_'// TRIM(nchar) //'.out'// TRIM(ncharcpu)
        open(unit=11, file=nomfich, status='old', form='unformatted')
        read(11)
        read(11) nvarh
        read(11)
        read(11)
        read(11)
        read(11) gamma

        allocate(ngridfile(1:ncpu + nboundary, 1:nlevelmax))
        allocate(ngridlevel(1:ncpu, 1:nlevelmax))
        if (nboundary > 0) then
            allocate(ngridbound(1:nboundary, 1:nlevelmax))
        end if

        ! Read grid numbers
        read(10) ngridlevel
        ngridfile(1:ncpu, 1:nlevelmax) = ngridlevel(1:ncpu, 1:nlevelmax)
        read(10)
        if (nboundary > 0) then
            do i = 1, 2
                read(10)
            end do
            read(10) ngridbound
            ngridfile(ncpu + 1:ncpu + nboundary, 1:nlevelmax) = ngridbound(1:nboundary, 1:nlevelmax)
        end if
        read(10)
        ! ROM: comment the single follwing line for old stuff
        read(10)
        if (TRIM(ordering) == 'bisection') then
            do i = 1, 5
                read(10)
            end do
        else
            read(10)
        end if
        read(10)
        read(10)
        read(10)

        ! Compute total number of grids
        ngridtot = 0
        do j = 1, ncpu
            do ilevel = 1, nlevelmax
                ngridtot = ngridtot + ngridfile(j, ilevel)
            end do
        end do

        ! Loop over levels
        istart = 0
        do ilevel = 1, lmax
            ! Loop over domains
            do j = 1, nboundary + ncpu
                ! Read AMR data
                if (ngridfile(j, ilevel) > 0) then
                    if (j == icpu .and. ilevel > levmaxlev(k)) levmaxlev(k) = ilevel
                    allocate(xdp(1:ngridfile(j, ilevel)))
                    allocate(xxdp(1:ngridfile(j, ilevel), 1:ndim))
                    allocate(vvdp(1:ngridfile(j, ilevel), 1:twotondim, 1:nvarh))
                    allocate(sdp(1:ngridfile(j, ilevel), 1:twotondim))
                    allocate(idp(1:ngridfile(j, ilevel)))
                    read(10) ! Skip grid index
                    read(10) ! Skip next index
                    read(10) ! Skip prev index
                    ! Read grid center
                    do ind = 1, ndim
                        read(10) xdp
                        do i = 1, ngridfile(j, ilevel)
                            xxdp(i, ind) = xdp(i)
                        end do
                    end do
                    read(10) ! Skip father index
                    do ind = 1, 2 * ndim
                        read(10) ! Skip nbor index
                    end do
                    ! Read son index
                    do ind = 1, twotondim
                        read(10) idp
                        do i = 1, ngridfile(j, ilevel)
                            sdp(i, ind) = idp(i)
                        end do
                    end do
                    ! Skip cpu map
                    do ind = 1, twotondim
                        read(10)
                    end do
                    ! Skip refinement map
                    do ind = 1, twotondim
                        read(10)
                    end do
                end if
                ! Read HYDRO data
                read(11)
                read(11)
                if (ngridfile(j, ilevel) > 0) then
                    ! Read hydro variables
                    do ind = 1, twotondim
                        do ivar = 1, nvarh
                            read(11) xdp
                            do i = 1, ngridfile(j, ilevel)

                                dl = 0.5d0 ** dble(lmax)

                                if (ndim == 3) then
                                    vvdp(i, ind, ivar) = xdp(i)
                                    dx = 0.5d0 ** ilevel
                                    iz = (ind - 1) / 4
                                    iy = (ind - 1 - 4 * iz) / 2
                                    ix = (ind - 1 - 2 * iy - 4 * iz)
                                    xc(1) = boxlen * (xxdp(i, 1) + (dble(ix) - 0.5D0) * dx - xbound(1))
                                    xc(2) = boxlen * (xxdp(i, 2) + (dble(iy) - 0.5D0) * dx - xbound(2))
                                    xc(3) = boxlen * (xxdp(i, 3) + (dble(iz) - 0.5D0) * dx - xbound(3))
                                    if (ivar == nvarh .and. j == icpu .and. ilevel >= lmin .and. (sdp(i, ind) == 0&
                                        & .or. ilevel == lmax) .and. (xmin <= xc(1) .and. xc(1) <= xmax) .and.&
                                        & (ymin <= xc(2) .and. xc(2) <= ymax) .and. (zmin <= xc(3) .and. xc(3) <= zmax)) then
                                    volume = volume + dx * dx * dx
                                    partmass = partmass + vvdp(i, ind, 1) * dx * dx * dx
                                    partm(k) = partm(k) + vvdp(i, ind, 1) * dx * dx * dx
                                    indexcell(k) = indexcell(k) + 1
                                    indexcelltot = indexcelltot + 1
                                end if
                            end if

                            if (ndim == 2) then
                                vvdp(i, ind, ivar) = xdp(i)
                                dx = 0.5d0 ** ilevel
                                iz = 0
                                iy = (ind - 1) / 2
                                ix = (ind - 1 - 2 * iy)
                                xc(1) = boxlen * (xxdp(i, 1) + (dble(ix) - 0.5D0) * dx - xbound(1))
                                xc(2) = boxlen * (xxdp(i, 2) + (dble(iy) - 0.5D0) * dx - xbound(2))
                                xc(3) = (zmin + zmax) / 2
                                if (ivar == nvarh .and. j == icpu .and. ilevel >= lmin .and. (sdp(i, ind) == 0&
                                    & .or. ilevel == lmax) .and. (xmin <= xc(1) .and. xc(1) <= xmax) .and.&
                                    & (ymin <= xc(2) .and. xc(2) <= ymax) .and. (zmin <= xc(3) .and. xc(3) <= zmax)) then
                                volume = volume + dx * dx * dl
                                partmass = partmass + vvdp(i, ind, 1) * dx * dx * dl
                                partm(k) = partm(k) + vvdp(i, ind, 1) * dx * dx * dl
                                indexcell(k) = indexcell(k) + 1
                                indexcelltot = indexcelltot + 1
                            end if
                        end if

                        if (ndim == 1) then
                            vvdp(i, ind, ivar) = xdp(i)
                            dx = 0.5d0 ** ilevel
                            iz = 0
                            iy = 0
                            ix = ind - 1
                            xc(1) = boxlen * (xxdp(i, 1) + (dble(ix) - 0.5D0) * dx - xbound(1))
                            xc(2) = (ymin + ymax) / 2
                            xc(3) = (zmin + zmax) / 2
                            if (ivar == nvarh .and. j == icpu .and. ilevel >= lmin .and. (sdp(i, ind) == 0&
                                & .or. ilevel == lmax) .and. (xmin <= xc(1) .and. xc(1) <= xmax) .and.&
                                & (ymin <= xc(2) .and. xc(2) <= ymax) .and. (zmin <= xc(3) .and. xc(3) <= zmax)) then
                            volume = volume + dx * dl * dl
                            partmass = partmass + vvdp(i, ind, 1) * dx * dl * dl
                            partm(k) = partm(k) + vvdp(i, ind, 1) * dx * dl * dl
                            indexcell(k) = indexcell(k) + 1
                            indexcelltot = indexcelltot + 1
                        end if
                    end if

                end do
            end do
        end do
        istart = istart + ngridfile(j, ilevel)
        deallocate(xdp)
        deallocate(xxdp)
        deallocate(vvdp)
        deallocate(sdp)
        deallocate(idp)
    end if
end do
end do
! End loop over levels

close(10)
close(11)

nfake(k) = nint(partm(k) / mdm)

deallocate(ngridfile, ngridlevel)
if (nboundary > 0) deallocate(ngridbound)
end do

write(*, *) 'AMR grid read'

write(*, *) 'Total gas mass = ', partmass

if (partmass / mdm > 1d9) then
    write(*, *) 'Too many gas particle'
    write(*, *) 'Stop'
    stop
end if

ndummypart = int(partmass / mdm) + 1

averdens = partmass / volume
partmass = partmass / dble(ndummypart)

nmin = 1
nmax = ndummypart

write(*, *) 'Number of gas dummy particles = ', ndummypart
write(*, *) 'Gas particle mass = ', partmass
write(*, *) 'Target gas particle mass = ', mdm


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Distribute gas particles
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

nmm = nmax - nmin + 1

allocate(xp(1:nmm, 1:ndim))
allocate(varp(1:nmm, 1:nvarh))
allocate(mleft(1:ncpu_read))
allocate(ilowd(1:ncpu_read))

pc = 0
partcount = 0
denspartcount = 0
masslefttot = 0d0
ilowdtot = 0

do k = 1, ncpu_read
    icpu = cpu_list(k)
    call title(icpu, ncharcpu)

    ! Open AMR file and skip header
    ipos = INDEX(repository,'output_')
    nchar = repository(ipos + 7:ipos + 13)
    nomfich = TRIM(repository) //'/amr_'// TRIM(nchar) //'.out'// TRIM(ncharcpu)
    open(unit=10, file=nomfich, status='old', form='unformatted')
    read(10)
    read(10) ndim
    read(10) nx, ny, nz
    read(10) nlevelmax
    read(10) ngridmax
    read(10) nboundary
    read(10) ngrid_current
    read(10) boxlen
    do i = 1, 13
        read(10)
    end do
    twotondim = 2 ** ndim
    xbound = (/ dble(nx / 2), dble(ny / 2), dble(nz / 2) /)
    allocate(ngridfile(1:ncpu + nboundary, 1:nlevelmax))
    allocate(ngridlevel(1:ncpu, 1:nlevelmax))
    if (nboundary > 0) allocate(ngridbound(1:nboundary, 1:nlevelmax))

    ! Open HYDRO file and skip headernpartlocal=int(rnpartlocal)
    nomfich = TRIM(repository) //'/hydro_'// TRIM(nchar) //'.out'// TRIM(ncharcpu)
    open(unit=11, file=nomfich, status='old', form='unformatted')
    read(11)
    read(11) nvarh
    read(11)
    read(11)
    read(11)
    read(11) gamma

    ! Read grid numbers
    read(10) ngridlevel
    ngridfile(1:ncpu, 1:nlevelmax) = ngridlevel(1:ncpu, 1:nlevelmax)
    read(10)
    if (nboundary > 0) then
        do i = 1, 2
            read(10)
        end do
        read(10) ngridbound
        ngridfile(ncpu + 1:ncpu + nboundary, 1:nlevelmax) = ngridbound(1:nboundary, 1:nlevelmax)
    end if
    read(10)
    ! ROM: comment the single follwing line for old stuff
    read(10)
    if (TRIM(ordering) == 'bisection') then
        do i = 1, 5
            read(10)
        end do
    else
        read(10)
    end if
    read(10)
    read(10)
    read(10)

    massleft = 0d0
    ilowd(k) = 0
    mleft(k) = 0d0
    locind(k) = 0

    ! Loop over levels
    istart = 0
    do ilevel = 1, lmax
        ! Loop over domains
        do j = 1, nboundary + ncpu
            ! Read AMR data
            if (ngridfile(j, ilevel) > 0) then
                allocate(xdp(1:ngridfile(j, ilevel)))
                allocate(xxdp(1:ngridfile(j, ilevel), 1:ndim))
                allocate(vvdp(1:ngridfile(j, ilevel), 1:twotondim, 1:nvarh))
                allocate(sdp(1:ngridfile(j, ilevel), 1:twotondim))
                allocate(idp(1:ngridfile(j, ilevel)))
                read(10) ! Skip grid index
                read(10) ! Skip next index
                read(10) ! Skip prev index
                ! Read grid center
                do ind = 1, ndim
                    read(10) xdp
                    do i = 1, ngridfile(j, ilevel)
                        xxdp(i, ind) = xdp(i)
                    end do
                end do
                read(10) ! Skip father index
                do ind = 1, 2 * ndim
                    read(10) ! Skip nbor index
                end do
                ! Read son index
                do ind = 1, twotondim
                    read(10) idp
                    do i = 1, ngridfile(j, ilevel)
                        sdp(i, ind) = idp(i)
                    end do
                end do
                ! Skip cpu map
                do ind = 1, twotondim
                    read(10)
                end do
                ! Skip refinement map
                do ind = 1, twotondim
                    read(10)
                end do
            end if
            ! Read HYDRO data
            read(11)
            read(11)
            if (ngridfile(j, ilevel) > 0) then
                ! Read hydro variables
                do ind = 1, twotondim
                    do ivar = 1, nvarh
                        read(11) xdp
                        do i = 1, ngridfile(j, ilevel)

                            dl = 0.5d0 ** dble(lmax)

                            if (ndim == 3) then
                                vvdp(i, ind, ivar) = xdp(i)
                                dx = 0.5d0 ** ilevel
                                iz = (ind - 1) / 4
                                iy = (ind - 1 - 4 * iz) / 2
                                ix = (ind - 1 - 2 * iy - 4 * iz)
                                xc(1) = boxlen * (xxdp(i, 1) + (dble(ix) - 0.5D0) * dx - xbound(1))
                                xc(2) = boxlen * (xxdp(i, 2) + (dble(iy) - 0.5D0) * dx - xbound(2))
                                xc(3) = boxlen * (xxdp(i, 3) + (dble(iz) - 0.5D0) * dx - xbound(3))
                                if (j == icpu .and. ivar == nvarh .and. (ilevel >= lmin) .and. (sdp(i, ind) == 0&
                                    & .or. ilevel == lmax) .and. (xmin <= xc(1) .and. xc(1) <= xmax) .and.&
                                    & (ymin <= xc(2) .and. xc(2) <= ymax) .and. (zmin <= xc(3) .and. xc(3) <= zmax)) then
                                rnpartlocal = (vvdp(i, ind, 1) * dx * dx * dx) / partmass
                                if (rnpartlocal >= 1) then
                                    npartlocal = int(rnpartlocal)
                                    massleft = massleft + (vvdp(i, ind, 1) * dx * dx * dx - dble(npartlocal) * partmass)
                                    jj = 1
                                    do while (jj <= npartlocal)
                                        jj = jj + 1
                                        locind(k) = locind(k) + 1
                                        partcount = partcount + 1
                                        if (nmin <= partcount .and. partcount <= nmax) then
                                            pc = pc + 1
                                            do l = 1, ndim
                                                call ranf(localseed, xx)
                                                xp(pc, l) = xx * boxlen * dx + xc(l) - boxlen * dx / 2
                                            end do
                                            do l = 1, nvarh
                                                varp(pc, l) = vvdp(i, ind, l)
                                            end do
                                            if (varp(pc, 1) > facdens * averdens) denspartcount = denspartcount + 1
                                        end if
                                    end do
                                end if
                                if (rnpartlocal < 1) then
                                    ilowd(k) = ilowd(k) + 1
                                    massleft = massleft + vvdp(i, ind, 1) * dx * dx * dx
                                    ilowdtot = ilowdtot + 1
                                    masslefttot = masslefttot + vvdp(i, ind, 1) * dx * dx * dx
                                end if
                            end if
                        end if

                        if (ndim == 2) then
                            vvdp(i, ind, ivar) = xdp(i)
                            dx = 0.5d0 ** ilevel
                            iz = 0
                            iy = (ind - 1) / 2
                            ix = (ind - 1 - 2 * iy)
                            xc(1) = boxlen * (xxdp(i, 1) + (dble(ix) - 0.5D0) * dx - xbound(1))
                            xc(2) = boxlen * (xxdp(i, 2) + (dble(iy) - 0.5D0) * dx - xbound(2))
                            xc(3) = (zmin + zmax) / 2
                            if (j == icpu .and. ivar == nvarh .and. (ilevel >= lmin) .and. (sdp(i, ind) == 0&
                                & .or. ilevel == lmax) .and. (xmin <= xc(1) .and. xc(1) <= xmax) .and.&
                                & (ymin <= xc(2) .and. xc(2) <= ymax) .and. (zmin <= xc(3) .and. xc(3) <= zmax)) then
                            rnpartlocal = (vvdp(i, ind, 1) * dx * dx * dl) / partmass
                            if (rnpartlocal >= 1) then
                                npartlocal = int(rnpartlocal)
                                massleft = massleft + (vvdp(i, ind, 1) * dx * dx * dl - dble(npartlocal) * partmass)
                                jj = 1
                                do while (jj <= npartlocal)
                                    jj = jj + 1
                                    locind(k) = locind(k) + 1
                                    partcount = partcount + 1
                                    if (nmin <= partcount .and. partcount <= nmax) then
                                        pc = pc + 1
                                        do l = 1, ndim
                                            call ranf(localseed, xx)
                                            xp(pc, l) = xx * boxlen * dx + xc(l) - boxlen * dx / 2
                                        end do
                                        do l = 1, nvarh
                                            varp(pc, l) = vvdp(i, ind, l)
                                        end do
                                        if (varp(pc, 1) > facdens * averdens) denspartcount = denspartcount + 1
                                    end if
                                end do
                            end if
                            if (rnpartlocal < 1) then
                                ilowd(k) = ilowd(k) + 1
                                massleft = massleft + vvdp(i, ind, 1) * dx * dx * dl
                                ilowdtot = ilowdtot + 1
                                masslefttot = masslefttot + vvdp(i, ind, 1) * dx * dx * dl
                            end if
                        end if
                    end if

                    if (ndim == 1) then
                        vvdp(i, ind, ivar) = xdp(i)
                        dx = 0.5d0 ** ilevel
                        iz = 0
                        iy = 0
                        ix = ind - 1
                        xc(1) = boxlen * (xxdp(i, 1) + (dble(ix) - 0.5D0) * dx - xbound(1))
                        xc(2) = (ymin + ymax) / 2
                        xc(3) = (zmin + zmax) / 2
                        if (j == icpu .and. ivar == nvarh .and. (ilevel >= lmin) .and. (sdp(i, ind) == 0&
                            & .or. ilevel == lmax) .and. (xmin <= xc(1) .and. xc(1) <= xmax) .and.&
                            & (ymin <= xc(2) .and. xc(2) <= ymax) .and. (zmin <= xc(3) .and. xc(3) <= zmax)) then
                        rnpartlocal = (vvdp(i, ind, 1) * dx * dl * dl) / partmass
                        if (rnpartlocal >= 1) then
                            npartlocal = int(rnpartlocal)
                            massleft = massleft + (vvdp(i, ind, 1) * dx * dl * dl - dble(npartlocal) * partmass)
                            jj = 1
                            do while (jj <= npartlocal)
                                jj = jj + 1
                                locind(k) = locind(k) + 1
                                partcount = partcount + 1
                                if (nmin <= partcount .and. partcount <= nmax) then
                                    pc = pc + 1
                                    do l = 1, ndim
                                        call ranf(localseed, xx)
                                        xp(pc, l) = xx * boxlen * dx + xc(l) - boxlen * dx / 2
                                    end do
                                    do l = 1, nvarh
                                        varp(pc, l) = vvdp(i, ind, l)
                                    end do
                                    if (varp(pc, 1) > facdens * averdens) denspartcount = denspartcount + 1
                                end if
                            end do
                        end if
                        if (rnpartlocal < 1) then
                            ilowd(k) = ilowd(k) + 1
                            massleft = massleft + vvdp(i, ind, 1) * dx * dl * dl
                            ilowdtot = ilowdtot + 1
                            masslefttot = masslefttot + vvdp(i, ind, 1) * dx * dl * dl
                        end if
                    end if
                end if

            end do
        end do
    end do
    istart = istart + ngridfile(j, ilevel)
    deallocate(xdp)
    deallocate(xxdp)
    deallocate(vvdp)
    deallocate(sdp)
    deallocate(idp)
end if
end do
end do

mleft(k) = massleft

close(10)
close(11)

deallocate(ngridfile, ngridlevel)
if (nboundary > 0) deallocate(ngridbound)
end do


write(*, *) 'gaspart3 I done'

if (partcount > ndummypart .or. pc > nmax - nmin + 1) WRITE(*,*)'ERROR! TOO MANY PARTICLES!'

ncpufull = 0

do k = 1, ncpu_read
    icpu = cpu_list(k)
    call title(icpu, ncharcpu)

    if (partcount < ndummypart .and. pc < nmax - nmin + 1) then

        ! Open AMR file and skip header
        ipos = INDEX(repository,'output_')
        nchar = repository(ipos + 7:ipos + 13)
        nomfich = TRIM(repository) //'/amr_'// TRIM(nchar) //'.out'// TRIM(ncharcpu)
        open(unit=10, file=nomfich, status='old', form='unformatted')
        read(10)
        read(10) ndim
        read(10) nx, ny, nz
        read(10) nlevelmax
        read(10) ngridmax
        read(10) nboundary
        read(10) ngrid_current
        read(10) boxlen
        do i = 1, 13
            read(10)
        end do
        twotondim = 2 ** ndim
        xbound = (/ dble(nx / 2), dble(ny / 2), dble(nz / 2) /)
        allocate(ngridfile(1:ncpu + nboundary, 1:nlevelmax))
        allocate(ngridlevel(1:ncpu, 1:nlevelmax))
        if (nboundary > 0) allocate(ngridbound(1:nboundary, 1:nlevelmax))

        ! Open HYDRO file and skip header
        nomfich = TRIM(repository) //'/hydro_'// TRIM(nchar) //'.out'// TRIM(ncharcpu)
        open(unit=11, file=nomfich, status='old', form='unformatted')
        read(11)
        read(11) nvarh
        read(11)
        read(11)
        read(11)
        read(11) gamma

        ! Read grid numbers
        read(10) ngridlevel
        ngridfile(1:ncpu, 1:nlevelmax) = ngridlevel(1:ncpu, 1:nlevelmax)
        read(10)
        if (nboundary > 0) then
            do i = 1, 2
                read(10)
            end do
            read(10) ngridbound
            ngridfile(ncpu + 1:ncpu + nboundary, 1:nlevelmax) = ngridbound(1:nboundary, 1:nlevelmax)
        end if
        read(10)
        ! ROM: comment the single follwing line for old stuff
        read(10)
        if (TRIM(ordering) == 'bisection') then
            do i = 1, 5
                read(10)
            end do
        else
            read(10)
        end if
        read(10)
        read(10)
        read(10)

        ! Loop over levels
        istart = 0
        do ilevel = 1, lmax
            ! Loop over domains
            do j = 1, nboundary + ncpu
                if (ngridfile(j, ilevel) > 0) then
                    allocate(xdp(1:ngridfile(j, ilevel)))
                    allocate(xxdp(1:ngridfile(j, ilevel), 1:ndim))
                    allocate(vvdp(1:ngridfile(j, ilevel), 1:twotondim, 1:nvarh))
                    allocate(sdp(1:ngridfile(j, ilevel), 1:twotondim))
                    allocate(idp(1:ngridfile(j, ilevel)))
                    read(10) ! Skip grid index
                    read(10) ! Skip next index
                    read(10) ! Skip prev index
                    ! Read grid center
                    do ind = 1, ndim
                        read(10) xdp
                        do i = 1, ngridfile(j, ilevel)
                            xxdp(i, ind) = xdp(i)
                        end do
                    end do
                    read(10) ! Skip father index
                    do ind = 1, 2 * ndim
                        read(10) ! Skip nbor index
                    end do
                    ! Read son index
                    do ind = 1, twotondim
                        read(10) idp
                        do i = 1, ngridfile(j, ilevel)
                            sdp(i, ind) = idp(i)
                        end do
                    end do
                    ! Skip cpu map
                    do ind = 1, twotondim
                        read(10)
                    end do
                    ! Skip refinement map
                    do ind = 1, twotondim
                        read(10)
                    end do
                end if
                ! Read HYDRO data
                read(11)
                read(11)
                if (ngridfile(j, ilevel) > 0) then
                    respart = ndummypart - partcount
                    respc = nmax - nmin + 1 - pc
                    do ind = 1, twotondim
                        do ivar = 1, nvarh
                            read(11) xdp
                            do i = 1, ngridfile(j, ilevel)

                                dl = 0.5d0 ** dble(lmax)

                                if (ndim == 3) then
                                    vvdp(i, ind, ivar) = xdp(i)
                                    dx = 0.5d0 ** ilevel
                                    iz = (ind - 1) / 4
                                    iy = (ind - 1 - 4 * iz) / 2
                                    ix = (ind - 1 - 2 * iy - 4 * iz)
                                    xc(1) = boxlen * (xxdp(i, 1) + (dble(ix) - 0.5D0) * dx - xbound(1))
                                    xc(2) = boxlen * (xxdp(i, 2) + (dble(iy) - 0.5D0) * dx - xbound(2))
                                    xc(3) = boxlen * (xxdp(i, 3) + (dble(iz) - 0.5D0) * dx - xbound(3))
                                    if (j == icpu .and. ivar == nvarh .and. (ilevel >= lmin) .and. (sdp(i, ind) == 0&
                                        & .or. ilevel == lmax) .and. (xmin <= xc(1) .and. xc(1) <= xmax) .and.&
                                        & (ymin <= xc(2) .and. xc(2) <= ymax) .and. (zmin <= xc(3) .and. xc(3) <= zmax)) then
                                    rnpartlocal = (vvdp(i, ind, 1) * dx * dx * dx) / partmass
                                    if (rnpartlocal < 1) then
                                        call poissdev(localseed, rnpartlocal, poisson)
                                        poisson = poisson + nint(masslefttot / dble(ilowdtot) / partmass)
                                        jj = 1
                                        do while (jj <= poisson .and. partcount < ndummypart)
                                            jj = jj + 1
                                            locind(k) = locind(k) + 1
                                            partcount = partcount + 1
                                            if (nmin <= partcount .and. partcount <= nmax) then
                                                pc = pc + 1
                                                do l = 1, ndim
                                                    call ranf(localseed, xx)
                                                    xp(pc, l) = xx * boxlen * dx + xc(l) - boxlen * dx / 2
                                                end do
                                                do kk = 1, nvarh
                                                    varp(pc, kk) = vvdp(i, ind, kk)
                                                end do
                                                if (varp(pc, 1) > facdens * averdens) denspartcount = denspartcount + 1
                                            end if
                                        end do
                                    end if
                                end if
                            end if

                            if (ndim == 2) then
                                vvdp(i, ind, ivar) = xdp(i)
                                dx = 0.5d0 ** ilevel
                                iz = 0
                                iy = (ind - 1) / 2
                                ix = (ind - 1 - 2 * iy)
                                xc(1) = boxlen * (xxdp(i, 1) + (dble(ix) - 0.5D0) * dx - xbound(1))
                                xc(2) = boxlen * (xxdp(i, 2) + (dble(iy) - 0.5D0) * dx - xbound(2))
                                xc(3) = (zmin + zmax) / 2
                                if (j == icpu .and. ivar == nvarh .and. (ilevel >= lmin) .and. (sdp(i, ind) == 0&
                                    & .or. ilevel == lmax) .and. (xmin <= xc(1) .and. xc(1) <= xmax) .and.&
                                    & (ymin <= xc(2) .and. xc(2) <= ymax) .and. (zmin <= xc(3) .and. xc(3) <= zmax)) then
                                rnpartlocal = (vvdp(i, ind, 1) * dx * dx * dl) / partmass
                                if (rnpartlocal < 1) then
                                    call poissdev(localseed, rnpartlocal, poisson)
                                    poisson = poisson + nint(masslefttot / dble(ilowdtot) / partmass)
                                    jj = 1
                                    do while (jj <= poisson .and. partcount < ndummypart)
                                        jj = jj + 1
                                        locind(k) = locind(k) + 1
                                        partcount = partcount + 1
                                        if (nmin <= partcount .and. partcount <= nmax) then
                                            pc = pc + 1
                                            do l = 1, ndim
                                                call ranf(localseed, xx)
                                                xp(pc, l) = xx * boxlen * dx + xc(l) - boxlen * dx / 2
                                            end do
                                            do kk = 1, nvarh
                                                varp(pc, kk) = vvdp(i, ind, kk)
                                            end do
                                            if (varp(pc, 1) > facdens * averdens) denspartcount = denspartcount + 1
                                        end if
                                    end do
                                end if
                            end if
                        end if

                        if (ndim == 1) then
                            vvdp(i, ind, ivar) = xdp(i)
                            dx = 0.5d0 ** ilevel
                            iz = 0
                            iy = 0
                            ix = (ind - 1)
                            xc(1) = boxlen * (xxdp(i, 1) + (dble(ix) - 0.5D0) * dx - xbound(1))
                            xc(2) = (ymin + ymax) / 2
                            xc(3) = (zmin + zmax) / 2
                            if (j == icpu .and. ivar == nvarh .and. (ilevel >= lmin) .and. (sdp(i, ind) == 0&
                                & .or. ilevel == lmax) .and. (xmin <= xc(1) .and. xc(1) <= xmax) .and.&
                                & (ymin <= xc(2) .and. xc(2) <= ymax) .and. (zmin <= xc(3) .and. xc(3) <= zmax)) then
                            rnpartlocal = (vvdp(i, ind, 1) * dx * dl * dl) / partmass
                            if (rnpartlocal < 1) then
                                call poissdev(localseed, rnpartlocal, poisson)
                                poisson = poisson + nint(masslefttot / dble(ilowdtot) / partmass)
                                jj = 1
                                do while (jj <= poisson .and. partcount < ndummypart)
                                    jj = jj + 1
                                    locind(k) = locind(k) + 1
                                    partcount = partcount + 1
                                    if (nmin <= partcount .and. partcount <= nmax) then
                                        pc = pc + 1
                                        do l = 1, ndim
                                            call ranf(localseed, xx)
                                            xp(pc, l) = xx * boxlen * dx + xc(l) - boxlen * dx / 2
                                        end do
                                        do kk = 1, nvarh
                                            varp(pc, kk) = vvdp(i, ind, kk)
                                        end do
                                        if (varp(pc, 1) > facdens * averdens) denspartcount = denspartcount + 1
                                    end if
                                end do
                            end if
                        end if
                    end if

                end do
            end do
        end do
        istart = istart + ngridfile(j, ilevel)
        deallocate(xdp)
        deallocate(xxdp)
        deallocate(vvdp)
        deallocate(sdp)
        deallocate(idp)
    end if
end do
end do

! End loop over levels
close(10)
close(11)

deallocate(ngridfile, ngridlevel)
if (nboundary > 0) deallocate(ngridbound)

end if

if (indexcell(k) > 0) ncpufull = ncpufull + 1

end do

write(*, *) 'gaspart3 II done'

if (partcount < ndummypart .or. pc < nmax - nmin + 1) then

    respart = int((ndummypart - partcount)) + 1

    allocate(rindex(1:respart))
    allocate(flagcell(1:indexcelltot))

    flagcell(1:indexcelltot) = 0

    do i = 1, respart
        call ranf(localseed, xx)
        rindex(i) = int(xx * dble(indexcelltot - 1)) + 1
        flagcell(rindex(i)) = flagcell(rindex(i)) + 1
    end do

    indexcelltmp = 0
    icell = 0

    do k = 1, ncpu_read

        icpu = cpu_list(k)
        call title(icpu, ncharcpu)

        ! Open AMR file and skip header
        ipos = INDEX(repository,'output_')
        nchar = repository(ipos + 7:ipos + 13)
        nomfich = TRIM(repository) //'/amr_'// TRIM(nchar) //'.out'// TRIM(ncharcpu)
        open(unit=10, file=nomfich, status='old', form='unformatted')
        read(10)
        read(10) ndim
        read(10) nx, ny, nz
        read(10) nlevelmax
        read(10) ngridmax
        read(10) nboundary
        read(10) ngrid_current
        read(10) boxlen
        do i = 1, 13
            read(10)
        end do
        twotondim = 2 ** ndim
        xbound = (/ dble(nx / 2), dble(ny / 2), dble(nz / 2) /)
        allocate(ngridfile(1:ncpu + nboundary, 1:nlevelmax))
        allocate(ngridlevel(1:ncpu, 1:nlevelmax))
        if (nboundary > 0) allocate(ngridbound(1:nboundary, 1:nlevelmax))

        ! Open HYDRO file and skip header
        nomfich = TRIM(repository) //'/hydro_'// TRIM(nchar) //'.out'// TRIM(ncharcpu)
        open(unit=11, file=nomfich, status='old', form='unformatted')
        read(11)
        read(11) nvarh
        read(11)
        read(11)
        read(11)
        read(11) gamma

        ! Read grid numbers
        read(10) ngridlevel
        ngridfile(1:ncpu, 1:nlevelmax) = ngridlevel
        read(10)
        if (nboundary > 0) then
            do i = 1, 2
                read(10)
            end do
            read(10) ngridbound
            ngridfile(ncpu + 1:ncpu + nboundary, 1:nlevelmax) = ngridbound
        end if
        read(10)
        ! ROM: comment the single follwing line for old stuff
        read(10)
        if (TRIM(ordering) == 'bisection') then
            do i = 1, 5
                read(10)
            end do
        else
            read(10)
        end if
        read(10)
        read(10)
        read(10)

        ! Compute total number of grids
        ngridtot = 0
        do j = 1, ncpu
            do ilevel = 1, nlevelmax
                ngridtot = ngridtot + ngridfile(j, ilevel)
            end do
        end do

        ! Loop over levels
        istart = 0
        ! icell=0
        do ilevel = 1, lmax
            ! Loop over domains
            do j = 1, nboundary + ncpu
                ! Read AMR data
                if (ngridfile(j, ilevel) > 0) then
                    allocate(xdp(1:ngridfile(j, ilevel)))
                    allocate(xxdp(1:ngridfile(j, ilevel), 1:ndim))
                    allocate(vvdp(1:ngridfile(j, ilevel), 1:twotondim, 1:nvarh))
                    allocate(sdp(1:ngridfile(j, ilevel), 1:twotondim))
                    allocate(idp(1:ngridfile(j, ilevel)))
                    read(10) ! Skip grid index
                    read(10) ! Skip next index
                    read(10) ! Skip prev index
                    ! Read grid center
                    do ind = 1, ndim
                        read(10) xdp
                        do i = 1, ngridfile(j, ilevel)
                            xxdp(i, ind) = xdp(i)
                        end do
                    end do
                    read(10) ! Skip father index
                    do ind = 1, 2 * ndim
                        read(10) ! Skip nbor index
                    end do
                    ! Read son index
                    do ind = 1, twotondim
                        read(10) idp
                        do i = 1, ngridfile(j, ilevel)
                            sdp(i, ind) = idp(i)
                        end do
                    end do
                    ! Skip cpu map
                    do ind = 1, twotondim
                        read(10)
                    end do
                    ! Skip refinement map
                    do ind = 1, twotondim
                        read(10)
                    end do
                end if
                ! Read HYDRO data
                read(11)
                read(11)
                if (ngridfile(j, ilevel) > 0) then
                    ! Read hydro variables
                    do ind = 1, twotondim
                        do ivar = 1, nvarh
                            read(11) xdp
                            do i = 1, ngridfile(j, ilevel)

                                dl = 0.5d0 ** dble(lmax)

                                if (ndim == 3) then
                                    vvdp(i, ind, ivar) = xdp(i)
                                    dx = 0.5d0 ** ilevel
                                    iz = (ind - 1) / 4
                                    iy = (ind - 1 - 4 * iz) / 2
                                    ix = (ind - 1 - 2 * iy - 4 * iz)
                                    xc(1) = boxlen * (xxdp(i, 1) + (dble(ix) - 0.5D0) * dx - xbound(1))
                                    xc(2) = boxlen * (xxdp(i, 2) + (dble(iy) - 0.5D0) * dx - xbound(2))
                                    xc(3) = boxlen * (xxdp(i, 3) + (dble(iz) - 0.5D0) * dx - xbound(3))
                                    if (j == icpu .and. ivar == nvarh &
                                        & .and. (ilevel >= lmin) &
                                        & .and. (ilevel == lmax &
                                        & .or. sdp(i, ind) == 0) &
                                        & .and. (xmin <= xc(1) &
                                        & .and. xc(1) <= xmax) &
                                        & .and. (ymin <= xc(2) &
                                        & .and. xc(2) <= ymax) &
                                        & .and. (zmin <= xc(3) &
                                        & .and. xc(3) <= zmax)) then
                                    icell = icell + 1
                                    if (flagcell(icell) >= 1 .and. nmin <= partcount .and. partcount < nmax) then
                                        do ii = 1, flagcell(icell)
                                            locind(k) = locind(k) + 1
                                            partcount = partcount + 1
                                            pc = pc + 1
                                            do l = 1, ndim
                                                call ranf(localseed, xx)
                                                xp(pc, l) = xx * boxlen * dx + xc(l) - boxlen * dx / 2
                                            end do
                                            do kk = 1, nvarh
                                                varp(pc, kk) = vvdp(i, ind, kk)
                                            end do
                                            if (varp(pc, 1) > facdens * averdens) denspartcount = denspartcount + 1
                                        end do
                                    end if
                                end if
                            end if

                            if (ndim == 2) then
                                vvdp(i, ind, ivar) = xdp(i)
                                dx = 0.5d0 ** ilevel
                                iz = 0
                                iy = (ind - 1) / 2
                                ix = (ind - 1 - 2 * iy)
                                xc(1) = boxlen * (xxdp(i, 1) + (dble(ix) - 0.5D0) * dx - xbound(1))
                                xc(2) = boxlen * (xxdp(i, 2) + (dble(iy) - 0.5D0) * dx - xbound(2))
                                xc(3) = (zmin + zmax) / 2
                                if (j == icpu .and. ivar == nvarh .and. &
                                    & (ilevel >= lmin) .and. (ilevel == lmax &
                                    & .or. sdp(i, ind) == 0) .and. &
                                    & (xmin <= xc(1) .and. &
                                    & xc(1) <= xmax) .and. &
                                    & (ymin <= xc(2) .and. &
                                    & xc(2) <= ymax) .and. &
                                    & (zmin <= xc(3) .and. &
                                    & xc(3) <= zmax)) then
                                icell = icell + 1
                                if (flagcell(icell) >= 1 .and. nmin <= partcount .and. partcount < nmax) then
                                    do ii = 1, flagcell(icell)
                                        locind(k) = locind(k) + 1
                                        partcount = partcount + 1
                                        pc = pc + 1
                                        do l = 1, ndim
                                            call ranf(localseed, xx)
                                            xp(pc, l) = xx * boxlen * dx + xc(l) - boxlen * dx / 2
                                        end do
                                        do kk = 1, nvarh
                                            varp(pc, kk) = vvdp(i, ind, kk)
                                        end do
                                        if (varp(pc, 1) > facdens * averdens) denspartcount = denspartcount + 1
                                    end do
                                end if
                            end if
                        end if

                        if (ndim == 1) then
                            vvdp(i, ind, ivar) = xdp(i)
                            dx = 0.5d0 ** ilevel
                            iz = 0
                            iy = 0
                            ix = (ind - 1)
                            xc(1) = boxlen * (xxdp(i, 1) + (dble(ix) - 0.5D0) * dx - xbound(1))
                            xc(2) = (ymin + ymax) / 2
                            xc(3) = (zmin + zmax) / 2
                            if (j == icpu .and. ivar == nvarh .and. (ilevel >= lmin) .and.&
                                & (ilevel == lmax .or. sdp(i, ind) == 0) .and.&
                                & (xmin <= xc(1) .and. xc(1) <= xmax) .and.&
                                & (ymin <= xc(2) .and. xc(2) <= ymax) .and.&
                                & (zmin <= xc(3) .and. xc(3) <= zmax)) then
                            icell = icell + 1
                            if (flagcell(icell) >= 1 .and. nmin <= partcount .and.&
                                & partcount < nmax) then
                            do ii = 1, flagcell(icell)
                                locind(k) = locind(k) + 1
                                partcount = partcount + 1
                                pc = pc + 1
                                do l = 1, ndim
                                    call ranf(localseed, xx)
                                    xp(pc, l) = xx * boxlen * dx + xc(l) - boxlen * dx / 2
                                end do
                                do kk = 1, nvarh
                                    varp(pc, kk) = vvdp(i, ind, kk)
                                end do
                                if (varp(pc, 1) > facdens * averdens) &
                                    & denspartcount = denspartcount + 1
                            end do
                        end if
                    end if
                end if

            end do
        end do
    end do
    indexcelltmp = indexcelltmp + indexcell(k)
    istart = istart + ngridfile(j, ilevel)
    deallocate(xdp)
    deallocate(xxdp)
    deallocate(vvdp)
    deallocate(sdp)
    deallocate(idp)
end if
end do
end do
! End loop over levels

close(10)
close(11)

deallocate(ngridfile, ngridlevel)
if (nboundary > 0) deallocate(ngridbound)
end do
deallocate(rindex)
deallocate(flagcell)
end if

write(*, *) 'gaspart3 III done'
! write(*,*)'indexcell', indexcelltot

return

end subroutine gaspart3

! =======================================================================
! =======================================================================
! =======================================================================

subroutine readpart(ncpu, ncpu_read, cpu_list, ndim, repository, metal, star, sink,&
        & lmin, lmax, xmin, xmax, ymin, ymax, zmin, zmax, nmin, nmax, npart_actual, ndm_actual,&
        & nstar_actual, fampart, tagpart)

    implicit none

    integer :: ncpu, ndim, lmax, lmin
    integer :: ipos, i, j, k
    integer :: nmin, nmax, ncpu_read, icpu
    integer, dimension(1:ncpu) :: cpu_list
    real(kind=8) :: xmin, xmax, ymin, ymax, zmin, zmax

    character(LEN=128) :: repository, nomfich
    character(LEN=5) :: ncharcpu, nchar
    real(kind=8), dimension(:), allocatable :: xdp
    integer, dimension(:), allocatable :: idp
    integer(1), dimension(:), allocatable :: ishort
    real(kind=8), dimension(1:3) :: xcc

    integer :: npart_tot, nsink_tot, nstar_tot, npart, nsink, ipa
    integer :: npart_actual, ndm_actual, nstar_actual, npartnow
    integer, dimension(:), allocatable :: levpart, idpart
    integer(1), dimension(:), allocatable :: fampart, tagpart
    real(KIND=8), dimension(:,:), allocatable :: xpart, vpart
    real(KIND=8), dimension(:), allocatable :: mpart, age, met

    logical :: metal, star, sink

    npart_tot = 0
    nsink_tot = 0
    do k = 1, ncpu_read
        icpu = cpu_list(k)
        ! Read number of particles from the Part file
        ipos = INDEX(repository,'output_')
        nchar = repository(ipos + 7:ipos + 13)
        call title(icpu, ncharcpu)
        nomfich = TRIM(repository) //'/part_'// TRIM(nchar) //'.out'// TRIM(ncharcpu)
        open(unit=11, file=nomfich, status='old', form='unformatted')
        read(11)
        read(11)
        read(11) npart
        read(11)
        read(11) nstar_tot
        read(11)
        read(11)
        read(11) nsink
        close(11)
        nsink_tot = nsink_tot + nsink
        npart_tot = npart_tot + npart
    end do

    ! -----------------------------------------------
    ! Initializing Part quantities
    ! -----------------------------------------------

    allocate(xpart(1:npart_tot, 1:ndim))
    allocate(vpart(1:npart_tot, 1:ndim))
    allocate(mpart(1:npart_tot))
    allocate(idpart(1:npart_tot))
    allocate(levpart(1:npart_tot))
    allocate(fampart(1:npart_tot))
    allocate(tagpart(1:npart_tot))
    if (nstar_tot > 0) then
        allocate(age(1:npart_tot))
        if (metal) then
            allocate(met(1:npart_tot))
        end if
    end if

    ! -----------------------------------------------
    ! -----------------------------------------------

    ! ----------------------------------------------------------
    ! Read up variable from the AMR grid and from the Part file
    ! ----------------------------------------------------------

    npartnow = 0

    do k = 1, ncpu_read
        icpu = cpu_list(k)
        call title(icpu, ncharcpu)
        ipos = INDEX(repository,'output_')
        nchar = repository(ipos + 7:ipos + 13)
        ! Open Part file and skip header
        nomfich = TRIM(repository) //'/part_'// TRIM(nchar) //'.out'// TRIM(ncharcpu)
        open(unit=12, file=nomfich, status='old', form='unformatted')
        read(12)
        read(12)
        read(12) npart
        read(12)
        read(12)
        read(12)
        read(12)
        read(12)
        ! Read Part position
        allocate(xdp(1:npart))
        allocate(idp(1:npart))
        allocate(ishort(1:npart))
        do i = 1, ndim
            read(12) xdp
            xpart(npartnow + 1:npartnow + npart, i) = xdp(1:npart)
        end do
        ! Read Part velocity
        do i = 1, ndim
            read(12) xdp
            vpart(npartnow + 1:npartnow + npart, i) = xdp(1:npart)
        end do
        ! Read Part mass, id and level
        read(12) xdp
        mpart(npartnow + 1:npartnow + npart) = xdp(1:npart)
        read(12) idp
        idpart(npartnow + 1:npartnow + npart) = idp(1:npart)
        read(12) idp
        levpart(npartnow + 1:npartnow + npart) = idp(1:npart)
        read(12) ishort
        fampart(npartnow + 1:npartnow + npart) = ishort
        read(12) ishort
        tagpart(npartnow + 1:npartnow + npart) = ishort
        if (sink .or. star) then
            read(12) xdp
            age(npartnow + 1:npartnow + npart) = xdp(1:npart)
            if (metal) then
                read(12) xdp
                met(npartnow + 1:npartnow + npart) = xdp(1:npart)
            end if
        end if
        close(12)
        npartnow = npartnow + npart
        deallocate(xdp)
        deallocate(idp)

    end do
    ! End loop over cpu

    ! -----------------------------------------------------------------------
    ! Counting Particles an writing output arrays
    ! -----------------------------------------------------------------------

    npart_actual = 0
    nstar_actual = 0
    ndm_actual = 0
    do i = 1, npart_tot
        if (ndim == 3) then
            xcc(1) = xpart(i, 1)
            xcc(2) = xpart(i, 2)
            xcc(3) = xpart(i, 3)
        end if
        if (ndim == 2) then
            xcc(1) = xpart(i, 1)
            xcc(2) = xpart(i, 2)
            xcc(3) = (zmin + zmax) / 2
        end if
        if (ndim == 1) then
            xcc(1) = xpart(i, 1)
            xcc(2) = (ymin + ymax) / 2
            xcc(3) = (zmin + zmax) / 2
        end if
        ! if(lmin<=levpart(i).and.levpart(i)<=lmax.and.xmin<=xcc(1).and.xcc(1)<=xmax&
        if (         xmin <= xcc(1) .and. xcc(1) <= xmax&
            & .and. ymin <= xcc(2) .and. xcc(2) <= ymax&
            & .and. zmin <= xcc(3) .and. xcc(3) <= zmax&
            & .and. idpart(i) > 0) then
        npart_actual = npart_actual + 1 ! NOTE: only particles with ID>0. TO SELECT ALSO SINK PARTICLES USE READPART2
    end if

    if ((.not. star) .and. (.not. sink)) then
        if (         xmin <= xcc(1) .and. xcc(1) <= xmax&
            & .and. ymin <= xcc(2) .and. xcc(2) <= ymax&
            & .and. zmin <= xcc(3) .and. xcc(3) <= zmax) then
        ndm_actual = ndm_actual + 1
    end if
end if

if ((.not. star) .and. (sink)) then
    if (idpart(i) > 0&
        & .and. xmin <= xcc(1) .and. xcc(1) <= xmax&
        & .and. ymin <= xcc(2) .and. xcc(2) <= ymax&
        & .and. zmin <= xcc(3) .and. xcc(3) <= zmax) then
    ndm_actual = ndm_actual + 1
end if
end if

if (star .or. sink) then
    if (age(i) == 0d0 .and. idpart(i) > 0&
        & .and. xmin <= xcc(1) .and. xcc(1) <= xmax&
        & .and. ymin <= xcc(2) .and. xcc(2) <= ymax&
        & .and. zmin <= xcc(3) .and. xcc(3) <= zmax) then
    ndm_actual = ndm_actual + 1
end if
end if

if (star .or. sink) then
    if (age(i) /= 0d0 .and. idpart(i) > 0&
        & .and. xmin <= xcc(1) .and. xcc(1) <= xmax&
        & .and. ymin <= xcc(2) .and. xcc(2) <= ymax&
        & .and. zmin <= xcc(3) .and. xcc(3) <= zmax) then
    nstar_actual = nstar_actual + 1
end if
end if
end do

allocate(xout(1:npart_actual, 1:ndim))
allocate(vout(1:npart_actual, 1:ndim))
allocate(mout(1:npart_actual))
allocate(idout(1:npart_actual))
allocate(ageout(1:npart_actual))
allocate(metout(1:npart_actual))
ageout = 0d0
metout = 0d0

ipa = 0
do i = 1, npart_tot
    if (ndim == 3) then
        xcc(1) = xpart(i, 1)
        xcc(2) = xpart(i, 2)
        xcc(3) = xpart(i, 3)
    end if
    if (ndim == 2) then
        xcc(1) = xpart(i, 1)
        xcc(2) = xpart(i, 2)
        xcc(3) = (zmin + zmax) / 2
    end if
    if (ndim == 1) then
        xcc(1) = xpart(i, 1)
        xcc(2) = (ymin + ymax) / 2
        xcc(3) = (zmin + zmax) / 2
    end if
    if (idpart(i) > 0&
        & .and. xmin <= xcc(1) .and. xcc(1) <= xmax&
        & .and. ymin <= xcc(2) .and. xcc(2) <= ymax&
        & .and. zmin <= xcc(3) .and. xcc(3) <= zmax) then
    ipa = ipa + 1
    do j = 1, ndim
        xout(ipa, j) = xpart(i, j)
        vout(ipa, j) = vpart(i, j)
    end do
    mout(ipa) = mpart(i)
    idout(ipa) = idpart(i)
    if (nstar_actual > 0) then
        ageout(ipa) = age(i)
        if (metal) then
            metout(ipa) = met(i)
        end if
    end if
end if
end do

deallocate(xpart, vpart, mpart, idpart, levpart)
if (nstar_tot > 0) then
    deallocate(age)
    if (metal) deallocate(met)
end if

return

end subroutine readpart

! =======================================================================
! =======================================================================
! =======================================================================


subroutine readpart2(ncpu, ncpu_read, cpu_list, ndim, repository, metal, star, sink,&
        & lmin, lmax, xmin, xmax, ymin, ymax, zmin, zmax, nmin, nmax, npart_actual, ndm_actual,&
        & nstar_actual)

    ! DAVIDE MARTIZZI 2010

    implicit none

    integer :: ncpu, ndim, lmax, lmin
    integer :: ipos, i, j, k
    integer :: nmin, nmax, ncpu_read, icpu
    integer, dimension(1:ncpu) :: cpu_list
    real(kind=8) :: xmin, xmax, ymin, ymax, zmin, zmax

    character(LEN=128) :: repository, nomfich
    character(LEN=5) :: ncharcpu, nchar
    real(kind=8), dimension(:), allocatable :: xdp
    integer, dimension(:)  , allocatable :: idp
    real(kind=8), dimension(1:3) :: xcc

    integer :: npart_tot, nsink_tot, nstar_tot, npart, nsink, ipa
    integer :: npart_actual, ndm_actual, nstar_actual, npartnow
    integer, dimension(:), allocatable :: levpart, idpart
    real(KIND=8), dimension(:,:), allocatable :: xpart, vpart
    real(KIND=8), dimension(:), allocatable :: mpart, age, met

    logical :: metal, star, sink


    npart_tot = 0
    nsink_tot = 0
    do k = 1, ncpu_read
        icpu = cpu_list(k)
        ! Read number of particles from the Part file
        ipos = INDEX(repository,'output_')
        nchar = repository(ipos + 7:ipos + 13)
        call title(icpu, ncharcpu)
        nomfich = TRIM(repository) //'/part_'// TRIM(nchar) //'.out'// TRIM(ncharcpu)
        open(unit=11, file=nomfich, status='old', form='unformatted')
        read(11)
        read(11)
        read(11) npart
        read(11)
        read(11) nstar_tot
        read(11)
        read(11)
        read(11) nsink
        close(11)
        nsink_tot = nsink_tot + nsink
        npart_tot = npart_tot + npart
    end do

    ! -----------------------------------------------
    ! Initializing Part quantities
    ! -----------------------------------------------

    allocate(xpart(1:npart_tot, 1:ndim))
    allocate(vpart(1:npart_tot, 1:ndim))
    allocate(mpart(1:npart_tot))
    allocate(idpart(1:npart_tot))
    allocate(levpart(1:npart_tot))
    if (nstar_tot > 0) then
        allocate(age(1:npart_tot))
        if (metal) then
            allocate(met(1:npart_tot))
        end if
    end if

    ! -----------------------------------------------
    ! -----------------------------------------------

    ! ----------------------------------------------------------
    ! Read up variable from the AMR grid and from the Part file
    ! ----------------------------------------------------------

    npartnow = 0

    do k = 1, ncpu_read
        icpu = cpu_list(k)
        call title(icpu, ncharcpu)
        ipos = INDEX(repository,'output_')
        nchar = repository(ipos + 7:ipos + 13)
        ! Open Part file and skip header
        nomfich = TRIM(repository) //'/part_'// TRIM(nchar) //'.out'// TRIM(ncharcpu)
        open(unit=12, file=nomfich, status='old', form='unformatted')
        read(12)
        read(12)
        read(12) npart
        read(12)
        read(12)
        read(12)
        read(12)
        read(12)
        ! Read Part position
        allocate(xdp(1:npart))
        allocate(idp(1:npart))
        do i = 1, ndim
            read(12) xdp
            xpart(npartnow + 1:npartnow + npart, i) = xdp(1:npart)
        end do
        ! Read Part velocity
        do i = 1, ndim
            read(12) xdp
            vpart(npartnow + 1:npartnow + npart, i) = xdp(1:npart)
        end do
        ! Read Part mass, id and level
        read(12) xdp
        mpart(npartnow + 1:npartnow + npart) = xdp(1:npart)
        read(12) idp
        idpart(npartnow + 1:npartnow + npart) = idp(1:npart)
        read(12) idp
        levpart(npartnow + 1:npartnow + npart) = idp(1:npart)
        if (sink .or. star) then
            read(12) xdp
            age(npartnow + 1:npartnow + npart) = xdp(1:npart)
            if (metal) then
                read(12) xdp
                met(npartnow + 1:npartnow + npart) = xdp(1:npart)
            end if
        end if
        close(12)
        npartnow = npartnow + npart
        deallocate(xdp)
        deallocate(idp)

    end do
    ! End loop over cpu

    ! -----------------------------------------------------------------------
    ! Counting Particles an writing output arrays
    ! -----------------------------------------------------------------------

    npart_actual = 0
    nstar_actual = 0
    ndm_actual = 0
    do i = 1, npart_tot
        if (ndim == 3) then
            xcc(1) = xpart(i, 1)
            xcc(2) = xpart(i, 2)
            xcc(3) = xpart(i, 3)
        end if
        if (ndim == 2) then
            xcc(1) = xpart(i, 1)
            xcc(2) = xpart(i, 2)
            xcc(3) = (zmin + zmax) / 2
        end if
        if (ndim == 1) then
            xcc(1) = xpart(i, 1)
            xcc(2) = (ymin + ymax) / 2
            xcc(3) = (zmin + zmax) / 2
        end if

        if (lmin <= levpart(i) .and. levpart(i) <= lmax .and. xmin <= xcc(1) .and. xcc(1) <= xmax &
            & .and. ymin <= xcc(2) .and. xcc(2) <= ymax .and. zmin <= xcc(3) .and. xcc(3) <= zmax)&
            npart_actual = npart_actual + 1 ! NOTE: SELECT ALSO SINK PARTICLES!!!

        if ((.not. star) .and. (.not. sink)) then
            if (lmin <= levpart(i) .and. levpart(i) <= lmax .and. xmin <= xcc(1)&
                & .and. xcc(1) <= xmax .and. ymin <= xcc(2) .and. xcc(2) <= ymax .and. zmin <= xcc(3)&
                & .and. xcc(3) <= zmax) ndm_actual = ndm_actual + 1
        end if

        if ((.not. star) .and. (sink)) then
            if (idpart(i) > 0 .and. lmin <= levpart(i) .and. levpart(i) <= lmax .and.&
                & xmin <= xcc(1) .and. xcc(1) <= xmax .and. ymin <= xcc(2) .and. xcc(2) <= ymax&
                & .and. zmin <= xcc(3) .and. xcc(3) <= zmax) ndm_actual = ndm_actual + 1
        end if

        if (star .or. sink) then
            if (age(i) == 0d0 .and. idpart(i) > 0 .and. lmin <= levpart(i) .and. levpart(i) <= lmax .and.&
                & xmin <= xcc(1) .and. xcc(1) <= xmax .and. ymin <= xcc(2) .and. xcc(2) <= ymax&
                & .and. zmin <= xcc(3) .and. xcc(3) <= zmax) ndm_actual = ndm_actual + 1
        end if

        if (star .or. sink) then
            if (age(i) /= 0d0 .and. idpart(i) > 0 .and. lmin <= levpart(i) .and. levpart(i) <= lmax .and.&
                & xmin <= xcc(1) .and. xcc(1) <= xmax .and. ymin <= xcc(2) .and. xcc(2) <= ymax&
                & .and. zmin <= xcc(3) .and. xcc(3) <= zmax) nstar_actual = nstar_actual + 1
        end if
    end do

    allocate(xout(1:npart_actual, 1:ndim))
    allocate(vout(1:npart_actual, 1:ndim))
    allocate(mout(1:npart_actual))
    allocate(idout(1:npart_actual))
    allocate(ageout(1:npart_actual))
    allocate(metout(1:npart_actual))
    ageout = 0d0
    metout = 0d0

    ipa = 0
    do i = 1, npart_tot
        if (ndim == 3) then
            xcc(1) = xpart(i, 1)
            xcc(2) = xpart(i, 2)
            xcc(3) = xpart(i, 3)
        end if
        if (ndim == 2) then
            xcc(1) = xpart(i, 1)
            xcc(2) = xpart(i, 2)
            xcc(3) = (zmin + zmax) / 2
        end if
        if (ndim == 1) then
            xcc(1) = xpart(i, 1)
            xcc(2) = (ymin + ymax) / 2
            xcc(3) = (zmin + zmax) / 2
        end if
        if ((lmin <= levpart(i) .and. levpart(i) <= lmax) .and. (xmin <= xcc(1) .and. xcc(1) <= xmax)&
            & .and. (ymin <= xcc(2) .and. xcc(2) <= ymax) .and. (zmin <= xcc(3) .and. xcc(3) <= zmax)) then
        ipa = ipa + 1
        do j = 1, ndim
            xout(ipa, j) = xpart(i, j)
            vout(ipa, j) = vpart(i, j)
        end do
        mout(ipa) = mpart(i)
        idout(ipa) = idpart(i)
        if (nstar_actual > 0) then
            ageout(ipa) = age(i)
            if (metal) then
                metout(ipa) = met(i)
            end if
        end if
    end if
end do

deallocate(xpart, vpart, mpart, idpart, levpart)
if (nstar_tot > 0) then
    deallocate(age)
    if (metal) deallocate(met)
end if

return

end subroutine readpart2

! =======================================================================
! =======================================================================
! =======================================================================

end module io_ramses
