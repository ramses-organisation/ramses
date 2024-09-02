program icrefine
! --------------------------------------------------------------------------
! Ce programme calcule la carte de densite surfacique projetee
! des particules de matiere noire d'une simulation RAMSES.
! Version F90 par R. Teyssier le 01/04/01.
! --------------------------------------------------------------------------
implicit none
integer :: ncpu, ndim, npart, ngrid, n, i, j, k, icpu, ipos, nstar, nstart, inull
integer :: ncpu2, npart2, ndim2, levelmin, levelmax, ilevel, ndark, ismooth
integer :: nx=0, ny=0, ix, iy, iz, ixp1, iyp1, idim, jdim, ncpu_read, smt=0, rsm=1
real(KIND=8) :: mtot, ddx, ddy, dex, dey, t, soft, poty, mass, btime, unit_l, aexp, unit_t
real(KIND=8) :: metal=- 1
real(KIND=8) :: xmin=0, xmax=1, ymin=0, ymax=1, zmin=0, zmax=1, r, xc=0.5, yc=0.5, zc=0.5, rad=- 1
integer :: imin, imax, jmin, jmax, kmin, kmax, lmin, ipart
real(KIND=8) :: xxmin, xxmax, yymin, yymax, dy, deltax, fakeage
real(KIND=4), dimension(:,:), allocatable :: toto
real(KIND=8), dimension(:,:), allocatable :: map
real(KIND=8), dimension(:)  , allocatable :: x
real(KIND=8), dimension(:)  , allocatable :: y
real(KIND=8), dimension(:)  , allocatable :: z
real(KIND=8), dimension(:)  , allocatable :: vx
real(KIND=8), dimension(:)  , allocatable :: vy
real(KIND=8), dimension(:)  , allocatable :: vz
real(KIND=8), dimension(:)  , allocatable :: m
real(KIND=8), dimension(:)  , allocatable :: temp
real(KIND=8), dimension(:)  , allocatable :: met
real(KIND=8), dimension(:)  , allocatable :: bt
integer , allocatable, dimension(:) :: temp2
integer , allocatable, dimension(:) :: id
integer , allocatable, dimension(:) :: idpart
real , allocatable, dimension(:,:,:) :: imark
character(LEN=1) :: proj='z'
character(LEN=5) :: nchar, ncharcpu
character(LEN=80) :: ordering, format_grille
character(LEN=80) :: GMGM
character(LEN=128) :: nomfich, repository, filetype='bin', grafic
logical :: ok, ok_part, periodic=.false., star=.false., okerode=.false.
logical :: gid=.false., fid=.false.
integer :: impi, ndom, bit_length, maxdom, maxid, idd
integer, dimension(1:8) :: idom, jdom, kdom, cpu_min, cpu_max
real(KIND=8), dimension(1:8) :: bounding_min, bounding_max
real(KIND=8) :: dkey, order_min, dmax, vfact
real(kind=8), dimension(:), allocatable :: bound_key
logical, dimension(:), allocatable :: cpu_read
integer, dimension(:), allocatable :: cpu_list
integer(kind=4) :: np1, np2, np3
real :: dx, dx2, x1o, x2o, x3o, astart, omegam, omegav, h0, x1or, x2or, x3or, dxor, omegak

call read_params

! -----------------------------------------------
! Lecture du fichier particules au format RAMSES
! -----------------------------------------------
ipos = INDEX(repository,'output_')
nchar = repository(ipos + 7:ipos + 13)
nomfich = TRIM(repository) //'/part_'// TRIM(nchar) //'.out00001'
inquire(file = nomfich, exist = ok) ! verify input file
if ( .not. ok ) then
    print *, TRIM(nomfich) //' not found.'
    stop
end if

nomfich = TRIM(repository) //'/info_'// TRIM(nchar) //'.txt'
inquire(file = nomfich, exist = ok) ! verify input file
if ( .not. ok ) then
    print *, TRIM(nomfich) //' not found.'
    stop
end if
open(unit=10, file=nomfich, form='formatted', status='old')
read(10,'("ncpu        =",I11)') ncpu
read(10,'("ndim        =",I11)') ndim
read(10,'("levelmin    =",I11)') levelmin
read(10,'("levelmax    =",I11)') levelmax
read(10,*)
read(10,*)
read(10,*)
write(*, *) ncpu, ndim, levelmin, levelmax

read(10,*)
read(10,'("time        =",E23.15)') t
read(10,'("aexp        =",E23.15)') aexp
read(10,*)
read(10,*)
read(10,*)
read(10,*)
read(10,*)
read(10,'("unit_l      =",E23.15)') unit_l
read(10,*)
read(10,'("unit_t      =",E23.15)') unit_t

read(10,*)
read(10,'("ordering type=",A80)'), ordering
read(10,*)
write(*, '(" ordering type=",A20)') , TRIM(ordering)
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

if (rad > 0) then
    xmin = xc - rad
    xmax = xc + rad
    ymin = yc - rad
    ymax = yc + rad
    zmin = zc - rad
    zmax = zc + rad
end if

if (TRIM(ordering) == 'hilbert') then

    dmax = max(xmax - xmin, ymax - ymin, zmax - zmin)
    do ilevel = 1, levelmax
        deltax = 0.5d0 ** ilevel
        if (deltax < dmax) exit
    end do
    lmin = ilevel
    bit_length = lmin - 1
    maxdom = 2 ** bit_length
    imin = 0; imax = 0; jmin = 0; jmax = 0; kmin = 0; kmax = 0
    if (bit_length > 0) then
        imin = int(xmin * dble(maxdom))
        imax = imin + 1
        jmin = int(ymin * dble(maxdom))
        jmax = jmin + 1
        kmin = int(zmin * dble(maxdom))
        kmax = kmin + 1
    end if

    dkey = (dble(2 ** (levelmax + 1) / dble(maxdom))) ** ndim
    ndom = 1
    if (bit_length > 0) ndom = 8
    idom(1) = imin; idom(2) = imax
    idom(3) = imin; idom(4) = imax
    idom(5) = imin; idom(6) = imax
    idom(7) = imin; idom(8) = imax
    jdom(1) = jmin; jdom(2) = jmin
    jdom(3) = jmax; jdom(4) = jmax
    jdom(5) = jmin; jdom(6) = jmin
    jdom(7) = jmax; jdom(8) = jmax
    kdom(1) = kmin; kdom(2) = kmin
    kdom(3) = kmin; kdom(4) = kmin
    kdom(5) = kmax; kdom(6) = kmax
    kdom(7) = kmax; kdom(8) = kmax

    do i = 1, ndom
        if (bit_length > 0) then
            call hilbert3d(idom(i), jdom(i), kdom(i), order_min, bit_length, 1)
        else
            order_min = 0.0d0
        end if
        bounding_min(i) = (order_min) * dkey
        bounding_max(i) = (order_min + 1.0D0) * dkey
    end do
    cpu_min = 0; cpu_max = 0
    do impi = 1, ncpu
        do i = 1, ndom
            if (   bound_key(impi - 1) <= bounding_min(i) .and.&
                & bound_key(impi  ) > bounding_min(i)) then
            cpu_min(i) = impi
        end if
        if (   bound_key(impi - 1) < bounding_max(i) .and.&
            & bound_key(impi  ) >= bounding_max(i)) then
        cpu_max(i) = impi
    end if
end do
end do

ncpu_read = 0
do i = 1, ndom
    do j = cpu_min(i), cpu_max(i)
        if (.not. cpu_read(j)) then
            ncpu_read = ncpu_read + 1
            cpu_list(ncpu_read) = j
            cpu_read(j) = .true.
        end if
    end do
end do
else
ncpu_read = ncpu
do j = 1, ncpu
    cpu_list(j) = j
end do
end  if

npart = 0
do k = 1, ncpu_read
    write(*, *) 'CPU=', k
    icpu = cpu_list(k)
    call title(icpu, ncharcpu)
    nomfich = TRIM(repository) //'/part_'// TRIM(nchar) //'.out'// TRIM(ncharcpu)
    open(unit=1, file=nomfich, status='old', form='unformatted')
    read(1) ncpu2
    read(1) ndim2
    read(1) npart2
    read(1)
    read(1) nstar
    close(1)
    npart = npart + npart2
end do
write(*, *) npart,' particles in the region'
allocate(m(1:npart))
allocate(x(1:npart))
allocate(y(1:npart))
allocate(z(1:npart))
allocate(vx(1:npart))
allocate(vy(1:npart))
allocate(vz(1:npart))
allocate(id(1:npart))
if (nstar > 0) then
    allocate(bt(1:npart))
end if

! -----------------------------------------------
! Compute projected mass using CIC smoothing
! ----------------------------------------------
mtot = 0.0d0
nstart = 1
do k = 1, ncpu_read
    icpu = cpu_list(k)
    call title(icpu, ncharcpu)
    nomfich = TRIM(repository) //'/part_'// TRIM(nchar) //'.out'// TRIM(ncharcpu)
    open(unit=1, file=nomfich, status='old', form='unformatted')
    write(*, *) 'Processing file '// TRIM(nomfich)
    read(1) ncpu2
    read(1) ndim2
    read(1) npart2
    read(1)
    read(1)
    read(1)
    read(1)
    read(1)
    allocate(temp(1:npart2))
    allocate(temp2(1:npart2))
    ! Read positions
    read(1) temp
    x(nstart:nstart + npart2 - 1) = temp
    read(1) temp
    y(nstart:nstart + npart2 - 1) = temp
    read(1) temp
    z(nstart:nstart + npart2 - 1) = temp
    ! Read velocity
    read(1) temp
    vx(nstart:nstart + npart2 - 1) = temp
    read(1) temp
    vy(nstart:nstart + npart2 - 1) = temp
    read(1) temp
    vz(nstart:nstart + npart2 - 1) = temp
    ! Read mass
    read(1) temp
    m(nstart:nstart + npart2 - 1) = temp
    ! Read identity
    read(1) temp2
    id(nstart:nstart + npart2 - 1) = temp2
    ! Read level
    read(1) temp2
    if (nstar > 0) then
        ! Read BT
        read(1) temp
        bt(nstart:nstart + npart2 - 1) = temp
    end if
    ! ----------------------------
    nstart = nstart + npart2  ! Fill up the next set
    deallocate(temp)
    deallocate(temp2)
end do

40 format(3e16.8)
50 format(2I16)
! Outputs IDs of selected particles
ipart = 0
mass = 0.0
write(*, *) 'Getting IDs...'
open(18, file='partID.dat', form='formatted')
maxid = 0
do i = 1, npart ! To get maximum identity of the particle
    if (nstar == 0) then  ! Only DM particles
        btime = 0
    else
        btime = bt(i)
    end if
    if (btime == 0) then
        ok_part = (x(i) >= xmin .and. x(i) <= xmax .and. &
            &   y(i) >= ymin .and. y(i) <= ymax .and. &
            &   z(i) >= zmin .and. z(i) <= zmax)
        if (rad > 0) then
            r = (x(i) - xc) ** 2 + (y(i) - yc) ** 2 + (z(i) - zc) ** 2
            ok_part = (sqrt(r) <= rad)
        end if
        if (ok_part) then
            maxid = max(maxid, id(i))
            ipart = ipart + 1
            mass = mass + m(i)
        end if
    end if
end do
write(*, *) 'We have', ipart,' particles in selected region'
write(*, *) 'Total mass =', mass

30 format(i16)
write(18, 50) ipart, npart, maxid
do i = 1, npart  ! Start finding the IDs
    if (nstar == 0) then  ! Only DM particles
        btime = 0
    else
        btime = bt(i)
    end if
    if (btime == 0) then
        ok_part = (x(i) >= xmin .and. x(i) <= xmax .and. &
            &   y(i) >= ymin .and. y(i) <= ymax .and. &
            &   z(i) >= zmin .and. z(i) <= zmax)
        if (rad > 0) then
            r = (x(i) - xc) ** 2 + (y(i) - yc) ** 2 + (z(i) - zc) ** 2
            ok_part = sqrt(r) <= rad
        end if
        if (ok_part) then
            write(18, 30) id(i)   ! Write IDs
        end if
    end if
end do

contains

subroutine read_params

    implicit none

    integer       :: i, n

    character(len=4)   :: opt
    character(len=128) :: arg
    logical       :: bad, ok

    n = command_argument_count()
    if (n < 4) then
        print *, 'usage: geticref  -inp  input_dir'
        print *, '                 [-dir axis] '
        print *, '                 [-xc xc] '
        print *, '                 [-yc yc] '
        print *, '                 [-zc zc] '
        print *, '                 [-rad rad] '
        print *, 'ex: geticref -inp output_00001 -xc 0.5 -yc 0.5 -zc 0.5 -rad 0.1'
        stop
    end if

    do i = 1, n, 2
        call get_command_argument(i, opt)
        if (i == n) then
            print '("option ",a2," has no argument")', opt
            stop 2
        end if
        call get_command_argument(i + 1, arg)
        select case (opt)
            case ('-inp')
            repository = trim(arg)
            case ('-dir')
            proj = trim(arg)
            case ('-xmi')
            read (arg,*) xmin
            case ('-xma')
            read (arg,*) xmax
            case ('-ymi')
            read (arg,*) ymin
            case ('-yma')
            read (arg,*) ymax
            case ('-zmi')
            read (arg,*) zmin
            case ('-zma')
            read (arg,*) zmax
            case ('-xc')
            read (arg,*) xc
            case ('-yc')
            read (arg,*) yc
            case ('-zc')
            read (arg,*) zc
            case ('-rad')
            read (arg,*) rad
            case ('-per')
            read (arg,*) periodic
            case ('-gid')
            read (arg,*) gid
            case ('-fid')
            read (arg,*) fid
            case ('-smt')
            read (arg,*) smt
            case ('-met')
            read (arg,*) metal
            case ('-rsm')
            read (arg,*) rsm
            case ('-gfc')
            grafic = trim(arg)
            case ('-fil')
            filetype = trim(arg)
            case default
            print '("unknown option ",a2," ignored")', opt
        end select
    end do

    return

end subroutine read_params
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
function fy(a)
    implicit none
    ! Computes the integrand
    real(kind=8) :: fy
    real(kind=8) :: y, a

    y = omegam * (1d0 / a - 1d0) + omegav * (a * a - 1d0) + 1d0
    fy = 1d0 / y ** 1.5d0

    return
end function fy

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
function fpeebl(a)
    implicit none
    real(kind=8) :: fpeebl, a
    ! Computes the growth factor f=d\log D1/d\log a.
    real(kind=8) :: fact, y, eps

    eps = 1.0d-6
    y = omegam * (1d0 / a - 1d0) + omegav * (a * a - 1d0) + 1d0
    fact = rombint(eps, a, eps)
    fpeebl = (omegav * a * a - 0.5d0 * omegam / a) / y - 1d0 + a * fy(a) / fact
    return
end function fpeebl
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
function rombint(a, b, tol)
    implicit none
    real(kind=8) :: rombint
    !
    ! Rombint returns the integral from a to b of f(x)dx using Romberg
    ! integration. The method converges provided that f(x) is continuous
    ! in (a,b). The function f must be double precision and must be
    ! declared external in the calling routine.
    ! tol indicates the desired relative accuracy in the integral.
    !
    integer :: maxiter=16, maxj=5
    real(kind=8), dimension(100) :: g
    real(kind=8) :: a, b, tol, fourj
    real(kind=8) :: h, error, gmax, g0, g1
    integer :: nint, i, j, k, jmax

    h = 0.5d0 * (b - a)
    gmax = h * (fy(a) + fy(b))
    g(1) = gmax
    nint = 1
    error = 1.0d20
    i = 0
10  i = i + 1
    if (.not.  (i > maxiter .or. (i > 5 .and. abs(error) < tol))) then
        ! Calculate next trapezoidal rule approximation to integral.

        g0 = 0.0d0
        do k = 1, nint
            g0 = g0 + fy(a + (k + k - 1) * h)
        end do
        g0 = 0.5d0 * g(1) + h * g0
        h = 0.5d0 * h
        nint = nint + nint
        jmax = min(i, maxj)
        fourj = 1.0d0

        do j = 1, jmax
            ! Use Richardson extrapolation.
            fourj = 4.0d0 * fourj
            g1 = g0 + (g0 - g(j)) / (fourj - 1.0d0)
            g(j) = g0
            g0 = g1
        end do
        if (abs(g0) > tol) then
            error = 1.0d0 - gmax / g0
        else
            error = gmax
        end if
        gmax = g0
        g(jmax + 1) = g0
        go to 10
    end if
    rombint = g0
    if (i > maxiter .and. abs(error) > tol) &
        &    write(*, *) 'Rombint failed to converge; integral, error=', &
        &    rombint, error
    return
end function rombint
end program icrefine

! =======================================================================
subroutine title(n, nchar)
    ! =======================================================================
    implicit none
    integer :: n
    character(5) :: nchar

    character(1) :: nchar1
    character(2) :: nchar2
    character(3) :: nchar3
    character(4) :: nchar4
    character(5) :: nchar5

    if (n >= 10000) then
        write(nchar5, '(i5)') n
        nchar = nchar5
    elseif (n >= 1000) then
        write(nchar4, '(i4)') n
        nchar = '0'// nchar4
    elseif (n >= 100) then
        write(nchar3, '(i3)') n
        nchar = '00'// nchar3
    elseif (n >= 10) then
        write(nchar2, '(i2)') n
        nchar = '000'// nchar2
    else
        write(nchar1, '(i1)') n
        nchar = '0000'// nchar1
    end if

end subroutine title

! ================================================================
! ================================================================
! ================================================================
! ================================================================
subroutine hilbert3d(x, y, z, order, bit_length, npoint)
    implicit none

    integer     , INTENT(IN)                     :: bit_length, npoint
    integer     , INTENT(IN) , dimension(1:npoint) :: x, y, z
    real(kind=8), INTENT(OUT), dimension(1:npoint) :: order

    logical, dimension(0:3 * bit_length - 1) :: i_bit_mask
    logical, dimension(0:1 * bit_length - 1) :: x_bit_mask, y_bit_mask, z_bit_mask
    integer, dimension(0:7, 0:1, 0:11) :: state_diagram
    integer :: i, ip, cstate, nstate, b0, b1, b2, sdigit, hdigit

    if (bit_length > bit_size(bit_length)) then
        write(*, *) 'Maximum bit length=', bit_size(bit_length)
        write(*, *) 'stop in hilbert3d'
        stop
    end if

    state_diagram = RESHAPE( (/   1, 2, 3, 2, 4, 5, 3, 5,&
        &   0, 1, 3, 2, 7, 6, 4, 5,&
        &   2, 6, 0, 7, 8, 8, 0, 7,&
        &   0, 7, 1, 6, 3, 4, 2, 5,&
        &   0, 9, 10, 9, 1, 1, 11, 11,&
        &   0, 3, 7, 4, 1, 2, 6, 5,&
        &   6, 0, 6, 11, 9, 0, 9, 8,&
        &   2, 3, 1, 0, 5, 4, 6, 7,&
        &  11, 11, 0, 7, 5, 9, 0, 7,&
        &   4, 3, 5, 2, 7, 0, 6, 1,&
        &   4, 4, 8, 8, 0, 6, 10, 6,&
        &   6, 5, 1, 2, 7, 4, 0, 3,&
        &   5, 7, 5, 3, 1, 1, 11, 11,&
        &   4, 7, 3, 0, 5, 6, 2, 1,&
        &   6, 1, 6, 10, 9, 4, 9, 10,&
        &   6, 7, 5, 4, 1, 0, 2, 3,&
        &  10, 3, 1, 1, 10, 3, 5, 9,&
        &   2, 5, 3, 4, 1, 6, 0, 7,&
        &   4, 4, 8, 8, 2, 7, 2, 3,&
        &   2, 1, 5, 6, 3, 0, 4, 7,&
        &   7, 2, 11, 2, 7, 5, 8, 5,&
        &   4, 5, 7, 6, 3, 2, 0, 1,&
        &  10, 3, 2, 6, 10, 3, 4, 4,&
        &   6, 1, 7, 0, 5, 2, 4, 3 /), &
        & (/ 8 , 2, 12 /) )

    do ip = 1, npoint

        ! convert to binary
        do i = 0, bit_length - 1
            x_bit_mask(i) = btest(x(ip), i)
            y_bit_mask(i) = btest(y(ip), i)
            z_bit_mask(i) = btest(z(ip), i)
        end do

        ! interleave bits
        do i = 0, bit_length - 1
            i_bit_mask(3 * i + 2) = x_bit_mask(i)
            i_bit_mask(3 * i + 1) = y_bit_mask(i)
            i_bit_mask(3 * i  ) = z_bit_mask(i)
        end do

        ! build Hilbert ordering using state diagram
        cstate = 0
        do i = bit_length - 1, 0, - 1
            b2 = 0 ; if (i_bit_mask(3 * i + 2)) b2 = 1
            b1 = 0 ; if (i_bit_mask(3 * i + 1)) b1 = 1
            b0 = 0 ; if (i_bit_mask(3 * i  )) b0 = 1
            sdigit = b2 * 4 + b1 * 2 + b0
            nstate = state_diagram(sdigit, 0, cstate)
            hdigit = state_diagram(sdigit, 1, cstate)
            i_bit_mask(3 * i + 2) = btest(hdigit, 2)
            i_bit_mask(3 * i + 1) = btest(hdigit, 1)
            i_bit_mask(3 * i  ) = btest(hdigit, 0)
            cstate = nstate
        end do

        ! save Hilbert key as double precision real
        order(ip) = 0.
        do i = 0, 3 * bit_length - 1
            b0 = 0 ; if (i_bit_mask(i)) b0 = 1
            order(ip) = order(ip) + dble(b0) * dble(2) ** i
        end do

    end do

end subroutine hilbert3d
