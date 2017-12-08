program part2list
  !--------------------------------------------------------------------------
  ! This software extracts particle from a region.
  ! Version F90 by R. Teyssier le 01/04/01.
  ! Update for new tracers by C. Cadiou, M. Trebitsch, H. Choi, O. Snaith and R. Bieri
  !--------------------------------------------------------------------------
  use utils
  use iso_fortran_env
  implicit none
  integer :: ncpu, ndim, npart, i, j, k, icpu, ipos, n_frw, nstar
  integer :: ncpu2, npart2, ndim2, levelmin, levelmax, ilevel, iii
  integer :: nx = 0, ny = 0, nz = 0, ncpu_read
  real(real64) :: mtot, t, time, time_tot, time_simu, weight
  real(real64) :: xmin = 0, xmax = 1, ymin = 0, ymax = 1, zmin = 0, zmax = 1
  real(real64) :: aexp, omega_m, omega_l, omega_b, omega_k, h0, unit_l, unit_t, unit_d
  integer :: imin, imax, jmin, jmax, kmin, kmax, lmin
  real(real64) :: deltax
  real(real64), dimension(:), allocatable :: aexp_frw, hexp_frw, tau_frw, t_frw
  real(real64), dimension(:,:), allocatable :: x, xout
  real(real64), dimension(:)  , allocatable :: m, age
  character(len=5) :: nchar, ncharcpu
  character(len=80) :: ordering
  character(len=80) :: gmgm
  character(len=128) :: filename, repository, outfich
  logical :: ok, periodic = .false., star = .false., ageweight = .false.
  logical :: ok_part
  integer :: impi, ndom, bit_length, maxdom
  integer, dimension(1:8) :: idom, jdom, kdom, cpu_min, cpu_max
  real(real64), dimension(1:8) :: bounding_min, bounding_max
  real(real64) :: dkey, order_min, dmax
  real(real64), dimension(:), allocatable :: bound_key
  logical, dimension(:), allocatable :: cpu_read
  integer, dimension(:), allocatable :: cpu_list
  integer(int8), dimension(:), allocatable :: family, tag

  integer :: unit_out, unit_in
  call read_params

  !-----------------------------------------------
  ! Reading particle files in RAMSES format
  !-----------------------------------------------
  ipos = INDEX(repository,'output_')
  nchar = repository(ipos+7:ipos+13)
  filename = TRIM(repository)//'/part_'//TRIM(nchar)//'.out00001'
  inquire(file = filename, exist = ok) ! verify input file
  if ( .not. ok ) then
     print *, TRIM(filename)//' not found.'
     stop
  end if

  filename = TRIM(repository)//'/info_'//TRIM(nchar)//'.txt'
  inquire(file = filename, exist = ok) ! verify input file
  if ( .not. ok ) then
     print *, TRIM(filename)//' not found.'
     stop
  end if
  open(newunit=unit_in, file=filename, form='formatted', status='old')
  ! read(unit_in,'("ncpu         = ", I11)') ncpu
  read(unit_in,'(A13, I11)') GMGM, ncpu
  ! read(unit_in,'("ndim         = ", I11)') ndim
  read(unit_in,'(A13, I11)') GMGM, ndim
  ! read(unit_in,'("levelmin     = ", I11)') levelmin
  read(unit_in,'(A13, I11)') GMGM, levelmin
  ! read(unit_in,'("levelmax     = ", I11)') levelmax
  read(unit_in,'(A13, I11)') GMGM, levelmax
  read(unit_in,*)
  read(unit_in,*)
  read(unit_in,*)

  read(unit_in,*)
  ! read(unit_in,'("time         = ", E23.15)') t
  read(unit_in,'(A13, E23.15)') GMGM, t

  !  read(unit_in,'("aexp         = ", E23.15)') aexp
  read(unit_in,'(A13, E23.15)') GMGM, aexp

  !  read(unit_in,'("H0           = ", E23.15)') h0
  read(unit_in,'(A13, E23.15)') GMGM, h0

  ! read(unit_in,'("omega_m      = ", E23.15)') omega_m
  read(unit_in,'(A13, E23.15)') GMGM, omega_m

  ! read(unit_in,'("omega_l      = ", E23.15)') omega_l
  read(unit_in,'(A13, E23.15)') GMGM, omega_l

  ! read(unit_in,'("omega_k      = ", E23.15)') omega_k
  read(unit_in,'(A13, E23.15)') GMGM, omega_k

  ! read(unit_in,'("omega_b      = ", E23.15)') omega_b
  read(unit_in,'(A13, E23.15)') GMGM, omega_b

  ! read(unit_in,'("unit_l       = ", E23.15)') unit_l
  read(unit_in,'(A13, E23.15)') GMGM, unit_l

  ! read(unit_in,'("unit_d       = ", E23.15)') unit_d
  read(unit_in,'(A13, E23.15)') GMGM, unit_d

  ! read(unit_in,'("unit_t       = ", E23.15)') unit_t
  read(unit_in,'(A13, E23.15)') GMGM, unit_t

  read(unit_in,*)

  ! read(unit_in,'("ordering type = ", A80)'), ordering
  read(unit_in,'(A14, A80)') GMGM, ordering
  write(*,'(" ordering type=", A20)') TRIM(ordering)
  read(unit_in,*)
  allocate(cpu_list(1:ncpu))
  if (TRIM(ordering) .eq.'hilbert') then
     allocate(bound_key(0:ncpu))
     allocate(cpu_read(1:ncpu))
     cpu_read = .false.
     do impi = 1, ncpu
        read(unit_in,'(I8, 1X, E23.15, 1X, E23.15)') i, bound_key(impi-1), bound_key(impi)
     end do
  end if
  close(unit_in)


  !-----------------------
  ! Cosmological model
  !-----------------------
  ! Allocate look-up tables
  n_frw = 1000
  allocate(aexp_frw(0:n_frw), hexp_frw(0:n_frw))
  allocate(tau_frw(0:n_frw), t_frw(0:n_frw))

  ! Compute Friedman model look up table
  write(*,*)'Computing Friedman model'
  call friedman(dble(omega_m), dble(omega_l), dble(omega_k), &
       & 1.d-6, 1.d-3, aexp_frw, hexp_frw, tau_frw, t_frw, n_frw, time_tot)

  ! Find neighboring expansion factors
  i = 1
  do while (aexp_frw(i) > aexp .and. i < n_frw)
     i = i+1
  end do
  ! Interploate time
  time_simu = t_frw(i)*(aexp-aexp_frw(i-1))/(aexp_frw(i)-aexp_frw(i-1))+ &
       & t_frw(i-1)*(aexp-aexp_frw(i))/(aexp_frw(i-1)-aexp_frw(i))
  write(*,*)'Age simu=',(time_tot+time_simu)/(h0*1d5/3.08d24)/(365.*24.*3600.*1d9)

  !-----------------------
  ! Map parameters
  !-----------------------
  if (nx == 0) then
     nx = 2**levelmin
  end if
  if (ny == 0) then
     ny = nx
  end if
  if (nz == 0) then
     nz = nx
  end if
  write(*,*)'time=', t

  if (TRIM(ordering) .eq.'hilbert') then

     dmax = max(xmax-xmin, ymax-ymin, zmax-zmin)
     do ilevel = 1, levelmax
        deltax = 0.5d0**ilevel
        if (deltax .lt. dmax) exit
     end do
     lmin = ilevel
     bit_length = lmin-1
     maxdom = 2**bit_length
     imin = 0; imax = 0; jmin = 0; jmax = 0; kmin = 0; kmax = 0
     if (bit_length > 0) then
        imin = int(xmin*dble(maxdom))
        imax = imin+1
        jmin = int(ymin*dble(maxdom))
        jmax = jmin+1
        kmin = int(zmin*dble(maxdom))
        kmax = kmin+1
     end if

     dkey = (dble(2**(levelmax+1)/dble(maxdom)))**ndim
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
        bounding_min(i) = (order_min)*dkey
        bounding_max(i) = (order_min+1.0D0)*dkey
     end do
     cpu_min = 0; cpu_max = 0
     do impi = 1, ncpu
        do i = 1, ndom
           if (   bound_key(impi-1) .le. bounding_min(i) .and.&
                & bound_key(impi  ) .gt. bounding_min(i)) then
              cpu_min(i) = impi
           end if
           if (   bound_key(impi-1) .lt. bounding_max(i) .and.&
                & bound_key(impi  ) .ge. bounding_max(i)) then
              cpu_max(i) = impi
           end if
        end do
     end do

     ncpu_read = 0
     do i = 1, ndom
        do j = cpu_min(i), cpu_max(i)
           if (.not. cpu_read(j)) then
              ncpu_read = ncpu_read+1
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
     icpu = cpu_list(k)
     call title(icpu, ncharcpu)
     filename = TRIM(repository)//'/part_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
     open(unit=1, file=filename, status='old', form='unformatted')
     read(1) ncpu2
     read(1) ndim2
     read(1) npart2
     read(1)
     read(1) nstar
     close(1)
     npart = npart+npart2
  end do
  write(*,*)'Found ', npart,' particles'
  if (nstar > 0) then
     if (star) then
        write(*,*)'Keeping star particles.'
     else
        write(*,*)'Discard star particles.'
     end if
  end if

  allocate(xout(1:npart, 1:ndim))
  j = 0
  !-----------------------------------------------
  ! Compute projected mass using CIC smoothing
  !----------------------------------------------
  mtot = 0.0d0
  do k = 1, ncpu_read
     icpu = cpu_list(k)
     call title(icpu, ncharcpu)
     filename = TRIM(repository)//'/part_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
     open(unit=1, file=filename, status='old', form='unformatted')
     write(*,*)'Processing file '//TRIM(filename)
     read(1) ncpu2
     read(1) ndim2
     read(1) npart2
     read(1)
     read(1)
     read(1)
     read(1)
     read(1)
     allocate(m(1:npart2))
     if (nstar > 0) allocate(age(1:npart2))
     allocate(x(1:npart2, 1:ndim2))
     allocate(family(1:npart2))
     allocate(tag(1:npart2))
     ! Read position
     do i = 1, ndim
        read(1) m
        x(1:npart2, i) = m
     end do
     ! Skip velocity
     do i = 1, ndim
        read(1) ! m
     end do
     ! Read mass
     read(1) m
     if (nstar > 0) then
        read(1) ! Skip identity
        read(1) ! Skip level
        read(1) family
        read(1) tag
        read(1) age
     end if

     close(1)
     if (periodic) then
        do i = 1, npart2
           ok_part = (x(i, 1) >= xmin .and. x(i, 1) <= xmax .and. &
                &   x(i, 2) >= ymin .and. x(i, 2) <= ymax .and. &
                &   x(i, 3) >= zmin .and. x(i, 3) <= zmax)

           if (nstar > 0) then
              if (star) then
                 ! Keep only stars
                 ok_part = ok_part .and. (family(i) == 2)
                 if (ageweight) then
                    iii = 1
                    do while (tau_frw(iii) > age(i) .and. iii < n_frw)
                       iii = iii+1
                    end do
                    ! Interpolate time
                    time = t_frw(iii)*(age(i)-tau_frw(iii-1))/(tau_frw(iii)-tau_frw(iii-1))+ &
                         & t_frw(iii-1)*(age(i)-tau_frw(iii))/(tau_frw(iii-1)-tau_frw(iii))
                    time = (time_simu-time)/(h0*1d5/3.08d24)/(365.*24.*3600.*1d9)
                    if (time > 0.01) then
                       weight = (time/0.01)**(-0.7)
                    end if
                 end if
              else
                 ! Keep only DM
                 ok_part = ok_part .and. (family(i) == 1)
              end if
           end if

           if (ok_part) then
              mtot = mtot + m(i)
              j = j + 1
              xout(j, :) = x(i, :)
           end if

        end do
     else
        do i = 1, npart2
           weight = 1.0
           ok_part = (x(i, 1) >= xmin .and. x(i, 1) <= xmax .and. &
                &   x(i, 2) >= ymin .and. x(i, 2) <= ymax .and. &
                &   x(i, 3) >= zmin .and. x(i, 3) <= zmax)

           if (nstar > 0) then
              if (star) then
                 ok_part = ok_part .and. (family(i) == 2)
                 if (ageweight) then
                    iii = 1
                    do while (tau_frw(iii) > age(i) .and. iii < n_frw)
                       iii = iii+1
                    end do
                    ! Interpolate time
                    time = t_frw(iii)*(age(i)-tau_frw(iii-1))/(tau_frw(iii)-tau_frw(iii-1))+ &
                         & t_frw(iii-1)*(age(i)-tau_frw(iii))/(tau_frw(iii-1)-tau_frw(iii))
                    time = (time_simu-time)/(h0*1d5/3.08d24)/(365.*24.*3600.*1d9)
                    if (time > 0.01) then
                       weight = (time/0.01)**(-0.7)
                    end if
                 end if
              else
                 ok_part = ok_part .and. (family(i) == 1)
              end if
           end if

           if (ok_part) then
              mtot = mtot + m(i)
              j = j + 1
              xout(j, :) = x(i, :)
           end if
        end do
     end if
     deallocate(x, m, family, tag)
     if (nstar > 0) deallocate(age)
  end do
  npart = j
  write(*,"('Total mass=', es14.5, ' in ', i10, ' particles')") mtot, npart

  ! Output file
  filename = TRIM(outfich)
  write(*,*)'Writing output to '//TRIM(filename)
  open(newunit=unit_out, file=filename, form='unformatted')
  write(unit_out) npart

  write(unit_out) xout(:npart, 1)
  write(unit_out) xout(:npart, 2)
  write(unit_out) xout(:npart, 3)
  close(unit_out)

contains

  subroutine read_params

    implicit none

    integer       :: i, n
    character(len=4)   :: opt
    character(len=128) :: arg

    n = command_argument_count()
    if (n < 4) then
       print *, 'usage: part2list  -inp  input_dir'
       print *, '                  -out  output_file'
       print *, '                 [-xmi xmin] '
       print *, '                 [-xma xmax] '
       print *, '                 [-ymi ymin] '
       print *, '                 [-yma ymax] '
       print *, '                 [-zmi zmin] '
       print *, '                 [-zma zmax] '
       print *, '                 [-nx  nx  ] '
       print *, '                 [-ny  ny  ] '
       print *, '                 [-nz  nz  ] '
       print *, '                 [-per flag] '
       print *, '                 [-str flag] '
       print *, 'ex: part2list -inp output_00001 -out list.dat'// &
            &   ' -xmi 0.1 -xma 0.7'
       stop
    end if

    do i = 1, n, 2
       call get_command_argument(i, opt)
       if (i == n) then
          print '("option ", a2," has no argument")', opt
          stop 2
       end if
       call get_command_argument(i+1, arg)
       select case (opt)
       case ('-age')
          read (arg,*) ageweight
       case ('-inp')
          repository = trim(arg)
       case ('-out')
          outfich = trim(arg)
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
       case ('-nx')
          read (arg,*) nx
       case ('-ny')
          read (arg,*) ny
       case ('-nz')
          read (arg,*) nz
       case ('-per')
          read (arg,*) periodic
       case ('-str')
          read (arg,*) star
       case default
          print '("unknown option ", a2," ignored")', opt
       end select
    end do

    return

  end subroutine read_params

end program part2list
