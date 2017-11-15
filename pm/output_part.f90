module dump_utils
  use iso_fortran_env
  implicit none

  private

  interface generic_dump
     module procedure logicaldump
     module procedure int8dump
     module procedure int16dump
     module procedure int32dump
     module procedure int64dump
     module procedure real32dump
     module procedure real64dump
  end interface generic_dump
  public :: generic_dump
contains

  subroutine logicaldump(varname, ivar, data, unit_out, dump_info, unit_info)
    character(len=*), intent(in) :: varname
    logical, intent(in), dimension(:) :: data
    integer, intent(inout) :: ivar
    integer, intent(in) :: unit_out, unit_info
    logical, intent(in) :: dump_info

    character(len=1) :: kind

    write(unit_out) data
    kind = '?'
    if (dump_info) call dump_var_info(varname, ivar, kind, unit_info)
    ivar = ivar + 1
  end subroutine logicaldump

  subroutine int8dump(varname, ivar, data, unit_out, dump_info, unit_info)
    character(len=*), intent(in) :: varname
    integer(int8), intent(in), dimension(:) :: data
    integer, intent(inout) :: ivar
    integer, intent(in) :: unit_out, unit_info
    logical, intent(in) :: dump_info

    character(len=1) :: kind

    write(unit_out) data
    kind = 'b'
    if (dump_info) call dump_var_info(varname, ivar, kind, unit_info)
    ivar = ivar + 1
  end subroutine int8dump

  subroutine int16dump(varname, ivar, data, unit_out, dump_info, unit_info)
    character(len=*), intent(in) :: varname
    integer(int16), intent(in), dimension(:) :: data
    integer, intent(inout) :: ivar
    integer, intent(in) :: unit_out, unit_info
    logical, intent(in) :: dump_info

    character(len=1) :: kind

    write(unit_out) data
    kind = 'h'
    if (dump_info) call dump_var_info(varname, ivar, kind, unit_info)
    ivar = ivar + 1
  end subroutine int16dump

  subroutine int32dump(varname, ivar, data, unit_out, dump_info, unit_info)
    character(len=*), intent(in) :: varname
    integer(int32), intent(in), dimension(:) :: data
    integer, intent(inout) :: ivar
    integer, intent(in) :: unit_out, unit_info
    logical, intent(in) :: dump_info

    character(len=1) :: kind

    write(unit_out) data
    kind = 'i'
    if (dump_info) call dump_var_info(varname, ivar, kind, unit_info)
    ivar = ivar + 1
  end subroutine int32dump

  subroutine int64dump(varname, ivar, data, unit_out, dump_info, unit_info)
    character(len=*), intent(in) :: varname
    integer(int64), intent(in), dimension(:) :: data
    integer, intent(inout) :: ivar
    integer, intent(in) :: unit_out, unit_info
    logical, intent(in) :: dump_info

    character(len=1) :: kind

    write(unit_out) data
    kind = 'q'
    if (dump_info) call dump_var_info(varname, ivar, kind, unit_info)
    ivar = ivar + 1
  end subroutine int64dump

  subroutine real32dump(varname, ivar, data, unit_out, dump_info, unit_info)
    character(len=*), intent(in) :: varname
    real(real32), intent(in), dimension(:) :: data
    integer, intent(inout) :: ivar
    integer, intent(in) :: unit_out, unit_info
    logical, intent(in) :: dump_info

    character(len=1) :: kind

    write(unit_out) data
    kind = 'f'
    if (dump_info) call dump_var_info(varname, ivar, kind, unit_info)
    ivar = ivar + 1
  end subroutine real32dump

  subroutine real64dump(varname, ivar, data, unit_out, dump_info, unit_info)
    character(len=*), intent(in) :: varname
    real(real64), intent(in), dimension(:) :: data
    integer, intent(inout) :: ivar
    integer, intent(in) :: unit_out, unit_info
    logical, intent(in) :: dump_info

    character(len=1) :: kind

    write(unit_out) data
    kind = 'd'
    if (dump_info) call dump_var_info(varname, ivar, kind, unit_info)
    ivar = ivar + 1
  end subroutine real64dump

  subroutine dump_var_info(varname, ivar, kind, unit_info)
    character(len=*), intent(in) :: varname, kind
    integer, intent(in) :: unit_info
    integer, intent(in) :: ivar

    write(unit_info, '(I2,", ", a, ", ", a)') ivar, trim(varname), trim(kind)
  end subroutine dump_var_info

end module dump_utils


subroutine backup_part(filename, file_desc)
  use amr_commons
  use pm_commons
  use dump_utils, only : generic_dump
  use iso_fortran_env
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer :: dummy_io, info2
  integer, parameter :: tag = 1122
#endif
  character(len=80) :: filename, file_desc

  integer :: i, idim, unit_out, ipart
  character(len=80) :: fileloc
  character(len=5) :: nchar
  real(dp), allocatable, dimension(:) :: xdp
  integer(i8b), allocatable, dimension(:) :: ii8
  integer, allocatable, dimension(:) :: ll
  integer(int8), allocatable, dimension(:) :: ii1

  integer :: unit_info, ivar
  logical :: dump_info

  character(len=1), dimension(1:3) :: dim_keys = (/"x", "y", "z"/)

  if (verbose) write(*,*) 'Entering backup_part'

  ! Set ivar to 1 for first variable
  ivar = 1

  ! Wait for the token
#ifndef WITHOUTMPI
  if (IOGROUPSIZE > 0) then
     if (mod(myid-1, IOGROUPSIZE) /= 0) then
        call MPI_RECV(dummy_io, 1, MPI_INTEGER, myid-1-1, tag, &
             & MPI_COMM_WORLD, MPI_STATUS_IGNORE, info2)
     end if
  end if
#endif


  call title(myid, nchar)
  fileloc = TRIM(filename) // TRIM(nchar)
  open(newunit=unit_out, file=TRIM(fileloc), form='unformatted')
  dump_info = .false.
  if (myid == 1) then
     open(newunit=unit_info, file=trim(file_desc), form='formatted')
     write(unit_info, '("# version: ", i2)') 1
     write(unit_info, '("# ivar, variable_name, variable_type")')
     dump_info = .true.
  end if

  rewind(unit_out)
  ! Write header
  write(unit_out) ncpu
  write(unit_out) ndim
  write(unit_out) npart
  write(unit_out) localseed
  write(unit_out) nstar_tot
  write(unit_out) mstar_tot
  write(unit_out) mstar_lost
  write(unit_out) nsink
  ! Write position
  allocate(xdp(1:npart))
  do idim = 1, ndim
     ipart = 0
     do i = 1, npartmax
        if (levelp(i) > 0) then
           ipart = ipart+1
           xdp(ipart) = xp(i, idim)
        end if
     end do
     call generic_dump("position_"//dim_keys(idim), ivar, xdp, unit_out, dump_info, unit_info)
  end do
  ! Write velocity
  do  idim = 1, ndim
     ipart = 0
     do i = 1, npartmax
        if (levelp(i) > 0) then
           ipart = ipart+1
           xdp(ipart) = vp(i, idim)
        end if
     end do
     call generic_dump("velocity_"//dim_keys(idim), ivar, xdp, unit_out, dump_info, unit_info)
  end do
  ! Write mass
  ipart = 0
  do i = 1, npartmax
     if (levelp(i) > 0) then
        ipart = ipart+1
        xdp(ipart) = mp(i)
     end if
  end do
  call generic_dump("mass", ivar, xdp, unit_out, dump_info, unit_info)
  deallocate(xdp)
  ! Write identity
  allocate(ii8(1:npart))
  ipart = 0
  do i = 1, npartmax
     if (levelp(i) > 0) then
        ipart = ipart+1
        ii8(ipart) = idp(i)
     end if
  end do
  call generic_dump("identity", ivar, ii8, unit_out, dump_info, unit_info)
  deallocate(ii8)

  ! Write level
  allocate(ll(1:npart))
  ipart = 0
  do i = 1, npartmax
     if (levelp(i) > 0) then
        ipart = ipart+1
        ll(ipart) = levelp(i)
     end if
  end do
  call generic_dump("levelp", ivar, ll, unit_out, dump_info, unit_info)

  deallocate(ll)

  ! Write family
  allocate(ii1(1:npart))
  ipart = 0
  do i = 1, npartmax
     if (levelp(i) > 0) then
        ipart = ipart+1
        ii1(ipart) = int(typep(i)%family, 1)
     end if
  end do
  call generic_dump("family", ivar, ii1, unit_out, dump_info, unit_info)

  ! Write tag
  ipart = 0
  do i = 1, npartmax
     if (levelp(i) > 0) then
        ipart = ipart+1
        ii1(ipart) = int(typep(i)%tag, 1)
     end if
  end do
  call generic_dump("tag", ivar, ii1, unit_out, dump_info, unit_info)
  deallocate(ii1)

#ifdef OUTPUT_PARTICLE_POTENTIAL
  ! Write potential (added by AP)
  allocate(xdp(1:npart))
  ipart = 0
  do i = 1, npartmax
     if (levelp(i) > 0) then
        ipart = ipart+1
        xdp(ipart) = ptcl_phi(i)
     end if
  end do
  call generic_dump("potential", ivar, xdp, unit_out, dump_info, unit_info)

  deallocate(xdp)
#endif

  ! Write birth epoch
  if (star .or. sink) then
     allocate(xdp(1:npart))
     ipart = 0
     do i = 1, npartmax
        if (levelp(i) > 0) then
           ipart = ipart+1
           xdp(ipart) = tp(i)
        end if
     end do
     call generic_dump("birth_time", ivar, xdp, unit_out, dump_info, unit_info)
     ! Write metallicity
     if (metal) then
        ipart = 0
        do i = 1, npartmax
           if (levelp(i) > 0) then
              ipart = ipart+1
              xdp(ipart) = zp(i)
           end if
        end do
        call generic_dump("metallicity", ivar, xdp, unit_out, dump_info, unit_info)
     end if
     deallocate(xdp)
  end if

  close(unit_out)
  if (myid == 1) close(unit_info)

  ! Send the token
#ifndef WITHOUTMPI
  if (IOGROUPSIZE > 0) then
     if (mod(myid, IOGROUPSIZE) /= 0 .and. (myid .lt. ncpu)) then
        dummy_io = 1
        call MPI_SEND(dummy_io, 1, MPI_INTEGER, myid-1+1, tag, &
             & MPI_COMM_WORLD, info2)
     end if
  end if
#endif

contains

end subroutine backup_part
