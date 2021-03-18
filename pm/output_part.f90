subroutine backup_part(filename, filename_desc)
  use amr_commons
  use pm_commons
  use dump_utils, only : generic_dump, dump_header_info, dim_keys
  use iso_fortran_env
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer :: dummy_io, info2
  integer, parameter :: tag = 1122
#endif
  character(len=80) :: filename, filename_desc

  integer :: i, idim, unit_out, ipart
  character(len=80) :: fileloc
  character(len=5) :: nchar
  real(dp), allocatable, dimension(:) :: xdp
  integer(i8b), allocatable, dimension(:) :: ii8
  integer, allocatable, dimension(:) :: ll
  integer(int8), allocatable, dimension(:) :: ii1

  integer :: unit_info, ivar
  logical :: dump_info

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
  if (myid == 1) then
     open(newunit=unit_info, file=trim(filename_desc), form='formatted')
     call dump_header_info(unit_info)
     dump_info = .true.
  else
     dump_info = .false.
  end if

  rewind(unit_out)
  ! Write header
  write(unit_out) ncpu
  write(unit_out) ndim
  write(unit_out) npart
  if (MC_tracer) then
     write(unit_out) localseed, tracer_seed
  else
     write(unit_out) localseed
  end if
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

  if (MC_tracer) then
     ! Dump particle pointer
     allocate(ll(1:npart))
     ! Get the idp of the stars on which tracers are attached
     ipart = 0
     do i = 1, npartmax
        if (levelp(i) > 0) then
           ipart = ipart + 1
           ! For star tracers, store the id of the star instead of local index
           if (is_star_tracer(typep(i))) then
              ll(ipart) = idp(partp(i))
           else ! store the relative location
              ll(ipart) = partp(i)
           end if
        end if
     end do

     call generic_dump("partp", ivar, ll, unit_out, dump_info, unit_info)
     deallocate(ll)
  end if

  !------------!
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

end subroutine backup_part
