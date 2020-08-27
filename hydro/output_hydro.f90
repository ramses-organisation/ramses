subroutine backup_hydro(filename, filename_desc)
  use amr_commons
  use hydro_commons
  use dump_utils, only : dump_header_info, generic_dump, dim_keys
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer :: dummy_io, info2
#endif

  character(len=80), intent(in) :: filename, filename_desc

  integer :: i, ivar, ncache, ind, ilevel, igrid, iskip, istart, ibound
  integer :: unit_out, unit_info
  integer, allocatable, dimension(:) :: ind_grid
  real(dp), allocatable, dimension(:) :: xdp
  character(LEN = 5) :: nchar
  character(LEN = 80) :: fileloc
  integer, parameter :: tag = 1121
#if NENER > 0
  integer :: irad
#endif
  logical :: dump_info_flag
  integer :: info_var_count
  character(len=100) :: field_name

  if (verbose) write(*,*)'Entering backup_hydro'

  call title(myid, nchar)
  fileloc = TRIM(filename)//TRIM(nchar)

  ! Wait for the token
#ifndef WITHOUTMPI
  if (IOGROUPSIZE > 0) then
     if (mod(myid-1, IOGROUPSIZE) /= 0) then
        call MPI_RECV(dummy_io, 1, MPI_INTEGER, myid-1-1, tag,&
             & MPI_COMM_WORLD, MPI_STATUS_IGNORE, info2)
     end if
  end if
#endif

  open(newunit=unit_out, file=fileloc, form='unformatted')

  if (myid == 1) then
     open(newunit=unit_info, file=filename_desc, form='formatted')
     call dump_header_info(unit_info)
     info_var_count = 1
     dump_info_flag = .true.
  else
     dump_info_flag = .false.
  end if

  write(unit_out) ncpu
  if(strict_equilibrium>0)then
     write(unit_out) nvar+2
  else
     write(unit_out) nvar
  endif
  write(unit_out) ndim
  write(unit_out) nlevelmax
  write(unit_out) nboundary
  write(unit_out) gamma
  do ilevel = 1, nlevelmax
     do ibound = 1, nboundary+ncpu
        if (ibound <= ncpu) then
           ncache = numbl(ibound, ilevel)
           istart = headl(ibound, ilevel)
        else
           ncache = numbb(ibound-ncpu, ilevel)
           istart = headb(ibound-ncpu, ilevel)
        end if
        write(unit_out) ilevel
        write(unit_out) ncache
        if (ncache > 0) then
           allocate(ind_grid(1:ncache), xdp(1:ncache))
           ! Loop over level grids
           igrid = istart
           do i = 1, ncache
              ind_grid(i) = igrid
              igrid = next(igrid)
           end do
           ! Loop over cells
           do ind = 1, twotondim
              iskip = ncoarse+(ind-1)*ngridmax
              do ivar = 1, ndim+1
                 if (ivar == 1) then
                    ! Write density
                    do i = 1, ncache
                       xdp(i) = uold(ind_grid(i)+iskip, 1)
                    end do
                    field_name = 'density'
                 else if (ivar >= 2 .and. ivar <= ndim+1) then
                    ! Write velocity field
                    do i = 1, ncache
                       xdp(i) = uold(ind_grid(i)+iskip, ivar)/max(uold(ind_grid(i)+iskip, 1), smallr)
                    end do
                    field_name = 'velocity_' // dim_keys(ivar - 1)
                 end if
                 call generic_dump(field_name, info_var_count, xdp, unit_out, dump_info_flag, unit_info)
              end do
#if NENER > 0
              ! Write non-thermal pressures
              do ivar = ndim+3, ndim+2+nener
                 do i = 1, ncache
                    xdp(i) = (gamma_rad(ivar-ndim-2)-1d0)*uold(ind_grid(i)+iskip, ivar)
                 end do
                 write(field_name, '("non_thermal_energy_", i0.2)') ivar-3
                 call generic_dump(field_name, info_var_count, xdp, unit_out, dump_info_flag, unit_info)
              end do
#endif
              ! Write thermal pressure
              do i = 1, ncache
                 xdp(i) = uold(ind_grid(i)+iskip, ndim+2)
                 xdp(i) = xdp(i)-0.5d0*uold(ind_grid(i)+iskip, 2)**2/max(uold(ind_grid(i)+iskip, 1), smallr)
#if NDIM > 1
                 xdp(i) = xdp(i)-0.5d0*uold(ind_grid(i)+iskip, 3)**2/max(uold(ind_grid(i)+iskip, 1), smallr)
#endif
#if NDIM > 2
                 xdp(i) = xdp(i)-0.5d0*uold(ind_grid(i)+iskip, 4)**2/max(uold(ind_grid(i)+iskip, 1), smallr)
#endif
#if NENER > 0
                 do irad = 1, nener
                    xdp(i) = xdp(i)-uold(ind_grid(i)+iskip, ndim+2+irad)
                 end do
#endif
                 xdp(i) = (gamma-1d0)*xdp(i)
              end do
              field_name = 'pressure'
              call generic_dump(field_name, info_var_count, xdp, unit_out, dump_info_flag, unit_info)
#if NVAR > NDIM+2+NENER
              ! Write passive scalars
              do ivar = ndim+3+nener, nvar
                 do i = 1, ncache
                    xdp(i) = uold(ind_grid(i)+iskip, ivar)/max(uold(ind_grid(i)+iskip, 1), smallr)
                 end do
                 if (metal .and. imetal == ivar) then
                    field_name = 'metallicity'
                 else
                    write(field_name, '("scalar_", i0.2)') ivar - ndim - 3 - nener
                 end if
                 call generic_dump(field_name, info_var_count, xdp, unit_out, dump_info_flag, unit_info)
              end do
#endif
              if(strict_equilibrium>0)then
                 do i = 1, ncache
                    xdp(i) = rho_eq(ind_grid(i)+iskip)
                 end do
                 field_name = 'equilibrium_density'
                 call generic_dump(field_name, info_var_count, xdp, unit_out, dump_info_flag, unit_info)
                 do i = 1, ncache
                    xdp(i) = p_eq(ind_grid(i)+iskip)
                 end do
                 field_name = 'equilibrium_pressure'
                 call generic_dump(field_name, info_var_count, xdp, unit_out, dump_info_flag, unit_info)
              endif
              ! We did one output, deactivate dumping of variables
              dump_info_flag = .false.
           end do
           deallocate(ind_grid, xdp)

        end if
     end do
  end do
  close(unit_out)

  if (myid == 1) close(unit_info)
  ! Send the token
#ifndef WITHOUTMPI
  if (IOGROUPSIZE > 0) then
     if (mod(myid, IOGROUPSIZE) /= 0 .and.(myid .lt. ncpu)) then
        dummy_io = 1
        call MPI_SEND(dummy_io, 1, MPI_INTEGER, myid-1+1, tag, &
             & MPI_COMM_WORLD, info2)
     end if
  end if
#endif


end subroutine backup_hydro
