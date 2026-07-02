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
     write(unit_out) nvar_all+2
  else
     write(unit_out) nvar_all
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
              ! Write density
              ! (always first because we need to to convert from/to primitive variables)
              field_name = 'density'
              call gather_conservative_from_uold(ind_grid, iskip, 1, xdp, ncache)
              call generic_dump(field_name, info_var_count, xdp, unit_out, dump_info_flag, unit_info)
              ! Write velocity field
              do ivar = 2, neul-1
                 if(write_conservative)then
                    field_name = 'momentum_' // dim_keys(ivar - 1)
                    call gather_conservative_from_uold(ind_grid, iskip, ivar, xdp, ncache)
                 else
                    field_name = 'velocity_' // dim_keys(ivar - 1)
                    call gather_primitive_from_uold(ind_grid, iskip, ivar, xdp, ncache)
                 end if
                 call generic_dump(field_name, info_var_count, xdp, unit_out, dump_info_flag, unit_info)
              end do
#ifdef SOLVERmhd
              ! Write left B field
              ! (before thermal pressure because we need it to convert between total energy and pressure)
              do ivar = 6, 8
                 field_name = 'B_' // dim_keys(ivar - 6 + 1) // '_left'
                 call gather_conservative_from_uold(ind_grid, iskip, ivar, xdp, ncache)
                 call generic_dump(field_name, info_var_count, xdp, unit_out, dump_info_flag, unit_info)
              end do
              ! Write right B field
              do ivar = nvar+1, nvar+3
                 field_name = 'B_' // dim_keys(ivar - (nvar+1) + 1) // '_right'
                 call gather_conservative_from_uold(ind_grid, iskip, ivar, xdp, ncache)
                 call generic_dump(field_name, info_var_count, xdp, unit_out, dump_info_flag, unit_info)
              end do
#endif
#if NENER > 0
              ! Write non-thermal pressures
              ! (before thermal pressure because we need it to convert between total energy and pressure)
              do ivar = nhydro+1, nhydro+nener
                 if(write_conservative)then
                    write(field_name, '("non_thermal_energy_", i0.2)') ivar-nhydro
                    call gather_conservative_from_uold(ind_grid, iskip, ivar, xdp, ncache)
                 else
                    write(field_name, '("non_thermal_pressure_", i0.2)') ivar-nhydro
                    call gather_primitive_from_uold(ind_grid, iskip, ivar, xdp, ncache)
                 end if
                 call generic_dump(field_name, info_var_count, xdp, unit_out, dump_info_flag, unit_info)
              end do
#endif
              if(write_conservative) then
                 ! Write total energy as stored in uold
                 field_name = 'total_energy'
                 call gather_conservative_from_uold(ind_grid, iskip, neul, xdp, ncache)
              else
                 ! Write thermal pressure (after all other pressures or energies)
                 field_name = 'pressure'
                 call calc_thermal_pressure_from_total_energy(ind_grid, iskip, xdp, ncache)
              end if
              call generic_dump(field_name, info_var_count, xdp, unit_out, dump_info_flag, unit_info)
#if NVAR > NHYDRO+NENER
              ! Write passive scalars if any
              do ivar = nhydro+1+nener, nvar
                 if(write_conservative) then
                    if (metal .and. imetal == ivar) then
                       field_name = 'metal_density'
                    else
                       write(field_name, '("scalar_density_", i0.2)') ivar - nhydro - 1 - nener
                    end if
                    call gather_conservative_from_uold(ind_grid, iskip, ivar, xdp, ncache)
                 else
                    if (metal .and. imetal == ivar) then
                       field_name = 'metallicity'
                    else
                       write(field_name, '("scalar_", i0.2)') ivar - nhydro - 1 - nener
                    end if
                    call gather_primitive_from_uold(ind_grid, iskip, ivar, xdp, ncache)
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
!#####################################################################
!#####################################################################
!#####################################################################
subroutine gather_conservative_from_uold(ind_grid, iskip, ivar, xdp, ncache)
   use amr_parameters, only:dp
   use hydro_commons
   implicit none
   integer,intent(in)::ivar,iskip,ncache
   integer, dimension(1:ncache),intent(in)::ind_grid
   real(dp), dimension(1:ncache),intent(out)::xdp
   !-----------------------------------------------------------
   ! Gather the variable that is present in uold at index ivar
   ! (e.g. density, B_left, B_right)
   !-----------------------------------------------------------
   integer::i

   do i = 1, ncache
      xdp(i) = uold(ind_grid(i)+iskip, ivar)
   end do

end subroutine gather_conservative_from_uold
!#####################################################################
!#####################################################################
!#####################################################################
subroutine gather_primitive_from_uold(ind_grid, iskip, ivar, xdp, ncache)
   use amr_parameters, only:dp
   use hydro_commons
   implicit none
   integer,intent(in)::ivar,iskip,ncache
   integer, dimension(1:ncache),intent(in)::ind_grid
   real(dp), dimension(1:ncache),intent(out)::xdp
   !-----------------------------------------------------------
   ! Gather the primitive variable from the conservative variable
   ! that is present in uold at index ivar
   ! (e.g. velocity)
   !-----------------------------------------------------------
   integer::i

   select case(ivar)

#if NENER > 0
   case(nhydro+1:nhydro+nener)
      ! non-thermal pressures
      do i = 1, ncache
         xdp(i) = (gamma_rad(ivar-nhydro)-1d0)*uold(ind_grid(i)+iskip, ivar)
      end do
#endif

   case default
      do i = 1, ncache
         xdp(i) = uold(ind_grid(i)+iskip, ivar)/max(uold(ind_grid(i)+iskip, 1), smallr)
      end do

   end select

end subroutine gather_primitive_from_uold
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine calc_thermal_pressure_from_total_energy(ind_grid, iskip, pressure, ncache)
   use amr_parameters, only:dp
   use hydro_commons
   implicit none
   integer,intent(in)::iskip,ncache
   integer, dimension(1:ncache),intent(in)::ind_grid
   real(dp), dimension(1:ncache),intent(out)::pressure
   !--------------------------------------------------------------------------------------
   ! Calculate the thermal pressure from the total energy,
   ! which is stored in uold(:,neul)
   !--------------------------------------------------------------------------------------
   integer::i
   real(dp)::d,energy
#if NENER > 0
   integer :: irad
#endif
#ifdef SOLVERmhd
   real(dp) :: A, B, C
#endif

   do i = 1, ncache
      d = max(uold(ind_grid(i)+iskip, 1), smallr)
      ! total energy
      energy = uold(ind_grid(i)+iskip, neul)
      ! subtract kinetic energy
      energy = energy - 0.5d0*uold(ind_grid(i)+iskip, 2)**2/d
#if NDIM > 1 || SOLVERmhd
      energy = energy - 0.5d0*uold(ind_grid(i)+iskip, 3)**2/d
#endif
#if NDIM > 2 || SOLVERmhd
      energy = energy - 0.5d0*uold(ind_grid(i)+iskip, 4)**2/d
#endif
#ifdef SOLVERmhd
      ! subtract magnetic energy
      A = 0.5d0*(uold(ind_grid(i)+iskip, 6)+uold(ind_grid(i)+iskip, nvar+1))
      B = 0.5d0*(uold(ind_grid(i)+iskip, 7)+uold(ind_grid(i)+iskip, nvar+2))
      C = 0.5d0*(uold(ind_grid(i)+iskip, 8)+uold(ind_grid(i)+iskip, nvar+3))
      energy = energy - 0.5*(A**2+B**2+C**2)
#endif
#if NENER > 0
      ! subtract non-thermal energies
      do irad = 1, nener
         energy = energy-uold(ind_grid(i)+iskip, nhydro+irad)
      end do
#endif

      ! convert to pressure
      pressure(i) = (gamma-1d0)*energy
   end do

end subroutine calc_thermal_pressure_from_total_energy
