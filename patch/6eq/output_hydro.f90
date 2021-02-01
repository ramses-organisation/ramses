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
  logical :: inv
  integer :: info_var_count
  integer :: idim,imat
  character(len=100) :: field_name
  real(dp)::ekin,erad
  real(dp),dimension(1:nvector,1:nmat),save::ff,gg
  real(dp),dimension(1:nvector,1:npri),save::qq
  real(dp),dimension(1:nvector),save::dtot,pp,cc

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
  write(unit_out) nvar
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
              ! Calculate total density 
              dtot(1:ncache) = 0.0
              do imat = 1,nmat
                do i=1,ncache
                  dtot(i) = dtot(i) + uold(ind_grid(i)+iskip,nmat+imat)
                end do
              end do
              ! Write volume fractions
              do imat = 1,nmat
                 ivar = imat
                 do i = 1, ncache
                    xdp(i) = uold(ind_grid(i)+iskip,ivar)
                 end do
                 write(field_name, '("vol_frac_", i0.2)') imat
                 call generic_dump(field_name, info_var_count, xdp, unit_out, dump_info_flag, unit_info)
              end do
              ! Write true densities
              do imat = 1,nmat
                 ivar = nmat+imat
                 do i = 1, ncache
                    xdp(i) = uold(ind_grid(i)+iskip,ivar)/max(uold(ind_grid(i)+iskip,imat),smallf)
                 end do
                 write(field_name, '("true_dens_", i0.2)') imat
                 call generic_dump(field_name, info_var_count, xdp, unit_out, dump_info_flag, unit_info)
              end do
              do ivar = 2*nmat+1, 2*nmat+ndim
                do i = 1, ncache
                  xdp(i) = uold(ind_grid(i)+iskip, ivar)/max(dtot(i),smallr)
                end do
                field_name = 'velocity_' // dim_keys(ivar - 1)
                call generic_dump(field_name, info_var_count, xdp, unit_out, dump_info_flag, unit_info)
              end do
#if NENER > 0
              ! Write non-thermal pressures
              do ivar = 3*nmat+ndim+1,3*nmat+ndim+nener
                 do i = 1, ncache
                    xdp(i) = (gamma_rad(ivar-ndim-2)-1d0)*uold(ind_grid(i)+iskip, ivar)
                 end do
                 write(field_name, '("non_thermal_energy_", i0.2)') ivar-3
                 call generic_dump(field_name, info_var_count, xdp, unit_out, dump_info_flag, unit_info)
              end do
#endif
              ! Calculate individual internal + radiative energies
              inv=.false.
              do imat = 1,nmat
                do i = 1, ncache
                  ff(i,imat)   = uold(ind_grid(i)+iskip,imat)
                  gg(i,imat)   = uold(ind_grid(i)+iskip,imat+nmat)/max(ff(i,imat),smallf)
                  ekin=0.0
                  do idim=1,ndim
                    qq(i,idim) = uold(ind_grid(i)+iskip,2*nmat+idim)/max(dtot(i),smallr)
                    ekin       = ekin + 0.5d0*qq(i,idim)**2
                  end do
                  erad=00
#if NENER > 0
                  do irad = 1,nener
                    erad       = erad + uold(ind_grid(i)+iskip,3*nmat+ndim+irad)
                  end do
#endif
                  qq(i,ndim+nmat+imat) = uold(ind_grid(i)+iskip,2*nmat+ndim+imat)/max(ff(i,imat),smallf) - gg(i,imat)*ekin - erad
                end do
              end do

              ! Write thermal pressure
              do imat=1,nmat
                call eos(gg(:,imat),qq(:,ndim+nmat+imat),pp,cc,imat,inv,ncache)
                do i=1,ncache
                  xdp(i) = pp(i)       ! Pressure
                end do
                write(field_name, '("pressure_", i0.2)') imat
                call generic_dump(field_name, info_var_count, xdp, unit_out, dump_info_flag, unit_info)
              end do


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






