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
  integer :: idim,imat
  character(len=100) :: field_name
  real(dp)::dtot,ekin,erad
  real(dp),dimension(1:nvector,1:nmat),save::ff,gg,kk_mat
  real(dp),dimension(1:nvector,1:npri),save::qq
  real(dp),dimension(1:nvector),save::pp,cc,kk_hat

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
              do ivar = 1, ndim+1
                 if (ivar == 1) then
                    ! Write total density
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
                 dtot=max(uold(ind_grid(i)+iskip,1),smallr)
                 do imat=1,nmat
                    ff(1,imat)=uold(ind_grid(i)+iskip,imat+nener+npri)
                    gg(1,imat)=uold(ind_grid(i)+iskip,imat+nener+npri+nmat)
                 end do
                 qq(1,1)=dtot
                 ekin=0.0
                 do idim=1,ndim
                    qq(1,idim+1)=uold(ind_grid(i)+iskip,idim+1)/dtot
                    ekin=ekin+0.5d0*qq(1,idim+1)**2
                 end do
                 erad=00
#if NENER > 0
                 do irad = 1,nener
                    erad=erad+uold(ind_grid(i)+iskip,ndim+2+irad)
                 end do
#endif
                 qq(1,npri)=uold(ind_grid(i)+iskip,npri)-dtot*ekin-erad
                 call eos(ff,gg,qq,pp,cc,kk_mat,kk_hat,1)
                 xdp(i)=pp(1)       ! Pressure
              end do
              field_name = 'pressure'
              call generic_dump(field_name, info_var_count, xdp, unit_out, dump_info_flag, unit_info)
              ! Write volume fractions
              do imat = 1,nmat
                 ivar = ndim+2+nener+imat
                 do i = 1, ncache
                    xdp(i) = uold(ind_grid(i)+iskip,ivar)
                 end do
                 write(field_name, '("vol_frac_", i0.2)') imat
                 call generic_dump(field_name, info_var_count, xdp, unit_out, dump_info_flag, unit_info)
              end do
              ! Write true density
              do imat = 1,nmat
                 ivar = ndim+2+nener+nmat+imat
                 do i = 1, ncache
                    xdp(i) = uold(ind_grid(i)+iskip,ivar)
                 end do
                 write(field_name, '("true_dens_", i0.2)') imat
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
