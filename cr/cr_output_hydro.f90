! Outputting of cosmic-ray variables: separate cr_NNNNN.out mode.
!************************************************************************
subroutine cr_backup_hydro(filename, filename_desc)
!------------------------------------------------------------------------
  use amr_commons
  use cr_parameters
  use cr_hydro_commons
  use hydro_parameters, only: gamma
  use dump_utils, only : dump_header_info, generic_dump, dim_keys
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer :: dummy_io, info2
  integer, parameter :: tag = 1132
#endif
  character(len=*), intent(in) :: filename, filename_desc

  integer :: i, idim, igroup, ncache, ind, ilevel, igrid, iskip, istart, ibound
  integer :: unit_out, unit_info
  integer, allocatable, dimension(:) :: ind_grid
  real(dp), allocatable, dimension(:) :: xdp
  character(len=5) :: nchar
  character(len=80) :: fileloc
  logical :: dump_info_flag
  character(len=100) :: field_name
  integer :: info_var_count
!------------------------------------------------------------------------
  if (verbose) write(*,*)'Entering cr_backup_hydro'

  if (myid == 1) then
     open(newunit=unit_info, file=trim(filename_desc), form='formatted')
     call dump_header_info(unit_info)
     dump_info_flag = .true.
     info_var_count = 1
  else
     dump_info_flag = .false.
  end if

  ! Wait for the token
#ifndef WITHOUTMPI
  if (IOGROUPSIZE > 0) then
     if (mod(myid-1, IOGROUPSIZE) /= 0) then
        call MPI_RECV(dummy_io, 1, MPI_INTEGER, myid-1-1, tag,&
             & MPI_COMM_WORLD, MPI_STATUS_IGNORE, info2)
     end if
  end if
#endif

  call title(myid, nchar)
  fileloc = TRIM(filename)//TRIM(nchar)
  open(newunit=unit_out, file=fileloc, form='unformatted')
  write(unit_out) ncpu
  write(unit_out) ncrvar
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
              do igroup = 1, ncr_groups
                 ! CR energy density of group igroup
                 do i = 1, ncache
                    xdp(i) = cruold(ind_grid(i)+iskip, iCRu+(ndim+1)*(igroup-1))
                 end do
                 write(field_name, '("CRegy_", i0.2)') igroup
                 call generic_dump(field_name, info_var_count, xdp, unit_out, dump_info_flag, unit_info)
                 do idim = 1, ndim
                    ! CR flux components of group igroup
                    do i = 1, ncache
                       xdp(i) = cruold(ind_grid(i)+iskip, iCRu+(ndim+1)*(igroup-1)+idim)
                    end do
                    write(field_name, '("CRflx_", i0.2, "_", a1)') igroup, dim_keys(idim)
                    call generic_dump(field_name, info_var_count, xdp, unit_out, dump_info_flag, unit_info)
                 end do
              end do
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

end subroutine cr_backup_hydro
