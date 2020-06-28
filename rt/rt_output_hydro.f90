! Outputting of radiation variables and groups info
!

!************************************************************************
subroutine rt_backup_hydro(filename, filename_desc)

!------------------------------------------------------------------------
  use amr_commons
  use rt_hydro_commons
  use rt_parameters
  use dump_utils, only : dump_header_info, generic_dump, dim_keys
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer :: dummy_io, info2
  integer, parameter :: tag = 1131
#endif
  character(len=*), intent(in) :: filename, filename_desc

  character(len=80) :: filedir, rt_filename

  integer :: i, ivar, idim, ncache, ind, ilevel, igrid, iskip, istart, ibound
  integer :: unit_out, unit_info
  integer, allocatable, dimension(:) :: ind_grid
  real(dp), allocatable, dimension(:) :: xdp
  character(len=5) :: nchar, ncharcpu
  character(len=80) :: fileloc

  logical :: dump_info_flag
  character(len=100) :: field_name
  integer :: info_var_count
!------------------------------------------------------------------------
  if (verbose) write(*,*)'Entering backup_rt'

  if (myid == 1) then
     call title(ifout-1, nchar)
     if (IOGROUPSIZEREP > 0) then
        call title(((myid-1)/IOGROUPSIZEREP)+1, ncharcpu)
        filedir = 'output_'//TRIM(nchar)//'/group_'//TRIM(ncharcpu)//'/'
     else
        filedir = 'output_'//TRIM(nchar)//'/'
     end if

     rt_filename = TRIM(filedir)//'info_rt_'//TRIM(nchar)//'.txt'
     call output_rtInfo(rt_filename)

     ! File descriptor
     open(newunit=unit_info, file=trim(filename_desc), form='formatted')
     call dump_header_info(unit_info)
     dump_info_flag = .true.
     info_var_count = 1
  else
     dump_info_flag = .false.
  end if

  if (.not. rt) return
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
  write(unit_out) nrtvar
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
              do ivar = 1, nGroups
                 ! Store photon density in flux units
                 do i = 1, ncache
                    xdp(i) = rt_c*rtuold(ind_grid(i)+iskip, iGroups(ivar))
                 end do
                 write(field_name, '("photon_density_", i0.2)') ivar
                 call generic_dump(field_name, info_var_count, xdp, unit_out, dump_info_flag, unit_info)
                 do idim = 1, ndim
                    ! Store photon flux
                    do i = 1, ncache
                       xdp(i) = rtuold(ind_grid(i)+iskip, iGroups(ivar)+idim)
                    end do
                    write(field_name, '("photon_flux_", i0.2, "_", a)') ivar, dim_keys(idim)
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

end subroutine rt_backup_hydro

!************************************************************************
SUBROUTINE output_rtInfo(filename)

! Output rt information into info_rt_XXXXX.txt
!------------------------------------------------------------------------
  use amr_commons
  use rt_hydro_commons
  use rt_parameters
  use rt_cooling_module
  implicit none
  character(LEN = 80) :: filename
  integer :: ilun
  real(dp) :: scale_np, scale_pf
  character(LEN = 80) :: fileloc
!------------------------------------------------------------------------
  if (verbose) write(*,*)'Entering output_rtInfo'

  ilun = myid+103

  ! Conversion factor from user units to cgs units
  call rt_units(scale_np, scale_pf)

  ! Open file
  fileloc = TRIM(filename)
  open(unit=ilun, file=fileloc, form='formatted')

  ! Write run parameters
  write(ilun,'("nRTvar       = ", I11)') nRTvar
  write(ilun,'("nIons        = ", I11)') nIons
  write(ilun,'("nGroups      = ", I11)') nGroups
  write(ilun,'("iIons        = ", I11)') iIons
  write(ilun,*)

  ! Write cooling parameters
  write(ilun,'("X_fraction   = ", E23.15)') X
  write(ilun,'("Y_fraction   = ", E23.15)') Y
  write(ilun,*)

  ! Write physical parameters
  write(ilun,'("unit_np      = ", E23.15)') scale_np
  write(ilun,'("unit_pf      = ", E23.15)') scale_pf
  write(ilun,'("rt_c_frac    = ", E23.15)') rt_c_fraction
  write(ilun,*)

  ! Write polytropic parameters
  write(ilun,'("n_star       = ", E23.15)') n_star
  write(ilun,'("T2_star      = ", E23.15)') T2_star
  write(ilun,'("g_star       = ", E23.15)') g_star
  write(ilun,*)
  call write_group_props(.false., ilun)

  close(ilun)

end subroutine output_rtInfo

!************************************************************************
SUBROUTINE write_group_props(update, lun)

! Write photon group properties to file or std output.
! lun => File identifier (use 6 for std. output)
!------------------------------------------------------------------------
  use rt_parameters
  use amr_commons, only:myid
  implicit none
  logical :: update
  integer :: ip, lun
!------------------------------------------------------------------------
  if (myid .ne. 1) RETURN
  if (.not. update) then
     write(lun,*) 'Photon group properties------------------------------ '
  else
     write(lun,*) 'Photon properties have been changed to----------------- '
  end if
  write(lun, 901) groupL0(:)
  write(lun, 902) groupL1(:)
  write(lun, 903) spec2group(:)
  do ip = 1, nGroups
     write(lun, 907) ip
     write(lun, 904) group_egy(ip)
     write(lun, 905) group_csn(ip,:)
     write(lun, 906) group_cse(ip,:)
  end do
  write (lun,*) '-------------------------------------------------------'

901 format ('  groupL0  [eV]  = ', 20f12.3)
902 format ('  groupL1  [eV]  = ', 20f12.3)
903 format ('  spec2group     = ', 20I12)
904 format ('  egy      [eV]  = ', 1pe12.3)
905 format ('  csn    [cm^2]  = ', 20(1pe12.3))
906 format ('  cse    [cm^2]  = ', 20(1pe12.3))
907 format ('  ---Group', I2)

END SUBROUTINE write_group_props

!*************************************************************************
SUBROUTINE output_rt_stats

! Output and reset rt statistics. These are cooling statistics and
! star rt feedback statistics
!-------------------------------------------------------------------------
  use amr_commons
  use rt_parameters
  use constants, only: yr2sec
  use mpi_mod
  implicit none
  integer(8) :: max_all, tot_all, cells_all, loopCodes_tot
  integer(8) :: loopCodes_all(20)
  real(dp) :: step_nPhot_all, step_nStar_all, step_mStar_all
  real(dp) :: scale_l, scale_t, scale_d, scale_v, scale_nh, scale_T2
#ifndef WITHOUTMPI
  integer :: info
#endif
!-------------------------------------------------------------------------
  call units(scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2)
  ! Cooling statistics:
  if (neq_chem .and. rt_output_coolstats) then
     cells_all = 0 ; tot_all = 0 ; max_all = 0 ; loopCodes_all = 0
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(n_cool_cells,         cells_all,     1, &
          MPI_INTEGER,          MPI_SUM, MPI_COMM_WORLD, info)
     call MPI_ALLREDUCE(tot_cool_loopcnt,     tot_all,       1, &
          MPI_INTEGER,          MPI_SUM, MPI_COMM_WORLD, info)
     call MPI_ALLREDUCE(max_cool_loopcnt,     max_all,       1, &
          MPI_INTEGER,          MPI_MAX, MPI_COMM_WORLD, info)
     call MPI_ALLREDUCE(loopCodes,            loopCodes_all, 20, &
          MPI_INTEGER,          MPI_MAX, MPI_COMM_WORLD, info)
     n_cool_cells     = cells_all ; tot_cool_loopcnt = tot_all
     max_cool_loopcnt = max_all   ; loopCodes        = loopCodes_all
#endif
     if (myid .eq. 1) then
        if (n_cool_cells .eq. 0) n_cool_cells = 1
        write(*, 111) dble(tot_cool_loopcnt)/n_cool_cells, max_cool_loopcnt, rt_advect
        loopCodes_tot = SUM(loopCodes)
        if (loopCodes_tot .gt. 0) then
           write(*, 112) dble(loopCodes(1:7))/dble(loopCodes_tot)
        else
           write(*, 112) dble(loopCodes(1:7))
        end if
     end if
     max_cool_loopcnt = 0; tot_cool_loopcnt = 0; n_cool_cells = 0; loopCodes(:) = 0
  end if ! output_coolstats
111 format(' Coolstats: Avg. # loops = ', f21.6, ', max. # loops = ', I10, ', rt_adv = ', L)
112 format(' Subcycling codes [Np, Fp, T, Npd, Td, xHII, xe, xH2]% = ', 8(f7.3, ''))

  ! Stellar rt feedback statistics:
  if (showSEDstats .and. rt_star) then
     step_nPhot_all = 0d0 ; step_nStar_all = 0d0 ; step_mStar_all = 0d0
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(step_nPhot,           step_nPhot_all,  1,        &
          MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
     call MPI_ALLREDUCE(step_nStar,           step_nStar_all,  1,        &
          MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
     call MPI_ALLREDUCE(step_mStar,           step_mStar_all,  1,        &
          MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
     step_nPhot  = step_nPhot_all
     step_nStar  = step_nStar_all
     step_mStar  = step_mStar_all
#endif
     tot_nPhot = tot_nPhot + step_nPhot
     if (myid .eq. 1)                                                     &
          write(*, 113) step_nPhot, tot_nPhot, step_nStar/dtnew(levelmin)&
          , step_mStar/dtnew(levelmin), dtnew(levelmin)*scale_t/yr2sec
     step_nPhot = 0d0 ; step_nStar = 0d0 ; step_mStar = 0d0
  end if
113 format(' SED feedback(phot/step/1d50, phot/tot/1d50, *, */Msun , dt[yr])= '  &
                                                             , 10(1pe9.2))
END SUBROUTINE output_rt_stats
