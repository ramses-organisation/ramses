subroutine read_nimhd_params(nml_ok)
  use amr_commons
  use nimhd_parameters
  implicit none
  logical, intent(inout) ::nml_ok
  !--------------------------------------------------
  ! Namelist definitions
  !--------------------------------------------------
  namelist/nonidealmhd_params/nambipolar,gammaAD &
        & ,nmagdiffu,etaMD,ntestDADM,use_x3d &
        & ,coefad, nminitimestep, coefalfven,coefdtohm

  ! Checks on non-ideal MHD parameters
  rewind(1)
  read(1,NML=nonidealmhd_params,END=109)

  if(nambipolar.or.nmagdiffu)then
    use_nonideal_mhd = .true.
  else
    use_nonideal_mhd = .false.
  endif

  if(.not.use_nonideal_mhd) return
  
  if(myid==1) then
    write(*,*)'!!!!!!!!!!!!!!!  Non Ideal MHD   !!!!!!!!!!!!!!!!'
    write(*,*)'Non ideal MHD parameters'
    write(*,*)'Making a test ? (Yes=1 No=0)',ntestDADM
    if(nambipolar) then
      write(*,*)'Ambipolar diffusion switched ON'
      write(*,*)'Ambipolar diffusion coefficient',gammaAD
      write(*,*)'Ambipolar diffusion time coefficient',coefad
      write(*,*)'Ionisation coefficient',coefionis
      if(nminitimestep.eq.1) then
        write(*,*)'Mini time step switched ON'
        write(*,*)'Mini time step coefficient',coefalfven
      else
        write(*,*)'Mini time step switched OFF'
      endif
    else
      write(*,*)'Ambipolar diffusion switched OFF'
    endif
  
    if(nmagdiffu)then
      write(*,*)'Magnetic diffusion switched ON : multiple time stepping'
      write(*,*)'Magnetic diffusion coefficient',etaMD
      write(*,*)'Magnetic diffusion  time coefficient',coefohm
    else
      write(*,*)'Magnetic diffusion switched OFF'
    endif

  endif

  call read_resistivities

  109 continue

end subroutine read_nimhd_params