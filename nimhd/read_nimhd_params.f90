subroutine read_nimhd_params(nml_ok)
  use amr_commons
  use nimhd_parameters
  implicit none
  logical, intent(inout) ::nml_ok
  !--------------------------------------------------
  ! Namelist definitions
  !--------------------------------------------------
  namelist/nonidealmhd_params/nambipolar,gammaAD &
        & ,nmagdiffu,etaMD &
        & ,coefad,coefohm, nminitimestep, coefalfven,coefdtohm,nimhdheating_in_flux,nimhdheating_source_term

  ! Checks on non-ideal MHD parameters
  rewind(1)
  read(1,NML=nonidealmhd_params,END=109)

  if(nambipolar.or.nmagdiffu)then
    use_nonideal_mhd = .true.
    if(nimhdheating_in_flux.eqv..false.) nimhdheating_source_term = .true.
    if(nimhdheating_source_term.eqv..true.)  nimhdheating_in_flux     = .false.
  else
    use_nonideal_mhd = .false.
  endif

  if(.not.use_nonideal_mhd) return

  if(myid==1) then
    write(*,*)'!!!!!!!!!!!!!!!  Non Ideal MHD   !!!!!!!!!!!!!!!!'
    write(*,*)'Non ideal MHD parameters'
    if(nambipolar) then
      write(*,*)'Ambipolar diffusion switched ON'
      if(nminitimestep) then
        write(*,*)'Mini time step switched ON'
        write(*,*)'Mini time step coefficient',coefalfven
      else
        write(*,*)'Mini time step switched OFF'
      endif
    else
      write(*,*)'Ambipolar diffusion switched OFF'
    endif

    if(nmagdiffu)then
      write(*,*)'Magnetic diffusion switched ON'
    else
      write(*,*)'Magnetic diffusion switched OFF'
    endif

  endif

  109 continue

end subroutine read_nimhd_params
