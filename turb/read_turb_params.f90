subroutine read_turb_params(nml_ok)
  use amr_commons
  use turb_commons
  implicit none
  logical, intent(inout) ::nml_ok

  !--------------------------------------------------
  ! Namelist definitions
  !--------------------------------------------------
  namelist/turb_params/turb, turb_seed, turb_type, instant_turb, comp_frac,&
       & forcing_power_spectrum, turb_T, turb_Ndt, turb_rms, turb_min_rho,&
       & turb_exact_rms, initial_turb, initial_turb_vrms,&
       & initial_turb_spectrum, initial_turb_comp_frac

  !--------------------------------------------------
  ! Read namelist; check variables that have been loaded
  !--------------------------------------------------

  ! Read namelist file
  rewind(1)
  read(1,NML=turb_params,END=87)

  if (.NOT. (turb.OR.initial_turb)) return

  ! turb_type=3 was never a driving mode: it applied the field once as an
  ! initial velocity. Translate it onto the dedicated parameters.
  if (turb .AND. turb_type == 3) then
     if (myid==1) then
        write (*,*) "WARNING: turb_type=3 is deprecated."
        write (*,*) "Use turb=.false., initial_turb=.true. and"
        write (*,*) "initial_turb_vrms=sqrt(ndim)*turb_rms instead."
     end if
     turb = .FALSE.
     initial_turb = .TRUE.
     initial_turb_vrms = sqrt(real(ndim,dp)) * turb_rms
  end if

  ! The initial field defaults to the same spectrum and compressive fraction as
  ! the driving
  if (initial_turb_spectrum == '') initial_turb_spectrum = forcing_power_spectrum
  if (initial_turb_comp_frac < 0.0_dp) initial_turb_comp_frac = comp_frac

  if (initial_turb_comp_frac > 1.0_dp) then
     write (*,*) "Invalid initial compressive fraction selected! (0.0 to 1.0)"
     nml_ok = .FALSE.
  end if

  if (initial_turb .AND. initial_turb_vrms <= 0.0_dp) then
     write (*,*) "Initial turbulent velocity dispersion must be > 0.0!"
     nml_ok = .FALSE.
  end if

  if (turb .AND. (turb_type < 1 .OR. turb_type > 2)) then
     write (*,*) "Invalid turbulence type selected! (1 or 2)"
     nml_ok = .FALSE.
  end if

  if (turb .AND. turb_type==2 .AND. .NOT.turb_exact_rms) then
     write (*,*) "Turb_type 2: Setting turb_exact_rms=.true."
     turb_exact_rms = .TRUE.
  end if

  if (comp_frac < 0.0_dp .OR. comp_frac > 1.0_dp) then
     write (*,*) "Invalid compressive fraction selected! (0.0 to 1.0)"
     nml_ok = .FALSE.
  end if

  if (turb_T <= 0.0_dp) then
     write (*,*) "Turbulent autocorrelation time must be > 0!"
     nml_ok = .FALSE.
  end if

  if (turb_Ndt <= 0) then
     write (*,*) "Number of timesteps per autocorrelation time must be > 0!"
     nml_ok = .FALSE.
  end if

  if (turb .AND. turb_rms <= 0.0_dp) then
     write (*,*) "Turbulent forcing rms acceleration must be > 0.0!"
     nml_ok = .FALSE.
  end if

87 continue

end subroutine read_turb_params
