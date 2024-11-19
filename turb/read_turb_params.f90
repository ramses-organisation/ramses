#if USE_TURB==1
subroutine read_turb_params(nml_ok)
  use amr_commons
  use turb_commons
  implicit none
  logical, intent(inout) ::nml_ok

  !--------------------------------------------------
  ! Namelist definitions
  !--------------------------------------------------
  namelist/turb_params/turb, turb_seed, turb_type, instant_turb, comp_frac,&
       & forcing_power_spectrum, turb_T, turb_Ndt, turb_rms, turb_min_rho

  !--------------------------------------------------
  ! Read namelist; check variables that have been loaded
  !--------------------------------------------------

  ! Read namelist file
  rewind(1)
  read(1,NML=turb_params,END=87)

  if (.NOT. turb) return

  if (turb_type < 1 .OR. turb_type > 3) then
     write (*,*) "Invalid turbulence type selected! (1 to 3)"
     nml_ok = .FALSE.
  end if

  ! BUG: upon restart, turb_type 2 gives the wrong rms.
  if (turb_type == 2) then
     write (*,*) "Turbulence type 2 is bugged. Please select 1 instead."
     nml_ok = .FALSE.
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

  if (turb_rms <= 0.0_dp) then
     write (*,*) "Turbulent forcing rms acceleration must be > 0.0!"
     nml_ok = .FALSE.
  end if

87 continue

end subroutine read_turb_params
#endif
