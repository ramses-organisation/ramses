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

  ! turb_type 2 used to give the wrong rms upon restart: init_turb interpolated
  ! between the two stored fields, and on a restart t falls part way through the
  ! interval, so the static field came out as a blend of two independent fields,
  ! whose rms is below turb_rms. init_turb now takes a single field for type 2
  ! and never interpolates it. Covered by the turb/driving-fixed test, whose
  ! bottom panel checks that the rms stays constant across a restart.

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

  ! instant_turb evolves the field for several autocorrelation times at startup
  ! so that an evolving run (turb_type=1) begins already saturated rather than
  ! climbing from 14% of the target amplitude over the first turb_T. For fixed
  ! (2) and decaying (3) turbulence a single draw is renormalised onto the
  ! requested amplitude in init_turb, so the spin-up (500 FFTs) achieves nothing
  ! and is switched off here.
  if (turb_type /= 1 .AND. instant_turb) then
     if (myid==1) write (*,*) &
          & "instant_turb only affects turb_type=1; disabling it for this run."
     instant_turb = .FALSE.
  end if

87 continue

end subroutine read_turb_params
