subroutine read_turb_params(nml_ok)
  use amr_commons
  use turb_commons
  implicit none
  logical, intent(inout) ::nml_ok

  !--------------------------------------------------
  ! Namelist definitions
  !--------------------------------------------------
  namelist/turb_params/turb, driven_turb, turb_seed, turb_type, instant_turb, comp_frac,&
       & turb_evolving, forcing_power_spectrum, turb_T, turb_Ndt, turb_rms,&
       & turb_min_rho, turb_exact_rms, initial_turb, initial_turb_vrms,&
       & initial_turb_spectrum, initial_turb_comp_frac

  !--------------------------------------------------
  ! Read namelist; check variables that have been loaded
  !--------------------------------------------------

  ! Read namelist file
  rewind(1)
  read(1,NML=turb_params,END=87)

  ! Namelist parameter turb has been replaced by driven_turb
  if (turb) then
    driven_turb = .true.
  end if

  if (.NOT. (driven_turb.OR.initial_turb)) return

  ! turb_type is deprecated: 1 and 2 are now the boolean turb_evolving, and 3
  ! was never a driving mode at all - it applied the field once as an initial
  ! velocity, which is what initial_turb does.
  if (turb_type /= -1) then
     select case (turb_type)
     case (1)
        if (myid==1) write (*,*) &
             & "WARNING: turb_type=1 is deprecated, use turb_evolving=.true."
        turb_evolving = .TRUE.
     case (2)
        if (myid==1) write (*,*) &
             & "WARNING: turb_type=2 is deprecated, use turb_evolving=.false."
        turb_evolving = .FALSE.
     case (3)
        if (myid==1) then
           write (*,*) "WARNING: turb_type=3 is deprecated."
           write (*,*) "Use driven_turb=.false., initial_turb=.true. and"
           write (*,*) "initial_turb_vrms=sqrt(ndim)*turb_rms instead."
        end if
        if (driven_turb) then
           driven_turb = .FALSE.
           initial_turb = .TRUE.
           initial_turb_vrms = sqrt(real(ndim,dp)) * turb_rms
           ! turb_type=3 drew from the driving spectrum, so a legacy namelist
           ! expresses its intent through those parameters
           initial_turb_spectrum = forcing_power_spectrum
           initial_turb_comp_frac = comp_frac
        end if
     case default
        write (*,*) "Invalid turbulence type selected! (1 to 3)"
        nml_ok = .FALSE.
     end select
  end if

  if (initial_turb_comp_frac < 0.0_dp .OR. initial_turb_comp_frac > 1.0_dp) then
     write (*,*) "Invalid initial compressive fraction selected! (0.0 to 1.0)"
     nml_ok = .FALSE.
  end if

  if (initial_turb .AND. initial_turb_vrms <= 0.0_dp) then
     write (*,*) "Initial turbulent velocity dispersion must be > 0.0!"
     nml_ok = .FALSE.
  end if

  ! A static field would otherwise carry the scatter of its single draw as a
  ! systematic amplitude bias for the whole run
  if (driven_turb .AND. .NOT.turb_evolving .AND. .NOT.turb_exact_rms) then
     if (myid==1) write (*,*) &
          & "Non-evolving driving: setting turb_exact_rms=.true."
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

  if (driven_turb .AND. turb_rms <= 0.0_dp) then
     write (*,*) "Turbulent forcing rms acceleration must be > 0.0!"
     nml_ok = .FALSE.
  end if

87 continue

end subroutine read_turb_params
