module nimhd_parameters
  use amr_parameters

  logical :: nambipolar = .false. ! flag to activate ambipolar diffusion
  logical :: nmagdiffu  = .false. ! flag to activate magnetic  diffusion
  logical :: use_nonideal_mhd = .false.     ! true if any of the non-ideal MHD effects is used

  logical :: nimhdheating_in_flux = .true. ! flag to activate the non-ideal mhd energy fluxes
  logical :: nimhdheating_source_term = .false. ! flag to activate the non-ideal mhd heating as a source term

  ! Fixed resistivity coefficients (in code units)
  real(dp):: gammaAD=1.0d0            ! ambipolar diffusion coefficient (beta = 1/(gammaAD*rho))
  real(dp):: etaMD=1d0                ! Ohmic (magnetic) diffusion coefficient

  ! timestep regulation
  real(dp):: coefad = 0.1d0    ! CFL condition for ambipolar diffusion
  real(dp):: coefohm = 0.05d0  ! CFL condition for ohmic dissipation
  logical :: nimhd_dt_cap = .false.  ! flag to activate timestep limitation hack (artificailly increase it)
  real(dp):: frac_dt_cap_ad = 1d-10  ! If capping dt: maximal ratio between ambipolar diffusion timestep and ideal MHD timestep.
  real(dp):: frac_dt_cap_ohm = 1d-10 ! same as frac_dt_cap_ad but for Ohmic dissipation.

end module nimhd_parameters
