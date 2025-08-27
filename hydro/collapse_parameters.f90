module collapse_parameters
  use amr_parameters,only:dp

  ! Cloud parameters
  real(dp)::delta_rho=0.0         ! m=2 density perturbation amplitude
  real(dp)::alpha_dense_core=0.54 ! thermal-to-gravitational energy ratio
  real(dp)::beta_dense_core=0.08  ! rotational-to-gravitational energy ratio
  real(dp)::crit_dense_core=0.0   ! 1/mu for Bfield strength
  real(dp)::theta_mag=0.0         ! angle in degrees for rotation misalignment between Bfield and rotation
  real(dp)::mass_c=1.0            ! mass of the cloud in solar masses
  real(dp)::Mach=0.0              ! Mach number of the cloud
end module collapse_parameters
