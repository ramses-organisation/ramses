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

  ! derived cloud parameters
  real(dp)::r0      ! radius
  real(dp)::d0      ! density
  real(dp)::omega0  ! rotation
  real(dp)::p0      ! pressure
  real(dp)::B0      ! vertical magnetic field

end module collapse_parameters

module collapse_commons
  use amr_parameters,only:dp

  ! Initial turbulence from file
  real(dp)::vx_tot,vy_tot,vz_tot,v_rms
  integer::count_vrms
  integer::n_size
  real(dp),dimension(1:3,1:100,1:100,1:100)::q_idl

end module collapse_commons
