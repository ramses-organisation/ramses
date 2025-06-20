module collapse_parameters
  use amr_parameters,only:dp

  ! Cloud parameters
  real(dp)::delta_rho=0.0
  real(dp)::alpha_dense_core=0.54
  real(dp)::beta_dense_core=0.08
  real(dp)::crit_dense_core=0.0 ! 1/mu for Bfield strength

end module collapse_parameters
