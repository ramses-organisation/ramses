module hydro_commons
  use amr_parameters
  use hydro_parameters
  real(dp),allocatable,dimension(:,:)::uold,unew ! State vector and its update
  real(dp),allocatable,dimension(:)::pstarold,pstarnew ! Stellar momentum and its update
  real(dp),allocatable,dimension(:)::divu,enew ! Non conservative variables
  real(dp),allocatable,dimension(:)::rho_eq,p_eq ! Strict hydrostatic equilibrium
  real(dp)::mass_tot=0.0D0,mass_tot_0=0.0D0
  real(dp)::e_tot_0=0d0,e_tot=0d0
  real(dp)::mom_tot_0=0d0,mom_tot=0d0
  real(dp)::lor_max=1d0
end module hydro_commons

module const
  use amr_parameters
  real(dp)::bigreal = 1.0e+30
  real(dp)::zero = 0.0
  real(dp)::one = 1.0
  real(dp)::two = 2.0
  real(dp)::three = 3.0
  real(dp)::four = 4.0
  real(dp)::two3rd = 0.6666666666666667
  real(dp)::half = 0.5
  real(dp)::third = 0.33333333333333333
  real(dp)::forth = 0.25
  real(dp)::sixth = 0.16666666666666667
end module const
