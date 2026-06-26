module hydro_commons
  use amr_parameters
  use hydro_parameters
  real(dp),allocatable,dimension(:,:)::uold,unew ! State vector and its update
  real(dp),allocatable,dimension(:,:)::fluxes ! Mass flux on the faces of the cells
  real(dp),allocatable,dimension(:)::pstarold,pstarnew ! Stellar momentum and its update
  real(dp),allocatable,dimension(:)::divu,enew ! Non conservative variables
  real(dp),allocatable,dimension(:)::rho_eq,p_eq ! Strict hydrostatic equilibrium
  real(dp)::mass_tot=0,mass_tot_0=0
  real(dp)::ana_xmi,ana_xma,ana_ymi,ana_yma,ana_zmi,ana_zma
  integer::nbins

  ! EOS tables
  real(dp),allocatable,dimension(:,:)::Rho_eos,Ener_eos,Temp_eos,P_eos,Cs_eos,S_eos,eint_eos
  real(dp),allocatable,dimension(:,:)::xH_eos, xH2_eos, xHe_eos,xHep_eos,Cv_eos,Dc_eos

end module hydro_commons

module const
  use amr_parameters
  real(dp)::bigreal = 1.0d+30
  real(dp)::zero = 0
  real(dp)::one = 1
  real(dp)::two = 2
  real(dp)::three = 3
  real(dp)::four = 4
  real(dp)::two3rd = 2/3d0
  real(dp)::half = 1/2d0
  real(dp)::third = 1/3d0
  real(dp)::forth = 1/4d0
  real(dp)::sixth = 1/6d0
end module const
