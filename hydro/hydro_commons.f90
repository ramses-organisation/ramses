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

#if USE_FLD==1
  integer,allocatable,dimension(:)::liste_ind
  integer::nb_ind
  real(dp)   ::dt_imp                            ! Implicit timestep               
  logical,allocatable,dimension(:)::in_sink !false -> true if cell within sink radius
#endif
end module hydro_commons

module const
  use amr_parameters
  real(dp)::bigreal = 1.0d+30
  real(dp)::zero = 0
  real(dp)::one = 1
  real(dp)::two = 2
  real(dp)::three = 3
  real(dp)::four = 4
  real(dp)::five = 5
  real(dp)::six = 6
  real(dp)::seven = 7
  real(dp)::eight = 8
  real(dp)::nine = 9
  real(dp)::ten = 10
  real(dp)::two3rd = 2/3d0
  real(dp)::half = 1/2d0
  real(dp)::third = 1/3d0
  real(dp)::forth = 1/4d0
  real(dp)::sixth = 1/6d0
end module const

#if USE_FLD==1
! Units
module units_commons
  use amr_parameters, only : dp
  real(dp):: scale_E0,scale_m,P_cal,C_cal
end module units_commons
#endif
