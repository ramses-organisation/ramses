module fld_commons
  use amr_parameters, only : dp
  integer,allocatable,dimension(:)::liste_ind
  integer::nb_ind
  real(dp)   ::dt_imp                            ! Implicit timestep               
  logical,allocatable,dimension(:)::in_sink !false -> true if cell within sink radius

end module fld_commons

! Units
module units_commons
  use amr_parameters, only : dp
  real(dp):: scale_E0,scale_kappa,scale_m,P_cal,C_cal
end module units_commons
