module cr_hydro_commons
  use amr_parameters
  use cr_parameters
  real(dp),allocatable,dimension(:,:)::cruold,crunew ! CR state vector and its update
end module cr_hydro_commons
