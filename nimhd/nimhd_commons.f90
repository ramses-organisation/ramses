module nimhd_commons
   use amr_parameters, ONLY:dp
   implicit none

   ! for RAMSES
   real(dp),allocatable,dimension(:,:,:,:,:)::resistivite_chimie ! resistivites chimie
   real(dp),allocatable,dimension(:,:)::resistivite_chimie_res   ! resistivites chimie (TC: used for 1d table, to generalize)

end module nimhd_commons
