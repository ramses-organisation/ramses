module nimhd_commons
   use amr_parameters, ONLY:dp
   implicit none

   ! for RAMSES
   real(dp),allocatable,dimension(:,:,:,:,:)::resistivite_chimie ! resistivites chimie

end module nimhd_commons
