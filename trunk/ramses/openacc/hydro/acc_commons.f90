module acc_commons
  use amr_parameters
  integer,parameter::ncube=1000
  real(dp),allocatable,dimension(:)::uloc_tot,gloc_tot
  real(dp),allocatable,dimension(:)::divu_acc,enew_acc
   
  real(dp),dimension(:,:,:,:,:),allocatable::qin 
  real(dp),dimension(:,:,:,:),allocatable::cin
  real(dp),dimension(:,:,:,:,:,:),allocatable::dq
  real(dp),dimension(:,:,:,:,:,:),allocatable::qm
  real(dp),dimension(:,:,:,:,:,:),allocatable::qp
  
  integer,dimension(1:ncube)::ucount0,gcount0,igrid_acc
  real(dp)::oneontwotondim
  
  integer,dimension(0:7,0:1,0:11)::state_diagram
  integer,dimension(1:8,1:6)::ggg_acc,hhh_acc
  
  real(dp) :: dx_acc, dtdx, dtdy, dtdz
  
  logical::copy_amr
 
end module acc_commons