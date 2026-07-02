!############################################################
!############################################################
!############################################################
!############################################################
subroutine cr_boundana(x,u,dx,ibound,ncell)
  use amr_parameters, ONLY: dp,ndim,nvector
  use cr_parameters, ONLY: ncr_groups,ncrvar,Ecr_idx,smallcr
  implicit none
  integer ::ibound                          ! Index of boundary region
  integer ::ncell                           ! Number of active cells
  real(dp)::dx                              ! Cell size
  real(dp),dimension(1:nvector,1:ncrvar)::u ! CR conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x   ! Cell center position.
  !================================================================
  ! Analytic CR boundary state; u holds per-group [energy, ndim fluxes],
  ! ibound = boundary-region index. Default: energy = 0., flux = 0.
  !================================================================
  integer::i,igrp,icrE

  do i=1,ncell
     do igrp=1,ncr_groups
        icrE=Ecr_idx(igrp)
        u(i,icrE)=0d0
        u(i,icrE+1:icrE+ndim)=0d0
     end do
  end do

  ! Add here, if you wish, some user-defined CR boundary conditions
  ! ........

end subroutine cr_boundana
