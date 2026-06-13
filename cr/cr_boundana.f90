!############################################################
!############################################################
!############################################################
!############################################################
subroutine cr_boundana(x,u,dx,ibound,ncell)
  use amr_parameters, ONLY: dp,ndim,nvector
  use cr_parameters, ONLY: ncr,ncrvars,iCRu,smallcr
  implicit none
  integer ::ibound                          ! Index of boundary region
  integer ::ncell                           ! Number of active cells
  real(dp)::dx                              ! Cell size
  real(dp),dimension(1:nvector,1:ncrvars)::u ! CR conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x   ! Cell center position.
  !================================================================
  ! Analytic boundary state for the separated CR variables.
  ! Positions are in user units: x(i,1:3) are in [0,boxlen]**ndim.
  ! u(i,1:ncrvars) is the CR conservative vector for group g:
  !   energy at iCRu+(ndim+1)*(g-1), the ndim fluxes in the slots after it.
  ! ibound is the namelist boundary-region index.
  !
  ! Default: impose the CR energy floor smallcr and zero flux. Test-specific
  ! imposed boundaries (e.g. the Jiang & Oh streaming waves) override this in
  ! the patch version of this file selected via PATCH=.
  !================================================================
  integer::i,igrp,icrE

  ! Default imposed CR boundary: energy floor, zero flux.
  do i=1,ncell
     do igrp=1,ncr
        icrE=iCRu+(ndim+1)*(igrp-1)
        u(i,icrE)=smallcr
        u(i,icrE+1:icrE+ndim)=0d0
     end do
  end do

  ! Add here, if you wish, some user-defined CR boundary conditions
  ! ........

end subroutine cr_boundana
