!================================================================
!================================================================
!================================================================
!================================================================
subroutine cr_condinit(x,u,dx,nn,ilevel)
  use amr_parameters
  use cr_parameters
  implicit none
  integer ::nn                              ! Number of cells
  integer:: ilevel                          ! Refinement level
  real(dp)::dx                              ! Cell size
  real(dp),dimension(1:nvector,1:ncrvars)::u ! CR conservative variables
  real(dp),dimension(1:nvector,1:ndim  )::x ! Cell center position.
  !================================================================
  ! This routine generates CR initial conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! u(i,1:ncrvars) is the CR conservative vector for group g:
  !   energy at iCRu+(ndim+1)*(g-1), the ndim fluxes in the slots after it.
  !   (iCRu=1; CR indexing never references nvar.)
  ! u(:,:) are in user units.
  !
  ! Default: set every CR energy to the smallcr floor and every CR flux to
  ! zero. Test-specific CR initial states (e.g. the Jiang & Oh Gaussian
  ! pulse) override this in the patch version of this file selected via
  ! PATCH=.
  !================================================================
  integer::i,igrp,icrE

  ! Built-in default CR initial condition: energy floor, zero flux.
  do i=1,nn
     do igrp=1,ncr
        icrE=iCRu+(ndim+1)*(igrp-1)
        u(i,icrE)=smallcr
        u(i,icrE+1:icrE+ndim)=0d0
     end do
  end do

  ! Add here, if you wish, some user-defined CR initial conditions
  ! ........

end subroutine cr_condinit
