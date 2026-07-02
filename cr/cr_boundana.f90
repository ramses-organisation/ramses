!############################################################
!############################################################
!############################################################
!############################################################
subroutine cr_boundana(x,u,dx,ibound,ncell)
  use amr_parameters, ONLY: dp,ndim,nvector
  use cr_parameters, ONLY: ncrvar
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

  u(1:ncell,1:ncrvar)=0d0

  ! Add here, if you wish, some user-defined CR boundary conditions
  ! ........

end subroutine cr_boundana
