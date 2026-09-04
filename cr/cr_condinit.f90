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
  real(dp),dimension(1:nvector,1:ncrvar)::u ! CR conservative variables
  real(dp),dimension(1:nvector,1:ndim  )::x ! Cell center position.
  !================================================================
  ! Generate CR initial conditions; u holds per-group [energy, ndim fluxes].
  ! Default: energy = 0., flux = 0; test ICs override via PATCH.
  !================================================================

  u(1:nn,1:ncrvar)=0d0

  ! Add here, if you wish, some user-defined CR initial conditions
  ! ........

end subroutine cr_condinit
