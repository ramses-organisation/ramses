!############################################################
!############################################################
!############################################################
!############################################################
subroutine cr_boundana(x,u,dx,ibound,ncell)
  use amr_parameters, ONLY: dp,ndim,nvector,boxlen
  use amr_commons,    ONLY: t
  use cr_parameters,  ONLY: ncrvar,Ecr_idx,gamma_cr,cr_boundana_kind,cr_boundary_u
  implicit none
  integer ::ibound                          ! Index of boundary region
  integer ::ncell                           ! Number of active cells
  real(dp)::dx                              ! Cell size
  real(dp),dimension(1:nvector,1:ncrvar)::u ! CR conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x   ! Cell center position.
  !================================================================
  ! CR state imposed on a bound_type=3 boundary; u holds per-group
  ! [energy, ndim fluxes]. Default: the namelist state cr_boundary_u(ibound,:).
  ! cr_boundana_kind selects a time-dependent analytic override.
  !================================================================
  integer::i

  do i=1,ncell
     u(i,1:ncrvar)=cr_boundary_u(ibound,1:ncrvar)
  end do

  select case(trim(cr_boundana_kind))

  case('jiang_411_triangular')
     ! Streaming wave entering along x: E_cr = 2 + gamma_cr*t - |x-boxlen/2|,
     ! F_cr = +/- gamma_cr*E_cr, pointing inward on each side.
     do i=1,ncell
        u(i,Ecr_idx(1))=2d0+gamma_cr(1)*t-abs(x(i,1)-boxlen*0.5d0)
        if(x(i,1)<boxlen*0.5d0)then
           u(i,Ecr_idx(1)+1)=-gamma_cr(1)*u(i,Ecr_idx(1))
        else
           u(i,Ecr_idx(1)+1)= gamma_cr(1)*u(i,Ecr_idx(1))
        endif
        if(ncrvar>2) u(i,Ecr_idx(1)+2:Ecr_idx(1)+ncrvar-1)=0d0
     end do

  end select

end subroutine cr_boundana
