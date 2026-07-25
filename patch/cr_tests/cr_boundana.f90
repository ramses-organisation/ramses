!############################################################
!############################################################
!############################################################
!############################################################
subroutine cr_boundana(x,u,dx,ibound,ncell)
  use amr_parameters, ONLY: dp,ndim,nvector,boxlen
  use amr_commons,    ONLY: t
  use cr_parameters,  ONLY: ncr_groups,ncrvar,Ecr_idx,cr_efloor,gamma_cr,jiang_test,cr_boundary_u
  implicit none
  integer ::ibound                          ! Index of boundary region
  integer ::ncell                           ! Number of active cells
  real(dp)::dx                              ! Cell size
  real(dp),dimension(1:nvector,1:ncrvar)::u ! CR conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x   ! Cell center position.
  !================================================================
  ! Jiang & Oh (2018) CR test suite -- imposed (bound_type=3) CR boundary; sets
  ! only u(1:ncell,1:ncrvar). Default = floor + zero flux, test branches override.
  !================================================================
  integer::i,igrp,icrE

  ! Default imposed CR boundary: energy floor, zero flux.
  do i=1,ncell
     do igrp=1,ncr_groups
        icrE=Ecr_idx(igrp)
        u(i,icrE)=cr_efloor
        u(i,icrE+1:icrE+ndim)=0d0
     end do
  end do

  ! Test-specific analytic boundary: 411_triangular streaming wave entering
  ! along x:  E_cr = 2 + gamma_cr*t - |x - boxlen/2|,  F_cr = +/- gamma_cr*E_cr
  ! (flux points inward on each side). Higher flux components vanish.
  if(trim(jiang_test)=='411_triangular')then
     do i=1,ncell
        u(i,Ecr_idx(1))=2d0+gamma_cr(1)*t-abs(x(i,1)-boxlen*0.5d0)
        if(x(i,1)<boxlen*0.5d0)then
           u(i,Ecr_idx(1)+1)=-gamma_cr(1)*u(i,Ecr_idx(1))
        else
           u(i,Ecr_idx(1)+1)= gamma_cr(1)*u(i,Ecr_idx(1))
        endif
        if(ncrvar>2) u(i,Ecr_idx(1)+2:Ecr_idx(1)+ncrvar-1)=0d0
     end do
  endif

  ! Two-pressure shock (tp-*): imposed CR boundary held at the namelist
  ! cr_boundary_u(ibound,:) -- energy at index 1, CR flux at 2:ncrvar.
  if(trim(jiang_test)=='tp-nostream' .or. trim(jiang_test)=='tp-stream-va075' &
       & .or. trim(jiang_test)=='tp-stream-va15')then
     do i=1,ncell
        u(i,1:ncrvar)=cr_boundary_u(ibound,1:ncrvar)
     end do
  endif

end subroutine cr_boundana
