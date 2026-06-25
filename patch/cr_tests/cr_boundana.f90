!############################################################
!############################################################
!############################################################
!############################################################
subroutine cr_boundana(x,u,dx,ibound,ncell)
  use amr_parameters, ONLY: dp,ndim,nvector,boxlen
  use amr_commons,    ONLY: t
  use cr_parameters,  ONLY: ncr,ncrvars,iCRu,smallcr,gamma_cr,jiang_test,crmom_bound
  implicit none
  integer ::ibound                          ! Index of boundary region
  integer ::ncell                           ! Number of active cells
  real(dp)::dx                              ! Cell size
  real(dp),dimension(1:nvector,1:ncrvars)::u ! CR conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x   ! Cell center position.
  !================================================================
  ! Jiang & Oh (2018) cosmic-ray test suite -- separated-CR analytic boundary.
  !
  ! Patch override of cr/cr_boundana.f90: sets ONLY the CR conservative
  ! variables u(1:ncell,1:ncrvars) in an imposed (bound_type=3) boundary
  ! region. ibound is the namelist boundary-region index.
  !
  ! Default: impose the CR energy floor smallcr and zero flux (identical to the
  ! base cr_boundana). The 411_triangular test overrides this with a
  ! time-dependent analytic streaming wave; all other jiang_test values keep
  ! the default, so this patch is inert for them.
  !
  ! Ported from ramses_cral patch/jiang_tests/boundana.f90, applying the central
  ! transformation: cral wrote the embedded u(:,icrU...) with icrU=nvar+4; here
  ! we write the separated CR buffer u(:,iCRu...) with iCRu=1. gamma_cr/t/boxlen
  ! and the wave formula are byte-faithful to cral.
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

  ! Test-specific analytic boundary: 411_triangular streaming wave entering
  ! along x:  E_cr = 2 + gamma_cr*t - |x - boxlen/2|,  F_cr = +/- gamma_cr*E_cr
  ! (flux points inward on each side). Higher flux components vanish.
  if(trim(jiang_test)=='411_triangular')then
     do i=1,ncell
        u(i,iCRu)=2d0+gamma_cr(1)*t-abs(x(i,1)-boxlen*0.5d0)
        if(x(i,1)<boxlen*0.5d0)then
           u(i,iCRu+1)=-gamma_cr(1)*u(i,iCRu)
        else
           u(i,iCRu+1)= gamma_cr(1)*u(i,iCRu)
        endif
        if(ncrvars>2) u(i,iCRu+2:iCRu+ncrvars-1)=0d0
     end do
  endif

  ! Two-pressure shock (tp_*): imposed CR boundary held at the namelist
  ! crmom_bound(ibound,:) -- energy at index 1, CR flux at 2:ncrvars -- the
  ! per-boundary analogue of crmom_region (cral's crmom_bound, separated here).
  if(trim(jiang_test)=='tp_nostream' .or. trim(jiang_test)=='tp_stream_va075' &
       & .or. trim(jiang_test)=='tp_stream_va15')then
     do i=1,ncell
        u(i,1:ncrvars)=crmom_bound(ibound,1:ncrvars)
     end do
  endif

end subroutine cr_boundana
