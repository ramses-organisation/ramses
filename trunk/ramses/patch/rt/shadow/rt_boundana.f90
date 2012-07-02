!############################################################
!############################################################
!############################################################
!############################################################
subroutine rt_boundana(x,u,dx,ibound,ncell)
  use amr_parameters, ONLY: dp,ndim,nvector
  use rt_parameters, ONLY: nrtvar,rt_boundary_var,nPacs,iPac,rt_c
  implicit none
  integer ::ibound                          ! Index of boundary region
  integer ::ncell                           ! Number of active cells
  real(dp)::dx                              ! Cell size
  real(dp),dimension(1:nvector,1:nrtvar)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x   ! Cell center position.
  !================================================================
  ! This routine generates boundary conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): N, U(i,2:ndim+1): Fx,Fy,Fz.
  ! U is in user units.
  ! ibound is the index of the boundary region defined in the namelist.
  !================================================================
  integer::ivar,i,iP
  real(dp)::scale_Np, scale_Fp
  real(dp),parameter::Fpx(3)=(/ 0.447d6, 0.494d6, 0.059d6 /)

  call rt_units(scale_Np, scale_Fp)
  do ip=1,nPacs
     do i=1,ncell
        ! Photon density:
        u(i,iPac(ip))   = Fpx(ip) / scale_Fp /rt_c
        ! Photon x-flux:
        u(i,iPac(ip)+1) = Fpx(ip) / scale_Fp
        ! Photon y-flux:
        u(i,iPac(ip)+2) = 0
     end do
  end do

end subroutine rt_boundana
