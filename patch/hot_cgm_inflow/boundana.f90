!############################################################
!############################################################
!############################################################
!############################################################
subroutine boundana(x,u,dx,ibound,ncell)
  use amr_parameters, ONLY: dp,ndim,nvector,boxlen
  use hydro_parameters, ONLY: nvar,nvar_all,boundary_var,gamma
  use cgm_commons
  use constants, only: kB, mH
  implicit none
  integer ::ibound                        ! Index of boundary region
  integer ::ncell                         ! Number of active cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar_all)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates analytical boundary conditions for the
  ! hot rotating CGM inflow problem.
  !
  ! Positions are in user units: x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector:
  ! U(i,1): d, U(i,2:ndim+1): d.u,d.v,d.w and U(i,ndim+2): E.
  !================================================================
  integer::i,ivar
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp)::xx,yy,zz,rr,RR_cyl,theta
  real(dp)::rho_cgs,T_cgs,P_cgs,v_phi,v_r
  real(dp)::vx,vy,vz,cos_phi,sin_phi
  real(dp)::x_center,y_center,z_center
  real(dp)::r_soft,rho,P,ekin

  ! Get unit conversions
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Box center
  x_center = boxlen / 2.0d0
  y_center = boxlen / 2.0d0
  z_center = boxlen / 2.0d0

  ! Softening radius (code units)
  r_soft = cgm_r_soft / scale_l

  do i=1,ncell
     ! Position relative to box center
     xx = x(i,1) - x_center
     yy = x(i,2) - y_center
     zz = x(i,3) - z_center

     ! Spherical and cylindrical radii
     rr = sqrt(xx**2 + yy**2 + zz**2 + r_soft**2)
     RR_cyl = sqrt(xx**2 + yy**2 + r_soft**2)

     ! Polar angle
     theta = acos(zz / max(rr, 1d-10))

     ! Get analytical CGM state
     call get_cgm_state(rr*scale_l, RR_cyl*scale_l, theta, &
          rho_cgs, T_cgs, P_cgs, v_r, v_phi)

     ! Convert to code units
     rho = rho_cgs / scale_d
     P = P_cgs / (scale_d * scale_v**2)

     ! Velocity components
     if (RR_cyl * scale_l > cgm_r_soft) then
        cos_phi = xx / RR_cyl
        sin_phi = yy / RR_cyl
     else
        cos_phi = 1.0d0
        sin_phi = 0.0d0
     endif

     vx = -v_phi * sin_phi / scale_v + v_r * (xx/rr) / scale_v
     vy =  v_phi * cos_phi / scale_v + v_r * (yy/rr) / scale_v
     vz =  v_r * (zz/rr) / scale_v

     ! Fill conservative variables
     ! Density
     u(i,1) = rho

     ! Momentum
     u(i,2) = rho * vx
#if NDIM>1
     u(i,3) = rho * vy
#endif
#if NDIM>2
     u(i,4) = rho * vz
#endif

     ! Total energy = kinetic + thermal
     ekin = 0.5d0 * rho * (vx**2 + vy**2 + vz**2)
     u(i,ndim+2) = ekin + P / (gamma - 1.0d0)

  end do

  ! Passive scalars (if any)
#if NVAR>NDIM+2
  do ivar=ndim+3,nvar_all
     do i=1,ncell
        u(i,ivar) = 0.0d0
     end do
  end do
#endif

end subroutine boundana


!============================================================
! Compute analytical state at reference cell position
! This is used for non-reflecting boundary conditions
!============================================================
subroutine boundana_ref(x,u,dx,ibound,ncell)
  use amr_parameters, ONLY: dp,ndim,nvector,boxlen
  use hydro_parameters, ONLY: nvar,nvar_all,gamma
  use cgm_commons
  use constants, only: kB, mH
  implicit none
  integer ::ibound
  integer ::ncell
  real(dp)::dx
  real(dp),dimension(1:nvector,1:nvar_all)::u
  real(dp),dimension(1:nvector,1:ndim)::x
  !================================================================
  ! Same as boundana but for reference cell positions
  ! Used to compute perturbations for non-reflecting BCs
  !================================================================

  ! Just call boundana - same computation
  call boundana(x,u,dx,ibound,ncell)

end subroutine boundana_ref
