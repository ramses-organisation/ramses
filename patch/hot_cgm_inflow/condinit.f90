!================================================================
!================================================================
!================================================================
!================================================================
subroutine condinit(x,u,dx,nn)
  use amr_commons
  use hydro_parameters
  use cgm_commons
  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! Hot rotating CGM inflow initial conditions
  ! Based on Stern et al. (2024) MNRAS 530, 1711
  !
  ! Setup: Hot (~10^6 K) rotating CGM in quasi-hydrostatic equilibrium
  ! with an isothermal gravitational potential (v_c = const).
  ! The CGM has net angular momentum and slowly inflows due to cooling.
  !================================================================
  real(dp),dimension(1:nvector,1:nvar),save::q   ! Primitive variables
  integer::i,ivar
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp)::xx,yy,zz,rr,RR_cyl,theta
  real(dp)::rho_cgs,T_cgs,P_cgs,v_phi,v_r
  real(dp)::vx,vy,vz,cos_phi,sin_phi
  real(dp)::x_center,y_center,z_center
  real(dp)::r_soft,n_H,tcool,tff,cs

  ! Get unit conversions
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Box center
  x_center = boxlen / 2.0d0
  y_center = boxlen / 2.0d0
  z_center = boxlen / 2.0d0

  ! Softening radius to avoid singularity at center (in code units)
  r_soft = cgm_r_soft / scale_l  ! Convert from cm to code units

  do i=1,nn
     ! Position relative to box center (code units)
     xx = x(i,1) - x_center
     yy = x(i,2) - y_center
     zz = x(i,3) - z_center

     ! Spherical and cylindrical radii
     rr = sqrt(xx**2 + yy**2 + zz**2 + r_soft**2)
     RR_cyl = sqrt(xx**2 + yy**2 + r_soft**2)

     ! Polar angle from rotation axis (z-axis)
     theta = acos(zz / max(rr, 1d-10))

     ! Get analytical CGM state at this position
     call get_cgm_state(rr*scale_l, RR_cyl*scale_l, theta, &
          rho_cgs, T_cgs, P_cgs, v_r, v_phi)

     ! Convert to code units
     ! Density
     q(i,1) = rho_cgs / scale_d

     ! Velocities: convert v_phi (azimuthal) and v_r (radial) to Cartesian
     ! v_phi is in the x-y plane, perpendicular to R_cyl
     if (RR_cyl * scale_l > cgm_r_soft) then
        cos_phi = xx / (RR_cyl)
        sin_phi = yy / (RR_cyl)
     else
        cos_phi = 1.0d0
        sin_phi = 0.0d0
     endif

     ! Azimuthal velocity: v_phi points in (-sin_phi, cos_phi, 0) direction
     ! Radial velocity: v_r points in (xx, yy, zz)/rr direction
     vx = -v_phi * sin_phi / scale_v + v_r * (xx/rr) / scale_v
     vy =  v_phi * cos_phi / scale_v + v_r * (yy/rr) / scale_v
     vz =  v_r * (zz/rr) / scale_v

     q(i,2) = vx
     q(i,3) = vy
     q(i,4) = vz

     ! Pressure (code units: P / (scale_d * scale_v^2))
     q(i,neul) = P_cgs / (scale_d * scale_v**2)

  end do

  ! Initialize passive scalars to zero if present
#if NVAR>NHYDRO+NENER
  do ivar=nhydro+1+nener,nvar
     q(1:nn,ivar) = 0.0d0
  end do
#endif

  ! Convert primitive to conservative variables
  ! density -> density
  u(1:nn,1)=q(1:nn,1)
  ! velocity -> momentum
  u(1:nn,2)=q(1:nn,1)*q(1:nn,2)
#if NDIM>1
  u(1:nn,3)=q(1:nn,1)*q(1:nn,3)
#endif
#if NDIM>2
  u(1:nn,4)=q(1:nn,1)*q(1:nn,4)
#endif
  ! kinetic energy
  u(1:nn,neul)=0.0d0
  u(1:nn,neul)=u(1:nn,neul)+0.5d0*q(1:nn,1)*q(1:nn,2)**2
#if NDIM>1
  u(1:nn,neul)=u(1:nn,neul)+0.5d0*q(1:nn,1)*q(1:nn,3)**2
#endif
#if NDIM>2
  u(1:nn,neul)=u(1:nn,neul)+0.5d0*q(1:nn,1)*q(1:nn,4)**2
#endif
  ! thermal pressure -> total fluid energy
  u(1:nn,neul)=u(1:nn,neul)+q(1:nn,neul)/(gamma-1.0d0)
#if NENER>0
  ! radiative pressure -> radiative energy
  ! radiative energy -> total fluid energy
  do ivar=1,nener
     u(1:nn,nhydro+ivar)=q(1:nn,nhydro+ivar)/(gamma_rad(ivar)-1.0d0)
     u(1:nn,neul)=u(1:nn,neul)+u(1:nn,nhydro+ivar)
  enddo
#endif
#if NVAR>NHYDRO+NENER
  ! passive scalars
  do ivar=nhydro+1+nener,nvar
     u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)
  end do
#endif

end subroutine condinit


!================================================================
! Compute analytical CGM state at given position
!================================================================
subroutine get_cgm_state(r_cgs, R_cyl_cgs, theta, rho, T, P, v_r, v_phi)
  use amr_parameters, only: dp
  use cgm_commons
  use constants, only: kB, mH
  implicit none

  real(dp), intent(in)  :: r_cgs      ! Spherical radius [cm]
  real(dp), intent(in)  :: R_cyl_cgs  ! Cylindrical radius [cm]
  real(dp), intent(in)  :: theta      ! Polar angle [rad]
  real(dp), intent(out) :: rho        ! Density [g/cm^3]
  real(dp), intent(out) :: T          ! Temperature [K]
  real(dp), intent(out) :: P          ! Pressure [erg/cm^3]
  real(dp), intent(out) :: v_r        ! Radial velocity [cm/s]
  real(dp), intent(out) :: v_phi      ! Azimuthal velocity [cm/s]

  real(dp) :: r_kpc, n_H, tcool, tff, cs2, mu
  real(dp) :: Rc_max_cm, v_c_cgs, Lambda_cgs
  real(dp) :: eps  ! Rotation correction parameter

  ! Mean molecular weight for fully ionized gas
  mu = cgm_mu

  ! Convert parameters to CGS
  v_c_cgs = cgm_v_c * 1d5           ! km/s -> cm/s
  Rc_max_cm = cgm_Rc_max * kpc2cm   ! kpc -> cm
  Lambda_cgs = cgm_Lambda           ! Already in erg cm^3/s

  ! Radius in kpc for scaling relations
  r_kpc = r_cgs / kpc2cm

  ! Sound speed squared for virial temperature (Eq. 10 in Stern+24)
  ! c_s^2 = (10/9) * v_c^2
  cs2 = (10.0d0/9.0d0) * v_c_cgs**2

  ! Temperature: T = (10/9) * mu * m_H * v_c^2 / k_B
  ! T ~ 2.0 x 10^6 K for v_c = 200 km/s
  T = mu * mH * cs2 / kB

  ! Hydrogen number density from Eq. 10 (scaled by Mdot and Lambda)
  ! n_H = 0.8 x 10^-3 * r_10^(-1.5) * v_c,200 * Mdot_1^0.5 * Lambda_-22^-0.5 cm^-3
  ! We use the relation that sets n_H for a given Mdot
  n_H = cgm_n_H_0 * (r_kpc / 10.0d0)**(-1.5d0)

  ! Ensure minimum density
  n_H = max(n_H, cgm_n_H_floor)

  ! Mass density
  rho = n_H * mu * mH

  ! Pressure from ideal gas law
  P = n_H * kB * T

  ! Cooling time (Eq. 4): t_cool = (3/2) P / (n_H^2 Lambda)
  tcool = 1.5d0 * P / (n_H**2 * Lambda_cgs)

  ! Free-fall time: t_ff = sqrt(2) * r / v_c
  tff = sqrt(2.0d0) * r_cgs / v_c_cgs

  ! Radial inflow velocity (Eq. 10): v_r = -r / t_cool
  v_r = -r_cgs / tcool

  ! Limit inflow velocity to subsonic
  v_r = max(v_r, -0.5d0 * sqrt(cs2))

  ! Azimuthal velocity from angular momentum conservation
  ! At large r: v_phi = v_c * (R_c,max / r) * sin(theta)
  ! This gives j = v_phi * R_cyl = v_c * R_c,max * sin^2(theta)
  if (r_cgs > Rc_max_cm) then
     v_phi = v_c_cgs * Rc_max_cm / r_cgs * sin(theta)
  else
     ! Inside circularization radius: Keplerian
     v_phi = v_c_cgs * sin(theta)
  endif

  ! Apply rotation corrections to density and temperature (Eq. 41)
  ! Only significant when r ~ R_c,max
  eps = (Rc_max_cm / r_cgs)**2
  if (eps < 1.0d0) then
     ! Density increases toward midplane
     rho = rho * (1.0d0 + eps * (11.0d0/4.0d0 * sin(theta)**2 - 35.0d0/24.0d0))
     ! Temperature decreases toward midplane
     T = T * (1.0d0 - eps * (2.0d0 * sin(theta)**2 - 5.0d0/6.0d0))
     ! Recalculate pressure for consistency
     P = rho / (mu * mH) * kB * T
  endif

contains
  real(dp) function kpc2cm_local()
    kpc2cm_local = 3.0856776d21
  end function kpc2cm_local

end subroutine get_cgm_state
