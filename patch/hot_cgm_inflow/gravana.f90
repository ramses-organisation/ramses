!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine gravana(x,f,dx,ncell)
  use amr_parameters
  use poisson_parameters
  use cgm_commons
  use constants

  implicit none
  integer ::ncell                         ! Size of input arrays
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:ndim)::f ! Gravitational acceleration
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine computes the gravitational acceleration for the
  ! hot rotating CGM inflow problem.
  !
  ! Isothermal potential: Phi = v_c^2 * ln(r)
  ! Acceleration: g = -v_c^2 / r * r_hat
  !
  ! gravity_type = 4: Isothermal halo (v_c = const)
  !   gravity_params(1) = v_c [km/s]
  !   gravity_params(2) = softening radius [kpc]
  !   gravity_params(3:5) = center position [code units]
  !================================================================
  integer::i
  real(dp)::rx,ry,rz,rr,rr2
  real(dp)::v_c_cgs,v_c_code,r_soft,r_soft_code
  real(dp)::x_center,y_center,z_center
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp)::gmass,emass,xmass,ymass,zmass
  real(dp)::a1,a2,z0,f_max

  ! Get unit conversions
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Initialize to zero
  f(1:ncell,1:ndim) = 0.0d0

  !================================================================
  ! Standard RAMSES gravity types (1-3)
  !================================================================

  ! Constant vector
  if(gravity_type==1)then
     do i=1,ncell
        f(i,1) = gravity_params(1)
#if NDIM>1
        f(i,2) = gravity_params(2)
#endif
#if NDIM>2
        f(i,3) = gravity_params(3)
#endif
     end do
  end if

  ! Point mass
  if(gravity_type==2)then
     gmass=gravity_params(1) ! GM
     emass=gravity_params(2) ! Softening length
     xmass=gravity_params(3) ! Point mass coordinates
     ymass=gravity_params(4)
     zmass=gravity_params(5)
     do i=1,ncell
        rx=0.0d0; ry=0.0d0; rz=0.0d0
        rx=x(i,1)-xmass
#if NDIM>1
        ry=x(i,2)-ymass
#endif
#if NDIM>2
        rz=x(i,3)-zmass
#endif
        rr=sqrt(rx**2+ry**2+rz**2+emass**2)
        f(i,1)=-gmass*rx/rr**3
#if NDIM>1
        f(i,2)=-gmass*ry/rr**3
#endif
#if NDIM>2
        f(i,3)=-gmass*rz/rr**3
#endif
     end do
  end if

  ! Vertical galactic field (Kuijken & Gilmore 1989)
  if(gravity_type==3)then
     a1 = gravity_params(1)
     a2 = gravity_params(2)
     z0 = gravity_params(3)
     a1 = a1 * kpc2cm / Myr2sec**2 / scale_l * scale_t**2
     a2 = a2 / Myr2sec**2 * scale_t**2
     z0 = z0 * pc2cm / scale_l
     do i=1,ncell
        rz=x(i,ndim)-0.5d0*boxlen
        f(i,ndim)=-a1*rz/(rz**2+z0**2)**0.5 - a2*rz
     end do
  end if

  !================================================================
  ! Isothermal halo for CGM simulations (gravity_type = 4)
  ! g = -v_c^2 / r * r_hat
  !================================================================
  if(gravity_type==4)then
     ! Parameters
     v_c_cgs = gravity_params(1) * 1.0d5  ! km/s -> cm/s
     r_soft = gravity_params(2) * kpc2cm  ! kpc -> cm

     ! Center position (code units, default to box center)
     x_center = gravity_params(3)
     y_center = gravity_params(4)
     z_center = gravity_params(5)

     ! If center not specified, use box center
     if (x_center <= 0.0d0) x_center = boxlen / 2.0d0
     if (y_center <= 0.0d0) y_center = boxlen / 2.0d0
     if (z_center <= 0.0d0) z_center = boxlen / 2.0d0

     ! Convert v_c to code units: v_code = v_cgs / scale_v
     v_c_code = v_c_cgs / scale_v

     ! Softening in code units
     r_soft_code = r_soft / scale_l

     do i=1,ncell
        ! Position relative to center
        rx = x(i,1) - x_center
#if NDIM>1
        ry = x(i,2) - y_center
#endif
#if NDIM>2
        rz = x(i,3) - z_center
#endif

        ! Radius with softening
        rr2 = rx**2 + r_soft_code**2
#if NDIM>1
        rr2 = rr2 + ry**2
#endif
#if NDIM>2
        rr2 = rr2 + rz**2
#endif
        rr = sqrt(rr2)

        ! Gravitational acceleration: g = -v_c^2 / r * r_hat
        ! In code units where G=1, this comes from Phi = v_c^2 * ln(r)
        f(i,1) = -v_c_code**2 * rx / rr2
#if NDIM>1
        f(i,2) = -v_c_code**2 * ry / rr2
#endif
#if NDIM>2
        f(i,3) = -v_c_code**2 * rz / rr2
#endif
     end do
  end if

end subroutine gravana

!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine phi_ana(rr,pp,ngrid)
  use amr_commons
  use poisson_commons
  use constants, only: twopi
  implicit none
  integer::ngrid
  real(dp),dimension(1:nvector)::rr,pp
  ! -------------------------------------------------------------------
  ! Analytical potential for boundary conditions
  ! -------------------------------------------------------------------

  integer :: i
  real(dp):: fourpi

  fourpi=2*twopi

  do i=1,ngrid
#if NDIM==1
     pp(i)=multipole(1)*fourpi/2*rr(i)
#elif NDIM==2
     pp(i)=multipole(1)*2*log(rr(i))
#elif NDIM==3
     pp(i)=-multipole(1)/rr(i)
#endif
  end do
end subroutine phi_ana
