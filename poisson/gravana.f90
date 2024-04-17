!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine gravana(x,f,dx,ncell)
  use amr_parameters
  use poisson_parameters
  use poisson_commons, only: multipole
  use constants

  implicit none
  integer ::ncell                         ! Size of input arrays
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:ndim)::f ! Gravitational acceleration
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine computes the acceleration using analytical models.
  ! x(i,1:ndim) are cell center position in [0,boxlen] (user units).
  ! f(i,1:ndim) is the gravitational acceleration in user units.
  ! Only if there is no self-gravity
  !================================================================
  integer::idim,i
  real(dp)::gmass,emass,xmass,ymass,zmass,rr,rx,ry,rz
  real(dp):: a1,a2,z0,sigma,f_max
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2

  ! Constant vector
  if(gravity_type==1)then
     do idim=1,ndim
        do i=1,ncell
           f(i,idim)=gravity_params(idim)
        end do
     end do
  end if

  ! Point mass
  if(gravity_type==2)then
     gmass=gravity_params(1) ! GM
     emass=dx
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

  if(gravity_type==3)then
     ! vertical galactic gravitational field
     ! Kuijken & Gilmore 1989 taken from Joung & MacLow (2006)
     ! g = -a1 z / sqrt(z^2+z0^2) - a2 z
     a1 = gravity_params(1) ! Star potential coefficient in kpc Myr-2
     a2 = gravity_params(2) ! DM potential coefficient in Myr-2
     z0 = gravity_params(3) ! Scale height in pc
    ! standard values are: a1 = 1.42d-3, a2 = 5.49d-4, z0 = 0.18d3 pc

    ! The gravitational field is given by
    ! g = -a1 z / sqrt(z^2+z0^2) - a2 z
    ! rho = [(a1 / z0) ( (z/z0)^2 + 1)^(-3/2) + a2] / (4piG)

    ! convert to code units
     call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
     a1 = a1 * kpc2cm / Myr2sec**2 / scale_l * scale_t**2
     a2 = a2 / Myr2sec**2 * scale_t**2
     z0 = z0 * pc2cm / scale_l

     do i=1,ncell
        ! the last dimension is vertical (1D -> x, 2D -> y, 3D -> z)
        rz=x(i,ndim)-0.5d0*boxlen
        f(i,ndim)=-a1*rz/(rz**2+z0**2)**0.5 - a2*rz
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
  ! This routine set up boundary conditions for fine levels.
  ! -------------------------------------------------------------------

  integer :: i
  real(dp):: fourpi

  fourpi=2*twopi

#if NDIM==1
  do i=1,ngrid
     pp(i)=multipole(1)*fourpi/2*rr(i)
  end do
#endif
#if NDIM==2
  do i=1,ngrid
     pp(i)=multipole(1)*2*log(rr(i))
  end do
#endif
#if NDIM==3
  do i=1,ngrid
     pp(i)=-multipole(1)/rr(i)
  end do
#endif
end subroutine phi_ana
