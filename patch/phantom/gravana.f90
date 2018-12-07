!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine gravana(x,force,dx,ncell)
  use amr_parameters
  use poisson_parameters
  use poisson_commons
  use mond_commons
  implicit none
  integer ::ncell                         ! Size of input arrays
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:ndim)::force ! Gravitational acceleration
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine computes the acceleration using analytical models.
  ! x(i,1:ndim) are cell center position in [0,boxlen] (user units).
  ! f(i,1:ndim) is the gravitational acceleration in user units.
  !================================================================
  integer::idim,i
  real(dp)::gmass,emass,xmass,ymass,zmass,rr,rx,ry,rz,constant

  real(dp),dimension(1:nvector)::r ! Cell center position.

  ! Self-gravity
  ! TODO: Check boundary conditions
  if (gravity_type==0) then

     ! Compute position relative to center of mass
     r = 0.0
     do idim=1,ndim
        do i=1,ncell
           x(i,idim) = x(i,idim) - multipole(idim+1)/multipole(1)
           r(i) = r(i) + x(i,idim)**2
        end do
     enddo
     do i=1,ncell
         r(i)=sqrt(r(i))
     end do

     if (connected_Mond) then
        ! Compute MONDian acceleration a_i = -sqrt(G*M*a0)*x_i/r^2  (a ~ 1/r)
        constant = -sqrt(multipole(1)*a0)
        do idim=1,ndim
           do i=1,ncell
              force(i,idim) = const * x(i,idim) / r(i)**2
           end do
        enddo
     else
        ! Compute Newtonian acceleration a_i = -G*M*x_i/r^3  (a ~ 1/r^2)
        constant = -multipole(1)
        do idim=1,ndim
           do i=1,ncell
              force(i,idim) = constant * x(i,idim) / r(i)**3
           end do
        enddo
     endif

  endif

  ! Constant vector
  if(gravity_type==1)then
     do idim=1,ndim
        do i=1,ncell
           force(i,idim)=gravity_params(idim)
        end do
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
        force(i,1)=-gmass*rx/rr**3
#if NDIM>1
        force(i,2)=-gmass*ry/rr**3
#endif
#if NDIM>2
        force(i,3)=-gmass*rz/rr**3
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
  use mond_commons
  implicit none
  integer::ngrid
  real(dp),dimension(1:nvector)::rr,pp
  ! -------------------------------------------------------------------
  ! This routine set up boundary conditions for fine levels.
  ! -------------------------------------------------------------------

  integer :: i
  real(dp):: fourpi
  real(dp)::constant

#if NDIM==1
  fourpi=4d0*ACOS(-1.0D0)
  constant = multipole(1)*fourpi/2d0
  do i=1,ngrid
     pp(i)=constant*rr(i)
  end do
#endif
#if NDIM==2
  constant = multipole(1)*2d0
  do i=1,ngrid
     pp(i)=constant*log(rr(i))
  end do
#endif
#if NDIM==3
  if (connected_Mond) then
     ! Compute the MONDian potential, phi = sqrt(G M a0) log(r)
     ! This approximation requires:
     ! (a) being in the deep MOND-limit
     ! (b) r being large enough [this is also the case for Newtonian dynamics]
     constant = sqrt(multipole(1)*a0)
     do i=1,ngrid
        pp(i)=constant*log(rr(i))
     end do
  else
     ! Compute the Newtonian potential, phi = -G M/r
     constant = -multipole(1)
     do i=1,ngrid
        pp(i)=constant/rr(i)
     end do
  endif
#endif
end subroutine phi_ana
