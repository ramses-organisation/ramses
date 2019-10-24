!================================================================
!================================================================
!================================================================
!================================================================
subroutine condinit(x,u,dx,nn)
  use amr_parameters
  use hydro_parameters
  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar+3)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:4): d.u,d.v,d.w, U(i,5): E, U(i,6:8): Bleft,
  ! U(i,nvar+1:nvar+3): Bright
  ! Q is the primitive variable vector. Conventions are here:
  ! Q(i,1): d, Q(i,2:4):u,v,w, Q(i,5): P, Q(i,6:8): Bleft,
  ! Q(i,nvar+1:nvar+3): Bright
  ! If nvar > 8, remaining variables (9:nvar) are treated as passive
  ! scalars in the hydro solver.
  ! U(:,:) and Q(:,:) are in user units.
  !================================================================
  integer::ivar,i
  real(dp),dimension(1:nvector,1:nvar+3),save::q   ! Primitive variables
  real(dp)::xc,xr,xl,yl,yr,yc,al,ar,B0,pi

  pi=ACOS(-1.0_dp)
  B0=1.0_dp/sqrt(4.0_dp*pi)
  do i=1,nn

     xl=x(i,1)-0.5_dp*dx
     xr=x(i,1)+0.5_dp*dx
     xc=x(i,1)
     yl=x(i,2)-0.5_dp*dx
     yr=x(i,2)+0.5_dp*dx
     yc=x(i,2)

     q(i,1)=25.0_dp/(36.0_dp*pi)
     q(i,2)=-sin(2.0_dp*pi*yc)
     q(i,3)=+sin(2.0_dp*pi*xc)
     q(i,4)=0.0_dp
     q(i,5)=5.0_dp/(12.0_dp*pi)

     Ar = B0*(cos(4.0_dp*pi*xl)/(4.0_dp*pi)+cos(2.0_dp*pi*yr)/(2.0_dp*pi))
     Al = B0*(cos(4.0_dp*pi*xl)/(4.0_dp*pi)+cos(2.0_dp*pi*yl)/(2.0_dp*pi))
     q(i,6)=(Ar-Al)/dx
     Ar = B0*(cos(4.0_dp*pi*xr)/(4.0_dp*pi)+cos(2.0_dp*pi*yr)/(2.0_dp*pi))
     Al = B0*(cos(4.0_dp*pi*xr)/(4.0_dp*pi)+cos(2.0_dp*pi*yl)/(2.0_dp*pi))
     q(i,nvar+1)=(Ar-Al)/dx
     Ar = B0*(cos(4.0_dp*pi*xr)/(4.0_dp*pi)+cos(2.0_dp*pi*yl)/(2.0_dp*pi))
     Al = B0*(cos(4.0_dp*pi*xl)/(4.0_dp*pi)+cos(2.0_dp*pi*yl)/(2.0_dp*pi))
     q(i,7)=(Al-Ar)/dx
     Ar = B0*(cos(4.0_dp*pi*xr)/(4.0_dp*pi)+cos(2.0_dp*pi*yr)/(2.0_dp*pi))
     Al = B0*(cos(4.0_dp*pi*xl)/(4.0_dp*pi)+cos(2.0_dp*pi*yr)/(2.0_dp*pi))
     q(i,nvar+2)=(Al-Ar)/dx

     q(i,8)=0.0_dp
     q(i,nvar+3)=0.0_dp
  end do


  ! Convert primitive to conservative variables
  ! density -> density
  u(1:nn,1)=q(1:nn,1)
  ! velocity -> momentum
  u(1:nn,2)=q(1:nn,1)*q(1:nn,2)
  u(1:nn,3)=q(1:nn,1)*q(1:nn,3)
  u(1:nn,4)=q(1:nn,1)*q(1:nn,4)
  ! kinetic energy
  u(1:nn,5)=0.0_dp
  u(1:nn,5)=u(1:nn,5)+0.5_dp*q(1:nn,1)*q(1:nn,2)**2
  u(1:nn,5)=u(1:nn,5)+0.5_dp*q(1:nn,1)*q(1:nn,3)**2
  u(1:nn,5)=u(1:nn,5)+0.5_dp*q(1:nn,1)*q(1:nn,4)**2
  ! pressure -> total fluid energy
  u(1:nn,5)=u(1:nn,5)+q(1:nn,5)/(gamma-1.0_dp)
  ! magnetic energy -> total fluid energy
  u(1:nn,5)=u(1:nn,5)+0.125_dp*(q(1:nn,6)+q(1:nn,nvar+1))**2
  u(1:nn,5)=u(1:nn,5)+0.125_dp*(q(1:nn,7)+q(1:nn,nvar+2))**2
  u(1:nn,5)=u(1:nn,5)+0.125_dp*(q(1:nn,8)+q(1:nn,nvar+3))**2
  u(1:nn,6:8)=q(1:nn,6:8)
  u(1:nn,nvar+1:nvar+3)=q(1:nn,nvar+1:nvar+3)
  ! passive scalars
  do ivar=9,nvar
     u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)
  end do

end subroutine condinit
!================================================================
!================================================================
!================================================================
!================================================================
subroutine velana(x,v,dx,t,ncell)
  use amr_parameters
  use hydro_parameters
  implicit none
  integer ::ncell                         ! Size of input arrays
  real(dp)::dx                            ! Cell size
  real(dp)::t                             ! Current time
  real(dp),dimension(1:nvector,1:3)::v    ! Velocity field
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine computes the user defined velocity fields.
  ! x(i,1:ndim) are cell center position in [0,boxlen] (user units).
  ! v(i,1:3) is the imposed 3-velocity in user units.
  !================================================================
  integer::i
  real(dp)::xx,yy,zz,vx,vy,vz,rr,tt,omega,aa,twopi

  ! Add here, if you wish, some user-defined initial conditions
  aa=1.0_dp
  twopi=2.0_dp*ACOS(-1.0_dp)
  do i=1,ncell

     xx=x(i,1)
     yy=x(i,2)
#if NDIM>2
     zz=x(i,3)
#endif

     ! ABC
     vx=aa*(cos(twopi*yy)+sin(twopi*zz))
     vy=aa*(sin(twopi*xx)+cos(twopi*zz))
     vz=aa*(cos(twopi*xx)+sin(twopi*yy))

     v(i,1)=vx
     v(i,2)=vy
     v(i,3)=vz

  end do


end subroutine velana
