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
  real(dp)::xc,xr,xl,yl,yr,yc,al,ar,r0,a0,twopi,rr,ss,tt

  ! density
  q(1:nn,1)=1.0

  ! velocity
  q(1:nn,2:4)=0.0
  R0=1.0
  do i=1,nn
     xc=x(i,1)-boxlen/2.0
     yc=x(i,2)-boxlen/2.0
     rr = sqrt(xc**2+yc**2)
     if(rr < R0)q(i,4)=1.0
  end do

  ! pressure
  q(1:nn,5)=1.0

  ! magnetic field
  call mag_screw(x,q,dx,nn)

  ! Convert primitive to conservative variables
  ! density -> density
  u(1:nn,1)=q(1:nn,1)
  ! velocity -> momentum
  u(1:nn,2)=q(1:nn,1)*q(1:nn,2)
  u(1:nn,3)=q(1:nn,1)*q(1:nn,3)
  u(1:nn,4)=q(1:nn,1)*q(1:nn,4)
  ! kinetic energy
  u(1:nn,5)=0.0d0
  u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,2)**2
  u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,3)**2
  u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,4)**2
  ! pressure -> total fluid energy
  u(1:nn,5)=u(1:nn,5)+q(1:nn,5)/(gamma-1.0d0)
  ! magnetic energy -> total fluid energy
  u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,6)+q(1:nn,nvar+1))**2
  u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,7)+q(1:nn,nvar+2))**2
  u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,8)+q(1:nn,nvar+3))**2
  u(1:nn,6:8)=q(1:nn,6:8)
  u(1:nn,nvar+1:nvar+3)=q(1:nn,nvar+1:nvar+3)
  ! passive scalars
#if NDIM > 8
  do ivar=9,nvar
     u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)
  end do
#endif

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
  aa=1.0
  twopi=2d0*ACOS(-1d0)
  do i=1,ncell

     xx=x(i,1)
#if NDIM > 1
     yy=x(i,2)
#endif
#if NDIM > 2
     zz=x(i,3)
#endif
     ! Ponomarenko
     xx=xx-boxlen/2.0
     yy=yy-boxlen/2.0
     rr=sqrt(xx**2+yy**2)
     if(yy>0)then
        tt=acos(xx/rr)
     else
        tt=-acos(xx/rr)+twopi
     endif
     if(rr<1.0)then
        omega=0.609711
        vz=0.792624
     else
        omega=0.0
        vz=0.0
     endif
     vx=-sin(tt)*rr*omega
     vy=+cos(tt)*rr*omega

     v(i,1)=vx
#if NDIM > 1
     v(i,2)=vy
#endif
#if NDIM > 2
     v(i,3)=vz
#endif
  end do


end subroutine velana
!================================================================
!================================================================
!================================================================
!================================================================
subroutine mag_screw(x,q,dx,nn)
  use amr_parameters
  use hydro_parameters, only: nvar
  use const
  implicit none

  real(dp),dimension(1:nvector,1:ndim)::x    ! Cell center position
  real(dp),dimension(1:nvector,1:nvar+3)::q  ! Primitive variables
  integer::nn                                ! Number of cells
  real(dp)::dx                               ! Cell size
  real(dp)::xx,yy,zz,A_screw
  integer::i,it,nticks
  real(dp)::dxmin,zfl,zceil,rtick,x_edge,y_edge
  real(dp),dimension(1:2,1:2,1:2)::A_mag     ! x/y,left/right,up/down
  real(dp)::dmax,dnfw

  ! Ponomarenko
  dxmin=boxlen*0.5d0**nlevelmax
  nticks=nint(dx/dxmin)

  do i=1,nn
    ! box-centered coordinates
    xx=x(i,1) - boxlen * 0.5d0
    yy=x(i,2) - boxlen * 0.5d0
    zz=x(i,3) - boxlen * 0.5d0

    zfl   = zz - 0.5d0*dx
    zceil = zz + 0.5d0*dx

    A_mag = 0.0

    ! A_(x,l)
    x_edge = xx + 0.5*(dxmin-dx)
    y_edge = yy - 0.5*dx
    do it=1,nticks
      ! floor
      A_mag(1,1,1) = A_mag(1,1,1) + A_screw(x_edge,y_edge,zfl,1)
      ! ceiling
      A_mag(1,1,2) = A_mag(1,1,2) + A_screw(x_edge,y_edge,zceil,1)
      ! advance along x
      x_edge=x_edge + dxmin
    end do

    ! A_(x,r)
    x_edge = xx + 0.5*(dxmin-dx)
    y_edge = yy + 0.5*dx
    do it=1,nticks
      ! floor
      A_mag(1,2,1) = A_mag(1,2,1) + A_screw(x_edge,y_edge,zfl,1)
      ! ceiling
      A_mag(1,2,2) = A_mag(1,2,2) + A_screw(x_edge,y_edge,zceil,1)
      ! advance along x
      x_edge=x_edge + dxmin
    end do

    ! A_(y,l)
    x_edge = xx - 0.5*dx
    y_edge = yy + 0.5*(dxmin-dx)
    do it=1,nticks
      ! floor
      A_mag(2,1,1) = A_mag(2,1,1) + A_screw(x_edge,y_edge,zfl,2)
      ! ceiling
      A_mag(2,1,2) = A_mag(2,1,2) + A_screw(x_edge,y_edge,zceil,2)
      ! advance along y
      y_edge=y_edge + dxmin
    end do

    ! A_(y,r)
    x_edge = xx + 0.5*dx
    y_edge = yy + 0.5*(dxmin-dx)
    do it=1,nticks
      ! floor
      A_mag(2,2,1) = A_mag(2,2,1) + A_screw(x_edge,y_edge,zfl,2)
      ! ceiling
      A_mag(2,2,2) = A_mag(2,2,2) + A_screw(x_edge,y_edge,zceil,2)
      ! advance along y
      y_edge=y_edge + dxmin
    end do

    ! average value
    A_mag = A_mag / DBLE(nticks)

    ! B left
    ! X direction
    q(i,6)= - (A_mag(2,1,2)-A_mag(2,1,1)) / dx ! [-d/dz A_y]
    ! Y direction
    q(i,7)=   (A_mag(1,1,2)-A_mag(1,1,1)) / dx ! [ d/dz A_x]
    ! Z direction  [d/dx A_y - d/dy A_x]
    q(i,8)=   ( A_mag(2,2,1)-A_mag(2,1,1) - A_mag(1,2,1)+A_mag(1,1,1) ) / dx

    ! B right
    ! X direction
    q(i,nvar+1)= - (A_mag(2,2,2)-A_mag(2,2,1)) / dx ! [-d/dz A_y]
    ! Y direction
    q(i,nvar+2)=   (A_mag(1,2,2)-A_mag(1,2,1)) / dx ! [ d/dz A_x]
    ! Z direction  [d/dx A_y - d/dy A_x]
    q(i,nvar+3)=   ( A_mag(2,2,2)-A_mag(2,1,2) - A_mag(1,2,2)+A_mag(1,1,2) ) / dx
  end do

end subroutine mag_screw

function A_screw(xx,yy,zz,dir)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none
  real(kind=8)::xx,yy,zz,A_screw
  real(kind=8)::A,rr,tt,twopi
  integer::dir

  twopi=2d0*acos(-1d0)

  rr=sqrt(xx**2+yy**2)
  if(yy>0)then
     tt=acos(xx/rr)
  else
     tt=-acos(xx/rr)+twopi
  endif
  A_screw=0d0
  if(rr<2.0)then
     A_screw=1d-1*(2d0-rr)*rr*cos(tt)*sin(twopi*(zz+boxlen/2)/boxlen)
     if(dir==1)A_screw=A_screw*cos(tt)
     if(dir==2)A_screw=A_screw*sin(tt)
  endif

  return
end function A_screw
