!############################################################
!############################################################
!############################################################
!############################################################
subroutine boundana(x,u,dx,ibound,ncell)
  use amr_parameters, ONLY: dp,ndim,nvector
  use hydro_parameters, ONLY: nvar,boundary_var,gamma
  implicit none
  integer ::ibound                        ! Index of boundary region
  integer ::ncell                         ! Number of active cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates boundary conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:ndim+1): d.u,d.v,d.w and U(i,ndim+2): E.
  ! U is in user units.
  ! ibound is the index of the boundary region defined in the namelist.
  !================================================================
  integer::ivar,i
  real(dp)::r,rmin,rmax,vx,vy,vz,lor,d,h,p,boxlen
  real(dp)::x0,y0

  ! Add here, if you wish, some user-defined boudary conditions
  !Boundary conditions for a 3D jet coming from z=0
  ! the box is 0<z<20 by -8<x,y<8. Initially the jet
  ! is restricted to z<1 ,r<1 with inflow boundary conditions.
  !The rest of the left boundary is zero gradient
  x0=0.
  y0=0.

  do i=1,ncell! ou nn?
     r=sqrt((x(i,1)-x0)**2+(x(i,2)-y0)**2)
     if (r<1.) then
           !inflowing jet
           vx=0.99d0
           vy=0.0d0
           vz=0.0d0
           d =.1d0
           p =0.01d0
           h=1.0d0+gamma/(gamma-1.0d0)*p/d
           lor=(1.-(vx**2+vy**2+vz**2))**(-1d0/2d0)

           U(i,1) = d*lor
           U(i,2) = lor**2*d*h*vx
           U(i,3) = lor**2*d*h*vy
           U(i,4) = lor**2*d*h*vz
           U(i,5) = lor**2*d*h-p

        endif
  end do










end subroutine boundana
