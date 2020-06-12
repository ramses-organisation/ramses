!############################################################
!############################################################
!############################################################
!############################################################
subroutine boundana(x,u,dx,ibound,nn)
  use amr_parameters
  use hydro_parameters
  use poisson_parameters
  implicit none
  integer ::ibound                        ! Index of boundary region
  integer ::nn                            ! Number of active cells
  real(dp)::dx                            ! Cell size
#ifdef SOLVERmhd
  real(dp),dimension(1:nvector,1:nvar+3)::u ! Conservative variables
#else
  real(dp),dimension(1:nvector,1:nvar)::u ! Conservative variables
#endif
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates boundary conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:ndim+1): d.u,d.v,d.w and U(i,ndim+2): E.
  ! If MHD, then:
  ! U(i,1): d, U(i,2:4): d.u,d.v,d.w, U(i,5): E,
  ! U(i,6:8): Bleft, U(i,nvar+1:nvar+3): Bright
  ! U is in user units.
  ! ibound is the index of the boundary region defined in the namelist.
  !================================================================
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! OLD version BOUNDANA
!   integer::ivar,i
!
! #ifdef SOLVERmhd
!   do ivar=1,nvar+3
! #else
!   do ivar=1,nvar
! #endif
!      do i=1,ncell
!         u(i,ivar)=boundary_var(ibound,ivar)
!      end do
!   end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer::i,j,id,iu,iv,iw,ip
  real(dp)::x1,x2,g,T_tmp
  real(dp)::rho1,rho2,rho3
  real(dp)::p1,p2,p3,pfactor
  real(dp)::gammainit1,gammainit2,gammainit3
  integer::ivar
  real(dp),dimension(1:nvector,1:nvar),save::q   ! Primitive variables

  id=1; iu=2; iv=3; iw=4; ip=ndim+2
  x1=x_center(1)
  x2=x_center(2)
  
  rho1=d_region(1)
  rho2=d_region(2)
  rho3=d_region(3)
  
  p1=p_region(1)
  p2=p_region(2)
  p3=p_region(3)

  gammainit1=gamma_region(1)
  gammainit2=gamma_region(2)
  gammainit3=gamma_region(3)
  ! gammainit1=gamma
  ! gammainit2=gamma
  ! gammainit3=gamma
  
  g = abs(gravity_params(1))

  do i=1,nn
    if(x(i,1) .lt. x1)then
      T_tmp = 1 - ((gammainit1-1.0)/gammainit1)*g*(rho1/p1)*(x(i,1))
      q(i,id)=rho1*((T_tmp)**(1.0/(gammainit1-1.0)))
      q(i,ip)=p1*((T_tmp)**(gammainit1/(gammainit1-1.0)))
    else if ((x(i,1) .ge. x1) .and. (x(i,1) .lt. x2)) then
      T_tmp = 1 - ((gammainit2-1.0)/gammainit2)*g*(rho2/p2)*((x(i,1)-x1))
      q(i,id)=rho2*((T_tmp)**(1.0/(gammainit2-1.0)))
      q(i,ip)=p2*((T_tmp)**(gammainit2/(gammainit2-1.0)))
    else 
      T_tmp = 1 - ((gammainit3-1.0)/gammainit3)*g*(rho3/p3)*((x(i,1)-x2))
      q(i,id)=rho3*((T_tmp)**(1.0/(gammainit3-1.0)))
      q(i,ip)=p3*((T_tmp)**(gammainit3/(gammainit3-1.0)))
    endif
    q(i,iu)=0.0d0
    if(ndim>1)q(i,iv)=0.0d0
    if(ndim>2)q(i,iw)=0.0d0
  end do

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
  u(1:nn,ndim+2)=0.0d0
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,2)**2
#if NDIM>1
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,3)**2
#endif
#if NDIM>2
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,4)**2
#endif

  ! pressure -> total fluid energy
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+q(1:nn,ndim+2)/(gamma-1.0d0)

  ! passive scalars
  do ivar=ndim+3,nvar
     u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)
  end do
  
end subroutine boundana
