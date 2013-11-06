!================================================================
!================================================================
!================================================================
!================================================================
subroutine condinit(x,u,dx,nn)
  use amr_parameters
  use hydro_parameters
  use poisson_parameters
  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:ndim+1): d.u,d.v,d.w and U(i,ndim+2): E.
  ! Q is the primitive variable vector. Conventions are here:
  ! Q(i,1): d, Q(i,2:ndim+1):u,v,w and Q(i,ndim+2): P.
  ! If nvar >= ndim+3, remaining variables are treated as passive
  ! scalars in the hydro solver.
  ! U(:,:) and Q(:,:) are in user units.
  !================================================================
  integer::i,j,id,iu,iv,iw,ip
  real(dp)::rho1,rho2,v1,v2,p1,p2,u1,u2,w1,w2
  integer::ivar
  real(dp),dimension(1:nvector,1:nvar),save::q   ! Primitive variables
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2,pi
  real(dp)::ramp_width,ramp,f
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  pi=3.14159265d0

  id=1; iu=2; iv=3; iw=4; ip=ndim+2

  rho1=d_region(1)/scale_d
  rho2=d_region(2)/scale_d
  v1=v_region(1)/scale_v
  v2=v_region(2)/scale_v
  u1=u_region(1)/scale_v
  u2=u_region(2)/scale_v
  w1=w_region(1)/scale_v
  w2=w_region(2)/scale_v
  p1=p_region(1)*scale_l**2*scale_d/(scale_t**2)
  p2=p_region(2)*scale_l**2*scale_d/(scale_t**2)
  
  ramp_width=0.2/scale_l
!  write(*,*) rho1,p1,rho1*p1/scale_T2

!  write(*,*) rho1,rho2,v1,v2,p1,p2,boxlen
  do i=1,nn
     if(x(i,2)>boxlen/2.)then 
        q(i,id)=rho2
        !q(i,iu)=u1*cos(2.*pi*(10.*x(i,1)))*exp(-abs(x(i,2)-boxlen/2.0)/(boxlen/10.0)) 

#if NDIM>2  
        q(i,iw)=0.0
#endif
        q(i,ip)=p1
     else 
        q(i,id)=rho1
        !q(i,iu)=u2*cos(2.*pi*(10.*x(i,1)))*exp(-abs(x(i,2)-boxlen/2.0)/(boxlen/10.0)) 

#if NDIM>2  
        q(i,iw)=0.0
#endif
        q(i,ip)=p2
     endif
     ramp=1./(1.+exp(-2/ramp_width*(x(i,2)-boxlen/2.0)))
     f=cos(2*pi*10.*x(i,1))*cos(2*pi*3.*x(i,1))*exp(-10.*abs(x(i,2)-0.5))
     q(i,iu)=f*(u1+ramp*(u2-u1))
     q(i,iv)=v1+ramp*(v2-v1)
  enddo
 
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

end subroutine condinit
