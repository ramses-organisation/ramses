! ---------------------------------------------------------------
!  UNSPLIT     Unsplit second order Godunov integrator for
!              polytropic relativistic gas dynamics using MUSCL-HANCOCK  with
!              hll and hllc riemann solver
!              with various slope limiters (moncen crashes for the moment).
!              This routine is inspired by the mhd/umuscl scheme
!
!
!  inputs/outputs
!  uin         => (const)  input state
!  gravin      => (const)  input gravitational acceleration
!  iu1,iu2     => (const)  first and last index of input array,
!  ju1,ju2     => (const)  cell centered,
!  ku1,ku2     => (const)  including buffer cells.
!  flux       <=  (modify) return fluxes in the 3 coord directions
!  if1,if2     => (const)  first and last index of output array,
!  jf1,jf2     => (const)  edge centered,
!  kf1,kf2     => (const)  for active cells only.
!  dx,dy,dz    => (const)  (dx,dy,dz)
!  dt          => (const)  time step
!  ngrid       => (const)  number of sub-grids
!  ndim        => (const)  number of dimensions
!
!  uin = (\rho, \rho u, \rho v, \rho w, Etot)
!  Note that here we have 3 components for v  whatever ndim.
!
!  This routine was written by Astrid Lamberts and Sebastien Fromang
! ----------------------------------------------------------------
!
! ----------------------------------------------------------------
subroutine unsplit(uin,gravin,flux,tmp,dx,dy,dz,dt,ngrid)
  use amr_parameters
  use const
  use hydro_parameters
  implicit none

  integer ::ngrid
  real(dp)::dx,dy,dz,dt

  ! Input states
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::uin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim)::gravin

  ! Output fluxes
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:nvar,1:ndim)::flux
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:2   ,1:ndim)::tmp

  ! Primitive variables
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar),save::qin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2       ),save::cin

  ! Slopes
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim),save::dq

  ! Left and right state arrays
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim),save::qm
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim),save::qp

  ! Intermediate fluxes
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar),save::fx
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:2   ),save::tx

  ! Local scalar variables
  integer::i,j,k,l,ivar
  integer::ilo,ihi,jlo,jhi,klo,khi

  ilo=MIN(1,iu1+2); ihi=MAX(1,iu2-2)
  jlo=MIN(1,ju1+2); jhi=MAX(1,ju2-2)
  klo=MIN(1,ku1+2); khi=MAX(1,ku2-2)


  ! Translate to primative variables, compute sound speeds
  call ctoprim(uin,qin,  gravin,dt,ngrid)

  ! Compute TVD slopes
  call uslope(qin,dq,dx,dt,ngrid)

  ! Compute 3D traced-states in all three directions
  if(scheme=='muscl')then
#if NDIM==1
     call trace1d(qin,dq,qm,qp,dx      ,dt,ngrid)
#endif
#if NDIM==2
     call trace2d(qin,dq,qm,qp,dx,dy   ,dt,ngrid)
#endif
#if NDIM==3
     call trace3d(qin,dq,qm,qp,dx,dy,dz,dt,ngrid)
#endif

   else if (scheme=='plmde')then
!#if NDIM==1
     call tracex  (qin,dq,cin,qm,qp,dx      ,dt,ngrid)
!#endif
#if NDIM>1
     call tracexy (qin,dq,cin,qm,qp,dx,dy   ,dt,ngrid)
#endif
#if NDIM>2
     call tracexyz(qin,dq,cin,qm,qp,dx,dy,dz,dt,ngrid)
#endif
   else
      write(*,*),'no valid scheme'
      stop

   endif




  ! Solve for 1D flux in X direction
  call cmpflxm(qm,iu1+1,iu2+1,ju1  ,ju2  ,ku1  ,ku2  , &
       &       qp,iu1  ,iu2  ,ju1  ,ju2  ,ku1  ,ku2  , &
       &          if1  ,if2  ,jlo  ,jhi  ,klo  ,khi  , 2,3,4,fx,tx,ngrid)



  ! Save flux in output array
  do i=if1,if2
  do j=jlo,jhi
  do k=klo,khi
     do ivar=1,nvar
        do l=1,ngrid
           flux(l,i,j,k,ivar,1)=fx(l,i,j,k,ivar)*dt/dx
        end do
     end do
     do ivar=1,2
        do l=1,ngrid
           tmp (l,i,j,k,ivar,1)=tx(l,i,j,k,ivar)*dt/dx
        end do
     end do
  end do
  end do
  end do


#if NDIM>1
  ! Solve for 1D flux in Y direction
  call cmpflxm(qm,iu1  ,iu2  ,ju1+1,ju2+1,ku1  ,ku2  , &
       &       qp,iu1  ,iu2  ,ju1  ,ju2  ,ku1  ,ku2  , &
       &          ilo  ,ihi  ,jf1  ,jf2  ,klo  ,khi  , 3,2,4,fx,tx,ngrid)
  ! Save flux in output array
  do i=ilo,ihi
  do j=jf1,jf2
  do k=klo,khi
     do ivar=1,nvar
        do l=1,ngrid
           flux(l,i,j,k,ivar,2)=fx(l,i,j,k,ivar)*dt/dy
        end do
     end do
     do ivar=1,2
        do l=1,ngrid
           tmp (l,i,j,k,ivar,2)=tx(l,i,j,k,ivar)*dt/dy
        end do
     end do
  end do
  end do
  end do
#endif

  ! Solve for 1D flux in Z direction
#if NDIM>2
  call cmpflxm(qm,iu1  ,iu2  ,ju1  ,ju2  ,ku1+1,ku2+1, &
       &       qp,iu1  ,iu2  ,ju1  ,ju2  ,ku1  ,ku2  , &
       &          ilo  ,ihi  ,jlo  ,jhi  ,kf1  ,kf2  , 4,2,3,fx,tx,ngrid)
  ! Save flux in output array
  do i=ilo,ihi
  do j=jlo,jhi
  do k=kf1,kf2
     do ivar=1,nvar
        do l=1,ngrid
           flux(l,i,j,k,ivar,3)=fx(l,i,j,k,ivar)*dt/dz
        end do
     end do
     do ivar=1,2
        do l=1,ngrid
           tmp (l,i,j,k,ivar,3)=tx(l,i,j,k,ivar)*dt/dz
        end do
     end do
  end do
  end do
  end do
#endif

end subroutine unsplit
!###########################################################
!###########################################################
!###########################################################
!###########################################################
!  Subroutine TRACE1D
!
!> This subroutine assumed that one considers a 1D problem
!! with 3 velocities (vx,vy,vz)

subroutine trace1d(q,dq,qm,qp,dx,dt,ngrid)

  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer :: ngrid
  real(dp) :: dx,dt

  real(dp), dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar       ) :: q
  real(dp), dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim) :: dq
  real(dp), dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim) :: qm
  real(dp), dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim) :: qp

  ! Local variables
  integer :: i, j, k, l, n
  integer :: ilo,ihi,jlo,jhi,klo,khi
  integer :: ir,iu,iv,iw,ip
  real(dp)    :: dtdx,smallp
  real(dp)    :: r,u,v,w,p,a
  real(dp)    :: drx,dux,dvx,dwx,dpx,dax
  real(dp)    :: sr0,su0=0.,sv0=0.,sw0=0.,sp0=0.,sa0=0.
  real(dp)    :: h,kappa,chi,f,ni,cs2,velsq,entho,lorsq
  real(dp)    :: a11,a12,a15,a22,a25,a32,a33,a35,a42,a44,a45,a52,a55
  real(dp)    :: lor,tau,ka

  smallp = smallr*smallc**2/gamma

  ilo=min(1,iu1+1); ihi=max(1,iu2-1)
  jlo=min(1,ju1+1); jhi=max(1,ju2-1)
  klo=min(1,ku1+1); khi=max(1,ku2-1)

  ir = 1; iu = 2; iv = 3 ; iw = 4 ; ip = 5

  dtdx = dt/dx
  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           do l= 1,ngrid

              ! Cell centered values
              r = q(l,i,j,k,ir)
              u = q(l,i,j,k,iu)
              v = q(l,i,j,k,iv)
              w = q(l,i,j,k,iw)
              p = q(l,i,j,k,ip)

              velsq=u*u+v*v+w*w ; lor=(1.0d0-velsq)**(-1./2.) ; lorsq=lor*lor

              ! TVD slopes in X direction
              drx = half*dq(l,i,j,k,ir,1)
              dux = half*dq(l,i,j,k,iu,1)
              dvx = half*dq(l,i,j,k,iv,1)
              dwx = half*dq(l,i,j,k,iw,1)
              dpx = half*dq(l,i,j,k,ip,1)

              ! Source terms (including transverse derivatives)
              tau=p/r
              if (eos .eq. 'TM') then
                 ka=9d0/4d0*tau**2+1
                 h=5d0/2d0*tau+sqrt(ka)
                 chi=-5d0/2d0*tau/r-9d0/4d0/sqrt(ka)*tau**2/r
                 kappa=5d0/2d0/r+9d0/4d0/sqrt(ka)*p/r**2
              else
                 entho = gamma/(gamma-one)
                 h     = 1.0d0+entho*p/r   !enthalpy
                 kappa = entho*1d0/r      !dh/dp
                 chi   = -entho*p/r**2     !dh/dr
              endif
              f=lorsq*(-h+h*r*kappa+r*chi*velsq) !normalizing factor

              a11= u*f
              a12= h*r*lorsq*(r*kappa-one)
              a15= u*(-r*lorsq*kappa+r*lorsq*velsq*kappa+1)

              a22=u*lorsq*(-h+r*chi+r*h*kappa)
              a25= (r*lorsq*kappa*h-r*h*kappa*lorsq*u*u-h-h*lorsq*(v*v+w*w)+lorsq*r*chi*(v*v+w*w))/(r*lorsq*h)

              a32= v*r*chi
              a33= u*f
              a35= -(u*v*(-h+r*chi+r*h*kappa))/(r*h)

              a42= w*r*chi
              a44= u*f
              a45= -(u*w*(-h+r*chi+r*h*kappa))/(r*h)

              a52= -r*r*h*lorsq*chi
              a55= lorsq*u*(-h+r*chi+h*r*kappa)

              sr0 = -a11*drx-a12*dux                -a15*dpx
              su0 =         -a22*dux                -a25*dpx
              sv0 =         -a32*dux-a33*dvx        -a35*dpx
              sw0 =         -a42*dux        -a44*dwx-a45*dpx
              sp0 =         -a52*dux                -a55*dpx

              sr0 = sr0/f*dtdx
              su0 = su0/f*dtdx
              sv0 = sv0/f*dtdx
              sw0 = sw0/f*dtdx
              sp0 = sp0/f*dtdx

              ! Cell-centered predicted states
              r = r + sr0
              u = u + su0
              v = v + sv0
              w = w + sw0
              p = p + sp0

              ! Right state at left interface
              qp(l,i,j,k,ir,1) = r - drx
              qp(l,i,j,k,iu,1) = u - dux
              qp(l,i,j,k,iv,1) = v - dvx
              qp(l,i,j,k,iw,1) = w - dwx
              qp(l,i,j,k,ip,1) = p - dpx
             if (qp(l,i,j,k,ir,1)<0.) write (*,*) 'd negative',i,qp(l,i,j,k,ir,1)
             if (qp(l,i,j,k,ip,1)<0.) write (*,*) 'p negative',i,qp(l,i,j,k,ip,1)
              qp(l,i,j,k,ir,1) = max(smallr, qp(l,i,j,k,ir,1))
              qp(l,i,j,k,ip,1) = max(smallp, qp(l,i,j,k,ip,1))

              ! Left state at right interface
              qm(l,i,j,k,ir,1) = r + drx
              qm(l,i,j,k,iu,1) = u + dux
              qm(l,i,j,k,iv,1) = v + dvx
              qm(l,i,j,k,iw,1) = w + dwx
              qm(l,i,j,k,ip,1) = p + dpx
             if (qm(l,i,j,k,ir,1)<0.) write (*,*) 'd negative',i,qm(l,i,j,k,ir,1)
             if (qm(l,i,j,k,ip,1)<0.) write (*,*) 'p negative',i,qm(l,i,j,k,ip,1)
              qm(l,i,j,k,ir,1) = max(smallr, qm(l,i,j,k,ir,1))
              qm(l,i,j,k,ip,1) = max(smallp, qm(l,i,j,k,ip,1))

              velsq=qp(l,i,j,k,iu,1)*qp(l,i,j,k,iu,1)+qp(l,i,j,k,iv,1)*qp(l,i,j,k,iv,1)+qp(l,i,j,k,iw,1)*qp(l,i,j,k,iw,1)
              if (velsq>1d0) then
                 qp(l,i,j,k,ir:ip,1)=q(l,i,j,k,ir:ip)
                 qm(l,i,j,k,ir:ip,1)=q(l,i,j,k,ir:ip)
              endif
              velsq=qm(l,i,j,k,iu,1)**2+qm(l,i,j,k,iv,1)**2+qm(l,i,j,k,iw,1)**2
              if (velsq>1d0) then
                 qp(l,i,j,k,ir:ip,1)=q(l,i,j,k,ir:ip)
                 qm(l,i,j,k,ir:ip,1)=q(l,i,j,k,ir:ip)
              endif

              ! Passive scalars
              do n = 6, nvar
                 a   = q(l,i,j,k,n)       ! Cell centered values
                 u   = q(l,i,j,k,iu)
                 dax = dq(l,i,j,k,n,1)    ! TVD slopes
                 sa0 = -u*dax             ! Source terms
                 qp(l,i,j,k,n,1) = a - half*dax + sa0*dtdx*half*lor   ! Right state
                 qm(l,i,j,k,n,1) = a + half*dax + sa0*dtdx*half*lor   ! Left state
              end do
           end do
        end do
     end do
  end do


end subroutine trace1d



!###########################################################
!###########################################################
!###########################################################
!###########################################################
#if NDIM>1
subroutine trace2d(q,dq,qm,qp,dx,dy,dt,ngrid)

  use amr_parameters
  use hydro_parameters
  use const

  implicit none

  real(dp)    :: dx,dy,dt
  integer :: ngrid

  real(dp), dimension(1:nvector,iu1:iu2  ,ju1:ju2  ,ku1:ku2  ,1:nvar       ) :: q
  real(dp), dimension(1:nvector,iu1:iu2  ,ju1:ju2  ,ku1:ku2  ,1:nvar,1:ndim) :: dq
  real(dp), dimension(1:nvector,iu1:iu2  ,ju1:ju2  ,ku1:ku2  ,1:nvar,1:ndim) :: qm
  real(dp), dimension(1:nvector,iu1:iu2  ,ju1:ju2  ,ku1:ku2  ,1:nvar,1:ndim) :: qp

  ! Local variables
  integer :: i,j,k,l,n
  integer :: ilo,ihi,jlo,jhi,klo,khi
  integer :: ir,iu,iv,iw,ip
  real(dp)    :: dtdx,dtdy,xc,xL,xR,yL,yR,smallp
  real(dp)    :: sr0,su0,sv0,sw0,sp0,sa0
  real(dp)    :: r,u,v,w,p ,a
  real(dp)    :: drx,dux,dvx,dwx,dpx, dax
  real(dp)    :: dry,duy,dvy,dwy,dpy, day
  real(dp)    :: lor,lorsq,h,kappa,chi,f,ni,cs2,vtot,velsq,entho,ka,tau
  real(dp)    :: a11,a12,a15,a22,a25,a32,a33,a35,a42,a44,a45,a52,a55
  real(dp)    :: b11,b13,b15,b22,b23,b25,b33,b35,b43,b44,b45,b53,b55

  smallp = smallr*smallc**2/gamma

  ilo=min(1,iu1+1); ihi=max(1,iu2-1)
  jlo=min(1,ju1+1); jhi=max(1,ju2-1)
  klo=min(1,ku1+1); khi=max(1,ku2-1)

  ir=1; iu=2; iv=3; iw=4 ; ip=5


  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           do l = 1, ngrid

              dtdx = dt/dx
              dtdy = dt/dy

              ! Cell centered values
              r = q(l,i,j,k,ir)
              u = q(l,i,j,k,iu)
              v = q(l,i,j,k,iv)
              w = q(l,i,j,k,iw)
              p = q(l,i,j,k,ip)

              velsq=u**2+v**2+w**2 ; lorsq=1d0/(1d0-velsq) ; lor=sqrt(lorsq)
              vtot=sqrt(u**2+v**2+w**2)

              ! Cell centered TVD slopes in X direction
              drx = half * dq(l,i,j,k,ir,1)
              dux = half * dq(l,i,j,k,iu,1)
              dvx = half * dq(l,i,j,k,iv,1)
              dwx = half * dq(l,i,j,k,iw,1)
              dpx = half * dq(l,i,j,k,ip,1)

              ! Cell centered TVD slopes in Y direction
              dry = half * dq(l,i,j,k,ir,2)
              duy = half * dq(l,i,j,k,iu,2)
              dvy = half * dq(l,i,j,k,iv,2)
              dwy = half * dq(l,i,j,k,iw,2)
              dpy = half * dq(l,i,j,k,ip,2)


              ! Source terms (including transverse derivatives)

              if (eos .eq. 'TM') then
                 tau=p/r
                 ka=9d0/4d0*tau**2+1
                 h=5d0/2d0*tau+sqrt(ka)
                 chi=-5d0/2d0*tau/r-9d0/4d0/sqrt(ka)*tau**2/r
                 kappa=5d0/2d0/r+9d0/4d0/sqrt(ka)*p/r**2
              else
                 entho = gamma/(gamma-one)
                 h     = 1.0d0+entho*p/r   !enthalpy
                 kappa = entho*1d0/r      !dh/dp
                 chi   = -entho*p/r**2     !dh/dr
              endif



              f=lorsq*(-h+h*r*kappa+r*chi*velsq)

              a11= u*f
              a12= h*r*lorsq*(r*kappa-1)
              a15= u*(-r*lorsq*kappa+r*lorsq*velsq*kappa+1)

              a22= u*lorsq*(-h+r*chi+r*h*kappa)
              a25= (r*lorsq*kappa*h-r*h*kappa*lorsq*u**2-h-h*lorsq*(v**2+w**2)+lorsq*r*chi*(v**2+w**2))/(r*lorsq*h)

              a32= v*r*chi
              a33= u*f
              a35= -(u*v*(-h+r*chi+r*h*kappa))/(r*h)

              a42= w*r*chi
              a44= u*f
              a45= -(u*w*(-h+r*chi+r*h*kappa))/(r*h)

              a52= -r**2*h*lorsq*chi
              a55= lorsq*u*(-h+r*chi+h*r*kappa)


              b11= v*f
              b13= h*r*lorsq*(r*kappa-1)
              b15= v*(-r*lorsq*kappa+r*lorsq*velsq*kappa+1)

              b22= v*f
              b23= u*r*chi
              b25=-(v*u*(-h+r*chi+r*h*kappa))/(r*h)

              b33= v*lorsq*(-h+r*chi+r*h*kappa)
              b35= (r*lorsq*kappa*h-r*h*kappa*lorsq*v**2-h-h*lorsq*(u**2+w**2)+lorsq*r*chi*(u**2+w**2))/(r*lorsq*h)

              b43= w*r*chi
              b44= v*f
              b45=-(v*w*(-h+r*chi+r*h*kappa))/(r*h)

              b53=-r**2*h*lorsq*chi
              b55= lorsq*v*(-h+r*chi+h*r*kappa)


              ! Source terms (including transverse derivatives)
              sr0 = (-a11*drx-a12*dux                -a15*dpx)*dtdx/f + (-b11*dry        -b13*dvy        -b15*dpy)*dtdy/f
              su0 = (        -a22*dux                -a25*dpx)*dtdx/f + (        -b22*duy-b23*dvy        -b25*dpy)*dtdy/f
              sv0 = (        -a32*dux-a33*dvx        -a35*dpx)*dtdx/f + (                -b33*dvy        -b35*dpy)*dtdy/f
              sw0 = (        -a42*dux        -a44*dwx-a45*dpx)*dtdx/f + (                -b43*dvy-b44*dwy-b45*dpy)*dtdy/f
              sp0 = (        -a52*dux                -a55*dpx)*dtdx/f + (                -b53*dvy        -b55*dpy)*dtdy/f


              ! Cell-centered predicted states
              r = r + sr0
              u = u + su0
              v = v + sv0
              w = w + sw0
              p = p + sp0

              ! Right state at left interface
              qp(l,i,j,k,ir,1) = r - drx
              qp(l,i,j,k,iu,1) = u - dux
              qp(l,i,j,k,iv,1) = v - dvx
              qp(l,i,j,k,iw,1) = w - dwx
              qp(l,i,j,k,ip,1) = p - dpx
              if (qp(l,i,j,k,ir,1)<0.) write (*,*) 'd negative',i,qp(l,i,j,k,ir,1)
             if (qp(l,i,j,k,ip,1)<0.) write (*,*) 'p negative',i,qp(l,i,j,k,ip,1)

              qp(l,i,j,k,ir,1) = max(smallr, qp(l,i,j,k,ir,1))
              qp(l,i,j,k,ip,1) = max(smallp, qp(l,i,j,k,ip,1))

              ! Left state at right interface
              qm(l,i,j,k,ir,1) = r + drx
              qm(l,i,j,k,iu,1) = u + dux
              qm(l,i,j,k,iv,1) = v + dvx
              qm(l,i,j,k,iw,1) = w + dwx
              qm(l,i,j,k,ip,1) = p + dpx
              if (qm(l,i,j,k,ir,1)<0.) write (*,*) 'd negative',i,qm(l,i,j,k,ir,1)
              if (qm(l,i,j,k,ip,1)<0.) write (*,*) 'p negative',i,qm(l,i,j,k,ip,1)

              qm(l,i,j,k,ir,1) = max(smallr, qm(l,i,j,k,ir,1))
              qm(l,i,j,k,ip,1) = max(smallp, qm(l,i,j,k,ip,1))


              ! Top state at bottom interface
              qp(l,i,j,k,ir,2) = r - dry
              qp(l,i,j,k,iu,2) = u - duy
              qp(l,i,j,k,iv,2) = v - dvy
              qp(l,i,j,k,iw,2) = w - dwy
              qp(l,i,j,k,ip,2) = p - dpy
              if (qp(l,i,j,k,ir,2)<0.) write (*,*) 'd negative',i,qp(l,i,j,k,ir,2)
             if (qp(l,i,j,k,ip,2)<0.) write (*,*) 'p negative',i,qp(l,i,j,k,ip,2)

              qp(l,i,j,k,ir,2) = max(smallr, qp(l,i,j,k,ir,2))
              qp(l,i,j,k,ip,2) = max(smallp, qp(l,i,j,k,ip,2))


              ! Bottom state at top interface
              qm(l,i,j,k,ir,2) = r + dry
              qm(l,i,j,k,iu,2) = u + duy
              qm(l,i,j,k,iv,2) = v + dvy
              qm(l,i,j,k,iw,2) = w + dwy
              qm(l,i,j,k,ip,2) = p + dpy
              if (qm(l,i,j,k,ir,2)<0.) write (*,*) 'd negative',i,qm(l,i,j,k,ir,2)
             if (qm(l,i,j,k,ip,2)<0.) write (*,*) 'p negative',i,qm(l,i,j,k,ip,2)

              qm(l,i,j,k,ir,2) = max(smallr, qm(l,i,j,k,ir,2))
              qm(l,i,j,k,ip,2) = max(smallp, qm(l,i,j,k,ip,2))

              velsq=qp(l,i,j,k,iu,1)**2+qp(l,i,j,k,iv,1)**2+qp(l,i,j,k,iw,1)**2
              if ((velsq>1d0) .or.(qp(l,i,j,k,ir,1)<0d0).or.(qp(l,i,j,k,ip,1)<0d0)) then
                 qp(l,i,j,k,ir:ip,1)=q(l,i,j,k,ir:ip)
                 qm(l,i,j,k,ir:ip,1)=q(l,i,j,k,ir:ip)
                 qp(l,i,j,k,ir:ip,2)=q(l,i,j,k,ir:ip)
                 qm(l,i,j,k,ir:ip,2)=q(l,i,j,k,ir:ip)
              endif
              velsq=qm(l,i,j,k,iu,1)**2+qm(l,i,j,k,iv,1)**2+qm(l,i,j,k,iw,1)**2
              if ((velsq>1d0).or.(qm(l,i,j,k,ir,1)<0d0).or.(qm(l,i,j,k,ip,1)<0d0)) then
                 qp(l,i,j,k,ir:ip,1)=q(l,i,j,k,ir:ip)
                 qm(l,i,j,k,ir:ip,1)=q(l,i,j,k,ir:ip)
                 qp(l,i,j,k,ir:ip,2)=q(l,i,j,k,ir:ip)
                 qm(l,i,j,k,ir:ip,2)=q(l,i,j,k,ir:ip)
              endif
              velsq=qp(l,i,j,k,iu,2)**2+qp(l,i,j,k,iv,2)**2+qp(l,i,j,k,iw,2)**2
              if ((velsq>1d0) .or.(qp(l,i,j,k,ir,2)<0d0).or.(qp(l,i,j,k,ip,2)<0d0)) then
                 qp(l,i,j,k,ir:ip,1)=q(l,i,j,k,ir:ip)
                 qm(l,i,j,k,ir:ip,1)=q(l,i,j,k,ir:ip)
                 qp(l,i,j,k,ir:ip,2)=q(l,i,j,k,ir:ip)
                 qm(l,i,j,k,ir:ip,2)=q(l,i,j,k,ir:ip)
              endif
              velsq=qm(l,i,j,k,iu,2)**2+qm(l,i,j,k,iv,2)**2+qm(l,i,j,k,iw,2)**2
              if ((velsq>1d0)  .or.(qm(l,i,j,k,ir,2)<0d0).or.(qm(l,i,j,k,ip,2)<0d0)) then
                 qp(l,i,j,k,ir:ip,1)=q(l,i,j,k,ir:ip)
                 qm(l,i,j,k,ir:ip,1)=q(l,i,j,k,ir:ip)
                 qp(l,i,j,k,ir:ip,2)=q(l,i,j,k,ir:ip)
                 qm(l,i,j,k,ir:ip,2)=q(l,i,j,k,ir:ip)
              endif
              ! passive scalars
              do n = 6, nvar
                 a   = q(l,i,j,k,n)       ! Cell centered values
                 u   = q(l,i,j,k,iu)
                 v   = q(l,i,j,k,iv)
                 dax = dq(l,i,j,k,n,1)    ! TVD slopes
                 day = dq(l,i,j,k,n,2)
                 sa0 = -u*dax-v*day       ! Source terms
                 qp(l,i,j,k,n,1) = a - half*dax + sa0*dtdx*half*lor   ! Right state
                 qm(l,i,j,k,n,1) = a + half*dax + sa0*dtdx*half*lor   ! Left state
                 qp(l,i,j,k,n,2) = a - half*day + sa0*dtdy*half*lor   ! Top state
                 qm(l,i,j,k,n,2) = a + half*day + sa0*dtdy*half*lor   ! Bottom state
              end do
           end do
        end do
     end do
  end do

end subroutine trace2d
#endif
!###########################################################
!###########################################################
!###########################################################
!###########################################################
#if NDIM==3
! Trace in 3 dimensions.
subroutine  trace3d(q,dq,qm,qp,dx,dy,dz,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use const

  implicit none

  real(dp) :: dx,dy,dz,dt
  integer :: ngrid

  real(dp), dimension(1:nvector,iu1:iu2  ,ju1:ju2  ,ku1:ku2  ,1:nvar       ) :: q
  real(dp), dimension(1:nvector,iu1:iu2  ,ju1:ju2  ,ku1:ku2  ,1:nvar,1:ndim) :: dq
  real(dp), dimension(1:nvector,iu1:iu2  ,ju1:ju2  ,ku1:ku2  ,1:nvar,1:ndim) :: qm
  real(dp), dimension(1:nvector,iu1:iu2  ,ju1:ju2  ,ku1:ku2  ,1:nvar,1:ndim) :: qp

  ! Declare local variables

  integer :: i,j,k,l,n
  integer :: ilo,ihi,jlo,jhi,klo,khi
  integer :: ir,iu,iv,iw,ip
  real(dp)    :: dtdx,dtdy,dtdz,xc,xL,xR,yL,yR,smallp
  real(dp)    :: lor,lorsq,h,kappa,chi,f,ni,cs2,vtot,velsq,entho,tau,ka
  real(dp)    :: r,u,v,w,p
  real(dp)    :: drx,dux,dvx,dwx,dpx
  real(dp)    :: dry,duy,dvy,dwy,dpy
  real(dp)    :: drz,duz,dvz,dwz,dpz
  real(dp)    :: sr0,su0,sv0,sw0,sp0
  real(dp)    :: a, dax,day,daz, sa0
  real(dp)    :: a11,a12,a15,a22,a25,a32,a33,a35,a42,a44,a45,a52,a55
  real(dp)    :: b11,b13,b15,b22,b23,b25,b33,b35,b43,b44,b45,b53,b55
  real(dp)    :: c11,c14,c15,c22,c24,c25,c33,c34,c35,c44,c45,c54,c55
  real(dp)::vmax=0d0

  smallp = smallr*smallc**2/gamma

  ir=1; iu=2; iv=3; iw=4; ip=5

  ilo=min(1,iu1+1); ihi=max(1,iu2-1)
  jlo=min(1,ju1+1); jhi=max(1,ju2-1)
  klo=min(1,ku1+1); khi=max(1,ku2-1)


  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           do l = 1, ngrid

              dtdx = dt/dx
              dtdy = dt/dy
              dtdz = dt/dz



              ! Cell centered values
              r =    q(l,i,j,k,ir)
              u =    q(l,i,j,k,iu)
              v =    q(l,i,j,k,iv)
              w =    q(l,i,j,k,iw)
              p =    q(l,i,j,k,ip)

              ! Cell centered TVD slopes in X direction
              drx = half * dq(l,i,j,k,ir,1)
              dux = half * dq(l,i,j,k,iu,1)
              dvx = half * dq(l,i,j,k,iv,1)
              dwx = half * dq(l,i,j,k,iw,1)
              dpx = half * dq(l,i,j,k,ip,1)

              ! Cell centered TVD slopes in Y direction
              dry = half * dq(l,i,j,k,ir,2)
              duy = half * dq(l,i,j,k,iu,2)
              dvy = half * dq(l,i,j,k,iv,2)
              dwy = half * dq(l,i,j,k,iw,2)
              dpy = half * dq(l,i,j,k,ip,2)

              ! Cell centered TVD slopes in Y direction
              drz = half * dq(l,i,j,k,ir,3)
              duz = half * dq(l,i,j,k,iu,3)
              dvz = half * dq(l,i,j,k,iv,3)
              dwz = half * dq(l,i,j,k,iw,3)
              dpz = half * dq(l,i,j,k,ip,3)

              ! Source terms (including transverse derivatives)
              velsq=u**2+v**2+w**2 ; lorsq=1d0/(1d0-velsq) ; lor=sqrt(lorsq)
              vtot=sqrt(u**2+v**2+w**2)

              if (eos .eq. 'TM') then
                 tau=p/r
                 ka=9d0/4d0*tau**2+1
                 h=5d0/2d0*tau+sqrt(ka)
                 chi=-5d0/2d0*tau/r-9d0/4d0/sqrt(ka)*tau**2/r
                 kappa=5d0/2d0/r+9d0/4d0/sqrt(ka)*p/r**2
              else
                 entho = gamma/(gamma-one)
                 h     = 1.0d0+entho*p/r   !enthalpy
                 kappa = entho*1d0/r      !dh/dp
                 chi   = -entho*p/r**2     !dh/dr
              endif

              f=lorsq*(-h+h*r*kappa+r*chi*velsq)

              a11= u*f
              a12= h*r*lorsq*(r*kappa-1)
              a15= u*(-r*lorsq*kappa+r*lorsq*velsq*kappa+1)

              a22= u*lorsq*(-h+r*chi+r*h*kappa)
              a25= (r*lorsq*kappa*h-r*h*kappa*lorsq*u**2-h-h*lorsq*(v**2+w**2)+lorsq*r*chi*(v**2+w**2))/(r*lorsq*h)

              a32= v*r*chi
              a33= u*f
              a35= -(u*v*(-h+r*chi+r*h*kappa))/(r*h)

              a42= w*r*chi
              a44= u*f
              a45= -(u*w*(-h+r*chi+r*h*kappa))/(r*h)

              a52= -r**2*h*lorsq*chi
              a55= lorsq*u*(-h+r*chi+h*r*kappa)


              b11= v*f
              b13= h*r*lorsq*(r*kappa-1)
              b15= v*(-r*lorsq*kappa+r*lorsq*velsq*kappa+1)

              b22= v*f
              b23= u*r*chi
              b25=-(v*u*(-h+r*chi+r*h*kappa))/(r*h)

              b33= v*lorsq*(-h+r*chi+r*h*kappa)
              b35= (r*lorsq*kappa*h-r*h*kappa*lorsq*v**2-h-h*lorsq*(u**2+w**2)+lorsq*r*chi*(u**2+w**2))/(r*lorsq*h)

              b43= w*r*chi
              b44= v*f
              b45=-(v*w*(-h+r*chi+r*h*kappa))/(r*h)

              b53=-r**2*h*lorsq*chi
              b55= lorsq*v*(-h+r*chi+h*r*kappa)



              c11= w*f
              c14= h*r*lorsq*(r*kappa-1)
              c15= w*(-r*lorsq*kappa+r*lorsq*velsq*kappa+1)

              c22= w*f
              c24= u*r*chi
              c25=-(u*w*(-h+r*chi+r*h*kappa))/(r*h)

              c33= w*f
              c34= v*r*chi
              c35=-(v*w*(-h+r*chi+r*h*kappa))/(r*h)

              c44= w*lorsq*(r*h*kappa-h+r*chi)
              c45=(r*lorsq*kappa*h-r*h*kappa*lorsq*w**2-h-h*lorsq*(u**2+v**2)+lorsq*r*chi*(u**2+v**2))/(r*lorsq*h)

              c54=-r**2*h*lorsq*chi
              c55= lorsq*w*(-h+r*chi+h*r*kappa)


              ! Source terms (including transverse derivatives)
              sr0 = (-a11*drx-a12*dux                -a15*dpx)*dtdx/f + (-b11*dry        -b13*dvy        -b15*dpy)*dtdy/f + (-c11*drz                -c14*dwz-c15*dpz)*dtdz/f
              su0 = (        -a22*dux                -a25*dpx)*dtdx/f + (        -b22*duy-b23*dvy        -b25*dpy)*dtdy/f + (        -c22*duz        -c24*dwz-c25*dpz)*dtdz/f
              sv0 = (        -a32*dux-a33*dvx        -a35*dpx)*dtdx/f + (                -b33*dvy        -b35*dpy)*dtdy/f + (                -c33*dvz-c34*dwz-c35*dpz)*dtdz/f
              sw0 = (        -a42*dux        -a44*dwx-a45*dpx)*dtdx/f + (                -b43*dvy-b44*dwy-b45*dpy)*dtdy/f + (                        -c44*dwz-c45*dpz)*dtdz/f
              sp0 = (        -a52*dux                -a55*dpx)*dtdx/f + (                -b53*dvy        -b55*dpy)*dtdy/f + (                        -c54*dwz-c55*dpz)*dtdz/f


              !Update in time the  primitive variables
              r = r + sr0
              u = u + su0
              v = v + sv0
              w = w + sw0
              p = p + sp0

              ! Face averaged right state at left interface
              qp(l,i,j,k,ir,1) = r - drx
              qp(l,i,j,k,iu,1) = u - dux
              qp(l,i,j,k,iv,1) = v - dvx
              qp(l,i,j,k,iw,1) = w - dwx
              qp(l,i,j,k,ip,1) = p - dpx
              qp(l,i,j,k,ir,1) = max(smallr, qp(l,i,j,k,ir,1))
              qp(l,i,j,k,ip,1) = max(smallp, qp(l,i,j,k,ip,1))

              ! Face averaged left state at right interface
              qm(l,i,j,k,ir,1) = r + drx
              qm(l,i,j,k,iu,1) = u + dux
              qm(l,i,j,k,iv,1) = v + dvx
              qm(l,i,j,k,iw,1) = w + dwx
              qm(l,i,j,k,ip,1) = p + dpx
              qm(l,i,j,k,ir,1) = max(smallr, qm(l,i,j,k,ir,1))
              qm(l,i,j,k,ip,1) = max(smallp, qm(l,i,j,k,ip,1))

              ! Face averaged top state at bottom interface
              qp(l,i,j,k,ir,2) = r - dry
              qp(l,i,j,k,iu,2) = u - duy
              qp(l,i,j,k,iv,2) = v - dvy
              qp(l,i,j,k,iw,2) = w - dwy
              qp(l,i,j,k,ip,2) = p - dpy
              qp(l,i,j,k,ir,2) = max(smallr, qp(l,i,j,k,ir,2))
              qp(l,i,j,k,ip,2) = max(smallp, qp(l,i,j,k,ip,2))

              ! Face averaged bottom state at top interface
              qm(l,i,j,k,ir,2) = r + dry
              qm(l,i,j,k,iu,2) = u + duy
              qm(l,i,j,k,iv,2) = v + dvy
              qm(l,i,j,k,iw,2) = w + dwy
              qm(l,i,j,k,ip,2) = p + dpy
              qm(l,i,j,k,ir,2) = max(smallr, qm(l,i,j,k,ir,2))
              qm(l,i,j,k,ip,2) = max(smallp, qm(l,i,j,k,ip,2))

              ! Face averaged front state at back interface
              qp(l,i,j,k,ir,3) = r - drz
              qp(l,i,j,k,iu,3) = u - duz
              qp(l,i,j,k,iv,3) = v - dvz
              qp(l,i,j,k,iw,3) = w - dwz
              qp(l,i,j,k,ip,3) = p - dpz
              qp(l,i,j,k,ir,3) = max(smallr, qp(l,i,j,k,ir,3))
              qp(l,i,j,k,ip,3) = max(smallp, qp(l,i,j,k,ip,3))

              ! Face averaged back state at front interface
              qm(l,i,j,k,ir,3) = r + drz
              qm(l,i,j,k,iu,3) = u + duz
              qm(l,i,j,k,iv,3) = v + dvz
              qm(l,i,j,k,iw,3) = w + dwz
              qm(l,i,j,k,ip,3) = p + dpz
              qm(l,i,j,k,ir,3) = max(smallr, qm(l,i,j,k,ir,3))
              qm(l,i,j,k,ip,3) = max(smallp, qm(l,i,j,k,ip,3))


              velsq=qp(l,i,j,k,iu,1)**2+qp(l,i,j,k,iv,1)**2+qp(l,i,j,k,iw,1)**2
              if ((velsq>1d0) .or.(qp(l,i,j,k,ir,1)<0d0).or.(qp(l,i,j,k,ip,1)<0d0)) then
                 qp(l,i,j,k,ir:ip,1)=q(l,i,j,k,ir:ip)
                 qm(l,i,j,k,ir:ip,1)=q(l,i,j,k,ir:ip)
                 qp(l,i,j,k,ir:ip,2)=q(l,i,j,k,ir:ip)
                 qm(l,i,j,k,ir:ip,2)=q(l,i,j,k,ir:ip)
                 qp(l,i,j,k,ir:ip,3)=q(l,i,j,k,ir:ip)
                 qm(l,i,j,k,ir:ip,3)=q(l,i,j,k,ir:ip)
              endif
              velsq=qp(l,i,j,k,iu,1)**2+qp(l,i,j,k,iv,1)**2+qp(l,i,j,k,iw,1)**2
              if (sqrt(velsq) >vmax) then
                 vmax=sqrt(velsq)
              endif

              velsq=qm(l,i,j,k,iu,1)**2+qm(l,i,j,k,iv,1)**2+qm(l,i,j,k,iw,1)**2
              if ((velsq>1d0).or.(qm(l,i,j,k,ir,1)<0d0).or.(qm(l,i,j,k,ip,1)<0d0)) then
                 qp(l,i,j,k,ir:ip,1)=q(l,i,j,k,ir:ip)
                 qm(l,i,j,k,ir:ip,1)=q(l,i,j,k,ir:ip)
                 qp(l,i,j,k,ir:ip,2)=q(l,i,j,k,ir:ip)
                 qm(l,i,j,k,ir:ip,2)=q(l,i,j,k,ir:ip)
                 qp(l,i,j,k,ir:ip,3)=q(l,i,j,k,ir:ip)
                 qm(l,i,j,k,ir:ip,3)=q(l,i,j,k,ir:ip)
              endif
              velsq=qm(l,i,j,k,iu,1)**2+qm(l,i,j,k,iv,1)**2+qm(l,i,j,k,iw,1)**2
              if (sqrt(velsq) >vmax) then
                 vmax=sqrt(velsq)

              endif

              velsq=qp(l,i,j,k,iu,2)**2+qp(l,i,j,k,iv,2)**2+qp(l,i,j,k,iw,2)**2
              if ((velsq>1d0) .or.(qp(l,i,j,k,ir,2)<0d0).or.(qp(l,i,j,k,ip,2)<0d0)) then
                 qp(l,i,j,k,ir:ip,1)=q(l,i,j,k,ir:ip)
                 qm(l,i,j,k,ir:ip,1)=q(l,i,j,k,ir:ip)
                 qp(l,i,j,k,ir:ip,2)=q(l,i,j,k,ir:ip)
                 qm(l,i,j,k,ir:ip,2)=q(l,i,j,k,ir:ip)
                 qp(l,i,j,k,ir:ip,3)=q(l,i,j,k,ir:ip)
                 qm(l,i,j,k,ir:ip,3)=q(l,i,j,k,ir:ip)
              endif
              velsq=qp(l,i,j,k,iu,2)**2+qp(l,i,j,k,iv,2)**2+qp(l,i,j,k,iw,2)**2
              if (sqrt(velsq) >vmax) then
                 vmax=sqrt(velsq)
              endif

              velsq=qm(l,i,j,k,iu,2)**2+qm(l,i,j,k,iv,2)**2+qm(l,i,j,k,iw,2)**2
              if ((velsq>1d0)  .or.(qm(l,i,j,k,ir,2)<0d0).or.(qm(l,i,j,k,ip,2)<0d0)) then
                 qp(l,i,j,k,ir:ip,1)=q(l,i,j,k,ir:ip)
                 qm(l,i,j,k,ir:ip,1)=q(l,i,j,k,ir:ip)
                 qp(l,i,j,k,ir:ip,2)=q(l,i,j,k,ir:ip)
                 qm(l,i,j,k,ir:ip,2)=q(l,i,j,k,ir:ip)
                 qp(l,i,j,k,ir:ip,3)=q(l,i,j,k,ir:ip)
                 qm(l,i,j,k,ir:ip,3)=q(l,i,j,k,ir:ip)
              endif
              velsq=qm(l,i,j,k,iu,2)**2+qm(l,i,j,k,iv,2)**2+qm(l,i,j,k,iw,2)**2
              if (sqrt(velsq) >vmax) then
                 vmax=sqrt(velsq)
              endif

              velsq=qp(l,i,j,k,iu,3)**2+qp(l,i,j,k,iv,3)**2+qp(l,i,j,k,iw,3)**2
              if ((velsq>1d0) .or.(qp(l,i,j,k,ir,2)<0d0).or.(qp(l,i,j,k,ip,2)<0d0)) then
                 qp(l,i,j,k,ir:ip,1)=q(l,i,j,k,ir:ip)
                 qm(l,i,j,k,ir:ip,1)=q(l,i,j,k,ir:ip)
                 qp(l,i,j,k,ir:ip,2)=q(l,i,j,k,ir:ip)
                 qm(l,i,j,k,ir:ip,2)=q(l,i,j,k,ir:ip)
                 qp(l,i,j,k,ir:ip,3)=q(l,i,j,k,ir:ip)
                 qm(l,i,j,k,ir:ip,3)=q(l,i,j,k,ir:ip)
              endif
              velsq=qp(l,i,j,k,iu,3)**2+qp(l,i,j,k,iv,3)**2+qp(l,i,j,k,iw,3)**2
              if (sqrt(velsq) >vmax) then
                 vmax=sqrt(velsq)
              endif

              velsq=qm(l,i,j,k,iu,3)**2+qm(l,i,j,k,iv,3)**2+qm(l,i,j,k,iw,3)**2
              if ((velsq>1d0)  .or.(qm(l,i,j,k,ir,3)<0d0).or.(qm(l,i,j,k,ip,3)<0d0)) then
                 qp(l,i,j,k,ir:ip,1)=q(l,i,j,k,ir:ip)
                 qm(l,i,j,k,ir:ip,1)=q(l,i,j,k,ir:ip)
                 qp(l,i,j,k,ir:ip,2)=q(l,i,j,k,ir:ip)
                 qm(l,i,j,k,ir:ip,2)=q(l,i,j,k,ir:ip)
                 qp(l,i,j,k,ir:ip,3)=q(l,i,j,k,ir:ip)
                 qm(l,i,j,k,ir:ip,3)=q(l,i,j,k,ir:ip)
              endif
              velsq=qm(l,i,j,k,iu,3)**2+qm(l,i,j,k,iv,3)**2+qm(l,i,j,k,iw,3)**2
              if (sqrt(velsq) >vmax) then
                 vmax=sqrt(velsq)
              endif


              ! Passive scalars
              do n = 6, nvar
                 a   = q(l,i,j,k,n)       ! Cell centered values
                 u   = q(l,i,j,k,iu)
                 v   = q(l,i,j,k,iv)
                 w   = q(l,i,j,k,iw)
                 dax = dq(l,i,j,k,n,1)    ! TVD slopes
                 day = dq(l,i,j,k,n,2)
                 daz = dq(l,i,j,k,n,3)
                 sa0 = -u*dax-v*day-w*daz     ! Source terms
                 qp(l,i,j,k,n,1) = a - half*dax + sa0*dtdx*half*lor  ! Right state
                 qm(l,i,j,k,n,1) = a + half*dax + sa0*dtdx*half*lor  ! Left state
                 qp(l,i,j,k,n,2) = a - half*day + sa0*dtdy*half*lor  ! Bottom state
                 qm(l,i,j,k,n,2) = a + half*day + sa0*dtdy*half*lor  ! Upper state
                 qp(l,i,j,k,n,3) = a - half*daz + sa0*dtdz*half*lor  ! Front state
                 qm(l,i,j,k,n,3) = a + half*daz + sa0*dtdz*half*lor  ! Back state
              end do
           end do
        end do
     end do
  end do


end subroutine trace3d
#endif
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmpflxm(qm,im1,im2,jm1,jm2,km1,km2, &
     &             qp,ip1,ip2,jp1,jp2,kp1,kp2, &
     &                ilo,ihi,jlo,jhi,klo,khi, ln,lt1,lt2, &
     &            flx,tmp,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer ::ngrid
  integer ::ln,lt1,lt2
  integer ::im1,im2,jm1,jm2,km1,km2
  integer ::ip1,ip2,jp1,jp2,kp1,kp2
  integer ::ilo,ihi,jlo,jhi,klo,khi
  real(dp),dimension(1:nvector,im1:im2,jm1:jm2,km1:km2,1:nvar,1:ndim)::qm
  real(dp),dimension(1:nvector,ip1:ip2,jp1:jp2,kp1:kp2,1:nvar,1:ndim)::qp
  real(dp),dimension(1:nvector,ip1:ip2,jp1:jp2,kp1:kp2,1:nvar)::flx
  real(dp),dimension(1:nvector,ip1:ip2,jp1:jp2,kp1:kp2,1:2)::tmp

  ! local variables
  integer ::i, j, k, n, l, idim, xdim
  real(dp)::entho
  real(dp),dimension(1:nvar),save::qleft,qright,qgdnv,fgdnv

  entho=one/(gamma-one)
  xdim=ln-1

  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           do l = 1, ngrid

              ! Left state
              qleft (1) = qm(l,i,j,k,1  ,xdim) ! Mass density
              qleft (2) = qm(l,i,j,k,5  ,xdim) ! Pressure
              qleft (3) = qm(l,i,j,k,ln ,xdim) ! Normal velocity
              qleft (4) = qm(l,i,j,k,lt1,xdim) ! Tangential velocity 1
              qleft (5) = qm(l,i,j,k,lt2,xdim) ! Tangential velocity 2

              ! Right state
              qright(1) = qp(l,i,j,k,1  ,xdim) ! Mass density
              qright(2) = qp(l,i,j,k,5  ,xdim) ! Pressure
              qright(3) = qp(l,i,j,k,ln ,xdim) ! Normal velocity
              qright(4) = qp(l,i,j,k,lt1,xdim) ! Tangential velocity 1
              qright(5) = qp(l,i,j,k,lt2,xdim) ! Tangential velocity 2

              ! Passive scalars
              do n = 6, nvar
                 qleft (n) = qm(l,i,j,k,n,xdim)
                 qright(n) = qp(l,i,j,k,n,xdim)
              end do


              ! Solve 1D Riemann problem
              if (riemann.eq.'hll')then
                 call riemann_hll     (qleft,qright,fgdnv)
              else if (riemann.eq.'hllc')then
                 call riemann_hllc     (qleft,qright,fgdnv)
              else
                 write(*,*)'unknown Riemann solver'
                 stop
              end if


              ! Output fluxes
              flx(l,i,j,k,1  ) = fgdnv(1)  ! Mass density
              flx(l,i,j,k,5  ) = fgdnv(2)  ! Total energy
              flx(l,i,j,k,ln ) = fgdnv(3)  ! Normal momentum
              flx(l,i,j,k,lt1) = fgdnv(4)  ! Transverse momentum 1
              flx(l,i,j,k,lt2) = fgdnv(5)  ! Transverse momentum 2

              ! Passive scalars
              do n = 6, nvar
                 flx(l,i,j,k,n) = fgdnv(n)
              end do

              ! Normal velocity estimate
              tmp(l,i,j,k,1) = half*(qleft(3)+qright(3))
              ! Internal energy flux
              if(fgdnv(1)>zero)then
                 tmp(l,i,j,k,2) = qleft (2)/qleft (1)*entho*fgdnv(1)
              else
                 tmp(l,i,j,k,2) = qright(2)/qright(1)*entho*fgdnv(1)
              end if

           end do
        end do
     end do
  end do

end subroutine cmpflxm
!###########################################################
!###########################################################
!###########################################################
!###########################################################


!  Subroutine CTOPRIM
!
!> Converts conservative variables to primitive in rhd
! Modified to include special relativity. The equation of state is an extension of
!classical EOS , and may fail for adiabatic indexes >2.


subroutine ctoprim(uin,q,gravin,dt,ngrid)
  use hydro_parameters
  use amr_parameters
  use amr_commons
  implicit none

  integer ::ngrid
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::uin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim)::gravin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q
!  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::c

  integer :: i,j,k,l,n,idim
  real(dp) :: dt,time
  real(dp) :: lor,entho ! Lorentz factor
  real(dp) :: D,M,E,Mx,My,Mz,u2,Xsi,R
  real(dp)::half=0.5
  real(dp)::small_bigD=1d-12
  real(dp)::D2,M2,E2,tau,h
! a voir
  real(dp)::eken
  real(dp)::eint, smalle, dtxhalf, oneonrho,rho

  if (verbose) write (*,*) 'Entering ctoprim...'

  ! convert to primitive variables
  do k = ku1, ku2
     do j = ju1, ju2
        do i = iu1, iu2
           do l = 1 ,ngrid

             ! Compute density
             D = uin(l,i,j,k,1)
             ! Compute momentum
             Mx=uin(l,i,j,k,2) ; My=uin(l,i,j,k,3) ; Mz=uin(l,i,j,k,4)
             M = sqrt(Mx**2+My**2+Mz**2)
             ! Compute total energy
             E = uin(l,i,j,k,5)
             !Method from Mignone,McKinney,2007. Same as BS2011 except one uses E'=U-D and u^2=Lor^2*v^2
             if (D<0) then
               write(*,*),'D<0 in ctoprim'
                D=small_bigD
                uin(l,i,j,k,1)=D
             endif
             if (E<0) then
                E=sqrt(M**2+d**2+1d-8)
               write(*,*),'E<0 in ctoprim'
                uin(l,i,j,k,5)=E
             endif

             if (E**2<M**2+D**2) then

               write (*,*) 'Switch...ctoprim'
                E=sqrt(M**2+d**2+1d-8)
                uin(l,i,j,k,5)=E
             endif

             if ( M .eq. 0) then
                q(l,i,j,k,1) = D
                q(l,i,j,k,2) = 0d0
                q(l,i,j,k,3) = 0d0
                q(l,i,j,k,4) = 0d0
                if (eos .eq. 'TM') then
                   q(l,i,j,k,5) =(E**2-D**2)/3d0/E
                else
                   q(l,i,j,k,5)=(E-D)*(gamma-1d0)
                endif
                lor=1d0
             else

                call Newton_Raphson_Mignone(D,M,E,gamma,R)

                ! Compute the Lorentz factor
                u2  = M**2.0d0/(R**2.0d0-M**2.0d0)
                lor = (1.0d0+u2)**(1d0/2d0)

                ! Compute the density
                q(l,i,j,k,1) = D/lor

                ! compute velocities
                q(l,i,j,k,2) = Mx/R!max(1.0e-10,Mx/R)
                q(l,i,j,k,3) = My/R!max(1.0e-10,My/R)
                q(l,i,j,k,4) = Mz/R!max(1.0e-10,Mz/R)

                ! Compute pressure
                Xsi=((R-D)-u2/(lor+1d0)*D)/lor**2
                if (eos .eq. 'TM') then
                   rho=q(l,i,j,k,1)
                   q(l,i,j,k,5)=(2d0*xsi*(xsi+2d0*rho))/(5d0*(xsi+rho)+sqrt(9d0*xsi**2+18d0*rho*xsi+25d0*rho**2))
                else
                   q(l,i,j,k,5)=(gamma-1d0)/gamma*Xsi
                endif
             endif
             if ((q(l,i,j,k,1)<0d0).or.(q(l,i,j,k,5)<0d0).or.E<0d0) then
                write(*,*) 'negative pressure or density'
                stop
             endif

             ! Passive scalar
             do n = 6, nvar
                oneonrho = 1d0/q(l,i,j,k,1)
                q(l,i,j,k,n) = uin(l,i,j,k,n)*oneonrho/lor
             end do
             tau=q(l,i,j,k,5)/q(l,i,j,k,1)
             h=5d0/2d0*tau+3d0/2d0*sqrt(tau**2+4d0/9d0)
             D2= q(l,i,j,k,1)*lor
             M2=lor**2*q(l,i,j,k,1)*h*sqrt(q(l,i,j,k,2)**2+q(l,i,j,k,3)**2+q(l,i,j,k,4)**2)
             E2=lor**2*q(l,i,j,k,1)*h-q(l,i,j,k,5)
!             if ((abs(D-d2)/D2 .gt. 1e-5 )) write(*,*),D,D2,'D'
!             if ((abs(Mx-M2)/M2 .gt. 1e-5))  write (*,*)M,M2,'M'
!             if ((abs(E-E2)/E2 .gt. 1e-5))  write(*,*),E,E2,'E'


        end do
      end do
   end do
  end do


end subroutine ctoprim

!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine uslope(q,dq,dx,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer::ngrid
  real(dp)::dx,dt
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::dq

  ! local arrays
  integer::i, j, k, l, n
  real(dp)::dsgn, dlim, dcen, dlft, drgt, slop
  real(dp)::dfll,dflm,dflr,dfml,dfmm,dfmr,dfrl,dfrm,dfrr
  real(dp)::dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl
  real(dp)::dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm
  real(dp)::dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr
  real(dp)::vmin,vmax,dfx,dfy,dfz,dff
  integer::ilo,ihi,jlo,jhi,klo,khi

  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)

  if(slope_type==0)then
     dq=zero
     return
  end if

#if NDIM==1
  do n = 1, nvar
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              if(slope_type==1.or.slope_type==2.or.slope_type==3)then  ! minmod or average
                 do l = 1, ngrid
                    dlft = MIN(slope_type,2)*(q(l,i  ,j,k,n) - q(l,i-1,j,k,n))
                    drgt = MIN(slope_type,2)*(q(l,i+1,j,k,n) - q(l,i  ,j,k,n))
                    dcen = half*(dlft+drgt)/MIN(slope_type,2)
                    dsgn = sign(one, dcen)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,1) = dsgn*min(dlim,abs(dcen))
                 end do
              else if(slope_type==4)then ! superbee
                 do l = 1, ngrid
                    dcen = q(l,i,j,k,2)*dt/dx
                    dlft = two/(one+dcen)*(q(l,i,j,k,n)-q(l,i-1,j,k,n))
                    drgt = two/(one-dcen)*(q(l,i+1,j,k,n)-q(l,i,j,k,n))
                    dcen = half*(q(l,i+1,j,k,n)-q(l,i-1,j,k,n))
                    dsgn = sign(one, dlft)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,1) = dsgn*dlim !min(dlim,abs(dcen))
                 end do
              else if(slope_type==5)then ! ultrabee
                 if(n==1)then
                    do l = 1, ngrid
                       dcen = q(l,i,j,k,2)*dt/dx
                       if(dcen>=0)then
                          dlft = two/(zero+dcen+1d-10)*(q(l,i,j,k,n)-q(l,i-1,j,k,n))
                          drgt = two/(one -dcen      )*(q(l,i+1,j,k,n)-q(l,i,j,k,n))
                       else
                          dlft = two/(one +dcen      )*(q(l,i,j,k,n)-q(l,i-1,j,k,n))
                          drgt = two/(zero-dcen+1d-10)*(q(l,i+1,j,k,n)-q(l,i,j,k,n))
                       endif
                       dsgn = sign(one, dlft)
                       slop = min(abs(dlft),abs(drgt))
                       dlim = slop
                       dcen = half*(q(l,i+1,j,k,n)-q(l,i-1,j,k,n))
                       if((dlft*drgt)<=zero)dlim=zero
                       dq(l,i,j,k,n,1) = dsgn*dlim !min(dlim,abs(dcen))
                    end do
                 else
                    do l = 1, ngrid
                       dq(l,i,j,k,n,1) = 0.0
                    end do
                 end if
              else if(slope_type==6)then ! unstable
                 if(n==1)then
                    do l = 1, ngrid
                       dlft = (q(l,i,j,k,n)-q(l,i-1,j,k,n))
                       drgt = (q(l,i+1,j,k,n)-q(l,i,j,k,n))
                       slop = 0.5*(dlft+drgt)
                       dlim = slop
                       dq(l,i,j,k,n,1) = dlim
                    end do
                 else
                    do l = 1, ngrid
                       dq(l,i,j,k,n,1) = 0.0
                    end do
                 end if
              else
                 write(*,*)'Unknown slope type'
                 stop
              end if
           end do
        end do
     end do
  end do
#endif

#if NDIM==2
  if(slope_type==1.or.slope_type==2)then  ! minmod or average
     do n = 1, nvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 ! slopes in first coordinate direction
                 do l = 1, ngrid
                    dlft = slope_type*(q(l,i  ,j,k,n) - q(l,i-1,j,k,n))
                    drgt = slope_type*(q(l,i+1,j,k,n) - q(l,i  ,j,k,n))
                    dcen = half*(dlft+drgt)/slope_type
                    dsgn = sign(one, dcen)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,1) = dsgn*min(dlim,abs(dcen))
                 end do
                 ! slopes in second coordinate direction
                 do l = 1, ngrid
                    dlft = slope_type*(q(l,i,j  ,k,n) - q(l,i,j-1,k,n))
                    drgt = slope_type*(q(l,i,j+1,k,n) - q(l,i,j  ,k,n))
                    dcen = half*(dlft+drgt)/slope_type
                    dsgn = sign(one,dcen)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,2) = dsgn*min(dlim,abs(dcen))
                 end do
              end do
           end do
        end do
     end do
  else if(slope_type==3)then ! positivity preserving 2d unsplit slope
     do n = 1, nvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 do l = 1, ngrid
                    dfll = q(l,i-1,j-1,k,n)-q(l,i,j,k,n)
                    dflm = q(l,i-1,j  ,k,n)-q(l,i,j,k,n)
                    dflr = q(l,i-1,j+1,k,n)-q(l,i,j,k,n)
                    dfml = q(l,i  ,j-1,k,n)-q(l,i,j,k,n)
                    dfmm = q(l,i  ,j  ,k,n)-q(l,i,j,k,n)
                    dfmr = q(l,i  ,j+1,k,n)-q(l,i,j,k,n)
                    dfrl = q(l,i+1,j-1,k,n)-q(l,i,j,k,n)
                    dfrm = q(l,i+1,j  ,k,n)-q(l,i,j,k,n)
                    dfrr = q(l,i+1,j+1,k,n)-q(l,i,j,k,n)

                    vmin = min(dfll,dflm,dflr,dfml,dfmm,dfmr,dfrl,dfrm,dfrr)
                    vmax = max(dfll,dflm,dflr,dfml,dfmm,dfmr,dfrl,dfrm,dfrr)

                    dfx  = half*(q(l,i+1,j,k,n)-q(l,i-1,j,k,n))
                    dfy  = half*(q(l,i,j+1,k,n)-q(l,i,j-1,k,n))
                    dff  = half*(abs(dfx)+abs(dfy))

                    if(dff>zero)then
                       slop = min(one,min(abs(vmin),abs(vmax))/dff)
                    else
                       slop = one
                    endif

                    dlim = slop

                    dq(l,i,j,k,n,1) = dlim*dfx
                    dq(l,i,j,k,n,2) = dlim*dfy

                 end do
              end do
           end do
        end do
     end do
  else
     write(*,*)'Unknown slope type'
     stop
  endif
#endif

#if NDIM==3
  if(slope_type==1)then  ! minmod
     do n = 1, nvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 ! slopes in first coordinate direction
                 do l = 1, ngrid
                    dlft = q(l,i  ,j,k,n) - q(l,i-1,j,k,n)
                    drgt = q(l,i+1,j,k,n) - q(l,i  ,j,k,n)
                    if((dlft*drgt)<=zero) then
                       dq(l,i,j,k,n,1) = zero
                    else if(dlft>0) then
                       dq(l,i,j,k,n,1) = min(dlft,drgt)
                    else
                       dq(l,i,j,k,n,1) = max(dlft,drgt)
                    end if
                 end do
                 ! slopes in second coordinate direction
                 do l = 1, ngrid
                    dlft = q(l,i,j  ,k,n) - q(l,i,j-1,k,n)
                    drgt = q(l,i,j+1,k,n) - q(l,i,j  ,k,n)
                    if((dlft*drgt)<=zero) then
                       dq(l,i,j,k,n,2) = zero
                    else if(dlft>0) then
                       dq(l,i,j,k,n,2) = min(dlft,drgt)
                    else
                       dq(l,i,j,k,n,2) = max(dlft,drgt)
                    end if
                 end do
                 ! slopes in third coordinate direction
                 do l = 1, ngrid
                    dlft = q(l,i,j,k  ,n) - q(l,i,j,k-1,n)
                    drgt = q(l,i,j,k+1,n) - q(l,i,j,k  ,n)
                    if((dlft*drgt)<=zero) then
                       dq(l,i,j,k,n,3) = zero
                    else if(dlft>0) then
                       dq(l,i,j,k,n,3) = min(dlft,drgt)
                    else
                       dq(l,i,j,k,n,3) = max(dlft,drgt)
                    end if
                 end do
              end do
           end do
        end do
     end do
  else if(slope_type==2)then ! moncen
     do n = 1, nvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 ! slopes in first coordinate direction
                 do l = 1, ngrid
                    dlft = slope_type*(q(l,i  ,j,k,n) - q(l,i-1,j,k,n))
                    drgt = slope_type*(q(l,i+1,j,k,n) - q(l,i  ,j,k,n))
                    dcen = half*(dlft+drgt)/slope_type
                    dsgn = sign(one, dcen)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,1) = dsgn*min(dlim,abs(dcen))
                 end do
                 ! slopes in second coordinate direction
                 do l = 1, ngrid
                    dlft = slope_type*(q(l,i,j  ,k,n) - q(l,i,j-1,k,n))
                    drgt = slope_type*(q(l,i,j+1,k,n) - q(l,i,j  ,k,n))
                    dcen = half*(dlft+drgt)/slope_type
                    dsgn = sign(one,dcen)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,2) = dsgn*min(dlim,abs(dcen))
                 end do
                 ! slopes in third coordinate direction
                 do l = 1, ngrid
                    dlft = slope_type*(q(l,i,j,k  ,n) - q(l,i,j,k-1,n))
                    drgt = slope_type*(q(l,i,j,k+1,n) - q(l,i,j,k  ,n))
                    dcen = half*(dlft+drgt)/slope_type
                    dsgn = sign(one,dcen)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,3) = dsgn*min(dlim,abs(dcen))
                 end do
              end do
           end do
        end do
     end do
  else if(slope_type==3)then ! positivity preserving 3d unsplit slope
     do n = 1, nvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 do l = 1, ngrid
                    dflll = q(l,i-1,j-1,k-1,n)-q(l,i,j,k,n)
                    dflml = q(l,i-1,j  ,k-1,n)-q(l,i,j,k,n)
                    dflrl = q(l,i-1,j+1,k-1,n)-q(l,i,j,k,n)
                    dfmll = q(l,i  ,j-1,k-1,n)-q(l,i,j,k,n)
                    dfmml = q(l,i  ,j  ,k-1,n)-q(l,i,j,k,n)
                    dfmrl = q(l,i  ,j+1,k-1,n)-q(l,i,j,k,n)
                    dfrll = q(l,i+1,j-1,k-1,n)-q(l,i,j,k,n)
                    dfrml = q(l,i+1,j  ,k-1,n)-q(l,i,j,k,n)
                    dfrrl = q(l,i+1,j+1,k-1,n)-q(l,i,j,k,n)

                    dfllm = q(l,i-1,j-1,k  ,n)-q(l,i,j,k,n)
                    dflmm = q(l,i-1,j  ,k  ,n)-q(l,i,j,k,n)
                    dflrm = q(l,i-1,j+1,k  ,n)-q(l,i,j,k,n)
                    dfmlm = q(l,i  ,j-1,k  ,n)-q(l,i,j,k,n)
                    dfmmm = q(l,i  ,j  ,k  ,n)-q(l,i,j,k,n)
                    dfmrm = q(l,i  ,j+1,k  ,n)-q(l,i,j,k,n)
                    dfrlm = q(l,i+1,j-1,k  ,n)-q(l,i,j,k,n)
                    dfrmm = q(l,i+1,j  ,k  ,n)-q(l,i,j,k,n)
                    dfrrm = q(l,i+1,j+1,k  ,n)-q(l,i,j,k,n)

                    dfllr = q(l,i-1,j-1,k+1,n)-q(l,i,j,k,n)
                    dflmr = q(l,i-1,j  ,k+1,n)-q(l,i,j,k,n)
                    dflrr = q(l,i-1,j+1,k+1,n)-q(l,i,j,k,n)
                    dfmlr = q(l,i  ,j-1,k+1,n)-q(l,i,j,k,n)
                    dfmmr = q(l,i  ,j  ,k+1,n)-q(l,i,j,k,n)
                    dfmrr = q(l,i  ,j+1,k+1,n)-q(l,i,j,k,n)
                    dfrlr = q(l,i+1,j-1,k+1,n)-q(l,i,j,k,n)
                    dfrmr = q(l,i+1,j  ,k+1,n)-q(l,i,j,k,n)
                    dfrrr = q(l,i+1,j+1,k+1,n)-q(l,i,j,k,n)

                    vmin = min(dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl, &
                         &     dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm, &
                         &     dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr)
                    vmax = max(dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl, &
                         &     dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm, &
                         &     dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr)

                    dfx  = half*(q(l,i+1,j,k,n)-q(l,i-1,j,k,n))
                    dfy  = half*(q(l,i,j+1,k,n)-q(l,i,j-1,k,n))
                    dfz  = half*(q(l,i,j,k+1,n)-q(l,i,j,k-1,n))
                    dff  = half*(abs(dfx)+abs(dfy)+abs(dfz))

                    if(dff>zero)then
                       slop = min(one,min(abs(vmin),abs(vmax))/dff)
                    else
                       slop = one
                    endif

                    dlim = slop

                    dq(l,i,j,k,n,1) = dlim*dfx
                    dq(l,i,j,k,n,2) = dlim*dfy
                    dq(l,i,j,k,n,3) = dlim*dfz

                 end do
              end do
           end do
        end do
     end do
  else
     write(*,*)'Unknown slope type'
     stop
  endif
#endif

end subroutine uslope



!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine Newton_Raphson_Mignone(D,M,E,gamma,R)
  use amr_parameters
  implicit none
  real(dp)::x,R,v0
  real(dp)::f_Mignone,f_prim_Mignone
  real(dp)::epsilon,Delta
  real(dp)::D,E,Eprim,M,gamma

  !initial guess
  Delta = 16.0d0*E**2-12d0*M**2
  R = (4d0*E+sqrt(Delta))/6d0
  !Switch to prime variables
  R=R-D ; Eprim=E-D
  !NR loop
  epsilon=1.
  do while (abs(epsilon)>1d-10)
     epsilon=f_Mignone(R,D,M,Eprim,gamma)/f_prim_Mignone(R,D,M,Eprim,gamma)/R
     R=R*(1d0-epsilon)
  enddo


! Go back to Q instead of Qprime
R=R+D

return
end subroutine Newton_Raphson_Mignone
!###########################################################
!###########################################################
!###########################################################
!###########################################################
function f_Mignone(R,D,M,Eprim,gamma)
  use amr_parameters
  use hydro_parameters,only:eos
  implicit none
  real(dp)::R,D,M,Eprim,gamma
  real(dp)::f_Mignone
  real(dp)::u2,lor,Xsi,P,rho
  u2=M**2/((R+D)**2-M**2) ; lor=(1+u2)**(0.5)
  Xsi=(R-u2/(lor+1)*D)/lor**2
  if (eos .eq. 'TM') then
     rho=D/lor
     P=(2d0*xsi*(xsi+2d0*rho))/(5d0*(xsi+rho)+sqrt(9d0*xsi**2+18d0*rho*xsi+25d0*rho**2))
  else
     P=(gamma-1.0d0)/gamma*Xsi
  endif
  f_Mignone  = R-P-Eprim
  return
end function f_Mignone
!###########################################################
!###########################################################
!###########################################################
!###########################################################
function f_prim_Mignone(R,D,M,Eprim,gamma)
  use amr_parameters
  use hydro_parameters,only:eos
  implicit none
  real(dp)::R,D,M,Eprim,gamma
  real(dp)::f_prim_Mignone
  real(dp)::u2,lor,dpdR,dpdxsi,rho,xsi,dpdrho,dv2dR,dxsidR,drhodR,P
  u2=M**2/((R+D)**2-M**2) ; lor=(1+u2)**(0.5)

  if (eos .eq. 'TM') then
     Xsi=(R-u2/(lor+1)*D)/lor**2
     rho=D/lor
     P=(2d0*xsi*(xsi+2d0*rho))/(5d0*(xsi+rho)+sqrt(9d0*xsi**2+18d0*rho*xsi+25d0*rho**2))

     dpdxsi=(2d0*xsi+2d0*rho-5d0*P)/(5d0*rho+5d0*xsi-8d0*P)
     dpdrho=(2d0*xsi-5d0*P)/(5d0*rho+5d0*xsi-8d0*P)
     dv2dR=-2d0*M**2/(R+D)**3

     dxsidR=1d0/lor**2-lor/2d0*(D+2d0*lor*xsi)*dv2dR
     drhodR=D*lor/2d0*dv2dR

     dpdR=dpdxsi*dxsidR+dpdrho*drhodR
  else

     dpdR=(gamma-1d0)/gamma
     dpdR=dpdR*(1d0+M**2/(R+D)**2*(1d0-D*lor/(R+D)))
  endif
  f_prim_Mignone=1d0-dpdR
  return
end function f_prim_Mignone
