! ---------------------------------------------------------------
!  UNSPLIT     Unsplit second order Godunov integrator for
!              polytropic gas dynamics using either
!              MUSCL-HANCOCK scheme or Collela's PLMDE scheme
!              with various slope limiters.
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
! ----------------------------------------------------------------
subroutine unsplit(uin,gravin,pin,flux,tmp,dx,dy,dz,dt,ngrid)
  use amr_parameters
  use const
  use hydro_parameters
  implicit none

  integer ::ngrid
  real(dp)::dx,dy,dz,dt
  real(dp)::dtdx

  ! Input states
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::pin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::uin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim)::gravin

  ! Output fluxes
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:nvar,1:ndim)::flux
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:2   ,1:ndim)::tmp

  ! Primitive variables
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar),save::qin

  ! Slopes
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim),save::dq

  ! Left and right state arrays
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim),save::qm
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim),save::qp

  ! Velocity divergence
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2)::divu

  ! Local scalar variables
  integer::i,j,k,l,ivar
  integer::ilo,ihi,jlo,jhi,klo,khi

  ilo=MIN(1,iu1+2); ihi=MAX(1,iu2-2)
  jlo=MIN(1,ju1+2); jhi=MAX(1,ju2-2)
  klo=MIN(1,ku1+2); khi=MAX(1,ku2-2)
  dtdx = dt/dx

  ! Translate to primitive variables, compute sound speeds
  call ctoprim(uin,qin,gravin,dt,ngrid)

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
  endif
  if(scheme=='plmde')then
#if NDIM==1
     call tracex  (qin,dq,qm,qp,dx      ,dt,ngrid)
#endif
#if NDIM==2
     call tracexy (qin,dq,qm,qp,dx,dy   ,dt,ngrid)
#endif
#if NDIM==3
     call tracexyz(qin,dq,qm,qp,dx,dy,dz,dt,ngrid)
#endif
  endif

  ! Solve for 1D flux in X direction
  call cmpflxm(qm,iu1+1,iu2+1,ju1  ,ju2  ,ku1  ,ku2  , &
       &       qp,iu1  ,iu2  ,ju1  ,ju2  ,ku1  ,ku2  , &
       &          if1  ,if2  ,jlo  ,jhi  ,klo  ,khi  , &
       &       2,3,4,flux,tmp,1,dtdx,ngrid)

  ! Solve for 1D flux in Y direction
#if NDIM>1
  call cmpflxm(qm,iu1  ,iu2  ,ju1+1,ju2+1,ku1  ,ku2  , &
       &       qp,iu1  ,iu2  ,ju1  ,ju2  ,ku1  ,ku2  , &
       &          ilo  ,ihi  ,jf1  ,jf2  ,klo  ,khi  , &
       &       3,2,4,flux,tmp,2,dtdx,ngrid)
#endif

  ! Solve for 1D flux in Z direction
#if NDIM>2
  call cmpflxm(qm,iu1  ,iu2  ,ju1  ,ju2  ,ku1+1,ku2+1, &
       &       qp,iu1  ,iu2  ,ju1  ,ju2  ,ku1  ,ku2  , &
       &          ilo  ,ihi  ,jlo  ,jhi  ,kf1  ,kf2  , &
       &       4,2,3,flux,tmp,3,dtdx,ngrid)
#endif

  if(difmag>0.0)then
    call cmpdivu(qin,divu,dx,dy,dz,ngrid)
    call consup(uin,flux,divu,dt,ngrid)
  endif

end subroutine unsplit
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine trace1d(q,dq,qm,qp,dx,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer ::ngrid
  real(dp)::dx, dt

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::dq
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qm
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qp

  ! Local variables
  integer ::i, j, k, l
  integer ::ilo,ihi,jlo,jhi,klo,khi
  integer ::ir, iu, ip
  real(dp)::dtdx
  real(dp)::r, u, p
  real(dp)::drx, dux, dpx
  real(dp)::sr0, su0, sp0
#if NENER>0
  integer::irad
  real(dp),dimension(1:nener)::e, dex, se0
#endif
#if NVAR > NHYDRO + NENER
  integer::n
  real(dp)::a, dax, sa0
#endif

  dtdx = dt/dx

  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)
  ir=1; iu=2; ip=3

  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           do l = 1, ngrid

              ! Cell centered values
              r   =  q(l,i,j,k,ir)
              u   =  q(l,i,j,k,iu)
              p   =  q(l,i,j,k,ip)
#if NENER>0
              do irad=1,nener
                 e(irad) = q(l,i,j,k,ip+irad)
              end do
#endif
              ! TVD slopes in X direction
              drx = dq(l,i,j,k,ir,1)
              dux = dq(l,i,j,k,iu,1)
              dpx = dq(l,i,j,k,ip,1)
#if NENER>0
              do irad=1,nener
                 dex(irad) = dq(l,i,j,k,ip+irad,1)
              end do
#endif

              ! Source terms (including transverse derivatives)
              sr0 = -u*drx - (dux)*r
              sp0 = -u*dpx - (dux)*gamma*p
              su0 = -u*dux - (dpx)/r
#if NENER>0
              do irad=1,nener
                 su0 = su0 - (dex(irad))/r
                 se0(irad) = -u*dex(irad) &
                      & - (dux)*gamma_rad(irad)*e(irad)
              end do
#endif

              ! Right state
              qp(l,i,j,k,ir,1) = r - half*drx + sr0*dtdx*half
              qp(l,i,j,k,iu,1) = u - half*dux + su0*dtdx*half
              qp(l,i,j,k,ip,1) = p - half*dpx + sp0*dtdx*half
!              qp(l,i,j,k,ir,1) = max(smallr, qp(l,i,j,k,ir,1))
              if(qp(l,i,j,k,ir,1)<smallr)qp(l,i,j,k,ir,1)=r
#if NENER>0
              do irad=1,nener
                 qp(l,i,j,k,ip+irad,1) = e(irad) - half*dex(irad) + se0(irad)*dtdx*half
              end do
#endif

              ! Left state
              qm(l,i,j,k,ir,1) = r + half*drx + sr0*dtdx*half
              qm(l,i,j,k,iu,1) = u + half*dux + su0*dtdx*half
              qm(l,i,j,k,ip,1) = p + half*dpx + sp0*dtdx*half
!              qm(l,i,j,k,ir,1) = max(smallr, qm(l,i,j,k,ir,1))
              if(qm(l,i,j,k,ir,1)<smallr)qm(l,i,j,k,ir,1)=r
#if NENER>0
              do irad=1,nener
                 qm(l,i,j,k,ip+irad,1) = e(irad) + half*dex(irad) + se0(irad)*dtdx*half
              end do
#endif

           end do
        end do
     end do
  end do

#if NVAR > NHYDRO + NENER
  ! Passive scalars
  do n = ndim+nener+3, nvar
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              do l = 1, ngrid
                 a   = q(l,i,j,k,n)       ! Cell centered values
                 u   = q(l,i,j,k,iu)
                 dax = dq(l,i,j,k,n,1)    ! TVD slopes
                 sa0 = -u*dax             ! Source terms
                 qp(l,i,j,k,n,1) = a - half*dax + sa0*dtdx*half   ! Right state
                 qm(l,i,j,k,n,1) = a + half*dax + sa0*dtdx*half   ! Left state
              end do
           end do
        end do
     end do
  end do
#endif

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

  integer ::ngrid
  real(dp)::dx, dy, dt

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::dq
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qm
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qp

  ! declare local variables
  integer ::i, j, k, l
  integer ::ilo,ihi,jlo,jhi,klo,khi
  integer ::ir, iu, iv, ip
  real(dp)::dtdx, dtdy
  real(dp)::r, u, v, p
  real(dp)::drx, dux, dvx, dpx
  real(dp)::dry, duy, dvy, dpy
  real(dp)::sr0, su0, sv0, sp0
#if NENER>0
  integer ::irad
  real(dp),dimension(1:nener)::e, dex, dey, se0
#endif
#if NVAR > NHYDRO + NENER
  integer ::n
  real(dp)::a, dax, day, sa0
#endif

  dtdx = dt/dx
  dtdy = dt/dy
  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)
  ir=1; iu=2; iv=3; ip=4

  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           do l = 1, ngrid

              ! Cell centered values
              r   =  q(l,i,j,k,ir)
              u   =  q(l,i,j,k,iu)
              v   =  q(l,i,j,k,iv)
              p   =  q(l,i,j,k,ip)
#if NENER>0
              do irad=1,nener
                 e(irad) = q(l,i,j,k,ip+irad)
              end do
#endif

              ! TVD slopes in all directions
              drx = dq(l,i,j,k,ir,1)
              dux = dq(l,i,j,k,iu,1)
              dvx = dq(l,i,j,k,iv,1)
              dpx = dq(l,i,j,k,ip,1)
#if NENER>0
              do irad=1,nener
                 dex(irad) = dq(l,i,j,k,ip+irad,1)
              end do
#endif

              dry = dq(l,i,j,k,ir,2)
              duy = dq(l,i,j,k,iu,2)
              dvy = dq(l,i,j,k,iv,2)
              dpy = dq(l,i,j,k,ip,2)
#if NENER>0
              do irad=1,nener
                 dey(irad) = dq(l,i,j,k,ip+irad,2)
              end do
#endif

              ! source terms (with transverse derivatives)
              sr0 = -u*drx-v*dry - (dux+dvy)*r
              sp0 = -u*dpx-v*dpy - (dux+dvy)*gamma*p
              su0 = -u*dux-v*duy - (dpx    )/r
              sv0 = -u*dvx-v*dvy - (dpy    )/r
#if NENER>0
              do irad=1,nener
                 su0 = su0 - (dex(irad))/r
                 sv0 = sv0 - (dey(irad))/r
                 se0(irad) = -u*dex(irad)-v*dey(irad) &
                      & - (dux+dvy)*gamma_rad(irad)*e(irad)
              end do
#endif

              ! Right state at left interface
              qp(l,i,j,k,ir,1) = r - half*drx + sr0*dtdx*half
              qp(l,i,j,k,iu,1) = u - half*dux + su0*dtdx*half
              qp(l,i,j,k,iv,1) = v - half*dvx + sv0*dtdx*half
              qp(l,i,j,k,ip,1) = p - half*dpx + sp0*dtdx*half
!              qp(l,i,j,k,ir,1) = max(smallr, qp(l,i,j,k,ir,1))
              if(qp(l,i,j,k,ir,1)<smallr)qp(l,i,j,k,ir,1)=r
#if NENER>0
              do irad=1,nener
                 qp(l,i,j,k,ip+irad,1) = e(irad) - half*dex(irad) + se0(irad)*dtdx*half
              end do
#endif

              ! Left state at right interface
              qm(l,i,j,k,ir,1) = r + half*drx + sr0*dtdx*half
              qm(l,i,j,k,iu,1) = u + half*dux + su0*dtdx*half
              qm(l,i,j,k,iv,1) = v + half*dvx + sv0*dtdx*half
              qm(l,i,j,k,ip,1) = p + half*dpx + sp0*dtdx*half
!              qm(l,i,j,k,ir,1) = max(smallr, qm(l,i,j,k,ir,1))
              if(qm(l,i,j,k,ir,1)<smallr)qm(l,i,j,k,ir,1)=r
#if NENER>0
              do irad=1,nener
                 qm(l,i,j,k,ip+irad,1) = e(irad) + half*dex(irad) + se0(irad)*dtdx*half
              end do
#endif

              ! Top state at bottom interface
              qp(l,i,j,k,ir,2) = r - half*dry + sr0*dtdy*half
              qp(l,i,j,k,iu,2) = u - half*duy + su0*dtdy*half
              qp(l,i,j,k,iv,2) = v - half*dvy + sv0*dtdy*half
              qp(l,i,j,k,ip,2) = p - half*dpy + sp0*dtdy*half
!              qp(l,i,j,k,ir,2) = max(smallr, qp(l,i,j,k,ir,2))
              if(qp(l,i,j,k,ir,2)<smallr)qp(l,i,j,k,ir,2)=r
#if NENER>0
              do irad=1,nener
                 qp(l,i,j,k,ip+irad,2) = e(irad) - half*dey(irad) + se0(irad)*dtdy*half
              end do
#endif

              ! Bottom state at top interface
              qm(l,i,j,k,ir,2) = r + half*dry + sr0*dtdy*half
              qm(l,i,j,k,iu,2) = u + half*duy + su0*dtdy*half
              qm(l,i,j,k,iv,2) = v + half*dvy + sv0*dtdy*half
              qm(l,i,j,k,ip,2) = p + half*dpy + sp0*dtdy*half
!              qm(l,i,j,k,ir,2) = max(smallr, qm(l,i,j,k,ir,2))
              if(qm(l,i,j,k,ir,2)<smallr)qm(l,i,j,k,ir,2)=r
#if NENER>0
              do irad=1,nener
                 qm(l,i,j,k,ip+irad,2) = e(irad) + half*dey(irad) + se0(irad)*dtdy*half
              end do
#endif

           end do
        end do
     end do
  end do

#if NVAR > NHYDRO + NENER
  ! passive scalars
  do n = ndim+nener+3, nvar
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              do l = 1, ngrid
                 a   = q(l,i,j,k,n)       ! Cell centered values
                 u   = q(l,i,j,k,iu)
                 v   = q(l,i,j,k,iv)
                 dax = dq(l,i,j,k,n,1)    ! TVD slopes
                 day = dq(l,i,j,k,n,2)
                 sa0 = -u*dax-v*day       ! Source terms
                 qp(l,i,j,k,n,1) = a - half*dax + sa0*dtdx*half   ! Right state
                 qm(l,i,j,k,n,1) = a + half*dax + sa0*dtdx*half   ! Left state
                 qp(l,i,j,k,n,2) = a - half*day + sa0*dtdy*half   ! Top state
                 qm(l,i,j,k,n,2) = a + half*day + sa0*dtdy*half   ! Bottom state
              end do
           end do
        end do
     end do
  end do
#endif

end subroutine trace2d
#endif
!###########################################################
!###########################################################
!###########################################################
!###########################################################
#if NDIM>2
subroutine trace3d(q,dq,qm,qp,dx,dy,dz,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer ::ngrid
  real(dp)::dx, dy, dz, dt

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::dq
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qm
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qp

  ! declare local variables
  integer ::i, j, k, l
  integer ::ilo,ihi,jlo,jhi,klo,khi
  integer ::ir, iu, iv, iw, ip
  real(dp)::dtdx, dtdy, dtdz
  real(dp)::r, u, v, w, p
  real(dp)::drx, dux, dvx, dwx, dpx
  real(dp)::dry, duy, dvy, dwy, dpy
  real(dp)::drz, duz, dvz, dwz, dpz
  real(dp)::sr0, su0, sv0, sw0, sp0
#if NENER>0
  integer ::irad
  real(dp),dimension(1:nener)::e, dex, dey, dez, se0
#endif
#if NVAR > NHYDRO + NENER
  integer ::n
  real(dp)::a, dax, day, daz, sa0
#endif

  dtdx = dt/dx
  dtdy = dt/dy
  dtdz = dt/dz
  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)
  ir=1; iu=2; iv=3; iw=4; ip=5

  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           do l = 1, ngrid

              ! Cell centered values
              r   =  q(l,i,j,k,ir)
              u   =  q(l,i,j,k,iu)
              v   =  q(l,i,j,k,iv)
              w   =  q(l,i,j,k,iw)
              p   =  q(l,i,j,k,ip)
#if NENER>0
              do irad=1,nener
                 e(irad) = q(l,i,j,k,ip+irad)
              end do
#endif

              ! TVD slopes in all 3 directions
              drx = dq(l,i,j,k,ir,1)
              dpx = dq(l,i,j,k,ip,1)
              dux = dq(l,i,j,k,iu,1)
              dvx = dq(l,i,j,k,iv,1)
              dwx = dq(l,i,j,k,iw,1)
#if NENER>0
              do irad=1,nener
                 dex(irad) = dq(l,i,j,k,ip+irad,1)
              end do
#endif

              dry = dq(l,i,j,k,ir,2)
              dpy = dq(l,i,j,k,ip,2)
              duy = dq(l,i,j,k,iu,2)
              dvy = dq(l,i,j,k,iv,2)
              dwy = dq(l,i,j,k,iw,2)
#if NENER>0
              do irad=1,nener
                 dey(irad) = dq(l,i,j,k,ip+irad,2)
              end do
#endif

              drz = dq(l,i,j,k,ir,3)
              dpz = dq(l,i,j,k,ip,3)
              duz = dq(l,i,j,k,iu,3)
              dvz = dq(l,i,j,k,iv,3)
              dwz = dq(l,i,j,k,iw,3)
#if NENER>0
              do irad=1,nener
                 dez(irad) = dq(l,i,j,k,ip+irad,3)
              end do
#endif

              ! Source terms (including transverse derivatives)
              sr0 = -u*drx-v*dry-w*drz - (dux+dvy+dwz)*r
              sp0 = -u*dpx-v*dpy-w*dpz - (dux+dvy+dwz)*gamma*p
              su0 = -u*dux-v*duy-w*duz - (dpx        )/r
              sv0 = -u*dvx-v*dvy-w*dvz - (dpy        )/r
              sw0 = -u*dwx-v*dwy-w*dwz - (dpz        )/r
#if NENER>0
              do irad=1,nener
                 su0 = su0 - (dex(irad))/r
                 sv0 = sv0 - (dey(irad))/r
                 sw0 = sw0 - (dez(irad))/r
                 se0(irad) = -u*dex(irad)-v*dey(irad)-w*dez(irad) &
                      & - (dux+dvy+dwz)*gamma_rad(irad)*e(irad)
              end do
#endif

              ! Right state at left interface
              qp(l,i,j,k,ir,1) = r - half*drx + sr0*dtdx*half
              qp(l,i,j,k,ip,1) = p - half*dpx + sp0*dtdx*half
              qp(l,i,j,k,iu,1) = u - half*dux + su0*dtdx*half
              qp(l,i,j,k,iv,1) = v - half*dvx + sv0*dtdx*half
              qp(l,i,j,k,iw,1) = w - half*dwx + sw0*dtdx*half
!              qp(l,i,j,k,ir,1) = max(smallr, qp(l,i,j,k,ir,1))
              if(qp(l,i,j,k,ir,1)<smallr)qp(l,i,j,k,ir,1)=r
#if NENER>0
              do irad=1,nener
                 qp(l,i,j,k,ip+irad,1) = e(irad) - half*dex(irad) + se0(irad)*dtdx*half
              end do
#endif

              ! Left state at left interface
              qm(l,i,j,k,ir,1) = r + half*drx + sr0*dtdx*half
              qm(l,i,j,k,ip,1) = p + half*dpx + sp0*dtdx*half
              qm(l,i,j,k,iu,1) = u + half*dux + su0*dtdx*half
              qm(l,i,j,k,iv,1) = v + half*dvx + sv0*dtdx*half
              qm(l,i,j,k,iw,1) = w + half*dwx + sw0*dtdx*half
!              qm(l,i,j,k,ir,1) = max(smallr, qm(l,i,j,k,ir,1))
              if(qm(l,i,j,k,ir,1)<smallr)qm(l,i,j,k,ir,1)=r
#if NENER>0
              do irad=1,nener
                 qm(l,i,j,k,ip+irad,1) = e(irad) + half*dex(irad) + se0(irad)*dtdx*half
              end do
#endif

              ! Top state at bottom interface
              qp(l,i,j,k,ir,2) = r - half*dry + sr0*dtdy*half
              qp(l,i,j,k,ip,2) = p - half*dpy + sp0*dtdy*half
              qp(l,i,j,k,iu,2) = u - half*duy + su0*dtdy*half
              qp(l,i,j,k,iv,2) = v - half*dvy + sv0*dtdy*half
              qp(l,i,j,k,iw,2) = w - half*dwy + sw0*dtdy*half
!              qp(l,i,j,k,ir,2) = max(smallr, qp(l,i,j,k,ir,2))
              if(qp(l,i,j,k,ir,2)<smallr)qp(l,i,j,k,ir,2)=r
#if NENER>0
              do irad=1,nener
                 qp(l,i,j,k,ip+irad,2) = e(irad) - half*dey(irad) + se0(irad)*dtdy*half
              end do
#endif

              ! Bottom state at top interface
              qm(l,i,j,k,ir,2) = r + half*dry + sr0*dtdy*half
              qm(l,i,j,k,ip,2) = p + half*dpy + sp0*dtdy*half
              qm(l,i,j,k,iu,2) = u + half*duy + su0*dtdy*half
              qm(l,i,j,k,iv,2) = v + half*dvy + sv0*dtdy*half
              qm(l,i,j,k,iw,2) = w + half*dwy + sw0*dtdy*half
!              qm(l,i,j,k,ir,2) = max(smallr, qm(l,i,j,k,ir,2))
              if(qm(l,i,j,k,ir,2)<smallr)qm(l,i,j,k,ir,2)=r
#if NENER>0
              do irad=1,nener
                 qm(l,i,j,k,ip+irad,2) = e(irad) + half*dey(irad) + se0(irad)*dtdy*half
              end do
#endif

              ! Back state at front interface
              qp(l,i,j,k,ir,3) = r - half*drz + sr0*dtdz*half
              qp(l,i,j,k,ip,3) = p - half*dpz + sp0*dtdz*half
              qp(l,i,j,k,iu,3) = u - half*duz + su0*dtdz*half
              qp(l,i,j,k,iv,3) = v - half*dvz + sv0*dtdz*half
              qp(l,i,j,k,iw,3) = w - half*dwz + sw0*dtdz*half
!              qp(l,i,j,k,ir,3) = max(smallr, qp(l,i,j,k,ir,3))
              if(qp(l,i,j,k,ir,3)<smallr)qp(l,i,j,k,ir,3)=r
#if NENER>0
              do irad=1,nener
                 qp(l,i,j,k,ip+irad,3) = e(irad) - half*dez(irad) + se0(irad)*dtdz*half
              end do
#endif

              ! Front state at back interface
              qm(l,i,j,k,ir,3) = r + half*drz + sr0*dtdz*half
              qm(l,i,j,k,ip,3) = p + half*dpz + sp0*dtdz*half
              qm(l,i,j,k,iu,3) = u + half*duz + su0*dtdz*half
              qm(l,i,j,k,iv,3) = v + half*dvz + sv0*dtdz*half
              qm(l,i,j,k,iw,3) = w + half*dwz + sw0*dtdz*half
!              qm(l,i,j,k,ir,3) = max(smallr, qm(l,i,j,k,ir,3))
              if(qm(l,i,j,k,ir,3)<smallr)qm(l,i,j,k,ir,3)=r
#if NENER>0
              do irad=1,nener
                 qm(l,i,j,k,ip+irad,3) = e(irad) + half*dez(irad) + se0(irad)*dtdz*half
              end do
#endif

           end do
        end do
     end do
  end do

#if NVAR > NHYDRO + NENER
  ! Passive scalars
  do n = ndim+nener+3, nvar
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              do l = 1, ngrid
                 a   = q(l,i,j,k,n)       ! Cell centered values
                 u   = q(l,i,j,k,iu)
                 v   = q(l,i,j,k,iv)
                 w   = q(l,i,j,k,iw)
                 dax = dq(l,i,j,k,n,1)    ! TVD slopes
                 day = dq(l,i,j,k,n,2)
                 daz = dq(l,i,j,k,n,3)
                 sa0 = -u*dax-v*day-w*daz     ! Source terms
                 qp(l,i,j,k,n,1) = a - half*dax + sa0*dtdx*half  ! Right state
                 qm(l,i,j,k,n,1) = a + half*dax + sa0*dtdx*half  ! Left state
                 qp(l,i,j,k,n,2) = a - half*day + sa0*dtdy*half  ! Bottom state
                 qm(l,i,j,k,n,2) = a + half*day + sa0*dtdy*half  ! Upper state
                 qp(l,i,j,k,n,3) = a - half*daz + sa0*dtdz*half  ! Front state
                 qm(l,i,j,k,n,3) = a + half*daz + sa0*dtdz*half  ! Back state
              end do
           end do
        end do
     end do
  end do
#endif

end subroutine trace3d
#endif
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmpflxm(qm,im1,im2,jm1,jm2,km1,km2, &
     &             qp,ip1,ip2,jp1,jp2,kp1,kp2, &
     &                ilo,ihi,jlo,jhi,klo,khi, ln,lt1,lt2, &
     &             flux,tmp,idim,dtdx,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  real(dp)::dtdx
  integer ::idim,ngrid
  integer ::ln,lt1,lt2
  integer ::im1,im2,jm1,jm2,km1,km2
  integer ::ip1,ip2,jp1,jp2,kp1,kp2
  integer ::ilo,ihi,jlo,jhi,klo,khi
  real(dp),dimension(1:nvector,im1:im2,jm1:jm2,km1:km2,1:nvar,1:ndim)::qm
  real(dp),dimension(1:nvector,ip1:ip2,jp1:jp2,kp1:kp2,1:nvar,1:ndim)::qp
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:nvar,1:ndim)::flux
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:2,   1:ndim)::tmp

  ! local variables
  integer ::i, j, k, l, xdim
  real(dp)::entho
  real(dp),dimension(1:nvector,1:nvar),save::qleft,qright
  real(dp),dimension(1:nvector,1:nvar+1),save::fgdnv
#if NVAR > NHYDRO
  integer ::n
#endif

  entho=one/(gamma-one)
  xdim=ln-1

  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi

           ! Mass density
           do l = 1, ngrid
              qleft (l,1) = qm(l,i,j,k,1,xdim)
              qright(l,1) = qp(l,i,j,k,1,xdim)
           end do

           ! Normal velocity
           do l = 1, ngrid
              qleft (l,2) = qm(l,i,j,k,ln,xdim)
              qright(l,2) = qp(l,i,j,k,ln,xdim)
           end do

           ! Pressure
           do l = 1, ngrid
              qleft (l,3) = qm(l,i,j,k,neul,xdim)
              qright(l,3) = qp(l,i,j,k,neul,xdim)
           end do

           ! Tangential velocity 1
#if NDIM>1
           do l = 1, ngrid
              qleft (l,4) = qm(l,i,j,k,lt1,xdim)
              qright(l,4) = qp(l,i,j,k,lt1,xdim)
           end do
#endif
           ! Tangential velocity 2
#if NDIM>2
           do l = 1, ngrid
              qleft (l,5) = qm(l,i,j,k,lt2,xdim)
              qright(l,5) = qp(l,i,j,k,lt2,xdim)
           end do
#endif
#if NVAR > NHYDRO
           ! Other advected quantities
           do n = nhydro+1, nvar
              do l = 1, ngrid
                 qleft (l,n) = qm(l,i,j,k,n,xdim)
                 qright(l,n) = qp(l,i,j,k,n,xdim)
              end do
           end do
#endif
           ! Solve Riemann problem
           if(riemann.eq.'acoustic')then
              call riemann_acoustic(qleft,qright,fgdnv,ngrid)
           else if (riemann.eq.'exact')then
              call riemann_approx  (qleft,qright,fgdnv,ngrid)
           else if (riemann.eq.'llf')then
              call riemann_llf     (qleft,qright,fgdnv,ngrid)
           else if (riemann.eq.'hllc')then
              call riemann_hllc    (qleft,qright,fgdnv,ngrid)
           else if (riemann.eq.'hll')then
              call riemann_hll     (qleft,qright,fgdnv,ngrid)
           else
              write(*,*)'unknown Riemann solver'
              stop
           end if

           ! Compute fluxes

           ! Mass density
           do l = 1, ngrid
              flux(l,i,j,k,1,idim) = fgdnv(l,1) * dtdx
           end do

           ! Normal momentum
           do l = 1, ngrid
              flux(l,i,j,k,ln,idim) = fgdnv(l,2) * dtdx
           end do

           ! Transverse momentum 1
#if NDIM>1
           do l = 1, ngrid
              flux(l,i,j,k,lt1,idim) = fgdnv(l,4) * dtdx
           end do
#endif
           ! Transverse momentum 2
#if NDIM>2
           do l = 1, ngrid
              flux(l,i,j,k,lt2,idim) = fgdnv(l,5) * dtdx
           end do
#endif
           ! Total energy
           do l = 1, ngrid
              flux(l,i,j,k,neul,idim) = fgdnv(l,3) * dtdx
           end do

#if NVAR > NHYDRO
           ! Other advected quantities
           do n = nhydro+1, nvar
              do l = 1, ngrid
                 flux(l,i,j,k,n,idim) = fgdnv(l,n) * dtdx
              end do
           end do
#endif
           ! Normal velocity
           do l = 1, ngrid
              tmp(l,i,j,k,1,idim) = half*(qleft(l,2)+qright(l,2)) * dtdx
           end do
           ! Internal energy flux
           do l = 1,ngrid
              tmp(l,i,j,k,2,idim) = fgdnv(l,nvar+1) * dtdx
           end do

        end do
     end do
  end do

end subroutine cmpflxm
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine ctoprim(uin,q,gravin,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none
  integer, intent(in)::ngrid
  real(dp),intent(in)::dt
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar),intent(in)::uin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim),intent(in)::gravin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar),intent(out)::q
  !---------------------------------------------------------
  ! Translate Conservative variables uin to PRIMitive variables q.
  !---------------------------------------------------------
  integer ::i, j, k, l
  real(dp) ::rho_grid,vx_grid,vy_grid,vz_grid
  real(dp)::eint, smalle, dtxhalf, oneoverrho
  real(dp)::eken, erad
#if NVAR > NHYDRO + NENER
  integer ::n
#endif
#if NENER>0
  integer ::irad
#endif

  smalle = smallc**2/gamma/(gamma-one)
  dtxhalf = dt*half

  ! Convert to primitive variable
  do k = ku1, ku2
     do j = ju1, ju2
        do i = iu1, iu2
           do l = 1, ngrid

              ! Compute density
              ! We keep it in a local scalar to avoid store-load dependencies of q(l,i,j,k,1)
              rho_grid = max(uin(l,i,j,k,1),smallr)

              ! Compute velocities
              oneoverrho = one/rho_grid
              vx_grid = uin(l,i,j,k,2)*oneoverrho
#if NDIM>1
              vy_grid = uin(l,i,j,k,3)*oneoverrho
#endif
#if NDIM>2
              vz_grid = uin(l,i,j,k,4)*oneoverrho
#endif

              ! Compute specific kinetic energy
              eken = half*vx_grid*vx_grid
#if NDIM>1
              eken = eken + half*vy_grid*vy_grid
#endif
#if NDIM>2
              eken = eken + half*vz_grid*vz_grid
#endif
              ! Compute non-thermal pressure
              erad = zero
#if NENER>0
              do irad = 1,nener
                 q(l,i,j,k,nhydro+irad) = (gamma_rad(irad)-one)*uin(l,i,j,k,nhydro+irad)
                 erad = erad+uin(l,i,j,k,nhydro+irad)*oneoverrho
              enddo
#endif
              ! Compute thermal pressure
              eint = MAX(uin(l,i,j,k,neul)*oneoverrho-eken-erad,smalle)
              q(l,i,j,k,neul) = (gamma-one)*rho_grid*eint

              ! Now, we store the density
              q(l,i,j,k,1) = rho_grid

              ! Store velocity and apply gravity predictor step
              q(l,i,j,k,2) = vx_grid + gravin(l,i,j,k,1)*dtxhalf
#if NDIM>1
              q(l,i,j,k,3) = vy_grid + gravin(l,i,j,k,2)*dtxhalf
#endif
#if NDIM>2
              q(l,i,j,k,4) = vz_grid + gravin(l,i,j,k,3)*dtxhalf
#endif

           end do
        end do
     end do
  end do

#if NVAR > NHYDRO + NENER
  ! Passive scalar
  do n = ndim+nener+3, nvar
     do k = ku1, ku2
        do j = ju1, ju2
           do i = iu1, iu2
              do l = 1, ngrid
                 oneoverrho = one/q(l,i,j,k,1)
                 q(l,i,j,k,n) = uin(l,i,j,k,n)*oneoverrho
              end do
           end do
        end do
     end do
  end do
#endif

end subroutine ctoprim
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine uslope(q,dq,dx,dt,ngrid)
  use amr_parameters, only:dp,nvector,ndim
  use hydro_parameters, only:nvar,slope_type,iu1,iu2,ju1,ju2,ku1,ku2
  use const
  use slope_types
  implicit none

  integer,intent(in)::ngrid
  real(dp),intent(in)::dx,dt
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar),intent(in)::q
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim),intent(out)::dq

  ! local arrays
  integer::i, j, k, l, n
  real(dp)::slope_type_real,dtdx
  integer::ilo,ihi,jlo,jhi,klo,khi

  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)

  slope_type_real = REAL(slope_type, kind=dp)
  dtdx=dt/dx

  do n = 1, nvar
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              if(slope_type==0)then
                 dq(:,i,j,k,n,:) = zero
#if NDIM==1
              else if(slope_type==1.or.slope_type==2.or.slope_type==3)then  ! minmod or average
#elif NDIM==2
              else if(slope_type==1.or.slope_type==2)then  ! minmod or average
#else
              else if(slope_type==1)then ! minmod
#endif
                 call calc_uslope_minmod_average(q,dq,i,j,k,n,ngrid,slope_type_real)
#if NDIM==3
              else if(slope_type==2)then
                 ! moncen
                 call calc_uslope_moncen(q,dq,i,j,k,n,ngrid)
#endif
#if NDIM>1
              else if(slope_type==3)then
                 ! positivity preserving unsplit slope (2D or 3D)
                 call calc_uslope_positivity_preserving(q,dq,i,j,k,n,ngrid)
#endif
#if NDIM==1
              else if(slope_type==4)then
                 ! superbee (only 1D)
                 call calc_uslope_superbee(q,dq,i,j,k,n,ngrid,dtdx)

              else if(slope_type==5)then
                 ! ultrabee (only 1D)
                 if(n==1)then
                    call calc_uslope_ultrabee(q,dq,i,j,k,n,ngrid,dtdx)
                 else
                    dq(:,i,j,k,n,:) = zero
                 endif

              else if(slope_type==6)then
                 ! unstable (only 1D)
                 if(n==1)then
                    call calc_uslope_unstable(q,dq,i,j,k,n,ngrid)
                 else
                    dq(:,i,j,k,n,:) = zero
                 endif
#endif
              else if(slope_type==7)then
                 ! van Leer
                 call calc_uslope_vanLeer(q,dq,i,j,k,n,ngrid)

              else if(slope_type==8)then
                 ! generalized moncen/minmod parameterisation (van Leer 1979)
                 call calc_uslope_vanLeer_bis(q,dq,i,j,k,n,ngrid)

              else
                 write(*,*)'Unknown slope type',dx,dt
                 call clean_stop
              endif

           end do
        end do
     end do
  end do


end subroutine uslope
