! ---------------------------------------------------------------
!  UNSPLIT     Unsplit second order Godunov integrator for
!              polytropic gas dynamics using either
!              MUSCL-HANCOCK scheme or Collela''s PLMDE scheme
!              with various slope limiters.
!
!  inputs/outputs
!  uin         => (const)  input state
!  gravin      => (const)  input gravitational acceleration
!  iu10,iu20     => (const)  first and last index of input array,
!  ju10,ju20     => (const)  cell centered,    
!  ku10,ku20     => (const)  including buffer cells.
!  flux       <=  (modify) return fluxes in the 3 coord directions
!  if1,if2     => (const)  first and last index of output array,
!  jf1,jf2     => (const)  edge centered,
!  kf1,kf2     => (const)  for active cells only.
!  dx,dy,dz    => (const)  (dx,dy,dz)
!  dt          => (const)  time step
!  ngrid       => (const)  number of sub-grids
!  ndim        => (const)  number of dimensions
! ----------------------------------------------------------------
subroutine unsplit_gpu_2d(uin,gravin,flux,tmp,dx,nxp,dt)
  use amr_parameters
  use const             
  use hydro_parameters
  use acc_parameters
  implicit none 

  integer ::nxp
  real(dp)::dx,dy,dz,dt

  ! Input states
  real(dp),dimension(-1:nxp+2,-1:nxp+2,-1:nxp+2,1:nvar)::uin 
  real(dp),dimension(-1:nxp+2,-1:nxp+2,-1:nxp+2,1:ndim)::gravin 

  ! Output fluxes
  real(dp),dimension(1:nxp+1,1:nxp+1,1:nxp+1,1:nvar,1:ndim)::flux
  real(dp),dimension(1:nxp+1,1:nxp+1,1:nxp+1,1:2   ,1:ndim)::tmp 

  ! Primitive variables
  real(dp),dimension(-1:nxp+2,-1:nxp+2,-1:nxp+2,1:nvar)::qin 
  real(dp),dimension(-1:nxp+2,-1:nxp+2,-1:nxp+2       )::cin

  ! Slopes
  real(dp),dimension(-1:nxp+2,-1:nxp+2,-1:nxp+2,1:nvar,1:ndim)::dq

  ! Left and right state arrays
  real(dp),dimension(-1:nxp+2,-1:nxp+2,-1:nxp+2,1:nvar,1:ndim)::qm
  real(dp),dimension(-1:nxp+2,-1:nxp+2,-1:nxp+2,1:nvar,1:ndim)::qp
  
  ! Intermediate fluxes
  real(dp),dimension(-1:nxp+2,-1:nxp+2,-1:nxp+2,1:nvar)::fx
  real(dp),dimension(-1:nxp+2,-1:nxp+2,-1:nxp+2,1:2   )::tx

!!!!!! REMEMBER DIVU !!!!!!!!!!!
  ! Velocity divergence
  !real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2)::divu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Local scalar variables
  integer::i,j,k,l,ivar
  integer::ilo,ihi,jlo,jhi,klo,khi

! these variable/settings are needed for compatibility with the non GPU code
  iu10=-1
  iu20=nxp+2
  ju10=-1
  ju20=nxp+2
  ku10=-1
  ku20=nxp+2
  if10=1
  if20=nxp+1
  jf10=1
  jf20=nxp+1
  kf10=1
  kf20=nxp+1

  dy=dx
  dz=dx
!!!!!!!!!!!!!!!!!!!!!!

  ilo=MIN(1,iu10+2); ihi=MAX(1,iu20-2)
  jlo=MIN(1,ju10+2); jhi=MAX(1,ju20-2)
  klo=MIN(1,ku10+2); khi=MAX(1,ku20-2)

#define NTPB 128

!$acc data pcopyin(uin,gravin) create(qin,dq,cin,qm,qp,fx,tx) &
!$acc& pcopyout(flux,tmp)

  ! Translate to primative variables, compute sound speeds  
  call ctoprim_gpu_2d(uin,qin,cin,gravin,dt,nxp)

  ! Compute TVD slopes
  call uslope_gpu_2d(qin,dq,dx,dt,nxp)

  ! Compute 3D traced-states in all three directions
  if(scheme=='muscl')then
#if NDIM==1
     call trace1d_gpu_2d(qin,dq,qm,qp,dx      ,dt,nxp)
#endif
#if NDIM==2
     call trace2d_gpu_2d(qin,dq,qm,qp,dx,dy   ,dt,nxp)
#endif
#if NDIM==3
     call trace3d_gpu_2d(qin,dq,qm,qp,dx,dy,dz,dt,nxp)
#endif
  endif
  if(scheme=='plmde')then
#if NDIM==1
     call tracex_gpu_2d(qin,dq,cin,qm,qp,dx      ,dt,nxp)
#endif
#if NDIM==2
     call tracexy_gpu_2d(qin,dq,cin,qm,qp,dx,dy   ,dt,nxp)
#endif
#if NDIM==3
     call tracexyz_gpu_2d(qin,dq,cin,qm,qp,dx,dy,dz,dt,nxp)
#endif
  endif

  ! Solve for 1D flux in X direction
  fx=0.0
  call cmpflxm_gpu_2d(qm,iu10+1,iu20+1,ju10  ,ju20  ,ku10  ,ku20  , &
       &       qp,iu10  ,iu20  ,ju10  ,ju20  ,ku10  ,ku20  , &
       &          if10  ,if20  ,jlo  ,jhi  ,klo  ,khi  , 2,3,4,fx,tx,nxp)
  ! Save flux in output array
!$acc parallel loop collapse(3) gang worker vector vector_length(NTPB)
  do i=if10,if20
  do j=jlo,jhi
  do k=klo,khi
     do ivar=1,nvar
           flux(i,j,k,ivar,1)=fx(i,j,k,ivar)*dt/dx
     end do
     do ivar=1,2
           tmp (i,j,k,ivar,1)=tx(i,j,k,ivar)*dt/dx
     end do
  end do
  end do
  end do
!$acc end parallel loop

  ! Solve for 1D flux in Y direction
#if NDIM>1
  call cmpflxm_gpu_2d(qm,iu10  ,iu20  ,ju10+1,ju20+1,ku10  ,ku20  , &
       &       qp,iu10  ,iu20  ,ju10  ,ju20  ,ku10  ,ku20  , &
       &          ilo  ,ihi  ,jf10  ,jf20  ,klo  ,khi  , 3,2,4,fx,tx,nxp)
  ! Save flux in output array
!$acc parallel loop collapse(3) gang worker vector vector_length(NTPB)
  do i=ilo,ihi
  do j=jf10,jf20
  do k=klo,khi
     do ivar=1,nvar
           flux(i,j,k,ivar,2)=fx(i,j,k,ivar)*dt/dy
     end do
     do ivar=1,2
           tmp (i,j,k,ivar,2)=tx(i,j,k,ivar)*dt/dy
     end do
  end do
  end do
  end do
!$acc end parallel loop
#endif

  ! Solve for 1D flux in Z direction
#if NDIM>2
  call cmpflxm_gpu_2d(qm,iu10  ,iu20  ,ju10  ,ju20  ,ku10+1,ku20+1, &
       &       qp,iu10  ,iu20  ,ju10  ,ju20  ,ku10  ,ku20  , &
       &          ilo  ,ihi  ,jlo  ,jhi  ,kf10  ,kf20  , 4,2,3,fx,tx,nxp)
  ! Save flux in output array
!$acc parallel loop collapse(3) gang worker vector vector_length(NTPB)
  do i=ilo,ihi
  do j=jlo,jhi
  do k=kf10,kf20
     do ivar=1,nvar
        flux(i,j,k,ivar,3)=fx(i,j,k,ivar)*dt/dz
     end do
     do ivar=1,2
        tmp (i,j,k,ivar,3)=tx(i,j,k,ivar)*dt/dz
     end do
  end do
  end do
  end do
!$acc end parallel loop
#endif
WRITE(*,*)"<>>>>>>>>>>><<<<<><><><><><><><><><><"

!$acc end data
WRITE(*,*)"ooooooooooooooooooooooooooooooooooooo"

end subroutine unsplit_gpu_2d
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine trace1d_gpu_2d(q,dq,qm,qp,dx,dt,nxp)
  use amr_parameters
  use hydro_parameters
  use acc_parameters
  use const
  implicit none

  integer ::nxp
  real(dp)::dx, dt

  real(dp),dimension(-1:nxp+2,-1:nxp+2,-1:nxp+2,1:nvar)::q  
  real(dp),dimension(-1:nxp+2,-1:nxp+2,-1:nxp+2,1:nvar,1:ndim)::dq 
  real(dp),dimension(-1:nxp+2,-1:nxp+2,-1:nxp+2,1:nvar,1:ndim)::qm 
  real(dp),dimension(-1:nxp+2,-1:nxp+2,-1:nxp+2,1:nvar,1:ndim)::qp 

  ! Local variables
  integer ::i, j, k, l, n
  integer ::ilo,ihi,jlo,jhi,klo,khi
  integer ::ir, iu, ip
  real(dp)::dtdx
  real(dp)::r, u, p, a
  real(dp)::drx, dux, dpx, dax
  real(dp)::sr0, su0, sp0, sa0
  
  dtdx = dt/dx

  ilo=MIN(1,iu10+1); ihi=MAX(1,iu20-1)
  jlo=MIN(1,ju10+1); jhi=MAX(1,ju20-1)
  klo=MIN(1,ku10+1); khi=MAX(1,ku20-1)
  ir=1; iu=2; ip=3

  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi

              ! Cell centered values
              r   =  q(i,j,k,ir)
              u   =  q(i,j,k,iu)
              p   =  q(i,j,k,ip)

              ! TVD slopes in X direction
              drx = dq(i,j,k,ir,1)
              dux = dq(i,j,k,iu,1)
              dpx = dq(i,j,k,ip,1)
              
              ! Source terms (including transverse derivatives)
              sr0 = -u*drx - (dux)*r
              su0 = -u*dux - (dpx)/r
              sp0 = -u*dpx - (dux)*gamma*p

              ! Right state
              qp(i,j,k,ir,1) = r - half*drx + sr0*dtdx*half
              qp(i,j,k,iu,1) = u - half*dux + su0*dtdx*half
              qp(i,j,k,ip,1) = p - half*dpx + sp0*dtdx*half
              qp(i,j,k,ir,1) = max(smallr, qp(i,j,k,ir,1))

              ! Left state
              qm(i,j,k,ir,1) = r + half*drx + sr0*dtdx*half
              qm(i,j,k,iu,1) = u + half*dux + su0*dtdx*half
              qm(i,j,k,ip,1) = p + half*dpx + sp0*dtdx*half
              qm(i,j,k,ir,1) = max(smallr, qm(i,j,k,ir,1))

        end do
     end do
  end do

  ! Passive scalars
  do n = ndim+3, nvar
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
                 a   = q(i,j,k,n)       ! Cell centered values
                 u   = q(i,j,k,iu)
                 dax = dq(i,j,k,n,1)    ! TVD slopes
                 sa0 = -u*dax             ! Source terms
                 qp(i,j,k,n,1) = a - half*dax + sa0*dtdx*half   ! Right state
                 qm(i,j,k,n,1) = a + half*dax + sa0*dtdx*half   ! Left state
           end do
        end do
     end do
  end do

end subroutine trace1d_gpu_2d
!###########################################################
!###########################################################
!###########################################################
!###########################################################
#if NDIM>1
subroutine trace2d_gpu_2d(q,dq,qm,qp,dx,dy,dt,nxp)
  use amr_parameters
  use hydro_parameters
  use acc_parameters
  use const
  implicit none

  integer ::nxp
  real(dp)::dx, dy, dt

  real(dp),dimension(-1:nxp+2,-1:nxp+2,-1:nxp+2,1:nvar)::q  
  real(dp),dimension(-1:nxp+2,-1:nxp+2,-1:nxp+2,1:nvar,1:ndim)::dq 
  real(dp),dimension(-1:nxp+2,-1:nxp+2,-1:nxp+2,1:nvar,1:ndim)::qm 
  real(dp),dimension(-1:nxp+2,-1:nxp+2,-1:nxp+2,1:nvar,1:ndim)::qp 

  ! declare local variables
  integer ::i, j, k, l, n
  integer ::ilo,ihi,jlo,jhi,klo,khi
  integer ::ir, iu, iv, ip
  real(dp)::dtdx, dtdy
  real(dp)::r, u, v, p, a
  real(dp)::drx, dux, dvx, dpx, dax
  real(dp)::dry, duy, dvy, dpy, day
  real(dp)::sr0, su0, sv0, sp0, sa0
  
  dtdx = dt/dx
  dtdy = dt/dy
  ilo=MIN(1,iu10+1); ihi=MAX(1,iu20-1)
  jlo=MIN(1,ju10+1); jhi=MAX(1,ju20-1)
  klo=MIN(1,ku10+1); khi=MAX(1,ku20-1)
  ir=1; iu=2; iv=3; ip=4

  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi

              ! Cell centered values
              r   =  q(i,j,k,ir)
              u   =  q(i,j,k,iu)
              v   =  q(i,j,k,iv)
              p   =  q(i,j,k,ip)

              ! TVD slopes in all directions
              drx = dq(i,j,k,ir,1)
              dux = dq(i,j,k,iu,1)
              dvx = dq(i,j,k,iv,1)
              dpx = dq(i,j,k,ip,1)
              
              dry = dq(i,j,k,ir,2)
              duy = dq(i,j,k,iu,2)
              dvy = dq(i,j,k,iv,2)
              dpy = dq(i,j,k,ip,2)
              
              ! source terms (with transverse derivatives)
              sr0 = -u*drx-v*dry - (dux+dvy)*r
              su0 = -u*dux-v*duy - (dpx    )/r
              sv0 = -u*dvx-v*dvy - (dpy    )/r
              sp0 = -u*dpx-v*dpy - (dux+dvy)*gamma*p

              ! Right state at left interface
              qp(i,j,k,ir,1) = r - half*drx + sr0*dtdx*half
              qp(i,j,k,iu,1) = u - half*dux + su0*dtdx*half
              qp(i,j,k,iv,1) = v - half*dvx + sv0*dtdx*half
              qp(i,j,k,ip,1) = p - half*dpx + sp0*dtdx*half
              qp(i,j,k,ir,1) = max(smallr, qp(i,j,k,ir,1))

              ! Left state at right interface
              qm(i,j,k,ir,1) = r + half*drx + sr0*dtdx*half
              qm(i,j,k,iu,1) = u + half*dux + su0*dtdx*half
              qm(i,j,k,iv,1) = v + half*dvx + sv0*dtdx*half
              qm(i,j,k,ip,1) = p + half*dpx + sp0*dtdx*half
              qm(i,j,k,ir,1) = max(smallr, qm(i,j,k,ir,1))

              ! Top state at bottom interface
              qp(i,j,k,ir,2) = r - half*dry + sr0*dtdy*half
              qp(i,j,k,iu,2) = u - half*duy + su0*dtdy*half
              qp(i,j,k,iv,2) = v - half*dvy + sv0*dtdy*half
              qp(i,j,k,ip,2) = p - half*dpy + sp0*dtdy*half
              qp(i,j,k,ir,2) = max(smallr, qp(i,j,k,ir,2))

              ! Bottom state at top interface
              qm(i,j,k,ir,2) = r + half*dry + sr0*dtdy*half
              qm(i,j,k,iu,2) = u + half*duy + su0*dtdy*half
              qm(i,j,k,iv,2) = v + half*dvy + sv0*dtdy*half
              qm(i,j,k,ip,2) = p + half*dpy + sp0*dtdy*half
              qm(i,j,k,ir,2) = max(smallr, qm(i,j,k,ir,2))

        end do
     end do
  end do

  ! passive scalars
  do n = ndim+3, nvar
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
                 a   = q(i,j,k,n)       ! Cell centered values
                 u   = q(i,j,k,iu)
                 v   = q(i,j,k,iv)
                 dax = dq(i,j,k,n,1)    ! TVD slopes
                 day = dq(i,j,k,n,2)
                 sa0 = -u*dax-v*day       ! Source terms
                 qp(i,j,k,n,1) = a - half*dax + sa0*dtdx*half   ! Right state
                 qm(i,j,k,n,1) = a + half*dax + sa0*dtdx*half   ! Left state
                 qp(i,j,k,n,2) = a - half*day + sa0*dtdy*half   ! Top state
                 qm(i,j,k,n,2) = a + half*day + sa0*dtdy*half   ! Bottom state
           end do
        end do
     end do
  end do

end subroutine trace2d_gpu_2d
#endif
!###########################################################
!###########################################################
!###########################################################
!###########################################################
#if NDIM>2
subroutine trace3d_gpu_2d(q,dq,qm,qp,dx,dy,dz,dt,nxp)
  use amr_parameters
  use hydro_parameters
  use acc_parameters
  use const
  implicit none

  INTEGER, INTENT(in) :: nxp
  REAL(dp), INTENT(in) :: dx, dy, dz, dt

  REAL(dp),DIMENSION(-1:nxp+2,-1:nxp+2,-1:nxp+2,1:nvar), &
       INTENT(in) :: q  
  REAL(dp),DIMENSION(-1:nxp+2,-1:nxp+2,-1:nxp+2,1:nvar,1:ndim), &
       INTENT(in) :: dq 
  REAL(dp),DIMENSION(-1:nxp+2,-1:nxp+2,-1:nxp+2,1:nvar,1:ndim), &
       INTENT(out) :: qm 
  REAL(dp),DIMENSION(-1:nxp+2,-1:nxp+2,-1:nxp+2,1:nvar,1:ndim), &
       INTENT(out) :: qp 

  ! declare local variables
  integer ::i, j, k, l, n
  integer ::ilo,ihi,jlo,jhi,klo,khi
!  integer ::ir, iu, iv, iw, ip
  real(dp)::dtdx, dtdy, dtdz
  real(dp)::r, u, v, w, p, a
  real(dp)::drx, dux, dvx, dwx, dpx, dax
  real(dp)::dry, duy, dvy, dwy, dpy, day
  real(dp)::drz, duz, dvz, dwz, dpz, daz
  real(dp)::sr0, su0, sv0, sw0, sp0, sa0

!!$ AH: These do not change, so use parameters to help the compiler
!!$ and avoid possible use of private variables/registers
  INTEGER, PARAMETER :: ir=1, iu=2, iv=3, iw=4, ip=5
 
! !$omp acc_data acc_copyin(q,dq,dx,dy,dz,dt) acc_copyout(qm,qp)
!$acc data present(q,dq,qm,qp) copyin(dx,dy,dz,dt)

  dtdx = dt/dx
  dtdy = dt/dy
  dtdz = dt/dz
  ilo=MIN(1,iu10+1); ihi=MAX(1,iu20-1)
  jlo=MIN(1,ju10+1); jhi=MAX(1,ju20-1)
  klo=MIN(1,ku10+1); khi=MAX(1,ku20-1)
!  ir=1; iu=2; iv=3; iw=4; ip=5

! !$omp acc_region_loop num_pes(2:NTPB)
!$acc parallel loop collapse(3) vector_length(NTPB)
  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           
           ! Cell centered values
           r   =  q(i,j,k,ir)
           u   =  q(i,j,k,iu)
           v   =  q(i,j,k,iv)
           w   =  q(i,j,k,iw)
           p   =  q(i,j,k,ip)
           
           ! TVD slopes in all 3 directions
           drx = dq(i,j,k,ir,1)
           dpx = dq(i,j,k,ip,1)
           dux = dq(i,j,k,iu,1)
           dvx = dq(i,j,k,iv,1)
           dwx = dq(i,j,k,iw,1)
           
           dry = dq(i,j,k,ir,2)
           dpy = dq(i,j,k,ip,2)
           duy = dq(i,j,k,iu,2)
           dvy = dq(i,j,k,iv,2)
           dwy = dq(i,j,k,iw,2)
           
           drz = dq(i,j,k,ir,3)
           dpz = dq(i,j,k,ip,3)
           duz = dq(i,j,k,iu,3)
           dvz = dq(i,j,k,iv,3)
           dwz = dq(i,j,k,iw,3)
           
           ! Source terms (including transverse derivatives)
           sr0 = -u*drx-v*dry-w*drz - (dux+dvy+dwz)*r
           sp0 = -u*dpx-v*dpy-w*dpz - (dux+dvy+dwz)*gamma*p
           su0 = -u*dux-v*duy-w*duz - (dpx        )/r
           sv0 = -u*dvx-v*dvy-w*dvz - (dpy        )/r
           sw0 = -u*dwx-v*dwy-w*dwz - (dpz        )/r
           
           ! Right state at left interface
           qp(i,j,k,ir,1) = r - half*drx + sr0*dtdx*half
           qp(i,j,k,ip,1) = p - half*dpx + sp0*dtdx*half
           qp(i,j,k,iu,1) = u - half*dux + su0*dtdx*half
           qp(i,j,k,iv,1) = v - half*dvx + sv0*dtdx*half
           qp(i,j,k,iw,1) = w - half*dwx + sw0*dtdx*half
           qp(i,j,k,ir,1) = max(smallr, qp(i,j,k,ir,1))
           
           ! Left state at left interface
           qm(i,j,k,ir,1) = r + half*drx + sr0*dtdx*half
           qm(i,j,k,ip,1) = p + half*dpx + sp0*dtdx*half
           qm(i,j,k,iu,1) = u + half*dux + su0*dtdx*half
           qm(i,j,k,iv,1) = v + half*dvx + sv0*dtdx*half
           qm(i,j,k,iw,1) = w + half*dwx + sw0*dtdx*half
           qm(i,j,k,ir,1) = max(smallr, qm(i,j,k,ir,1))
           
           ! Top state at bottom interface
           qp(i,j,k,ir,2) = r - half*dry + sr0*dtdy*half
           qp(i,j,k,ip,2) = p - half*dpy + sp0*dtdy*half
           qp(i,j,k,iu,2) = u - half*duy + su0*dtdy*half
           qp(i,j,k,iv,2) = v - half*dvy + sv0*dtdy*half
           qp(i,j,k,iw,2) = w - half*dwy + sw0*dtdy*half
           qp(i,j,k,ir,2) = max(smallr, qp(i,j,k,ir,2))
           
           ! Bottom state at top interface
           qm(i,j,k,ir,2) = r + half*dry + sr0*dtdy*half
           qm(i,j,k,ip,2) = p + half*dpy + sp0*dtdy*half
           qm(i,j,k,iu,2) = u + half*duy + su0*dtdy*half
           qm(i,j,k,iv,2) = v + half*dvy + sv0*dtdy*half
           qm(i,j,k,iw,2) = w + half*dwy + sw0*dtdy*half
           qm(i,j,k,ir,2) = max(smallr, qm(i,j,k,ir,2))
           
           ! Back state at front interface
           qp(i,j,k,ir,3) = r - half*drz + sr0*dtdz*half
           qp(i,j,k,ip,3) = p - half*dpz + sp0*dtdz*half
           qp(i,j,k,iu,3) = u - half*duz + su0*dtdz*half
           qp(i,j,k,iv,3) = v - half*dvz + sv0*dtdz*half
           qp(i,j,k,iw,3) = w - half*dwz + sw0*dtdz*half
           qp(i,j,k,ir,3) = max(smallr, qp(i,j,k,ir,3))
           
           ! Front state at back interface
           qm(i,j,k,ir,3) = r + half*drz + sr0*dtdz*half
           qm(i,j,k,ip,3) = p + half*dpz + sp0*dtdz*half
           qm(i,j,k,iu,3) = u + half*duz + su0*dtdz*half
           qm(i,j,k,iv,3) = v + half*dvz + sv0*dtdz*half
           qm(i,j,k,iw,3) = w + half*dwz + sw0*dtdz*half
           qm(i,j,k,ir,3) = max(smallr, qm(i,j,k,ir,3))
           
        end do
     end do
  end do
!$acc end parallel loop
! !$omp end acc_region_loop

!!$ AH: I think the restriction on n >= (ndim+3=6) means that the
!!$     two loopnests are independent

!!$ AH2: As NVAR=NDIM+2 for the builds we are considering, this
!!$      loopnest is never executed. The CCE realises this, so
!!$      there is no loopmark for this nest.

  ! Passive scalars
! !$omp acc_region_loop num_pes(2:NTPB)
!$acc parallel loop vector_length(NTPB)
  do n = ndim+3, nvar
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              a   = q(i,j,k,n)       ! Cell centered values
              u   = q(i,j,k,iu)
              v   = q(i,j,k,iv)
              w   = q(i,j,k,iw)
              dax = dq(i,j,k,n,1)    ! TVD slopes
              day = dq(i,j,k,n,2)
              daz = dq(i,j,k,n,3)
              sa0 = -u*dax-v*day-w*daz     ! Source terms
              qp(i,j,k,n,1) = a - half*dax + sa0*dtdx*half  ! Right state
              qm(i,j,k,n,1) = a + half*dax + sa0*dtdx*half  ! Left state
              qp(i,j,k,n,2) = a - half*day + sa0*dtdy*half  ! Bottom state
              qm(i,j,k,n,2) = a + half*day + sa0*dtdy*half  ! Upper state
              qp(i,j,k,n,3) = a - half*daz + sa0*dtdz*half  ! Front state
              qm(i,j,k,n,3) = a + half*daz + sa0*dtdz*half  ! Back state
           end do
        end do
     end do
  end do
!$acc end parallel loop
! !$omp end acc_region_loop

!$acc end data
! !$omp end acc_data

end subroutine trace3d_gpu_2d
#endif
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmpflxm_gpu_2d(qm,im1,im2,jm1,jm2,km1,km2, &
     &             qp,ip1,ip2,jp1,jp2,kp1,kp2, &
     &                ilo,ihi,jlo,jhi,klo,khi, ln,lt1,lt2, &
     &            flx,tmp,nxp)
  use amr_parameters
  use hydro_parameters
  use acc_parameters
  use const
  implicit none

  integer ::nxp
  integer ::ln,lt1,lt2
  integer ::im1,im2,jm1,jm2,km1,km2
  integer ::ip1,ip2,jp1,jp2,kp1,kp2
  integer ::ilo,ihi,jlo,jhi,klo,khi
  real(dp),dimension(im1:im2,jm1:jm2,km1:km2,1:nvar,1:ndim)::qm
  real(dp),dimension(ip1:ip2,jp1:jp2,kp1:kp2,1:nvar,1:ndim)::qp 
  real(dp),dimension(ip1:ip2,jp1:jp2,kp1:kp2,1:nvar)::flx
  real(dp),dimension(ip1:ip2,jp1:jp2,kp1:kp2,1:2)::tmp
  
  ! local variables
  integer ::i, j, k, n, l, idim, xdim
  real(dp)::entho
  real(dp),dimension(1:nvar)::qleft,qright,qgdnv,fgdnv

! Defining Collapse is dangerous. I think there is a bug.
#undef Collapse
#undef Nofuse

#ifdef Collapse
  INTEGER :: di,dj,dk,ijk,ij
#endif

  entho=one/(gamma-one)
  xdim=ln-1

#ifdef Collapse
  di = ihi-ilo
  dj = jhi-jlo
  dk = khi-klo
#endif

!$acc parallel loop present(qm,qp,flx,tmp) private(qleft,qright,qgdnv,fgdnv) &
!$acc& vector_length(NTPB) collapse(3) gang worker vector
!dir$ noblocking
  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi

           ! Mass density
           qleft (1) = qm(i,j,k,1,xdim)
           qright(1) = qp(i,j,k,1,xdim)
           
           ! Normal velocity
           qleft (2) = qm(i,j,k,ln,xdim)
           qright(2) = qp(i,j,k,ln,xdim)
           
           ! Pressure
           qleft (3) = qm(i,j,k,ndim+2,xdim)
           qright(3) = qp(i,j,k,ndim+2,xdim)
           
           ! Tangential velocity 1
#if NDIM>1
           qleft (4) = qm(i,j,k,lt1,xdim)
           qright(4) = qp(i,j,k,lt1,xdim)
#endif
           ! Tangential velocity 2
#if NDIM>2
           qleft (5) = qm(i,j,k,lt2,xdim)
           qright(5) = qp(i,j,k,lt2,xdim)
#endif           
           ! Other advected quantities
           do n = ndim+3, nvar
              qleft (n) = qm(i,j,k,n,xdim)
              qright(n) = qp(i,j,k,n,xdim)
           end do

!! CLAU AH VERSION          
           ! Solve Riemann problem
!           if(riemann.eq.'acoustic')then
!              call riemann_acoustic(qleft,qright,qgdnv,fgdnv,ngrid)
!           else if (riemann.eq.'exact')then
!              call riemann_approx  (qleft,qright,qgdnv,fgdnv,ngrid)
!           else if (riemann.eq.'llf')then
!              call riemann_llf     (qleft,qright,qgdnv,fgdnv,ngrid)
              call riemann_llf_scalar(qleft,qright,qgdnv,fgdnv)
!           else if (riemann.eq.'hllc')then
!              call riemann_hllc    (qleft,qright,qgdnv,fgdnv,ngrid)
!           else if (riemann.eq.'hll')then
!              call riemann_hll     (qleft,qright,qgdnv,fgdnv,ngrid)
!           else
!              write(*,*)'unknown Riemann solver'
!              stop
!           end if
           
           ! Compute fluxes
           
           ! Mass density
           flx(i,j,k,1) = fgdnv(1)
           
           ! Normal momentum
           flx(i,j,k,ln) = fgdnv(2)

           ! Transverse momentum 1
#if NDIM>1
           flx(i,j,k,lt1) = fgdnv(4)
#endif
           ! Transverse momentum 2
#if NDIM>2
           flx(i,j,k,lt2) = fgdnv(5)
#endif           
           ! Total energy
           flx(i,j,k,ndim+2) = fgdnv(3)

           ! Other advected quantities
           do n = ndim+3, nvar
              flx(i,j,k,n) = fgdnv(n)
           end do

           ! Temporary Godunov states
           tmp(i,j,k,1) = qgdnv(2)   ! Normal velocity
           if(qgdnv(2)>zero)then       ! Internal energy flux
              tmp(i,j,k,2) = qleft (3)*qgdnv(2)*entho
           else
              tmp(i,j,k,2) = qright(3)*qgdnv(2)*entho
           end if

        end do
     end do
  end do
!$acc end parallel loop

end subroutine cmpflxm_gpu_2d
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine ctoprim_gpu_2d(uin,q,c,gravin,dt,nxp)
  use amr_parameters
  use hydro_parameters
  use acc_parameters
  use const
  implicit none

  integer ::nxp
  real(dp)::dt
  real(dp),dimension(-1:nxp+2,-1:nxp+2,-1:nxp+2,1:nvar)::uin
  real(dp),dimension(-1:nxp+2,-1:nxp+2,-1:nxp+2,1:ndim)::gravin
  real(dp),dimension(-1:nxp+2,-1:nxp+2,-1:nxp+2,1:nvar)::q  
  real(dp),dimension(-1:nxp+2,-1:nxp+2,-1:nxp+2)::c  

  integer ::i, j, k, l, n, idim
  real(dp)::eint, smalle, dtxhalf, oneonrho
  real(dp)::eken

!$acc data present(uin,q,c,gravin)

  smalle = smallc**2/gamma/(gamma-one)
  dtxhalf = dt*half

  ! Convert to primitive variable
!$acc parallel loop collapse(3) gang worker vector vector_length(NTPB)
  do k = ku10, ku20
     do j = ju10, ju20
        do i = iu10, iu20

           ! Compute density
           q(i,j,k,1) = max(uin(i,j,k,1),smallr)
           
           ! Compute velocities
           oneonrho = 1.d0/q(i,j,k,1)
           q(i,j,k,2) = uin(i,j,k,2)*oneonrho
#if NDIM>1
           q(i,j,k,3) = uin(i,j,k,3)*oneonrho
#endif
#if NDIM>2
           q(i,j,k,4) = uin(i,j,k,4)*oneonrho
#endif
           
           ! Compute specific kinetic energy
           eken = half*q(i,j,k,2)*q(i,j,k,2)
#if NDIM>1
           eken = eken + half*q(i,j,k,3)*q(i,j,k,3)
#endif
#if NDIM>2
           eken = eken + half*q(i,j,k,4)*q(i,j,k,4)
#endif
           
           ! Compute pressure
           eint = MAX(uin(i,j,k,ndim+2)*oneonrho-eken,smalle)
           q(i,j,k,ndim+2)=(gamma-one)*q(i,j,k,1)*eint
           
           ! Compute sound speed
           c(i,j,k)=sqrt(gamma*q(i,j,k,ndim+2)*oneonrho)
           
           ! Gravity predictor step
           q(i,j,k,2) = q(i,j,k,2) + gravin(i,j,k,1)*dtxhalf
#if NDIM>1
           q(i,j,k,3) = q(i,j,k,3) + gravin(i,j,k,2)*dtxhalf
#endif
#if NDIM>2
           q(i,j,k,4) = q(i,j,k,4) + gravin(i,j,k,3)*dtxhalf
#endif
           
        end do
     end do
  end do
!$acc end parallel loop

  ! Passive scalar
!$acc parallel loop vector_length(NTPB)
  do n = ndim+3, nvar
     do k = ku10, ku20
        do j = ju10, ju20
           do i = iu10, iu20
              oneonrho = 1.d0/q(i,j,k,1)
              q(i,j,k,n) = uin(i,j,k,n)*oneonrho
           end do
        end do
     end do
  end do
!$acc end parallel loop 

!$acc end data

end subroutine ctoprim_gpu_2d
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine uslope_gpu_2d(q,dq,dx,dt,nxp)
  use amr_parameters
  use hydro_parameters
  use acc_parameters
  use const
  implicit none

  integer::nxp
  real(dp)::dx,dt
  real(dp),dimension(-1:nxp+2,-1:nxp+2,-1:nxp+2,1:nvar)::q 
  real(dp),dimension(-1:nxp+2,-1:nxp+2,-1:nxp+2,1:nvar,1:ndim)::dq

  ! local arrays
  integer::i, j, k, l, n
  real(dp)::dsgn, dlim, dcen, dlft, drgt, slop
  real(dp)::dfll,dflm,dflr,dfmdfmm,dfmr,dfrdfrm,dfrr
  real(dp)::dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl
  real(dp)::dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm
  real(dp)::dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr
  real(dp)::vmin,vmax,dfx,dfy,dfz,dff
  integer::ilo,ihi,jlo,jhi,klo,khi
  
  ilo=MIN(1,iu10+1); ihi=MAX(1,iu20-1)
  jlo=MIN(1,ju10+1); jhi=MAX(1,ju20-1)
  klo=MIN(1,ku10+1); khi=MAX(1,ku20-1)

  if(slope_type==0)then
     dq=zero
     return
  end if

!$acc data present(q,dq)
#if NDIM==1
  do n = 1, nvar
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              if(slope_type==1.or.slope_type==2.or.slope_type==3)then  ! minmod or average
                    dlft = MIN(slope_type,2)*(q(i  ,j,k,n) - q(i-1,j,k,n))
                    drgt = MIN(slope_type,2)*(q(i+1,j,k,n) - q(i  ,j,k,n))
                    dcen = half*(dlft+drgt)/MIN(slope_type,2)
                    dsgn = sign(one, dcen)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(i,j,k,n,1) = dsgn*min(dlim,abs(dcen))
              else if(slope_type==4)then ! superbee
                    dcen = q(i,j,k,2)*dt/dx
                    dlft = two/(one+dcen)*(q(i,j,k,n)-q(i-1,j,k,n))
                    drgt = two/(one-dcen)*(q(i+1,j,k,n)-q(i,j,k,n))
                    dcen = half*(q(i+1,j,k,n)-q(i-1,j,k,n))
                    dsgn = sign(one, dlft)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(i,j,k,n,1) = dsgn*dlim !min(dlim,abs(dcen))
              else if(slope_type==5)then ! ultrabee
                 if(n==1)then
                       dcen = q(i,j,k,2)*dt/dx
                       if(dcen>=0)then
                          dlft = two/(zero+dcen+1d-10)*(q(i,j,k,n)-q(i-1,j,k,n))
                          drgt = two/(one -dcen      )*(q(i+1,j,k,n)-q(i,j,k,n))
                       else
                          dlft = two/(one +dcen      )*(q(i,j,k,n)-q(i-1,j,k,n))
                          drgt = two/(zero-dcen+1d-10)*(q(i+1,j,k,n)-q(i,j,k,n))
                       endif
                       dsgn = sign(one, dlft)
                       slop = min(abs(dlft),abs(drgt))
                       dlim = slop
                       dcen = half*(q(i+1,j,k,n)-q(i-1,j,k,n))
                       if((dlft*drgt)<=zero)dlim=zero
                       dq(i,j,k,n,1) = dsgn*dlim !min(dlim,abs(dcen))
                 else
                       dq(i,j,k,n,1) = 0.0
                 end if
              else if(slope_type==6)then ! unstable
                 if(n==1)then
                       dlft = (q(i,j,k,n)-q(i-1,j,k,n))
                       drgt = (q(i+1,j,k,n)-q(i,j,k,n))
                       slop = 0.5*(dlft+drgt)
                       dlim = slop
                       dq(i,j,k,n,1) = dlim
                 else
                       dq(i,j,k,n,1) = 0.0
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
                    dlft = slope_type*(q(i  ,j,k,n) - q(i-1,j,k,n))
                    drgt = slope_type*(q(i+1,j,k,n) - q(i  ,j,k,n))
                    dcen = half*(dlft+drgt)/slope_type
                    dsgn = sign(one, dcen)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(i,j,k,n,1) = dsgn*min(dlim,abs(dcen))
                 ! slopes in second coordinate direction
                    dlft = slope_type*(q(i,j  ,k,n) - q(i,j-1,k,n))
                    drgt = slope_type*(q(i,j+1,k,n) - q(i,j  ,k,n))
                    dcen = half*(dlft+drgt)/slope_type
                    dsgn = sign(one,dcen)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(i,j,k,n,2) = dsgn*min(dlim,abs(dcen))
              end do
           end do
        end do
     end do
  else if(slope_type==3)then ! positivity preserving 2d unsplit slope
     do n = 1, nvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                    dfll = q(i-1,j-1,k,n)-q(i,j,k,n)
                    dflm = q(i-1,j  ,k,n)-q(i,j,k,n)
                    dflr = q(i-1,j+1,k,n)-q(i,j,k,n)
                    dfml = q(i  ,j-1,k,n)-q(i,j,k,n)
                    dfmm = q(i  ,j  ,k,n)-q(i,j,k,n)
                    dfmr = q(i  ,j+1,k,n)-q(i,j,k,n)
                    dfrl = q(i+1,j-1,k,n)-q(i,j,k,n)
                    dfrm = q(i+1,j  ,k,n)-q(i,j,k,n)
                    dfrr = q(i+1,j+1,k,n)-q(i,j,k,n)
                    
                    vmin = min(dfldflm,dflr,dfmdfmm,dfmr,dfrdfrm,dfrr)
                    vmax = max(dfldflm,dflr,dfmdfmm,dfmr,dfrdfrm,dfrr)
                    
                    dfx  = half*(q(i+1,j,k,n)-q(i-1,j,k,n))
                    dfy  = half*(q(i,j+1,k,n)-q(i,j-1,k,n))
                    dff  = half*(abs(dfx)+abs(dfy))
                    
                    if(dff>zero)then
                       slop = min(one,min(abs(vmin),abs(vmax))/dff)
                    else
                       slop = one
                    endif
                    
                    dlim = slop
                    
                    dq(i,j,k,n,1) = dlim*dfx
                    dq(i,j,k,n,2) = dlim*dfy

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
!$acc parallel loop collapse(4) gang worker vector vector_length(NTPB)
     do n = 1, nvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 ! slopes in first coordinate direction
                 dlft = q(i  ,j,k,n) - q(i-1,j,k,n)
                 drgt = q(i+1,j,k,n) - q(i  ,j,k,n)
                 if((dlft*drgt)<=zero) then
                    dq(i,j,k,n,1) = zero
                 else if(dlft>0) then
                    dq(i,j,k,n,1) = min(dlft,drgt)
                 else
                    dq(i,j,k,n,1) = max(dlft,drgt)
                 end if
                 ! slopes in second coordinate direction
                 dlft = q(i,j  ,k,n) - q(i,j-1,k,n)
                 drgt = q(i,j+1,k,n) - q(i,j  ,k,n)
                 if((dlft*drgt)<=zero) then
                    dq(i,j,k,n,2) = zero
                 else if(dlft>0) then
                    dq(i,j,k,n,2) = min(dlft,drgt)
                 else
                    dq(i,j,k,n,2) = max(dlft,drgt)
                 end if
                 ! slopes in third coordinate direction
                 dlft = q(i,j,k  ,n) - q(i,j,k-1,n)
                 drgt = q(i,j,k+1,n) - q(i,j,k  ,n)
                 if((dlft*drgt)<=zero) then
                    dq(i,j,k,n,3) = zero
                 else if(dlft>0) then
                    dq(i,j,k,n,3) = min(dlft,drgt)
                 else
                    dq(i,j,k,n,3) = max(dlft,drgt)
                 end if
              end do
           end do
        end do
     end do
!$acc end parallel loop
  else if(slope_type==2)then ! moncen
     do n = 1, nvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 ! slopes in first coordinate direction
                 dlft = slope_type*(q(i  ,j,k,n) - q(i-1,j,k,n))
                 drgt = slope_type*(q(i+1,j,k,n) - q(i  ,j,k,n))
                 dcen = half*(dlft+drgt)/slope_type
                 dsgn = sign(one, dcen)
                 slop = min(abs(dlft),abs(drgt))
                 dlim = slop
                 if((dlft*drgt)<=zero)dlim=zero
                 dq(i,j,k,n,1) = dsgn*min(dlim,abs(dcen))
                 ! slopes in second coordinate direction
                 dlft = slope_type*(q(i,j  ,k,n) - q(i,j-1,k,n))
                 drgt = slope_type*(q(i,j+1,k,n) - q(i,j  ,k,n))
                 dcen = half*(dlft+drgt)/slope_type
                 dsgn = sign(one,dcen)
                 slop = min(abs(dlft),abs(drgt))
                 dlim = slop
                 if((dlft*drgt)<=zero)dlim=zero
                 dq(i,j,k,n,2) = dsgn*min(dlim,abs(dcen))
                 ! slopes in third coordinate direction
                 dlft = slope_type*(q(i,j,k  ,n) - q(i,j,k-1,n))
                 drgt = slope_type*(q(i,j,k+1,n) - q(i,j,k  ,n))
                 dcen = half*(dlft+drgt)/slope_type
                 dsgn = sign(one,dcen)
                 slop = min(abs(dlft),abs(drgt))
                 dlim = slop
                 if((dlft*drgt)<=zero)dlim=zero
                 dq(i,j,k,n,3) = dsgn*min(dlim,abs(dcen))
              end do
           end do
        end do
     end do
  else if(slope_type==3)then ! positivity preserving 3d unsplit slope
     do n = 1, nvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 dflll = q(i-1,j-1,k-1,n)-q(i,j,k,n)
                 dflml = q(i-1,j  ,k-1,n)-q(i,j,k,n)
                 dflrl = q(i-1,j+1,k-1,n)-q(i,j,k,n)
                 dfmll = q(i  ,j-1,k-1,n)-q(i,j,k,n)
                 dfmml = q(i  ,j  ,k-1,n)-q(i,j,k,n)
                 dfmrl = q(i  ,j+1,k-1,n)-q(i,j,k,n)
                 dfrll = q(i+1,j-1,k-1,n)-q(i,j,k,n)
                 dfrml = q(i+1,j  ,k-1,n)-q(i,j,k,n)
                 dfrrl = q(i+1,j+1,k-1,n)-q(i,j,k,n)
                 
                 dfllm = q(i-1,j-1,k  ,n)-q(i,j,k,n)
                 dflmm = q(i-1,j  ,k  ,n)-q(i,j,k,n)
                 dflrm = q(i-1,j+1,k  ,n)-q(i,j,k,n)
                 dfmlm = q(i  ,j-1,k  ,n)-q(i,j,k,n)
                 dfmmm = q(i  ,j  ,k  ,n)-q(i,j,k,n)
                 dfmrm = q(i  ,j+1,k  ,n)-q(i,j,k,n)
                 dfrlm = q(i+1,j-1,k  ,n)-q(i,j,k,n)
                 dfrmm = q(i+1,j  ,k  ,n)-q(i,j,k,n)
                 dfrrm = q(i+1,j+1,k  ,n)-q(i,j,k,n)
                 
                 dfllr = q(i-1,j-1,k+1,n)-q(i,j,k,n)
                 dflmr = q(i-1,j  ,k+1,n)-q(i,j,k,n)
                 dflrr = q(i-1,j+1,k+1,n)-q(i,j,k,n)
                 dfmlr = q(i  ,j-1,k+1,n)-q(i,j,k,n)
                 dfmmr = q(i  ,j  ,k+1,n)-q(i,j,k,n)
                 dfmrr = q(i  ,j+1,k+1,n)-q(i,j,k,n)
                 dfrlr = q(i+1,j-1,k+1,n)-q(i,j,k,n)
                 dfrmr = q(i+1,j  ,k+1,n)-q(i,j,k,n)
                 dfrrr = q(i+1,j+1,k+1,n)-q(i,j,k,n)
                 
                 
                 vmin = min(dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl, &
                      &     dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm, &
                      &     dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr)
                 vmax = max(dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl, &
                      &     dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm, &
                      &     dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr)
                 
                 dfx  = half*(q(i+1,j,k,n)-q(i-1,j,k,n))
                 dfy  = half*(q(i,j+1,k,n)-q(i,j-1,k,n))
                 dfz  = half*(q(i,j,k+1,n)-q(i,j,k-1,n))
                 dff  = half*(abs(dfx)+abs(dfy)+abs(dfz))
                 
                 if(dff>zero)then
                    slop = min(one,min(abs(vmin),abs(vmax))/dff)
                 else
                    slop = one
                 endif
                 
                 dlim = slop
                 
                 dq(i,j,k,n,1) = dlim*dfx
                 dq(i,j,k,n,2) = dlim*dfy
                 dq(i,j,k,n,3) = dlim*dfz
                 
              end do
           end do
        end do
     end do
  else
     write(*,*)'Unknown slope type'
     stop
  endif     
#endif

!$acc end data
  
end subroutine uslope_gpu_2d
