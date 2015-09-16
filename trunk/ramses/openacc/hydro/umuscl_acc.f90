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
subroutine unsplit_gpu(uin,gravin,dx,nxp,dt,ncube,tot_cube,flux,tmp)
  use amr_parameters
  use const             
  use hydro_parameters
  use acc_parameters
  use acc_commons, only:dtdx,dtdy,dtdz,qin,cin,dq,qm,qp
  implicit none 

  integer ::nxp,ncube,tot_cube
  real(dp)::dx,dy,dz,dt

  ! Input states
  real(dp),dimension(-1:nxp+2,-1:nxp+2,-1:nxp+2,1:nvar,1:tot_cube)::uin 
  real(dp),dimension(-1:nxp+2,-1:nxp+2,-1:nxp+2,1:ndim,1:tot_cube)::gravin 

  ! Output fluxes
  real(dp),dimension(1:nxp+1,1:nxp+1,1:nxp+1,1:nvar,1:ndim,1:ncube)::flux
  real(dp),dimension(1:nxp+1,1:nxp+1,1:nxp+1,1:2   ,1:ndim,1:ncube)::tmp 

!!!!!! REMEMBER DIVU !!!!!!!!!!!
  ! Velocity divergence
  !real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2)::divu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Local scalar variables
  integer::i,j,k,l,ivar,nc
  integer::i0,j0,k0
  
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


#define NTPB 128
  
!$acc data present(uin,gravin,tmp,flux,qin,cin,dq,qm,qp)
  
  ! Translate to primative variables, compute sound speeds  
  call ctoprim_gpu_2d(uin,gravin,dt,nxp,ncube)
  
  ! Compute TVD slopes
  call uslope_gpu_2d(dx,dt,nxp,ncube)
  
  ! Compute 3D traced-states in all three directions
  if(scheme=='muscl')then
#if NDIM==1
     call trace1d_gpu_2d(qin,dq,qm,qp,dx      ,dt,nxp)
#endif
#if NDIM==2
     call trace2d_gpu_2d(qin,dq,qm,qp,dx,dy   ,dt,nxp)
#endif
#if NDIM==3
     call trace3d_gpu_2d(nxp,ncube)! async(3)
#endif
  endif
  if(scheme=='plmde')then
#if NDIM==1
     call tracex_gpu_2d(qin,dq,cin,qm,qp,dx,dt,nxp)
#endif
#if NDIM==2
     call tracexy_gpu_2d(qin,dq,cin,qm,qp,dx,dy,dt,nxp)
#endif
#if NDIM==3
     call tracexyz_gpu_2d(qin,dq,cin,qm,qp,dx,dy,dz,dt,nxp,ncube)
#endif
  endif
  
  ! Solve for 1D flux in X direction
  i0 = 1; j0 = 0; k0 = 0;
  call cmpflxm_gpu_2d(flux,tmp,dtdx,i0,j0,k0,1,2,3,4,nxp,ncube) 
  
  ! Solve for 1D flux in Y direction
  i0 = 0; j0 = 1; k0 = 0;
  call cmpflxm_gpu_2d(flux,tmp,dtdy,i0,j0,k0,2,3,2,4,nxp,ncube)
  
  ! Solve for 1D flux in Z direction
  i0 = 0; j0 = 0; k0 = 1;
  call cmpflxm_gpu_2d(flux,tmp,dtdz,i0,j0,k0,3,4,2,3,nxp,ncube) 
  
  
  ! DIVU STUFF, SEE HOST CODE
  

!$acc end data

end subroutine unsplit_gpu
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmpflxm_gpu_2d(flux,tmp,dtdx,i0,j0,k0,idim,ln,lt1,lt2,nxp,ncube)
  use amr_parameters
  use hydro_parameters
  use acc_parameters
  use acc_commons, only: qm,qp
  use const
  implicit none

  integer ::nxp,ncube,idim
  integer ::ln,lt1,lt2
  integer ::i0,j0,k0
  real(dp)::dtdx
  
  real(dp),dimension(1:nxp+1,1:nxp+1,1:nxp+1,1:nvar,1:ndim,1:ncube)::flux
  real(dp),dimension(1:nxp+1,1:nxp+1,1:nxp+1,1:2   ,1:ndim,1:ncube)::tmp
  
  ! local variables
  integer ::i, j, k, n, l, nc
  real(dp),dimension(1:nvar)::qleft,qright,qgdnv,fgdnv

!$acc data present(qm,qp,dtdx,flux,tmp)

!$acc parallel loop collapse(4) gang vector private(qleft,qright,qgdnv,fgdnv) !async(3)
  do nc=1,ncube
   do k = 1, nxp+k0
      do j = 1, nxp+j0
         do i = 1, nxp+i0

            ! Mass density
            qleft (1) = qm(i-i0,j-j0,k-k0,1,idim,nc)
            qright(1) = qp(i,j,k,1,idim,nc)

            ! Normal velocity
            qleft (2) = qm(i-i0,j-j0,k-k0,ln,idim,nc)
            qright(2) = qp(i,j,k,ln,idim,nc)

            ! Pressure
            qleft (3) = qm(i-i0,j-j0,k-k0,ndim+2,idim,nc)
            qright(3) = qp(i,j,k,ndim+2,idim,nc)

            ! Tangential velocity 1
#if NDIM>1
            qleft (4) = qm(i-i0,j-j0,k-k0,lt1,idim,nc)
            qright(4) = qp(i,j,k,lt1,idim,nc)
#endif
            ! Tangential velocity 2
#if NDIM>2
            qleft (5) = qm(i-i0,j-j0,k-k0,lt2,idim,nc)
            qright(5) = qp(i,j,k,lt2,idim,nc)
#endif           
            ! Other advected quantities
            do n = ndim+3, nvar
               qleft (n) = qm(i-i0,j-j0,k-k0,n,idim,nc)
               qright(n) = qp(i,j,k,n,idim,nc)
            end do

 !! CLAU AH VERSION          
            ! Solve Riemann problem
 !           if(riemann.eq.'acoustic')then
 !              call riemann_acoustic(qleft,qright,fgdnv,ngrid)
 !           else if (riemann.eq.'exact')then
 !              call riemann_approx  (qleft,qright,fgdnv,ngrid)
 !           else if (riemann.eq.'llf')then
 !              call riemann_llf     (qleft,qright,fgdnv,ngrid)
               call riemann_llf_scalar(qleft,qright,qgdnv,fgdnv)
 !           else if (riemann.eq.'hllc')then
 !              call riemann_hllc    (qleft,qright,fgdnv,ngrid)
 !           else if (riemann.eq.'hll')then
 !              call riemann_hll     (qleft,qright,fgdnv,ngrid)
 !           else
 !              write(*,*)'unknown Riemann solver'
 !              stop
 !           end if

            ! Compute fluxes

            ! Mass density
            flux(i,j,k,1,idim,nc) = fgdnv(1)*dtdx

            ! Normal momentum
            flux(i,j,k,ln,idim,nc) = fgdnv(2)*dtdx

            ! Transverse momentum 1
#if NDIM>1
 ! prima riga in cui fallisce 5000
            flux(i,j,k,lt1,idim,nc) = fgdnv(4)*dtdx
#endif
            ! Transverse momentum 2
#if NDIM>2
            flux(i,j,k,lt2,idim,nc) = fgdnv(5)*dtdx
#endif           
            ! Total energy
            flux(i,j,k,ndim+2,idim,nc) = fgdnv(3)*dtdx

            ! Other advected quantities
            do n = ndim+3, nvar
               flux(i,j,k,n,idim,nc) = fgdnv(n)*dtdx
            end do

            ! Temporary Godunov states
            tmp(i,j,k,1,idim,nc) = qgdnv(2)*dtdx   ! Normal velocity
            if(qgdnv(2)>zero)then       ! Internal energy flux
               tmp(i,j,k,2,idim,nc) = (qleft(3)*qgdnv(2)/(gamma-one))*dtdx
            else
               tmp(i,j,k,2,idim,nc) = (qright(3)*qgdnv(2)/(gamma-one))*dtdx
            end if
         end do
      end do
   end do
  end do
  
!$acc end parallel loop
  
!$acc end data

end subroutine cmpflxm_gpu_2d
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

#if NVAR > NDIM + 2
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
#endif

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
subroutine trace3d_gpu_2d(nxp,ncube)
  use amr_parameters
  use hydro_parameters
  use acc_parameters
  use acc_commons, only: dtdx,dtdy,dtdz,qin,dq,qm,qp
  use const
  implicit none

  integer, intent(in) :: nxp,ncube

  ! declare local variables
  integer, parameter :: ir=1, iu=2, iv=3, iw=4, ip=5
  integer ::i, j, k, l, n, nc
  real(dp)::r, u, v, w, p, a
  real(dp)::drx, dux, dvx, dwx, dpx, dax
  real(dp)::dry, duy, dvy, dwy, dpy, day
  real(dp)::drz, duz, dvz, dwz, dpz, daz
  real(dp)::sr0, su0, sv0, sw0, sp0, sa0

  
!$acc data present(qin,dq,qm,qp,dtdx,dtdy,dtdz)
  
!$acc parallel loop collapse(4) gang vector !async(3)
  do nc=1,ncube
   do k = 0, nxp+1
      do j = 0, nxp+1
         do i = 0, nxp+1

            ! Cell centered values
            r   =  qin(i,j,k,ir,nc)
            u   =  qin(i,j,k,iu,nc)
            v   =  qin(i,j,k,iv,nc)
            w   =  qin(i,j,k,iw,nc)
            p   =  qin(i,j,k,ip,nc)

            ! TVD slopes in all 3 directions
            drx = dq(i,j,k,ir,1,nc)
            dpx = dq(i,j,k,ip,1,nc)
            dux = dq(i,j,k,iu,1,nc)
            dvx = dq(i,j,k,iv,1,nc)
            dwx = dq(i,j,k,iw,1,nc)

            dry = dq(i,j,k,ir,2,nc)
            dpy = dq(i,j,k,ip,2,nc)
            duy = dq(i,j,k,iu,2,nc)
            dvy = dq(i,j,k,iv,2,nc)
            dwy = dq(i,j,k,iw,2,nc)

            drz = dq(i,j,k,ir,3,nc)
            dpz = dq(i,j,k,ip,3,nc)
            duz = dq(i,j,k,iu,3,nc)
            dvz = dq(i,j,k,iv,3,nc)
            dwz = dq(i,j,k,iw,3,nc)

            ! Source terms (including transverse derivatives)
            sr0 = -u*drx-v*dry-w*drz - (dux+dvy+dwz)*r
            sp0 = -u*dpx-v*dpy-w*dpz - (dux+dvy+dwz)*gamma*p
            su0 = -u*dux-v*duy-w*duz - (dpx        )/r
            sv0 = -u*dvx-v*dvy-w*dvz - (dpy        )/r
            sw0 = -u*dwx-v*dwy-w*dwz - (dpz        )/r

            ! Right state at left interface
            qp(i,j,k,ir,1,nc) = r - half*drx + sr0*dtdx*half
            qp(i,j,k,ip,1,nc) = p - half*dpx + sp0*dtdx*half
            qp(i,j,k,iu,1,nc) = u - half*dux + su0*dtdx*half
            qp(i,j,k,iv,1,nc) = v - half*dvx + sv0*dtdx*half
            qp(i,j,k,iw,1,nc) = w - half*dwx + sw0*dtdx*half
            qp(i,j,k,ir,1,nc) = max(smallr, qp(i,j,k,ir,1,nc))

            ! Left state at left interface
            qm(i,j,k,ir,1,nc) = r + half*drx + sr0*dtdx*half
            qm(i,j,k,ip,1,nc) = p + half*dpx + sp0*dtdx*half
            qm(i,j,k,iu,1,nc) = u + half*dux + su0*dtdx*half
            qm(i,j,k,iv,1,nc) = v + half*dvx + sv0*dtdx*half
            qm(i,j,k,iw,1,nc) = w + half*dwx + sw0*dtdx*half
            qm(i,j,k,ir,1,nc) = max(smallr, qm(i,j,k,ir,1,nc))

            ! Top state at bottom interface
            qp(i,j,k,ir,2,nc) = r - half*dry + sr0*dtdy*half
            qp(i,j,k,ip,2,nc) = p - half*dpy + sp0*dtdy*half
            qp(i,j,k,iu,2,nc) = u - half*duy + su0*dtdy*half
            qp(i,j,k,iv,2,nc) = v - half*dvy + sv0*dtdy*half
            qp(i,j,k,iw,2,nc) = w - half*dwy + sw0*dtdy*half
            qp(i,j,k,ir,2,nc) = max(smallr, qp(i,j,k,ir,2,nc))

            ! Bottom state at top interface
            qm(i,j,k,ir,2,nc) = r + half*dry + sr0*dtdy*half
            qm(i,j,k,ip,2,nc) = p + half*dpy + sp0*dtdy*half
            qm(i,j,k,iu,2,nc) = u + half*duy + su0*dtdy*half
            qm(i,j,k,iv,2,nc) = v + half*dvy + sv0*dtdy*half
            qm(i,j,k,iw,2,nc) = w + half*dwy + sw0*dtdy*half
            qm(i,j,k,ir,2,nc) = max(smallr, qm(i,j,k,ir,2,nc))

            ! Back state at front interface
            qp(i,j,k,ir,3,nc) = r - half*drz + sr0*dtdz*half
            qp(i,j,k,ip,3,nc) = p - half*dpz + sp0*dtdz*half
            qp(i,j,k,iu,3,nc) = u - half*duz + su0*dtdz*half
            qp(i,j,k,iv,3,nc) = v - half*dvz + sv0*dtdz*half
            qp(i,j,k,iw,3,nc) = w - half*dwz + sw0*dtdz*half
            qp(i,j,k,ir,3,nc) = max(smallr, qp(i,j,k,ir,3,nc))

            ! Front state at back interface
            qm(i,j,k,ir,3,nc) = r + half*drz + sr0*dtdz*half
            qm(i,j,k,ip,3,nc) = p + half*dpz + sp0*dtdz*half
            qm(i,j,k,iu,3,nc) = u + half*duz + su0*dtdz*half
            qm(i,j,k,iv,3,nc) = v + half*dvz + sv0*dtdz*half
            qm(i,j,k,iw,3,nc) = w + half*dwz + sw0*dtdz*half
            qm(i,j,k,ir,3,nc) = max(smallr, qm(i,j,k,ir,3,nc))

         end do
      end do
   end do
  end do
!$acc end parallel loop
  
  
#if NVAR > NDIM + 2
  ! Passive scalars
!$acc parallel loop collapse(5) gang vector
  do nc=1,ncube
   do n = ndim+3, nvar
      do k = 0, nxp+1
         do j = 0, nxp+1
            do i = 0, nxp+1
               a   = qin(i,j,k,n,nc)       ! Cell centered values
               u   = qin(i,j,k,iu,nc)
               v   = qin(i,j,k,iv,nc)
               w   = qin(i,j,k,iw,nc)
               dax = dq(i,j,k,n,1,nc)    ! TVD slopes
               day = dq(i,j,k,n,2,nc)
               daz = dq(i,j,k,n,3,nc)
               sa0 = -u*dax-v*day-w*daz     ! Source terms
               qp(i,j,k,n,1,nc) = a - half*dax + sa0*dtdx*half  ! Right state
               qm(i,j,k,n,1,nc) = a + half*dax + sa0*dtdx*half  ! Left state
               qp(i,j,k,n,2,nc) = a - half*day + sa0*dtdy*half  ! Bottom state
               qm(i,j,k,n,2,nc) = a + half*day + sa0*dtdy*half  ! Upper state
               qp(i,j,k,n,3,nc) = a - half*daz + sa0*dtdz*half  ! Front state
               qm(i,j,k,n,3,nc) = a + half*daz + sa0*dtdz*half  ! Back state
            end do
         end do
      end do
   end do
  end do
!$acc end parallel loop
#endif

!$acc end data

end subroutine trace3d_gpu_2d
#endif
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine ctoprim_gpu_2d(uin,gravin,dt,nxp,ncube)
  use amr_parameters
  use hydro_parameters
  use acc_commons, only:qin,cin
  use const
  implicit none

  integer ::nxp,ncube
  real(dp)::dt
  real(dp),dimension(-1:nxp+2,-1:nxp+2,-1:nxp+2,1:nvar,1:ncube)::uin
  real(dp),dimension(-1:nxp+2,-1:nxp+2,-1:nxp+2,1:ndim,1:ncube)::gravin

  integer ::i, j, k, l, n, idim, nc
  real(dp)::eint, smalle, dtxhalf, oneonrho
  real(dp)::eken

!$acc data present(uin,qin,cin,gravin)

  smalle = smallc**2/gamma/(gamma-one)
  dtxhalf = dt*half

  ! Convert to primitive variable
!$acc parallel loop collapse(4) gang vector !async(3)
  do nc=1,ncube
   do k = -1, nxp+2
      do j = -1, nxp+2
         do i = -1, nxp+2

            ! Compute density
            qin(i,j,k,1,nc) = max(uin(i,j,k,1,nc),smallr)

            ! Compute velocities
            oneonrho = 1.d0/qin(i,j,k,1,nc)
            qin(i,j,k,2,nc) = uin(i,j,k,2,nc)*oneonrho
#if NDIM>1
            qin(i,j,k,3,nc) = uin(i,j,k,3,nc)*oneonrho
#endif
#if NDIM>2
            qin(i,j,k,4,nc) = uin(i,j,k,4,nc)*oneonrho
#endif

            ! Compute specific kinetic energy
            eken = half*qin(i,j,k,2,nc)*qin(i,j,k,2,nc)
#if NDIM>1
            eken = eken + half*qin(i,j,k,3,nc)*qin(i,j,k,3,nc)
#endif
#if NDIM>2
            eken = eken + half*qin(i,j,k,4,nc)*qin(i,j,k,4,nc)
#endif

            ! Compute pressure
            eint = MAX(uin(i,j,k,ndim+2,nc)*oneonrho-eken,smalle)
            qin(i,j,k,ndim+2,nc)=(gamma-one)*qin(i,j,k,1,nc)*eint

            ! Compute sound speed
            cin(i,j,k,nc)=sqrt(gamma*qin(i,j,k,ndim+2,nc)*oneonrho)

            ! Gravity predictor step
            qin(i,j,k,2,nc) = qin(i,j,k,2,nc) + gravin(i,j,k,1,nc)*dtxhalf
#if NDIM>1
            qin(i,j,k,3,nc) = qin(i,j,k,3,nc) + gravin(i,j,k,2,nc)*dtxhalf
#endif
#if NDIM>2
            qin(i,j,k,4,nc) = qin(i,j,k,4,nc) + gravin(i,j,k,3,nc)*dtxhalf
#endif
         end do
      end do
   end do
  end do
!$acc end parallel loop

#if NVAR > NDIM + 2
  ! Passive scalar
!$acc parallel loop collapse(5) gang vector
  do nc=1,ncube
   do n = ndim+3, nvar
      do k = -1, nxp+2
         do j = -1, nxp+2
            do i = -1, nxp+2
               oneonrho = 1.d0/qin(i,j,k,1,nc)
               qin(i,j,k,n,nc) = uin(i,j,k,n,nc)*oneonrho
            end do
         end do
      end do
   end do
  end do
!$acc end parallel loop 
#endif

!$acc end data

end subroutine ctoprim_gpu_2d
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine uslope_gpu_2d(dx,dt,nxp,ncube)
  use amr_parameters
  use hydro_parameters
  use acc_commons, only: qin,dq
  use const
  implicit none

  integer::nxp,ncube
  real(dp)::dx,dt

  ! local arrays
  integer::i, j, k, l, n, nc, ivar, idim
  real(dp)::dsgn, dlim, dcen, dlft, drgt, slop
  real(dp)::dfll,dflm,dflr,dfmdfmm,dfmr,dfrdfrm,dfrr
  real(dp)::dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl
  real(dp)::dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm
  real(dp)::dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr
  real(dp)::vmin,vmax,dfx,dfy,dfz,dff

  if(slope_type==0)then
     !$acc parallel loop collapse(6) gang vector present(dq) copyin(zero)
     do n=1,ncube
     do idim=1,ndim
     do ivar=1,nvar
     do k=-1,nxp+2
     do j=-1,nxp+2
     do i=-1,nxp+2
         dq(i,j,k,ivar,idim,n)=zero
     end do
     end do
     end do
     end do
     end do
     end do
     !$acc end parallel loop
     return
  end if

!$acc data present(qin,dq)
#if NDIM==1
  do n = 1, nvar
     do k = 0, nxp+1
        do j = 0, nxp+1
           do i = 0, nxp+1
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
        do k = 0, nxp+1
           do j = 0, nxp+1
              do i = 0, nxp+1
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
        do k = 0, nxp+1
           do j = 0, nxp+1
              do i = 0, nxp+1
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
!$acc parallel loop collapse(5) gang vector !async(3)
  do nc=1,ncube
     do n = 1, nvar
        do k = 0, nxp+1
           do j = 0, nxp+1
              do i = 0, nxp+1
                 ! slopes in first coordinate direction
                 dlft = qin(i  ,j,k,n,nc) - qin(i-1,j,k,n,nc)
                 drgt = qin(i+1,j,k,n,nc) - qin(i  ,j,k,n,nc)
                 if((dlft*drgt)<=zero) then
                    dq(i,j,k,n,1,nc) = zero
                 else if(dlft>0) then
                    dq(i,j,k,n,1,nc) = min(dlft,drgt)
                 else
                    dq(i,j,k,n,1,nc) = max(dlft,drgt)
                 end if
                 ! slopes in second coordinate direction
                 dlft = qin(i,j  ,k,n,nc) - qin(i,j-1,k,n,nc)
                 drgt = qin(i,j+1,k,n,nc) - qin(i,j  ,k,n,nc)
                 if((dlft*drgt)<=zero) then
                    dq(i,j,k,n,2,nc) = zero
                 else if(dlft>0) then
                    dq(i,j,k,n,2,nc) = min(dlft,drgt)
                 else
                    dq(i,j,k,n,2,nc) = max(dlft,drgt)
                 end if
                 ! slopes in third coordinate direction
                 dlft = qin(i,j,k  ,n,nc) - qin(i,j,k-1,n,nc)
                 drgt = qin(i,j,k+1,n,nc) - qin(i,j,k  ,n,nc)
                 if((dlft*drgt)<=zero) then
                    dq(i,j,k,n,3,nc) = zero
                 else if(dlft>0) then
                    dq(i,j,k,n,3,nc) = min(dlft,drgt)
                 else
                    dq(i,j,k,n,3,nc) = max(dlft,drgt)
                 end if
              end do
           end do
        end do
     end do
  end do
!$acc end parallel loop
  
  else if(slope_type==2)then ! moncen
!$acc parallel loop collapse(5) gang vector !async(3)
     do nc=1,ncube
     do n = 1, nvar
        do k = 0, nxp+1
           do j = 0, nxp+1
              do i = 0, nxp+1
                 ! slopes in first coordinate direction
                 dlft = slope_type*(qin(i  ,j,k,n,nc) - qin(i-1,j,k,n,nc))
                 drgt = slope_type*(qin(i+1,j,k,n,nc) - qin(i  ,j,k,n,nc))
                 dcen = half*(dlft+drgt)/slope_type
                 dsgn = sign(one, dcen)
                 slop = min(abs(dlft),abs(drgt))
                 dlim = slop
                 if((dlft*drgt)<=zero)dlim=zero
                 dq(i,j,k,n,1,nc) = dsgn*min(dlim,abs(dcen))
                 ! slopes in second coordinate direction
                 dlft = slope_type*(qin(i,j  ,k,n,nc) - qin(i,j-1,k,n,nc))
                 drgt = slope_type*(qin(i,j+1,k,n,nc) - qin(i,j  ,k,n,nc))
                 dcen = half*(dlft+drgt)/slope_type
                 dsgn = sign(one,dcen)
                 slop = min(abs(dlft),abs(drgt))
                 dlim = slop
                 if((dlft*drgt)<=zero)dlim=zero
                 dq(i,j,k,n,2,nc) = dsgn*min(dlim,abs(dcen))
                 ! slopes in third coordinate direction
                 dlft = slope_type*(qin(i,j,k  ,n,nc) - qin(i,j,k-1,n,nc))
                 drgt = slope_type*(qin(i,j,k+1,n,nc) - qin(i,j,k  ,n,nc))
                 dcen = half*(dlft+drgt)/slope_type
                 dsgn = sign(one,dcen)
                 slop = min(abs(dlft),abs(drgt))
                 dlim = slop
                 if((dlft*drgt)<=zero)dlim=zero
                 dq(i,j,k,n,3,nc) = dsgn*min(dlim,abs(dcen))
              end do
           end do
        end do
     end do
     end do
!$acc end parallel loop
     
  else if(slope_type==3)then ! positivity preserving 3d unsplit slope
!$acc parallel loop collapse(5) gang vector !async(3)
     do nc=1,ncube
     do n = 1, nvar
        do k = 0, nxp+1
           do j = 0, nxp+1
              do i = 0, nxp+1
                 dflll = qin(i-1,j-1,k-1,n,nc)-qin(i,j,k,n,nc)
                 dflml = qin(i-1,j  ,k-1,n,nc)-qin(i,j,k,n,nc)
                 dflrl = qin(i-1,j+1,k-1,n,nc)-qin(i,j,k,n,nc)
                 dfmll = qin(i  ,j-1,k-1,n,nc)-qin(i,j,k,n,nc)
                 dfmml = qin(i  ,j  ,k-1,n,nc)-qin(i,j,k,n,nc)
                 dfmrl = qin(i  ,j+1,k-1,n,nc)-qin(i,j,k,n,nc)
                 dfrll = qin(i+1,j-1,k-1,n,nc)-qin(i,j,k,n,nc)
                 dfrml = qin(i+1,j  ,k-1,n,nc)-qin(i,j,k,n,nc)
                 dfrrl = qin(i+1,j+1,k-1,n,nc)-qin(i,j,k,n,nc)
                 
                 dfllm = qin(i-1,j-1,k  ,n,nc)-qin(i,j,k,n,nc)
                 dflmm = qin(i-1,j  ,k  ,n,nc)-qin(i,j,k,n,nc)
                 dflrm = qin(i-1,j+1,k  ,n,nc)-qin(i,j,k,n,nc)
                 dfmlm = qin(i  ,j-1,k  ,n,nc)-qin(i,j,k,n,nc)
                 dfmmm = qin(i  ,j  ,k  ,n,nc)-qin(i,j,k,n,nc)
                 dfmrm = qin(i  ,j+1,k  ,n,nc)-qin(i,j,k,n,nc)
                 dfrlm = qin(i+1,j-1,k  ,n,nc)-qin(i,j,k,n,nc)
                 dfrmm = qin(i+1,j  ,k  ,n,nc)-qin(i,j,k,n,nc)
                 dfrrm = qin(i+1,j+1,k  ,n,nc)-qin(i,j,k,n,nc)
                 
                 dfllr = qin(i-1,j-1,k+1,n,nc)-qin(i,j,k,n,nc)
                 dflmr = qin(i-1,j  ,k+1,n,nc)-qin(i,j,k,n,nc)
                 dflrr = qin(i-1,j+1,k+1,n,nc)-qin(i,j,k,n,nc)
                 dfmlr = qin(i  ,j-1,k+1,n,nc)-qin(i,j,k,n,nc)
                 dfmmr = qin(i  ,j  ,k+1,n,nc)-qin(i,j,k,n,nc)
                 dfmrr = qin(i  ,j+1,k+1,n,nc)-qin(i,j,k,n,nc)
                 dfrlr = qin(i+1,j-1,k+1,n,nc)-qin(i,j,k,n,nc)
                 dfrmr = qin(i+1,j  ,k+1,n,nc)-qin(i,j,k,n,nc)
                 dfrrr = qin(i+1,j+1,k+1,n,nc)-qin(i,j,k,n,nc)
                 
                 
                 vmin = min(dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl, &
                      &     dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm, &
                      &     dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr)
                 vmax = max(dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl, &
                      &     dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm, &
                      &     dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr)
                 
                 dfx  = half*(qin(i+1,j,k,n,nc)-qin(i-1,j,k,n,nc))
                 dfy  = half*(qin(i,j+1,k,n,nc)-qin(i,j-1,k,n,nc))
                 dfz  = half*(qin(i,j,k+1,n,nc)-qin(i,j,k-1,n,nc))
                 dff  = half*(abs(dfx)+abs(dfy)+abs(dfz))
                 
                 if(dff>zero)then
                    slop = min(one,min(abs(vmin),abs(vmax))/dff)
                 else
                    slop = one
                 endif
                 
                 dlim = slop
                 
                 dq(i,j,k,n,1,nc) = dlim*dfx
                 dq(i,j,k,n,2,nc) = dlim*dfy
                 dq(i,j,k,n,3,nc) = dlim*dfz
                 
              end do
           end do
        end do
     end do
     end do
!$acc end parallel loop
  else
     write(*,*)'Unknown slope type'
     stop
  endif     
#endif

!$acc end data
  
end subroutine uslope_gpu_2d
