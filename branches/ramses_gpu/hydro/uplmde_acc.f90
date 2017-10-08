!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine tracex_gpu_2d(q,dq,c,qm,qp,dx,dt,ngrid,nxp)
  use amr_parameters
  use hydro_parameters
  use acc_parameters
  use const
  implicit none

  integer::ngrid,nxp
  real(dp)::dx, dt  

  real(dp),dimension(-1:nxp+2,-1:nxp+2,-1:nxp+2,1:nvar)::q  
  real(dp),dimension(-1:nxp+2,-1:nxp+2,-1:nxp+2,1:nvar,1:ndim)::dq 
  real(dp),dimension(-1:nxp+2,-1:nxp+2,-1:nxp+2,1:nvar,1:ndim)::qm 
  real(dp),dimension(-1:nxp+2,-1:nxp+2,-1:nxp+2,1:nvar,1:ndim)::qp 
  real(dp),dimension(-1:nxp+2,-1:nxp+2,-1:nxp+2)::c  

  ! Local variables
  integer ::ilo,ihi,jlo,jhi,klo,khi
  integer ::i, j, k, l, n
  integer ::ir, iu, ip
  real(dp)::dtdx,project_out
  real(dp)::cc, ccc, csq, r, u, p, a
  real(dp)::drx, dux, dpx, dax
  real(dp)::alpham, alphap, alpha0r
  real(dp)::spminus, spplus, spzero
  real(dp)::apright, amright, azrright, azaright
  real(dp)::apleft,  amleft,  azrleft,  azaleft
  
  dtdx = dt/dx
  ilo=MIN(1,iu10+1); ihi=MAX(1,iu20-1)
  jlo=MIN(1,ju10+1); jhi=MAX(1,ju20-1)
  klo=MIN(1,ku10+1); khi=MAX(1,ku20-1)
  ir=1; iu=2; ip=3
  project_out=one !zero

  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           do l = 1, ngrid

              ! Cell centered values
              cc  = c (i,j,k)
              r   = q (i,j,k,ir)
              u   = q (i,j,k,iu)
              p   = q (i,j,k,ip)
              csq = gamma*p/r

              ! TVD slopes in X direction
              drx = dq(i,j,k,ir,1)
              dux = dq(i,j,k,iu,1)
              dpx = dq(i,j,k,ip,1)
              
              ! Supersonic fix for high-velocity gradients
              ccc = cc
              if(ABS(dux) > three*cc)ccc=zero

              ! Characteristic analysis along X direction
              alpham  = half*(dpx/csq - dux*r/cc)
              alphap  = half*(dpx/csq + dux*r/cc)
              alpha0r = drx - dpx/csq

              ! Right state
              spminus = (u-ccc)*dtdx
              spplus  = (u+ccc)*dtdx
              spzero  = (u    )*dtdx
              if((u+ccc)>zero)spplus =-project_out
              if((u-ccc)>zero)spminus=-project_out
              if( u     >zero)spzero =-project_out

              apright  = half*(-one-spplus )*alphap
              amright  = half*(-one-spminus)*alpham 
              azrright = half*(-one-spzero )*alpha0r

              qp(i,j,k,ir,1) = r + (apright+amright+azrright)
              qp(i,j,k,iu,1) = u + (apright-amright         )*cc/r
              qp(i,j,k,ip,1) = p + (apright+amright         )*csq
              qp(i,j,k,ir,1) = max(smallr,qp(i,j,k,ir,1))

              ! Left state
              spminus = (u-ccc)*dtdx
              spplus  = (u+ccc)*dtdx
              spzero  = (u    )*dtdx
              if((u+ccc)<=zero)spplus =+project_out
              if((u-ccc)<=zero)spminus=+project_out
              if( u     <=zero)spzero =+project_out

              apleft   = half*(+one-spplus )*alphap
              amleft   = half*(+one-spminus)*alpham
              azrleft  = half*(+one-spzero )*alpha0r

              qm(i,j,k,ir,1) = r + (apleft+amleft+azrleft)
              qm(i,j,k,iu,1) = u + (apleft-amleft        )*cc/r
              qm(i,j,k,ip,1) = p + (apleft+amleft        )*csq
              qm(i,j,k,ir,1) = max(smallr, qm(i,j,k,ir,1))

           end do
        end do
     end do
  end do

#if NVAR > NDIM + 2
  do n = ndim+3, nvar
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              do l = 1, ngrid
                 a   =  q(i,j,k,n)    ! Cell centered values
                 u   =  q(i,j,k,iu)
                 dax = dq(i,j,k,n,1)  ! TVD slopes
                 
                 ! Right state
                 spzero=(u    )*dtdx
                 if(u>zero)spzero=-project_out
                 azaright = half*(-one-spzero )*dax
                 qp(i,j,k,n,1) = a + azaright

                 ! Left state
                 spzero=(u    )*dtdx
                 if(u<=zero)spzero=+project_out
                 azaleft = half*(+one-spzero)*dax
                 qm(i,j,k,n,1) = a + azaleft
              end do
           end do
        end do
     end do
  end do
#endif

end subroutine tracex_gpu_2d
!###########################################################
!###########################################################
!###########################################################
!###########################################################
#if NDIM>1
subroutine tracexy_gpu_2d(q,dq,c,qm,qp,dx,dy,dt,ngrid,nxp)
  use amr_parameters
  use hydro_parameters
  use acc_parameters
  use const
  implicit none

  integer ::ngrid,nxp
  real(dp)::dx, dy, dt

  real(dp),dimension(-1:nxp+2,-1:nxp+2,-1:nxp+2,1:nvar)::q  
  real(dp),dimension(-1:nxp+2,-1:nxp+2,-1:nxp+2,1:nvar,1:ndim)::dq 
  real(dp),dimension(-1:nxp+2,-1:nxp+2,-1:nxp+2,1:nvar,1:ndim)::qm 
  real(dp),dimension(-1:nxp+2,-1:nxp+2,-1:nxp+2,1:nvar,1:ndim)::qp 
  real(dp),dimension(-1:nxp+2,-1:nxp+2,-1:nxp+2)::c  

  ! declare local variables
  integer ::ilo,ihi,jlo,jhi,klo,khi
  integer ::i, j, k, l, n
  integer ::ir, iu, iv, ip
  real(dp)::dtdx,dtdy,project_out
  real(dp)::cc, ccc, csq, r, u, v, p, a
  real(dp)::drx, dux, dvx, dpx, dax
  real(dp)::dry, duy, dvy, dpy, day
  real(dp)::alpham, alphap, alpha0r, alpha0u, alpha0v
  real(dp)::spminus, spplus, spzero
  real(dp)::apright, amright, azrright, azuright, azvright, azaright
  real(dp)::apleft,  amleft,  azrleft,  azuleft,  azvleft,  azaleft
  real(dp)::srx,sux,svx,spx,sax
  real(dp)::sry,suy,svy,spy,say
    
  dtdx = dt/dx; dtdy = dt/dy
  ilo=MIN(1,iu10+1); ihi=MAX(1,iu20-1)
  jlo=MIN(1,ju10+1); jhi=MAX(1,ju20-1)
  klo=MIN(1,ku10+1); khi=MAX(1,ku20-1)
  ir=1; iu=2; iv=3; ip=4
  project_out=one !zero

  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           do l = 1, ngrid

              ! cell centered values
              cc  = c  (i,j,k)
              r   = q  (i,j,k,ir)
              u   = q  (i,j,k,iu)
              v   = q  (i,j,k,iv)
              p   = q  (i,j,k,ip)
              csq = gamma*p/r

              ! TVD slopes in X and Y directions
              drx = dq(i,j,k,ir,1)
              dux = dq(i,j,k,iu,1)
              dvx = dq(i,j,k,iv,1)
              dpx = dq(i,j,k,ip,1)
              
              dry = dq(i,j,k,ir,2)
              duy = dq(i,j,k,iu,2)
              dvy = dq(i,j,k,iv,2)
              dpy = dq(i,j,k,ip,2)
              
              ! Transverse derivatives
              srx = half*dtdy*(-v*dry - (dvy)*r      )
              sux = half*dtdy*(-v*duy                )
              svx = half*dtdy*(-v*dvy - (dpy)/r      )
              spx = half*dtdy*(-v*dpy - (dvy)*gamma*p)

              sry = half*dtdx*(-u*drx - (dux)*r      )
              suy = half*dtdx*(-u*dux - (dpx)/r      )
              svy = half*dtdx*(-u*dvx                )
              spy = half*dtdx*(-u*dpx - (dux)*gamma*p)

              ! Characteristic analysis along X direction
              alpham  = half*(dpx/csq - dux*r/cc)
              alphap  = half*(dpx/csq + dux*r/cc)
              alpha0r = drx - dpx/csq
              alpha0v = dvx

              ! Supersonic fix for high-velocity gradients
              ccc = cc
              if(ABS(dux) > three*cc)ccc=zero

              ! Right state
              spminus = (u-ccc)*dtdx
              spplus  = (u+ccc)*dtdx
              spzero  = (u    )*dtdx
              if((u+ccc)>zero)spplus =-project_out
              if((u-ccc)>zero)spminus=-project_out
              if( u     >zero)spzero =-project_out

              apright  = half*(-one-spplus )*alphap 
              amright  = half*(-one-spminus)*alpham 
              azrright = half*(-one-spzero )*alpha0r
              azvright = half*(-one-spzero )*alpha0v

              qp(i,j,k,ir,1) = r + (apright+amright+azrright)     +srx
              qp(i,j,k,iu,1) = u + (apright-amright         )*cc/r+sux
              qp(i,j,k,ip,1) = p + (apright+amright         )*csq +spx
              qp(i,j,k,iv,1) = v + (                azvright)     +svx
              qp(i,j,k,ir,1) = max(smallr,qp(i,j,k,ir,1))

              ! Left state
              spminus = (u-ccc)*dtdx
              spplus  = (u+ccc)*dtdx
              spzero  = (u    )*dtdx
              if((u+ccc)<=zero)spplus =+project_out
              if((u-ccc)<=zero)spminus=+project_out
              if( u     <=zero)spzero =+project_out

              apleft   = half*(+one-spplus )*alphap 
              amleft   = half*(+one-spminus)*alpham 
              azrleft  = half*(+one-spzero )*alpha0r
              azvleft  = half*(+one-spzero )*alpha0v
              
              qm(i,j,k,ir,1) = r + (apleft+amleft+azrleft)     +srx
              qm(i,j,k,iu,1) = u + (apleft-amleft        )*cc/r+sux
              qm(i,j,k,ip,1) = p + (apleft+amleft        )*csq +spx
              qm(i,j,k,iv,1) = v + (              azvleft)     +svx
              qm(i,j,k,ir,1) = max(smallr, qm(i,j,k,ir,1))

              ! Characteristic analysis along Y direction
              alpham  = half*(dpy/csq - dvy*r/cc)
              alphap  = half*(dpy/csq + dvy*r/cc)
              alpha0r = dry - dpy/csq
              alpha0u = duy

              ! Supersonic fix for high-velocity gradients
              ccc = cc
              if(ABS(dvy) > three*cc)ccc=zero

              ! Top state
              spminus = (v-ccc)*dtdy
              spplus  = (v+ccc)*dtdy
              spzero  = (v    )*dtdy
              if((v+ccc)>zero)spplus =-project_out
              if((v-ccc)>zero)spminus=-project_out
              if( v     >zero)spzero =-project_out

              apright  = half*(-one-spplus )*alphap
              amright  = half*(-one-spminus)*alpham 
              azrright = half*(-one-spzero )*alpha0r
              azuright = half*(-one-spzero )*alpha0u

              qp(i,j,k,ir,2) = r + (apright+amright+azrright)     +sry
              qp(i,j,k,iv,2) = v + (apright-amright         )*cc/r+svy
              qp(i,j,k,ip,2) = p + (apright+amright         )*csq +spy
              qp(i,j,k,iu,2) = u + (                azuright)     +suy
              qp(i,j,k,ir,2) = max(smallr,qp(i,j,k,ir,2))

              ! Bottom state
              spminus = (v-ccc)*dtdy
              spplus  = (v+ccc)*dtdy
              spzero  = (v    )*dtdy
              if((v+ccc)<=zero)spplus =+project_out
              if((v-ccc)<=zero)spminus=+project_out
              if( v     <=zero)spzero =+project_out

              apleft   = half*(+one-spplus )*alphap
              amleft   = half*(+one-spminus)*alpham
              azrleft  = half*(+one-spzero )*alpha0r
              azuleft  = half*(+one-spzero )*alpha0u

              qm(i,j,k,ir,2) = r + (apleft+amleft+azrleft)     +sry
              qm(i,j,k,iv,2) = v + (apleft-amleft        )*cc/r+svy
              qm(i,j,k,ip,2) = p + (apleft+amleft        )*csq +spy
              qm(i,j,k,iu,2) = u + (              azuleft)     +suy
              qm(i,j,k,ir,2) = max(smallr, qm(i,j,k,ir,2))

           end do
        end do
     end do
  end do

#if NVAR > NDIM + 2
  ! Passive scalars
  do n = ndim+3, nvar
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              do l = 1, ngrid
                 a = q(i,j,k,n)     ! Cell centered values
                 u = q(i,j,k,iu)
                 v = q(i,j,k,iv)
                 dax = dq(i,j,k,n,1)    ! TVD slopes
                 day = dq(i,j,k,n,2)
                 sax = half*dtdy*(-v*day) ! Transverse
                 say = half*dtdx*(-u*dax) ! derivatives

                 ! Right state
                 spzero=(u    )*dtdx
                 if(u>zero)spzero=-project_out
                 azaright = half*(-one-spzero )*dax
                 qp(i,j,k,n,1) = a + azaright + sax
                 
                 ! Left state
                 spzero=(u    )*dtdx
                 if(u<=zero)spzero=+project_out
                 azaleft = half*(+one-spzero )*dax
                 qm(i,j,k,n,1) = a + azaleft + sax
                 
                 ! Top state
                 spzero=(v    )*dtdy
                 if(v>zero)spzero=-project_out
                 azaright = half*(-one-spzero )*day
                 qp(i,j,k,n,2) = a + azaright + say
                 
                 ! Bottom state
                 spzero=(v    )*dtdy
                 if(v<=zero)spzero=+project_out
                 azaleft = half*(+one-spzero )*day
                 qm(i,j,k,n,2) = a + azaleft + say
                 
              end do
           end do
        end do
     end do
  end do
#endif

end subroutine tracexy_gpu_2d
#endif
!###########################################################
!###########################################################
!###########################################################
!###########################################################
#if NDIM>2
subroutine tracexyz_gpu_2d(q,dq,c,qm,qp,dx,dy,dz,dt,ngrid,nxp)
  use amr_parameters
  use hydro_parameters
  use acc_parameters
  use const
  implicit none

  integer ::ngrid,nxp
  real(dp)::dx,dy,dz, dt

  real(dp),dimension(-1:nxp+2,-1:nxp+2,-1:nxp+2,1:nvar)::q  
  real(dp),dimension(-1:nxp+2,-1:nxp+2,-1:nxp+2,1:nvar,1:ndim)::dq 
  real(dp),dimension(-1:nxp+2,-1:nxp+2,-1:nxp+2,1:nvar,1:ndim)::qm 
  real(dp),dimension(-1:nxp+2,-1:nxp+2,-1:nxp+2,1:nvar,1:ndim)::qp 
  real(dp),dimension(-1:nxp+2,-1:nxp+2,-1:nxp+2)::c  

  ! declare local variables
  integer ::ilo,ihi,jlo,jhi,klo,khi
  integer ::i, j, k, l, n
  integer ::ir, iu, iv, iw, ip
  real(dp)::dtdx,dtdy,dtdz,project_out
  real(dp)::cc, ccc, csq, r, u, v, w, p, a
  real(dp)::drx, dux, dvx, dwx, dpx, dax
  real(dp)::dry, duy, dvy, dwy, dpy, day
  real(dp)::drz, duz, dvz, dwz, dpz, daz
  real(dp)::alpham, alphap, alpha0r, alpha0u, alpha0v, alpha0w
  real(dp)::spminus, spplus, spzero
  real(dp)::apright, amright, azrright, azuright, azvright, azwright, azaright
  real(dp)::apleft,  amleft,  azrleft,  azuleft,  azvleft,  azwleft,  azaleft
  real(dp)::srx,sux,svx,swx,spx,sax
  real(dp)::sry,suy,svy,swy,spy,say
  real(dp)::srz,suz,svz,swz,spz,saz
    
  dtdx = dt/dx; dtdy = dt/dy; dtdz = dt/dz
  ilo=MIN(1,iu10+1); ihi=MAX(1,iu20-1)
  jlo=MIN(1,ju10+1); jhi=MAX(1,ju20-1)
  klo=MIN(1,ku10+1); khi=MAX(1,ku20-1)
  ir=1; iu=2; iv=3; iw=4; ip=5
  project_out=one !zero
  
  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           do l = 1, ngrid
  
              ! Cell centered values
              cc  = c  (i,j,k)
              r   = q  (i,j,k,ir)
              u   = q  (i,j,k,iu)
              v   = q  (i,j,k,iv)
              w   = q  (i,j,k,iw)
              p   = q  (i,j,k,ip)
              csq = gamma*p/r

              ! TVD slopes in all 3 directions
              drx = dq(i,j,k,ir,1)
              dux = dq(i,j,k,iu,1)
              dvx = dq(i,j,k,iv,1)
              dwx = dq(i,j,k,iw,1)
              dpx = dq(i,j,k,ip,1)
              
              dry = dq(i,j,k,ir,2)
              duy = dq(i,j,k,iu,2)
              dvy = dq(i,j,k,iv,2)
              dwy = dq(i,j,k,iw,2)
              dpy = dq(i,j,k,ip,2)
              
              drz = dq(i,j,k,ir,3)
              duz = dq(i,j,k,iu,3)
              dvz = dq(i,j,k,iv,3)
              dwz = dq(i,j,k,iw,3)
              dpz = dq(i,j,k,ip,3)

              ! Transverse derivatives
              srx = half*dtdx*(-v*dry-w*drz - (dvy+dwz)*r      )
              spx = half*dtdx*(-v*dpy-w*dpz - (dvy+dwz)*gamma*p)
              sux = half*dtdx*(-v*duy-w*duz                    )
              svx = half*dtdx*(-v*dvy-w*dvz - (dpy)/r          )
              swx = half*dtdx*(-v*dwy-w*dwz - (dpz)/r          )

              sry = half*dtdx*(-u*drx-w*drz - (dux+dwz)*r      )
              spy = half*dtdx*(-u*dpx-w*dpz - (dux+dwz)*gamma*p)
              suy = half*dtdx*(-u*dux-w*duz - (dpx)/r          )
              svy = half*dtdx*(-u*dvx-w*dvz                    )
              swy = half*dtdx*(-u*dwx-w*dwz - (dpz)/r          )

              srz = half*dtdx*(-v*dry-u*drx - (dvy+dux)*r      )
              spz = half*dtdx*(-v*dpy-u*dpx - (dvy+dux)*gamma*p)
              suz = half*dtdx*(-v*duy-u*dux - (dpx)/r          )
              svz = half*dtdx*(-v*dvy-u*dvx - (dpy)/r          )
              swz = half*dtdx*(-v*dwy-u*dwx                    )

              ! Characteristic analysis along X direction
              alpham  = half*(dpx/csq - dux*r/cc)
              alphap  = half*(dpx/csq + dux*r/cc)
              alpha0r = drx - dpx/csq
              alpha0v = dvx
              alpha0w = dwx

              ! Supersonic fix for high-velocity gradients
              ccc = cc
              if(ABS(dux) > three*cc)ccc=zero

              ! Right state
              spminus = (u-ccc)*dtdx
              spplus  = (u+ccc)*dtdx
              spzero  = (u    )*dtdx
              if((u+ccc)>zero)spplus =-project_out
              if((u-ccc)>zero)spminus=-project_out
              if( u     >zero)spzero =-project_out

              apright  = half*(-one-spplus )*alphap
              amright  = half*(-one-spminus)*alpham
              azrright = half*(-one-spzero )*alpha0r
              azvright = half*(-one-spzero )*alpha0v
              azwright = half*(-one-spzero )*alpha0w

              qp(i,j,k,ir,1) = r + (apright+amright+azrright)     +srx
              qp(i,j,k,iu,1) = u + (apright-amright         )*cc/r+sux
              qp(i,j,k,ip,1) = p + (apright+amright         )*csq +spx
              qp(i,j,k,iv,1) = v + (                azvright)     +svx
              qp(i,j,k,iw,1) = w + (                azwright)     +swx
              qp(i,j,k,ir,1) = max(smallr,qp(i,j,k,ir,1))

              ! Left state
              spminus = (u-ccc)*dtdx
              spplus  = (u+ccc)*dtdx
              spzero  = (u    )*dtdx
              if((u+ccc)<=zero)spplus =+project_out
              if((u-ccc)<=zero)spminus=+project_out
              if( u     <=zero)spzero =+project_out

              apleft   = half*(+one-spplus )*alphap
              amleft   = half*(+one-spminus)*alpham
              azrleft  = half*(+one-spzero )*alpha0r
              azvleft  = half*(+one-spzero )*alpha0v
              azwleft  = half*(+one-spzero )*alpha0w

              qm(i,j,k,ir,1) = r + (apleft+amleft+azrleft)     +srx
              qm(i,j,k,iu,1) = u + (apleft-amleft        )*cc/r+sux
              qm(i,j,k,ip,1) = p + (apleft+amleft        )*csq +spx
              qm(i,j,k,iv,1) = v + (              azvleft)     +svx
              qm(i,j,k,iw,1) = w + (              azwleft)     +swx
              qm(i,j,k,ir,1) = max(smallr, qm(i,j,k,ir,1))

              ! Characteristic analysis along Y direction
              alpham  = half*(dpy/csq - dvy*r/cc)
              alphap  = half*(dpy/csq + dvy*r/cc)
              alpha0r = dry - dpy/csq
              alpha0u = duy
              alpha0w = dwy

              ! Supersonic fix for high-velocity gradients
              ccc = cc
              if(ABS(dvy) > three*cc)ccc=zero

              ! Top state
              spminus = (v-ccc)*dtdy
              spplus  = (v+ccc)*dtdy
              spzero  = (v    )*dtdy
              if((v+ccc)>zero)spplus =-project_out
              if((v-ccc)>zero)spminus=-project_out
              if( v     >zero)spzero =-project_out

              apright  = half*(-one-spplus )*alphap 
              amright  = half*(-one-spminus)*alpham 
              azrright = half*(-one-spzero )*alpha0r
              azuright = half*(-one-spzero )*alpha0u
              azwright = half*(-one-spzero )*alpha0w

              qp(i,j,k,ir,2) = r + (apright+amright+azrright)     +sry
              qp(i,j,k,iv,2) = v + (apright-amright         )*cc/r+svy
              qp(i,j,k,ip,2) = p + (apright+amright         )*csq +spy
              qp(i,j,k,iu,2) = u + (                azuright)     +suy
              qp(i,j,k,iw,2) = w + (                azwright)     +swy
              qp(i,j,k,ir,2) = max(smallr,qp(i,j,k,ir,2))

              ! Bottom state
              spminus = (v-ccc)*dtdy
              spplus  = (v+ccc)*dtdy
              spzero  = (v    )*dtdy
              if((v+ccc)<=zero)spplus =+project_out
              if((v-ccc)<=zero)spminus=+project_out
              if( v     <=zero)spzero =+project_out

              apleft   = half*(+one-spplus )*alphap 
              amleft   = half*(+one-spminus)*alpham 
              azrleft  = half*(+one-spzero )*alpha0r
              azuleft  = half*(+one-spzero )*alpha0u
              azwleft  = half*(+one-spzero )*alpha0w

              qm(i,j,k,ir,2) = r + (apleft+amleft+azrleft)     +sry
              qm(i,j,k,iv,2) = v + (apleft-amleft        )*cc/r+svy
              qm(i,j,k,ip,2) = p + (apleft+amleft        )*csq +spy
              qm(i,j,k,iu,2) = u + (              azuleft)     +suy
              qm(i,j,k,iw,2) = w + (              azwleft)     +swy
              qm(i,j,k,ir,2) = max(smallr, qm(i,j,k,ir,2))

              ! Characteristic analysis along Z direction
              alpham  = half*(dpz/csq - dwz*r/cc)
              alphap  = half*(dpz/csq + dwz*r/cc)
              alpha0r = drz - dpz/csq
              alpha0u = duz
              alpha0v = dvz

              ! Supersonic fix for high-velocity gradients
              ccc = cc
              if(ABS(dwz) > three*cc)ccc=zero

              ! Front state
              spminus = (w-ccc)*dtdz
              spplus  = (w+ccc)*dtdz
              spzero  = (w    )*dtdz
              if((w+ccc)>zero)spplus =-project_out
              if((w-ccc)>zero)spminus=-project_out
              if( w     >zero)spzero =-project_out

              apright  = half*(-one-spplus )*alphap 
              amright  = half*(-one-spminus)*alpham 
              azrright = half*(-one-spzero )*alpha0r
              azuright = half*(-one-spzero )*alpha0u
              azvright = half*(-one-spzero )*alpha0v

              qp(i,j,k,ir,3) = r + (apright+amright+azrright)     +srz
              qp(i,j,k,iw,3) = w + (apright-amright         )*cc/r+swz
              qp(i,j,k,ip,3) = p + (apright+amright         )*csq +spz
              qp(i,j,k,iu,3) = u + (                azuright)     +suz
              qp(i,j,k,iv,3) = v + (                azvright)     +svz
              qp(i,j,k,ir,3) = max(smallr,qp(i,j,k,ir,3))

              ! Back state
              spminus = (w-ccc)*dtdz
              spplus  = (w+ccc)*dtdz
              spzero  = (w    )*dtdz
              if((w+ccc)<=zero)spplus =+project_out
              if((w-ccc)<=zero)spminus=+project_out
              if( w     <=zero)spzero =+project_out

              apleft   = half*(+one-spplus )*alphap 
              amleft   = half*(+one-spminus)*alpham 
              azrleft  = half*(+one-spzero )*alpha0r
              azuleft  = half*(+one-spzero )*alpha0u
              azvleft  = half*(+one-spzero )*alpha0v

              qm(i,j,k,ir,3) = r + (apleft+amleft+azrleft)     +srz
              qm(i,j,k,iw,3) = w + (apleft-amleft        )*cc/r+swz
              qm(i,j,k,ip,3) = p + (apleft+amleft        )*csq +spz
              qm(i,j,k,iu,3) = u + (              azuleft)     +suz
              qm(i,j,k,iv,3) = v + (              azvleft)     +svz
              qm(i,j,k,ir,3) = max(smallr, qm(i,j,k,ir,3))
           end do

        end do
     end do
  end do

#if NVAR > NDIM + 2
  ! Passive scalars
  do n = ndim+3, nvar
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              do l = 1, ngrid
                 a   =  q(i,j,k,n)    ! Cell centered values
                 u   =  q(i,j,k,iu)
                 v   =  q(i,j,k,iv)
                 w   =  q(i,j,k,iw)
                 dax = dq(i,j,k,n,1)  ! TVD slopes
                 day = dq(i,j,k,n,2)
                 daz = dq(i,j,k,n,3)
                 sax = half*dtdx*(-v*day-w*daz) ! Transverse
                 say = half*dtdx*(-u*dax-w*daz) ! derivatives
                 saz = half*dtdx*(-v*day-u*dax) ! 

                 
                 ! Right state
                 spzero = (u    )*dtdx
                 if(u>zero)spzero=-project_out
                 azaright = half*(-one-spzero )*dax
                 qp(i,j,k,n,1) = a + azaright + sax
                 
                 ! Left state
                 spzero = (u    )*dtdx
                 if(u<=zero)spzero=+project_out
                 azaleft = half*(+one-spzero )*dax
                 qm(i,j,k,n,1) = a + azaleft + sax

                 ! Top state
                 spzero = (v    )*dtdy
                 if(v>zero)spzero=-project_out
                 azaright = half*(-one-spzero )*day
                 qp(i,j,k,n,2) = a + azaright + say
                 
                 ! Bottom state
                 spzero = (v    )*dtdy
                 if(v<=zero)spzero=+project_out
                 azaleft = half*(+one-spzero )*day
                 qm(i,j,k,n,2) = a + azaleft + say

                 ! Front state
                 spzero = (w    )*dtdy
                 if(w>zero)spzero=-project_out
                 azaright = half*(-one-spzero )*daz
                 qp(i,j,k,n,3) = a + azaright + saz
                 
                 ! Back state
                 spzero = (w    )*dtdy
                 if(w<=zero)spzero=+project_out
                 azaleft = half*(+one-spzero )*daz
                 qm(i,j,k,n,3) = a + azaleft + saz
              end do
           end do
        end do
     end do
  end do
#endif

end subroutine tracexyz_gpu_2d
#endif
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmpdivu_gpu_2d(q,div,dx,dy,dz,ngrid,nxp)
  use amr_parameters
  use hydro_parameters
  use acc_parameters
  use const
  implicit none

  integer ::ngrid,nxp
  real(dp)::dx, dy, dz
  real(dp),dimension(-1:nxp+2,-1:nxp+2,-1:nxp+2,1:nvar)::q  
  real(dp),dimension(1:nxp+1,1:nxp+1,1:nxp+1)::div

  integer::i, j, k, l
  real(dp)::factorx, factory, factorz
  real(dp),dimension(1:nvector)::ux, vy, wz

  factorx=half**(ndim-1)/dx
  factory=half**(ndim-1)/dy
  factorz=half**(ndim-1)/dz

  do k = kf10, kf20
     do j = jf10, jf20
        do i = if10, if20

           ux = zero; vy=zero; wz=zero

           if(ndim>0)then
              do l=1, ngrid
                 ux(l)=ux(l)+factorx*(q(i  ,j  ,k  ,2) - q(i-1,j  ,k  ,2))
              end do
           end if

#if NDIM>1           
           if(ndim>1)then
              do l=1, ngrid
                 ux(l)=ux(l)+factorx*(q(i  ,j-1,k  ,2) - q(i-1,j-1,k  ,2))
                 vy(l)=vy(l)+factory*(q(i  ,j  ,k  ,3) - q(i  ,j-1,k  ,3)+&
                      &               q(i-1,j  ,k  ,3) - q(i-1,j-1,k  ,3))
              end do
           end if
#endif
#if NDIM>2
           if(ndim>2)then
              do l=1, ngrid
                 ux(l)=ux(l)+factorx*(q(i  ,j  ,k-1,2) - q(i-1,j  ,k-1,2)+&
                      &               q(i  ,j-1,k-1,2) - q(i-1,j-1,k-1,2))
                 vy(l)=vy(l)+factory*(q(i  ,j  ,k-1,3) - q(i  ,j-1,k-1,3)+&
                      &               q(i-1,j  ,k-1,3) - q(i-1,j-1,k-1,3))
                 wz(l)=wz(l)+factorz*(q(i  ,j  ,k  ,4) - q(i  ,j  ,k-1,4)+&
                      &               q(i  ,j-1,k  ,4) - q(i  ,j-1,k-1,4)+&
                      &               q(i-1,j  ,k  ,4) - q(i-1,j  ,k-1,4)+&
                      &               q(i-1,j-1,k  ,4) - q(i-1,j-1,k-1,4))
              end do
           end if
#endif
           do l=1,ngrid
              div(i,j,k) = ux(l) + vy(l) + wz(l)
           end do

        end do
     end do
  end do

end subroutine cmpdivu_gpu_2d
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine consup_gpu_2d(uin,flux,div,dt,ngrid,nxp)
  use amr_parameters
  use hydro_parameters
  use acc_parameters
  use const
  implicit none

  integer ::ngrid,nxp
  real(dp)::dt
  real(dp),dimension(-1:nxp+2,-1:nxp+2,-1:nxp+2,1:nvar)::uin 
  real(dp),dimension(1:nxp+1,1:nxp+1,1:nxp+1,1:nvar,1:ndim)::flux
  real(dp),dimension(1:nxp+1,1:nxp+1,1:nxp+1)::div 

  integer:: i, j, k, l, n
  real(dp)::factor
  real(dp),dimension(1:nvector),save:: div1

  factor=half**(ndim-1)

  ! Add diffusive flux where flow is compressing
  do n = 1, nvar

     do k = kf10, MAX(kf1,ku20-2)
        do j = jf10, MAX(jf1, ju20-2) 
           do i = if10, if20

              div1 = zero
              do l = 1, ngrid
                 div1(l) = factor*div(i,j,k)
              end do
#if NDIM>1
              if(ndim>1)then
                 do l = 1, ngrid
                    div1(l)=div1(l)+factor*div(i,j+1,k)
                 end do
              end if
#endif
#if NDIM>2
              if(ndim>2)then
                 do l = 1, ngrid
                    div1(l)=div1(l)+factor*(div(i,j,k+1)+div(i,j+1,k+1))
                    div1(l)=difmag*min(zero,div1(l))
                 end do
              end if
#endif
              do l = 1, ngrid
                 div1(l) = difmag*min(zero,div1(l))
              end do
              do l = 1, ngrid
                 flux(i,j,k,n,1) = flux(i,j,k,n,1) + &
                      &  dt*div1(l)*(uin(i,j,k,n) - uin(i-1,j,k,n))
              end do

           end do
        end do
     end do

#if NDIM>1
     if(ndim>1)then
     do k = kf10, MAX(kf1,ku20-2)
        do j = jf10, jf20
           do i = iu10+2, iu20-2

              div1 = zero
              do l = 1, ngrid
                 div1(l)=div1(l)+factor*(div(i,j,k ) + div(i+1,j,k))
              end do
              if(ndim>2)then
                 do l = 1, ngrid
                    div1(l)=div1(l)+factor*(div(i,j,k+1) + div(i+1,j,k+1))
                 end do
              end if
              do l = 1, ngrid
                 div1(l) = difmag*min(zero,div1(l))
              end do
              do l = 1, ngrid
                 flux(i,j,k,n,2) = flux(i,j,k,n,2) + &
                      &  dt*div1(l)*(uin(i,j,k,n) - uin(i,j-1,k,n))
              end do

           end do
        end do
     end do
     end if
#endif

#if NDIM>2
     if(ndim>2)then
     do k = kf10, kf20
        do j = ju10+2, ju20-2 
           do i = iu10+2, iu20-2 

              do l = 1, ngrid
                 div1(l)=factor*(div(i,j  ,k) + div(i+1,j  ,k) &
                      &        + div(i,j+1,k) + div(i+1,j+1,k))
              end do
              do l = 1, ngrid
                 div1(l) = difmag*min(zero,div1(l))
              end do
              do l = 1, ngrid
                 flux(i,j,k,n,3) = flux(i,j,k,n,3) + &
                      &  dt*div1(l)*(uin(i,j,k,n) - uin(i,j,k-1,n))
              end do
           end do
        end do
     end do
     end if
#endif

  end do

end subroutine consup_gpu_2d
