!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine tracex(q,dq,c, qm,qp,dx,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer::ngrid
  real(dp)::dx, dt

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::dq
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qm
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qp
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::c

  ! Local variables
  integer ::ilo,ihi,jlo,jhi,klo,khi
  integer ::i, j, k, l, n
  integer ::ir, iu, ip
  real(dp)::dtdx,project_out
  real(dp)::cc, ccc, csq, r, u, p, a,b,v,w,eta
  real(dp)::drx, dux, dpx, dax, dvx, dwx
  real(dp)::alpham, alphap, alpha0x,alpha0y,alpha0z
  real(dp)::spminus, spplus, spzero
  real(dp)::apright, amright, azrxright, azryright, azrzright
  real(dp)::apleft,  amleft,  azrxleft, azryleft, azrzleft
  real(dp)::ni,h,lor,velsq,smallp, rp , rm,entho,vtot2,coef,vperp2


  dtdx = dt/dx
  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)
  ir=1; iu=2; ip=5
  project_out=one !zero
  ni = 1./(gamma-1.)
  smallp = 1.0e-10

!!!!
  entho = gamma/(gamma-one)

  if (eos .eq. 'TM') then
     write(*,*),'plmde does not work with TM eos, switch to a MUSCL scheme of change the EOS'
     stop
  endif

  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           do l = 1, ngrid



              ! Cell centered values
              r   = q (l,i,j,k,1)
              u   = q (l,i,j,k,2)
              v   = q (l,i,j,k,3)
              w   = q (l,i,j,k,4)
              p   = q (l,i,j,k,5)
              h   = one+P/r*gamma/(gamma-one)
              csq = (gamma*p/(r*h) )
              cc  =sqrt(csq) ! relativistic sound speed


              eta =(one-u**2-csq*(v**2+w**2))**(half)
              b= h*csq
              vtot2= u**2+v**2+w**2
              lor=(one-vtot2)**(-half)
              a= cc*eta/(lor*r)

              ! TVD slopes in X direction
              drx = dq(l,i,j,k,1,1)
              dux = dq(l,i,j,k,2,1)
              dvx = dq(l,i,j,k,3,1)
              dwx = dq(l,i,j,k,4,1)
              dpx = dq(l,i,j,k,5,1)


              ! Supersonic fix for high-velocity gradients
              ccc = cc

              if(ABS(dux) > three*cc)ccc=zero
              if(ABS(dvx) > three*cc)ccc=zero
              if(ABS(dwx) > three*cc)ccc=zero

              ! Characteristic analysis along X direction
              alpham =                   (-half/a)*dux               +half/b*dpx
              alpha0x = drx                                          -one/b*dpx
              alpha0y =             u*v/(one-u**2)*dux +dvx         + v/(one-u**2)/lor**2/r/h*dpx
              alpha0z =             u*w/(one-u**2)*dux       + dwx  + w/(one-u**2)/lor**2/r/h*dpx
              alphap =                    (half/a)*dux              +half/b*dpx

              !eigenvalues
              coef= (one-vtot2)*(one-u**2-ccc**2*(v**2+w**2))
              spminus = (u*(one-ccc**2)-ccc*sqrt(coef))/(one-vtot2*ccc**2)*dtdx
              spplus  = (u*(one-ccc**2)+ccc*sqrt(coef))/(one-vtot2*ccc**2)*dtdx
              spzero  = (u    )*dtdx



              ! Right state


              if(((u*(one-ccc**2)+ccc*sqrt(coef))/(one-vtot2*ccc**2))>zero)spplus =-project_out
              if(((u*(one-ccc**2)-ccc*sqrt(coef))/(one-vtot2*ccc**2))>zero)spminus=-project_out
              if( u     >zero)spzero =-project_out


              apright  = half*(-one-spplus )*alphap
              amright  = half*(-one-spminus)*alpham
              azrxright = half*(-one-spzero )*alpha0x
              azryright = half*(-one-spzero )*alpha0y
              azrzright = half*(-one-spzero )*alpha0z


              rm=u*(one-csq)-cc*eta/lor
              rp=u*(one-csq)+cc*eta/lor

              rm=cc*rm/(r*(u*cc+lor*eta))
              rp=cc*rp/(r*(u*cc-lor*eta))


              qp(l,i,j,k,1,1) = r + (apright+ azrxright+amright)
              qp(l,i,j,k,2,1) = u + (-amright+apright)*a
              qp(l,i,j,k,3,1) = v + (amright*v*rm + azryright + apright*v*rp)
              qp(l,i,j,k,4,1) = w + (amright*w*rm + azrzright + apright*w*rp)
              qp(l,i,j,k,5,1) = p + ( amright+apright)*b

              qp(l,i,j,k,ir,1) = max(smallr, qp(l,i,j,k,ir,1))
              qp(l,i,j,k,ip,1) = max(smallp, qp(l,i,j,k,ip,1))


              ! Left state

              spminus = (u*(one-ccc**2)-ccc*sqrt(coef))/(one-vtot2*ccc**2)*dtdx
              spplus  = (u*(one-ccc**2)+ccc*sqrt(coef))/(one-vtot2*ccc**2)*dtdx
              spzero  = (u    )*dtdx


              if( (u*(one-ccc**2)+ccc*sqrt(coef))/(one-vtot2*ccc**2)   <=zero)spplus =+project_out
              if( (u*(one-ccc**2)-ccc*sqrt(coef))/(one-vtot2*ccc**2)   <=zero)spminus =+project_out
              if( u     <=zero)spzero =+project_out


              apleft   = half*(+one-spplus )*alphap
              amleft   = half*(+one-spminus)*alpham
              azrxleft = half*(+one-spzero )*alpha0x
              azryleft = half*(+one-spzero )*alpha0y
              azrzleft = half*(+one-spzero )*alpha0z


              qm(l,i,j,k,1,1) = r + (apleft+azrxleft+amleft)
              qm(l,i,j,k,2,1) = u + (-amleft+apleft)*a
              qm(l,i,j,k,3,1) = v + (amleft*v*rm + azryleft + apleft*v*rp)
              qm(l,i,j,k,4,1) = w + (w*rm*amleft + azrzleft + apleft*w*rp)
              qm(l,i,j,k,5,1) = p + ( amleft+apleft)*b

              qm(l,i,j,k,ir,1) = max(smallr, qm(l,i,j,k,ir,1))
              qm(l,i,j,k,ip,1) = max(smallp, qm(l,i,j,k,ip,1))

              velsq=qp(l,i,j,k,2,1)*qp(l,i,j,k,2,1)+qp(l,i,j,k,3,1)*qp(l,i,j,k,3,1)+qp(l,i,j,k,4,1)*qp(l,i,j,k,4,1)
              if (velsq>1d0) then
                 qp(l,i,j,k,ir:ip,1)=q(l,i,j,k,ir:ip)
                 qm(l,i,j,k,ir:ip,1)=q(l,i,j,k,ir:ip)
              endif


              velsq=qm(l,i,j,k,2,1)*qm(l,i,j,k,2,1)+qm(l,i,j,k,3,1)*qm(l,i,j,k,3,1)+qm(l,i,j,k,4,1)*qm(l,i,j,k,4,1)
              if (velsq>1d0) then
                 qp(l,i,j,k,:,1)=q(l,i,j,k,:)
                 qm(l,i,j,k,:,1)=q(l,i,j,k,:)
              endif

              ! Passive scalars
              ! simple advection
              do n = 6, nvar
                 a = q(l,i,j,k,n)     ! Cell centered values
                 u = q(l,i,j,k,2)
                 dax = dq(l,i,j,k,n,1)    ! TVD slopes

                 ! Right state
                 spzero=(u    )*dtdx
                 if(u>zero)spzero=-project_out
                 azrxright = half*(-one-spzero )*dax
                 qp(l,i,j,k,n,1) = a + azrxright

                 ! Left state
                 spzero=(u    )*dtdx
                 if(u<=zero)spzero=+project_out
                 azrxleft = half*(+one-spzero )*dax
                 qm(l,i,j,k,n,1) = a + azrxleft

              end do
           end do
        end do
     end do
  end do


end subroutine tracex
!###########################################################
!###########################################################
!###########################################################
!###########################################################
#if NDIM>1
subroutine tracexy(q,dq,c,qm,qp,dx,dy,dt,ngrid)
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
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::c

  ! declare local variables
  integer ::ilo,ihi,jlo,jhi,klo,khi
  integer ::i, j, k, l, n
  integer ::ir, iu, iv, IW, ip
  real(dp)::dtdx,dtdy,project_out
  real(dp)::cc, ccc, csq, r, u, v, p, w, h,  a
  real(dp)::alpham, alphap, alpha0x, alpha0y, alpha0z
  real(dp)::spminus, spplus, spzero
  real(dp)::dry, duy, dvy, dwy, dpy, day,say, DAX, SAX
  real(dp)::apright, amright, azrxright, azryright, azrzright, azaright
  real(dp)::apleft,  amleft,  azrxleft,  azryleft,  azrzleft,  azaleft
  real(dp)::eta,b,vtot,lor,vtot2,coef,rm,rp,smallp,vels
  real(dp):: D,Mx,My,E,M,Mz,velsq



  dtdx = dt/dx; dtdy = dt/dy
  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)
  ir=1; iu=2; iv=3; iw=4;  ip=5
  project_out=one !zero

  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           do l = 1, ngrid


              ! cell centered values
              r   = q (l,i,j,k,1)
              u   = q (l,i,j,k,2)
              v   = q (l,i,j,k,3)
              w   = q (l,i,j,k,4)
              p   = q (l,i,j,k,5)
              h   = one+P/r*gamma/(gamma-one)
              csq = (gamma*p/(r*h) )
              cc  =sqrt(csq) ! relativistic sound speed

              eta =(one-v**2-csq*(u**2+w**2))**(half)
              b= h*csq
              vtot2= u**2+v**2+w**2
              lor=(one-vtot2)**(-1./2.)
              a= cc*eta/(lor*r)


              !Slopes in Y direction
              dry = dq(l,i,j,k,ir,2)
              duy = dq(l,i,j,k,iu,2)
              dvy = dq(l,i,j,k,iv,2)
              dwy = dq(l,i,k,k,iw,2)
              dpy = dq(l,i,j,k,ip,2)



              ! Supersonic fix for high-velocity gradients
              ccc = cc

              if(ABS(duy) > three*cc)ccc=zero
              if(ABS(dvy) > three*cc)ccc=zero
              if(ABS(dwy) > three*cc)ccc=zero



              ! Characteristic analysis along Y direction
              alpham =                 (-half/a)*dvy                          +half/b*dpy
              alpha0x =      duy +u*v/(one-v**2)*dvy        + u/(one-v**2)/lor**2/r/h*dpy
              alpha0y = dry                                                   - one/b*dpy
              alpha0z =          v*w/(one-v**2)*dvy  + dwy  + w/(one-v**2)/lor**2/r/h*dpy
              alphap  =                (half/a)*dvy                           +half/b*dpy



              ! Top state
              !eigenvalues
              coef= (one-vtot2)*(one-v**2-ccc**2*(u**2+w**2))
              spminus = (v*(one-ccc**2)-ccc*sqrt(coef))/(one-vtot2*ccc**2)*dtdy
              spplus  = (v*(one-ccc**2)+ccc*sqrt(coef))/(one-vtot2*ccc**2)*dtdy
              spzero  = (v    )*dtdy


              if(((v*(one-ccc**2)+ccc*sqrt(coef))/(one-vtot2*ccc**2))>zero)spplus =-project_out
              if(((v*(one-ccc**2)-ccc*sqrt(coef))/(one-vtot2*ccc**2))>zero)spminus=-project_out
              if( v     >zero)spzero =-project_out


              apright   = half*(-one-spplus )*alphap
              amright   = half*(-one-spminus)*alpham
              azrxright = half*(-one-spzero )*alpha0x
              azryright = half*(-one-spzero )*alpha0y
              azrzright = half*(-one-spzero )*alpha0z

              rm=v*(one-csq)-cc*eta/lor
              rp=v*(one-csq)+cc*eta/lor

              rm=cc*rm/(r*(v*cc+lor*eta))
              rp=cc*rp/(r*(v*cc-lor*eta))


              qp(l,i,j,k,1,2) = r + (apright+ azryright+amright)
              qp(l,i,j,k,2,2) = u + (amright*u*rm+ azrxright + apright*u*rp)
              qp(l,i,j,k,3,2) = v + (-amright + apright)*a
              qp(l,i,j,k,4,2) = w + (amright*w*rm + azrzright + apright*w*rp)
              qp(l,i,j,k,5,2) = p + ( amright+apright)*b

              qp(l,i,j,k,ir,2) = max(smallr, qp(l,i,j,k,ir,2))
              qp(l,i,j,k,ip,2) = max(smallp, qp(l,i,j,k,ip,2))


              ! Bottom state

              !eigenvalues

              spminus = (v*(one-ccc**2)-ccc*sqrt(coef))/(one-vtot2*ccc**2)*dtdy
              spplus  = (v*(one-ccc**2)+ccc*sqrt(coef))/(one-vtot2*ccc**2)*dtdy
              spzero  = (v    )*dtdy


              if(((v*(one-ccc**2)+ccc*sqrt(coef))/(one-vtot2*ccc**2))<=zero)spplus =+project_out
              if(((v*(one-ccc**2)-ccc*sqrt(coef))/(one-vtot2*ccc**2))<=zero)spminus=+project_out
              if( v     <zero)spzero =+project_out

              apleft   = half*(+one-spplus )*alphap
              amleft   = half*(+one-spminus)*alpham
              azrxleft = half*(+one-spzero )*alpha0x
              azryleft = half*(+one-spzero )*alpha0y
              azrzleft = half*(+one-spzero )*alpha0z

              qm(l,i,j,k,1,2) = r + (apleft+ azryleft+amleft)
              qm(l,i,j,k,2,2) = u + (amleft*u*rm+ azrxleft + apleft*u*rp)
              qm(l,i,j,k,3,2) = v + (-amleft + apleft)*a
              qm(l,i,j,k,4,2) = w + (amleft*w*rm + azrzleft + apleft*w*rp)
              qm(l,i,j,k,5,2) = p + ( amleft+apleft)*b

              qm(l,i,j,k,ir,2) = max(smallr, qm(l,i,j,k,ir,2))
              qm(l,i,j,k,ip,2) = max(smallp, qm(l,i,j,k,ip,2))



              velsq=qp(l,i,j,k,2,2)*qp(l,i,j,k,2,2)+qp(l,i,j,k,3,2)*qp(l,i,j,k,3,2)+qp(l,i,j,k,4,2)*qp(l,i,j,k,4,2)
              if (velsq>1d0) then
                 qp(l,i,j,k,ir:ip,2)=q(l,i,j,k,ir:ip)
                 qm(l,i,j,k,ir:ip,2)=q(l,i,j,k,ir:ip)
              endif


              velsq=qm(l,i,j,k,2,2)*qm(l,i,j,k,2,2)+qm(l,i,j,k,3,2)*qm(l,i,j,k,3,2)+qm(l,i,j,k,4,2)*qm(l,i,j,k,4,2)
              if (velsq>1d0) then
                 qp(l,i,j,k,:,2)=q(l,i,j,k,:)
                 qm(l,i,j,k,:,2)=q(l,i,j,k,:)
              endif

           end do
        end do
     end do
  end do

  ! Passive scalars
  do n = ndim+3, nvar
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              do l = 1, ngrid
                 a = q(l,i,j,k,n)     ! Cell centered values
                 u = q(l,i,j,k,iu)
                 v = q(l,i,j,k,iv)
                 dax = dq(l,i,j,k,n,1)    ! TVD slopes
                 day = dq(l,i,j,k,n,2)
                 sax = half*dtdy*(-v*day) ! Transverse
                 say = half*dtdx*(-u*dax) ! derivatives

                 ! Right state
                 spzero=(u    )*dtdx
                 if(u>zero)spzero=-project_out
                 azaright = half*(-one-spzero )*dax
                 qp(l,i,j,k,n,1) = a + azaright + sax

                 ! Left state
                 spzero=(u    )*dtdx
                 if(u<=zero)spzero=+project_out
                 azaleft = half*(+one-spzero )*dax
                 qm(l,i,j,k,n,1) = a + azaleft + sax

                 ! Top state
                 spzero=(v    )*dtdy
                 if(v>zero)spzero=-project_out
                 azaright = half*(-one-spzero )*day
                 qp(l,i,j,k,n,2) = a + azaright + say

                 ! Bottom state
                 spzero=(v    )*dtdy
                 if(v<=zero)spzero=+project_out
                 azaleft = half*(+one-spzero )*day
                 qm(l,i,j,k,n,2) = a + azaleft + say

              end do
           end do
        end do
     end do
  end do

end subroutine tracexy
#endif
!###########################################################
!###########################################################
!###########################################################
!###########################################################
#if NDIM>2
subroutine tracexyz(q,dq,c,qm,qp,dx,dy,dz,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer ::ngrid
  real(dp)::dx,dy,dz, dt

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::dq
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qm
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qp
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::c

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
  real(dp)::eta,b,vtot,lor,vtot2,coef,rm,rp,smallp,velsq,h,alpha0x,alpha0y,alpha0z
  real(dp)::azrzright,azrxright,azryright

  dtdx = dt/dx; dtdy = dt/dy; dtdz = dt/dz
  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)
  ir=1; iu=2; iv=3; iw=4; ip=5
  project_out=one !zero

  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           do l = 1, ngrid

              ! Cell centered values
              cc  = c  (l,i,j,k)
              r   = q  (l,i,j,k,ir)
              u   = q  (l,i,j,k,iu)
              v   = q  (l,i,j,k,iv)
              w   = q  (l,i,j,k,iw)
              p   = q  (l,i,j,k,ip)

              h   = one+P/r*gamma/(gamma-one)
              csq = (gamma*p/(r*h) )
              cc  =sqrt(csq) ! relativistic sound speed

              eta =(one-v**2-csq*(u**2+w**2))**(half)
              b= h*csq
              vtot2= u**2+v**2+w**2
              lor=(one-vtot2)**(-1./2.)
              a= cc*eta/(lor*r)

              ! TVD slopes in all Z directions

              drz = dq(l,i,j,k,ir,3)
              duz = dq(l,i,j,k,iu,3)
              dvz = dq(l,i,j,k,iv,3)
              dwz = dq(l,i,j,k,iw,3)
              dpz = dq(l,i,j,k,ip,3)


              ! Supersonic fix for high-velocity gradients
              ! has to be checked
              ccc = cc

              if(ABS(duz) > three*cc)ccc=zero
              if(ABS(dvz) > three*cc)ccc=zero
              if(ABS(dwz) > three*cc)ccc=zero



              ! Characteristic analysis along Z direction
              alpham  =                   (-1./2./a)*dwz     +1./2./b*dpz
              alpha0x =      duz      + u*w/(1-w**2)*dwz   + u/(1-w**2)/lor**2/r/h*dpz
              alpha0y =          +dvz + v*w/(1-w**2)*dwz   + v/(1-w**2)/lor**2/r/h*dpz
              alpha0z = drz                                   -1./b*dpz
              alphap  =                    (1./2./a)*dwz      +1./2./b*dpz


              ! Front state
              !eigenvalues
              coef    = (one-vtot2)*(one-w**2-ccc**2*(u**2+v**2))
              spminus = (w*(one-ccc**2)-ccc*sqrt(coef))/(one-vtot2*ccc**2)*dtdz
              spplus  = (w*(one-ccc**2)+ccc*sqrt(coef))/(one-vtot2*ccc**2)*dtdz
              spzero  = (w    )*dtdz

              if(((w*(one-ccc**2)+ccc*sqrt(coef))/(one-vtot2*ccc**2))>zero)spplus =-project_out
              if(((w*(one-ccc**2)-ccc*sqrt(coef))/(one-vtot2*ccc**2))>zero)spminus=-project_out
              if( w     >zero)spzero =-project_out


              apright  = half*(-one-spplus )*alphap
              amright  = half*(-one-spminus)*alpham
              azrright = half*(-one-spzero )*alpha0r
              azvright = half*(-one-spzero )*alpha0v
              azwright = half*(-one-spzero )*alpha0w

              rm=w*(one-csq)-cc*eta/lor
              rp=w*(one-csq)+cc*eta/lor

              rm=cc*rm/(r*(w*cc+lor*eta))
              rp=cc*rp/(r*(w*cc-lor*eta))



              qp(l,i,j,k,1,3) = r + (apright+ azrzright+amright)
              qp(l,i,j,k,2,3) = u + (amright*u*rm+ azrxright + apright*u*rp)
              qp(l,i,j,k,3,3) = v + (amright*v*rm + azryright + apright*w*rp)
              qp(l,i,j,k,4,3) = w + (-amright + apright)*a
              qp(l,i,j,k,5,3) = p + ( amright+apright)*b

              qp(l,i,j,k,ir,3) = max(smallr, qp(l,i,j,k,ir,3))
              qp(l,i,j,k,ip,3) = max(smallp, qp(l,i,j,k,ip,3))


              ! Back state
              coef    = (one-vtot2)*(one-w**2-ccc**2*(u**2+v**2))
              spminus = (w*(one-ccc**2)-ccc*sqrt(coef))/(one-vtot2*ccc**2)*dtdz
              spplus  = (w*(one-ccc**2)+ccc*sqrt(coef))/(one-vtot2*ccc**2)*dtdz
              spzero  = (w    )*dtdz

              if(((w*(one-ccc**2)+ccc*sqrt(coef))/(one-vtot2*ccc**2))<=zero)spplus =+project_out
              if(((w*(one-ccc**2)-ccc*sqrt(coef))/(one-vtot2*ccc**2))<=zero)spminus=+project_out
              if( w     <zero)spzero =+project_out


              apleft   = half*(+one-spplus )*alphap
              amleft   = half*(+one-spminus)*alpham
              azrleft  = half*(+one-spzero )*alpha0r
              azvleft  = half*(+one-spzero )*alpha0v
              azwleft  = half*(+one-spzero )*alpha0w


              qp(l,i,j,k,1,3) = r - (apright+ azrzright+amright)
              qp(l,i,j,k,2,3) = u - (amright*u*rm+ azrxright + apright*u*rp)
              qp(l,i,j,k,3,3) = v - (amright*v*rm + azryright + apright*w*rp)
              qp(l,i,j,k,4,3) = w - (-amright + apright)*a
              qp(l,i,j,k,5,3) = p - ( amright+apright)*b

              qp(l,i,j,k,ir,3) = max(smallr, qp(l,i,j,k,ir,3))
              qp(l,i,j,k,ip,3) = max(smallp, qp(l,i,j,k,ip,3))

           end do
        end do
     end do
  end do

  ! Passive scalars
  do n = ndim+3, nvar
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              do l = 1, ngrid
                 a   =  q(l,i,j,k,n)    ! Cell centered values
                 u   =  q(l,i,j,k,iu)
                 v   =  q(l,i,j,k,iv)
                 w   =  q(l,i,j,k,iw)
                 dax = dq(l,i,j,k,n,1)  ! TVD slopes
                 day = dq(l,i,j,k,n,2)
                 daz = dq(l,i,j,k,n,3)
                 sax = half*dtdx*(-v*day-w*daz) ! Transverse
                 say = half*dtdx*(-u*dax-w*daz) ! derivatives
                 saz = half*dtdx*(-v*day-u*dax) !


                 ! Right state
                 spzero = (u    )*dtdx
                 if(u>zero)spzero=-project_out
                 azaright = half*(-one-spzero )*dax
                 qp(l,i,j,k,n,1) = a + azaright + sax

                 ! Left state
                 spzero = (u    )*dtdx
                 if(u<=zero)spzero=+project_out
                 azaleft = half*(+one-spzero )*dax
                 qm(l,i,j,k,n,1) = a + azaleft + sax

                 ! Top state
                 spzero = (v    )*dtdy
                 if(v>zero)spzero=-project_out
                 azaright = half*(-one-spzero )*day
                 qp(l,i,j,k,n,2) = a + azaright + say

                 ! Bottom state
                 spzero = (v    )*dtdy
                 if(v<=zero)spzero=+project_out
                 azaleft = half*(+one-spzero )*day
                 qm(l,i,j,k,n,2) = a + azaleft + say

                 ! Front state
                 spzero = (w    )*dtdy
                 if(w>zero)spzero=-project_out
                 azaright = half*(-one-spzero )*daz
                 qp(l,i,j,k,n,3) = a + azaright + saz

                 ! Back state
                 spzero = (w    )*dtdy
                 if(w<=zero)spzero=+project_out
                 azaleft = half*(+one-spzero )*daz
                 qm(l,i,j,k,n,3) = a + azaleft + saz
              end do
           end do
        end do
     end do
  end do

end subroutine tracexyz
#endif
