
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmpdt(uu,gg,dx,dt,ncell)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer::ncell
  real(dp)::dx,dt
  real(dp),dimension(1:nvector,1:nvar)::uu ! conservative
  real(dp),dimension(1:nvector,1:ndim)::gg
  real(dp),dimension(1:nvector,1:nvar)::q

  real(dp) ::dtcell,smallp
  integer  ::k,idim
  real(dp) ::velx,vely,velz
  real(dp) :: lor,entho ! Lorentz factor
  real(dp) :: D,M,E,Mx,My,Mz,u2,Xsi,R
  real(dp) ::rho,p,vpar,vx,vy,vz
  integer::a,b,c

!  smallp = smallc**2/gamma
  dt=courant_factor*dx/smallc
  !convert to primitive variables

  call ctoprimbis(uu,ncell,q)

  do  k=1,ncell
  !compute fastest signal speed (x dir)
     call find_speed_info((/q(k,1),q(k,5),q(k,2),q(k,3),q(k,4)/),velx)

#if NDIM > 1
     !compute fastest signal speed (y dir)
     call find_speed_info((/q(k,1),q(k,5),q(k,3),q(k,2),q(k,4)/),vely)
#endif

#if NDIM == 3
  !compute fastest signal speed (z dir)
     call find_speed_info((/q(k,1),q(k,5),q(k,4),q(k,2),q(k,3)/),velz)
#endif


#if NDIM == 1
     dt=min(dt,dx/velx)
#endif

#if NDIM == 2
     dt=min(dt,dx*courant_factor/(velx+vely))
#endif

#if NDIM == 3
     dt=min(dt,dx/(velx+vely+velz) )
#endif


  end do


end subroutine cmpdt
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine hydro_refine(ug,um,ud,ok,nn)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none
  ! dummy arguments
  integer nn
  real(dp),dimension(1:nvector,1:nvar)::ug
  real(dp),dimension(1:nvector,1:nvar)::um
  real(dp),dimension(1:nvector,1:nvar)::ud
  real(dp),dimension(1:nvector,1:nvar)::qg
  real(dp),dimension(1:nvector,1:nvar)::qm
  real(dp),dimension(1:nvector,1:nvar)::qd
  real(dp),dimension(1:nvar)          ::qvarg,qvarm,qvard
  logical ::ok(1:nvector)



  integer::k,idim
  real(dp),dimension(1:nvector),save::eking,ekinm,ekind
  real(dp)::dg,dm,dd,pg,pm,pd,vg,vm,vd,cg,cm,cd,error,lorg,lorm,lord


     !convert to primitive variables
  call ctoprimbis(ug,nn,qg)
  call ctoprimbis(um,nn,qm)
  call ctoprimbis(ud,nn,qd)

  ug=qg
  um=qm
  ud=qd

  ! Compute errors
  if(err_grad_d >= 0.)then
     do k=1,nn
        dg=ug(k,1); dm=um(k,1); dd=ud(k,1)
        error=2.0d0*MAX( &
             & ABS((dd-dm)/(dd+dm+floor_d)) , &
             & ABS((dm-dg)/(dm+dg+floor_d)) )
        ok(k) = ok(k) .or. error > err_grad_d
     end do
  end if

  if(err_grad_p >= 0.)then
     do k=1,nn
        pg=ug(k,5); pm=um(k,5); pd=ud(k,5)
        error=2.0d0*MAX( &
             & ABS((pd-pm)/(pd+pm+floor_p)), &
             & ABS((pm-pg)/(pm+pg+floor_p)) )
        ok(k) = ok(k) .or. error > err_grad_p
     end do
  end if


  if(err_grad_u >= 0.)then
        do k=1,nn

           ! compute signal velocity
           qvarg(1)=ug(k,1); qvarm(1)=um(k,1) ; qvard(1)=ud(k,1)
           qvarg(2)=ug(k,5); qvarm(2)=um(k,5) ; qvard(2)=ud(k,5)

           ! x direction
           qvarg(3)=ug(k,2) ; qvarm(3)=um(k,2)  ; qvard(3)=ud(k,2)
           qvarg(4)=ug(k,3) ; qvarm(4)=um(k,3)  ; qvard(3)=ud(k,3)
           qvarg(5)=ug(k,4) ; qvarm(5)=um(k,4)  ; qvard(3)=ud(k,4)
           call find_speed_info(qvarg,vg)
           call find_speed_info(qvarm,vm)
           call find_speed_info(qvard,vd)
           error=2.0d0*MAX( &
                & ABS((ud(k,2)-um(k,2))/(ABS(vd)+ABS(vm)+floor_u)) , &
                & ABS((ug(k,2)-um(k,2))/(ABS(vg)+ABS(vm)+floor_u)) )
           ok(k) = ok(k) .or. error > err_grad_u



           ! y direction
           qvarg(3)=ug(k,3) ; qvarm(3)=um(k,3)  ; qvard(3)=ud(k,3)
           qvarg(4)=ug(k,2) ; qvarm(4)=um(k,2)  ; qvard(3)=ud(k,2)
           qvarg(5)=ug(k,4) ; qvarm(5)=um(k,4)  ; qvard(3)=ud(k,4)
           call find_speed_info(qvarg,vg)
           call find_speed_info(qvarm,vm)
           call find_speed_info(qvard,vd)

           error=2.0d0*MAX( &
                & ABS((ud(k,3)-um(k,3))/(ABS(vd)+ABS(vm)+floor_u)) , &
                & ABS((ug(k,3)-um(k,3))/(ABS(vg)+ABS(vm)+floor_u)) )
           ok(k) = ok(k) .or. error > err_grad_u



           ! z direction
           qvarg(3)=ug(k,4) ; qvarm(3)=um(k,4)  ; qvard(3)=ud(k,4)
           qvarm(4)=ug(k,3) ; qvarm(4)=um(k,3)  ; qvard(3)=ud(k,3)
           qvard(5)=ug(k,2) ; qvarm(5)=um(k,2)  ; qvard(3)=ud(k,2)
           call find_speed_info(qvarg,vg)
           call find_speed_info(qvarm,vm)
           call find_speed_info(qvard,vd)

           error=2.0d0*MAX( &
                & ABS((ud(k,4)-um(k,4))/(ABS(vd)+ABS(vm)+floor_u)) , &
                & ABS((ug(k,4)-um(k,4))/(ABS(vg)+ABS(vm)+floor_u)) )
           ok(k) = ok(k) .or. error > err_grad_u


        end do

  end if


  if(err_grad_lor >= 0.)then

        do k=1,nn
           lorg=(1.-ug(k,2)**2-ug(k,3)**2-ug(k,4)**4)**(-1./2.)
           lorm=(1.-um(k,2)**2-um(k,3)**2-um(k,4)**4)**(-1./2.)
           lord=(1.-ud(k,2)**2-ud(k,3)**2-ud(k,4)**4)**(-1./2.)
           error=2.0d0*MAX( &
                & ABS((lord-lorm)/(lord +lorm )) , &
                & ABS((lorg-lorm)/(lorg+lorm  )))
           ok(k) = ok(k) .or. error > err_grad_lor

        end do

  end if



end subroutine hydro_refine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine riemann_hllc(qleft,qright,fgdnv)

  use amr_parameters
  use const
  use hydro_parameters

  implicit none

  real(dp), dimension(1:nvar) :: qleft,qright,fgdnv
  real(dp), dimension(1:nvar) :: fleft,fright
  real(dp), dimension(1:nvar) :: uleft,uright
  real(dp),dimension(1:nvar)  :: usl,usr,fhll,uhll
  real(dp)                    :: vleft,vright,cfleft,cfright,sl,sr
  real(dp)                    :: lleftp,lleftm,lrightp,lrightm,factorl,factorr
  real(dp)                    :: vp1l,vp1r,vp2l,vp2r,vnl,vnr,lorl,lorr,vplsq,vprsq,vlsq,vrsq
  real(dp)                    :: den,ps,sstar,ovs,a,b,c,quad,srsl
  integer                     :: i
  real(dp)                    :: h,lor




  call find_mhd_flux(qleft ,uleft ,fleft )
  call find_mhd_flux(qright,uright,fright)

!--- Step 1. ------------------------------------------------------------------
! Compute the max and min wave speeds (used in Mignone )
!
  call find_speed_fast(qleft ,cfleft )
  call find_speed_fast(qright,cfright)


  vnl  =qleft (3)
  vnr  =qright(3)
  vp1l =qleft (4)
  vp1r =qright(4)
  vp2l =qleft (5)
  vp2r =qright(5)

  vlsq   =vnl**2+vp1l**2+vp2l**2
  vrsq   =vnr**2+vp1r**2+vp2r**2
  vplsq  =vp1l**2+vp2l**2
  vprsq  =vp1l**2+vp2l**2

  ! lorentz factors
  lorl=(1.0d0-vlsq)**(-1./2.)
  lorr=(1.0d0-vrsq)**(-1./2.)


!From Del Zanna et al, 2002.Case with transverse velocity component
  !wave speeds

  factorl = cfleft *(sqrt(1.0d0-vnl**2-cfleft**2 *vplsq))/lorl
  factorr = cfright*(sqrt(1.0d0-vnr**2-cfright**2*vprsq))/lorr

  lleftp =(vnl*(1.0d0- cfleft**2) +factorl)/(1.0d0- cfleft**2*vlsq)
  lleftm =(vnl*(1.0d0- cfleft**2) -factorl)/(1.0d0- cfleft**2*vlsq)

  lrightp=(vnr*(1.0d0- cfright**2)+factorr)/(1.0d0- cfright**2*vrsq)
  lrightm=(vnr*(1.0d0- cfright**2)-factorr)/(1.0d0- cfright**2*vrsq)

  sr=max(lleftp,lrightp)
  sl=min( lleftm,lrightm)

!--- Step 2. ------------------------------------------------------------------
! Compute L/R fluxes according to Mignone 2 (done above)

!--- Step 3. ------------------------------------------------------------------
! Compute HLL flux using Mignone Eq 11 (necessary for computing lmdas (Eq 18)
! Compute HLL conserved quantities using Mignone eq 9

  ovs = 1.0 / ( sr - sl )
  srsl = sr*sl
    do i=1, nvar
     fhll(i)=(sr * fleft(i) -  sl * fright(i) + srsl * (uright(i) - uleft(i))) *ovs
     uhll(i)=(sr * uright(i) - sl * uleft(i)  + fleft(i) - fright(i) ) * ovs
  enddo

!--- Step 4. ------------------------------------------------------------------
! Compute contact wave speed using larger root from Mignone Eq 18
! Physical root is the root with the minus sign

 ! quadratic formUla calcUlation

  a = fhll(2)
  b = -(uhll(2) +fhll(3))
  c = uhll(3)

  if (b >0) then
     quad = -0.5*(b +sqrt(b*b - 4.0*a*c))
  else
     quad = -0.5*(b -sqrt(b*b - 4.0*a*c))
  endif
  sstar = c/quad

!--- Step 5. ------------------------------------------------------------------
! Determine intercell flux according to Mignone 13


  if( sl >= 0.0) then ! Fl
     ! intercell flux is left flux
     fgdnv = fleft
  else if  (sstar >= 0.0)  then ! Fls
     ! Mignone 2006 Eq 48
     ps = -Fhll(2)*sstar + Fhll(3)

   ! now calcUlate Usl with Mignone Eq 16
    den = 1.0 / (sl - sstar);
    usl(1) =  uleft(1) * (sl - qleft(3))                                 * den
    usl(3) = (uleft(3) * (sl - qleft(3) )+ps - qleft(2))                 * den
    usl(4) =  uleft(4) * (sl - qleft(3))                                 * den
    usl(5) =  uleft(5) * (sl - qleft(3))                                 * den
    usl(2) = (uleft(2) * (sl - qleft(3)) +ps *sstar - qleft(2)*qleft(3)) * den

    ! now calcUlate Fsr using Mignone Eq 14
    fgdnv= sl *(usl -uleft) +fleft

  else if( sr >= 0.0) then ! Frs

    ! Mignone 2006 Eq 48
    ps = -fhll(2)*sstar + Fhll(3)

    ! now calcUlate Usr with Mignone Eq 16
    den = 1.0 / (sr - sstar)

    usr(1) =  uright(1) * (sr - qright(3))* den
    usr(3) = (uright(3) * (sr - qright(3))+ps - qright(2))* den
    usr(4) =  uright(4) * (sr - qright(3))* den
    usr(5) =  uright(5) * (sr - qright(3))* den
    usr(2) = (uright(2) * (sr - qright(3)) +ps *sstar - qright(2)*qright(3))* den

    !now calcUlate Fsr using Mignone Eq 14

    fgdnv= sr *(usr -uright) +fright

  else
   ! intercell flux is right flux
    fgdnv =fright
 endif
 !passive scalars
 do i=6,nvar  ! il faut un facteur de lorentz????
    if(fgdnv(1)>0)then
       fgdnv(i) = fgdnv(1)*qleft (i)
    else
       fgdnv(i) = fgdnv(1)*qright(i)
    endif
 enddo


  return
end subroutine riemann_hllc

!###########################################################
!###########################################################
!###########################################################
!###########################################################
!  Subroutine HLL
!
!> HLL Riemann solver.
!< q=(d,p,v_n,b_n,v_p1,b_p1,v_p2,b_p2)

subroutine riemann_hll(qleft,qright,fgdnv)

  use amr_parameters
  use const
  use hydro_parameters

  implicit none

  real(dp), dimension(1:nvar) :: qleft,qright,fgdnv
  real(dp), dimension(1:nvar) :: fleft,fright
  real(dp), dimension(1:nvar) :: uleft,uright,udiff,ugdnv
  real(dp)                    :: vleft,vright,cfleft,cfright,sl,sr
  real(dp)                    :: vlsq,vrsq,vplsq,vprsq,factorl,factorr
  real(dp)                    :: lleftp,lleftm,lrightp,lrightm,vp1l,vp1r,vp2l,vp2r,vnl,vnr,lorl,lorr


  call find_mhd_flux(qleft,uleft ,fleft )
  call find_mhd_flux(qright,uright,fright)

!sound speed
  call find_speed_fast(qleft ,cfleft )
  call find_speed_fast(qright,cfright)

  vnl  =qleft (3)
  vnr  =qright(3)
  vp1l =qleft (4)
  vp1r =qright(4)
  vp2l =qleft (5)
  vp2r =qright(5)

  vlsq   =vnl**2+vp1l**2+vp2l**2
  vrsq   =vnr**2+vp1r**2+vp2r**2
  vplsq  =vp1l**2+vp2l**2
  vprsq  =vp1l**2+vp2l**2


  ! lorentz factors
  lorl=(1.0d0-vlsq)**(-1./2.)
  lorr=(1.0d0-vrsq)**(-1./2.)

  !From Del Zanna et al, 2002.Case with transverse velocity component
  !wave speeds

  factorl = cfleft *(sqrt(1.0d0-vnl**2-cfleft**2 *vplsq))/lorl
  factorr = cfright*(sqrt(1.0d0-vnl**2-cfright**2*vprsq))/lorr

  lleftp =(vnl*(1.0d0- cfleft**2) +factorl)/(1.0d0- cfleft**2*vlsq)
  lleftm =(vnl*(1.0d0- cfleft**2) -factorl)/(1.0d0- cfleft**2*vlsq)

  lrightp=(vnr*(1.0d0- cfright**2)+factorr)/(1.0d0- cfright**2*vrsq)
  lrightm=(vnr*(1.0d0- cfright**2)-factorr)/(1.0d0- cfright**2*vrsq)


  sr=max(zero, lleftp,lrightp)
  sl=max(zero, -lleftm,-lrightm)

  !Compute HLL fluxes

  fgdnv = (sr*fleft+sl*fright-sr*sl*(uright-uleft))/(sr+sl)
  return
end subroutine riemann_hll
!###########################################################
!###########################################################
!###########################################################
!###########################################################

SUBROUTINE find_speed_info(qvar,vel_info)
  USE amr_parameters
  USE const
  USE hydro_parameters
  !! calculate the fastest velocity at which information is exchanged
  !! at the interface
  !! the structure of qvar is : rho, Pressure, Vnormal, Vpar1,Vpar2
  IMPLICIT NONE
  integer                 :: i,ib,iv
  real(dp), dimension(1:nvar) :: qvar
  real(dp)                    :: vel_info
  real(dp)                    :: d,p,u,v,w,velsq,vperpsq,factor
  real(dp)                    :: ein,enth,cs,lor,tau

  d=qvar(1); p=qvar(2); u=qvar(3) ; v=qvar(4); w=qvar(5)
  velsq=u*u+v*v+w*w ; vperpsq=v*v+w*w
  lor=sqrt(one/(one-velsq))
  !modification for special relativity
  if (eos .eq.'TM') then
     tau=p/d
     enth=5d0/2d0*tau+3d0/2d0*sqrt(tau**2+4d0/9d0)
     cs=sqrt(tau/3d0/enth*(5d0*enth-8d0*tau)/(enth-tau))
  else
     ein = d+p/(gamma-1) ; enth= ein + p
     cs  = sqrt(gamma*p/enth)
  endif
  factor = cs*(sqrt(1.0d0-u*u-cs*cs*vperpsq))/lor
  vel_info=max(abs(u*(1.0d0-cs*cs)+factor),abs(u*(1.0d0-cs*cs)-factor))/(1.0d0-cs*cs*velsq)

  return


END SUBROUTINE find_speed_info


!###########################################################
!###########################################################
!###########################################################
!###########################################################
!  Subroutine FIND_MHD_FLUX
!
!> Compute the 1d mhd fluxes. The structure of qvar is : rho, pressure,
!! vnormal, vtransverse1,  vtransverse2 (+passive scalars)
! Follows equation 6 from Del Zanna, Bucciantini, 2002,A&A,390,1177

subroutine find_mhd_flux(qvar,cvar,ff)

  use amr_parameters
  use const
  use hydro_parameters

  implicit none

  integer                 :: i,ivar,ib
  real(dp), dimension(1:nvar) :: qvar,cvar,ff
  real(dp)                    :: etot,d,u,v,w,p,entho,ein,enth,tau
  real(dp)                :: lor ! Lorentz factor

  ! local variables

  d=qvar(1); p=qvar(2); u=qvar(3);  v=qvar(4); w=qvar(5)
  lor  = (1-(u**2+v**2+w**2))**(-1./2.)
  entho = one/(gamma-one)

  if (eos .eq. 'TM') then
     tau=p/d
     enth=d*(5d0/2d0*tau+3d0/2d0*sqrt(tau**2+4d0/9d0))
  else
     ein  = p*entho+d
     enth = ein + p !relativistic enthalpy
  endif
  ! compute conservative variables
  cvar(1) = d*lor
  cvar(2) = enth*lor**2-p
  cvar(3) = enth*lor**2*u
  cvar(4) = enth*lor**2*v
  cvar(5) = enth*lor**2*w
  do ivar = 6,nvar
     cvar(ivar) = d*qvar(ivar)*lor
  end do

  ! compute fluxes
  ff(1) = d*u*lor
  ff(2) = enth*lor**2*u
  ff(3) = enth*lor**2*u**2+p
  ff(4) = enth*lor**2*u*v
  ff(5) = enth*lor**2*u*w
  do ivar =6,nvar
     ff(ivar) = d*u*qvar(ivar)*lor
  end do
  return

end subroutine find_mhd_flux



!###########################################################
!###########################################################
!###########################################################
!###########################################################
!  Subroutine FIND_SPEED_FAST
!
!> Calculate the sound speed
!! The structure of qvar is : rho, pressure, vnormal,vtransverse1, vtransverse2,

subroutine find_speed_fast(qvar,vel_info)

  use amr_parameters
  use const
  use hydro_parameters

  implicit none

  real(dp), dimension(1:nvar) :: qvar
  real(dp)                    :: vel_info,cs
  real(dp)                    :: d,p,u,v,w,enth,lor,ein,tau

  d=qvar(1); p=qvar(2); u=qvar(3); v=qvar(4); w=qvar(5)
  if (eos .eq. 'TM') then
     tau=p/d
     enth=5d0/2d0*tau+3d0/2d0*sqrt(tau**2+4d0/9d0)
     cs=sqrt(tau/3d0/enth*(5d0*enth-8d0*tau)/(enth-tau))
  else
     ein = d+p/(gamma-1)
     enth= ein + p
     cs  = sqrt(gamma*p/(enth))
  endif
  !In hydro case, just returns the sound speed
  vel_info = cs
  return

end subroutine find_speed_fast


subroutine ctoprimbis(uu,n,q)
use amr_parameters
use const
use hydro_parameters

implicit none

  real(dp),dimension(1:nvector,1:nvar)::uu
  real(dp),dimension(1:nvector,1:nvar)::q

  integer::n,k,idim
  real(dp) :: lor,entho ! Lorentz factor
  real(dp) :: D,M,E,Mx,My,Mz,u2,Xsi,R,rho
  real(dp) :: small_bigD=1e-12
  real(dp)::D2,M2,E2,tau,h

  do k=1,n
     !convert to primitive variables

     ! Compute density
     D = uu(k,1)
     ! Compute momentum
     Mx=uu(k,2) ; My=uu(k,3) ; Mz=uu(k,4)
     M = sqrt(Mx**2+My**2+Mz**2)
     ! Compute total energy
     E = uu(k,5)
     !Method from Mignone,McKinney,2007. Same as BS2011 except one uses E'=U-D and u^2=Lor^2*v^2

     if (D<0) then
       write(*,*),'D<0 in ctoprimbis'
        D=small_bigD
        uu(k,1)=D
     endif
     if (E<0) then
        E=sqrt(M**2+d**2+1d-8)
       write(*,*),'E<0 in ctoprimbis'
        uu(k,5)=E
     endif

     if (E**2<M**2+D**2) then

       write (*,*) 'Switch...ctoprimbis'
        E=sqrt(M**2+d**2+1d-8)
        uu(k,5)=E
     endif

     if ( M .eq. 0) then
        q(k,1) = D
        q(k,2) = 0d0
        q(k,3) = 0d0
        q(k,4) = 0d0
        if (eos .eq. 'TM') then
           q(k,5) =(E**2-D**2)/3d0/E
        else
           q(k,5)=(E-D)*(gamma-1d0)
        endif
        lor=1d0

     else

        call Newton_Raphson_Mignone(D,M,E,gamma,R)
        ! Compute the Lorentz factor
        u2  = M**2.0d0/(R**2.0d0-M**2.0d0)
        lor = (1.0d0+u2)**(1d0/2d0)

        ! Compute the density
        q(k,1) = D/lor

        ! compute velocities
        q(k,2) = Mx/R
        q(k,3) = My/R
        q(k,4) = Mz/R

        ! Compute pressure
        Xsi=((R-D)-u2/(lor+1d0)*D)/lor**2
        if (eos .eq. 'TM') then
           rho=q(k,1)
           q(k,5)=(2d0*xsi*(xsi+2d0*rho))/(5d0*(xsi+rho)+sqrt(9d0*xsi**2+18d0*rho*xsi+25d0*rho**2))
        else
           q(k,5)=(gamma-1d0)/gamma*Xsi
        endif
     endif
     if ((q(k,1)<0d0).or.(q(k,5)<0d0).or.E<0d0) then
        write(*,*) 'negative pressure or density ctoprim 2'
        stop
     endif
     ! Passive scalar
     do idim = 6, nvar
        q(k,idim) = uu(k,idim)/q(k,1)/lor
     end do
     tau=q(k,5)/q(k,1)
     h=5d0/2d0*tau+3d0/2d0*sqrt(tau**2+4d0/9d0)
     enddo



return


end subroutine ctoprimbis
