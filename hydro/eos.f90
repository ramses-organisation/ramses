subroutine barotropic_eos_temperature(nH, temperature)
   use amr_parameters
   implicit none
   !--------------------------------------------------------------
   ! This routine selects the chosen EOS and calculates the
   ! temperature T[in Kelvin]/mu from the density nH[in H/cc]
   !--------------------------------------------------------------
   real(dp), intent(in) ::nH
   real(dp), intent(out)::temperature
   real(dp)::factor1,factor2,factor3

   SELECT CASE (barotropic_eos_form)
   CASE ('isothermal')
      temperature = T2_eos
   CASE ('polytrope')
      temperature = T2_eos*(nH/polytrope_n(1))**(polytrope_index(1)-1.0d0)
   CASE ('double_polytrope')
      temperature = T2_eos * (1 + (nH/polytrope_n(1))**(polytrope_index(1)-1.0d0))
   CASE ('2nd_collapse')
      ! Machida & Inutsuka 2006, Marchand et al. 2016
      factor1 = sqrt(1 + (nH/polytrope_n(1))**(2*polytrope_index(1)))
      factor2 = (1 + (nH/polytrope_n(2)))**polytrope_index(2)
      factor3 = (1 + (nH/polytrope_n(3)))**polytrope_index(3)
      temperature = T2_eos * factor1 * factor2 * factor3
   CASE ('custom')
      ! WRITE YOUR FAVORITE EOS HERE
      if(nH<polytrope_n(1))then
         temperature = T2_eos
      else
         temperature = T2_eos * (nH/polytrope_n(1))**(polytrope_index(1)-1.0d0)
      endif
   CASE DEFAULT
     write(*,*)'unknown barotropic eos form'
     call clean_stop
   END SELECT

end subroutine barotropic_eos_temperature
!################################################################
!################################################################
!################################################################
!################################################################
subroutine temperature_eos(rho_temp,Enint_temp,Teos,ht)
  use amr_parameters      ,only:dp,mu_gas
  use hydro_commons       !,only:gamma
  use constants           ,only:kB,mH
  use cooling_module, only:barotropic_eos
  implicit none
  !--------------------------------------------------------------
  ! This routine computes the temperature from the density and 
  ! internal volumic energy. Inputs/output are in code units.
  !--------------------------------------------------------------
  real(dp), intent(in) :: Enint_temp,rho_temp
  integer , intent(out):: ht 
  real(dp), intent(out):: Teos
  real(dp)::rho,Enint
  real(dp)::scale_nH,scale_T2,scale_t,scale_v,scale_d,scale_l

  integer::i_t,i_r,i
  real(dp)::logr,tt,uu,y1,y2,y3,y4
  real(dp):: le,lr
  real(dp):: dd1,dd2,de1,de2
  integer :: ir,ie
  real(dp):: xx,drho,dener

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

if(eos)then
  if (enint_temp ==0.d0) then
     teos=0.d0
  else
     ht=0

     rho   = rho_temp * scale_d             
     Enint = Enint_temp * scale_d*scale_v**2

     drho  = (rhomax-rhomin)/float(nRho)
     lr = 0.5d0 + (log10(rho )- rhomin )/drho

     if (lr .ge. Nrho) then
        write(*,*)'pb 1'
        stop
     endif
     ir = floor(lr)

     dEner = ( Emax  -   emin)/float(nEnergy)
     le = 0.5d0 + (log10(Enint) - emin - log10(rho) )/dEner

     if ((le .ge. nEnergy) ) then
        write(*,*)'pb 2'
        stop
     endif

     ie = floor(le)
     if  (ir < 1 .or. ie < 1 ) then
        write(*,*) 'inter_tp hors limite ir,ie,rho,enint = ',ir,ie,rho ,Enint 
        ir=1.0d0
        ie=1.0d0
        stop
     endif

     dd1 = lr - float(ir)
     dd2 = le - float(ie)

     de1 = 1.0d0 - dd1
     de2 = 1.0d0 - dd2

     Teos=0.d0

     Teos = Teos + de1*de2*Temp_eos(ir  ,ie  )
     Teos = Teos + dd1*de2*Temp_eos(ir+1,ie  )
     Teos = Teos + de1*dd2*Temp_eos(ir  ,ie+1)
     Teos = Teos + dd1*dd2*Temp_eos(ir+1,ie+1)

     Teos = Teos ! give T in K

     xx =  Temp_eos(ir,ie)*Temp_eos(ir+1,ie)*Temp_eos(ir,ie+1)*Temp_eos(ir+1,ie+1) 

     if (xx .eq. 0.0d0 ) then
        ht=1
        !     write(*,*) '**************** Pb_eos ****************'
        !     write(*,*) 'ir,ie,i =',ir,ie,i
        !     write(*,*) rho ,Enint 
        !     write(*,*) Temp_eos(ir,ie  ), Temp_eos(ir+1,ie  ) 
        !     write(*,*) Temp_eos(ir,ie+1), Temp_eos(ir+1,ie+1)
        !     write(*,*) Temp_eos(ir,ie+2), Temp_eos(ir-1,ie+1)
        !     write(*,*) 0.25*(Temp_eos(ir,ie+2) + Temp_eos(ir-1,ie+1) +   Temp_eos(ir+1,ie+1) +  Temp_eos(ir,ie))
        !     stop
     endif
  endif
else 
  rho   = rho_temp*scale_d
  Enint = Enint_temp*scale_d*scale_v**2 
  if (barotropic_eos) then
      call barotropic_eos_temperature(rho_temp, Teos)
      Teos = Teos*mu_gas
  else
      Teos = Enint/(rho*kB/(mu_gas*mH*(gamma-1.0d0)))
  endif
  ht=1
endif

  return

end subroutine temperature_eos
!################################################################
!################################################################
!################################################################
!################################################################
subroutine enerint_eos(rho_temp,temp_temp,Eeos)
  use amr_parameters      ,only:dp,mu_gas
  use hydro_commons       !,only:gamma
  use constants           ,only:kB,mH
  implicit none
  !--------------------------------------------------------------
  ! This routine computes the internal volumic energy from  
  ! the density and the temperature. Inputs/output are in code units.
  !--------------------------------------------------------------
  real(dp), intent(in) :: temp_temp,rho_temp
  real(dp), intent(out):: Eeos
  real(dp)::rho,temp
  real(dp)::scale_nH,scale_T2,scale_t,scale_v,scale_d,scale_l

  integer::i_t,i_r,i
  real(dp)::logr,tt,uu,y1,y2,y3,y4
  real(dp):: le,lr
  real(dp):: dd1,dd2,de1,de2
  integer :: ir,ie
  real(dp):: xx,drho,dtemp

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

 if(eos)then
     rho = rho_temp * scale_d
     temp   = temp_temp
     
     drho  = (rhomax-rhomin)/float(nRho)
     lr = 0.5d0 + (log10(rho )- rhomin )/drho

  if (lr .ge. Nrho) then
     write(*,*)'pb 11'
     stop
  endif
  ir = floor(lr)

  dtemp = ( log10(Tmax)  -   log10(Tmin))/float(ntemp)
  le = 1.0d0 + (log10(temp) - log10(Tmin))/dtemp

  if ((le .ge. ntemp) ) then
     write(*,*)'pb 22'
     stop
  endif

  ie = floor(le)
  if  (ir < 1 .or. ie < 1 ) then
     write(*,*) 'inter_ener hors limite ir,ie,rho,enint cooling= ',ir,ie,rho ,temp
     ir=1.0d0
     ie=1.0d0
     stop
  endif


  dd1 = lr - float(ir)
  dd2 = le - float(ie)

  de1 = 1.0d0 - dd1
  de2 = 1.0d0 - dd2

  Eeos=0.d0

  Eeos = Eeos + de1*de2*eint_eos(ir  ,ie  )
  Eeos = Eeos + dd1*de2*eint_eos(ir+1,ie  )
  Eeos = Eeos + de1*dd2*eint_eos(ir  ,ie+1)
  Eeos = Eeos + dd1*dd2*eint_eos(ir+1,ie+1)

  Eeos = Eeos / (scale_d*scale_v**2) ! give energy in code units

  xx =  eint_eos(ir,ie)*eint_eos(ir+1,ie)*eint_eos(ir,ie+1)*eint_eos(ir+1,ie+1) 

  if (xx .eq. 0.0d0 ) then
     !     write(*,*) '**************** Pb_eos ****************'
     !     write(*,*) 'ir,ie,i =',ir,ie,i
     !     write(*,*) rho ,Enint 
     !     write(*,*) eint_eos(ir,ie  ), eint_eos(ir+1,ie  ) 
     !     write(*,*) eint_eos(ir,ie+1), eint_eos(ir+1,ie+1)
     !     write(*,*) eint_eos(ir,ie+2), eint_eos(ir-1,ie+1)
     !     write(*,*) 0.25*(eint_eos(ir,ie+2) + eint_eos(ir-1,ie+1) +   eint_eos(ir+1,ie+1) +  eint_eos(ir,ie))
     !     stop
  endif

else
  rho  = rho_temp * scale_d
  temp = temp_temp

  Eeos = rho*kB/(mu_gas*mH*(gamma-1.0))*temp/(scale_d*scale_v**2)
endif

  return

end subroutine enerint_eos
!################################################################
!################################################################
!################################################################
!################################################################
subroutine soundspeed_eos(rho_temp,Enint_temp,Cseos)
  use amr_parameters      ,only:dp
  use hydro_commons       !,only:gamma
  implicit none
  !--------------------------------------------------------------
  ! This routine computes the sound speed from the internal volumic energy 
  ! and the temperature. Inputs/output are in code units.
  !--------------------------------------------------------------
  real(dp), intent(in) :: Enint_temp,rho_temp
  real(dp), intent(out):: Cseos

  integer::i_t,i_r,i
  real(dp):: Enint,rho
  real(dp)::logr,tt,uu,y1,y2,y3,y4
  real(dp):: le,lr
  real(dp):: dd1,dd2,de1,de2
  integer :: ir,ie
  real(dp):: xx,drho,dener

  real(dp)::scale_nH,scale_T2,scale_t,scale_v,scale_d,scale_l

  if(eos)then
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
   rho = rho_temp * scale_d
  Enint   = Enint_temp * scale_d*scale_v**2

  drho  = (rhomax-rhomin)/float(nRho)
  lr = 0.50d0 + (log10(rho )- rhomin )/drho

  if (lr .ge. Nrho) then
     write(*,*)'intermaxrho',rho ,rho_eos(nRho,1)
     stop
  endif
  ir = floor(lr)

  dEner = ( Emax  -   emin)/float(nEnergy)
  le = 0.50d0 + (log10(Enint) - emin - log10(rho) )/dEner

  if ((le .ge. nEnergy) ) then
     write(*,*)'intermaxE_s',Enint ,Ener_eos(ir,nEnergy)
     stop
  endif

  ie = floor(le)
  if  (ir < 1 .or. ie < 1 ) then
     write(*,*) 'inter_Cs hors limite i,ir,ie = ',ir,ie,rho ,Enint
     ir=1.0d0
     ie=1.0d0
     stop
  endif

  dd1 = lr - float(ir)
  dd2 = le - float(ie)

  de1 = 1.0d0 - dd1
  de2 = 1.0d0 - dd2

  Cseos=0.d0


  Cseos = Cseos + de1*de2*Cs_eos(ir  ,ie  )
  Cseos = Cseos + dd1*de2*Cs_eos(ir+1,ie  )
  Cseos = Cseos + de1*dd2*Cs_eos(ir  ,ie+1)
  Cseos = Cseos + dd1*dd2*Cs_eos(ir+1,ie+1)

  Cseos = Cseos /scale_v ! give Cs in code units

  xx =  Cs_eos(ir,ie)*Cs_eos(ir+1,ie)*Cs_eos(ir,ie+1)*Cs_eos(ir+1,ie+1) 

  if (xx .eq. 0.0d0 ) then
     write(*,*) '**************** Cs_eos ****************'
     write(*,*) 'ir,ie,i =',ir,ie,i
     write(*,*) rho ,Enint 
     write(*,*) Cs_eos(ir,ie  ), Cs_eos(ir+1,ie  ) 
     write(*,*) Cs_eos(ir,ie+1), Cs_eos(ir+1,ie+1)
     write(*,*) Cs_eos(ir,ie+2), Cs_eos(ir-1,ie+1)
     write(*,*) 0.25*(Cs_eos(ir,ie+2) + Cs_eos(ir-1,ie+1) +   Cs_eos(ir+1,ie+1) +  Cs_eos(ir,ie))
     stop
  endif
else
  Cseos = sqrt(gamma*(gamma-1.d0)*Enint_temp/rho_temp)
endif

  return

end subroutine soundspeed_eos
!################################################################
!################################################################
!################################################################
!################################################################
function cmp_Cv_eos(rho,Enint)
  use amr_parameters      ,only:dp,mu_gas
  use hydro_commons       !,only:gamma
  use constants           ,only:kB,mH
  implicit none
  !--------------------------------------------------------------
  ! This function computes the Cv from the density and 
  ! internal volumic energy. Inputs/output are in code units.
  !--------------------------------------------------------------
  real(dp)   :: rho,Enint
  real(dp)   :: cmp_Cv_eos

  real(dp)::scale_nH,scale_T2,scale_t,scale_v,scale_d,scale_l

  real(dp)   :: drho,dener,xx
  real(dp)   :: dd1,dd2,de1,de2
  real(dp)   :: lr,le
  integer    :: ir,ie,i

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  
 
if(eos)then
     drho  = (rhomax-rhomin)/float(nRho)
  lr = 0.5d0 + (log10(rho )- rhomin )/drho

  if (lr .ge. Nrho) then
     write(*,*)'intermaxrho',rho ,rho_eos(nRho,1)
     stop
  endif
  ir = floor(lr)

  dEner = ( Emax  -   emin)/float(nEnergy)
  le = 0.5d0 + (log10(Enint) - emin - log10(rho) )/dEner

  if ((le .ge. nEnergy) ) then
     write(*,*)'intermaxE2',Enint ,Ener_eos(ir,nEnergy)
     stop
  endif

  ie = floor(le)
  if  (ir < 1 .or. ie < 1 ) then
!!$      write(*,*) 'inter_cv hors limite i,ir,ie = ',ir,ie,rho ,Enint
     ir=1.0d0
     ie=1.0d0
     !      stop
  endif

  dd1 = lr - float(ir)
  dd2 = le - float(ie)

  de1 = 1.0d0 - dd1
  de2 = 1.0d0 - dd2

  cmp_Cv_eos=0.d0

  cmp_Cv_eos = cmp_Cv_eos + de1*de2*Cv_eos(ir  ,ie  )
  cmp_Cv_eos = cmp_Cv_eos + dd1*de2*Cv_eos(ir+1,ie  )
  cmp_Cv_eos = cmp_Cv_eos + de1*dd2*Cv_eos(ir  ,ie+1)
  cmp_Cv_eos = cmp_Cv_eos + dd1*dd2*Cv_eos(ir+1,ie+1)

  cmp_Cv_eos = cmp_Cv_eos

  xx =  Cv_eos(ir,ie)*Cv_eos(ir+1,ie)*Cv_eos(ir,ie+1)*Cv_eos(ir+1,ie+1)

  if (xx .eq. 0.0d0 ) then
     write(*,*) '**************** Pb_eos ****************'
     write(*,*) 'ir,ie,i =',ir,ie,i
     write(*,*) rho ,Enint
     write(*,*) Cv_eos(ir,ie  ), Cv_eos(ir+1,ie  )
     write(*,*) Cv_eos(ir,ie+1), Cv_eos(ir+1,ie+1)
     write(*,*) Cv_eos(ir,ie+2), Cv_eos(ir-1,ie+1)
     write(*,*) 0.25*(Cv_eos(ir,ie+2) + Cv_eos(ir-1,ie+1) +   Cv_eos(ir+1,ie+1) +  Cv_eos(ir,ie))
     stop
  endif

else
  cmp_Cv_eos = rho*kB/(mu_gas*mH*(gamma-1.0d0))/scale_v**2
endif

end function cmp_Cv_eos
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_eos
  use amr_commons
  use hydro_commons
  use constants

   implicit none

   integer::ht !j,irad,ht
   real(dp)::scale_nH,scale_T2,scale_t,scale_v,scale_d,scale_l

  integer::ii,jj,kk,ee,hh,gg,ie,ir,k,it
  real(dp)::xx,yy,vv,ww,zz
  real(dp)::dtemp1,Temp_new2,epsilon_n,eint_old,T0,temp_new,d_loc,eint_new

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

     !--------------------------------
     ! Read eos tables
     !--------------------------------
!      open(14,file='verif.dat')
     open(10,file='tab_eos.dat',status='old',form='unformatted')
     read(10) nRho,nEnergy
     read(10) rhomin,rhomax,emin,Emax,yHe
     
     allocate(Rho_eos(nRho,nEnergy),Ener_eos(nRho,nEnergy),Temp_eos(nRho,nEnergy),P_eos(nRho,nEnergy))
     allocate(  Cs_eos(nRho,nEnergy),S_eos(nRho,nEnergy),  xH_eos(nRho,nEnergy), xH2_eos(nRho,nEnergy)                  )
     allocate(xHe_eos(nRho,nEnergy),xHep_eos(nRho,nEnergy),Cv_eos(nRho,nEnergy)                                       )
     !inversion de la table eos
     nTemp=nEnergy
     allocate(eint_eos(nRho,nTemp))
     
     read(10)  rho_eos
     read(10) Ener_eos
     read(10) Temp_eos
     read(10)    P_eos
     read(10)    S_eos
     read(10)   Cs_eos
     read(10)   xH_eos
     read(10)  xH2_eos
     read(10)  xHe_eos
     read(10) xHep_eos
     close(10)
     
     rho_eos(:,:) = log10(rho_eos(:,:))
     ener_eos(:,:) = log10(ener_eos(:,:))
     
     do k=1,5
        ii=0
        jj=0
        kk=0
        hh=0
        ee=0
        gg=0
        do ir=2,nRho-1
           do ie=2,nEnergy-1
              if (P_eos(ir,ie) .eq. 0.0d0) then
                 ii = ii+1
                 xx = P_eos(ir,ie+1) * P_eos(ir,ie-1) *  P_eos(ir-1,ie) * P_eos(ir+1,ie)
                 yy = P_eos(ir+1,ie+1) * P_eos(ir+1,ie-1) *  P_eos(ir-1,ie-1) * P_eos(ir-1,ie+1)
                 if(ie > 2 .and. ie < nEnergy-1 .and. ir > 2 .and. ir < nRho-1)then
                    ww = P_eos(ir,ie+2) * P_eos(ir,ie-2) *  P_eos(ir-2,ie) * P_eos(ir+2,ie)
                 else
                    ww = 0.0_dp
                 endif
                 if(ie > 3 .and. ie < nEnergy-2 .and. ir > 3 .and. ir < nRho-2)then
                    zz = P_eos(ir+3,ie+3) * P_eos(ir-3,ie-3) *  P_eos(ir-3,ie+3) * P_eos(ir+3,ie-3)
                 else
                    zz = 0.0_dp
                 endif
                 if (xx .ne. 0.) then
                    P_eos(ir,ie) = 0.25d0*(P_eos(ir,ie+1) + P_eos(ir,ie-1) + P_eos(ir-1,ie) + P_eos(ir+1,ie))
                    jj=jj+1              
                 else if (yy .ne. 0. .and. k > 0) then
                    P_eos(ir,ie) = 0.25d0*(P_eos(ir+1,ie+1) + P_eos(ir+1,ie-1) + P_eos(ir-1,ie+1)+P_eos(ir-1,ie-1))
                    kk=kk+1
                 else if (ww .ne. 0 .and. k > 1) then
                    ee = ee +1
                    P_eos(ir,ie) = 0.25d0*(P_eos(ir,ie+2) + P_eos(ir,ie-2) + P_eos(ir-2,ie) + P_eos(ir+2,ie))
                 else if (zz .ne. 0 .and. k > 2) then
                    hh=hh+1
                    P_eos(ir,ie) = 0.25d0*(P_eos(ir+3,ie+3) + P_eos(ir+3,ie-3) + P_eos(ir-3,ie+3)+P_eos(ir-3,ie-3))
                 else 
                    gg=gg+1
                 endif
              endif
           enddo
        end do
        if (myid == 1) print*, "on bouche les trous P_eos", ii,jj,kk,ee,hh,gg, "iter", k
     end do
     
     do k=1,5
        ii=0
        jj=0
        kk=0
        hh=0
        ee=0
        gg=0
        do ir=2,nRho-1
           do ie=2,nEnergy-1
              if (Cs_eos(ir,ie) .eq. 0.0d0) then           
                 ii = ii+1
                 xx = Cs_eos(ir,ie+1) * Cs_eos(ir,ie-1) *  Cs_eos(ir-1,ie) * Cs_eos(ir+1,ie)
                 yy = Cs_eos(ir+1,ie+1) * Cs_eos(ir+1,ie-1) *  Cs_eos(ir-1,ie-1) * Cs_eos(ir-1,ie+1)
                 if(ie > 2 .and. ie < nEnergy-1 .and. ir > 2 .and. ir < nRho-1)then
                    ww = Cs_eos(ir,ie+2) * Cs_eos(ir,ie-2) *  Cs_eos(ir-2,ie) * Cs_eos(ir+2,ie)
                 else
                    ww = 0.0_dp
                 endif
                 if(ie > 3 .and. ie < nEnergy-2 .and. ir > 3 .and. ir < nRho-2)then
                    zz = Cs_eos(ir+3,ie+3) * Cs_eos(ir-3,ie-3) *  Cs_eos(ir-3,ie+3) * Cs_eos(ir+3,ie-3)
                 else
                    zz = 0.0_dp
                 endif
                 if (xx .ne. 0.) then
                    Cs_eos(ir,ie) = 0.25d0*(Cs_eos(ir,ie+1) + Cs_eos(ir,ie-1) + Cs_eos(ir-1,ie) + Cs_eos(ir+1,ie))
                    jj=jj+1              
                 else if (yy .ne. 0. .and. k > 0) then
                    Cs_eos(ir,ie) = 0.25d0*(Cs_eos(ir+1,ie+1) + Cs_eos(ir+1,ie-1) + Cs_eos(ir-1,ie+1)+Cs_eos(ir-1,ie-1))
                    kk=kk+1
                 else if (ww .ne. 0 .and. k > 1) then
                    ee = ee +1
                    Cs_eos(ir,ie) = 0.25d0*(Cs_eos(ir,ie+2) + Cs_eos(ir,ie-2) + Cs_eos(ir-2,ie) + Cs_eos(ir+2,ie))
                 else if (zz .ne. 0 .and. k > 2) then
                    hh=hh+1
                    Cs_eos(ir,ie) = 0.25d0*(Cs_eos(ir+3,ie+3) + Cs_eos(ir+3,ie-3) + Cs_eos(ir-3,ie+3)+Cs_eos(ir-3,ie-3))
                 else 
                    gg=gg+1
                 endif
              endif
           enddo
        end do
        if (myid == 1) print*, "on bouche les trous Cs_eos", ii,jj,kk,ee,hh,gg, "iter", k
     end do
     
     do k=1,5
        ii=0
        jj=0
        kk=0
        hh=0
        ee=0
        gg=0
        do ir=2,nRho-1
           do ie=2,nEnergy-1
              if (Temp_eos(ir,ie) .eq. 0.0d0) then           
                 ii = ii+1
                 xx = Temp_eos(ir,ie+1) * Temp_eos(ir,ie-1) *  Temp_eos(ir-1,ie) * Temp_eos(ir+1,ie)
                 yy = Temp_eos(ir+1,ie+1) * Temp_eos(ir+1,ie-1) *  Temp_eos(ir-1,ie-1) * Temp_eos(ir-1,ie+1)
                 if(ie > 2 .and. ie < nEnergy-1 .and. ir > 2 .and. ir < nRho-1)then
                    ww = Temp_eos(ir,ie+2) * Temp_eos(ir,ie-2) *  Temp_eos(ir-2,ie) * Temp_eos(ir+2,ie)
                 else
                    ww = 0.0_dp
                 endif
                 if(ie > 3 .and. ie < nEnergy-2 .and. ir > 3 .and. ir < nRho-2)then
                    zz = Temp_eos(ir+3,ie+3) * Temp_eos(ir-3,ie-3) *  Temp_eos(ir-3,ie+3) * Temp_eos(ir+3,ie-3)
                 else
                    zz = 0.0_dp
                 endif
                 if (xx .ne. 0.) then
                    Temp_eos(ir,ie) = 0.25d0*(Temp_eos(ir,ie+1)+Temp_eos(ir,ie-1)+Temp_eos(ir-1,ie)+Temp_eos(ir+1,ie))
                    jj=jj+1              
                 else if (yy .ne. 0. .and. k > 0) then
                    Temp_eos(ir,ie) = 0.25d0*(Temp_eos(ir+1,ie+1)+Temp_eos(ir+1,ie-1)+Temp_eos(ir-1,ie+1)+Temp_eos(ir-1,ie-1))
                    kk=kk+1
                 else if (ww .ne. 0 .and. k > 1) then
                    ee = ee +1
                    Temp_eos(ir,ie) = 0.25d0*(Temp_eos(ir,ie+2)+Temp_eos(ir,ie-2)+Temp_eos(ir-2,ie)+Temp_eos(ir+2,ie))
                 else if (zz .ne. 0 .and. k > 2) then
                    hh=hh+1
                    Temp_eos(ir,ie) = 0.25d0*(Temp_eos(ir+3,ie+3)+Temp_eos(ir+3,ie-3)+Temp_eos(ir-3,ie+3)+Temp_eos(ir-3,ie-3))
                 else 
                    gg=gg+1
                 endif
              endif
           enddo
        end do
        
        if (myid == 1) print*, "on bouche les trous Temp_eos", ii,jj,kk,ee,hh,gg, "iter", k
     end do
     
     Tmin=3.0d0
     Tmax=1.0d5
     dtemp1 =(log10(Tmax) - log10(Tmin))/ntemp
     eint_eos(:,:)=0.0d0
     do ir=2,nRho-1
        do it=1,ntemp
           d_loc = (10.**rho_eos(ir,1))
           T0 = 10.**(log10(Tmin) + (it-1.0d0)*dtemp1)
           
           eint_old = d_loc*kb*T0/(mu_gas*mh*(gamma-1.0d0))
           if (it >1) then
              eint_old = max(d_loc*kb*T0/(mu_gas*mh*(gamma-1.0d0)),eint_eos(ir,it-1))
           end if
           
           epsilon_n = 1.0d0
           
           do ii=1,1000
              call temperature_eos(d_loc/scale_d,eint_old/(scale_d*scale_v**2),temp_new,ht)
              if (ht==1) then
                 eint_old=0.d0
                 exit
              end if
              call temperature_eos(d_loc/scale_d,eint_old*1.001_dp/(scale_d*scale_v**2),temp_new2,ht)
              if (ht==1) then
                 eint_old=0.d0
                 exit
              end if
              
              if(abs(Temp_new2-Temp_new) .ne.0)then
                 eint_new = eint_old - (Temp_new-T0)/((Temp_new2-Temp_new)/(0.001*eint_old))
              else
                 eint_new = eint_old
              endif
              epsilon_n = abs(eint_new - eint_old)/eint_old
              eint_old = eint_new
              if  (abs(epsilon_n) .lt. 1.d-4) then
                 exit
              else if (ii==1000) then
                 print*, "newton for e(rho,T) did not converge at ", log10(d_loc),log10(T0)
              end if
           end do
           eint_eos(ir,it) = eint_old 
        enddo
     enddo
     
     Cv_eos(:,:)=0.0d0
     
     do  ir=2,nRho-1
        do  ie=2,nEnergy-1
           d_loc = (10.**rho_eos(ir,1))
           T0 = 10.**(ener_eos(ir,ie))
           
           call temperature_eos(d_loc/scale_d,(T0-0.001_dp*T0)/(scale_d*scale_v**2),temp_new,ht)
           call temperature_eos(d_loc/scale_d,(T0+0.001_dp*T0)/(scale_d*scale_v**2),temp_new2,ht)
           
           if((temp_new2-temp_new) .ne. 0.0_dp)then
              Cv_eos(ir,ie)=(0.002_dp*T0)/(temp_new2-temp_new)
           else
              Cv_eos(ir,ie) = 1.0_dp
           endif
        end do
     end do
     Cv_eos(:,nEnergy)=Cv_eos(:,nEnergy-1)
     
     if (myid == 1) print*, "ok pour le bouchage"

   return

end subroutine init_eos