MODULE coolrates_module
  use amr_parameters,only:dp
  use rt_parameters,only:nIons,isH2,isHe,ixHI,ixHII,ixHeII,ixHeIII
  implicit none
  
CONTAINS

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION comp_Alpha_H2(T,Z)

! Returns creation rate of h2 on dust [cm^3 s-1] (Draine and Bertoldi 1996)
! plus gas phase rate for low Z on H- assuming equilibrium abundances for H-
! as explained in the Appendix of McKee and Krumolz (2012)
! T           => Temperature [K]
! Z           => Metallicity in Solar units
!-------------------------------------------------------------------------
  implicit none  
  real(dp),intent(in)::T,Z
  real(dp)::comp_Alpha_H2, lambda
!-------------------------------------------------------------------------
  comp_Alpha_H2 =  Z * 6.0d-18*(T**0.5)/ &
       & (1.0+0.4*(T/100.)**0.5+0.2*(T/100.)+0.08*(T/100.)**2) !&
!       & + 8.0d-19*((T/1000.0)**0.88)
END FUNCTION comp_Alpha_H2

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION comp_AlphaA_HII(T)

! Returns case A rec. coefficient [cm3 s-1] for HII (Hui&Gnedin'97)
! T           => Temperature [K]
!-------------------------------------------------------------------------
  implicit none  
  real(dp),intent(in)::T
  real(dp)::comp_AlphaA_HII, lambda
!-------------------------------------------------------------------------
  lambda=  315614./T                                ! 2.d0 * 157807.d0 / T
  comp_AlphaA_HII =  1.269d-13 * lambda**1.503 &
                     / ( ( 1.d0+(lambda/0.522)**0.47 )**1.923 )
END FUNCTION comp_AlphaA_HII

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION comp_dAlphaA_dT_HII(T)

! Returns Temperature derivative of the case A rec. coeff. of HII 
! T           => Temperature [K]
!-------------------------------------------------------------------------
  implicit none  
  real(dp),intent(in)::T
  real(dp)::comp_dAlphaA_dT_HII, lambda, f
!-------------------------------------------------------------------------
  lambda=  315614./T
  f= 1.d0+(lambda/0.522)**0.47
  comp_dAlphaA_dT_HII = 1.269d-13 * lambda**1.503 / f**1.923             &
                        / T * ( 0.90381*(f-1.)/f - 1.503 )
END FUNCTION comp_dAlphaA_dT_HII

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION comp_AlphaA_HeII(T)

! Returns case A rec. coefficient [cm3 s-1] for HeII (Hui&Gnedin'97)
! T           => Temperature [K]
!-------------------------------------------------------------------------
  implicit none  
  real(dp),intent(in)::T
  real(dp)::comp_AlphaA_HeII, lambda
!-------------------------------------------------------------------------
  lambda=  570670./T
  comp_AlphaA_HeII =  3.d-14 * lambda**0.654
END FUNCTION comp_AlphaA_HeII

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION comp_AlphaA_HeIII(T)

! Returns case A rec. coefficient [cm3 s-1] for HeIII (Hui&Gnedin'97)
! T           => Temperature [K]
!-------------------------------------------------------------------------
  implicit none  
  real(dp),intent(in)::T
  real(dp)::comp_AlphaA_HeIII, lambda
!-------------------------------------------------------------------------
  lambda=  1263030./T
  comp_AlphaA_HeIII =  2.538d-13 * lambda**1.503                         &
                     / ( ( 1.d0+(lambda/0.522)**0.47 )**1.923 )
END FUNCTION comp_AlphaA_HeIII

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION comp_AlphaB_HII(T)

! Returns case B rec. coefficient [cm3 s-1] for HII (Hui&Gnedin'97)
! T           => Temperature [K]
!-------------------------------------------------------------------------
  implicit none  
  real(dp),intent(in)::T
  real(dp)::comp_AlphaB_HII, lambda
!-------------------------------------------------------------------------
  lambda = 315614./T
  comp_AlphaB_HII =  &
       2.753d-14 * lambda**1.5 / ( (1.d0+(lambda/2.74)**0.407)**2.242 )
END FUNCTION comp_AlphaB_HII

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION comp_AlphaB_HeII(T)

! Returns case B rec. coefficient [cm3 s-1] for HeII (Hui&Gnedin'97)
! T           => Temperature [K]
!-------------------------------------------------------------------------
  implicit none  
  real(dp),intent(in)::T
  real(dp)::comp_AlphaB_HeII, lambda
!-------------------------------------------------------------------------
  lambda = 570670./T
  comp_AlphaB_HeII = 1.26d-14 * lambda**0.75
END FUNCTION comp_AlphaB_HeII

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION comp_AlphaB_HeIII(T)

! Returns case B rec. coefficient [cm3 s-1] for HeIII (Hui&Gnedin'97)
! T           => Temperature [K]
!-------------------------------------------------------------------------
  implicit none  
  real(dp),intent(in)::T
  real(dp)::comp_AlphaB_HeIII, lambda
!-------------------------------------------------------------------------
  lambda=  1263030./T
  comp_AlphaB_HeIII =  5.506d-14 * lambda**1.5                           &
                     / ( ( 1.d0+(lambda/2.74)**0.407 )**2.242 )
END FUNCTION comp_AlphaB_HeIII

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION comp_dAlphaB_dT_HII(T)

! Returns temperature derivative of the case B recombination rate
! Gnedin 1997
! T           => Temperature [K]
!-------------------------------------------------------------------------
  implicit none  
  real(dp),intent(in)::T
  real(dp)::comp_dAlphaB_dT_HII, lambda, f
!-------------------------------------------------------------------------
  lambda = 315614./T
  f= 1.d0+(lambda/2.74)**0.407
  comp_dAlphaB_dT_HII = 2.753d-14 * lambda**1.5 / f**2.242               &
                  / T * ( 0.912494*(f-1.)/f - 1.5 )
END FUNCTION comp_dAlphaB_dT_HII

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION comp_Beta_H2HI(T)

! Returns collisional rate [cm3 s-1] of H2 and HI (Abel 1997->Dove&Mandy 1986)
! T           => Temperature [K]
!-------------------------------------------------------------------------
  implicit none  
  real(dp),intent(in)::T
  real(dp)::comp_Beta_H2HI, kbT
!-------------------------------------------------------------------------
	kbT=8.618d-5*T !eV
  comp_Beta_H2HI = &
                3.324d-9*(kbT**2.012)*exp(-4.463d0/kbT)/(1.d0+0.2472d0*kbT)**3.512
END FUNCTION comp_Beta_H2HI

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION comp_Beta_H2H2(T)

! Returns collisional rate [cm3 s-1] of H2 and H2 (Martin&Keogh&Mandy 1998y)
! T           => Temperature [K]
!-------------------------------------------------------------------------
  implicit none  
  real(dp),intent(in)::T
  real(dp)::comp_Beta_H2H2, kbT
!-------------------------------------------------------------------------
	kbT=8.618*T !eV
  comp_Beta_H2H2 = &
                2.519d-5*(kbt**4.1881)*exp(-0.1731d0/kbT)/(1.d0+2.1347d0*kbT)**5.6881
END FUNCTION comp_Beta_H2H2

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION comp_Beta_HI(T)

! Returns collisional ionization rate [cm3 s-1] of HI (Maselli&'03)
! T           => Temperature [K]
!-------------------------------------------------------------------------
  implicit none  
  real(dp),intent(in)::T
  real(dp)::comp_Beta_HI, T5
!-------------------------------------------------------------------------
  T5 = T/1.d5
  comp_Beta_HI = &
                5.85d-11 * sqrt(T) / (1.d0+sqrt(T5)) * exp(-157809.1d0/T)
END FUNCTION comp_Beta_HI

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION comp_dBeta_dT_HI(T)

! Returns T derivative of collisional ionization rate [cm3 s-1] of HI 
! (Maselli&'03)
! T           => Temperature [K]
!-------------------------------------------------------------------------
  implicit none  
  real(dp),intent(in)::T
  real(dp)::f, T5, comp_dBeta_dT_HI
!-------------------------------------------------------------------------
  T5 = T/1.d5
  f = 1.+sqrt(T5)
  comp_dBeta_dT_HI = 5.85d-11 * sqrt(T) / f * exp(-157809.1d0/T)         &
                * (0.5/T + 157809.1d0/T**2 - .5/(sqrt(1.d5*T)*f))
END FUNCTION comp_dBeta_dT_HI

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION comp_Beta_HeI(T)

! Returns collisional ionization rate [cm3 s-1] of HeI (Maselli&'03)
! T           => Temperature [K]
!-------------------------------------------------------------------------
  implicit none  
  real(dp),intent(in)::T
  real(dp)::comp_Beta_HeI, T5
!-------------------------------------------------------------------------
  T5 = T/1.d5
  comp_Beta_HeI = &
                2.38d-11 * sqrt(T) / (1.d0+sqrt(T5)) * exp(-285335.4d0/T)
END FUNCTION comp_Beta_HeI

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION comp_Beta_HeII(T)

! Returns collisional ionization rate [cm3 s-1] of HeII (Maselli&'03)
! T           => Temperature [K]
!-------------------------------------------------------------------------
  implicit none  
  real(dp),intent(in)::T
  real(dp)::comp_Beta_HeII, T5
!-------------------------------------------------------------------------
  T5 = T/1.d5
  comp_Beta_HeII = &
                5.68d-12 * sqrt(T) / (1.d0+sqrt(T5)) * exp(-631515.d0/T)
END FUNCTION comp_Beta_HeII

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION comp_collExrate_HI(TK)

! Gives Collisional Excitation rate coefficient for the 1->2
! transition in hydrogen atom energies, in [cm3 s-1], according to
! Callaway, Unnikrishnan & Oza 1987.
! TK      => Temperatue in Kelvin
!-------------------------------------------------------------------------
  real(dp),intent(in)::TK
  real(dp)::comp_collExrate_HI
  real(dp),parameter ::kB    = 1.38062d-16       ! Boltzm.const. [erg K-1]
!-------------------------------------------------------------------------
  comp_collExrate_HI =                                                   &
       2.41d-6/sqrt(TK) * (TK/1.d4)**0.22 * exp(-1.63d-11/(kb*TK))
END FUNCTION comp_collExrate_HI

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
FUNCTION compCoolrate(T, ne, nH, nN, nI, aexp, dcooldT, rt_OTSA)

! Compute cooling rate in a cell
! T        => Cell emperature [K]
! nN       => Neutral abundances
! nI       => Ionized abundances
! nH       => Hydrogen number density [cm-3]  
! aexp     => Cosmic expansion
! dcooldT  <= Temperature derivative of the rate
! dcooldx  <= Ionized fraction derivative of the rate
! returns:    Cooling rate [erg s-1 cm-3]  
!-------------------------------------------------------------------------
  implicit none  
  real(dp)::T, ne, nH, aexp
  real(dp),dimension(nIons)::nN,nI
  logical::RT_OTSA
  real(dp)::compCoolrate, dcooldT!----------------------------------------
  real(dp)::T2, sqrtT, sqrt5T, fT5, sqrtT_fT5, f, Ta, fT3
  real(dp)::laHII, laHeII, laHeIII 
  real(dp)::ci_HI, ci_HeI, ci_HeII, dci_HI, dci_HeI, dci_HeII
  real(dp)::r_HII, r_HeII, r_HeIII, dr_HII, dr_HeII, dr_HeIII
  real(dp)::ce_HI, ce_HeI, ce_HeII, dce_HI, dce_HeI, dce_HeII
  real(dp)::bre, brefac, dbre, com, dcom, die, ddie
  real(dp)::lterleft,lterright,lter,lter_hi,lter_h2
  real(dp)::ltevleft,ltevright,ltev,ltev_hi,ltev_h2
  real(dp)::lowrleft,lowrright,lowr_hi,lowr_h2
  real(dp)::lowvleft_hi,lowvright_hi,lowv_hi,lowv_h2
  real(dp)::dlter,dlter_hi,dlter_h2,dltev,dltev_hi,dltev_h2
  real(dp)::dlowr_hi,dlowr_h2,dlowv_hi,dlowv_h2
  real(dp)::TT,fTT3,lowtot_hi,lowtot_h2,ltetot
  real(dp)::cool_hicol,cool_h2col,coolH2
  real(dp)::dcoolr_hi,dcoolv_hi,dcoolr_h2,dcoolv_h2,dcoolH2
  real(dp),parameter::kb=1.3806d-16        ! Boltzmann constant [ergs K-1]
!-------------------------------------------------------------------------
  T2       = T**2
  sqrtT    = sqrt(T)
  sqrt5T   = sqrt(1.d5*T)
  fT5      = 1.+sqrt(T/1.d5)

  ! Coll. Ionization Cooling from Cen 1992 (via Maselli et al 2003)*******
  sqrtT_fT5 = sqrtT/fT5
  f    = 0.5/T - 0.5/(sqrt5T+T)
  ci_HI   = 1.27d-21 * sqrtT_fT5 * exp(-157809.1/T)     * ne   * nN(ixHII)
  dci_HI  = ci_HI *   ( f + 157809.1/T2 )                             !dHI
  if(isHe) then
     ci_HeI  = 9.38d-22 * sqrtT_fT5 * exp(-285335.4/T)  * ne  * nN(ixHeII)
     ci_HeII = 4.95d-22 * sqrtT_fT5 * exp(-631515. /T)  * ne * nN(ixHeIII)
     dci_HeI = ci_HeI *  ( f + 285335.4/T2 )                         !dHeI
     dci_HeII= ci_HeII * ( f + 631515.0/T2 )                        !dHeII
  endif

  ! Recombination Cooling (Hui&Gnedin'97)*********************************
  laHII  = 315614./T                                             
  if(isHe) then
     laHeII  = 570670./T                                            
     laHeIII  = 1263030./T
  endif
  if(.not. rt_otsa) then ! Case A-----------------------------------------
     f       = 1.d0+(laHII/0.541)**0.502                         
     r_HII   = 1.778d-29 * laHII**1.965 / f**2.697 * T    * ne * nI(ixHII)
     dr_HII  = r_HII/T * ( - 0.965 + 1.35389*(f-1.)/f )              !dHII
 
     if(isHe) then
        r_HeII  = 3.d-14 * laHeII**0.654 * kb * T        * ne * nI(ixHeII)
        dr_HeII = r_HeII / T * 0.346                                !dHeII
        f       = 1.d0+(laHeIII/0.541)**0.502                      
        r_HeIII = 14.224d-29 * laHeIII**1.965 / f**2.697 *T*ne*nI(ixHeIII)
        dr_HeIII= r_HeIII/T * ( - 0.965 + 1.35389*(f-1.)/f )       !dHeIII
     endif
  else ! Case B-----------------------------------------------------------
     f      = 1.d0+(laHII/2.25)**0.376                        
     r_HII  = 3.435d-30 * laHII**1.97 / f**3.72 * T       * ne * nI(ixHII)
     dr_HII = r_HII/T * ( - 0.97 + 1.39827*(f-1.)/f )                !dHII

     if(isHe) then
        r_HeII  = 1.26d-14 * laHeII**0.75 * kb * T       * ne * nI(ixHeII)
        dr_HeII = r_HeII / T * 0.25                                 !dHeII
        f       = 1.d0+(laHeIII/2.25)**0.376                       
        r_HeIII = 27.48d-30 * laHeIII**1.97 / f**3.72 * T * ne*nI(ixHeIII)
        dr_HeIII =r_HeIII/T * ( - 0.97 + 1.39827*(f-1.)/f )        !dHeIII
     endif
  endif

  ! Collisional excitation cooling from Cen'92****************************
  ce_HI  = 7.5d-19 * exp(-118348./T) / fT5                * ne * nN(ixHII)
  f = 0.5/(sqrt5T+T)
  dce_HI  = ce_HI * ( 118348./T2 - f)                                !dHI
  if(isHe) then
     ce_HeI = 9.10d-27 * exp(-13179./T) / fT5 * T**(-0.1687)*ne*nN(ixHeII)
     ce_HeII = 5.54d-17 * exp(-473638./T) / fT5*T**(-0.397)*ne*nN(ixHeIII)
     dce_HeI  = ce_HeI * ( 13179./T2 - f - 0.1687/T)                 !dHeI
     dce_HeII  = ce_HeII * ( 473638./T2 - f - 0.397/T)              !dHeII
  endif

  ! Bremsstrahlung from Osterbrock & Ferland 2006*************************
  brefac = 1.42d-27 * 1.5 * sqrtT  * ne
  bre  = brefac * nI(ixHII)                               
  if(isHe) then
     bre = bre + brefac * ( nI(ixHeII) + 4. * nI(ixHeIII) )
  endif
  dbre = bre * 0.5 / T

  ! Compton Cooling from Haimann et al. 96, via Maselli et al.************
  Ta   = 2.727/aexp                          
  com  = 1.017d-37 * Ta**4 * (T-Ta)                                   * ne    
  dcom  = com / (T-Ta)   

  ! Dielectronic recombination cooling, from Black 1981*******************
  if(isHe) then
     f = 1.24D-13*T**(-1.5d0)*exp(-470000.d0/T)         * ne * nN(ixHeIII)
     die = f*(1.D0+0.3d0*exp(-94000.d0/T))
     ddie = die/T2*(564000.-1.5*T) - f*94000./T2
  endif

  ! H2 cooling Hallenbach McKee (1979) + Halle Combes (2012) in cgs*******
  if(isH2) then

!!$     ! LTE cooling
!!$     fT3=T/1.d3
!!$     lterleft=9.5d-22*fT3**3.76*exp(-1.0d0*(0.13d0/fT3)**3)/(1.0d0+0.12d0*fT3**2.1)
!!$     lterright=3.0d-24*exp(-0.51d0/fT3)
!!$     lter=lterright+lterleft
!!$     ltevleft=6.7d-19*exp(-5.86d0/fT3)
!!$     ltevright=1.6d-18*exp(-11.7d0/fT3)
!!$     ltev=ltevleft+ltevright
!!$     ltetot=lter+ltev
  
     ! Collisional cooling
     lowrleft=1.25*exp(-1.70d2/T)*2.35d-14
     lowrright=1.75*exp(-5.05d2/T)*6.97d-14
     lowr_hi=lowrleft*gamma_hi(T,2.d0)+lowrright*gamma_hi(T,3.d0)
     lowr_h2=lowrleft*gamma_h2(T,2.d0)+lowrright*gamma_h2(T,3.d0)
     lowvleft_hi=exp(-5860.0/T)*8.09d-13*1.0d-12*sqrt(T)*exp(-1.0d3/T)
     lowvright_hi=exp(-11720.0/T)*1.6d-12*1.6d-12*sqrt(T)*exp(-1.0*(4.0d2/T)**2)
     lowv_hi=lowvleft_hi+lowvright_hi
     lowv_h2=exp(-5860.0/T)*8.09d-13*1.4d-12*sqrt(T)*exp(-1.2d4/(T+1.2d3))
     lowtot_hi=lowr_hi+lowv_hi
     lowtot_h2=lowr_h2+lowv_h2

     coolH2 = nN(ixHI)*(nN(ixHII)*lowtot_hi+nN(ixHI)*lowtot_h2)
     
     ! Same thing for the temperature derivative
     TT=T*1.0001

!!$     ! LTE cooling
!!$     fTT3=TT/1.d3
!!$     lterleft=9.5d-22*fTT3**3.76*exp(-1.0d0*(0.13d0/fTT3)**3)/(1.0d0+0.12d0*fTT3**2.1)
!!$     lterright=3.0d-24*exp(-0.51d0/fTT3)
!!$     lter=lterright+lterleft
!!$     ltevleft=6.7d-19*exp(-5.86d0/fTT3)
!!$     ltevright=1.6d-18*exp(-11.7d0/fTT3)
!!$     ltev=ltevleft+ltevright
!!$     ltetot=lter+ltev
     
     ! Collisional cooling
     lowrleft=1.25*exp(-1.70d2/TT)*2.35d-14
     lowrright=1.75*exp(-5.05d2/TT)*6.97d-14
     lowr_hi=lowrleft*gamma_hi(TT,2.d0)+lowrright*gamma_hi(TT,3.d0)
     lowr_h2=lowrleft*gamma_h2(TT,2.d0)+lowrright*gamma_h2(TT,3.d0)
     lowvleft_hi=exp(-5860.0/TT)*8.09d-13*1.0d-12*sqrt(TT)*exp(-1.0d3/TT)
     lowvright_hi=exp(-11720.0/TT)*1.6d-12*1.6d-12*sqrt(TT)*exp(-1.0*(4.0d2/TT)**2)
     lowv_hi=lowvleft_hi+lowvright_hi
     lowv_h2=exp(-5860.0/TT)*8.09d-13*1.4d-12*sqrt(TT)*exp(-1.2d4/(TT+1.2d3))
     lowtot_hi=lowr_hi+lowv_hi
     lowtot_h2=lowr_h2+lowv_h2
     
     dcoolH2 = nN(ixHI)*(nN(ixHII)*lowtot_hi+nN(ixHI)*lowtot_h2)

     dcoolH2 = (dcoolH2 - coolH2) / (0.0001*T)

  endif

  ! Overall Cooling*******************************************************
  compCoolrate = ci_HI + r_HII + ce_HI + com + bre
  if(isHe) compCoolrate = compCoolrate &
                 +ci_HeI +r_HeII +ce_HeI +die +ci_HeII  +r_HeIII  +ce_HeII
  if(isH2) compCoolrate = compCoolrate + coolH2

  dCooldT = dci_HI + dr_HII + dce_HI + dcom + dbre
  if(isHe) dCooldT = dCooldT  &
            +dci_HeI +dr_HeII +dce_HeI +ddie +dci_HeII +dr_HeIII +dce_HeII
  if(isH2) dCooldT = dCooldT + dcoolH2

END FUNCTION compCoolrate

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION gamma_hi(T,J)
  !-------------------------------------------------------------------------
  implicit none  
  real(dp),intent(in)::T,J
  real(dp)::gamma_hi, T3, jconsts
  !-------------------------------------------------------------------------
  T3 = T/1.d3
  jconsts=0.33+0.9*exp(-1.0d0*((J-3.5d0)/0.9d0)**2)
  gamma_hi = &
       jconsts*(1.0d-11*sqrt(T3)/(1.0d0+60.0d0*T3**(-4)) + 1.0d-12*T3)
               
END FUNCTION gamma_hi
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION gamma_h2(T,J)
  !-------------------------------------------------------------------------
  implicit none  
  real(dp),intent(in)::T,J
  real(dp)::gamma_h2, T3, jconsts
  !----------------------------------------------------------
  T3 = T/1.d3
  jconsts=(0.276*J**2)*exp(-1.0d0*(J/3.18d0)**1.7)
  gamma_h2 = jconsts*(3.3d-12 +6.6d-12*T3)
               
END FUNCTION gamma_h2
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION dgamma_hi(T,J)
  !-------------------------------------------------------------------------
  implicit none  
  real(dp),intent(in)::T,J
  real(dp)::dgamma_hi, T3, jconsts
  !-------------------------------------------------------------------------
  T3 = T/1.d3
  jconsts=0.33+0.9*exp(-1.0d0*((J-3.5d0)/0.9d0)**2)
  dgamma_hi = &
       jconsts*(1.d-11*sqrt(T3)/(1.+60.*T3**(-4))*(1./(2.*T3)+ &
			240.*T3**(-5)/(1.+60.*T3**(-4)))+1.d-12)/1.d3
               
END FUNCTION dgamma_hi

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION dgamma_h2(J)
  !-------------------------------------------------------------------------
  implicit none  
  real(dp),intent(in)::J
  real(dp)::dgamma_h2, jconsts
  !-------------------------------------------------------------------------
  jconsts=0.276*(J**2)*exp(-1.0d0*(J/3.18d0)**1.7)
  dgamma_h2 = &
       jconsts*6.6d-12
              
END FUNCTION dgamma_h2
END MODULE coolrates_module       
