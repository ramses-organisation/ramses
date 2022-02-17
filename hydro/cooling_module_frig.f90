!! Cooling from frig version including modifications for extinction
! Authors: Valeska Valdiva, Patrick Hennebelle, Benjamin Godard
!=======================================================================
subroutine solve_cooling_frig(nH,T2,zsolar,boost,dt,deltaT2,ncell)
!=======================================================================
  implicit none
  ! BRIDGE FUNCTION WITH SAME INTERFACE AS solve_cooling
  ! Input/output variables to this function
  ! nH - hydrogen number density in PHYSICAL units
  ! T2 - temperature / mu in PHYSICAL units
  ! zsolar - Metallicity in solar units (Zphys / 0.02)
  ! boost - raditation boost - exp(-100*nH) if self_shielding=True
  ! dt - cooling timestep in seconds
  ! deltaT2 - temperature change in K/mu (??)
  ! ncell - number of elements in the vector
  integer::ncell
  real(kind=8)::dt
  real(kind=8),dimension(1:ncell)::nH,T2,deltaT2,zsolar,boost
  ! Input/output variables to analytic function calc_temp 
  real(kind=8)::NN,TT, dt_tot_unicode
  ! Temporary variables
  integer::i
  real(kind=8)::TT_ini, mu
  ! Units
  real(kind=8) :: scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  ! HARD-CODED mu TO MAKE TEMPERATURE AGREE WITH HENNEBELLE CODE
  ! TC: this is ugly, we can do better
  mu = 1.4
  scale_T2 = scale_T2 * mu
  do i=1,ncell
     NN = nH(i) ! NOTE!! THE CODE BELOW ASSUMES scale_nH=1 !!
                ! SO WE LEAVE THIS AS IT IS TO KEEP UNITS CONSISTENCY
     TT = T2(i) / scale_T2
     TT_ini = TT
     dt_tot_unicode = dt / scale_t
     call calc_temp(NN,TT,dt_tot_unicode)
     deltaT2(i) = (TT - TT_ini) * scale_T2
  end do
end subroutine solve_cooling_frig
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine  calc_temp(NN,TT,dt_tot_unicode)
    use amr_parameters
    use hydro_commons
    use constants, only:kB
    implicit none

    integer :: n,i,j,k,idim, iter,ii
    real(dp) :: dt, dt_tot, temps, dt_max
    real(dp) :: rho,temp,dt_tot_unicode
    real(dp) :: alpha_ct,mu
    real(dp) :: NN,TT, TTold, ref,ref2,dRefdT, eps, vardt,varrel, dTemp
    real(dp) :: rhoutot2
    real(dp) :: scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
    ! HARD-CODED mu TO MAKE TEMPERATURE AGREE WITH HENNEBELLE CODE
    ! TC: this is ugly, we can do better
    mu = 1.4
    !
    ! cgs units are used here
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
    scale_T2 = scale_T2 * mu

    if( TT .le. 0.) then
        TT = 50. / scale_T2
        return
    endif

    vardt = 10.**(1./10.); varrel = 0.2

    dt_tot = dt_tot_unicode * scale_t
    TT     = TT * scale_T2

    if (NN .le. smallr) then
        if( NN .le. 0)  write(*,*) 'prob dens',NN
        NN = smallr  
    endif

    alpha_ct = NN*kB/(gamma-1.)

    ! eps - a small offset of T to find gradient in T
    eps = 1d-5

    iter  = 0 ; temps = 0.
    do while ( temps < dt_tot)
        if (TT .lt.0) then
            write(*,*) 'prob Temp',TT, NN
            NN = max(NN,smallr)
            TT = min(4000./NN,8000.)  !2.*4000. / NN
        endif

        TTold = TT

        ! Calculate cooling rate
        !NN is assumed to be in cc and TT in Kelvin
        if (TT < 10035.d0) then
            call cooling_low(TT,NN,ref)
            call cooling_low(TT*(1d0+eps),NN,ref2)
        else
            call cooling_high(TT,NN,ref)
            call cooling_high(TT*(1d0+eps),NN,ref2)
        end if
        
        ! dT = T*(1+eps)-T = eps*T
        dRefdT = (ref2-ref)/(TT*eps)

        ! TODO STG - COPY THIS FUNCTION UP TO HERE, USE ref, drefdT TO 
        !            REPLACE rt_cmp_metals SOMEHOW

        if (iter == 0) then
            if (dRefDT .ne. 0.) then
                dt = abs(1.0E-1 * alpha_ct/dRefDT)
            else
                dt = 1.0E-1 * dt_tot
            endif
            dt_max = dt_tot - temps
            if (dt > 0.7*dt_max) dt = dt_max*(1.+1.0E-12)
        endif

        dTemp = ref/(alpha_ct/dt - dRefdT)

        eps = abs(dTemp/TT)
        if (eps > 0.2) dTemp = 0.2*TTold*dTemp/abs(dTemp)

        TT = TTold + dTemp
        if (TT < 0.) then
            write(*,*) 'Temperature negative !!!'
            write(*,*) 'TTold,TT   = ',TTold,TT
            write(*,*) 'rho   = ',rho
            TT = 100.  !*kelvin
        endif
        iter = iter + 1
        temps = temps + dt
        dt = vardt*varrel*dt/Max(vardt*eps, varrel)
        dt_max = dt_tot - temps
        if (dt > 0.7*dt_max) dt = dt_max*(1.+1.0E-12)
    enddo

    !!now convert temperature in code units
    TT = TT / scale_T2

end subroutine calc_temp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cooling_high(T,n,ref)
    use amr_parameters
    implicit none

    real(dp) :: T,n,ref
    real(dp) :: cold,hot,logT

    ! cooling rate based on Dopita and Sutherland

    logT=log10(T)

    if (logT .LT. 4.0) then
        cold=0.1343*logT**3-1.3906*logT**2+5.1554*logT-31.967
    else if (logT .LT. 4.25) then
        cold=12.64*logT-75.56
    else if (logT .LT. 4.35) then
        cold=-0.3*logT-20.565
    else if (logT .LT. 4.9) then
        cold=1.745*logT-29.463
    else if (logT .LT. 5.4) then
        cold=-20.9125
    else if (logT .LT. 5.9) then
        cold=-1.795*logT-11.219
    else if (logT .LT. 6.2) then
        cold=-21.8095
    else if (logT .LT. 6.7) then
        cold=-1.261*logT-13.991
    else
        cold=-22.44
    endif

    cold=-1.0*10.0**(cold)

    ref=(n**2)*(cold)

end subroutine cooling_high
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cooling_low(T,n,ref)
  use amr_parameters
  implicit none

  real(dp) :: T,n,ref
  real(dp) :: cold,hot,cold_cII,cold_o,cold_h,cold_cII_m,cold_o_m,cold_rec
  real(dp) :: param,G0,epsilon,bet,x,x_ana,ne   ! x is the ionisation rate

  ! cooling and heating function computed from the cooling of 
  ! chemical elements

  ! Carbon abondance 3.5 10-4, depletion 0.4

  !!! We compute the ionisation
  !!! We assume that if x is over 1.d-4 then it is dominated by oxygen
  !!! and that the minimal value is given by the carbon and is 
  !!! 3.5 1.d-4 * depletion * density
  
  !!! For the electrons due to hydrogen we take the formula
  !!! from  Wolfire et al. 2003 appendice C2.
  !!! The ionization rate is set to 1.d-16 G0'=GO/1.7
  !!! Z'd = 1 et phipah=0.5

  ne = 2.4d-3*((T/100d0)**0.25d0)/0.5d0 ! formula C15 of Wolfire et al. 2003

  ! Analytic ionisation in absence of photoionisation
  x_ana = ne / N   ! ionisation
  x_ana = min(x_ana,0.1d0)
  x_ana = max(x_ana,3.5d-4*0.4d0)
  x = x_ana ! (Different variables in case we need photoionised values later)

 ! NOTE - HERE WE USE THE NON-PHOTOIONISED RATES AS THIS MIGHT 
 !        BE TOO HIGH AT x=1
  cold_cII =  92. * 1.38E-16 * 2. * (2.8E-7* ((T/100.)**(-0.5))*x_ana + 8.E-10*((T/100.)**(0.07))) &
       * 3.5E-4 * 0.4 * exp(-92./ T)

  ! oxygen-prompted cooling
  ! abundance 8.6 10-4 depletion 0.8
  cold_o = 1.E-26 * sqrt(T) * (24. * exp(-228./ T) + 7. * exp(-326./ T) )
  ! take oxygen abundance into account
  cold_o = cold_o * 4.5E-4

  ! Hydrogen-prompted cooling
  ! formula from Spitzer 1978
  ! NOTE - RT function handles hydrogen cooling out of equilibrium
  cold_h = 0d0
  if (.not. rt) then
     cold_h = 7.3E-19 * x * exp(-118400./ T )
  endif

  ! cooling from metastables metal lines
  ! formulas from Hollenbach and McKee 1989 (ApJ 342, 306)

  ! Ionised carbon, 1 transition 2P 4P
  ! weight is 1
  ! 2P->4P :
  ! The excitation coefficients depends on the temperature above 10^4 K
  ! les expressions des coefficients d'excitation ont une dependance
  ! abondance 3.5 d-4 depletion 0.4

  cold_cII_m = 6.2d4 * 1.38d-16 * 1.d0 * &    !transition 2P->4P
             ( 2.3d-8* (T/10000.)**(-0.5) * x + 1.d-12 ) *exp(-6.2d4 / T) &
             * 3.5d-4 * 0.4

  if ( T .le. 1.d4 ) then
     cold_o_m = 2.3d4 * 1.38d-16 / 3.d0 * &
             ( 5.1d-9 * (T/10000.)**(0.57) * x + 1.d-12) *exp(-2.3d4/T)
  
     cold_o_m = cold_o_m + &
                4.9d4 * 1.38d-16 / 3.d0  * &
             ( 2.5d-9 * (T/10000.)**(0.57) * x + 1.d-12) *exp(-4.9d4/T)
  
     cold_o_m = cold_o_m + &
                2.6d4 * 1.38d-16 * 1.d0  * &
             ( 5.2d-9 * (T/10000.)**(0.57) * x + 1.d-12) *exp(-2.6d4/T)
  else
     cold_o_m = 2.3d4 * 1.38d-16 / 3.d0 * &
             ( 5.1d-9 * (T/10000.)**(0.17) * x + 1.d-12) *exp(-2.3d4/T)
  
     cold_o_m = cold_o_m + &
                4.9d4 * 1.38d-16 / 3.d0  * &
             ( 2.5d-9 * (T/10000.)**(0.13) * x + 1.d-12) *exp(-4.9d4/T)

     cold_o_m = cold_o_m + &
                2.6d4 * 1.38d-16 * 1.d0  * &
             ( 5.2d-9 * (T/10000.)**(0.15) * x + 1.d-12) *exp(-2.6d4/T)
  endif

  ! oxigen abundance 
  cold_o_m = cold_o_m * 4.5d-4

  ! sum of the cooling terms
  cold = cold_cII  + cold_h  + cold_o  + cold_o_m +  cold_cII_m

  !!!! Computation of the heating term
  !!! Heating on grains is taken into account
  !!! formula 1 et 2  of Wolfire et al. 1995

  !!!! G0 is the UV flux compared to the one given by Habing et Draine
  G0 = 1./1.7

  param = G0 * sqrt(T)/(n*x)
  epsilon = 4.9E-2 / (1. + (param/1925.)**0.73)
  epsilon  = epsilon + 3.7E-2 * (T/1.E4)**0.7 / (1. + (param/5.E3) )

  hot = 1.E-24 * epsilon

  ! for a UV flux of G0/1.7
  hot = hot * G0

  ! recombination cooling on positively charged grains
  bet = 0.74/(T**0.068)
  cold_rec = 4.65E-30*(T**0.94)*(param**bet)*x

  ref = hot*n - (n**2)*(cold + cold_rec)

end subroutine cooling_low
!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine hot_cold_2(T,n,ref,dRefDT,XH2)    
  use amr_parameters
#ifdef RT
  use rt_parameters,only: isH2,iIons
#endif
  implicit none
  
  real(dp), intent (in)    :: T,n
  real(dp), intent (out)   :: ref,dRefDT

  real(dp)                 :: G0,zeta_p,phi_pah,x_cp,x_o,eps
  real(dp)                 :: T_1,T_2,ne_1,ne_2,x_1,x_2
  real(dp)                 :: cold_cII_1,cold_o_1,cold_h_1,cold_rec_1,cold_1
  real(dp)                 :: cold_cII_2,cold_o_2,cold_h_2,cold_rec_2,cold_2
  real(dp)                 :: hot_ph_1,hot_ph_2,hot_cr,hot_1,hot_2
  real(dp)                 :: ref_1,ref_2

  real(dp)                 :: XH2
  real(dp)                 :: x_cO=0.,cold_mol_1=0.,cold_mol_2=0.

  ! ======================================================================================
  ! 1 - Definition of the local UV field 
  !                       local C+ and O fraction
  !                       local electronic fraction 
  !     UV field is expressed in Habing unit (1.274e-4 erg cm-2 s-1 sr-1)
  !     p_UV is a scaling parameter defined in the namelist
  !     the shielding coefficient is computed using Valeska's method
  ! ======================================================================================
  G0 = 1.0_dp * p_UV
  ! TC: can't we get the UV flux from RT?

  ! ---------------------------------------------------
  ! C+ and O fraction
  ! - We assume that all Carbon is in ionized form 
  !   with a depletion onto grains of 0.40 => 1.4e-4 (Wolfire et al. 2003, Tab. 1)
  ! - We assume that all Oxygen is in atomic  form 
  !   with a depletion onto grains of 0.37 => 3.2e-4 (Wolfire et al. 2003, Tab. 1)
  ! ---------------------------------------------------
  x_cp = 3.5e-4_dp * 0.40_dp
  x_o  = 8.6e-4_dp * 0.37_dp

#ifdef RT
  if(isH2) then
     !calculate CO abundance using a simple prescription and recalculate the x_cp abundance
     CALL compute_Cp_CO_NL(n,XH2,G0,x_o,x_cp,x_co)
  endif
#endif

  ! ---------------------------------------------------
  ! ionization and PAH recombination parameter
  ! set to use the formalism of Wolfire et al. (2003)
  ! ---------------------------------------------------
  zeta_p  = 2.0_dp
  phi_pah = 0.5_dp

  ! ---------------------------------------------------
  ! temperatures at which cooling & heating are estimated 
  ! ---------------------------------------------------
  T_1 = T
  eps = 1.0e-5_dp
  T_2 = T * (1.0_dp + eps)

  ! ---------------------------------------------------
  ! electron density and electronic fraction
  ! ---------------------------------------------------
  CALL ELEC_DENS(zeta_p, G0, T_1, phi_pah, x_cp, N, ne_1)
  CALL ELEC_DENS(zeta_p, G0, T_2, phi_pah, x_cp, N, ne_2)
  x_1 = ne_1 / N
  x_2 = ne_2 / N

  ! ======================================================================================
  ! 2 - set the cooling functions
  !     a - hyperfine transitions at small temperature of CII and OI
  !         metastable lines of CII and OI (NOT USED)
  !     b - H cooling by excitation of its electronic lines
  !     c - cooling by recombination on positively charged grains
  ! ======================================================================================

  CALL COOL_CP(T_1, x_cp, x_1, cold_cII_1)
  CALL COOL_CP(T_2, x_cp, x_2, cold_cII_2)

  CALL COOL_O(T_1, x_o, cold_o_1)
  CALL COOL_O(T_2, x_o, cold_o_2)

  cold_h_1 = 0d0
  cold_h_2 = 0d0
  !if rt this cooling is taken into account in the rt module
  if (.not. rt) then
     CALL COOL_H(T_1, x_1, cold_h_1)
     CALL COOL_H(T_2, x_2, cold_h_2)
  endif

  CALL COOL_REC(G0, T_1, phi_pah, x_1, N, cold_rec_1)
  CALL COOL_REC(G0, T_2, phi_pah, x_2, N, cold_rec_2)

#ifdef RT
  if(isH2) then
  !calculate molecular cooling - in principle no need for RT and isH2 but for the sake of coherence
  !it is required for now - can be changed
  CALL cool_goldsmith(T_1,N,cold_mol_1)
  CALL cool_goldsmith(T_2,N,cold_mol_2)
  endif
#endif

  ! Sum all cooling functions
  cold_1 = cold_cII_1  + cold_o_1 + cold_h_1 + cold_rec_1 + cold_mol_1
  cold_2 = cold_cII_2  + cold_o_2 + cold_h_2 + cold_rec_2 + cold_mol_2

  ! ======================================================================================
  ! 3 - set the heating rates
  !     a - photoelectric effect
  !     b - cosmic rays (for dark cloud cores, Goldsmith 2001, Eq. 3)
  ! ======================================================================================

  CALL HEAT_PH (G0, T_1, phi_pah, x_1, N, hot_ph_1)
  CALL HEAT_PH (G0, T_2, phi_pah, x_2, N, hot_ph_2)

  !corrected CR as value seems to be much higher than orif-ginally estimated (McCall+2003)
  hot_cr = 5.0E-27_dp

  ! Sum all heating functions
  hot_1 = hot_ph_1 + hot_cr
  hot_2 = hot_ph_2 + hot_cr

  ! ======================================================================================
  ! 4 - compute the net heating rate (in erg cm-3 s-1) and the output variables
  ! ======================================================================================
  ref_1 = n * hot_1 - n**2.0_dp * cold_1
  ref_2 = n * hot_2 - n**2.0_dp * cold_2

  ref    = ref_1
  dRefDT = (ref_2-ref_1) / T / eps

end subroutine hot_cold_2

SUBROUTINE ELEC_DENS(zeta_p, G0_p, T, phi, xcp, nH, ne)
    !---------------------------------------------------------------------------
    ! Compute the local electron density based on Wolfire et al. (2003, eq. C9) 
    ! and adding to this equation the abundance of C+
    ! input variables :
    !     zeta_p -> ionization parameter (in units of 1.3e-16 s-1)
    !     G0_p   -> local UV radiation field (in Habing unit)
    !     T      -> local temperature
    !     phi    -> PAH recombination factor
    !     xcp    -> C+ abundance (relative to nH)
    !     nH     -> proton density
    ! output variables :
    !     ne     -> electron density (in cm-3)
    !---------------------------------------------------------------------------
    USE amr_parameters, only:dp
    IMPLICIT none

    REAL (KIND=dp), INTENT(IN)  :: zeta_p,G0_p,T,phi,xcp,nH
    REAL (KIND=dp), INTENT(OUT) :: ne

    ! Careful. The formula of Wolfire is given as functions of parameters which 
    ! are scaling of those used in their standard model of the solar neighbourghood 
    ! ===> in this model G0 = 1.7 (in Habing units)
    ne = 2.4e-3_dp * (zeta_p)**0.5_dp * (T/100_dp)**0.25_dp * (G0_p/1.7)**0.5_dp / phi + xcp * nH

    ! limit the electron fraction to 0.1 to get closer to Wolfire results (Fig. 10 top panels)
    ne = min(ne,0.1*nH)

END SUBROUTINE ELEC_DENS

SUBROUTINE COOL_CP(T, xcp, xe, cool)
    !---------------------------------------------------------------------------
    ! Compute the cooling rate due to the hyperfine 
    ! (and metastable - NOT USED HERE) lines of C+
    ! input variables :
    !     T    -> local temperature
    !     xcp  -> C+ abundance (relative to nH)
    !     xe   -> e- abundance (relative to nH)
    ! output variables :
    !     cool -> cooling rate (erg cm3 s-1) 
    !---------------------------------------------------------------------------
    USE amr_parameters, only:dp
    IMPLICIT none

    REAL (KIND=dp), INTENT(IN)  :: T,xcp,xe
    REAL (KIND=dp), INTENT(OUT) :: cool

    ! Ref : Wolfire et al. (2003, Eq. C1 & C2)
    cool = ( 2.25e-23_dp + 1.0e-20_dp * (T/100.0_dp)**(-0.5_dp) * xe ) * exp(-92.0_dp / T) * xcp

END SUBROUTINE COOL_CP

SUBROUTINE COOL_O(T, xo, cool)
    !---------------------------------------------------------------------------
    ! Compute the cooling rate due to the hyperfine 
    ! (and metastable - NOT USED HERE) lines of O
    ! input variables :
    !     T    -> local temperature
    !     xo   -> O abundance (relative to nH)
    ! output variables :
    !     cool -> cooling rate (erg cm3 s-1) 
    !---------------------------------------------------------------------------
    USE amr_parameters, only:dp
    IMPLICIT none

    REAL (KIND=dp), INTENT(IN)  :: T,xo
    REAL (KIND=dp), INTENT(OUT) :: cool

    ! Ref : Wolfire et al. (2003, Eq. C3)
    cool = 7.81e-24_dp * (T/100.0_dp)**(0.4_dp) * exp(-228.0_dp/T) * xo

END SUBROUTINE COOL_O

SUBROUTINE COOL_H(T, xe, cool)
    !---------------------------------------------------------------------------
    ! Compute the cooling rate due to the electronic lines of H
    ! Ref : taken from Spitzer (1978)
    ! input variables :
    !     T    -> local temperature
    !     xe   -> e- abundance (relative to nH)
    ! output variables :
    !     cool -> cooling rate (erg cm3 s-1) 
    !---------------------------------------------------------------------------
    USE amr_parameters, only:dp
    IMPLICIT none

    REAL (KIND=dp), INTENT(IN)  :: T,xe
    REAL (KIND=dp), INTENT(OUT) :: cool

    cool = 7.3e-19_dp * xe * exp(-118400.0_dp / T)

END SUBROUTINE COOL_H

!TC: carefull, this is also defined in cooling_module.f90
SUBROUTINE COOL_REC(G0_p, T, phi, xe, nH, cool)
    !---------------------------------------------------------------------------
    ! Compute the cooling rate due to the recombination
    ! of electrons on positively charged PAH
    ! input variables :
    !     G0_p -> local UV radiation field (Draine's unit)
    !     T    -> local temperature
    !     phi  -> PAH recombination factor
    !     xe   -> e- abundance (relative to nH)
    !     nH   -> proton density
    ! output variables :
    !     cool -> cooling rate (erg cm3 s-1) 
    !---------------------------------------------------------------------------
    USE amr_parameters, only:dp
    IMPLICIT none

    REAL (KIND=dp), INTENT(IN)  :: G0_p,T,phi,xe,nH
    REAL (KIND=dp), INTENT(OUT) :: cool
    REAL (KIND=dp)              :: param,bet

    param = G0_p * sqrt(T) / (nH * xe * phi)
    bet   = 0.74_dp / (T**0.068_dp)

    ! Ref : Wolfire et al. (2003, Eq. 21)
    cool = 4.65e-30_dp * T**0.94_dp * param**bet * xe * phi

END SUBROUTINE COOL_REC

SUBROUTINE HEAT_PH(G0_p, T, phi, xe, nH, heat)
    !---------------------------------------------------------------------------
    ! Compute the heating rate due to the photoelectric effect
    ! Ref : Wolfire et al. (2003, Eqs. 19 & 20)
    ! input variables :
    !     G0_p -> local UV radiation field (in Habing units)
    !     T    -> local temperature
    !     phi  -> PAH recombination factor
    !     xe   -> e- abundance (relative to nH)
    !     nH   -> proton density
    ! output variables :
    !     heat -> cooling rate (erg cm3 s-1) 
    !---------------------------------------------------------------------------
    USE amr_parameters, only:dp
    IMPLICIT none

    REAL (KIND=dp), INTENT(IN)  :: G0_p,T,phi,xe,nH
    REAL (KIND=dp), INTENT(OUT) :: heat
    REAL (KIND=dp)              :: param,epsilon

    param = G0_p * sqrt(T) / (nH * xe * phi)

    epsilon = 4.9e-2_dp / ( 1.0_dp + (param / 1925.0_dp)**0.73_dp ) &
            + 3.7e-2_dp / ( 1.0_dp + (param / 5000.0_dp)          ) * (T / 1.0e4_dp)**0.7_dp

    heat = 1.3E-24 * epsilon * G0_p

END SUBROUTINE HEAT_PH
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Molecular cooling from table 2 of Goldsmith 2001
! Temp: temperature in K
! ndens: density in cc
! cool_mol_dust: molecular and dust cooling 
subroutine cool_goldsmith(Temp,ndens,cool_mol_dust)
  use amr_parameters, only:dp
  implicit none

  real(dp):: Temp,ndens
  real(dp)::cool_mol_dust
  integer :: i
  real(dp):: log_n,alpha,beta,cool_mol,cool_dust
  !Table from Goldsmith 2001
  real(dp),dimension(7),parameter:: logn_v= (/2.,2.5,3.,4.,5.,6.,7./)
  real(dp),dimension(7),parameter:: alpha_v=(/6.3e-26,3.2e-25,1.1e-24,5.6e-24,2.3e-23,4.9e-23,7.4e-23/)
  real(dp),dimension(7),parameter:: beta_v=(/1.4,1.8,2.4,2.7,3.,3.4,3.8/)
  !dust temperature is assumed to be 10 K 
  !need to be improved particularly if stellar radiative feedback is accounted for
  real(dp):: Tdust=10.

  !validity range 
  if(ndens <  100. .or. ndens > 1.e7) then 
     cool_mol=0.
     ! TC: this should be cool_mol_dust?
     !     What about the dust contribution below 100 H/cc?
     return
  endif

  log_n = log10(ndens)
  do i = 1, 6
     if (log_n .ge. logn_v(i) .and. log_n .le. logn_v(i+1)) exit
     ! TC: exit necesary?
  enddo
    
  alpha = alpha_v(i+1) * (log_n - logn_v(i)) + alpha_v(i) * (-log_n + logn_v(i+1))
  alpha = alpha / (logn_v(i+1)-logn_v(i))

  beta = beta_v(i) * (log_n - logn_v(i)) + beta_v(i) * (-log_n + logn_v(i+1))
  beta = beta / (logn_v(i+1)-logn_v(i))

  cool_mol = alpha * (Temp/10.)**beta / ndens**2 !division by n^2 because cooling works in cm^3 (rate)
                                                 !while Goldsmith expressed it in cm^-3

  !take into accound the cooling/heating through dust
  cool_dust = 2.e-33 * (Temp-Tdust) * sqrt(Temp/10.)

  cool_mol_dust = cool_mol + cool_dust
                                                                                                       
end subroutine cool_goldsmith
!###########################################################                                          
!###########################################################                                          
!###########################################################                                          
!########################################################### 
! Compute the C+/CO transition at equilibrium from Koyama & Inutuska 2000
! taken from Nelson & Langer 1997
subroutine compute_Cp_CO_NL(ndens,XH2,G0,xO,xCp,xCO)
  use amr_parameters, only:dp
  implicit none

  real(dp):: ndens, XH2, G0, xO, xCp, xCO
  real(dp),parameter::k0=5.d-16 !#cm^3 s^-1
  real(dp),parameter::k1=5.d-10
  ! real(dp):: xCp_0=3.e-4,xO=3.e-4
  real(dp):: xCp_0,beta,gamma_CO,Gamma_CHx

  xCp_0=xCp
  gamma_CO = 1.d-10 * G0 !#s^-1
  Gamma_CHx = 5.d-10 * G0 
    
  beta = k1*xO / (k1*xO+gamma_CO/(XH2*ndens))
  
  xCp = xCp_0 * gamma_CO / (gamma_CO + k0*beta*ndens)
    
  xCO = xCp_0 - xCp
    
end subroutine compute_Cp_CO_NL
