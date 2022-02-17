!! Metal cooling from frig version (Audit & Hennebelle 2005)
!! TC: some of these processes depend on the composition of the gas which
!!     can be obtained from RT
!! Authors: Sam Geen, Patrick Hennebelle, Benjamin Godard
!=========================================================================
! Calculates metal cooling rates using tables
SUBROUTINE rt_metal_cool(Tin,Nin,xin,mu,metal_tot,metal_prime,XH2,aexp)
  use amr_parameters,only:dp
  implicit none
  ! Taken from the equilibrium cooling_module of RAMSES
  ! Compute cooling enhancement due to metals
  ! Tin          => Temperature in Kelvin, divided by mu
  ! Nin          => Hydrogen number density (H/cc)
  ! xin          => Hydrogen ionisation fraction (0->1)
  ! mu           => Average mass per particle in terms of mH
  ! ne           => Electron number density
  ! metal_tot   <=  Metal cooling contribution to de/dt / (nH*ne) [erg s-1 cm-3]
  ! metal_prime <=  d(metal_tot)/dT2 / (nH*ne) [erg s-1 cm-3 K-1]
  real(dp),intent(in)::Tin,Nin,xin,mu,aexp
  real(dp)::T1,T2,nH,cool1,cool2,eps,XH2
  real(dp),intent(out)::metal_tot,metal_prime

  ! Set a reference temperature to calculate gradient wrt T
  eps = 1d-5 ! A small value
  T1 = Tin*mu
  T2 = Tin*(1+eps)*mu
  
  ! Call a function mixing the two cooling functions
  call rt_metal_cool_mashup(T1,Nin,xin,mu,cool1,XH2,aexp)
  call rt_metal_cool_mashup(T2,Nin,xin,mu,cool2,XH2,aexp)
  
  ! Don't cool below 10K to prevent bound errors, but allow heating
  if ((Tin*mu .gt. 10d0) .or. (cool1 .lt. 0d0)) then
     ! Calculate gradient and output
     metal_tot = cool1
     ! T2 = T*(1+eps), so T2-T == eps*T
     metal_prime = (cool2 - cool1) / (Tin * mu * eps)
     ! NOTE !!!! NEED TO MULTIPLY BY nH*ne AFTER THIS IS OVER!!!!
     ! EXCLAMATION MARK EXCLAMATION MARK
  else
     ! Prevent runaway cooling below 10K
     metal_tot = 0d0
     metal_prime = 0d0
  endif

END SUBROUTINE rt_metal_cool

!=========================================================================
! Mixes Patrick Hennebelle's cooling function with Alex Riching's table
! Uses a sigmoid function centred at x=0.1 with a +/-0.05 spread to switch
SUBROUTINE rt_metal_cool_mashup(T,N,x,mu,cool,XH2,aexp)
  use amr_parameters,only:dp
  use cooling_module,only:cmp_metals
  implicit none
  ! Taken from the equilibrium cooling_module of RAMSES
  ! Compute cooling enhancement due to metals
  ! T            => Temperature in Kelvin *with mu included*
  ! N            => Hydrogen number density (H/cc)
  ! x            => Hydrogen ionisation fraction (0->1)
  ! cool         <=  Metal cooling [erg s-1 cm-3]
  ! XH2          => fraction of H2
  real(dp),intent(in)::T,N,x,mu,XH2,aexp
  real(dp)::coolph,coolphoto,drefdt
  real(dp),intent(out)::cool
  real(dp),parameter::scaleup=1d30

  cool = 0d0
  coolph = 0d0
  coolphoto = 0d0

  ! Get equilibrium metal cooling
  if (T < 10035.d0) then
     ! Patrick's low-temperature cooling function
     !call cooling_low(T,N,coolph)
     !note dRefDT not used here
     call hot_cold_2(T,N,coolph,dRefDT,XH2)    

     cool = -coolph
  else
     call cmp_metals(T/mu,N,mu,cool,dRefdT,aexp)
     ! Leave this out for now until thinking about it. Thoughts:
     ! 1) The is in CIE, not NEQ (see Sutherland & Dopita, 1993)
     ! 2) This seems to include hydrogen, which we already have in RAMSES-RT
     ! 3) Pretty close agreement > 10^6 with rt_cmp_cooling anyway
     !call hot_cold_2(T,N,coolph,dRefdT)
     !cool = -coolph
  end if
  ! Handle photoionisation in range where it matters (based on Ferland 2003)
  ! Add a threshold in x to make sure neutral gas is actually treated properly
  ! This is because sometimes the multiplier truncates cooling for low values
  ! Sam Geens's advice
  if ((T .lt. 1d5).and.(T .gt. 5000d0) .and.(x .gt.1d-2) .and. (N .lt. 1.d5) ) then
     ! Prevent floating point underflows by scaling up
     cool = cool*scaleup
     call cool_ferlandlike(T/mu,N,coolphoto)
     ! If the cooling is above this value just use this anyway
     if (coolphoto*scaleup .gt.cool) then
        cool = cool*(1d0-x) + coolphoto*x*scaleup
     endif
     ! Scale back down again to the physical value
     cool = cool/scaleup
  endif

END SUBROUTINE rt_metal_cool_mashup

!=========================================================================
subroutine hot_cold_2(T,n,ref,dRefDT,XH2)    
  use amr_parameters,only:dp,rt,p_UV
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
     ! TC: other waz to get XH2?
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

!=========================================================================
SUBROUTINE cool_ferlandlike(T,N,cool)
  use amr_parameters,only:dp
  implicit none
  ! Linear fit to Ferland 2003
  ! Similar to Osterbrock may with a floor below 10^3 K (assume no photoequilibrium < 1000 K...)
  ! Modified to meet our neq_chem metal cooling peak in rt_cmp_metals
  ! Compute cooling enhancement due to metals
  ! T            => Temperature in Kelvin *with mu included*
  ! N            => Hydrogen number density (H/cc)
  ! x            => Hydrogen ionisation fraction (0->1)
  ! cool         <=  Metal cooling [erg s-1 cm-3]
  real(dp),intent(in)::T,N
  real(dp),intent(out)::cool
  ! First piece: flat cooling @ 3d-24
  real(dp),parameter::cool0=3d-24
  real(dp),parameter::T0=9000d0
  ! Second piece: linear fit to meet rt_cmp_metals @ 1e5
  real(dp),parameter::cool1=2.2d-22
  real(dp),parameter::T1=1d5
  if (T < T0) then
     cool = cool0
  else
     cool = (log10(T) - log10(T0)) * (log10(cool1)-log10(cool0)) / &
          & (log10(T1)-log10(T0)) + log10(cool0)
     cool = 10d0**cool
  end if
  cool = cool*N*N ! N*Ne (fully ionised)

END SUBROUTINE cool_ferlandlike
!=========================================================================

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
