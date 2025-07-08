! molecules_module.f90
module molecules_module
  use amr_parameters, only: dp
  implicit none

  private  ! everything is private by default
  public :: alpha_H2, beta_H2_umist, beta_H2, beta_H2_krome, alpha_CO, beta_CO, alpha_H2_prim, alpha_H2_dust
  public :: comp_SH2, comp_Sd

CONTAINS

FUNCTION alpha_H2_prim(T, xe, H2_cosmic_ray_ionization_rate, G0, xHI, xHII) result(rate)
  implicit none

  real(dp), intent(in) :: T, xe, H2_cosmic_ray_ionization_rate, G0
  real(dp), intent(in) :: xHI, xHII
  real(dp) :: rate
  real(dp) :: logT, lnTe
  real(dp) :: k1, k2, k5, k13, k14, k15, k_hm_cr, k_hm_gamma

  ! H- channel for H2 formation
  rate = 0.d0

  ! Primordial channel
  logT = log10(T)
  lnTe = log(T*8.621738d-5) ! K -> eV

  ! Creation and destruction channels of H- included with updated rates from Glover et al. 2010
  ! H + e- -> H- + gamma
  k1 = (10.d0**(-17.845d0 + 0.762d0*logT + 0.1523d0*(logT**2.d0) - 0.03274d0*(logT**3.d0)))
  if (T .ge. 6000.d0) then
     k1 = 10.d0**(-16.42d0 + 0.1998d0*(logT**2.d0) - 5.447d-3*(logT**4.d0) + 4.0415d-5*(logT**6.d0))
  end if

  ! H- + H -> H2 + e
  k2 = 4.0d-9*(max(T,300.d0)**(-0.17d0))

  ! H- + H+ -> H + H
  k5 = 2.4d-6/sqrt(T)*(1.d0 + T/20000.d0)

  ! H- + CR --> H + e-
  k_hm_cr = 1.28d-13 * (H2_cosmic_ray_ionization_rate / 1.d-16)

  ! H- + gamma --> H + e-
  k_hm_gamma = 5.9d-9 * G0

  ! H- + e -> H + e + e
  k13 = -1.801849334d1 + 2.36085220d0*lnTe - 2.82744300d-1*(lnTe**2.d0) &
                 +1.62331664d-2*(lnTe**3.d0)-3.36501203d-2*(lnTe**4.d0)+1.17832978d-2*(lnTe**5.) &
                 -1.65619470d-3*(lnTe**6.d0)+1.06827520d-4*(lnTe**7.d0)-2.63128581d-6*(lnTe**8.)
  k13 = exp(max(k13,-92.d0)) ! Max needed to prevent result from diverging at low temperature

  ! H- + H --> H + H + e-
  ! I think this reaction was broken in glover so I took the results from
  ! https://www.aanda.org/articles/aa/pdf/2016/02/aa27262-15.pdf Table A1
  ! Harley added the fudge factor for continuity
  k14 = 1.357772745525155d0 * 2.5634d-15 * (exp(lnTe)**1.78186d0) ! Note that T must be in eV for this reaction to make sense
  if (T .gt. 1160.d0) then
     k14 = -3.388464953d1 + 1.13944933d0*lnTe - 1.4210135d-1*(lnTe**2.d0) &
            + 8.4644554d-3*(lnTe**3.d0) - 1.4328641d-3*(lnTe**4.d0) + 2.0122503d-4*(lnTe**5.d0) &
            + 8.6639632d-5*(lnTe**6.d0) - 2.5850097d-5*(lnTe**7.d0) + 2.4555012d-6*(lnTe**8.d0) &
            - 8.0683825d-8*(lnTe**9.d0)
     k14 = exp(k14)
  end if

  ! H- + H+ --> H2+ + e- --> H + H (via recombinative dissociation)
  k15 = 6.9d-9 * (T**(-0.35d0))
  if (T .gt. 8000.d0) then
     k15 = 9.6d-7 * (T**(-0.9d0))
  end if

  rate = rate + k1*k2*xe/(k2 + k5*xHII + k_hm_cr + k_hm_gamma + k13*xe + k14*xHI + k15*xHII)

END FUNCTION alpha_H2_prim

FUNCTION alpha_H2_dust(T, dust_to_gas_mass_ratio_over_mw) result(rate)
  ! Formation on dust
  implicit none

  real(dp), intent(in) :: T, dust_to_gas_mass_ratio_over_mw
  real(dp) :: rate
  real(dp) :: clumping_factor, T2

  clumping_factor = 1.d0
  T2 = T / 100.d0
  rate = dust_to_gas_mass_ratio_over_mw * (3.5d-17) * clumping_factor * sqrt(min(T2,1.d2))

END FUNCTION alpha_H2_dust

FUNCTION alpha_H2(T, dust_to_gas_mass_ratio_over_mw, xe, H2_cosmic_ray_ionization_rate, G0, xHI, xHII, nH) result(rate)
  ! Creation rate of molecular hydrogen
  ! We consider both the primordial channel (via H-) as well
  ! as formation on dust
  implicit none
  real(dp), intent(in) :: T, dust_to_gas_mass_ratio_over_mw, xe
  real(dp), intent(in) :: H2_cosmic_ray_ionization_rate, G0, xHI
  real(dp), intent(in) :: xHII, nH
  real(dp) :: rate

  rate = 0.d0

  ! Formation rate on dust. Consider only HI
  rate = rate + alpha_H2_dust(T, dust_to_gas_mass_ratio_over_mw) * xHI * nH

  ! Primordial H- channel
  rate = rate + alpha_H2_prim(T, xe, H2_cosmic_ray_ionization_rate, G0, xHI, xHII) * xHI * nH

END FUNCTION alpha_H2

FUNCTION beta_H2_umist(T, nH, ne, nH2) result(rate)
  ! H2 destruction from umist
  implicit none

  real(dp), intent(in) :: T, nH, ne, nH2
  real(dp) :: rate
  real(dp) :: T_loc

  rate = 0.d0

  !H2 + H2 --> H2 + H + H 
  T_loc = max(min(T,41000.d0),2803.d0)
  rate = rate + (1.00d-8 * ((T_loc/300d0)**0.0d0) * exp(-84100.d0/T_loc) * nH2)

  !H2 + e- --> H + H + e-
  T_loc = max(min(T,41000.d0),3400.d0)
  rate = rate + (3.22d-9 * ((T_loc/300d0)**0.35d0) * exp(-102000.d0/T_loc) * ne)

  !H2 + H --> H + H + H
  T_loc = max(min(T,41000.d0),1833.d0)
  rate = rate + (4.67d-7 * ((T_loc/300d0)**(-1.d0)) * exp(-55000.d0/T_loc) * nH)

END FUNCTION beta_H2_umist

FUNCTION beta_H2_krome(T, nH, ne, nH2, nHe) result(rate)
  ! From bovino 2016
  ! https://www.aanda.org/articles/aa/pdf/2016/06/aa28158-16.pdf
  implicit none
  real(dp), intent(in) :: T, nH, ne, nH2, nHe
  real(dp) :: rate

  rate = 0.d0

  ! H2 + H --> H + H + H (k18)
  rate = rate + ((6.67d-12 * sqrt(T) * exp(-1.d0 * (1.d0 + (63593.d0/T)))) * nH)

  ! H2 + H2 --> H2 + H + H (k19)
  rate = rate + ((5.996d-30 * (T**4.1881d0) * ((1.d0 + (6.761d0 * T))**(-5.6881d0)) * exp(-54657.4d0/T)) * nH2)

  ! H2 + e- --> H + H + e- (k22)
  rate = rate + (4.38d-10 * (T**0.35d0) * exp(-102000.d0/T) * ne)

  ! H2 + He --> H + H + He
  rate = rate + ((10.d0**(-27.029d0 + (3.801d0 * log10(T)) - (29487.d0/T))) * nHe)

END FUNCTION beta_H2_krome

FUNCTION beta_H2(T, nH, xHI, xH2, xHe, ne, nHI, nH2, nHeI) result(rate)
  ! Returns the collisional dissociation rates of H2 for four different
  ! reactions [cm3s-1] from Glover & Abel (2008)
  ! http://mnras.oxfordjournals.org/content/388/4/1627.full.pdf 
  implicit none

  real(dp), intent(in):: T, nH, xHI, xH2, xHe
  real(dp), intent(in):: ne, nHI, nH2, nHeI
  real(dp):: rate
  real(dp):: T4, ncrH, ncrH2, ncrHe, invncr
  real(dp):: LTEfac, NLTEfac
  real(dp):: k8, k9, k9L, k10, k10L, k11, k11L
  real(dp):: lk8, lk9, lk10, lk11

  T4 = T / 1d4

  ! Critical number densities.  See eqns 15, 16, and 17
  ncrH =  10.d0**(3.d0 - 0.416d0*log10(T4) - 0.327d0*log10(T4)*log10(T4))
  ncrH2 = 10.d0**(4.845d0 - 1.3d0*log10(T4) + 1.62d0*log10(T4)*log10(T4))
  ncrHe = 10.d0**(5.0792d0*(1.d0 - 1.23d-5*(T - 2000.d0)))

  ! 1/ncr.  see eqn 14
  invncr = (xHI/ncrH) + (xH2/ncrH2) + (xHe/ncrHe)

  ! prefactors for LTE and NLTE collision rates  see eqn 13
  LTEfac = (nH*invncr)/(1.d0 + (nH*invncr))
  NLTEfac = 1.d0/(1.d0 + (nH*invncr))
  if (ne .gt. nHI) then
     LTEfac = 1.d0
     NLTEfac = 0.d0
  end if

  !reaction rates from the appendix:
  !H2 + e- --> H + H + e-
  k8 = 3.73d-9 * (T**0.1121d0) * exp(-99430.d0/T) ! Glover et al. (2010)

  !H2 + H --> H + H + H
  k9 = (6.67d-12)*sqrt(T)*exp(-1.d0*(1.d0 + (63593.d0/T)))
  k9L = (3.52d-9)*exp(-43900.d0/T)

  !H2 + H2 --> H2 + H + H
  k10 = ( (5.996d-30*(T**4.1881d0)) / ((1.d0 + 6.761d-6*T)**5.6881d0)) * exp(-54657.4d0/T)
  k10L = (1.3d-9)*exp(-53300.d0/T)

  !H2 + He --> H + H + He
  k11 = 10.d0**(-27.029d0 + (3.801d0*log10(T)) - (29487.d0/T))
  k11L = 10.d0**(-2.729d0 - (1.75d0*log10(T)) - (23474.d0/T))

  k8   = max(k8, 1d-40)
  k9   = max(k9, 1d-40)
  k10  = max(k10, 1d-40)
  k11  = max(k11, 1d-40)
  k9L  = max(k9L, 1d-40)
  k10L = max(k10L, 1d-40)
  k11L = max(k11L, 1d-40)

  !Log of all the rates
  lk8 = log10(k8)
  lk9 = (LTEfac*log10(k9L)) + (NLTEfac*log10(k9))
  lk10 = (LTEfac*log10(k10L)) + (NLTEfac*log10(k10))
  lk11 = (LTEfac*log10(k11L)) + (NLTEfac*log10(k11))

  rate = (ne*(10.d0**lk8)) + (nHI*(10.d0**lk9)) + (nH2*(10.d0**lk10)) + (nHeI*(10.0**lk11))
  rate = max(rate, 1d-40)

END FUNCTION beta_H2

FUNCTION alpha_CO(G0, xi_cr_H2, nCII, nH2, xO) result(rate)
  ! see glover 2012
  implicit none

  real(dp), intent(in):: G0, xi_cr_H2, nCII, nH2, xO
  real(dp):: rate
  real(dp):: k0, k1, gammaCHx_cr, gammaCHx, beta

  k0 = 5.d-16 ! cm^3 s^-1
  k1 = 5.d-10 ! Rate coefficient for the formation of CO from O + CHx

  ! Assuming CO formation is modulated by CH2+, https://home.strw.leidenuniv.nl/~ewine/photo/display_ch2+_65e31a07e69e64dbd64d37801983018f.html
  gammaCHx_cr = 8.88d-15 * (xi_cr_H2 / 1d-16) ! Cosmic rays
  gammaCHx = (1.41d-10 * G0) + gammaCHx_cr

  beta = k1 * xO/(k1*xO + gammaCHx/nH2)
  rate = k0 * nCII * nH2 * beta

END FUNCTION alpha_CO

FUNCTION beta_CO(G0, xi_cr_H2) result(rate)
  ! CO destruction
  implicit none

  real(dp), intent(in):: G0, xi_cr_H2
  real(dp):: rate
  real(dp):: gammaCO, gammaCO_cr

  gammaCO = 2.43d-10 * G0
  gammaCO_cr = 4.62d-15 * (xi_cr_H2 / 1d-16)

  rate = gammaCO + gammaCO_cr

END FUNCTION beta_CO

FUNCTION comp_Sd(nHI, nH2, dx_SS, Z) result(ss_factor)
  ! Returns the self shielding factor for dust
  ! see section 2.2 http://iopscience.iop.org/0004-637X/697/1/55/pdf/apj_697_1_55.pdf
  implicit none
  real(dp), intent(in)::nHI, nH2, dx_SS, Z
  real(dp):: ss_factor
  real(dp):: Sdeff, cNHI, cNH2

  Sdeff = 2.34d-21    ! dust cross section cm^2 ! Updated for bare-gr-s
  cNHI = nHI*dx_SS    ! HI column density
  cNH2 = nH2*dx_SS    ! H2 column density

  ss_factor = exp(-Sdeff*Z*(cNHI + (2.d0*cNH2)))
END FUNCTION comp_Sd

FUNCTION comp_SH2(nH2, dx_SS) result(ss_factor)
  ! Returns the self shielding factor for dust
  ! see section 2.2 http://iopscience.iop.org/0004-637X/697/1/55/pdf/apj_697_1_55.pdf
  implicit none
  real(dp), intent(in):: nH2, dx_SS
  real(dp):: ss_factor
  real(dp):: xfac, cNH2, wH2, Sa, Sb, Sc

  cNH2 = nH2*dx_SS  !H2 column density
  xfac = cNH2/(5.d14)
  wH2 = 0.2d0

  Sa = (1.d0 - wH2)/((1.d0 + xfac)*(1.d0 + xfac))
  Sb = wH2/sqrt(1.d0 + xfac)
  Sc = exp(-0.00085d0*sqrt(1.d0 + xfac))

  ss_factor = Sa + (Sb*Sc)
END FUNCTION comp_SH2

end module molecules_module
