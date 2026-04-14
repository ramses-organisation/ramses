module constants
  use amr_commons, ONLY: dp

  ! Numerical constants
  real(dp),parameter ::twopi        = 6.2831853d0

  real(dp),parameter ::mu_mol       = 1.2195d0

  ! Physical constants
  ! Source:
  ! * SI - SI Brochure (2018)
  ! * PCAD - http://www.astro.wisc.edu/~dolan/constants.html
  ! * NIST - National Institute of Standards and Technology
  ! * IAU - Internatonal Astronomical Union resolution
#if USE_FLD==0
  real(dp),parameter ::pi           = twopi/2d0
  real(dp),parameter ::hplanck      = 6.6260702d-27 ! Planck const. [erg s]; SI
  real(dp),parameter ::eV2erg       = 1.6021766d-12 ! Electronvolt [erg]; SI
  real(dp),parameter ::kB           = 1.3806490d-16 ! Boltzmann const. [erg K-1]; SI
  real(dp),parameter ::c_cgs        = 2.9979246d+10 ! Speed of light [cm s-1]; SI
  real(dp),parameter ::a_r          = 7.5657233d-15 ! Radiation density const. [erg cm-3 K-4]; SI (derived)
  real(dp),parameter ::mH           = 1.6605390d-24 ! H atom mass [g] = amu, i.e. atomic mass unit; NIST
  real(dp),parameter ::M_sun        = 1.9891000d+33 ! Solar Mass [g]; IAU
  real(dp),parameter ::L_sun        = 3.8280000d+33 ! Solar Lum [erg s-1]; IAU
  real(dp),parameter ::R_sun        = 6.9570000d+10 ! Solar Radius [cm]; IAU

  ! Conversion factors - distance
  ! IAU 2012 convention:
  ! 1 pc = 648000 AU / pi
  ! 1 AU = 14 959 787 070 000 cm
  ! real(dp),parameter ::pc2cm        = 3.0856776d+18
  ! real(dp),parameter ::kpc2cm       = 3.0856776d+21
  ! real(dp),parameter ::Mpc2cm       = 3.0856776d+24
  ! real(dp),parameter ::Gpc2cm       = 3.0856776d+27
#else 
  ! take the ramses-ism values for getting automatic tests cases passed
  real(dp),parameter ::pi           = acos(-1.0d0)!twopi/2d0
  real(dp),parameter ::hplanck      = 6.6262000d-27!6.6260702d-27 ! Planck const. [erg s]; SI
  real(dp),parameter ::eV2erg       = 1.6022000d-12!1.6021766d-12 ! Electronvolt [erg]; SI
  real(dp),parameter ::kB           = 1.3806200d-16 !1.3806490d-16 ! Boltzmann const. [erg K-1]; SI
  real(dp),parameter ::c_cgs        = 2.9979250d+10!2.9979246d+10 ! Speed of light [cm s-1]; SI (TC:used)
  real(kind=8),parameter ::clight   = 2.9979250d+10 !TC: two inconsistent values are used (c_cgs)
  real(dp),parameter ::a_r          = 7.5657d-15 !7.5657233d-15 ! Radiation density const. [erg cm-3 K-4]; SI (derived)
  real(dp),parameter:: aR=7.56591469318689378e-015_dp
  real(dp),parameter ::mH           = 1.6600000d-24 !1.6605390d-24 ! H atom mass [g] = amu, i.e. atomic mass unit; NIST
  real(dp),parameter::M_sun         = 1.9889e33_dp  !Sun mass in g !TC: todo get rid of this
  real(dp),parameter::R_sun         = 6.95508e10_dp !Sun radius in cm
  real(dp),parameter::L_sun         = 3.846e33_dp   !Sun luminosity in erg/s

  real(dp),parameter ::pc2cm        = 3.08d+18
  real(dp),parameter ::kpc2cm       = 3.08d+21
  real(dp),parameter ::Mpc2cm       = 3.08d+24
  real(dp),parameter ::Gpc2cm       = 3.08d+27
#endif
  real(dp),parameter ::factG_in_cgs = 6.6740800d-08 ! Gravitational const. [cm3 g-1 s-2]; NIST
  real(dp),parameter ::sigma_T      = 6.6524587d-25 ! Thomson scattering cross-section [cm2]; NIST
  real(dp),parameter ::rhoc         = 1.8800000d-29 ! Crit. density [g cm-3]

  ! Conversion factors - distance
  ! IAU 2012 convention:
  ! 1 pc = 648000 AU / pi
  ! 1 AU = 14 959 787 070 000 cm
  ! real(dp),parameter ::pc2cm        = 3.0856776d+18
  ! real(dp),parameter ::kpc2cm       = 3.0856776d+21
  ! real(dp),parameter ::Mpc2cm       = 3.0856776d+24
  ! real(dp),parameter ::Gpc2cm       = 3.0856776d+27

  ! Conversion factors - time
  ! Year definition follows IAU recommendation
  ! https://www.iau.org/publications/proceedings_rules/units/
  real(dp),parameter ::yr2sec       = 3.15576000d+07 ! Year [s]
  real(dp),parameter ::kyr2sec      = 3.15576000d+10 ! Kyr [s]
  real(dp),parameter ::Myr2sec      = 3.15576000d+13 ! Myr [s]
  real(dp),parameter ::Gyr2sec      = 3.15576000d+16 ! Gyr [s]
  real(dp),parameter ::sec2Gyr      = 3.16880878d-17 ! sec [Gyr]



end module constants
