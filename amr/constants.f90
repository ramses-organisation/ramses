module constants
  use amr_commons, ONLY: dp

  ! Numerical constants
  real(dp),parameter ::twopi        = 6.2831853d0
  real(dp),parameter ::pi           = twopi/2d0

  real(dp),parameter ::mu_mol       = 1.2195d0

  ! Physical constants
  ! Source:
  ! * Wiki - Wikipedia
  ! * IAU - Internatonal Astro Union resolution/recommendation
  ! * PCAD - http://www.astro.wisc.edu/~dolan/constants.html
  real(dp),parameter ::hplanck      = 6.6260702d-27 ! Planck const. [erg s]; arXiv:1401.8160
  real(dp),parameter ::eV2erg       = 1.6021772d-12 ! Electronvolt [erg]; PCAD
  real(dp),parameter ::kB           = 1.3806580d-16 ! Boltzm.const. [erg K-1]; PCAD
  real(dp),parameter ::rhoc         = 1.8800000d-29 ! Crit. density [g cm-3]
  real(dp),parameter ::mH           = 1.6605400d-24 ! H atom mass [g] = amu, i.e. atomic mass unit; PCAD
  real(dp),parameter ::factG_in_cgs = 6.6740831d-08 ! Grav.const. [cm3 kg-1 s-2]; Wiki
  real(dp),parameter ::a_r          = 7.5657000d-15 ! Rad.const. [erg cm-3 K-4]; Wiki
  real(dp),parameter ::M_sun        = 1.98910000d33 ! Solar Mass [g]; IAU
  real(dp),parameter ::sigma_T      = 6.6524587d-25 ! Thomson [cm2]; Wiki
  real(dp),parameter ::c_cgs        = 2.9979246d+10 ! Speed of light [cm s-1]; Wiki
  real(dp),parameter ::L_sun        = 3.8280000d+33 ! Solar Lum [erg s-1]; IAU

  ! Conversion factors - distance
  real(dp),parameter ::pc2cm        = 3.0856775d+18
  real(dp),parameter ::kpc2cm       = 3.0856775d+21
  real(dp),parameter ::Mpc2cm       = 3.0856775d+24
  real(dp),parameter ::Gpc2cm       = 3.0856775d+27

  ! Conversion factors - time
  ! Year definition follows IAU recommendation
  ! https://www.iau.org/publications/proceedings_rules/units/
  real(dp),parameter ::yr2sec       = 3.15576000d+07 ! Year [s]
  real(dp),parameter ::kyr2sec      = 3.15576000d+10 ! Kyr [s]
  real(dp),parameter ::Myr2sec      = 3.15576000d+13 ! Myr [s]
  real(dp),parameter ::Gyr2sec      = 3.15576000d+16 ! Gyr [s]
  real(dp),parameter ::sec2Gyr      = 3.16880878d-17 ! sec [Gyr]



end module constants
