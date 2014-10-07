subroutine units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  use amr_commons
  use hydro_commons
  use cooling_module
  real(dp)::scale_nH,scale_T2,scale_t,scale_v,scale_d,scale_l
  real(dp)::solar_mass, kpc, G_cgs
  !-----------------------------------------------------------------------
  ! Conversion factors from user units into cgs units
  ! For gravity runs, make sure that G=1 in user units.
  !-----------------------------------------------------------------------

  ! scale_d converts mass density from user units into g/cc
  !scale_d = 1.0
  solar_mass    = 1.9891D33 ! g
  kpc           = 3.085677581282D21 ! cm
  ! Mass unit: 10**9 solar masses
  scale_d       = 1.0D9 * solar_mass / kpc**3
  if(cosmo) scale_d = omega_m * rhoc *(h0/100.)**2 / aexp**3

  ! scale_t converts time from user units into seconds
  G_cgs         = 6.67428D-8
  scale_t       = 1.0D0 / sqrt(G_cgs * scale_d)
  if(cosmo) scale_t = aexp**2 / (h0*1d5/3.08d24)

  ! scale_l converts distance from user units into cm
  ! Length unit: kpc
  scale_l       = kpc
  if(cosmo) scale_l = aexp * boxlen_ini * 3.08d24 / (h0/100)

  ! scale_v convert velocity in user units into cm/s
  scale_v = scale_l / scale_t

  ! scale_T2 converts (P/rho) in user unit into (T/mu) in Kelvin
  scale_T2 = mH/kB * scale_v**2

  ! scale_nH converts rho in user units into nH in H/cc
  scale_nH = X/mH * scale_d

end subroutine units
