subroutine units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  use amr_commons
  use hydro_commons
  use constants, only: Mpc2cm, mH, kB, rhoc
  use cooling_module, only: X
#if USE_FLD==1
  use units_commons, only : scale_kappa
#endif
  implicit none

  real(dp)::scale_nH,scale_T2,scale_t,scale_v,scale_d,scale_l
  real(dp),parameter::Grav=6.67e-08_dp  
  real(dp) :: pc
  
  pc = 3.08e18_dp

  !-----------------------------------------------------------------------
  ! Conversion factors from user units into cgs units
  ! For gravity runs, make sure that G=1 in user units.
  !-----------------------------------------------------------------------

  ! scale_d converts mass density from user units into g/cc
  scale_d = mu_gas*mH

  ! scale_t converts time from user units into seconds
  !scale_t = units_time
  !if(cosmo) scale_t = aexp**2 / (h0*1d5/Mpc2cm)
  ! scale_t converts time from user units into seconds
  scale_t = 1.0_dp/sqrt(Grav*scale_d)

  ! scale_l converts distance from user units into cm
  !scale_l = units_length
  !if(cosmo) scale_l = aexp * boxlen_ini * Mpc2cm / (h0/100)
  !calculate the initial cloud radius
  scale_l = pc

  ! scale_v converts velocity in user units into cm/s
  scale_v = scale_l / scale_t


  ! scale_T2 converts (P/rho) in user unit into (T/mu) in Kelvin
  scale_T2 = mH/kB * scale_v**2

  !scale_T2 = mu_gas**2 * mH**2 * pc**2 * Grav / kb

  ! scale_nH converts rho in user units into nH in H/cc
  scale_nH = 1.0d0

end subroutine units
