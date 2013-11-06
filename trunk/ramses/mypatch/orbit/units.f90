subroutine units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  use amr_commons
  use hydro_commons
  use cooling_module
  real(dp)::scale_nH,scale_T2,scale_t,scale_v,scale_d,scale_l
  !-----------------------------------------------------------------------
  ! Conversion factors from user units into cgs units
  ! For gravity runs, make sure that G=1 in user units.
  !-----------------------------------------------------------------------
  real(kind=8)::unitscheck

  ! scale_d converts mass density from user units into g/cc
  scale_d = 3.0096099334241354d-14
  if(cosmo) scale_d = omega_m * rhoc *(h0/100.)**2 / aexp**3

  ! scale_t converts time from user units into seconds
  scale_t = 2.231293696821217d10
  if(cosmo) scale_t = aexp**2 / (h0*1d5/3.08d24)

  ! scale_l converts distance from user units into cm
  scale_l = 4.787136d15
  if(cosmo) scale_l = aexp * boxlen_ini * 3.08d24 / (h0/100)

  ! scale_v convert velocity in user units into cm/s
  scale_v = scale_l / scale_t

  ! scale_T2 converts (P/rho) in user unit into (T/mu) in Kelvin
  scale_T2 = mH/kB * scale_v**2

  ! scale_nH converts rho in user units into nH in H/cc
  scale_nH = 1./mH * scale_d

  !check for G=1
  unitscheck=scale_t**2.*scale_d*6.6738480d-8-1.
  if (abs(unitscheck)>1.d-5)then 
     print*,'units not ok for gravity run: ',unitscheck
     call clean_stop
  endif

end subroutine units
