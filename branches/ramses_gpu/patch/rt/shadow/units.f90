!************************************************************************
subroutine units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

! Conversion factors from user units into cgs units
! For gravity runs, make sure that G=1 in user units.
!
! scale_l   <=  length scale [cm]
! scale_t   <=  time scale [t]
! scale_d   <=  mass density scale [g cm-3]
! scale_v   <=  speed scale [cm s-1]
! scale_nH  <=  number density scale [cm-3]
! scale_T2  <=  Temperature scale [Kelvin]
!------------------------------------------------------------------------
  use amr_commons
  use hydro_commons
  use cooling_module
  implicit none
  real(dp)::scale_nH,scale_T2,scale_t,scale_v,scale_d,scale_l
!------------------------------------------------------------------------
  ! scale_d converts mass density from user units into g/cc
  scale_d = mH

  ! scale_t converts time from user units (Myr) into seconds
  scale_t = 3.1556926d13

  ! scale_l converts distance from user units (kpc) into cm
  ! since the box always has length 1 in UU, scale_l is the 
  ! length of the box in cgs
  scale_l = 3.08568025d21

  ! scale_v converts velocity in user units into cm/s
  scale_v = scale_l / scale_t

  ! scale_T2 converts (P/rho) in user unit into (T/mu) in Kelvin,
  ! using ideal gas equation
  scale_T2 = mH/kB * scale_v**2 !1.1080168e+08

  ! scale_nH converts rho in user units into nH in H/cc
  scale_nH = X/mH * scale_d != X

  !scale_p=1.58708e-8   !!!1.79751d-13

end subroutine units
