!================================================================
! CGM Commons Module
! Parameters for hot rotating CGM inflow simulations
! Based on Stern et al. (2024) MNRAS 530, 1711
!================================================================
module cgm_commons
  use amr_parameters, only: dp

  implicit none

  ! Physical constants (CGS)
  real(dp), parameter :: kpc2cm = 3.0856776d21      ! kpc to cm
  real(dp), parameter :: km2cm  = 1.0d5             ! km to cm
  real(dp), parameter :: Myr2s  = 3.15576d13        ! Myr to seconds
  real(dp), parameter :: Msun   = 1.989d33          ! Solar mass in g
  real(dp), parameter :: yr2s   = 3.15576d7         ! Year to seconds

  !================================================================
  ! CGM model parameters (set via namelist)
  !================================================================

  ! Circular velocity [km/s] - sets virial temperature
  ! T_vir ~ 2 x 10^6 K for v_c = 200 km/s
  real(dp) :: cgm_v_c = 200.0d0

  ! Maximum circularization radius [kpc]
  ! Where angular momentum support becomes dominant
  real(dp) :: cgm_Rc_max = 15.0d0

  ! Central hydrogen number density [cm^-3] at r = 10 kpc
  ! From Eq. 10: n_H = 0.8e-3 * (r/10kpc)^-1.5 * ...
  real(dp) :: cgm_n_H_0 = 8.0d-4

  ! Minimum hydrogen density [cm^-3]
  real(dp) :: cgm_n_H_floor = 1.0d-6

  ! Cooling function normalization [erg cm^3 s^-1]
  ! Lambda ~ 3e-23 for T ~ 10^6 K, Z ~ 0.3 Zsun
  real(dp) :: cgm_Lambda = 3.0d-23

  ! Mean molecular weight (fully ionized, primordial)
  real(dp) :: cgm_mu = 0.6d0

  ! CGM metallicity [Zsun]
  real(dp) :: cgm_Z = 0.3d0

  ! Softening radius [cm] to avoid singularity at center
  real(dp) :: cgm_r_soft = 1.0d20  ! ~0.3 kpc

  ! Cooling radius [kpc] - outer boundary of inflow region
  real(dp) :: cgm_r_cool = 100.0d0

  ! Mass inflow rate [Msun/yr] - for reference
  real(dp) :: cgm_Mdot = 1.0d0

  !================================================================
  ! Boundary condition parameters
  !================================================================

  ! Boundary type for CGM:
  ! 0 = standard RAMSES boundary (from namelist)
  ! 1 = non-reflecting characteristic BC
  ! 2 = sponge layer (relaxation toward analytical)
  integer :: cgm_bc_type = 1

  ! Sponge layer width [code units, fraction of box]
  real(dp) :: cgm_sponge_width = 0.1d0

  ! Sponge relaxation timescale [in units of local sound crossing time]
  real(dp) :: cgm_sponge_tau = 1.0d0

  !================================================================
  ! Derived quantities (computed at runtime)
  !================================================================

  ! Virial temperature [K]
  real(dp) :: cgm_T_vir

  ! Sound speed at virial temperature [cm/s]
  real(dp) :: cgm_cs_vir

  ! t_cool/t_ff at R_c,max
  real(dp) :: cgm_tcool_tff_Rc

end module cgm_commons


!================================================================
! Initialize CGM parameters
!================================================================
subroutine init_cgm_params()
  use cgm_commons
  use constants, only: kB, mH
  implicit none

  real(dp) :: v_c_cgs, Rc_max_cm, n_H_Rc, tcool_Rc, tff_Rc

  ! Convert to CGS
  v_c_cgs = cgm_v_c * km2cm
  Rc_max_cm = cgm_Rc_max * kpc2cm

  ! Virial temperature: T = (10/9) * mu * m_H * v_c^2 / k_B
  cgm_T_vir = (10.0d0/9.0d0) * cgm_mu * mH * v_c_cgs**2 / kB

  ! Sound speed: c_s = sqrt(10/9) * v_c
  cgm_cs_vir = sqrt(10.0d0/9.0d0) * v_c_cgs

  ! Density at R_c,max
  n_H_Rc = cgm_n_H_0 * (cgm_Rc_max / 10.0d0)**(-1.5d0)

  ! Cooling time at R_c,max
  ! t_cool = (3/2) * n * k * T / (n^2 * Lambda) = (3/2) * k * T / (n * Lambda)
  tcool_Rc = 1.5d0 * kB * cgm_T_vir / (n_H_Rc * cgm_Lambda)

  ! Free-fall time at R_c,max
  tff_Rc = sqrt(2.0d0) * Rc_max_cm / v_c_cgs

  ! t_cool / t_ff ratio
  cgm_tcool_tff_Rc = tcool_Rc / tff_Rc

end subroutine init_cgm_params


!================================================================
! Read CGM parameters from namelist
!================================================================
subroutine read_cgm_params(nml_ok)
  use cgm_commons
  implicit none
  logical, intent(out) :: nml_ok

  namelist /cgm_params/ cgm_v_c, cgm_Rc_max, cgm_n_H_0, cgm_n_H_floor, &
       cgm_Lambda, cgm_mu, cgm_Z, cgm_r_soft, cgm_r_cool, cgm_Mdot, &
       cgm_bc_type, cgm_sponge_width, cgm_sponge_tau

  ! Set default values (already set in module, but be explicit)
  cgm_v_c = 200.0d0
  cgm_Rc_max = 15.0d0
  cgm_n_H_0 = 8.0d-4
  cgm_n_H_floor = 1.0d-6
  cgm_Lambda = 3.0d-23
  cgm_mu = 0.6d0
  cgm_Z = 0.3d0
  cgm_r_soft = 1.0d20
  cgm_r_cool = 100.0d0
  cgm_Mdot = 1.0d0
  cgm_bc_type = 1
  cgm_sponge_width = 0.1d0
  cgm_sponge_tau = 1.0d0

  ! Read namelist file
  rewind(1)
  read(1, NML=cgm_params, END=101, ERR=101)
  nml_ok = .true.
  goto 102

101 nml_ok = .false.

102 continue

  ! Initialize derived quantities
  call init_cgm_params()

end subroutine read_cgm_params
