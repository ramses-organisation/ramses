module hydro_parameters
#ifdef grackle
  use grackle_parameters
#endif
  use amr_parameters

  ! Number of independant variables
#ifndef NENER
  integer,parameter::nener=0
#else
  integer,parameter::nener=NENER
#endif
#ifndef NVAR
  integer,parameter::nvar=8+nener
#else
  integer,parameter::nvar=NVAR
#endif

  ! Size of hydro kernel
  integer,parameter::iu1=-1
  integer,parameter::iu2=+4
  integer,parameter::ju1=(1-ndim/2)-1*(ndim/2) !integer division
  integer,parameter::ju2=(1-ndim/2)+4*(ndim/2)
  integer,parameter::ku1=(1-ndim/3)-1*(ndim/3)
  integer,parameter::ku2=(1-ndim/3)+4*(ndim/3)
  integer,parameter::if1=1
  integer,parameter::if2=3
  integer,parameter::jf1=1
  integer,parameter::jf2=(1-ndim/2)+3*(ndim/2)
  integer,parameter::kf1=1
  integer,parameter::kf2=(1-ndim/3)+3*(ndim/3)

  ! Imposed boundary condition variables
  real(dp),dimension(1:MAXBOUND,1:nvar+3)::boundary_var
  real(dp),dimension(1:MAXBOUND)::d_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::p_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::u_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::v_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::w_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::A_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::B_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::C_bound=0.0d0
#if NENER>0
  real(dp),dimension(1:MAXBOUND,1:NENER)::prad_bound=0.0
#endif
#if NVAR>8+NENER
  real(dp),dimension(1:MAXBOUND,1:NVAR-8-NENER)::var_bound=0.0
#endif

  ! Refinement parameters for hydro
  real(dp)::err_grad_d=-1.0  ! Density gradient
  real(dp)::err_grad_u=-1.0  ! Velocity gradient
  real(dp)::err_grad_p=-1.0  ! Pressure gradient
  real(dp)::err_grad_A=-1.0  ! Bx gradient
  real(dp)::err_grad_B=-1.0  ! By gradient
  real(dp)::err_grad_C=-1.0  ! Bz gradient
  real(dp)::err_grad_B2=-1.0 ! B L2 norm gradient
  real(dp)::floor_d=1d-10   ! Density floor
  real(dp)::floor_u=1d-10   ! Velocity floor
  real(dp)::floor_p=1d-10   ! Pressure floor
  real(dp)::floor_A=1d-10   ! Bx floor
  real(dp)::floor_B=1d-10   ! By floor
  real(dp)::floor_C=1d-10   ! Bz floor
  real(dp)::floor_b2=1d-10  ! B L2 norm floor
  real(dp)::mass_sph=0.0D0   ! mass_sph
#if NENER>0
  real(dp),dimension(1:NENER)::err_grad_prad=-1.0
#endif
#if NVAR>8+NENER
  real(dp),dimension(1:NVAR-8-NENER)::err_grad_var=-1.0
#endif
  real(dp),dimension(1:MAXLEVEL)::jeans_refine=-1.0

  ! Initial conditions hydro variables
  real(dp),dimension(1:MAXREGION)::d_region=0.
  real(dp),dimension(1:MAXREGION)::u_region=0.
  real(dp),dimension(1:MAXREGION)::v_region=0.
  real(dp),dimension(1:MAXREGION)::w_region=0.
  real(dp),dimension(1:MAXREGION)::p_region=0.
  real(dp),dimension(1:MAXREGION)::A_region=0.
  real(dp),dimension(1:MAXREGION)::B_region=0.
  real(dp),dimension(1:MAXREGION)::C_region=0.
#if NENER>0
  real(dp),dimension(1:MAXREGION,1:NENER)::prad_region=0.0
#endif
#if NVAR>8+NENER
  real(dp),dimension(1:MAXREGION,1:NVAR-8-NENER)::var_region=0.0
#endif

  ! Hydro solver parameters
  integer ::niter_riemann=10
  integer ::slope_type=1
  integer ::slope_mag_type=-1
  real(dp)::slope_theta=1.5d0
  real(dp)::gamma=1.4d0
  real(dp),dimension(1:512)::gamma_rad=1.33333333334d0
  real(dp)::courant_factor=0.5d0
  real(dp)::difmag=0.0d0
  real(dp)::smallc=1d-10
  real(dp)::smallr=1d-10
  real(dp)::eta_mag=0.0d0
  character(LEN=10)::scheme='muscl'
  character(LEN=10)::riemann='llf'
  character(LEN=10)::riemann2d='llf'
  logical ::allow_switch_solver=.false.   ! enable on the fly switching 1D riemann solver hll or hlld to llf to prevent numerical crash
  logical ::allow_switch_solver2D=.false. ! switching for 2D riemann solver hlld to llf (checks only minimum density, needed in cosmology)
  real(dp)::switch_solv_B=1d20            ! value of B_tot**2/P above which to switch solver
  real(dp)::switch_solv_dens=1d20         ! switch solver when density discontinuity is larger than this factor
  real(dp)::switch_solv_min_dens=1d-20    ! switch solver when density is smaller than this value [c.u.]
  integer ::ischeme=0
  integer ::iriemann=0
  integer ::iriemann2d=0

  ! Interpolation parameters
  integer ::interpol_var=0
  integer ::interpol_type=1
  integer ::interpol_mag_type=-1

  ! Passive variables index
  integer::imetal=9
  integer::idelay=9
  integer::ixion=9
  integer::ichem=9
  integer::ivirial1=9
  integer::ivirial2=9
  integer::inener=9

end module hydro_parameters
