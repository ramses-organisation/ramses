module hydro_parameters

#ifdef grackle
  use grackle_parameters
#endif
  use amr_parameters
  use collapse_parameters

  ! Number of independant variables

  ! Euler variables: density, velocity, pressure
#ifdef SOLVERhydro
  integer,parameter::neul=ndim+2
#else
  ! SOLVERmhd and SOLVERrhd: always 3 velocities when MHD
  integer,parameter::neul=5
#endif

#ifndef NHYDRO
#ifdef SOLVERmhd
  ! additional magnetic field (B left)
  integer,parameter::nhydro=neul+3 !8
#else
  integer,parameter::nhydro=neul
#endif
#else
  integer,parameter::nhydro=NHYDRO
#endif
  ! non-thermal energies
#ifndef NENER
  integer,parameter::nener=0
#else
  integer,parameter::nener=NENER
#endif
  integer,parameter::inener=nhydro+1
  ! total amount of variables
#ifndef NVAR
  integer,parameter::nvar=nhydro+nener
#else
  integer,parameter::nvar=NVAR
#endif
#ifdef SOLVERmhd
  ! additional vars to store magnetic field on the right cell face
  integer,parameter::nvar_all=nvar+3
#else
  integer,parameter::nvar_all=nvar
#endif

  ! Size of hydro kernel
  integer,parameter::iu1=-1
  integer,parameter::iu2=+4
  integer,parameter::ju1=(1-ndim/2)-1*(ndim/2)
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
  real(dp),dimension(1:MAXBOUND,1:nvar_all)::boundary_var
  real(dp),dimension(1:MAXBOUND)::d_bound=0
  real(dp),dimension(1:MAXBOUND)::p_bound=0
  real(dp),dimension(1:MAXBOUND)::u_bound=0
  real(dp),dimension(1:MAXBOUND)::v_bound=0
  real(dp),dimension(1:MAXBOUND)::w_bound=0
#ifdef SOLVERmhd
  real(dp),dimension(1:MAXBOUND)::A_bound=0
  real(dp),dimension(1:MAXBOUND)::B_bound=0
  real(dp),dimension(1:MAXBOUND)::C_bound=0
#endif
  ! TODO allow other variables in inflow:
#if NENER>0
  real(dp),dimension(1:MAXBOUND,1:NENER)::prad_bound=0
#endif
#if NVAR>NHYDRO+NENER
  real(dp),dimension(1:MAXBOUND,1:NVAR-NHYDRO-NENER)::var_bound=0
#endif
  ! Refinement parameters for hydro
  real(dp)::err_grad_d=-1.0d0  ! Density gradient
  real(dp)::err_grad_u=-1.0d0  ! Velocity gradient
  real(dp)::err_grad_p=-1.0d0  ! Pressure gradient
  real(dp)::floor_d=1d-10     ! Density floor
  real(dp)::floor_u=1d-10     ! Velocity floor
  real(dp)::floor_p=1d-10     ! Pressure floor
#ifdef SOLVERmhd
  real(dp)::err_grad_A=-1.0d0  ! Bx gradient
  real(dp)::err_grad_B=-1.0d0  ! By gradient
  real(dp)::err_grad_C=-1.0d0  ! Bz gradient
  real(dp)::err_grad_B2=-1.0d0 ! B L2 norm gradient
  real(dp)::floor_A=1d-10     ! Bx floor
  real(dp)::floor_B=1d-10     ! By floor
  real(dp)::floor_C=1d-10     ! Bz floor
  real(dp)::floor_b2=1d-10    ! B L2 norm floor
#endif
#ifdef SOLVERrhd
  real(dp)::err_grad_lor=-1.0  ! Lorentz factor gradient
#endif
  real(dp)::mass_sph=0.0d0     ! mass_sph
  ! TODO allow for discontinuity-based refine on non-standard hydro vars:
#if NENER>0
  real(dp),dimension(1:NENER)::err_grad_prad=-1
#endif
#if NVAR>NHYDRO+NENER
  real(dp),dimension(1:NVAR-NHYDRO-NENER)::err_grad_var=-1
#endif
  real(dp),dimension(1:MAXLEVEL)::jeans_refine=-1

  ! Initial conditions hydro variables
  real(dp),dimension(1:MAXREGION)::d_region=0
  real(dp),dimension(1:MAXREGION)::u_region=0
  real(dp),dimension(1:MAXREGION)::v_region=0
  real(dp),dimension(1:MAXREGION)::w_region=0
  real(dp),dimension(1:MAXREGION)::p_region=0
#ifdef SOLVERmhd
  real(dp),dimension(1:MAXREGION)::A_region=0
  real(dp),dimension(1:MAXREGION)::B_region=0
  real(dp),dimension(1:MAXREGION)::C_region=0
#endif
#if NENER>0
  real(dp),dimension(1:MAXREGION,1:NENER)::prad_region=0
#endif
#if NVAR>NHYDRO+NENER
  real(dp),dimension(1:MAXREGION,1:NVAR-NHYDRO-NENER)::var_region=0
#endif
  ! Hydro solver parameters
  integer ::niter_riemann=10
  integer ::slope_type=1
  real(dp)::gamma=1.4d0
  real(dp),dimension(1:512)::gamma_rad=1.33333333334d0
  real(dp)::courant_factor=0.5d0
  real(dp)::difmag=0
  real(dp)::smallc=1.0d-10
  real(dp)::smallr=1.0d-10
  character(LEN=10)::scheme='muscl'
  character(LEN=10)::riemann='llf'
#ifdef SOLVERmhd
  integer ::slope_mag_type=-1
  real(dp)::eta_mag=0
  character(LEN=10)::riemann2d='llf'
  logical ::allow_switch_solver=.false.   ! enable on the fly switching 1D riemann solver hll or hlld to llf to prevent numerical crash
  logical ::allow_switch_solver2D=.false. ! switching for 2D riemann solver hlld to llf (checks only minimum density, needed in cosmology)
  real(dp)::switch_solv_B=1d20            ! value of B_tot**2/P above which to switch solver
  real(dp)::switch_solv_dens=1d20         ! switch solver when density discontinuity is larger than this factor
  real(dp)::switch_solv_min_dens=1d-20    ! switch solver when density is smaller than this value [c.u.]
  integer ::ischeme=0
  integer ::iriemann=0
  integer ::iriemann2d=0
#endif
#ifdef SOLVERrhd
  character(LEN=10)::eos_rhd='constant'
#endif

  ! Interpolation parameters
  integer ::interpol_var=0
  integer ::interpol_type=1
#ifdef SOLVERmhd
  integer ::interpol_mag_type=-1
#endif

  ! Passive variables index indices for star formation recipes etc, in order of access in uold/unew
  ! (updated in read_hydro_params)
  integer,parameter::imetal=inener+nener
  integer::idelay=imetal
  integer::ivirial1=imetal
  integer::ivirial2=imetal
  integer::ixion=imetal
  integer::ichem=imetal

end module hydro_parameters
