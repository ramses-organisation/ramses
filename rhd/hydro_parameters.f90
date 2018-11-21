module hydro_parameters
  use amr_parameters

  ! Number of independant variables
#ifndef NENER
  integer,parameter::nener=0
#else
  integer,parameter::nener=NENER
#endif
#ifndef NVAR
  integer,parameter::nvar=5
#else
  integer,parameter::nvar=NVAR
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
  real(dp),dimension(1:MAXBOUND,1:nvar)::boundary_var
  real(dp),dimension(1:MAXBOUND)::d_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::p_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::u_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::v_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::w_bound=0.0d0

  ! Refinement parameters for hydro
  real(dp)::err_grad_d=-1.0  ! Density gradient
  real(dp)::err_grad_u=-1.0  ! Velocity gradient
  real(dp)::err_grad_p=-1.0  ! Pressure gradient
  real(dp)::err_grad_lor=-1.0  ! Lorentz factor gradient
  real(dp)::floor_d=1d-10   ! Density floor
  real(dp)::floor_u=1d-10   ! Velocity floor
  real(dp)::floor_p=1d-10   ! Pressure floor
  real(dp)::mass_sph=0.0D0   ! mass_sph
  real(dp),dimension(1:MAXLEVEL)::jeans_refine=-1.0

  ! Initial conditions hydro variables
  real(dp),dimension(1:MAXREGION)::d_region=0.
  real(dp),dimension(1:MAXREGION)::u_region=0.
  real(dp),dimension(1:MAXREGION)::v_region=0.
  real(dp),dimension(1:MAXREGION)::w_region=0.
  real(dp),dimension(1:MAXREGION)::p_region=0.

  ! Hydro solver parameters
  integer ::niter_riemann=10
  integer ::slope_type=1
  real(dp)::gamma=1.4d0
  character(LEN=10)::eos='constant'
  real(dp)::courant_factor=0.5d0
  real(dp)::difmag=0.0d0
  real(dp)::smallc=1d-10
  real(dp)::smallr=1d-10
  character(LEN=10)::scheme='muscl'
  character(LEN=10)::riemann='llf'

  ! Interpolation parameters
  integer ::interpol_var=0
  integer ::interpol_type=1

  ! Passive variables index
  integer::imetal=6
  integer::idelay=6
  integer::ixion=6
  integer::ichem=6
  integer::ivirial1=6
  integer::ivirial2=6
  integer::inener=6

end module hydro_parameters
