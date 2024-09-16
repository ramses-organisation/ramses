module hydro_parameters

#ifdef grackle
  use grackle_parameters
#endif
  use amr_parameters

  ! Number of independant variables
#ifndef NMAT
  integer,parameter::nmat=2
#else
  integer,parameter::nmat=NMAT
#endif
  integer,parameter::nener=0
#ifndef NVAR
  integer,parameter::nvar=ndim+3*nmat
  integer,parameter::npri=ndim+2*nmat
#else
  integer,parameter::nvar=NVAR
  integer,parameter::npri=NVAR-nmat
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
  real(dp),dimension(1:MAXBOUND)::d1_bound=0,d2_bound=0,d3_bound=0,d4_bound=0
  real(dp),dimension(1:MAXBOUND)::f1_bound=0,f2_bound=0,f3_bound=0,f4_bound=0
  real(dp),dimension(1:MAXBOUND)::p1_bound=0,p2_bound=0,p3_bound=0,p4_bound=0
  real(dp),dimension(1:MAXBOUND)::u_bound=0,v_bound=0,w_bound=0

  ! Refinement parameters for hydro
  real(dp)::err_grad_f=-1.0d0    ! Volume fraction
  real(dp)::err_grad_d=-1.0d0    ! Density gradient
  real(dp)::err_grad_u=-1.0d0    ! Velocity gradient
  real(dp)::err_grad_p=-1.0d0    ! Pressure gradient
  real(dp)::floor_f=1d-10     ! Fraction floor
  real(dp)::floor_d=1d-10     ! Density floor
  real(dp)::floor_u=1d-10     ! Velocity floor
  real(dp)::floor_p=1d-10     ! Pressure floor
  real(dp)::mass_sph=0.0d0     ! mass_sph

  ! Initial conditions hydro variables
  real(dp),dimension(1:MAXREGION)::d1_region=0,d2_region=0,d3_region=0,d4_region=0
  real(dp),dimension(1:MAXREGION)::f1_region=0,f2_region=0,f3_region=0,f4_region=0
  real(dp),dimension(1:MAXREGION)::p1_region=0,p2_region=0,p3_region=0,p4_region=0
  real(dp),dimension(1:MAXREGION)::u_region=0.,v_region=0.,w_region=0.
  real(dp),dimension(1:MAXREGION)::p_region=0.

  ! Multi-material properties
  integer ,dimension(1:100)::eos_type=0 ! Perfect gas
  real(dp),dimension(1:100,1:10)::eos_params=1.4
  character(LEN=128),dimension(1:100)::eos_file=' '

  ! Hydro solver parameters
  integer ::geom=1
  integer ::slope_type=1
  real(dp)::gamma=-999d0
  real(dp)::courant_factor=0.5d0
  real(dp)::difmag=0
  real(dp)::smallc=1.0d-10
  real(dp)::smallr=1.0d-10
  real(dp)::smallf=1.0d-10
  character(LEN=10)::scheme='muscl'
  character(LEN=20)::eos_name='mie-grueneisen'

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
