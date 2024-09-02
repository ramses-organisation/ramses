module hydro_parameters

#ifdef grackle
    use grackle_parameters
#endif
    use amr_parameters

    ! Number of independant variables
#ifndef NENER
    integer, parameter :: nener=0
#else
    integer, parameter :: nener=NENER
#endif
#ifndef NVAR
    integer, parameter :: nvar=ndim + 2 + nener
#else
    integer, parameter :: nvar=NVAR
#endif
    ! Size of hydro kernel
    integer, parameter :: iu1=- 1
    integer, parameter :: iu2=+ 4
    integer, parameter :: ju1=(1 - ndim / 2) - 1 * (ndim / 2)
    integer, parameter :: ju2=(1 - ndim / 2) + 4 * (ndim / 2)
    integer, parameter :: ku1=(1 - ndim / 3) - 1 * (ndim / 3)
    integer, parameter :: ku2=(1 - ndim / 3) + 4 * (ndim / 3)
    integer, parameter :: if1=1
    integer, parameter :: if2=3
    integer, parameter :: jf1=1
    integer, parameter :: jf2=(1 - ndim / 2) + 3 * (ndim / 2)
    integer, parameter :: kf1=1
    integer, parameter :: kf2=(1 - ndim / 3) + 3 * (ndim / 3)

    ! Imposed boundary condition variables
    real(dp), dimension(1:MAXBOUND, 1:nvar) :: boundary_var
    real(dp), dimension(1:MAXBOUND) :: d_bound=0
    real(dp), dimension(1:MAXBOUND) :: p_bound=0
    real(dp), dimension(1:MAXBOUND) :: u_bound=0
    real(dp), dimension(1:MAXBOUND) :: v_bound=0
    real(dp), dimension(1:MAXBOUND) :: w_bound=0
#if NENER > 0
    real(dp), dimension(1:MAXBOUND, 1:NENER) :: prad_bound=0
#endif
#if NVAR > NDIM + 2 + NENER
    real(dp), dimension(1:MAXBOUND, 1:NVAR - NDIM - 2 - NENER) :: var_bound=0
#endif
    ! Refinement parameters for hydro
    real(dp) :: err_grad_d=- 1.0d0  ! Density gradient
    real(dp) :: err_grad_u=- 1.0d0  ! Velocity gradient
    real(dp) :: err_grad_p=- 1.0d0  ! Pressure gradient
    real(dp) :: floor_d=1d-10     ! Density floor
    real(dp) :: floor_u=1d-10     ! Velocity floor
    real(dp) :: floor_p=1d-10     ! Pressure floor
    real(dp) :: mass_sph=0.0d0     ! mass_sph
#if NENER > 0
    real(dp), dimension(1:NENER) :: err_grad_prad=- 1
#endif
#if NVAR > NDIM + 2 + NENER
    real(dp), dimension(1:NVAR - NDIM - 2) :: err_grad_var=- 1
#endif
    real(dp), dimension(1:MAXLEVEL) :: jeans_refine=- 1

    ! Initial conditions hydro variables
    real(dp), dimension(1:MAXREGION) :: d_region=0
    real(dp), dimension(1:MAXREGION) :: u_region=0
    real(dp), dimension(1:MAXREGION) :: v_region=0
    real(dp), dimension(1:MAXREGION) :: w_region=0
    real(dp), dimension(1:MAXREGION) :: p_region=0
#if NENER > 0
    real(dp), dimension(1:MAXREGION, 1:NENER) :: prad_region=0
#endif
#if NVAR > NDIM + 2 + NENER
    real(dp), dimension(1:MAXREGION, 1:NVAR - NDIM - 2 - NENER) :: var_region=0
#endif
    ! Hydro solver parameters
    integer :: niter_riemann=10
    integer :: slope_type=1
    real(dp) :: slope_theta=1.5d0
    real(dp) :: gamma=1.4d0
    real(dp), dimension(1:512) :: gamma_rad=1.33333333334d0
    real(dp) :: courant_factor=0.5d0
    real(dp) :: difmag=0
    real(dp) :: smallc=1.0d-10
    real(dp) :: smallr=1.0d-10
    character(LEN=10) :: scheme='muscl'
    character(LEN=10) :: riemann='llf'

    ! Interpolation parameters
    integer :: interpol_var=0
    integer :: interpol_type=1

    ! Passive variables index
    integer :: imetal=6
    integer :: idelay=6
    integer :: ixion=6
    integer :: ichem=6
    integer :: ivirial1=6
    integer :: ivirial2=6
    integer :: inener=6

end module hydro_parameters
