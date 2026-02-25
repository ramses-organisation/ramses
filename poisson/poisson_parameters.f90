module poisson_parameters
  use amr_parameters

  ! Convergence criterion for Poisson solvers
  real(dp)::epsilon=1.0D-4

  ! Type of force computation
  integer ::gravity_type=0             ! deprecated
  logical ::self_gravity=.true.        ! activate self-gravity
  integer ::gravity_rho_ana_type=0     ! include an analytical density field in Poisson source term
  integer ::gravity_force_ana_type=0   ! add analytical gravitational force

  ! Gravity parameters
  real(dp),dimension(1:10)::gravity_params=0           ! deprecated
  real(dp),dimension(1:10)::gravity_rho_ana_params=0   ! parameters for analytical density field
  real(dp),dimension(1:10)::gravity_force_ana_params=0  ! parameters for analytical gravity force

  ! Maximum level for CIC dark matter interpolation
  integer :: cic_levelmax=0

  ! Min level for CG solver
  ! level < cg_levelmin uses fine multigrid
  ! level >=cg_levelmin uses conjugate gradient
  integer :: cg_levelmin=999

  ! Gauss-Seidel smoothing sweeps for fine multigrid
  integer, parameter :: ngs_fine   = 2
  integer, parameter :: ngs_coarse = 2

  ! Number of multigrid cycles for coarse levels *in safe mode*
  !   1 is the fastest,
  !   2 is slower but can give much better convergence in some cases
  integer, parameter :: ncycles_coarse_safe = 1

end module poisson_parameters
