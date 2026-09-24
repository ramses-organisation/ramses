module poisson_commons
  use amr_commons
  use poisson_parameters

  real(dp),allocatable,dimension(:)  ::phi,phi_old       ! Potential
  real(dp),allocatable,dimension(:)  ::rho               ! Density
  real(dp),allocatable,dimension(:,:)::f                 ! 3-force

  real(dp),allocatable,dimension(:)  ::rho_top   ! Density at last CIC level

  ! Multigrid lookup table for amr -> mg index mapping
  integer, allocatable, dimension(:) :: lookup_mg   ! Lookup table

  ! Communicator arrays for multigrid levels
#ifdef LIGHT_MPI_COMM
  type(communicator_light), allocatable, dimension(:,:) :: active_mg
  type(communicator_varoct), allocatable, dimension(:) :: emission_mg
#else
  type(communicator), allocatable, dimension(:,:) :: active_mg
  type(communicator), allocatable, dimension(:,:) :: emission_mg
#endif

  ! Send/recv Multigrid temporary communicator (light) used in build_parent_comms_mg subroutine
  type communicator_mg
    integer                       ::ngrid
    integer, dimension(:), pointer::igrid
  end type communicator_mg

  ! Minimum MG level
  integer :: levelmin_mg

  ! Multigrid safety switch for the order of the multigrid boundary reconstruction scheme
  ! (see sect. 4 of Guillet & Teyssier 2011)
  !   .false. = second-order (default)
  !   .true.  = first-order
  logical, allocatable, dimension(:) :: safe_mode
  ! Number of solves before safe_mode is reset to .false., per level
  ! (includes the solve that triggered it)
  integer, allocatable, dimension(:) :: safe_mode_countdown

  ! Multipole coefficients
  real(dp),dimension(1:ndim+1)::multipole

end module poisson_commons
