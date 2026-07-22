module turb_parameters
  use amr_parameters

  integer, parameter  :: ILP = selected_int_kind(r=15) ! integer long precision

  ! Precision string
#ifndef NPRE
  integer,parameter::cdp=kind((1.0E0, 1.0E0)) ! default
  character(len=16), parameter :: precision_str='SINGLE_PRECISION'
#else
#if NPRE==4
  integer,parameter::cdp=kind((1.0E0, 1.0E0)) ! complex*4
  character(len=16), parameter :: precision_str='SINGLE_PRECISION'
#else
  integer,parameter::cdp=kind((1.0D0, 1.0D0)) ! complex*8
  character(len=16), parameter :: precision_str='DOUBLE_PRECISION'
#endif
#endif

  ! Turbulence variables
  integer, parameter  :: TURB_GS=64                    ! Turbulent grid size
  integer, parameter  :: TGRID_X=TURB_GS-1             ! Limit of grid, x dimension
#if NDIM>1
  integer, parameter  :: TGRID_Y=TURB_GS-1             ! Limit of grid, x dimension
#else
  integer, parameter  :: TGRID_Y=0                     ! Limit of grid, x dimension
#endif
#if NDIM>2
  integer, parameter  :: TGRID_Z=TURB_GS-1             ! Limit of grid, x dimension
#else
  integer, parameter  :: TGRID_Z=0                     ! Limit of grid, x dimension
#endif
  real(dp), parameter :: turb_gs_real=real(TURB_GS,dp) ! real(TURB_GS, dp)

  logical  :: turb=.FALSE.        ! Use turbulence driving?
  logical  :: turb_evolving=.TRUE. ! Does the driving field evolve over time?
                                  ! .true.  = Ornstein-Uhlenbeck evolution
                                  ! .false. = one static field for the whole run
  integer  :: turb_type=-1        ! Deprecated, replaced by turb_evolving and
                                  ! initial_turb. -1 means "not specified"
  integer  :: turb_seed=-1        ! Turbulent seed (-1=random)
  logical  :: instant_turb=.TRUE. ! Generate initial turbulence before start?
  character (LEN=100) :: forcing_power_spectrum='parabolic'
                                  ! Power spectrum type of turbulent forcing

  real(dp) :: comp_frac=0.3333_dp ! Compressive fraction
  real(dp) :: turb_T=1.0_dp       ! Turbulent velocity autocorrelation time
  integer  :: turb_Ndt=100        ! Number of timesteps per autocorr. time
  real(dp) :: turb_rms=1.0_dp     ! rms turbulent forcing acceleration
  logical  :: turb_exact_rms=.FALSE. ! Set the forcing rms to exactly match turb_rms?

  real(dp) :: turb_min_rho=1d-50  ! Minimum density for turbulence

  ! Turbulent initial velocity field, independent of the driving above
  logical  :: initial_turb=.FALSE.       ! Add turbulent velocity to the initial conditions?
  real(dp) :: initial_turb_vrms=0.0_dp   ! Velocity dispersion |v|_rms of that field
  real(dp) :: initial_turb_comp_frac=0.3333_dp
                                  ! Compressive fraction of the initial velocity
                                  ! field. 1/3 spreads the power evenly over the
                                  ! one longitudinal and two transverse degrees
                                  ! of freedom
  character (LEN=100) :: initial_turb_spectrum='power_law'
                                  ! Power spectrum of the initial velocity field.
                                  ! Broadband by default: initial conditions want
                                  ! power on every scale, whereas the driving
                                  ! usually acts on a few large-scale modes

end module turb_parameters
