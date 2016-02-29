module pm_parameters
  use amr_parameters, ONLY: dp
  integer::nsinkmax=20000           ! Maximum number of sinks
  integer::npartmax=0               ! Maximum number of particles
  integer::npart=0                  ! Actual number of particles
  integer::nsink=0                  ! Actual number of sinks
  integer::iseed=0                  ! Seed for stochastic star formation
  integer::nstar_tot=0              ! Total number of star particle
  real(dp)::mstar_tot=0             ! Total star mass
  real(dp)::mstar_lost=0            ! Missing star mass


  ! More sink related parameters, can all be set in namelist file

  integer::ir_cloud=4                        ! Radius of cloud region in unit of grid spacing (i.e. the ACCRETION RADIUS)
  integer::ir_cloud_massive=3                ! Radius of massive cloud region in unit of grid spacing for PM sinks
  real(dp)::sink_soft=2.d0                   ! Sink grav softening length in dx at levelmax for "direct force" sinks
  real(dp)::mass_sink_direct_force=-1.d0     ! mass above which sinks are treated as "direct force" objects
  
  logical::create_sinks=.false.              ! turn formation of new sinks on

  real(dp)::merging_timescale=-1.d0          ! time during which sinks are considered for merging (only when 'timescale' is used),                                             ! used also as contraction timescale in creation
  real(dp)::cont_speed=0.

  character(LEN=15)::accretion_scheme='none' ! Sink accretion scheme; options: 'none', 'flux', 'bondi', 'threshold'
  logical::flux_accretion=.false.
  logical::threshold_accretion=.false.
  logical::bondi_accretion=.false.

  logical::nol_accretion=.false.             ! Leave angular momentum in the gas at accretion
  real(dp)::mass_sink_seed=0.0               ! Initial sink mass. If < 0, use the AGN feedback based recipe
  real(dp)::c_acc=-1.0                       ! "courant factor" for sink accretion time step control.
                                             ! gives fration of available gas that can be accreted in one timestep.

  logical::eddington_limit=.false.           ! Switch for Eddington limit for the smbh case
  logical::sink_drag=.false.                 ! Gas dragging sink
  logical::clump_core=.false.                ! Trims the clump (for star formation)
  logical::verbose_AGN=.false.               ! Controls print verbosity for the SMBH case
  logical::quench_AGN=.false.                ! Controls if quenching of AGN is present
  real(dp)::acc_sink_boost=1.0               ! Boost coefficient for accretion
  real(dp)::mass_merger_vel_check_AGN=-1.0   ! Threshold for velocity check in  merging; in Msun; default: don't check

  character(LEN=15)::feedback_scheme='energy' ! AGN feedback scheme; options: 'energy' or 'momentum'
  real(dp)::T2_min=1.d7                      ! Minimum temperature of the gas to trigger AGN blast; in K
  real(dp)::T2_max=1.d9                      ! Maximum allowed temperature of the AGN blast; in K
  real(dp)::T2_AGN=1.d12                     ! AGN blast temperature; in K

  real(dp)::v_max=5.d4                       ! Maximum allowed velocity of the AGN blast; in km/s
  real(dp)::v_AGN=1.d4                       ! AGN blast velocity; in km/s
  real(dp)::cone_opening=180.                ! Outflow cone opening angle; in deg

  real(dp)::mass_halo_AGN=1.d10              ! Minimum mass of the halo for sink creation
  real(dp)::mass_clump_AGN=1.d10             ! Minimum mass of the clump for sink creation

end module pm_parameters
