module pm_parameters
  use amr_parameters, ONLY: dp
  integer::nsinkmax=2000            ! Maximum number of sinks
  integer::npartmax=0               ! Maximum number of particles
  integer::npart=0                  ! Actual number of particles
  integer::nsink=0                  ! Actual number of sinks
  integer::iseed=0                  ! Seed for stochastic star formation
  integer::tseed=0                  ! Seed for MC tracers
  integer::nstar_tot=0              ! Total number of star particles
  real(dp)::mstar_tot=0             ! Total star mass
  real(dp)::mstar_lost=0            ! Missing star mass

  integer::ntracer_tot=0            ! Total number of tracers
  ! More sink related parameters, can all be set in namelist file

  integer::ir_cloud=4                        ! Radius of cloud region in unit of grid spacing (i.e. the ACCRETION RADIUS)
  integer::ir_cloud_massive=4                ! Radius of massive cloud region in unit of grid spacing for PM sinks
  real(dp)::sink_soft=2                      ! Sink grav softening length in dx at levelmax for "direct force" sinks
  real(dp)::mass_sink_direct_force=-1        ! mass above which sinks are treated as "direct force" objects
  integer::nlevelmax_sink=0                  ! HACK to put sinks at coarser level (for sims which are not fully refined)

  logical::create_sinks=.false.              ! turn formation of new sinks on
  logical::check_energies=.true.             ! when flagging clumps for sink formation, check whether their gravitational energy is dominant

  real(dp)::merging_timescale=-1             ! time during which sinks are considered for merging (only when 'timescale' is used),
                                             ! used also as contraction timescale in creation
  real(dp)::cont_speed=0                     ! Clump contraction rate

  character(LEN=15)::accretion_scheme='none' ! Sink accretion scheme; options: 'none', 'bondi', 'threshold'
  logical::threshold_accretion=.false.       ! NOT A NAMELIST PARAMETER
  logical::bondi_accretion=.false.           ! NOT A NAMELIST PARAMETER
  logical::bondi_use_vrel=.true.             ! Use v_rel^2 in the denominator of Bondi formula
  real(dp)::c_acc=0.75                       ! "courant factor" for sink accretion
                                             ! gives fraction of available gas that can be accreted in one timestep
  real(dp)::mass_sink_seed=0                 ! Initial sink mass
  real(dp)::mass_smbh_seed=0                 ! Initial SMBH mass
  real(dp)::mass_merger_vel_check=1d100      ! Threshold for velocity check in merging in M_sun; default: don't check

  logical::eddington_limit=.false.           ! Switch for Eddington limit for the smbh case
  logical::clump_core=.false.                ! Trims the clump (for star formation)
  logical::verbose_AGN=.false.               ! Controls print verbosity for the SMBH case
  real(dp)::acc_sink_boost=1                 ! Boost coefficient for accretion

  real(dp)::AGN_fbk_frac_ener=1              ! Fraction of AGN feedback released as thermal blast
  real(dp)::AGN_fbk_frac_mom=0               ! Fraction of AGN feedback released as momentum injection

  real(dp)::T2_min=1d7                      ! Minimum temperature of the gas to trigger AGN blast; in K
  real(dp)::T2_max=1d9                      ! Maximum allowed temperature of the AGN blast; in K
  real(dp)::T2_AGN=1d12                     ! AGN blast temperature; in K
  real(dp)::v_max=2000                      ! Maximum allowed velocity of the AGN jet; in km/s

  real(dp)::cone_opening=180d0              ! Outflow cone opening angle; in deg
  real(dp)::epsilon_kin=1                    ! Efficiency of kinetic feedback
  real(dp)::kin_mass_loading=100d0           ! Mass loading of the jet
  real(dp)::AGN_fbk_mode_switch_threshold=0.01d0 ! M_Bondi/M_Edd ratio to switch between feedback modes
                                                 ! if rate gt <value> is thermal, else is momentum;
                                                 ! if <value> le 0 then not active

  real(dp)::mass_halo_AGN=1d10              ! Minimum mass of the halo for sink creation
  real(dp)::mass_clump_AGN=1d10             ! Minimum mass of the clump for sink creation
  real(dp)::mass_star_AGN=0d0               ! Minimum mass of stars in the clump for sink creation

  real(dp)::boost_threshold_density=0.1d0   ! Accretion boost threshold for Bondi

  real(dp)::max_mass_nsc=1d15               ! Maximum mass of the Nuclear Star Cluster (msink)

  logical::sink_descent=.false.             ! Switch for the sink descent
  real(dp)::gamma_grad_descent=0.0d0        ! Step for the gradient descent
  real(dp)::fudge_graddescent=1.0d0         ! Fudge factor for the for the BB gradient descent

  character(LEN=15)::agn_acc_method='mass'
  character(LEN=15)::agn_inj_method='volume'

  type part_t
     ! We store these two things contiguously in memory
     ! because they are fetched at similar times
     integer(1) :: family
     integer(1) :: tag
  end type part_t

  ! MC Tracer
  character(LEN=1025) :: tracer_feed             ! Filename to read the tracer from
  character(LEN=  10) :: tracer_feed_fmt='ascii' ! Format of the input (ascii or binary)
  real(dp)::tracer_mass=-1.0                     ! Mass of the tracers, used for outputs and seed

  integer :: tracer_first_balance_levelmin = -1  ! Set to >0 to add more weight on level finer than this
  integer :: tracer_first_balance_part_per_cell = 0 ! Typical initial number of parts per cell

end module pm_parameters
