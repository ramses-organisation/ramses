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
  real(dp)::msink_direct=-1.d0               ! mass above which sinks are treated as "direct force" objectfs
  
  logical::create_sinks=.false.              ! turn formation of new sinks on

  character(LEN=15)::merging_scheme='none'   ! sink merging scheme. options: 'none,'timescale', 'FOF'
  real(dp)::merging_timescale=-1.d0          ! time during which sinks are considered for merging (only when 'timescale' is used), 
                                             ! used also as contraction timescale in creation
  real(dp)::cont_speed=0.

  character(LEN=15)::accretion_scheme='none' ! sink accretion scheme. options: 'none', 'flux', 'bondi', 'threshold'
  logical::flux_accretion=.false.
  logical::threshold_accretion=.false.
  logical::bondi_accretion=.false.

  logical::nol_accretion=.false.              ! Leave angular momentum in the gas at accretion
  real(dp)::sink_seedmass=5.4d-4             ! Initial mass sinks are created with in bondi or flux accretion case (in solar masses)
  real(dp)::c_acc=-1.0                       ! "courant factor" for sink accretion time step control.
                                             ! gives fration of available gas that can be accreted in one timestep.
  logical::diffuse_acczone=.false.           ! If set to true, the difmag parameter will only be activated inside the sink accretion zone.
                                             ! This can help to prevent crashes because of negative densities, especially when
                                             ! nol_accretion=.true.

end module pm_parameters
