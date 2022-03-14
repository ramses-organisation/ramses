module amr_parameters

  ! Define real types
  integer,parameter::sp=kind(1.0E0)
#ifndef NPRE
  integer,parameter::dp=kind(1.0E0) ! default
#else
#if NPRE==4
  integer,parameter::dp=kind(1.0E0) ! real*4
#else
  integer,parameter::dp=kind(1.0D0) ! real*8
#endif
#endif
#ifdef QUADHILBERT
  integer,parameter::qdp=kind(1.0_16) ! real*16
#else
  integer,parameter::qdp=kind(1.0_8) ! real*8
#endif
  integer,parameter::MAXOUT=1000
  integer,parameter::MAXLEVEL=100

  ! Define integer types (for particle IDs mostly)
  ! Warning: compiler needs to accept fortran standard 2003.
  ! Specific numbers for fortran kinds are, in principle, implementation
  ! dependent, so "i8b=8" with "integer(i8b)" is implementation-dependent.
  ! See portability note for non-gcc compilers: (boud 2016-11-29) -
  ! https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fFORTRAN_005fENV.html
  ! The selected_int_kind approach below is more likely to be portable:
  integer,parameter::i4b=selected_int_kind(9) ! since log(2*10^9)/log(2)=30.9
#ifndef LONGINT
  !integer,parameter::i8b=4  ! default long int are 4-byte int
  integer,parameter::i8b=selected_int_kind(9) ! since log(2*10^9)/log(2)=30.9
#else
  integer, parameter :: i8b=selected_int_kind(18) ! since log(2*10^18)/log(2)=60.8
  !integer,parameter::i8b=8  ! long int are 8-byte int
#endif /* LONGINT */

  ! Number of dimensions
#ifndef NDIM
  integer,parameter::ndim=1
#else
  integer,parameter::ndim=NDIM
#endif
  integer,parameter::twotondim=2**ndim
  integer,parameter::threetondim=3**ndim
  integer,parameter::twondim=2*ndim

  ! Vectorization parameter
#ifndef NVECTOR
  integer,parameter::nvector=500  ! Size of vector sweeps
#else
  integer,parameter::nvector=NVECTOR
#endif

  integer, parameter :: nstride = 65536

  ! Run control
  logical::verbose =.false.   ! Write everything
  logical::hydro   =.false.   ! Hydro activated
  logical::pic     =.false.   ! Particle In Cell activated
  logical::poisson =.false.   ! Poisson solver activated
  logical::cosmo   =.false.   ! Cosmology activated
  logical::star    =.false.   ! Star formation activated
  logical::sink    =.false.   ! Sink particles activated
  logical::rt      =.false.   ! Radiative transfer activated
  logical::debug   =.false.   ! Debug mode activated
  logical::static  =.false.   ! Static mode activated
  logical::static_dm=.false.  ! Static mode for dm only activated
  logical::static_gas=.false. ! Static mode for gas only activated
  logical::static_stars=.false.! Static mode for stars only activated
  logical::tracer  =.false.   ! Tracer particles activated
  logical::MC_tracer = .false.! Use Monte Carlo tracer particle (https://arxiv.org/abs/1810.11401)
  logical::lightcone=.false.  ! Enable lightcone generation
  logical::clumpfind=.false.  ! Enable clump finder
  logical::unbind=.false.     ! Enable particle unbinding for the clump finder
  logical::make_mergertree=.false. ! Make on the fly mergertrees
  logical::aton=.false.       ! Enable ATON coarse grid radiation transfer

  ! Mesh parameters
  integer::nx=1,ny=1,nz=1                  ! Number of coarse cells in each dimension
  integer::levelmin=1                      ! Full refinement up to levelmin
  integer::nlevelmax=1                     ! Maximum number of level
  integer::ngridmax=0                      ! Maximum number of grids
  integer,dimension(1:MAXLEVEL)::nexpand=1 ! Number of mesh expansion
  integer::nexpand_bound=1                 ! Number of mesh expansion for virtual boundaries
  real(dp)::boxlen=1.0D0                   ! Box length along x direction
  character(len=128)::ordering='hilbert'
  logical::cost_weighting=.true.           ! Activate load balancing according to cpu time
  ! Recursive bisection tree parameters
  integer::nbilevelmax=1                   ! Max steps of bisection partitioning
  integer::nbinodes=3                      ! Max number of internal nodes
  integer::nbileafnodes=2                  ! Max number of leaf (terminal) nodes
  real(dp)::bisec_tol=0.05d0               ! Tolerance for bisection load balancing

                                 ! Step parameters
  integer::nrestart=0            ! New run or backup file number
  integer::nrestart_quad=0       ! Restart with double precision Hilbert keys
  real(dp)::trestart=0           ! Restart time
  logical::restart_remap=.false. ! Force load balance on restart
  integer::nstepmax=1000000      ! Maximum number of time steps
  integer::ncontrol=1            ! Write control variables
  integer::nremap=0              ! Load balancing frequency (0: never)
  integer,allocatable,dimension(:)::remap_pscalar

  ! Output parameters
  integer::iout=1                ! Increment for output times
  integer::ifout=1               ! Increment for output files
  integer::iback=1               ! Increment for backup files
  integer::noutput=1             ! Total number of outputs
  integer::foutput=1000000       ! Frequency of outputs
  logical::gadget_output=.false. ! Output in gadget format
  logical::output_now=.false.    ! write output next step
  real(dp)::walltime_hrs=-1      ! Wallclock time for submitted job
  real(dp)::minutes_dump=1       ! Dump an output minutes before walltime ends

  ! Lightcone parameters
  real(dp)::thetay_cone=12.5d0
  real(dp)::thetaz_cone=12.5d0
  real(dp)::zmax_cone=2d0

  ! Cosmology and physical parameters
  real(dp)::boxlen_ini               ! Box size in h-1 Mpc
  real(dp)::omega_b=0.045d0          ! Omega Baryon
  real(dp)::omega_m=1                ! Omega Matter
  real(dp)::omega_l=0                ! Omega Lambda
  real(dp)::omega_k=0                ! Omega Curvature
  real(dp)::h0     =1                ! Hubble constant in km/s/Mpc
  real(dp)::aexp   =1                ! Current expansion factor
  real(dp)::hexp   =0                ! Current Hubble parameter
  real(dp)::texp   =0                ! Current proper time
  real(dp)::n_sink = -1              ! Sink particle density threshold in H/cc
  real(dp)::rho_sink = -1            ! Sink particle density threshold in g/cc
  real(dp)::d_sink = -1              ! Sink particle density threshold in user units
  real(dp)::m_star =-1               ! Star particle mass in units of mass_sph
  real(dp)::n_star =0.1d0            ! Star formation density threshold in H/cc
  real(dp)::eps_star=0               ! Star formation efficiency
  real(dp)::T2_star=0                ! Typical ISM polytropic temperature
  real(dp)::g_star =1.6d0            ! Typical ISM polytropic index
  real(dp)::jeans_ncells=-1          ! Jeans polytropic EOS
  real(dp)::del_star=200             ! Minimum overdensity to define ISM
  real(dp)::eta_sn =0                ! Supernova mass fraction
  real(dp)::eta_ssn=0.95d0           ! Single supernova ejected mass fraction (sf_imf=.true. only)
  real(dp)::yield  =0                ! Supernova yield
  real(dp)::f_ek   =1                ! Supernovae kinetic energy fraction (only between 0 and 1)
  real(dp)::rbubble=0                ! Supernovae superbubble radius in pc
  real(dp)::f_w    =0                ! Supernovae mass loading factor
  integer ::ndebris=1                ! Supernovae debris particle number
  real(dp)::mass_gmc=-1              ! Stochastic exploding GMC mass
  real(dp)::z_ave  =0                ! Average metal abundance
  real(dp)::B_ave  =0                ! Average magnetic field
  real(dp)::z_reion=8.5d0            ! Reionization redshift
  real(dp)::T2_start                 ! Starting gas temperature
  real(dp)::T2max=huge(1._dp)        ! Temperature ceiling for cooling_fine
  real(dp)::t_delay=1.0d1            ! Feedback time delay in Myr
  real(dp)::t_diss =20               ! Dissipation timescale for feedback
  real(dp)::t_sne =10                ! Supernova blast time
  real(dp)::J21    =0                ! UV flux at threshold in 10^21 units
  real(dp)::a_spec =1                ! Slope of the UV spectrum
  real(dp)::beta_fix=0               ! Pressure fix parameter
  real(dp)::kappa_IR=0               ! IR dust opacity
  real(dp)::ind_rsink=4              ! Number of cells defining the radius of the sphere where AGN feedback is active
  real(dp)::ir_eff=0.75d0            ! efficiency of the IR feedback (only when ir_feedback=.true.)
  real(dp)::sf_trelax=0              ! Relaxation time for star formation (cosmo=.false. only)
  real(dp)::sf_tdiss=0               ! Dissipation timescale for subgrid turbulence in units of turbulent crossing time
  integer::sf_model=3                ! Virial star formation model
  integer::nlevel_collapse=3         ! Number of levels to follow initial dark matter collapse (cosmo=.true. only)
  real(dp)::mass_star_max=120        ! Maximum mass of a star in solar mass
  real(dp)::mass_sne_min=10          ! Minimum mass of a single supernova in solar mass
  integer::momentum_feedback=0       ! Use supernovae momentum feedback if cooling radius not resolved
  integer::strict_equilibrium=0      ! Hydro scheme to preserve exactly hydrostatic equilibrium

  ! PIC dust parameters
  real(dp)::charge_to_mass=0.0       ! Charge to mass ratio for dust grains
  real(dp)::t_stop=0.0               ! Stopping time for dust grains
  real(dp)::stopping_rate=-1.0       ! When greater than or equal to zero, overrides t_stop for constant_t_stop==.true. Allows for zero drag.
  real(dp)::grain_size=0.0           ! Grain size parameter rho_d^i r_d/(rho_g l_0). May wan to get rid of t_stop.
  logical::boris=.false.             ! Activate boris pusher for PIC solver for grain dynamics
  logical::constant_t_stop=.false.    ! Dictates whether stopping time is constant t_stop, or uses grain_size, gas density, velocity, etc.
  logical::second_order=.false.      ! Only works for constant t-stop
  real(dp)::dust_to_gas=1.0          ! Dust-to-gas mass ratio.
  real(dp),dimension(1:3)::accel_gr=0 ! constant external grain force
  integer,dimension(1:MAXOUT)::trajectories=0 ! determines whether or not to output trajectories, which particles to output, and how many.
  logical :: supersonic_drag=.true.   ! if true, Epstein drag is used. If false, drag depends only on density and the sound speed.
  integer :: ndust=1                  ! Determines how many dust grains we has as a multiple of the resolution.
  real(dp):: ddex=0.0                 ! Determines how many decades the dust spectrum spans.
  real(dp):: charge_slope=0.0         ! Determines how the grain charge scales with grain size (power law option)

  logical ::self_shielding=.false.
  logical ::pressure_fix=.false.
  logical ::nordlund_fix=.true.
  logical ::cooling=.false.
  logical ::neq_chem=.false.            ! Non-equilbrium chemistry activated
  logical ::isothermal=.false.          ! Enable equation of state for gas (heating and cooling disabled if .true.)
  logical ::barotropic_eos=.false.      ! New keyword to replace the confusing name "isothermal"
  logical ::metal=.false.
  logical ::haardt_madau=.false.
  logical ::delayed_cooling=.false.
  logical ::smbh=.false.
  logical ::agn=.false.
  logical ::use_proper_time=.false.
  logical ::convert_birth_times=.false. ! Convert stellar birthtimes: conformal -> proper
  logical ::ir_feedback=.false.         ! Activate ir feedback from accreting sinks
  logical ::sf_virial=.false.           ! Activate SF Virial criterion
  logical ::sf_log_properties=.false.   ! Log in ascii files birth properties of stars and supernovae
  logical ::sf_imf=.false.              ! Activate IMF sampling for SN feedback when resolution allows it
  logical ::sf_compressive=.false.      ! Advect compressive and solenoidal turbulence terms separately

  ! EOS parameters
  character(len=20)::barotropic_eos_form='legacy'  !Type of barotropic EOS: choose from:
                                        !'isothermal': constant temperature T0
                                        !'polytrope': T = T0*(rho/rho0)**(gamma-1) or P ~ rho**gamma for ideal gas
                                        !'double_polytrope': isothermal with T0 below rho0 and polytropic with gamma above
                                        !'custom': for patching your own eos
                                        !'legacy': same as polytrop but using the old n_star, g_star and T2_star
  real(dp)::polytrope_rho=1.0d50        ! sets rho0 in EOS = density normalisation or knee-density, in g/cm3
  real(dp)::polytrope_rho_cu=1.0d50     ! rho0 in code units
  real(dp)::polytrope_index=1.0d0       ! sets gamma in EOS = polytropic index
  real(dp)::T_eos=10                    ! sets T0 in EOS: isothermal temperature or temperature normalisation, in K
  real(dp)::mu_gas=1d0                  ! molecular weight
  real(dp)::T2_eos=10                   ! = T/mu, used in the computations

  ! Output times
  real(dp),dimension(1:MAXOUT)::aout=1.1d0      ! Output expansion factors
  real(dp),dimension(1:MAXOUT)::tout=HUGE(1.0D0)! Output times

  ! Movie
  integer,parameter::NMOV=5
  integer::imovout=0             ! Increment for output times
  integer::imov=1                ! Initialize
  real(kind=8)::tstartmov=0,astartmov=0
  real(kind=8)::tendmov=0,aendmov=0
  real(kind=8),allocatable,dimension(:)::amovout,tmovout
  logical::movie=.false.
  integer::nw_frame=512 ! prev: nx_frame, width of frame in pixels
  integer::nh_frame=512 ! prev: ny_frame, height of frame in pixels
  integer::levelmax_frame=0
  real(kind=8),dimension(1:4*NMOV)::xcentre_frame=0d0
  real(kind=8),dimension(1:4*NMOV)::ycentre_frame=0d0
  real(kind=8),dimension(1:4*NMOV)::zcentre_frame=0d0
  real(kind=8),dimension(1:2*NMOV)::deltax_frame=0d0
  real(kind=8),dimension(1:2*NMOV)::deltay_frame=0d0
  real(kind=8),dimension(1:2*NMOV)::deltaz_frame=0d0
  real(kind=8),dimension(1:NMOV)::dtheta_camera=0d0
  real(kind=8),dimension(1:NMOV)::dphi_camera=0d0
  real(kind=8),dimension(1:NMOV)::theta_camera=0d0
  real(kind=8),dimension(1:NMOV)::phi_camera=0d0
  real(kind=8),dimension(1:NMOV)::tstart_theta_camera=0d0
  real(kind=8),dimension(1:NMOV)::tstart_phi_camera=0d0
  real(kind=8),dimension(1:NMOV)::tend_theta_camera=0d0
  real(kind=8),dimension(1:NMOV)::tend_phi_camera=0d0
  real(kind=8),dimension(1:NMOV)::focal_camera=0d0
  real(kind=8),dimension(1:NMOV)::dist_camera=0d0
  real(kind=8),dimension(1:NMOV)::ddist_camera=0d0
  real(kind=8),dimension(1:NMOV)::smooth_frame=1d0
  real(kind=8),dimension(1:NMOV)::varmin_frame=-1d60
  real(kind=8),dimension(1:NMOV)::varmax_frame=1d60
  integer,dimension(1:NMOV)::ivar_frame=0
  logical,dimension(1:NMOV)::perspective_camera=.false.
  logical,dimension(1:NMOV)::zoom_only_frame=.false.
  character(LEN=NMOV)::proj_axis='z' ! x->x, y->y, projection along z
  character(LEN=6),dimension(1:NMOV)::shader_frame='square'
  character(LEN=10),dimension(1:NMOV)::method_frame='mean_mass'
  character(len=10),dimension(1:50)::movie_vars_txt=''
  integer::n_movie_vars
  integer::i_mv_temp=-1,  i_mv_dens=-1,       i_mv_p=-1
  integer::i_mv_speed=-1, i_mv_metallicity=-1
  integer::i_mv_vx=-1,    i_mv_vy=-1,         i_mv_vz=-1
  integer::i_mv_dm=-1,    i_mv_stars=-1,      i_mv_lum=-1
  integer::i_mv_var=-1,   i_mv_xh2=-1,        i_mv_xhi=-1
  integer:: i_mv_xhii=-1, i_mv_xheii=-1,      i_mv_xheiii=-1
  integer::i_mv_fp=-1,    i_mv_pmag=-1
  integer,dimension(1:50)::movie_vars=-1
  integer,dimension(1:50)::movie_var_number=1


                                                 ! Refinement parameters for each level
  real(dp),dimension(1:MAXLEVEL)::m_refine =-1   ! Lagrangian threshold
  real(dp),dimension(1:MAXLEVEL)::r_refine =-1   ! Radius of refinement region
  real(dp),dimension(1:MAXLEVEL)::x_refine = 0   ! Center of refinement region
  real(dp),dimension(1:MAXLEVEL)::y_refine = 0   ! Center of refinement region
  real(dp),dimension(1:MAXLEVEL)::z_refine = 0   ! Center of refinement region
  real(dp),dimension(1:MAXLEVEL)::exp_refine = 2 ! Exponent for distance
  real(dp),dimension(1:MAXLEVEL)::a_refine = 1   ! Ellipticity (Y/X)
  real(dp),dimension(1:MAXLEVEL)::b_refine = 1   ! Ellipticity (Z/X)
  real(dp)::var_cut_refine=-1                    ! Threshold for variable-based refinement
  real(dp)::mass_cut_refine=-1                   ! Mass threshold for particle-based refinement
  integer::ivar_refine=-1                        ! Variable index for refinement
  logical::sink_refine=.false.                   ! Fully refine on sink particles

  ! Initial condition files for each level
  logical::multiple=.false.
  character(LEN=80),dimension(1:MAXLEVEL)::initfile=' '
  character(LEN=20)::filetype='ascii'

  ! Initial condition regions parameters
  integer,parameter::MAXREGION=100
  integer                           ::nregion=0
  character(LEN=10),dimension(1:MAXREGION)::region_type='square'
  real(dp),dimension(1:MAXREGION)   ::x_center=0
  real(dp),dimension(1:MAXREGION)   ::y_center=0
  real(dp),dimension(1:MAXREGION)   ::z_center=0
  real(dp),dimension(1:MAXREGION)   ::length_x=1d10
  real(dp),dimension(1:MAXREGION)   ::length_y=1d10
  real(dp),dimension(1:MAXREGION)   ::length_z=1d10
  real(dp),dimension(1:MAXREGION)   ::exp_region=2

  ! Boundary conditions parameters
  integer,parameter::MAXBOUND=100
  logical                           ::simple_boundary=.false.
  integer                           ::nboundary=0
  integer                           ::icoarse_min=0
  integer                           ::icoarse_max=0
  integer                           ::jcoarse_min=0
  integer                           ::jcoarse_max=0
  integer                           ::kcoarse_min=0
  integer                           ::kcoarse_max=0
  integer ,dimension(1:MAXBOUND)    ::boundary_type=0
  integer ,dimension(1:MAXBOUND)    ::ibound_min=0
  integer ,dimension(1:MAXBOUND)    ::ibound_max=0
  integer ,dimension(1:MAXBOUND)    ::jbound_min=0
  integer ,dimension(1:MAXBOUND)    ::jbound_max=0
  integer ,dimension(1:MAXBOUND)    ::kbound_min=0
  integer ,dimension(1:MAXBOUND)    ::kbound_max=0
  logical                           ::no_inflow=.false.

  ! Number of processes sharing one token
  ! Only one process can write at a time in an I/O group
  integer::IOGROUPSIZE=0           ! Main snapshot
  integer::IOGROUPSIZECONE=0       ! Lightcone
  integer::IOGROUPSIZEREP=0        ! Subfolder size
  logical::withoutmkdir=.false.    ! If true mkdir should be done before the run
  logical::print_when_io=.false.   ! If true print when IO
  logical::synchro_when_io=.false. ! If true synchronize when IO



end module amr_parameters
