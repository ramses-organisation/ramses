module radiation_parameters
  use amr_parameters
  use hydro_parameters,only:ngrp,nvar_bicg,nvar_trad

  ! DT adaptatif
  real(dp),allocatable,dimension(:,:)::rad_flux ! Flux entrant dans une cellule
  real(dp),allocatable,dimension(:,:)::frad     ! Radiative force  
  real(dp)::Tr_floor=10.0 ! Background radiation field temperature - WARNING: it affects the pressure_fix in set_uold.
  integer::ntp,nfr

  real(dp):: alpha_imp = 1.0d0	!0.0:explicite 0.5:CN 1.0:implicite
  real(dp):: robin = 1.0d0	!0.0:Von Neumann 1.0:Dirichlet

  ! Multigroup
  integer,parameter::Nomega=100     ! Number of points in the omega data to compute Q moment term

  real(dp)::Tray_min=0.5d0 ! Minimum temperature in the radiative energy
  real(dp):: eray_min      ! minimum rad energy inside frequency group
  real(dp):: deray_min     ! minimum rad energy derivative inside frequency group
  real(dp):: small_er=1.0d-30       ! minimum rad energy inside frequency group in code units
  
  real(dp) :: numin=1.0d5,numax=1.0d19 ! Overall frequency boudaries
  real(dp) :: frequency_upperlimit=1.0d35 ! High end frequency if 'extra_end_group = .true.

  integer::Ninv_art4=1000                               ! Number of points in tabulated arT4 function
  real(dp),dimension(:    ),allocatable::dEr_inv_art4   ! Radiative energy increment
  real(dp),dimension(:,:  ),allocatable::inverse_art4_T ! array for tabulated arT4 function dT regular
  real(dp),dimension(:,:,:),allocatable::inverse_art4_E ! array for tabulated arT4 function dE regular
  
  real(dp), dimension(:), allocatable :: nu_min_hz ! minimum freq of given group in Hz
  real(dp), dimension(:), allocatable :: nu_max_hz ! maximum freq of given group in Hz
  real(dp), dimension(:), allocatable :: nu_min_ev ! minimum freq of given group in eV
  real(dp), dimension(:), allocatable :: nu_max_ev ! maximum freq of given group in eV
  
  real(dp),dimension(0:Nomega):: f_array,w1_array,dw1_array,w2_array,dw2_array ! Arrays of omega terms for Q computation

  logical :: freqs_in_Hz=.true.      ! Frequency units in Hz if true; if not eV
  logical :: read_groups=.false.     ! Read group boundaries from file if true
  logical :: split_groups_log=.true. ! Automatic splitting of group in log if true; if not use regular splitting
  logical :: extra_end_group=.false. ! The last group holds frequencies numax -> frequency_upperlimit if true
  logical :: grey_rad_transfer=.true.! Default: grey radiation transfer
  logical :: stellar_photon=.false.  ! Stellar photons are treated as a separate group (igrp=1). No emission for this group (radiation_source=0)
  logical :: sinks_opt_thin=.false.  ! Optically-thin (min_optical_depth) sink volumes, relevant for hybrid RT

  ! Opacities
  character(len=12) :: opacity_type = 'grey'  ! 'grey' or 'multigroup'
  logical :: sublimation_kuiper=.false. ! Mimicks dust sublimation with decreasing d/g ratio, see Kuiper+10 ApJ 
  logical :: fit_semenov=.false.     ! Use a fit of the Semenov et al. (2003) Rosseland anf Planck dust opacities

  ! Radiation solver parameters
  real(dp)::epsilon_diff=1d-6                        ! CG iteration break criteria
  character(LEN=10)::fld_limiter='nolim'             ! Flux limiter (nolim, levermore or minerbo)
  integer::i_fld_limiter
  integer,parameter::i_fld_limiter_nolim=0
  integer,parameter::i_fld_limiter_minerbo=1
  integer,parameter::i_fld_limiter_levermore=2
  integer :: niter=0                                 ! Total number of iteration
  real(dp),dimension(1:10)::dtdiff_params=1d10       ! Conduction time step behaviour
  real(dp),dimension(1:10)::rosseland_params=1.0     ! Rosseland opacity coefficient's parameters
  real(dp),dimension(1:10)::planck_params=1.0        ! Planck opacity coefficient's parameters
  real(dp)::min_optical_depth=1.d-6        ! set the minimum optical depth in the cell (it may accelerate convergence in optically thin regions)

  ! Variables needed for BICG scheme
  real(dp),dimension(:,:,:,:),allocatable :: coeff_glob_left,coeff_glob_right
  real(dp),dimension(:,:,:  ),allocatable :: var_bicg,precond_bicg
  real(dp),dimension(:,:,:  ),allocatable :: mat_residual_glob
  real(dp),dimension(:,:    ),allocatable :: residual_glob
  real(dp),dimension(:,:    ),allocatable :: kappaR_bicg
  logical::block_diagonal_precond_bicg ! if .false. only diagonal, if .true. block diagonal
  integer :: i_rho,i_beta,i_y,i_pAp,i_s
  integer , dimension(1:nvar_bicg) :: ind_bicg
  real(dp), dimension(1:nvar_bicg) :: norm_bicg
  integer , dimension(1:nvar_trad) :: ind_trad
  real(dp), dimension(1:nvar_trad) :: norm_trad
  logical , dimension(1:nvar_trad) :: is_radiative_energy

  integer                                   :: irad_trans_model        !< Integer designating radiative transfer model: 0 = P1, 1 = M1
  integer, parameter                        :: irad_trans_model_p1 = 0 !< P1 radiative transfer model identifier
  integer, parameter                        :: irad_trans_model_m1 = 1 !< M1 radiative transfer model identifier
  integer                                   :: n_points                !< Number of points in the tabulated eigenvalues curve
  real(dp), dimension(:,:,:,:), allocatable :: valp                    !< Array to hold the tabulated eigenvalues as a function of \f$\theta\f$ and \f$\epsilon\f$
  real(dp)                                  :: valp_min=0.0_dp

  logical::store_matrix=.true.

#if USE_FLD==1 && NGRP == 1
  logical, parameter :: bicg_to_cg = .true.
#else
  logical, parameter :: bicg_to_cg = .false.
#endif

end module radiation_parameters
