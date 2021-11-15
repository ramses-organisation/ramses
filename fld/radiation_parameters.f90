module radiation_parameters
  use amr_parameters
  use hydro_parameters,only:ngrp,nvar_bicg,nvar_trad

  ! DT adaptatif
  real(dp),allocatable,dimension(:,:)::rad_flux ! Flux entrant dans une cellule
  real(dp),allocatable,dimension(:,:)::urad     ! Old values of Erg in NR iterations
  real(dp),allocatable,dimension(:,:)::frad     ! Radiative force  
  real(dp)::Tr_floor=10.0 ! Background radiation field temperature - WARNING: it affects the pressure_fix in set_uold.
  integer::ntp,nfr

  real(dp):: alpha_imp = 1.0d0	!0.0:explicite 0.5:CN 1.0:implicite
  real(dp):: robin = 1.0d0	!0.0:Von Neumann 1.0:Dirichlet

  ! Multigroup
  integer,parameter::Nomega=100     ! Number of points in the omega data to compute Q moment term

  real(dp),parameter:: aR=7.56591469318689378e-015_dp
  real(dp),parameter::Tray_min=0.5d0 ! Minimum temperature in the radiative energy
  real(dp),parameter:: eray_min=(aR)*Tray_min**4 ! minimum rad energy inside frequency group
  real(dp),parameter:: deray_min=(4.0d0*aR)*Tray_min**3 ! minimum rad energy derivative inside frequency group
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
  logical :: external_radiation_field=.false. ! Default: No external radiation background (@ Tr_floor)
  logical :: stellar_photon=.false.  ! Stellar photons are treated as a separate group (igrp=1). No emission for this group (radiation_source=0)
  logical :: sinks_opt_thin=.false.  ! Optically-thin (min_optical_depth) sink volumes, relevant for hybrid RT

  ! Opacities
  character(len=12) :: opacity_type = 'grey'  ! 'grey' or 'multigroup'
  logical :: sublimation_kuiper=.false. ! Mimicks dust sublimation with decreasing d/g ratio, see Kuiper+10 ApJ 

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

!##################################################################################################
!##################################################################################################
!##################################################################################################
!##################################################################################################

!  Function PLANCK_ANA:
!
!> Compute Planck average opacity.
!<
function planck_ana(dens,Tp,Tr,igroup,insink)

  use amr_commons
  use radiation_parameters

  implicit none

  integer ,intent(in)    :: igroup
  real(dp),intent(in)    :: dens,Tp,Tr
  logical, intent(in)    :: insink
  real(dp)               :: planck_ana,pi,Tgd
  real(dp)               :: Tevap ! if sublimation_kuiper
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  if(sublimation_kuiper) then
     !# RMR #### Sublimation of dust grains as in kuiper+10 ApJ ####
     !### The highest dust temperature is the evaporation temperature
     !### Opacities are taken at this temperature
     !### The input temperature is kept to compute the d/g ratio
     Tgd = Tp
     Tevap = 2000.0d0*dens**0.0195  !! Evaporation temperature
     if(Tp .gt. Tevap) Tgd = Tevap
  endif

  if(sublimation_kuiper) then !! No dust grains above Tevap
     planck_ana = planck_params(1)*(dens**planck_params(2))*(Tgd**planck_params(3))
  else
     planck_ana = planck_params(1)*(dens**planck_params(2))*(Tp**planck_params(3))
  endif  

  if(sublimation_kuiper) then
     !# RMR ## Sublimation mimicked by a d/g ratio that decreases as a arctan function centered on Tevap ##
     pi=acos(-1.0d0)
     planck_ana = planck_ana*(0.5d0 - 1./pi*atan(0.01d0*(Tp - Tevap) ) ) & !planck_ana contains the d/g ratio of 0.01
          +dens*0.01d0*(1.0d0-0.01d0*(0.5d0 - 1./pi*atan(0.01d0*(Tp - Tevap) ) )) !=> quasi full gas
  endif

  if (sinks_opt_thin .and. insink) planck_ana = min_optical_depth/(0.5D0**nlevelmax*boxlen*scale_l)

end function planck_ana

!##################################################################################################
!##################################################################################################
!##################################################################################################
!##################################################################################################

!  Function ROSSELAND_ANA:
!
!> Compute Rosseland mean opacity.
!<
function rosseland_ana(dens,Tp,Tr,igroup,insink)

  use amr_commons
  use radiation_parameters

  implicit none

  integer ,intent(in)    :: igroup
  real(dp),intent(in)    :: dens,Tp,Tr
  logical, intent(in)    :: insink
  real(dp)               :: rosseland_ana,pi,Tgd
  real(dp)               :: Tevap ! if sublimation_kuiper
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  if(sublimation_kuiper) then
     !# RMR #### Sublimation of dust grains as in kuiper+10 ApJ ####
     !### The highest dust temperature is the evaporation temperature
     !### Opacities are taken at this temperature
     !### The input temperature is kept to compute the d/g ratio
     Tgd = Tp
     Tevap = 2000.0d0*dens**0.0195  !! Evaporation temperature
     if(Tp .gt. Tevap) Tgd = Tevap
  endif

  if(sublimation_kuiper) then !! No dust grains above Tevap
     rosseland_ana = rosseland_params(1)*(dens**rosseland_params(2))*(Tgd**rosseland_params(3))
  else
     rosseland_ana = rosseland_params(1)*(dens**rosseland_params(2))*(Tp**rosseland_params(3))
  endif

  if(sublimation_kuiper) then
     !# RMR ## Sublimation mimicked by a d/g ratio that decreases as a arctan function centered on Tevap ##
     pi=acos(-1.0d0)
     rosseland_ana = rosseland_ana*(0.5d0 - 1./pi*atan(0.01d0*(Tp - Tevap) ) ) & !ross_ana contains the d/g ratio of 0.01
          +dens*0.01d0*(1.0d0-0.01d0*(0.5d0 - 1./pi*atan(0.01d0*(Tp - Tevap) ) )) !=> quasi full gas
  endif

  if (sinks_opt_thin .and. insink) rosseland_ana = min_optical_depth/(0.5D0**nlevelmax*boxlen*scale_l)

end function rosseland_ana

!##################################################################################################
!##################################################################################################
!##################################################################################################
!##################################################################################################

!  Function SCATTERING_ANA:
!
!> This routine computes the scattering opacity kappa_s*rho
!! as a function of density and temperature.
!! Units are supposed to be in cgs here (as in units.f90)
!<
function scattering_ana(dens,Tp,Tr,igroup)

  use amr_commons
  use const

  implicit none

  integer ,intent(in)    :: igroup
  real(dp),intent(in)    :: dens,Tp,Tr
  real(dp)               :: scattering_ana
  
  scattering_ana = zero

end function scattering_ana
