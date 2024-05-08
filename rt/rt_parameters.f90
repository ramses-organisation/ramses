module rt_parameters
  use hydro_parameters
  use constants

#ifdef NGROUPS
  integer,parameter::nGroups=NGROUPS          ! # of photon groups (set in Makefile)
#else
  integer,parameter::nGroups=1
#endif
  integer,parameter::nRTvar=nGroups*(1+ndim) ! # of RT variables (photon density and flux)

  real(dp)::rt_c=0., rt_c2=0.                ! RT constants in user units (set in init_rt)
  real(dp)::rt_c_cgs=c_cgs                                 !   Reduced lightspeed [cm s-1]
  real(dp),parameter::one_over_c_cgs=3.335640484668562d-11 !         Save some computation
  real(dp),allocatable,dimension(:,:)::lambda1,lambda4                   ! HLL eigenvalues
#ifndef NPRE
  real(dp),parameter::smallNp=1d-30                 !               Minimum photon density
#else
#if NPRE==4
  real(dp),parameter::smallNp=1d-30                 !               Minimum photon density
#else
  real(dp),parameter::smallNp=1d-50                 !               Minimum photon density
#endif
#endif
  ! Ion species---------------------------------------------------------------------------
#ifdef NIONS
  integer,parameter::nIons=NIONS                    !   # of ion species (set in Makefile)
#else
  integer,parameter::nIons=3                        !   HII and optionally HI, HeII, HeIII
#endif
  integer::nIonsUsed=0                              ! # species used (as opposed to alloc)
  integer::iIons=6                                  !    Starting index of ion states in U
  ! Ionization energies
  real(dp),dimension(nIons)::ionEvs                 !                       Set in rt_init
  real(dp),parameter::ionEv_HI    = 11.20d0
  real(dp),parameter::ionEv_HII   = 13.60d0
  real(dp),parameter::ionEv_HeII  = 24.59d0
  real(dp),parameter::ionEv_HeIII = 54.42d0
  logical::isHe=.true.                              !     He ionization fractions tracked?
  logical::isH2=.false.                             !                          H2 tracked?
  integer::ixHI=0, ixHII=0, ixHeII=0, ixHeIII=0     !      Indexes of ionization fractions
  logical::is_SED_single_Z=.false.                  !                   Tracks if single Z
                                                    !                 in SED interpolation

  ! RT_PARAMS namelist--------------------------------------------------------------------
  logical::rt_advect=.false.           ! Advection of photons?                           !
  logical::rt_smooth=.false.           ! Smooth the discrete RT update of op. splitting  !
  real(dp)::rt_Tconst=-1               ! If pos. use this value for all T-depend. rates  !
  logical::rt_isTconst=.false.         ! Const rates activated?                          !
  logical::rt_star=.false.             ! Activate radiation from star particles?         !
  logical::rt_AGN=.false.              ! Activate radiation from sink particles on central cloud !
  logical::rt_sink=.false.             ! Activate radiation from sinks
  real(dp)::rt_esc_frac=1d0            ! Escape fraction of light from stellar particles !
  logical::rt_is_init_xion=.false.     ! Initialize ionization from T profile?           !
  character(LEN=10)::rt_flux_scheme='glf'                                                !
  logical::rt_use_hll=.false.          ! Use hll flux (or the default glf)               !
  logical::rt_is_outflow_bound=.false. ! Make all boundaries=outflow for RT              !
  real(dp)::rt_courant_factor=0.8d0    ! Courant factor for RT timesteps                 !
  logical::rt_refine=.false.           ! Refine on RT-related conditions?                !
  real(dp)::rt_err_grad_n=-1.0         ! Photon number density gradient for refinement   !
  real(dp)::rt_floor_n=1d-10           ! Photon number density floor for refinement      !
  real(dp)::rt_err_grad_xHI=-1.0       ! Ionization state gradient for refinement        !
  real(dp)::rt_err_grad_xHII=-1.0      ! Ionization state gradient for refinement        !
  real(dp)::rt_refine_aexp=-1.0        ! Start a for RT gradient refinement              !
  real(dp)::rt_floor_xHI=1d-10         ! Ionization state floor for refinement           !
  real(dp)::rt_floor_xHII=1d-10        ! Ionization state floor for refinement           !
  real(dp)::rt_c_fraction=1d0          ! Actual lightspeed fraction for RT lightspeed    !
  integer::rt_nsubcycle=1              ! Maximum number of RT-steps during one hydro/    !
                                       ! gravity/etc timestep                            !
  logical::rt_otsa=.true.              ! Use on-the-spot approximation                   !
  logical::rt_isDiffuseUVsrc=.false.   ! UV emission from low-density cells              !
  real(dp)::rt_UVsrc_nHmax=-1d0        ! Density threshold for UV emission               !
  logical::upload_equilibrium_x=.false.! Enforce equilibrium xion when uploading         !
  integer::heat_unresolved_HII=0       ! Subgrid model heating unresolved HII regions    !
  integer::iHIIheat=6                  ! Var index for HII heating                       !
  logical::cosmic_rays=.false.         ! Include cosmic ray ionisation                   !

  character(LEN=128)::hll_evals_file=''! File HLL eigenvalues                            !
  character(LEN=128)::sed_dir=''       ! Dir containing stellar energy distributions     !
  character(LEN=128)::uv_file=''       ! File containing stellar energy distributions    !

  ! RT_GROUPS namelist--------------------------------------------------------------------
  integer::sedprops_update=-1                      ! Update sedprops from star populations
  ! negative: never update, 0:update on init, pos x: update every x coarse steps
  logical::SED_isEgy=.false. ! Integrate energy out of SEDs rather than photon count
  ! Grop props: avg and energy weigthed photoionization c-section (cm2), avg. energy (ev).
  ! Indexes nGroups, nIons stand for photon group vs species (e.g. 1=H, 2=He).
  integer,dimension(nGroups)::iGroups=1                          ! Start indices of groups
  real(dp),dimension(nGroups,nIons)::group_csn=0, group_cse=0    !    Cross sections (cm2)
  real(dp),dimension(nGroups)::group_egy=0                       !  Avg photon energy (ev)
  real(dp),dimension(nGroups)::group_egy_AGNfrac=0               !  Fraction of AGN energy
  real(dp),dimension(nGroups)::groupL0=13.60                     ! Wavelength lower limits
  real(dp),dimension(nGroups)::groupL1=0                         ! Wavelength upper limits
  integer,dimension(nIons)::spec2group=0                 !Ion -> group # in recombinations

  ! Imposed boundary condition variables
  real(dp),dimension(1:MAXBOUND,1:nrtvar)::rt_boundary_var
  real(dp),dimension(1:MAXBOUND)::rt_n_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::rt_u_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::rt_v_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::rt_w_bound=0.0d0

  ! Initial condition RT regions parameters----------------------------------------------
  integer                           ::rt_nregion=0
  character(LEN=10),dimension(1:MAXREGION)::rt_region_type='square'
  real(dp),dimension(1:MAXREGION)   ::rt_reg_x_center=0.
  real(dp),dimension(1:MAXREGION)   ::rt_reg_y_center=0.
  real(dp),dimension(1:MAXREGION)   ::rt_reg_z_center=0.
  real(dp),dimension(1:MAXREGION)   ::rt_reg_length_x=1.E10
  real(dp),dimension(1:MAXREGION)   ::rt_reg_length_y=1.E10
  real(dp),dimension(1:MAXREGION)   ::rt_reg_length_z=1.E10
  real(dp),dimension(1:MAXREGION)   ::rt_exp_region=2.0
  integer,dimension(1:MAXREGION)    ::rt_reg_group=1
  real(dp),dimension(1:MAXREGION)   ::rt_n_region=0.                      ! Photon density
  real(dp),dimension(1:MAXREGION)   ::rt_u_region=0.                         ! Photon flux
  real(dp),dimension(1:MAXREGION)   ::rt_v_region=0.                         ! Photon flux
  real(dp),dimension(1:MAXREGION)   ::rt_w_region=0.                         ! Photon flux

   ! RT source regions parameters----------------------------------------------------------
  integer                           ::rt_nsource=0
  character(LEN=10),dimension(1:MAXREGION)::rt_source_type='square'
  real(dp),dimension(1:MAXREGION)   ::rt_src_x_center=0.
  real(dp),dimension(1:MAXREGION)   ::rt_src_y_center=0.
  real(dp),dimension(1:MAXREGION)   ::rt_src_z_center=0.
  real(dp),dimension(1:MAXREGION)   ::rt_src_length_x=1.E10
  real(dp),dimension(1:MAXREGION)   ::rt_src_length_y=1.E10
  real(dp),dimension(1:MAXREGION)   ::rt_src_length_z=1.E10
  real(dp),dimension(1:MAXREGION)   ::rt_exp_source=2.0
  integer, dimension(1:MAXREGION)   ::rt_src_group=1
  real(dp),dimension(1:MAXREGION)   ::rt_n_source=0.                      ! Photon density
  real(dp),dimension(1:MAXREGION)   ::rt_u_source=0.                         ! Photon flux
  real(dp),dimension(1:MAXREGION)   ::rt_v_source=0.                         ! Photon flux
  real(dp),dimension(1:MAXREGION)   ::rt_w_source=0.                         ! Photon flux

  ! Indexing in flux_module
  integer,parameter::ifrt1=0                                                           ! 0
  integer,parameter::jfrt1=1-ndim/2                                               ! 0 or 1
  integer,parameter::kfrt1=1-ndim/3                                               ! 0 or 1

  ! Cooling statistics: avg loop # per cell, maximum loop #, # of cooling calls-----------
  logical::rt_output_coolstats=.false.    ! Output cooling statistics                     !
  integer*8::tot_cool_loopcnt=0,max_cool_loopcnt=0,n_cool_cells=0
  integer*8,dimension(20)::loopCodes=0

  ! SED statistics: Radiation emitted, total, last coarse step [#photons/10^50]-----------
  logical::showSEDstats=.true.
  real(dp)::tot_nPhot, step_nPhot, step_nStar, step_mStar

  ! H2 parameters ------------------------------------------------------------------------
  ! Self-shielding factor, see Nickerson, Teyssier, & Rosdahl (2018)
  ! Array to track which groups are in the Lyman-Werner band, 11.2 eV to 13.6 eV
  real(dp),dimension(1:NGROUPS)::ssh2 = 1d0, isLW = 0d0
end module rt_parameters
