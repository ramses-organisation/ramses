module rt_parameters
  use hydro_parameters

#ifdef NPACS
  integer,parameter::nPacs=NPACS                  ! # of photon packages (set in Makefile)
#else
  integer,parameter::nPacs=1
#endif
  integer,parameter::nRTvar=nPacs*(1+ndim)    ! # of RT variables (photon density and flux)

  real(dp)::rt_c=0., rt_c2=0.                ! RT constants in user units (set in init_rt)
  real(dp),parameter::c_cgs=2.9979250d+10                  ! Actual lightspeed in [cm s-1]
  real(dp)::rt_c_cgs=c_cgs                                     ! RT lightspeed in [cm s-1]
  real(dp),parameter::m_sun=1.9891d33               ! Solar mass [g], for SED calculations
  real(dp),parameter::eV_to_erg=1.6022d-12          !        eV to erg conversion constant
  real(dp), parameter:: Gyr2sec = 3.15569d+16       !       Gyr to sec conversion constant
  real(dp), parameter:: Myr2sec = 3.15569d+13       !       Myr to sec conversion constant
  real(dp), parameter:: sec2Gyr = 3.16888d-17       !       sec to Gyr conversion constant
  real(dp),parameter:: hp=6.6262d-27                !            Planck const   [erg sec ]
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
  integer,parameter::nIons=NIONS                      ! # of ion species (set in Makefile)
#else
  integer,parameter::nIons=3                        !                     Hii, Heii, Heiii
#endif
  integer::iIons=6                                  !    Starting index of ion states in U
  ! Ionization energies
  real(dp),parameter,dimension(nIons)::ionEvs=(/13.60d0, 24.59d0, 54.42d0/)

  ! RT_PARAMS namelist--------------------------------------------------------------------
  logical::rt_pp=.false.               ! Post-processing (true) or coupled (false)?      !
  logical::rt_cool_only=.false.        ! Don't propagate photons but use them in cooling !
  logical::output_iState=.false.       ! Set to true if initial snapshot is to be output !
  logical::rt_smooth=.false.           ! Smooth the discrete RT update of op. splitting  !
  integer::rt_coupling_method=2        ! 1==full coupling (not working), 2==pp coupling, !
                                       ! 3=adaptive c, 4=rt_dt                           !
  logical::rt_coupling_full=.false.    ! Full (slow and complicated) RT-hydro coupling?  !
  logical::rt_coupling_pp=.false.      ! RT p-processing after each coarse hydro timestep!
  logical::rt_coupling_adc=.false.     ! RT timestep = hydro timestep                    !
  logical::rt_coupling_dtrt=.false.    ! Hydro timestep = RT timestep                    !
  logical::rt_cooling=.true.           ! Activate cooling inside radiation step          !
  logical::rt_star=.false.             ! Activate radiation from star particles?         !
  logical::rt_transp_ISM=.false.       ! Make polytropic cells transparent to RT         !
  real(dp)::rt_esc_frac=1.d0           ! Escape fraction of light from stellar particles !
  logical::rt_is_init_xion=.false.     ! Initialize ionization from T profile?           !
  character(LEN=10)::rt_flux_scheme='glf'                                                !
  logical::rt_use_hll=.false.          ! Use hll flux (or the default glf)               !
  logical::rt_is_outflow_bound=.false. ! Make all boundaries=outflow for RT              !
  real(dp)::rt_max_subcycles=-1        ! Max # of rt-steps in a hydro-step (<1 => no max)!
  real(dp)::rt_courant_factor=0.8d0    ! Courant factor for RT timesteps                 !
  logical::rt_refine=.false.           ! Refine on RT-related conditions?                !
  real(dp)::rt_err_grad_n=-1.0         ! Photon number density gradient for refinement   !
  real(dp)::rt_floor_n=1.d-10          ! Photon number density floor for refinement      !
  real(dp)::rt_err_grad_xHI=-1.0       ! Ionization state gradient for refinement        !
  real(dp)::rt_err_grad_xHII=-1.0      ! Ionization state gradient for refinement        !
  real(dp)::rt_refine_aexp=-1.0        ! Start a for RT gradient refinement              !
  real(dp)::rt_floor_xHI=1.d-10        ! Ionization state floor for refinement           !
  real(dp)::rt_floor_xHII=1.d-10       ! Ionization state floor for refinement           !
  real(dp)::rt_c_fraction=1.d0         ! Actual lightspeed fraction for RT lightspeed    !
  logical::rt_otsa=.true.              ! Use on-the-spot approximation                   !
  logical::rt_UV_hom=.false.           ! Homogeneous UV in every cell?                   !
  real(dp)::rt_UV_nHSS=1d10            ! Self Shielding density threshold for hom UV     !
  logical::rt_isDiffuseUVsrc=.false.   ! UV emission from low-density cells              !
  real(dp)::rt_UVsrc_nHmax=-1.d0       ! Density threshold for UV emission               !x
  logical::rt_freeflow=.false.         ! True => disable photon interaction with gas     !
  logical::upload_equilibrium_x=.false.! Enforce equilibrium xion when uploading         !

  character(LEN=128)::hll_evals_file=''! File HLL eigenvalues                            !
  character(LEN=128)::sed_dir=''       ! Dir containing stellar energy distributions     !
  character(LEN=128)::uv_file=''       ! File containing stellar energy distributions    !

  ! RT_PACS namelist----------------------------------------------------------------------
  integer::sedprops_update=-1                      ! Update sedprops from star populations
  ! negative: never update, 0:update on init, pos x: update every x coarse steps
  ! Sedprops: avg and energy weigthed photoionization c-section (cm2), avg. energy (ev).
  ! Indexes nPacs, nIons stand for photon package vs species (e.g. 1=H, 2=He).
  integer,dimension(nPacs)::iPac=1                               !   Start indices of pacs
  real(dp),dimension(nPacs,nIons)::pac_csn=0, pac_egy=0          !          Cross sections
  real(dp),dimension(nPacs)::pacL0=13.60                         ! Wavelength lower limits
  real(dp),dimension(nPacs)::pacL1=0                             ! Wavelength upper limits
  integer,dimension(nIons)::spec2pac=0                     !Ion -> pac # in recombinations

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
  real(dp),dimension(1:MAXREGION)   ::rt_n_region=0.                      ! Photon density
  real(dp),dimension(1:MAXREGION)   ::rt_u_region=0.                         ! Photon flux
  real(dp),dimension(1:MAXREGION)   ::rt_v_region=0.                         ! Photon flux
  real(dp),dimension(1:MAXREGION)   ::rt_w_region=0.                         ! Photon flux
  real(dp),dimension(1:MAXREGION)   ::rt_xion_region=0.  ! Xion state (square regions only)

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
  integer,dimension(1:MAXREGION)    ::rt_src_pac=1
  real(dp),dimension(1:MAXREGION)   ::rt_n_source=0.                      ! Photon density
  real(dp),dimension(1:MAXREGION)   ::rt_u_source=0.                         ! Photon flux
  real(dp),dimension(1:MAXREGION)   ::rt_v_source=0.                         ! Photon flux
  real(dp),dimension(1:MAXREGION)   ::rt_w_source=0.                         ! Photon flux
  real(dp),dimension(1:MAXREGION)   ::rt_wind_source=0.                     ! Stellar wind

  ! Indexing in flux_module
  integer,parameter::ifrt1=0                                                           ! 0
  integer,parameter::jfrt1=1-ndim/2                                               ! 0 or 1
  integer,parameter::kfrt1=1-ndim/3                                               ! 0 or 1


  ! Cooling statistics: avg loop # per cell, maximum loop #, # of cooling calls-----------
  integer*8::tot_cool_loopcnt=0,max_cool_loopcnt=0,n_cool_cells=0
  integer*8,dimension(4)::loopCodes=0

  ! SED statistics: Radiation emitted, total, last coarse step [#photons/10^50]-----------
  logical::showSEDstats=.true.
  real(dp)::tot_nPhot, step_nPhot, step_nStar, step_mStar

  logical::inLastCoarseStep=.false.    ! .t. when doing last ilevel step in coarse step  !
  logical::doDump = .false.


end module rt_parameters
