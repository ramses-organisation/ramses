module clfind_commons
  use amr_commons, ONLY: qdp,dp
  use sparse_matrix

  integer::ntest,itest                                    !number of cells above threshold per CPU
  integer::ivar_clump=1
  integer::levelmax_clfind
  integer::npeaks,npeaks_tot,npeaks_max
  integer,allocatable,dimension(:)::npeaks_per_cpu
  integer,allocatable,dimension(:)::ipeak_start
  real(dp)::tot_mass
  real(dp)::relevance_threshold=2
  real(dp)::density_threshold=-1
  real(dp)::saddle_threshold=-1
  real(dp)::rho_clfind=-1
  real(dp)::n_clfind=-1
  real(dp)::mass_threshold=0
  real(dp)::age_cut_clfind=0
  logical::merge_unbound=.false.
  logical::clinfo=.false.
  logical::output_clump_field=.false. !output the clump ID field

  ! Peak communicator arrays
  integer,allocatable,dimension(:)::peak_send_cnt,peak_send_oft
  integer,allocatable,dimension(:)::peak_recv_cnt,peak_recv_oft
  integer::peak_recv_tot,peak_send_tot
  integer,allocatable,dimension(:)::peak_send_buf,peak_recv_buf

  ! Spare matrix for saddle points densities
  type(sparse_mat)::sparse_saddle_dens

  ! Hash table variables
  integer::nhash,hfree,hcollision
  integer,dimension(:),allocatable::gkey,nkey,hkey

  ! Peak-patch properties
  real(dp),allocatable,dimension(:,:)::clump_size,center_of_mass,clump_velocity
  real(dp),allocatable,dimension(:,:,:)::Icl_d_3by3,Icl_3by3
  real(dp),allocatable,dimension(:)::max_dens,min_dens,av_dens
  real(dp),allocatable,dimension(:)::thermal_support,kinetic_support,magnetic_support
  real(dp),allocatable,dimension(:)::halo_mass,clump_mass,clump_vol
  real(dp),allocatable,dimension(:)::clump_mass4,clump_star4
  real(dp),allocatable,dimension(:,:)::peak_pos
  real(dp),allocatable,dimension(:)::relevance
  real(dp),allocatable,dimension(:)::Psurf,MagPsurf,MagTsurf
  real(dp),allocatable,dimension(:)::grav_term, rad_term
  real(dp),allocatable,dimension(:)::clump_check
  real(dp),allocatable,dimension(:)::Icl,Icl_d,Icl_dd
  integer,allocatable,dimension(:)::peak_cell,peak_cell_level
  integer,allocatable,dimension(:)::n_cells,n_cells_halo,lev_peak,new_peak
  integer,allocatable,dimension(:)::occupied,occupied_all,ind_halo
  integer,allocatable,dimension(:)::form,form_all
  logical,allocatable,dimension(:)::contracting
  real(dp),allocatable,dimension(:,:)::table_properties

  ! Cell-above-the-threshold properties
  real(dp),allocatable,dimension(:)::denp ! Density of the cells
  integer,allocatable,dimension(:)::imaxp,icellp,levp,testp_sort ! Sort indices

  ! Prime numbers for hash table
  integer,dimension(0:30)::prime=(/2,3,7,13,23,53,97,193,389,769,1543,&
       & 3079,6151,12289,24593,49157,98317,196613,393241,786433,1572869, &
       & 3145739,6291469,12582917,25165843,50331653,100663319,201326611, &
       & 402653189,805306457,1610612741/)



!#####################################################################################
!#######################                        ######################################
!#######################   PARTICLE UNBINDING   ######################################
!#######################                        ######################################
!#####################################################################################


  !----------------------------
  ! Particle unbinding related
  !----------------------------

  logical :: particlebased_clump_output = .false.               ! write resulting clump properties based on particles after
                                                                ! unbinding, not default cell-based properties
  integer :: nunbound, nunbound_tot, candidates, candidates_tot ! counters
  integer :: mergelevel_max                                     ! deepest merging level
  integer, allocatable, dimension(:)  :: clmppart_first         ! first particle in particle linked list for each peak id
  integer, allocatable, dimension(:)  :: clmppart_last          ! last particle in particle linked list for each peak id
  integer, allocatable, dimension(:)  :: clmppart_next          ! next particle in particle linked list for each peak id
  integer, allocatable, dimension(:)  :: nclmppart              ! number of particle in particle linked list for each peak id

  integer, allocatable,dimension(:)   :: clmpidp                ! ID of peak particle is in
  real(dp),allocatable,dimension(:,:) :: clmp_vel_pb            ! particle based clump velocity
  real(dp),allocatable,dimension(:)   :: clmp_mass_pb           ! particle based clump mass
#ifdef UNBINDINGCOM
  real(dp),allocatable,dimension(:,:) :: clmp_com_pb            ! particle based center of mass
#endif
  logical                             :: periodical             ! if simulation is periodical



  !-----------
  ! mass bins
  !-----------

  integer :: nmassbins=50
  logical :: logbins=.true.

  real(dp),allocatable,dimension(:,:) :: cmp !cumulative mass profile
  real(dp),allocatable,dimension(:,:) :: cmp_distances !CMP distances





  !-----------------
  ! Potential stuff
  !-----------------

  logical   :: saddle_pot=.true. ! subtract the potential at the closest saddle of the CoM for unbinding
                                 ! =considering neighbours for the exclusive unbinding


  real(dp),allocatable,dimension(:)   :: phi_unb        ! gravitational potential phi
  real(dp)  :: rmin=0.0
  real(dp)  :: GravConst                                ! gravitational Constant. =factG elsewhere.
  real(dp),allocatable,dimension(:)   :: closest_border ! closest border of clump to the center of mass
                                                        ! stores relative distance squared in each direction
                                                        ! (x^2+y^2+z^2)

  !---------------------------
  ! Iterative unbinding stuff
  !---------------------------

  logical   :: iter_properties=.true.  ! whether to repeat the unbinding with updated clump properties
  real(dp)  :: conv_limit = 0.01       ! convergence factor. If the v_clump_old/v_clump_new < conv_limit,
                                       ! stop iterating for this clump.
  integer   :: repeat_max = 100        ! maximal number of loops per level
  logical   :: loop_again              ! if the loop needs to be redone (for a level)


  logical,allocatable,dimension(:)    :: to_iter           ! whether to repeat the clump properties search
                                                           ! on this peak ID
  logical,allocatable,dimension(:)    :: is_namegiver      ! whether peak is namegiver
  logical,allocatable,dimension(:)    :: contributes       ! whether the particle still contributes
                                                           ! to the clump properties
  integer,allocatable,dimension(:)    :: hasatleastoneptcl ! clump has at least 1 particle that contributes

  ! counters
  integer :: niterunbound, niterunbound_tot




!#####################################################################################
!#######################                        ######################################
!#######################   MERGER TREES         ######################################
!#######################                        ######################################
!#####################################################################################


  !----------------------
  ! Mergertree general
  !----------------------

  logical :: use_exclusive_mass = .true.  ! whether to use exclusive or inclusive halo mass definition for treemaking
  integer :: nmost_bound = 200            ! maximal number of most bound particles to track

  real(dp), allocatable, dimension(:,:) :: most_bound_energy      ! stores the energy of nmost_bound particles per peak
  integer,  allocatable, dimension(:,:) :: most_bound_pid         ! stores particle ID of nmost_bound particles per peak

  integer,  allocatable, dimension(:)   :: prog_id                ! global ID of progenitors
  integer,  allocatable, dimension(:)   :: prog_owner             ! CPU that owns the progenitor
  real(dp), allocatable, dimension(:)   :: prog_mass              ! list of progenitors masses
  integer(i8b),allocatable, dimension(:):: tracers_all            ! list of progenitor tracers (global particle IDs) of all progs
  integer,  allocatable, dimension(:)   :: tracers_loc_pid        ! list of progenitor tracers (local particle IDs)
  integer,  allocatable, dimension(:)   :: tracer_loc_progids_all ! list of progenitor IDs for tracers (local prog ID) of all progs
  integer,  allocatable, dimension(:)   :: tracer_loc_progids     ! list of progenitor IDs for tracers (local prog ID)
                                                                  ! only on this CPU
  integer(i8b),allocatable, dimension(:):: galaxy_tracers         ! list of active galaxy tracers
                                                                  ! (the absolutely most bound  particle of progenitor)
  integer,  allocatable, dimension(:)   :: main_prog              ! main progenitor of each descendant
  integer,  allocatable, dimension(:)   :: main_desc              ! main descendant of each progenitor

  integer :: progenitorcount = 0          ! count the number of clumps that will be progenitors
  integer :: progenitorcount_written = 0  ! count the number of progenitors for output
  integer :: nprogs = 0                   ! number of progenitors read in/to work with for creating tree
  integer :: prog_free = 1                ! first free progenitor local index
  integer :: ntracers = 0                 ! number of tracers on this CPU



  real(dp), allocatable, dimension(:)   :: clmp_mass_exclusive  ! exclusive clump mass, containing only bound particles
  integer,  allocatable, dimension(:)   :: prog_outputnr        ! snapshot number of progenitor
  ! real(dp), allocatable, dimension(:,:) :: clmp_vel_exclusive
  integer  :: killed_tot, appended_tot ! count killed or appended clumps that were too small
  real(dp) :: partm_common






  !-------------------------------
  ! Progenitor-Descendant matrix
  !-------------------------------

  type prog_desc_mat


    ! dimension 1:nprogs
    integer, dimension(:), allocatable :: first ! first descendant for prog
    integer, dimension(:), allocatable :: cnt   ! number of descendants/progenitors

    ! dimension 1:10*npeaks_max
    integer, dimension(:), allocatable :: ntrace  ! number of tracers for this prog/desc pair
    integer, dimension(:), allocatable :: clmp_id ! descendant or progenitor ID. Desc: global; prog: local
    integer, dimension(:), allocatable :: next    ! next descendant for prog

    integer :: mat_free_ind = 1 ! first free index of matrix

  end type prog_desc_mat

  type(prog_desc_mat) :: p2d_links,d2p_links ! sparse matrix for progenitor/descendants linking





  !--------------------------
  ! Multi-shapshot matching
  !--------------------------


  integer, dimension(:), allocatable :: pmprogs         ! Past Merged Progenitors for multi-snapshot matching
  integer(i8b), dimension(:), allocatable :: pmprogs_galaxy  ! Past Merged Progenitors' galaxy particles
  integer, dimension(:), allocatable :: pmprogs_t       ! Time at which past progenitors have been merged (= ifout-1 at merging time)
  integer, dimension(:), allocatable :: pmprogs_owner   ! Current owner Merged Progenitors' galaxy particles
  real(dp),dimension(:), allocatable :: pmprogs_mass    ! Mass of past merged progenitors


  integer :: npastprogs = 0           ! number of past progenitors stored
  integer :: npastprogs_max           ! max number for array loops/allocation
  integer :: pmprog_free = 1          ! First free index in pmprogs* arrays
  integer :: max_past_snapshots = 0   ! maximal number of snapshots to store








!#####################################################################################
!#######################                        ######################################
!#######################   MOCK GALAXIES        ######################################
!#######################                        ######################################
!#####################################################################################


  logical :: make_mock_galaxies = .true.  ! whether to make galaxies

  real(dp), allocatable, dimension(:) :: mpeak                      ! peak mass and expansion factor at time of mass peak
  real(dp), allocatable, dimension(:) :: prog_mpeak                 ! peak mass and expansion factor for progenitors
  real(dp), allocatable, dimension(:) :: pmprogs_mpeak              ! stellar mass of past merged progenitors

  integer, allocatable, dimension(:)  :: orphans_local_pid          ! local particle id of orphans
  integer, allocatable, dimension(:)  :: prog_galaxy_local_id       ! local particle id of progenitor galaxies



!====================================================================
! ADDITIONAL VARIABLES FOR MERGERTREE/UNBINDING DEBUGGING
!====================================================================


#ifdef MTREEDEBUG

  logical :: mtreedebug_no_matrix_dump_prog         = .false. ! don't dump matrices for progenitors
  logical :: mtreedebug_no_matrix_dump_desc         = .false. ! don't dump matrices for descendants
  logical :: mtreedebug_no_unbinding_dump           = .false. ! don't dump unbinding data
  logical :: mtreedebug_no_unbinding_particle_dump  = .true.  ! don't write individual particle IDs in unbinding dumps
  logical :: mtreedebug_no_progdata_dump            = .false. ! don't write progenitor data dump
  logical :: mtreedebug_no_progdata_particle_dump   = .false. ! don't write individual particle IDs in progenitor data dumps
  logical :: mtreedebug_no_pmprogdata_dump          = .false. ! don't write past merged progenitor data dump
  logical :: mtreedebug_no_pmprogdata_particle_dump = .false. ! don't write individual particle IDs in past merged progenitor data dumps
  logical :: mtreedebug_no_mostbound_lists          = .false. ! don't dump most_bound particle lists

#endif


contains


  !================================================================================
  subroutine calc_stellar_mass_params(alpha, gam, delta, loge, logM1, xi, scale_m)
  !================================================================================
    !-----------------------------------------
    ! This subroutine computes the parameters
    ! for the stellar mass at the current
    ! snapshot (current redshift)
    !-----------------------------------------
    use amr_commons, ONLY: dp,aexp
    use cooling_module, ONLY: rhoc
    implicit none
    real(dp), intent(out) :: alpha, gam, delta, loge, logM1, xi, scale_m

    real(dp)      :: M_Sol    = 1.998d33  ! solar mass in g
    logical, save :: seed_set = .false.   ! whether seed for random scatter is set
    real(dp)      :: junk1, junk2, junk3, junk4

    real(dp) :: M_10    =   11.514
    real(dp) :: M_1a    = -  1.793
    real(dp) :: M_1z    = -  0.251
    real(dp) :: e_0     = -  1.777
    real(dp) :: e_a     = -  0.006
    real(dp) :: e_z     =    0.000
    real(dp) :: e_a2    = -  0.119
    real(dp) :: alpha_0 = -  1.412
    real(dp) :: alpha_a =    0.731
    real(dp) :: delta_0 =    3.508
    real(dp) :: delta_a =    2.608
    real(dp) :: delta_z = -  0.043
    real(dp) :: gamma_0 =    0.316
    real(dp) :: gamma_a =    1.319
    real(dp) :: gamma_z =    0.279
    real(dp) :: xi_0    =    0.218
    real(dp) :: xi_a    = -  0.023

    real(dp) :: nu, z, scale_d, scale_l
    integer :: n, clock, i
    integer, dimension(:), allocatable:: seed

    call units(scale_l,junk1,scale_d,junk2,junk3,junk4)
    scale_m = scale_d * scale_l**3 / M_Sol ! get mass in units of M_Sol

    z = 1d0/aexp - 1

    nu = exp(-4d0*aexp**2)
    logM1 = M_10    + nu*(M_1a   *(aexp-1)  + M_1z*z)
    loge  = e_0     + nu*(e_a    *(aexp-1)  + e_z*z ) + e_a2*(aexp-1)
    alpha = alpha_0 + nu*(alpha_a*(aexp-1))
    delta = delta_0 + nu*(delta_a*(aexp-1)            + delta_z*z)
    gam   = gamma_0 + nu*(gamma_a*(aexp-1)            + gamma_z*z)


    !-------------------------
    ! get seed for scatter
    !-------------------------
    ! set unique seed for every task if necessary
    if (.not.seed_set) then
      call random_seed(size=n)
      allocate(seed(n))
      call system_clock(count=clock)
      seed = clock + myid + 53*[(n-i, i=1, n)]
      call random_seed(put=seed)
      deallocate(seed)
      seed_set = .true.
    endif

    xi = (xi_0 + xi_a * (aexp-1))
  end subroutine calc_stellar_mass_params




  !=============================================================================
  real(dp) function stellar_mass(m, alpha, gam, delta, loge, logM1, xi, scale_m)
  !=============================================================================
    !-----------------------------------------------------------
    ! Computes stellar mass given peak mass and expansion
    ! factor a at time when clump had peak mass using a
    ! parametric SHAM relation taken from
    ! Behroozi, Wechsler & Conroy 2013
    ! DOI:	                10.1088/0004-637X/770/1/57
    ! Bibliographic Code:	  2013ApJ...770...57B
    ! http://adsabs.harvard.edu/abs/2013ApJ...770...57B
    !-----------------------------------------------------------
    use amr_commons, ONLY: dp
    implicit none
    real(dp), intent(in) :: m ! mass
    real(dp), intent(in) :: alpha, gam, delta, loge, logM1, xi, scale_m ! parameters computed by subroutine

    real(dp) :: euler   = 2.7182818284590
    real(dp) :: log10_2 = 0.301029995663981 ! log_10(2)

    integer, parameter :: quadp = kind(1.0_8) ! need quad precision for some exponents...
    real(quadp) :: f0, fm, logMM1, sig, xi_scale

    call random_number(xi_scale)  ! temporarily store random number in xi
    call random_number(sig)       ! get sign
    if (sig > 0.5d0) xi_scale = -xi_scale


    ! get f(0)
    f0 = -log10_2 + delta*log10_2**gam/(1d0 + euler)
    ! get log M/M1
    logMM1 = log10(m*scale_m)-logM1
    ! get f(log M/M1)
    fm = -log10(10d0**(alpha*logMM1) + 1) + &
      delta*log10(1d0+exp(logMM1))**gam/(1d0 + exp(10**(-logMM1)))

    stellar_mass = 10d0**DBLE(loge+logM1 + fm - f0 + xi_scale*xi) ! explicitly convert back to double precision
  end function stellar_mass




end module clfind_commons
