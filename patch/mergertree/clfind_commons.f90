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
  real(dp)::relevance_threshold=2.0
  real(dp)::density_threshold=-1.0d0
  real(dp)::saddle_threshold=-1.d0
  real(dp)::rho_clfind=-1.d0
  real(dp)::n_clfind=-1.d0
  real(dp)::mass_threshold=0.d0
  real(dp)::age_cut_clfind=0.d0
  logical::merge_unbound=.false.
  logical::clinfo=.false.
  logical::unbind=.true. !##### NEW HERE

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
  real(dp),allocatable,dimension(:)::clump_mass4
  real(dp),allocatable,dimension(:,:)::peak_pos
  real(dp),allocatable,dimension(:)::relevance
  real(dp),allocatable,dimension(:)::Psurf,MagPsurf,MagTsurf
  real(dp),allocatable,dimension(:)::grav_term, rad_term
  real(dp),allocatable,dimension(:)::clump_check
  real(dp),allocatable,dimension(:)::Icl,Icl_d,Icl_dd
  integer,allocatable,dimension(:)::peak_cell,peak_cell_level
  integer,allocatable,dimension(:)::n_cells,n_cells_halo,lev_peak,new_peak
  integer,allocatable,dimension(:)::occupied,occupied_all,ind_halo
  logical,allocatable,dimension(:)::contracting
!  integer,allocatable,dimension(:)::form,form_all ! Tells whether a sink has to be formed within a clump.

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
  !Particle unbinding related
  !----------------------------

  logical :: unbinding_formatted_output=.false.         !write unformatted output by request

  integer :: nunbound, nunbound_tot, candidates, candidates_tot !counters
  integer :: mergelevel_max                             !deepest merging level
  integer, allocatable, dimension(:)  :: clmppart_first !first particle in particle linked list for each peak id
  integer, allocatable, dimension(:)  :: clmppart_last  !last particle in particle linked list for each peak id
  integer, allocatable, dimension(:)  :: clmppart_next  !next particle in particle linked list for each peak id
  integer, allocatable, dimension(:)  :: nclmppart      !number of particle in particle linked list for each peak id
  
  integer, allocatable,dimension(:)   :: clmpidp        ! ID of peak particle is in
  real(dp),allocatable,dimension(:,:) :: clmp_vel_pb    ! particle based center of mass, clump velocity
  real(dp),allocatable,dimension(:)   :: clmp_mass_pb   ! particle based clump mass
  logical                             :: periodical
  


  !-----------
  !mass bins
  !-----------

  integer :: nmassbins=50
  logical :: logbins=.true.

  real(dp),allocatable,dimension(:,:) :: cmp !cumulative mass profile 
  real(dp),allocatable,dimension(:,:) :: cmp_distances !CMP distances





  !----------------
  !Potential stuff
  !----------------

  logical   :: saddle_pot=.true. !subtract the potential at the closest saddle of the CoM for unbinding
                                 !=considering neighbours for the exclusive unbinding


  real(dp),allocatable,dimension(:)   :: phi_unb !gravitational potential phi
  real(dp)  :: rmin=0.0
  real(dp)  :: GravConst        !gravitational Constant. =factG elsewhere.
  real(dp),allocatable,dimension(:)   :: closest_border ! closest border of clump to the center of mass
                                                        ! stores relative distance squared in each direction
                                                        ! (x^2+y^2+z^2)
  
  !-------------------------
  !Repeated unbinding stuff
  !-------------------------

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

  logical :: make_mergertree = .true.   ! whether to make merger trees
  integer :: nmost_bound = 250          ! maximal number of most bound particles to track

  real(dp), allocatable, dimension(:,:) :: most_bound_energy    ! stores the energy of nmost_bound particles per peak
  integer,  allocatable, dimension(:,:) :: most_bound_pid       ! stores particle ID of nmost_bound particles per peak

  integer,  allocatable, dimension(:)   :: prog_id              ! global ID of progenitors
  integer,  allocatable, dimension(:)   :: prog_owner           ! CPU that owns the progenitor
  real(dp), allocatable, dimension(:)   :: prog_mass            ! list of progenitors masses
  integer,  allocatable, dimension(:)   :: tracers_all          ! list of progenitor tracers (global particle IDs) of all progs
  integer,  allocatable, dimension(:)   :: tracers_loc_pid        ! list of progenitor tracers (local particle IDs)
  integer,  allocatable, dimension(:)   :: tracer_loc_progids_all ! list of progenitor IDs for tracers (local prog ID) of all progs
  integer,  allocatable, dimension(:)   :: tracer_loc_progids     ! list of progenitor IDs for tracers (local prog ID)
                                                                ! only on this CPU
  integer,  allocatable, dimension(:)   :: galaxy_tracers       ! list of active galaxy tracers 
                                                                ! (the absolutely most bound  particle of progenitor) 
  integer,  allocatable, dimension(:)   :: main_prog            ! main progenitor of each descendant 
  integer,  allocatable, dimension(:)   :: main_desc            ! main descendant of each progenitor
                                                                ! (the absolutely most bound  particle of progenitor) 
  integer :: progenitorcount = 0 ! count the number of clumps that will be progenitors
  integer :: progenitorcount_written = 0 ! count the number of progenitors for output
  integer :: nprogs = 0          ! number of progenitors read in/to work with for creating tree
  integer :: prog_free = 1       ! first free progenitor local index
  integer :: ntracers = 0        ! number of tracers on this CPU



  real(dp), allocatable, dimension(:)   :: clmp_mass_exclusive
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
  integer, dimension(:), allocatable :: pmprogs_galaxy  ! Past Merged Progenitors' galaxy particles
  integer, dimension(:), allocatable :: pmprogs_t       ! Time at which past progenitors have been merged (= ifout-1 at merging time)
  integer, dimension(:), allocatable :: pmprogs_owner   ! Current owner Merged Progenitors' galaxy particles
  real(dp),dimension(:), allocatable :: pmprogs_mass    ! Mass of past merged progenitors


  integer :: npastprogs = 0           ! number of past progenitors stored
  integer :: npastprogs_max           ! max number for array loops/allocation
  integer :: pmprog_free = 1          ! First free index in pmprogs* arrays
  integer :: max_past_snapshots = 0   ! maximal number of snapshots to store








end module clfind_commons

