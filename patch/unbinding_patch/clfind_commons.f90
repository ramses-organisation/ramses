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
  
  integer, allocatable,dimension(:)   :: clmpidp                   ! ID of peak particle is in
  real(dp),allocatable,dimension(:,:) :: clmp_com_pb,clmp_vel_pb   ! particle based center of mass, clump velocity
  real(dp),allocatable,dimension(:)   :: clmp_mass_pb              ! particle based clump mass
  !real(dp),allocatable,dimension(:,:) :: periodicity_correction
  logical                             :: periodical
  


  !-----------
  !mass bins
  !-----------

  integer :: nmassbins=50
  logical :: logbins=.true.

  real(dp),allocatable,dimension(:,:) :: cmp !cumulative mass profile, 
  real(dp),allocatable,dimension(:,:) :: cmp_distances !CMP distances





  !----------------
  !Potential stuff
  !----------------

  !logical   :: mp_pot=.false.   !unbinding potential calc: consider enclosed mass as point source.
  logical   :: saddle_pot=.true.!subtract the potential at the closest saddle of the CoM for unbinding


  real(dp),allocatable,dimension(:)   :: phi_unb !gravitational potential phi
  real(dp)  :: rmin=0.0
                                !=considering neighbours for the exclusive unbinding
  real(dp)  :: GravConst        !gravitational Constant. =factG elsewhere.
  real(dp),allocatable,dimension(:)   :: closest_border ! closest border of clump to the center of mass
                                                        ! stores relative distance squared in each direction
                                                        ! (x^2+y^2+z^2)
  
  !-------------------------
  !Repeated unbinding stuff
  !-------------------------

  logical   :: iter_properties=.true.  ! whether to repeat the unbinding with updated clump properties
  real(dp)  :: conv_limit = 0.01          ! convergence factor. If the v_clump_old/v_clump_new < conv_limit,
                                        ! stop iterating for this clump. 
  integer   :: repeat_max = 100          ! maximal number of loops per level
  logical   :: loop_again               ! if the loop needs to be redone (for a level)


  logical,allocatable,dimension(:)    :: to_iter         ! whether to repeat the clump properties search
                                                         ! on this peak ID
  logical,allocatable,dimension(:)    :: contributes     ! whether the particle still contributes
                                                         ! to the clump properties

  real(dp),allocatable,dimension(:,:) :: oldcom, oldvel     ! store old values: centre of mass, bulk velocits
  real(dp),allocatable,dimension(:)   :: oldcmpd,oldm       ! particle furthest away from CoM
  integer,allocatable,dimension(:)    :: hasatleastoneptcl  ! clump has at least 1 particle that contributes


    ! debugging
  integer :: testpeak
  integer :: niterunbound, niterunbound_tot

end module clfind_commons

