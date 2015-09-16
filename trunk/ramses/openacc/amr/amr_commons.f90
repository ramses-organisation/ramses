module amr_commons
  use amr_parameters

  logical::output_done=.false.                  ! Output just performed
  logical::init=.false.                         ! Set up or run
  logical::balance=.false.                      ! Load balance or run
  logical::shrink=.false.                       ! Shrink mesh or run
  integer::nstep=0                              ! Time step
  integer::nstep_coarse=0                       ! Coarse step
  integer::nstep_coarse_old=0                   ! Old coarse step
  integer::nflag,ncreate,nkill                  ! Refinements
  integer::ncoarse                              ! nx.ny.nz
  integer::ngrid_current                        ! Actual number of octs

  real(dp)::emag_tot=0.0D0                      ! Total magnetic energy
  real(dp)::ekin_tot=0.0D0                      ! Total kinetic energy
  real(dp)::eint_tot=0.0D0                      ! Total internal energy
  real(dp)::epot_tot=0.0D0                      ! Total potential energy
  real(dp)::epot_tot_old=0.0D0                  ! Old potential energy
  real(dp)::epot_tot_int=0.0D0                  ! Time integrated potential
  real(dp)::const=0.0D0                         ! Energy conservation
  real(dp)::aexp_old=1.0D0                      ! Old expansion factor
  real(dp)::rho_tot=0.0D0                       ! Mean density in the box
  real(dp)::t=0.0D0                             ! Time variable

  ! executable identification
  CHARACTER(LEN=80)::builddate,patchdir
  CHARACTER(LEN=80)::gitrepo,gitbranch,githash

  ! Save namelist filename
  CHARACTER(LEN=80)::namelist_file

  ! MPI variables
  integer::ncpu,ndomain,myid,overload=1

  ! Friedman model variables
  integer::n_frw
  real(dp),allocatable,dimension(:)::aexp_frw,hexp_frw,tau_frw,t_frw

  ! Initial conditions parameters from grafic
  integer                  ::nlevelmax_part
  real(dp)                 ::aexp_ini=10.
  real(dp),dimension(1:MAXLEVEL)::dfact=1.0d0,astart
  real(dp),dimension(1:MAXLEVEL)::vfact
  real(dp),dimension(1:MAXLEVEL)::xoff1,xoff2,xoff3,dxini
  integer ,dimension(1:MAXLEVEL)::n1,n2,n3

  ! Level related arrays
  real(dp),dimension(1:MAXLEVEL)::dtold,dtnew ! Time step at each level
  real(dp),dimension(1:MAXLEVEL)::rho_max     ! Maximum density at each level
  integer ,dimension(1:MAXLEVEL)::nsubcycle=2 ! Subcycling at each level

  ! Pointers for each level linked list
  integer,allocatable,dimension(:,:)::headl
  integer,allocatable,dimension(:,:)::taill
  integer,allocatable,dimension(:,:)::numbl
  integer(i8b),allocatable,dimension(:,:)::numbtot

  ! Pointers for each level boundary linked list
  integer,allocatable,dimension(:,:)::headb
  integer,allocatable,dimension(:,:)::tailb
  integer,allocatable,dimension(:,:)::numbb

  ! Pointers for free memory grid linked list
  integer::headf,tailf,numbf,used_mem,used_mem_tot

  ! Tree arrays
  real(dp),allocatable,dimension(:,:)::xg      ! grids position
  integer ,allocatable,dimension(:,:)::nbor    ! neighboring father cells
  integer ,allocatable,dimension(:)  ::father  ! father cell
  integer ,allocatable,dimension(:)  ::next    ! next grid in list
  integer ,allocatable,dimension(:)  ::prev    ! previous grid in list
  integer ,allocatable,dimension(:)  ::son     ! sons grids
  integer ,allocatable,dimension(:)  ::flag1   ! flag for refine
  integer ,allocatable,dimension(:)  ::flag2   ! flag for expansion

  ! Global indexing
  integer ,allocatable,dimension(:)  ::cpu_map  ! domain decomposition
  integer ,allocatable,dimension(:)  ::cpu_map2 ! new domain decomposition

  ! Hilbert key
  real(qdp),allocatable,dimension(:)::hilbert_key
  real(qdp),allocatable,dimension(:)::bound_key,bound_key2
  real(qdp)                         ::order_all_min,order_all_max

  ! Recursive bisection                                                                               
  real(dp),allocatable,dimension(:)    ::bisec_wall         ! bisection wall positions                
  integer ,allocatable,dimension(:,:)  ::bisec_next         ! next 2 child cells in bisection         
  integer::bisec_root                                       ! root of bisection tree                  

  integer,allocatable,dimension(:)     ::bisec_indx         ! map from leaf cell id to cpu id         
  real(dp),allocatable,dimension(:,:)  ::bisec_cpubox_min   ! cpu domains boxes                       
  real(dp),allocatable,dimension(:,:)  ::bisec_cpubox_max
  real(dp),allocatable,dimension(:,:)  ::bisec_cpubox_min2  ! cpu domains boxes for new decomp        
  real(dp),allocatable,dimension(:,:)  ::bisec_cpubox_max2

  integer,allocatable,dimension(:)     ::bisec_cpu_load     ! CPU loads (for stats)                   
  integer,allocatable,dimension(:,:)   ::bisec_hist         ! histograms for load computation         
  integer,allocatable,dimension(:)     ::bisec_hist_bounds  ! histogram splitting boundaries          
  integer,allocatable,dimension(:)     ::new_hist_bounds
  integer,allocatable,dimension(:)     ::bisec_ind_cell     ! histo swap id -> cell id map (big)      
  integer,allocatable,dimension(:)     ::cell_level         ! store the level of the cells (big)      

  real(dp)::bisec_res                                       ! resolution parameters                   
  integer ::bisec_nres

  ! Communication structure
  type communicator
     integer                                ::ngrid
     integer                                ::npart
     integer     ,dimension(:)  ,allocatable::igrid
     integer     ,dimension(:,:),allocatable::f
     real(kind=8),dimension(:,:),allocatable::u
     integer(i8b),dimension(:,:),allocatable::fp
     real(kind=8),dimension(:,:),allocatable::up
#ifdef ATON
     real(kind=8),dimension(:,:),pointer::u_radiation
#endif
  end type communicator

  ! Active grid, emission and reception communicators
  type(communicator),allocatable,dimension(:)  ::active
  type(communicator),allocatable,dimension(:,:)::boundary
  type(communicator),allocatable,dimension(:,:)::emission
  type(communicator),allocatable,dimension(:,:)::reception

  ! OpenACC communication buffers
  integer,allocatable,dimension(:,:)           ::nemisspe_acc,nrecpe_acc
  integer,allocatable,dimension(:,:)           ::sum_nemisspe,sum_nrecpe

  ! Types for physical boundary conditions
  CHARACTER(LEN=20)::type_hydro  ='hydro'
  CHARACTER(LEN=20)::type_accel  ='accel'
  CHARACTER(LEN=20)::type_flag   ='flag'

  ! Units specified by the user in the UNITS_PARAMS namelist for non-cosmo runs.
  ! These values shouldn't be used directly. Instead call units() in amr/units.f90.
  real(dp)::units_density=1.0  ! [g/cm^3]
  real(dp)::units_time=1.0     ! [seconds]
  real(dp)::units_length=1.0   ! [cm]

  ! Synchro_hydro_fine subroutine (reduce copies CPU <--> GPU)
  integer,dimension(1:nvector),save::ind_grid_shfv,ind_cell_shfv
  real(dp),dimension(1:nvector),save::pp_shfv

  ! Some compilers may have issues with the reading and writing to many files.
  ! So, this parameter is user-specified and tunes the bug -it is added to ilun-
  integer::ilun_parameter = 500

#ifdef MPROFILING
  !MPROFILING optional
  ! New parts of code (in branches/ramses_gpu)
  real(dp)::acc_t_make_virtual_fine_dp,acc_t_make_virtual_reverse_dp, &
& acc_t_set_unew,acc_t_set_uold,acc_t_add_gravity_source_terms,acc_t_synchro_hydro_fine, &
& acc_t_make_boundary_phi,acc_t_make_boundary_mask,acc_t_force_fine,acc_t_save_phi_old, &
& acc_t_restrict_mask_coarse_reverse,acc_t_cmp_residual_mg_coarse,acc_t_gauss_seidel_mg_coarse, &
& acc_t_restrict_residual_coarse_reverse,acc_t_interpolate_and_correct_coarse, &
& acc_t_set_scan_flag_coarse,acc_t_make_fine_mask,acc_t_make_fine_bc_rhs,acc_t_make_virtual_mg_dp, &
& acc_t_make_reverse_mg_dp,acc_t_restrict_mask_fine_reverse,acc_t_cmp_residual_mg_fine, &
& acc_t_gauss_seidel_mg_fine,acc_t_restrict_residual_fine_reverse,acc_t_interpolate_and_correct_fine,&
& acc_t_set_scan_flag_fine,acc_t_make_initial_phi,acc_t_make_multipole_phi,acc_t_total_time, &
& acc_t_amr_step,acc_t_godunov_fine,acc_t_add_pdv_source_terms,acc_t_multigrid_fine,acc_t_cmp_residual_norm2_fine
#endif

end module amr_commons

