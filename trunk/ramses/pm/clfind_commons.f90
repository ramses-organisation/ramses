module clfind_commons
  use amr_commons, ONLY: qdp,dp
  use sparse_matrix

  integer::ivar_clump=1
  integer::npeaks,npeaks_tot,npeaks_max
  integer,allocatable,dimension(:)::npeaks_per_cpu
  integer,allocatable,dimension(:)::ipeak_start
  real(dp)::tot_mass
  real(dp)::relevance_threshold=1.5
  real(dp)::density_threshold=0.
  real(dp)::saddle_threshold=-1
  real(dp)::mass_threshold=0.
  logical::merge_unbound=.false.,clinfo=.false.
  logical::first_pass

  ! Peak communicator arrays
  integer,allocatable,dimension(:)::peak_send_cnt,peak_send_oft
  integer,allocatable,dimension(:)::peak_recv_cnt,peak_recv_oft
  integer::peak_recv_tot,peak_send_tot
  integer,allocatable,dimension(:)::peak_send_buf,peak_recv_buf

  ! Spare matrix for saddle points densities
  type(sparse_mat)::sparse_saddle_dens

  ! Hash table variables
  integer::nhash,hfree
  integer,dimension(:),allocatable::gkey,nkey,hkey

  ! Peak-patch properties
  real(dp),allocatable,dimension(:,:)::clump_size,center_of_mass,clump_momentum,clump_force
  real(dp),allocatable,dimension(:,:,:)::second_moments,Icl_d_3by3,Icl_3by3
  real(dp),allocatable,dimension(:)::min_dens,av_dens,phi_min
  real(dp),allocatable,dimension(:)::max_dens,e_kin_int,e_bind,e_thermal
  real(dp),allocatable,dimension(:)::halo_mass,clump_mass,clump_vol,clump_mass4
  real(dp),allocatable,dimension(:,:)::peak_pos,bulk_momentum
  real(dp),allocatable,dimension(:)::saddle_max
  real(dp),allocatable,dimension(:)::relevance
  real(dp),allocatable,dimension(:)::phi_ref
  real(dp),allocatable,dimension(:)::Psurf,v_therm,v_rms,m4
  real(dp),allocatable,dimension(:)::e_bind_iso,e_therm_iso,e_kin_iso,grav_term,clump_virial
  real(dp),allocatable,dimension(:)::peak_check,ball4_check,isodens_check,clump_check
  real(dp),allocatable,dimension(:)::Icl,Icl_d,Icl_dd
  integer,allocatable,dimension(:)::n_cells,lev_peak,minmatch,new_peak
  integer,allocatable,dimension(:)::occupied,occupied_all,ind_halo
  integer,allocatable,dimension(:)::form,form_all
  logical,allocatable,dimension(:)::contracting

  ! Cell-above-the-threshold properties
  real(dp),allocatable,dimension(:)::denp ! Density of the cells
  integer,allocatable,dimension(:)::imaxp,icellp,levp,testp_sort ! Sort indices

  ! Prime numbers for hash table
  integer,dimension(0:30)::prime=(/2,3,7,13,23,53,97,193,389,769,1543,&
       & 3079,6151,12289,24593,49157,98317,196613,393241,786433,1572869, &
       & 3145739,6291469,12582917,25165843,50331653,100663319,201326611, &
       & 402653189,805306457,1610612741/)

end module clfind_commons
