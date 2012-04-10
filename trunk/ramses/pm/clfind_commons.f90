module clfind_commons




  integer::nparts,nparts_tot,npeaks,npeaks_tot
  real(kind=8)::tot_mass,relevance_threshold,density_threshold,mass_threshold
  

  !Big array for saddlepoint values
  real(kind=8),allocatable,dimension(:,:)::saddle_dens_tot


  !peak_batch properties
  real(kind=8),allocatable,dimension(:,:)::clump_size_tot,center_of_mass_tot,clump_momentum_tot
  real(kind=8),allocatable,dimension(:,:,:)::second_moments,second_moments_tot
  real(kind=8),allocatable,dimension(:)::min_dens_tot,av_dens_tot,phi_min_tot
  real(kind=8),allocatable,dimension(:)::max_dens_tot,e_kin_int_tot,e_bind_tot,e_thermal_tot
  real(kind=8),allocatable,dimension(:)::e_kin_int_tot4,e_bind_tot4,e_thermal_tot4
  real(kind=8),allocatable,dimension(:)::clump_mass_tot,clump_vol_tot
  real(kind=8),allocatable,dimension(:,:)::peak_pos_tot
  real(kind=8),allocatable,dimension(:)::saddle_max_tot
  real(kind=8),allocatable,dimension(:)::relevance_tot
  integer,allocatable,dimension(:)::n_cells_tot,minmatch_tot

  integer,allocatable,dimension(:)::sort_index
  integer,allocatable,dimension(:)::occupied,occupied_all !tells wheter there is already a sink in a clump


end module clfind_commons
