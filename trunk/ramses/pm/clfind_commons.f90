module clfind_commons
  use amr_commons, ONLY: qdp

  integer::nparts,nparts_tot,npeaks,npeaks_tot
  real(kind=8)::tot_mass
  real(kind=8)::relevance_threshold=1.5
  real(kind=8)::density_threshold=0.
  real(kind=8)::mass_threshold=0.

  ! Big array for saddlepoint values
  real(kind=8),allocatable,dimension(:,:)::saddle_dens,saddle_dens_tot


  ! Peak patch properties
  real(kind=8),allocatable,dimension(:,:)::clump_size_tot,center_of_mass_tot,clump_momentum_tot
  real(kind=8),allocatable,dimension(:,:,:)::second_moments,second_moments_tot
  real(kind=8),allocatable,dimension(:)::min_dens_tot,av_dens_tot,phi_min_tot
  real(kind=8),allocatable,dimension(:)::max_dens_tot,e_kin_int_tot,e_bind_tot,e_thermal_tot
  real(kind=8),allocatable,dimension(:)::e_kin_int_tot4,e_bind_tot4,e_thermal_tot4
  real(kind=8),allocatable,dimension(:)::clump_mass_tot,clump_vol_tot
  real(kind=8),allocatable,dimension(:,:)::peak_pos_tot,bulk_momentum_tot
  real(kind=8),allocatable,dimension(:)::saddle_max_tot
  real(kind=8),allocatable,dimension(:)::relevance_tot
  real(kind=8),allocatable,dimension(:)::phi_ref, phi_ref_tot,phi_ref2_tot
  real(kind=8),allocatable,dimension(:)::Psurf,Psurf_tot,v_therm_tot,v_rms_tot,m4_tot
  real(kind=8),allocatable,dimension(:)::E_bind_iso_tot,E_therm_iso_tot,E_kin_iso_tot
  real(kind=8),allocatable,dimension(:)::peak_check,ball4_check,isodens_check,clump_check

  ! Test particles properties
  real(qdp),allocatable,dimension(:)::denp ! Density of the cell containing a test particle. Davide: used by the clump finder.
  integer,allocatable,dimension(:)::iglobalp,icellp,levp,testp_sort ! Used to sort test particles by density  
  integer,allocatable,dimension(:)::n_cells_tot,minmatch_tot,new_peak
  integer,allocatable,dimension(:)::sort_index
  integer,allocatable,dimension(:)::occupied,occupied_all ! Tells whether there is already a sink in a clump

end module clfind_commons
