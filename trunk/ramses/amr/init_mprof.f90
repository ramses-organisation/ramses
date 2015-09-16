subroutine init_mprof
  use amr_commons
  implicit none
#ifdef MPROFILING
  acc_t_make_virtual_fine_dp=0.0d0
  acc_t_make_virtual_reverse_dp=0.0d0
  acc_t_set_unew=0.0d0
  acc_t_set_uold=0.0d0
  acc_t_add_gravity_source_terms=0.0d0
  acc_t_synchro_hydro_fine=0.0d0
  acc_t_make_boundary_phi=0.0d0
  acc_t_make_boundary_mask=0.0d0
  acc_t_force_fine=0.0d0
  acc_t_save_phi_old=0.0d0
  acc_t_restrict_mask_coarse_reverse=0.0d0
  acc_t_cmp_residual_mg_coarse=0.0d0
  acc_t_gauss_seidel_mg_coarse=0.0d0
  acc_t_restrict_residual_coarse_reverse=0.0d0
  acc_t_interpolate_and_correct_coarse=0.0d0
  acc_t_set_scan_flag_coarse=0.0d0
  acc_t_make_fine_mask=0.0d0
  acc_t_make_fine_bc_rhs=0.0d0
  acc_t_make_virtual_mg_dp=0.0d0
  acc_t_make_reverse_mg_dp=0.0d0
  acc_t_restrict_mask_fine_reverse=0.0d0
  acc_t_cmp_residual_mg_fine=0.0d0
  acc_t_gauss_seidel_mg_fine=0.0d0
  acc_t_restrict_residual_fine_reverse=0.0d0
  acc_t_interpolate_and_correct_fine=0.0d0
  acc_t_set_scan_flag_fine=0.0d0
  acc_t_make_initial_phi=0.0d0
  acc_t_make_multipole_phi=0.0d0
  acc_t_amr_step=0.0d0
  acc_t_godunov_fine=0.0d0
  acc_t_add_pdv_source_terms=0.0d0
  acc_t_multigrid_fine=0.0d0
  acc_t_cmp_residual_norm2_fine=0.0d0
  acc_t_total_time=0.0d0
#endif 
end subroutine
