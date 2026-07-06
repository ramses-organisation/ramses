module dice_commons
  use amr_commons,ONLY:dp
  use hydro_commons,ONLY:mass_sph
  use constants,only:mu_mol

  ! particle data
  character(len=512)::ic_file, ic_format
  ! misc
  real(dp)::IG_rho         = 1.0D-5
  real(dp)::IG_T2          = 1.0D7
  real(dp)::IG_metal       = 0.01
  real(dp)::ic_scale_pos   = 1.0
  real(dp)::ic_scale_vel   = 1.0
  real(dp)::ic_scale_mass  = 1.0
  real(dp)::ic_scale_u     = 1.0
  real(dp)::ic_scale_age   = 1.0
  real(dp)::ic_scale_metal = 1.0
  real(dp)::ic_t_restart   = 0.0D0
  integer::ic_mask_ivar    = 0
  real(dp)::ic_mask_min    = 1d40
  real(dp)::ic_mask_max    = -1d40
  integer::ic_mask_ptype   = -1
  integer::ic_ifout        = 1
  integer::ic_nfile        = 1
  integer,dimension(1:6)::ic_skip_type        = -1
  integer,dimension(1:6)::cosmo_add_gas_index = -1
  real(dp),dimension(1:3)::ic_mag_const = (/ 0.0, 0.0, 0.0 /)
  real(dp),dimension(1:3)::ic_center    = (/ 0.0, 0.0, 0.0 /)
  character(len=4)::ic_head_name  = 'HEAD'
  character(len=4)::ic_pos_name   = 'POS '
  character(len=4)::ic_vel_name   = 'VEL '
  character(len=4)::ic_id_name    = 'ID  '
  character(len=4)::ic_mass_name  = 'MASS'
  character(len=4)::ic_u_name     = 'U   '
  character(len=4)::ic_metal_name = 'Z   '
  character(len=4)::ic_age_name   = 'AGE '
  ! Gadget units in cgs
  real(dp)::gadget_scale_l = 3.085677581282D21
  real(dp)::gadget_scale_v = 1.0D5
  real(dp)::gadget_scale_m = 1.9891D43
  real(dp)::gadget_scale_t = 1.0D6*365*24*3600
  real(dp),allocatable,dimension(:)::up
  real(dp),allocatable,dimension(:)::maskp
  logical::dice_init       = .false.
  logical::amr_struct      = .false.
#ifdef SOLVERmhd
  ! magnetic
  integer,parameter::MAXGAL= 32
  real(dp),dimension(1:MAXGAL)::ic_mag_center_x = 0.0
  real(dp),dimension(1:MAXGAL)::ic_mag_center_y = 0.0
  real(dp),dimension(1:MAXGAL)::ic_mag_center_z = 0.0
  real(dp),dimension(1:MAXGAL)::ic_mag_axis_x   = 0.0
  real(dp),dimension(1:MAXGAL)::ic_mag_axis_y   = 0.0
  real(dp),dimension(1:MAXGAL)::ic_mag_axis_z   = 1.0
  real(dp),dimension(1:MAXGAL)::ic_mag_scale_R  = 1.0
  real(dp),dimension(1:MAXGAL)::ic_mag_scale_H  = 1.0
  real(dp),dimension(1:MAXGAL)::ic_mag_scale_B  = 0.0
#endif

end module dice_commons
