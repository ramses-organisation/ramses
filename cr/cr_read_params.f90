!################################################################
!################################################################
subroutine read_cr_params(nml_ok)
  ! Read the &CR_PARAMS namelist block (cosmic-ray runtime parameters)
  ! from the parameter file already open on unit 1, and sanity-check it.
  use amr_commons
  use cr_parameters
  implicit none
  logical::nml_ok
  integer::ig

  namelist/cr_params/cr_advect,cr_HLLE,cr_isotropic_pressure &
       & ,cr_flux_correction,gamma_cr,cr_feedback &
       & ,cr_smallr_decouple,cr_efloor,cr_c_fraction,cr_nsubcycle &
       & ,cr_varvmax,cr_varvmax_fudge,cr_varvmax_vdvs &
       & ,Dcr,DCRmax,Dcr_perp_factor,cr_streaming_diffusion &
       & ,cr_streaming_heating,cr_v_alfven,cr_f_taucell &
       & ,cr_cooling,zeta_cr,cr_ne,cr_fneut &
       & ,cr_condinit_kind,cr_boundana_kind,err_grad_cr &
       & ,cr_region_u,cr_boundary_u &
       & ,cr_nregion,cr_region_type,cr_reg_x_center,cr_reg_y_center,cr_reg_z_center &
       & ,cr_reg_length_x,cr_reg_length_y,cr_reg_length_z,cr_exp_region,cr_reg_group

  if(myid==1)write(*,'(A50)')"Reading cr_params namelist ..."
  rewind(1)
  read(1,NML=cr_params,END=101)
101 continue                                 ! No harm if no cr namelist

  ! Energy index of each CR group in cruold/crunew (mirror rt_init iGroups)
  do ig=1,ncr_groups
     Ecr_idx(ig)=1+(ndim+1)*(ig-1)
  end do

  ! Sanity checks
  if(cr_nsubcycle<1)then
     if(myid==1)write(*,*)'Error in cr_params: cr_nsubcycle must be >= 1'
     nml_ok=.false.
  endif
  if(cr_c_fraction<=0d0 .or. cr_c_fraction>1d0)then
     if(myid==1)write(*,*)'Error in cr_params: cr_c_fraction must be in (0,1]'
     nml_ok=.false.
  endif

  if(.not.(trim(cr_condinit_kind)==''         .or. &
       &   trim(cr_condinit_kind)=='jiang_411'.or. &
       &   trim(cr_condinit_kind)=='jiang_412'.or. &
       &   trim(cr_condinit_kind)=='jiang_413'.or. &
       &   trim(cr_condinit_kind)=='jiang_414'.or. &
       &   trim(cr_condinit_kind)=='jiang_415'.or. &
       &   trim(cr_condinit_kind)=='jiang_415_donut'.or. &
       &   trim(cr_condinit_kind)=='jiang_411_triangular'.or. &
       &   trim(cr_condinit_kind)=='jiang_421'.or. &
       &   trim(cr_condinit_kind)=='jiang_422'.or. &
       &   trim(cr_condinit_kind)=='jiang_424'))then
     if(myid==1)write(*,*)'Error in cr_params: unknown cr_condinit_kind="'// &
          & trim(cr_condinit_kind)//'". Valid values: '// &
          & '"", jiang_411, jiang_412, jiang_413, jiang_414, jiang_415, '// &
          & 'jiang_415_donut, jiang_411_triangular, jiang_421, jiang_422, jiang_424'
     nml_ok=.false.
  endif

  if(.not.(trim(cr_boundana_kind)==''                    .or. &
       &   trim(cr_boundana_kind)=='jiang_411_triangular'.or. &
       &   trim(cr_boundana_kind)=='jiang_413'           .or. &
       &   trim(cr_boundana_kind)=='jiang_424'))then
     if(myid==1)write(*,*)'Error in cr_params: unknown cr_boundana_kind="'// &
          & trim(cr_boundana_kind)//'". Valid values: '// &
          & '"", jiang_411_triangular, jiang_413, jiang_424'
     nml_ok=.false.
  endif

end subroutine read_cr_params
