subroutine read_fld_params(namelist_unit,nml_ok)
   use amr_commons, only:myid
   use hydro_parameters
   use radiation_parameters
   use constants
   use units_commons
   use cloud_module
#if RT==1
   use rt_parameters,only:rt_protostar_m1
#endif
   implicit none
   integer,intent(in)::namelist_unit
   logical,intent(inout)::nml_ok
   integer::nml_err
   !------------------------------------------------------------------------
   ! Reads the parameters for the flux-limited diffusion (FLD) radiation.
   !------------------------------------------------------------------------
#if USE_FLD==1
   integer::irad,i,j
   real(dp)::radiation_source
   character(len=2):: rad_trans_model='m1'
#endif

   ! local variables for unit conversion
   real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v

   namelist/radiation_params/grey_rad_transfer,dtdiff_params,dt_control &
        & ,rosseland_params,planck_params,epsilon_diff,fld_limiter &
        & ,freqs_in_Hz,read_groups,split_groups_log,extra_end_group  &
        & ,numin,numax,Tr_floor,robin,rad_trans_model,min_optical_depth,rt_feedback,Tray_min &
        & ,PMS_evol,Hosokawa_track,energy_fix,facc_star,facc_star_lum,valp_min,store_matrix &
        & ,rt_protostar_fld,sublimation_kuiper,lum_injection &
        & ,sinks_opt_thin, fit_semenov

   rewind(1)
   if(FLD)read(1,NML=radiation_params)

  ! Conversion factor from user units to cgs units (to be done after read physics_params with units_density...)
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

#if USE_FLD==1 || USE_M_1==1
  ! Initialize multigroup
  allocate(nu_min_hz(1:ngrp),nu_max_hz(1:ngrp),nu_min_ev(1:ngrp),nu_max_ev(1:ngrp))
  call create_groups
  eray_min = (aR)*Tray_min**4
  deray_min = (4.0d0*aR)*Tray_min**3
  small_er = eray_min/(scale_d*Scale_v**2)
  call tabulate_art4
  call read_omegas
  if(myid==1 .and. grey_rad_transfer .and. ngrp .gt.1) then
     print*,'Warning: Grey Radiation Transfer with NRAD>1'
     call clean_stop
  endif
  scale_E0 = aR*(Tr_floor**4)
  P_cal = scale_E0 / (scale_d * scale_v**2)
  C_cal = c_cgs / scale_v
  is_radiative_energy = .false.
#endif

#if USE_FLD==1
  ! Set i_fld_limiter
  i_fld_limiter=i_fld_limiter_nolim
  if(fld_limiter=='levermore') i_fld_limiter=i_fld_limiter_levermore
  if(fld_limiter=='minerbo')  i_fld_limiter=i_fld_limiter_minerbo
  ! Index array for radiative variables and temperature
  ! Needed in M1 because temperature is stored in uold(:,nvar)
  do irad = 1,nvar_bicg
     ind_bicg (irad) = firstindex_er+irad
     norm_bicg(irad) = P_cal
  enddo
  ind_trad(1) = nvar
  norm_trad(1) = Tr_floor
  do irad = 2,nvar_trad
     ind_trad(irad) = firstindex_er-1+irad
     norm_trad(irad) = P_cal
     is_radiative_energy(irad) = .true.
  enddo

#if RT==1
  if(rt_protostar_fld .and. rt_protostar_m1) then
     write(*,*)'Wrong choice for rt_protostar method : choose one kind not both'
     call clean_stop
  end if
#endif
#endif

#if USE_M_1==1
  ! Set radiative transfer model
  select case(rad_trans_model)
  case('P1','p1')
     irad_trans_model = irad_trans_model_p1
  case('M1','m1')
     irad_trans_model = irad_trans_model_m1
  case default
     if(myid==1) write(*,*) 'unknown radiative transfer model: '//rad_trans_model
     call clean_stop
  end select
  call compute_valp
  ! Index array for radiative variables and temperature
  ! Needed in M1 because temperature is stored in uold(:,nvar)
  ind_bicg(1) = nvar
  norm_bicg(1) = Tr_floor
  do irad = 2,nvar_bicg
     ind_bicg(irad) = firstindex_er-1+irad
     norm_bicg(irad) = P_cal
  enddo
  do irad = ngrp+2,nvar_bicg
     norm_bicg(irad) = norm_bicg(irad)*C_cal
  enddo
  ind_trad=ind_bicg
  norm_trad=norm_bicg
  is_radiative_energy(2:ngrp+1) = .true.
#endif

#if USE_FLD==1
  ! Multigroup opacities initialization
  if(fld)then
     if((opacity_type == 'grey') .and. (ngrp > 1) .and. (.not.stellar_photon))then
        if(myid == 1)then
           write(*,*) 'WARNING! Trying to use grey opacity table with ngrp =',ngrp
           write(*,*) 'Switching to multigroup opacities'
        endif
        opacity_type = 'multigroup'
     endif
     call init_opacities
  end if
 
  if(PMS_evol .and. rt_feedback .and. Hosokawa_track)then
     open(101,file='Hosokawa_track.dat', status='old')
     read(101,*)nmdot_PMS,nm_PMS,ndata_PMS
     allocate(nb_ligne_PMS(nmdot_PMS))
     allocate(data_PMS(nmdot_PMS,nm_PMS,ndata_PMS))
     read(101,*)nb_ligne_PMS(1),nb_ligne_PMS(2),nb_ligne_PMS(3),nb_ligne_PMS(4),nb_ligne_PMS(5)
     do i=1,nmdot_PMS
        do j=1,nm_PMS
           read(101,*)data_PMS(i,j,1),data_PMS(i,j,2),data_PMS(i,j,3),data_PMS(i,j,4) ! mdot,mass,luminosity,radius
        end do
     end do
     close(101)
  end if
#endif

end subroutine read_fld_params