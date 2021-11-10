SUBROUTINE read_sink_feedback_params(nml_ok)

  ! Read FEEDBACK_PARAMS namelist
  !-------------------------------------------------------------------------
    use amr_parameters,only:dp
    use amr_commons
    use sink_feedback_module
    implicit none
    logical::nml_ok
    real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  !-------------------------------------------------------------------------
    namelist/feed_params/  Tsat, Vsat, sn_r_sat, sn_p, sn_e, sn_mass, &
          & FB_nsource, FB_on, FB_start, FB_end, FB_sourcetype, &
          & FB_pos_x, FB_pos_y, FB_pos_z, &
          & FB_mejecta, FB_energy, FB_thermal, &
          & FB_radius, Vdisp
    rewind(1)
    read(1,NML=feed_params,END=101)
  101 continue
  
        call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  
        !normalise the supernova quantities
        sn_p_ref = sn_p / (scale_d * scale_v * scale_l**3)
        sn_e_ref = sn_e / (scale_d * scale_v**2 * scale_l**3)
        sn_mass_ref = sn_mass / (scale_d * scale_l**3) !1 solar mass ejected
  
        !normalise Vsat which is assumed to be in KM/S
        Vsat = Vsat * 1.e5 / scale_v
  
        !normalise Vdisp which is assumed to be in KM/S
        Vdisp = Vdisp * 1.e5 / scale_v
    
END SUBROUTINE read_sink_feedback_params
