subroutine read_hydro_params(nml_ok)
  use amr_commons
  use hydro_commons
  use mpi_mod
  implicit none
  logical::nml_ok
  !--------------------------------------------------
  ! Local variables
  !--------------------------------------------------
  integer::i,idim,ifixed,nboundary_true=0
  integer ,dimension(1:MAXBOUND)::bound_type
  real(dp)::ek_bound
  logical :: dummy
#ifdef SOLVERmhd
  real(dp)::em_bound
#endif

  !--------------------------------------------------
  ! Namelist definitions
  !--------------------------------------------------

  ! Initial conditions parameters
  namelist/init_params/filetype,initfile,multiple,nregion,region_type &
       & ,x_center,y_center,z_center,aexp_ini &
       & ,length_x,length_y,length_z,exp_region &
       & ,d_region,u_region,v_region,w_region,p_region,gamma_region &
       & ,pert_r, pert_dx, heating_r, heating_dx &
#ifdef SOLVERmhd
       & ,A_region,B_region,C_region,B_ave &
#if NVAR>8+NENER
       & ,var_region &
#endif
#else
#if NVAR>NDIM+2+NENER
       & ,var_region &
#endif
#endif
#if NENER>0
       & ,prad_region &
#endif
       & ,omega_b

  ! Hydro parameters
  namelist/hydro_params/gamma,courant_factor,smallr,smallc &
       & ,niter_riemann,slope_type,difmag &
#if NENER>0
       & ,gamma_rad &
#endif
#ifdef SOLVERmhd
       & ,riemann2d,slope_mag_type,eta_mag &
#endif
       & ,pressure_fix,beta_fix,scheme,riemann &
       & ,strict_equilibrium

  ! Refinement parameters
  namelist/refine_params/x_refine,y_refine,z_refine,r_refine &
       & ,a_refine,b_refine,exp_refine,jeans_refine,mass_cut_refine &
       & ,m_refine,mass_sph,err_grad_d,err_grad_p,err_grad_u &
       & ,floor_d,floor_u,floor_p,ivar_refine,var_cut_refine &
#ifdef SOLVERmhd
       & ,err_grad_A,err_grad_B,err_grad_C,err_grad_B2 &
       & ,floor_A,floor_B,floor_C,floor_B2,interpol_mag_type &
#endif
       & ,interpol_var,interpol_type,sink_refine

  ! Boundary parameters
  namelist/boundary_params/nboundary,bound_type &
       & ,ibound_min,ibound_max,jbound_min,jbound_max &
       & ,kbound_min,kbound_max &
#if NENER>0
       & ,prad_bound &
#endif
#ifdef SOLVERmhd
#if NVAR>8+NENER
       & ,var_bound &
#endif
       & ,A_bound,B_bound,C_bound &
#else
#if NVAR>NDIM+2+NENER
       & ,var_bound &
#endif
#endif
       & ,d_bound,u_bound,v_bound,w_bound,p_bound,no_inflow

  ! Feedback parameters
  namelist/feedback_params/eta_sn,eta_ssn,yield,rbubble,f_ek,ndebris &
       & ,f_w,mass_gmc,kappa_IR,delayed_cooling,momentum_feedback &
       & ,ir_feedback,ir_eff,t_diss,t_sne,mass_star_max,mass_sne_min

  ! Cooling / basic chemistry parameters
  namelist/cooling_params/cooling,metal,isothermal,haardt_madau,J21 &
       & ,a_spec,self_shielding, z_ave,z_reion,ind_rsink,T2max,neq_chem

  ! Star formation parameters
  namelist/sf_params/m_star,n_star,T2_star,g_star,del_star &
       & ,eps_star,jeans_ncells,sf_virial,sf_trelax,sf_tdiss,sf_model&
       & ,sf_log_properties,sf_imf,sf_compressive

  ! Units parameters
  namelist/units_params/units_density,units_time,units_length

  ! Dummy namelist for physics
  namelist/physics_params/ dummy

#ifdef grackle
   namelist/grackle_params/use_grackle,grackle_with_radiative_cooling,grackle_primordial_chemistry,grackle_metal_cooling &
       & ,grackle_UVbackground,grackle_cmb_temperature_floor,grackle_h2_on_dust,grackle_photoelectric_heating &
       & ,grackle_use_volumetric_heating_rate,grackle_use_specific_heating_rate,grackle_three_body_rate,grackle_cie_cooling &
       & ,grackle_h2_optical_depth_approximation,grackle_ih2co,grackle_ipiht,grackle_NumberOfTemperatureBins,grackle_CaseBRecombination &
       & ,grackle_Compton_xray_heating,grackle_LWbackground_sawtooth_suppression,grackle_NumberOfDustTemperatureBins,grackle_use_radiative_transfer &
       & ,grackle_radiative_transfer_coupled_rate_solver,grackle_radiative_transfer_intermediate_step,grackle_radiative_transfer_hydrogen_only &
       & ,grackle_self_shielding_method,grackle_Gamma,grackle_photoelectric_heating_rate,grackle_HydrogenFractionByMass &
       & ,grackle_DeuteriumToHydrogenRatio,grackle_SolarMetalFractionByMass,grackle_TemperatureStart,grackle_TemperatureEnd &
       & ,grackle_DustTemperatureStart,grackle_DustTemperatureEnd,grackle_LWbackground_intensity,grackle_UVbackground_redshift_on &
       & ,grackle_UVbackground_redshift_off,grackle_UVbackground_redshift_fullon,grackle_UVbackground_redshift_drop &
       & ,grackle_cloudy_electron_fraction_factor,grackle_data_file
#endif

  ! Read namelist file
  rewind(1)
  read(1,NML=init_params,END=121)
  goto 122
121 write(*,*)' You need to set up namelist &INIT_PARAMS in parameter file'
  call clean_stop
122 rewind(1)

  ! Fail if physics params is found
  read(1, NML=physics_params, end=110)
  if (myid == 1) &
       write(*, *) 'ERROR: the namelist contains the old `physics_params`. It has been depreciated in favor of the sections feedback_, cooling_, sf_ and units_params.'
  call clean_stop
110 rewind(1)
  if(nlevelmax>levelmin)read(1,NML=refine_params)
  rewind(1)
  if(hydro)read(1,NML=hydro_params)
  rewind(1)
  read(1,NML=boundary_params,END=103)
  simple_boundary=.true.
  goto 104
103 simple_boundary=.false.
104 if(nboundary>MAXBOUND)then
    write(*,*) 'Error: nboundary>MAXBOUND'
    call clean_stop
  end if
  rewind(1)
  read(1,NML=feedback_params,END=105)
105 continue
  rewind(1)
  read(1,NML=cooling_params,END=106)
106 continue
  rewind(1)
  read(1,NML=sf_params,END=107)
107 continue
  rewind(1)
  read(1,NML=units_params,END=108)
108 continue
#ifdef grackle
  rewind(1)
  read(1,NML=grackle_params)
#endif
#ifdef ATON
  if(aton)call read_radiation_params(1)
#endif

! MHD sets enums for riemann solvers, hydro just does a string compare
#ifdef SOLVERmhd
  !------------------------------------------------
  ! set ischeme
  !------------------------------------------------
  SELECT CASE (scheme)
  CASE ('muscl')
    ischeme = 0
  CASE ('induction')
    ischeme = 1

  CASE DEFAULT
    write(*,*)'unknown scheme'
    call clean_stop
  END SELECT
  !------------------------------------------------
  ! set iriemann
  !------------------------------------------------
  SELECT CASE (riemann)
  CASE ('llf')
    iriemann = 0
  CASE ('roe')
    iriemann = 1
  CASE ('hll')
    iriemann = 2
  CASE ('hlld')
    iriemann = 3
  CASE ('upwind')
    iriemann = 4
  CASE ('hydro')
    iriemann = 5

  CASE DEFAULT
    write(*,*)'unknown riemann solver'
    call clean_stop
  END SELECT
  !------------------------------------------------
  ! set iriemann
  !------------------------------------------------
  SELECT CASE (riemann2d)
  CASE ('llf')
    iriemann2d = 0
  CASE ('roe')
    iriemann2d = 1
  CASE ('upwind')
    iriemann2d = 2
  CASE ('hll')
    iriemann2d = 3
  CASE ('hlla')
    iriemann2d = 4
  CASE ('hlld')
    iriemann2d = 5
  CASE DEFAULT
    write(*,*)'unknown 2D riemann solver'
    call clean_stop
  END SELECT

  !--------------------------------------------------
  ! Make sure virtual boundaries are expanded to
  ! account for staggered mesh representation
  !--------------------------------------------------
  nexpand_bound=2
#endif

  !--------------------------------------------------
  ! Check for dm only cosmo run
  !--------------------------------------------------
  if(.not.hydro)then
     omega_b = 0.0D0
  endif

  !--------------------------------------------------
  ! Check for star formation
  !--------------------------------------------------
  if(eps_star>0)then
     ! For historical reference:
     ! t_star=0.1635449d0*(n_star/0.1d0)**(-0.5d0)/eps_star
     star=.true.
     pic=.true.
  endif

  !--------------------------------------------------
  ! Check for metal
  !--------------------------------------------------
#ifdef SOLVERmhd
  if(metal.and.nvar<(ndim+6))then
#else
  if(metal.and.nvar<(ndim+3))then
#endif
     if(myid==1)write(*,*)'Error: metals need nvar >= ndim+3'
     if(myid==1)write(*,*)'Modify hydro_parameters.f90 and recompile'
     nml_ok=.false.
  endif

  !--------------------------------------------------
  ! Check whether illegally trying non-eq chemistry
  !--------------------------------------------------
#ifndef RT
  if(neq_chem) then
     if(myid==1)write(*,*) 'Error: non-equilibrium chemistry unavailable'
     if(myid==1)write(*,*) 'Recompile with RT=True (or -DRT)'
     nml_ok=.false.
  endif
#endif

  !--------------------------------------------------
  ! Check for non-thermal energies
  !--------------------------------------------------
#if NENER>0
#ifdef SOLVERmhd
  if(nvar<(8+nener))then
#else
  if(nvar<(ndim+2+nener))then
#endif
     if(myid==1)write(*,*)'Error: non-thermal energy need nvar >= ndim+2+nener'
     if(myid==1)write(*,*)'Modify NENER and recompile'
     nml_ok=.false.
  endif
#endif

  !--------------------------------------------------
  ! Check ind_rsink
  !--------------------------------------------------
  if(ind_rsink<=0.0d0)then
     if(myid==1)write(*,*)'Error in the namelist'
     if(myid==1)write(*,*)'Check ind_rsink'
     nml_ok=.false.
  end if

  !-------------------------------------------------
  ! This section deals with hydro boundary conditions
  !-------------------------------------------------
  if(simple_boundary.and.nboundary==0)then
     simple_boundary=.false.
  endif

  if (simple_boundary)then

     ! Compute new coarse grid boundaries
     do i=1,nboundary
        if(ibound_min(i)*ibound_max(i)==1.and.ndim>0.and.bound_type(i)>0)then
           nx=nx+1
           if(ibound_min(i)==-1)then
              icoarse_min=icoarse_min+1
              icoarse_max=icoarse_max+1
           end if
           nboundary_true=nboundary_true+1
        end if
     end do
     do i=1,nboundary
        if(jbound_min(i)*jbound_max(i)==1.and.ndim>1.and.bound_type(i)>0)then
           ny=ny+1
           if(jbound_min(i)==-1)then
              jcoarse_min=jcoarse_min+1
              jcoarse_max=jcoarse_max+1
           end if
           nboundary_true=nboundary_true+1
        end if
     end do
     do i=1,nboundary
        if(kbound_min(i)*kbound_max(i)==1.and.ndim>2.and.bound_type(i)>0)then
           nz=nz+1
           if(kbound_min(i)==-1)then
              kcoarse_min=kcoarse_min+1
              kcoarse_max=kcoarse_max+1
           end if
           nboundary_true=nboundary_true+1
        end if
     end do

     ! Compute boundary geometry
     do i=1,nboundary
        if(ibound_min(i)*ibound_max(i)==1.and.ndim>0.and.bound_type(i)>0)then
           if(ibound_min(i)==-1)then
              ibound_min(i)=icoarse_min+ibound_min(i)
              ibound_max(i)=icoarse_min+ibound_max(i)
              if(bound_type(i)==1)boundary_type(i)=1
              if(bound_type(i)==2)boundary_type(i)=11
              if(bound_type(i)==3)boundary_type(i)=21
           else
              ibound_min(i)=icoarse_max+ibound_min(i)
              ibound_max(i)=icoarse_max+ibound_max(i)
              if(bound_type(i)==1)boundary_type(i)=2
              if(bound_type(i)==2)boundary_type(i)=12
              if(bound_type(i)==3)boundary_type(i)=22
           end if
           if(ndim>1)jbound_min(i)=jcoarse_min+jbound_min(i)
           if(ndim>1)jbound_max(i)=jcoarse_max+jbound_max(i)
           if(ndim>2)kbound_min(i)=kcoarse_min+kbound_min(i)
           if(ndim>2)kbound_max(i)=kcoarse_max+kbound_max(i)
        else if(jbound_min(i)*jbound_max(i)==1.and.ndim>1.and.bound_type(i)>0)then
           ibound_min(i)=icoarse_min+ibound_min(i)
           ibound_max(i)=icoarse_max+ibound_max(i)
           if(jbound_min(i)==-1)then
              jbound_min(i)=jcoarse_min+jbound_min(i)
              jbound_max(i)=jcoarse_min+jbound_max(i)
              if(bound_type(i)==1)boundary_type(i)=3
              if(bound_type(i)==2)boundary_type(i)=13
              if(bound_type(i)==3)boundary_type(i)=23
           else
              jbound_min(i)=jcoarse_max+jbound_min(i)
              jbound_max(i)=jcoarse_max+jbound_max(i)
              if(bound_type(i)==1)boundary_type(i)=4
              if(bound_type(i)==2)boundary_type(i)=14
              if(bound_type(i)==3)boundary_type(i)=24
           end if
           if(ndim>2)kbound_min(i)=kcoarse_min+kbound_min(i)
           if(ndim>2)kbound_max(i)=kcoarse_max+kbound_max(i)
        else if(kbound_min(i)*kbound_max(i)==1.and.ndim>2.and.bound_type(i)>0)then
           ibound_min(i)=icoarse_min+ibound_min(i)
           ibound_max(i)=icoarse_max+ibound_max(i)
           jbound_min(i)=jcoarse_min+jbound_min(i)
           jbound_max(i)=jcoarse_max+jbound_max(i)
           if(kbound_min(i)==-1)then
              kbound_min(i)=kcoarse_min+kbound_min(i)
              kbound_max(i)=kcoarse_min+kbound_max(i)
              if(bound_type(i)==1)boundary_type(i)=5
              if(bound_type(i)==2)boundary_type(i)=15
              if(bound_type(i)==3)boundary_type(i)=25
           else
              kbound_min(i)=kcoarse_max+kbound_min(i)
              kbound_max(i)=kcoarse_max+kbound_max(i)
              if(bound_type(i)==1)boundary_type(i)=6
              if(bound_type(i)==2)boundary_type(i)=16
              if(bound_type(i)==3)boundary_type(i)=26
           end if
        end if
     end do
     do i=1,nboundary
        ! Check for errors
        if( (ibound_min(i)<0.or.ibound_max(i)>(nx-1)) .and. (ndim>0) .and.bound_type(i)>0 )then
           if(myid==1)write(*,*)'Error in the namelist'
           if(myid==1)write(*,*)'Check boundary conditions along X direction',i
           nml_ok=.false.
        end if
        if( (jbound_min(i)<0.or.jbound_max(i)>(ny-1)) .and. (ndim>1) .and.bound_type(i)>0)then
           if(myid==1)write(*,*)'Error in the namelist'
           if(myid==1)write(*,*)'Check boundary conditions along Y direction',i
           nml_ok=.false.
        end if
        if( (kbound_min(i)<0.or.kbound_max(i)>(nz-1)) .and. (ndim>2) .and.bound_type(i)>0)then
           if(myid==1)write(*,*)'Error in the namelist'
           if(myid==1)write(*,*)'Check boundary conditions along Z direction',i
           nml_ok=.false.
        end if
     end do
  end if
  nboundary=nboundary_true
  if(simple_boundary.and.nboundary==0)then
     simple_boundary=.false.
  endif

  !--------------------------------------------------
  ! Compute boundary conservative variables
  !--------------------------------------------------
  do i=1,nboundary
     boundary_var(i,1)=MAX(d_bound(i),smallr)
     boundary_var(i,2)=d_bound(i)*u_bound(i)
#if NDIM>1 || SOLVERmhd
     boundary_var(i,3)=d_bound(i)*v_bound(i)
#endif
#if NDIM>2 || SOLVERmhd
     boundary_var(i,4)=d_bound(i)*w_bound(i)
#endif
     ek_bound=0.0d0
     do idim=1,ndim
        ek_bound=ek_bound+0.5d0*boundary_var(i,idim+1)**2/boundary_var(i,1)
     end do
     boundary_var(i,ndim+2)=ek_bound+P_bound(i)/(gamma-1.0d0)
#ifdef SOLVERmhd
     boundary_var(i,6)=A_bound(i)
     boundary_var(i,7)=B_bound(i)
     boundary_var(i,8)=C_bound(i)
     boundary_var(i,nvar+1)=A_bound(i)
     boundary_var(i,nvar+2)=B_bound(i)
     boundary_var(i,nvar+3)=C_bound(i)
     ek_bound=0.5d0*d_bound(i)*(u_bound(i)**2+v_bound(i)**2+w_bound(i)**2)
     em_bound=0.5d0*(A_bound(i)**2+B_bound(i)**2+C_bound(i)**2)
     boundary_var(i,5)=ek_bound+em_bound+P_bound(i)/(gamma-1.0d0)
#endif
  end do

  !-----------------------------------
  ! Rearrange level dependent arrays
  !-----------------------------------
  do i=nlevelmax,levelmin,-1
     jeans_refine(i)=jeans_refine(i-levelmin+1)
  end do
  do i=1,levelmin-1
     jeans_refine(i)=-1
  end do

  !-----------------------------------
  ! Sort out passive variable indices
  !-----------------------------------
#ifdef SOLVERmhd
  ! Hard-coded variables are rho,v*ndim,P,B*ndim
  ! MHD only works in 3D, so ndim=3
  ifixed=8
#else
  ! Hard-coded variables are rho,v*ndim,P
  ifixed=ndim+2
#endif
  inener=ifixed+1
  imetal=inener+nener
  idelay=imetal
  if(metal)idelay=imetal+1
  ivirial1=idelay
  ivirial2=idelay
  if(delayed_cooling)then
     ivirial1=idelay+1
     ivirial2=idelay+1
  endif
  if(sf_virial)then
     if(sf_compressive) ivirial2=ivirial1+1
  endif
  ixion=ivirial2
  if(sf_virial)ixion=ivirial2+1
  ichem=ixion
  if(aton)ichem=ixion+1
  if(myid==1.and.hydro.and.(nvar>ndim+2)) then
     write(*,'(A50)')"__________________________________________________"
     write(*,*) 'Hydro var indices:'
#if NENER>0
     write(*,*) '   inener   = ',inener
#endif
     if(metal)           write(*,*) '   imetal   = ',imetal
     if(delayed_cooling) write(*,*) '   idelay   = ',idelay
     if(sf_virial)then
        write(*,*) '   ivirial1 = ',ivirial1
        if(sf_compressive) write(*,*) '   ivirial2 = ',ivirial2
     endif
     if(aton)            write(*,*) '   ixion    = ',ixion
#ifdef RT
     if(rt) write(*,*) '   iIons    = ',ichem
#endif
     write(*,'(A50)')"__________________________________________________"
  endif

  ! Last variable is ichem

#ifdef SOLVERmhd
  !-----------------------------------
  ! Set magnetic slope limiters
  !-----------------------------------
  if (slope_mag_type == -1) then
    slope_mag_type = slope_type
  endif
  if (interpol_mag_type == -1) then
    interpol_mag_type = interpol_type
  endif
#endif

end subroutine read_hydro_params
