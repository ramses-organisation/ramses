!*************************************************************************
SUBROUTINE rt_init

!  Initialize everything for radiative transfer
!  Some initialisation is also needed in case of non-equilibrium
!  chemistry, even if rt=.false.
!-------------------------------------------------------------------------
  use amr_commons
  use hydro_commons
  use rt_hydro_commons
  use rt_flux_module
  use rt_cooling_module,only:rt_isIRtrap,iIRtrapVar
  use rt_parameters
  use SED_module
  use UV_module
  implicit none
  integer:: i, nvar_count
!-------------------------------------------------------------------------
  if(verbose)write(*,*)'Entering init_rt'
  ! Count the number of variables and check if ok:
  nvar_count = ichem-1     ! # of non-rt vars: rho u v w p (z) (delay) (x)
  if(rt_isIRtrap) &
     iIRtrapVar = inener  ! Trapped rad. stored in nonthermal pressure var
  if(heat_unresolved_HII .eq. 2) then
     ! Using NENER for unresolved HII region heating
     iHIIheat=inener
     if(rt_isIRtrap) iHIIheat=inener+1
     if(nener.lt.iHIIheat-inener+1) then
        if(myid==1) then
           write(*,*) 'Need more NENER FOR HEATING HII REGIONS'
           write(*,*) 'STOPPING'
        endif
        call clean_stop
     endif
  endif
  iIons=nvar_count+1         !      Starting index of ionisation fractions
  nvar_count = nvar_count+NIONS  !                            # hydro vars

  if(nvar_count .gt. nvar) then
     if(myid==1) then
        write(*,*) 'rt_init(): Something wrong with NVAR.'
        write(*,*) 'Should have NVAR=2+ndim+1*metal+1*dcool+1*aton+IRtrap+nIons'
        write(*,*) 'Have NVAR=',nvar
        write(*,*) 'Should have NVAR=',nvar_count
        write(*,*) 'STOPPING!'
     endif
     call clean_stop
  endif

  if(rt_star .or. sedprops_update .ge. 0) &
     call init_SED_table    ! init stellar energy distribution properties
  if(rt .and. .not. hydro) then
     if(myid==1) then
        write(*,*) 'hydro must be turned on when running radiative transfer.'
        write(*,*) 'STOPPING!'
     endif
     call clean_stop
  endif
  if(rt_star) use_proper_time=.true.    ! Need proper birth time for stars
  if(rt) neq_chem=.true.        ! Equilibrium cooling doesn't work with RT

  ! To maximize efficiency, rt advection and rt timestepping is turned off
  ! until needed.
  if(rt .and. .not.rt_otsa) rt_advect=.true.
  if(rt .and. rt_nsource .gt. 0) rt_advect=.true.
  if(rt .and. rt_nregion .gt. 0) rt_advect=.true.
  if(rt .and. rt_AGN ) rt_advect=.true.
  ! UV propagation is checked in set_model
  ! Star feedback is checked in amr_step

  do i=1,nGroups  ! Starting indices in uold and unew of each photon group
     iGroups(i)=1+(ndim+1)*(i-1)
     if(nrestart.eq.0) then
        rtuold(:,iGroups(i))=smallNp
     endif
  end do
  if(trim(rt_flux_scheme).eq.'hll') rt_use_hll=.true.
  if(rt_use_hll) call read_hll_eigenvalues
  tot_cool_loopcnt=0 ; max_cool_loopcnt=0 ; n_cool_cells=0
  loopCodes=0
  tot_nPhot=0d0 ;  step_nPhot=0d0; step_nStar=0d0; step_mStar=0d0

END SUBROUTINE rt_init

!*************************************************************************
SUBROUTINE update_rt_c

! Update the speed of light for radiative transfer, in code units.
! This cannot be just a constant, since scale_v changes with time in
! cosmological simulations.
!-------------------------------------------------------------------------
  use rt_parameters
  use amr_commons
  implicit none
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
!-------------------------------------------------------------------------
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  rt_c=rt_c_cgs/scale_v
  rt_c2=rt_c**2
END SUBROUTINE update_rt_c

!*************************************************************************
SUBROUTINE adaptive_rt_c_update(ilevel, dt)

! Set the lightspeed such that RT can be done at ilevel in time dt in
! a single step.
!-------------------------------------------------------------------------
  use amr_parameters
  use rt_parameters
  use SED_module
  implicit none
  integer:: ilevel, nx_loc
  real(dp):: dt, scale, dx
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
!-------------------------------------------------------------------------
  ! Mesh spacing at ilevel
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5D0**ilevel*scale

  ! new lightspeed
  rt_c = dx/3d0/dt * rt_courant_factor
  rt_c2 = rt_c**2

  ! new ligtspeed in cgs
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  rt_c_cgs = rt_c*scale_v
  rt_c_fraction = rt_c_cgs/c_cgs

  call updateRTGroups_CoolConstants        ! These change as a consequence

END SUBROUTINE adaptive_rt_c_update


!*************************************************************************
SUBROUTINE read_rt_params(nml_ok)

! Read rt_params namelist
!-------------------------------------------------------------------------
  use amr_commons
  use rt_parameters
  use cooling_module, only:X, Y
  use rt_cooling_module
  use UV_module
  use SED_module
  implicit none
  logical::nml_ok
  integer::iCount
!-------------------------------------------------------------------------
  namelist/rt_params/rt_star, rt_esc_frac, rt_flux_scheme, rt_smooth     &
       & ,rt_is_outflow_bound, rt_TConst, rt_courant_factor              &
       & ,rt_c_fraction, rt_nsubcycle, rt_otsa, sedprops_update          &
       & ,sed_dir, uv_file, rt_UVsrc_nHmax, nUVgroups, nSEDgroups        &
       & ,SED_isEgy, rt_output_coolstats, hll_evals_file                 &
       & ,upload_equilibrium_x, X, Y, rt_is_init_xion                    &
       & ,rt_err_grad_n, rt_floor_n, rt_err_grad_xHII, rt_floor_xHII     &
       & ,rt_err_grad_xHI, rt_floor_xHI, rt_refine_aexp, is_mu_H2,isHe   &
       & ,isH2, rt_isIR, is_kIR_T, rt_T_rad, rt_vc, rt_pressBoost        &
       & ,rt_isoPress, rt_isIRtrap, iPEH_group, heat_unresolved_HII      &
       & ,cosmic_rays                                                    &
       ! RT regions (for initialization)                                 &
       & ,rt_nregion, rt_region_type                                     &
       & ,rt_reg_x_center, rt_reg_y_center, rt_reg_z_center              &
       & ,rt_reg_length_x, rt_reg_length_y, rt_reg_length_z              &
       & ,rt_exp_region, rt_reg_group                                    &
       & ,rt_n_region, rt_u_region, rt_v_region, rt_w_region             &
       ! RT source regions (for every timestep)                          &
       & ,rt_nsource, rt_source_type                                     &
       & ,rt_src_x_center, rt_src_y_center, rt_src_z_center              &
       & ,rt_src_length_x, rt_src_length_y, rt_src_length_z              &
       & ,rt_exp_source, rt_src_group                                    &
       & ,rt_n_source, rt_u_source, rt_v_source, rt_w_source             &
       ! RT boundary (for boundary conditions)                           &
       & ,rt_n_bound,rt_u_bound,rt_v_bound,rt_w_bound                    &
       & ,rt_AGN, rt_sink


  ! Set default initialisation of ionisation states:
  ! -Off if restarting, but can set to true (for postprocessing)
  ! -On otherwise, but can set to false (for specificic initialisation)
  if(nrestart .gt. 0) then
     rt_is_init_xion=.false.
  else
     rt_is_init_xion=.true.
  endif

  ! Read namelist file
  rewind(1)
  read(1,NML=rt_params,END=101)
101 continue                                   ! No harm if no rt namelist

  if(nGroups.le.0) rt=.false. ! No sense  doing rt if there are no photons
  if(.not. rt .and. .not. rt_star) sedprops_update=-1

  if(rt_err_grad_n .gt. 0. .or. rt_err_grad_xHII .gt. 0.                 &
       .or. rt_err_grad_xHI .gt. 0.) rt_refine=.true.

  rt_c_cgs = c_cgs * rt_c_fraction
  !call update_rt_c

  ! Trapped IR pressure closure as in Rosdahl & Teyssier 2015, eq 43:
  if(rt_isIRtrap) gamma_rad(1) = rt_c_fraction / 3d0 + 1d0

  if(rt_Tconst .ge. 0d0) rt_isTconst=.true.

  ! Set number of used ionisation fractions, indexes of ionization
  ! fractions, and ionization energies, and check if we have enough
  ! ionization variables (NIONS)
  if(rt .or. neq_chem) then
     iCount=0
     if(isH2) then
        iCount=iCount+1
        ixHI=iCount    ; ionEvs(ixHI)=ionEv_HI
     endif
     iCount=iCount+1    ; ixHII=iCount   ; ionEvs(ixHII)=ionEv_HII
     if(isHe) then
        iCount=iCount+1 ; ixHeII=iCount  ; ionEvs(ixHeII)=ionEv_HeII
        iCount=iCount+1 ; ixHeIII=iCount ; ionEvs(ixHeIII)=ionEv_HeIII
     endif
     nIonsUsed = iCount
     if(nIonsUsed .gt. NIONS) then
        if(myid==1) then
           write(*,*) 'Not enough variables for ionization fractions'
           write(*,*) 'Have NIONS=',NIONS
           write(*,*) 'STOPPING!'
        endif
        nml_ok=.false.
     endif
     if(nIonsUsed .lt. NIONS) then
        if(myid==1) then
           write(*,*) 'Too many variables for ionization fractions'
           write(*,*) 'Have NIONS=',NIONS
           write(*,*) 'Need NIONS=',nIonsUsed
           write(*,*) 'Probably no harm, so still continuing...'
        endif
     endif
     if(myid==1) then
        write(*,*) 'Number of ionization fractions is:',nIonsUsed
        write(*,*) 'The indexes are iHI, iHII, iHeII, iHeIII ='              &
             , ixHI, ixHII, ixHeII, ixHeIII
     endif
  endif

  if(rt_sink.and.(.not.stellar))then
     write(*,*) 'Enable stellar particles to use rt_sink'
     nml_ok=.false.
  endif

  call read_rt_groups(nml_ok)
END SUBROUTINE read_rt_params

!*************************************************************************
SUBROUTINE read_rt_groups(nml_ok)

! Read rt_groups namelist
!-------------------------------------------------------------------------
  use amr_commons
  use rt_parameters
  use rt_cooling_module
  use SED_module
  implicit none
  logical::nml_ok
  integer::i,igroup_HI=0, igroup_HII=0, igroup_HeII=0, igroup_HeIII=0
!-------------------------------------------------------------------------
  namelist/rt_groups/group_csn, group_cse, group_egy, spec2group         &
       & , groupL0, groupL1, kappaAbs, kappaSc, group_egy_AGNfrac
  if(myid==1) then
     write(*,'(" Working with ",I2," photon groups and  "                &
          & ,I2, " ion species")') nGroups, nIons
     write(*,*) ''
  endif

  if(nGroups .le. 0) then
     rt = .false.
     return
  endif
#if NGROUPS>0
  !  Use H2, HI, HeI, HeII ionization energies  as default group intervals
  groupL0(1:min(nGroups,nIons))=ionEvs(1:min(nGroups,nIons))! Lower bounds
  groupL1(1:min(nGroups,nIons-1))=ionEvs(2:min(nGroups+1,nIons)) !   Upper
  groupL1(min(nGroups,nIons))=0.                    ! Upper bound=infinity

  i=0
  if(isH2) then ! Set index for H2 dissociating group
     i=i+1 ; igroup_HI=i
  endif
  if(i .lt. nGroups) then ! Set index for HI ionizing group
     i=i+1 ; igroup_HII=i
  endif
  if(i .lt. nGroups .and. isHe) then ! Set index for HeI ionizing group
     i=i+1 ; igroup_HeII=i
  endif
  if(i .lt. nGroups .and. isHe) then ! Set index for HeII ionizing group
     i=i+1 ; igroup_HeIII=i
  endif

  ! Default groups are all blackbodies at 10^5 Kelvin:
  group_csn=0d0 ; group_cse=0d0 ; group_egy=0d0         ! Default all zero
  if(igroup_HI .gt. 0) then
     if(ixHI .gt. 0) then                                ! H2 dissociation
        group_csn(igroup_HI,ixHI)=2.1d-19
        group_cse(igroup_HI,ixHI)=2.1d-19
     endif
     group_egy(igroup_HI)=12.44
  endif
  if(igroup_HII .gt. 0) then
     if(ixHI .gt. 0) then                   ! H2 ionization by HI photons
        group_csn(igroup_HII,ixHI)=5.0d-18
        group_cse(igroup_HII,ixHI)=5.3d-18
     endif
     group_csn(igroup_HII,ixHII)=3.007d-18                ! HI ionization
     group_cse(igroup_HII,ixHII)=2.781d-18
     group_egy(igroup_HII)=18.85
  endif
  if(igroup_HeII .gt. 0) then
     if(ixHI .gt. 0) then                  ! H2 ionization by HeI photons
        group_csn(igroup_HeII,ixHI)=2.9d-18
        group_cse(igroup_HeII,ixHI)=2.8d-18
     endif
     group_csn(igroup_HeII,ixHII)=5.687d-19! HI ionization by HeI photons
     group_cse(igroup_HeII,ixHII)=5.042d-19
     if(ixHeII .gt. 0) then                              ! HeI ionization
        group_csn(igroup_HeII,ixHeII)=4.478d-18
        group_cse(igroup_HeII,ixHeII)=4.130d-18
     endif
     group_egy(igroup_HeII)=35.079
  endif
  if(igroup_HeIII .gt. 0) then
     if(ixHI .gt. 0) then                 ! H2 ionization by HeII photons
        group_csn(igroup_HeIII,ixHI)=4.1d-19
        group_cse(igroup_HeIII,ixHI)=4.1d-19
     endif
     group_csn(igroup_HeIII,ixHII)=7.889d-20  ! HI ioniz. by HeII photons
     group_cse(igroup_HeIII,ixHII)=7.456d-20
     if(ixHeII .gt. 0) then                  ! HeI ioniz. by HeII photons
        group_csn(igroup_HeIII,ixHeII)=1.197d-18
        group_cse(igroup_HeIII,ixHeII)=1.142d-18
     endif
     if(ixHeIII .gt. 0) then                            ! HeII ionization
        group_csn(igroup_HeIII,ixHeIII)=1.055d-18
        group_cse(igroup_HeIII,ixHeIII)=1.001d-18
     endif
     group_egy(igroup_HeIII)=65.666
  endif
#endif

  do i=1,min(nIons,nGroups)
     spec2group(i)=i                   ! Species contributions to groups
  end do

  ! Read namelist file
  rewind(1)
  read(1,NML=rt_groups,END=101)
101 continue              ! no harm if no rt namelist

  if(minval(group_egy) .le. 0d0 .and. myid==1) then
     print*,'========================================================='
     print*,'WARNING! Some photon groups have zero or negative energy!'
     print*,'This could have unwanted effects, so be careful!!!'
     print*,'========================================================='
  endif

  if(isH2) then
     do i=1,nGroups
        if((groupL0(i) .ge. 11.2) .and. (groupL1(i) .le. 13.6)          &
           .and. (groupL0(i) .le. 13.6) .and. (groupL1(i) .ge. 11.2))then
           ssh2(i) = 4d2 ! H2 self-shielding factor
           isLW(i) = 1d0 ! Index for LW groups
        endif
    enddo
  endif

  call updateRTGroups_CoolConstants
  call write_group_props(.false.,6)
END SUBROUTINE read_rt_groups

!************************************************************************
SUBROUTINE add_rt_sources(ilevel,dt)

! Inject radiation from RT source regions (from the RT namelist). Since
! the light sources are continuously emitting radiation, this is called
! continuously during code execution, rather than just during
! initialization.
!
! ilevel => amr level at which to inject the radiation
! dt     => timestep for injection (since injected values are per time)
!------------------------------------------------------------------------
  use amr_commons
  use rt_parameters
  use rt_hydro_commons
  implicit none
  integer::ilevel
  real(dp)::dt
  integer::i,igrid,ncache,iskip,ngrid
  integer::ind,idim,ivar,ix,iy,iz,nx_loc
  integer ,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp)::dx,dx_loc,scale
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp),dimension(1:nvector,1:nrtvar),save::uu
!------------------------------------------------------------------------
  call add_UV_background(ilevel)
  if(numbtot(1,ilevel)==0)return    ! no grids at this level
  if(rt_nsource .le. 0) return      ! no rt sources
  if(verbose)write(*,111)ilevel

  ! Mesh size at level ilevel in coarse cell units
  dx=0.5D0**ilevel
  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  ! Local constants
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  ncache=active(ilevel)%ngrid
  ! dx (and dx_loc=dx) are just equal to 1/nx (where 1 is the boxlength)
  ! Loop over grids by vector sweeps
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     ! Loop over cells
     do ind=1,twotondim
        ! Gather cell indices
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        ! Gather cell centre positions
        do idim=1,ndim
           do i=1,ngrid
              xx(i,idim)=xg(ind_grid(i),idim)+xc(ind,idim)
           end do
        end do
        ! Rescale position from code units to user units
        do idim=1,ndim
           do i=1,ngrid
              xx(i,idim)=(xx(i,idim)-skip_loc(idim))*scale
           end do
        end do
        ! Read the RT variables
        do ivar=1,nrtvar
           do i=1,ngrid
              uu(i,ivar)=rtunew(ind_cell(i),ivar)
           end do
        end do
        ! find injected values per cell
        call rt_sources_vsweep(xx,uu,dx_loc,dt,ngrid)
        ! Write the RT variables
        do ivar=1,nrtvar
           do i=1,ngrid
              rtunew(ind_cell(i),ivar)=uu(i,ivar)
           end do
        end do
     end do
     ! End loop over cells
  end do
  ! End loop over grids

111 format('   Entering add_rt_sources for level ',I2)

END SUBROUTINE add_rt_sources

!************************************************************************
SUBROUTINE add_UV_background(ilevel)

! Inject radiation from RT source regions (from the RT namelist). Since
! the light sources are continuously emitting radiation, this is called
! continuously during code execution, rather than just during
! initialization.
!
! ilevel => amr level at which to inject the radiation
!------------------------------------------------------------------------
  use UV_module, ONLY: UV_Nphot_cgs, nUVgroups, iUVgroups
  use amr_commons
  use rt_parameters
  use hydro_commons
  use rt_hydro_commons
  implicit none
  integer::ilevel,i,igrid,ncache,iskip,ngrid,j,ind,ic,ig
  integer ,dimension(1:nvector),save::ind_grid
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_np  &
            ,scale_fp,efactor,nH
!------------------------------------------------------------------------
  if(numbtot(1,ilevel)==0)return     ! no grids at this level
  if(.not. rt_isDiffuseUVsrc) return ! no propagated UV background
  if(verbose)write(*,111)ilevel

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  call rt_units(scale_np, scale_fp)

  ncache=active(ilevel)%ngrid
  ! Loop over grids by vector sweeps
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     ! Loop over cells
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ic=iskip+ind_grid(i) ! cell index
           ! Read the local gas density and inject the UV radiation
           nH=unew(ic,1)*scale_nH
           efactor = exp(-nH/rt_UVsrc_nHmax)
           do j=1,nUVgroups
              ig = iGroups(iUVgroups(j))
              rtunew(ic,ig) = max(rtunew(ic,ig)                     &
                                 ,UV_Nphot_cgs(j)/scale_np * efactor)
           end do
        end do
     end do
     ! End loop over cells
  end do
  ! End loop over grids

111 format('   Entering add_UV_background for level ',I2)

END SUBROUTINE add_UV_background

!************************************************************************
SUBROUTINE rt_sources_vsweep(x,uu,dx,dt,nn)

! Do a vector sweep, injecting RT source regions into cells, that is if
! they are in any of these regions.
!
! x      =>  ncells*ndim: positions of grid cells
! uu    <=  ncells*nrtvars: injected rt variables in each cell
! dx     =>  real cell width in code units
! dt     =>  real timestep length in code units
! nn     =>  int number of cells
!------------------------------------------------------------------------
  use amr_commons
  use rt_parameters
  implicit none
  integer ::nn
  real(dp)::dx,dt,dx_cgs,dt_cgs
  real(dp),dimension(1:nvector,1:nrtvar)::uu
  real(dp),dimension(1:nvector,1:ndim)::x
  integer::i,k,group_ind
  real(dp)::vol,r,xn,yn,zn,en
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::scale_np,scale_fp
!------------------------------------------------------------------------
  ! Initialize everything to zero
  !  uu=0.0d0
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  call rt_units(scale_np, scale_fp)
  dx_cgs=dx*scale_l
  dt_cgs=dt*scale_t
  ! Loop over RT regions
  do k=1,rt_nsource

     ! Find which photon group we should be contributing to
     if(rt_src_group(k) .le. 0 .or. rt_src_group(k) .gt. nGroups) cycle
     group_ind = iGroups(rt_src_group(k))
     ! For "square" regions only:
     if(rt_source_type(k) .eq. 'square')then
       ! Exponent of choosen norm
        en=rt_exp_source(k)
        do i=1,nn
           ! Compute position in normalized coordinates
           xn=0.0d0; yn=0.0d0; zn=0.0d0
           xn=2.0d0*abs(x(i,1)-rt_src_x_center(k))/rt_src_length_x(k)
#if NDIM>1
           yn=2.0d0*abs(x(i,2)-rt_src_y_center(k))/rt_src_length_y(k)
#endif
#if NDIM>2
           zn=2.0d0*abs(x(i,3)-rt_src_z_center(k))/rt_src_length_z(k)
#endif
           ! Compute cell "radius" relative to region center
           if(rt_exp_source(k)<10)then
              r=(xn**en+yn**en+zn**en)**(1.0/en)
           else
              r=max(xn,yn,zn)
           end if
           ! If cell lies within region, inject value
           if(r<1.0)then
              uu(i,group_ind) = rt_n_source(k)/rt_c_cgs/scale_Np
              ! The input flux is the fraction Fp/(c*Np) (Max 1 magnitude)
              uu(i,group_ind+1) =                                       &
                        rt_u_source(k) * rt_c * rt_n_source(k) / scale_Np
#if NDIM>1
              uu(i,group_ind+2) =                                       &
                        rt_v_source(k) * rt_c * rt_n_source(k) / scale_Np
#endif
#if NDIM>2
              uu(i,group_ind+3) =                                       &
                        rt_w_source(k) * rt_c * rt_n_source(k) / scale_Np
#endif
           end if
        end do
     end if

     ! For "point" regions only:
     if(rt_source_type(k) .eq. 'point')then
        ! Volume elements
        vol=dx_cgs**ndim
        ! Compute CIC weights relative to region center
        do i=1,nn
           xn=1.0; yn=1.0; zn=1.0
           xn=max(1.0-abs(x(i,1)-rt_src_x_center(k))/dx, 0.0_dp)
#if NDIM>1
           yn=max(1.0-abs(x(i,2)-rt_src_y_center(k))/dx, 0.0_dp)
#endif
#if NDIM>2
           zn=max(1.0-abs(x(i,3)-rt_src_z_center(k))/dx, 0.0_dp)
#endif
           r=xn*yn*zn
           if(r .gt. 0.) then
              ! If cell lies within CIC cloud, inject value.
              ! Photon input is in # per sec...need to convert to uu
              uu(i,group_ind)=uu(i,group_ind)                            &
                            + rt_n_source(k) / scale_Np * r / vol * dt_cgs
              uu(i,group_ind+1)=uu(i,group_ind+1) + rt_u_source(k) *rt_c &
                            * rt_n_source(k) / scale_Np * r / vol * dt_cgs
#if NDIM>1
              uu(i,group_ind+2)=uu(i,group_ind+2) + rt_v_source(k) *rt_c &
                            * rt_n_source(k) / scale_Np * r / vol * dt_cgs
#endif
#if NDIM>2
              uu(i,group_ind+3)=uu(i,group_ind+3) + rt_w_source(k) *rt_c &
                            * rt_n_source(k) / scale_Np * r / vol * dt_cgs
#endif
           endif
        end do
     end if

     ! For shell regions only:
     if(rt_source_type(k) .eq. 'shell')then
        ! An emitting spherical shell with center coordinates given,
        ! along with inner and outer radius (rt_src_length_x,z).
        ! Compute CIC weights relative to region center
        do i=1,nn
           xn=0.0; yn=0.0; zn=0.0
           xn=max(abs(x(i,1)-rt_src_x_center(k)), 0.0_dp)
#if NDIM>1
           yn=max(abs(x(i,2)-rt_src_y_center(k)), 0.0_dp)
#endif
#if NDIM>2
           zn=max(abs(x(i,3)-rt_src_z_center(k)), 0.0_dp)
#endif
           r=sqrt(xn**2+yn**2+zn**2)
           if(r .gt. rt_src_length_x(k) .and. &
                r .lt. rt_src_length_y(k)) then
              ! If cell lies within CIC cloud, inject value
              ! photon input is in # per sec...need to convert to uu
              uu(i,group_ind)=rt_n_source(k) / scale_np
           endif
        end do
     end if
  end do

  return
END SUBROUTINE rt_sources_vsweep
