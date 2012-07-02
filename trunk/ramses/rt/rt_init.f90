!*************************************************************************
SUBROUTINE rt_init

!  Initialize everything for radiative transfer
!-------------------------------------------------------------------------
  use amr_commons
  use hydro_commons
  use rt_hydro_commons
  use rt_flux_module
  use rt_cooling_module, only: update_UVrates
  use SED_module
  implicit none
  integer:: i, ilevel, ivar, nvar_count
!-------------------------------------------------------------------------
  if(verbose)write(*,*)'Entering init_rt'
  
  call update_rt_c
  call update_UVrates(aexp)
  if(star)call update_SED_Pacprops

  ! Count the number of variables and check if ok:
  nvar_count=ndim+2                        ! Density, momenta and pressure
  if(metal) nvar_count=nvar_count+1
  if(delayed_cooling) nvar_count=nvar_count+1
  iIons=nvar_count+1                       !        Starting inxex of xion
  nvar_count=nvar_count+nIons              !                   Xion states
  
  !if(myid==1) print*,'indexes=',nvar_count,nvar,iIons,i_BCT,i_Lya,nrtvar
  if(nvar_count .gt. nvar) then 
     if(myid==1) then 
        write(*,*) 'rt_init(): Something wrong with NVAR.'
        write(*,*) 'Should have NVAR=2+ndim+dcool+metal+'
        write(*,*) 'nIons+1*BCT+2+Lya. STOPPING!'
        write(*,*) 'Have nvar=',nvar
        write(*,*) 'Should have nvar>=',nvar_count
        write(*,*) 'STOPPING!'
     endif
     call clean_stop
  endif

  ! Update hydro variable to the initial ionized species
  var_region(1:rt_nregion,iIons-ndim-2)=rt_xion_region(1:rt_nregion)

  do i=1,nPacs  ! Starting indices in uold and unew of each photon package
     iPac(i)=1+(ndim+1)*(i-1)
     if(nrestart.eq.0) then 
        rtuold(:,iPac(i))=smallNp
     endif
  end do
  if(trim(rt_flux_scheme).eq.'hll') rt_use_hll=.true.
  if(rt_use_hll) call read_hll_eigenvalues

  if(rt_star .or. sedprops_update .ge. 0) then
     call init_SED_table    ! init stellar energy distribution properties
     if(sedprops_update .ge. 0) &
          call update_SED_Pacprops ! set RT props from st. distributions 
  endif

  tot_cool_loopcnt=0 ; max_cool_loopcnt=0 ; n_cool_cells=0
  loopCodes=0
  tot_nPhot=0.d0 ;  step_nPhot=0.d0; step_nStar=0.d0

END SUBROUTINE rt_init

!************************************************************************
SUBROUTINE update_rt_c

! Update the speed of light for radiative transfer, in code units.
! This cannot be just a constant, since scale_v changes with time in 
! cosmological simulations.
!------------------------------------------------------------------------
  use rt_parameters
  use amr_commons
  implicit none
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
!------------------------------------------------------------------------
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  rt_c=rt_c_cgs/scale_v
  rt_c2=rt_c**2
END SUBROUTINE update_rt_c

!************************************************************************
SUBROUTINE adaptive_rt_c_update(ilevel, dt)

! Set the lightspeed such that RT can be done at ilevel in time dt in 
! a single step.
!------------------------------------------------------------------------
  use amr_parameters
  use rt_parameters
  use SED_module
  implicit none
  integer:: ilevel, nx_loc
  real(dp):: dt, scale, dx
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
!------------------------------------------------------------------------
  ! Mesh spacing at ilevel
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5D0**ilevel*scale

  ! new lightspeed
  rt_c = dx/3.d0/dt * rt_courant_factor 
  rt_c2 = rt_c**2

  ! new ligtspeed in cgs
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  rt_c_cgs = rt_c*scale_v
  rt_c_fraction = rt_c_cgs/c_cgs

  call updateRTPac_CoolConstants           ! These change as a consequence

END SUBROUTINE adaptive_rt_c_update


!************************************************************************
SUBROUTINE read_rt_params(nml_ok)

! Read rt_params namelist
!------------------------------------------------------------------------
  use amr_commons
  use rt_parameters
  use rt_cooling_module
  use UV_module
  use SED_module
  implicit none
  logical::nml_ok
!------------------------------------------------------------------------
  namelist/rt_params/rt, rt_pp, rt_coupling_method, rt_star, rt_esc_frac&
       & ,rt_flux_scheme, rt_smooth, rt_is_outflow_bound                &
       & ,rt_max_subcycles, rt_courant_factor, rt_c_fraction, rt_otsa   &
       & ,sedprops_update, hll_evals_file, sed_dir, uv_file             &
       & ,rt_UVsrc_nHmax, nUVpacs, rt_freeflow                          &
       & ,upload_equilibrium_x, X, Y, rt_is_init_xion, rt_UV_nhSS       &
       & ,rt_err_grad_n, rt_floor_n, rt_err_grad_xHII, rt_floor_xHII    &
       & ,rt_err_grad_xHI, rt_floor_xHI, rt_refine_aexp                 &
       ! RT regions (for initialization)
       & ,rt_nregion, rt_region_type                                    &
       & ,rt_reg_x_center, rt_reg_y_center, rt_reg_z_center             &
       & ,rt_reg_length_x, rt_reg_length_y, rt_reg_length_z             &
       & ,rt_exp_region                                                 &
       & ,rt_n_region, rt_u_region, rt_v_region, rt_w_region            &
       & ,rt_xion_region                                                &
       ! RT source regions (for every timestep)
       & ,rt_nsource, rt_source_type                                    &
       & ,rt_src_x_center, rt_src_y_center, rt_src_z_center             &
       & ,rt_src_length_x, rt_src_length_y, rt_src_length_z             &
       & ,rt_exp_source, rt_src_pac                                     &
       & ,rt_n_source, rt_u_source, rt_v_source, rt_w_source            &
       ! RT boundary (for boundary conditions)
       & ,rt_n_bound,rt_u_bound,rt_v_bound,rt_w_bound
  ! Read namelist file
  rewind(1)
  read(1,NML=rt_params,END=101)
101 continue                                   ! No harm if no rt namelist

  rt_coupling_full=.false.  !
  rt_coupling_pp=.false.    !       These are all false if not rt or rt-pp
  rt_coupling_adc=.false.   !
  rt_coupling_dtrt=.false.  !

  if(nPacs.le.0) rt=.false. ! No sense in doing rt if there are no photons
  if(rt_pp) rt=.true.
  if(.not. rt .and. .not. rt_star) then
     rt_pp=.false.
     sedprops_update=-1
  endif
  !if(rt  .or. rt_star .and. .not. rt_pp) then
     select case(rt_coupling_method)
     case(:1)
        rt_coupling_full=.true.
     case(2)
        rt_coupling_pp=.true.
        if(myid==1 .and. (rt_star .or. rt)) write(*,*) &
             'I suggest you use rt_max_subcycle with this RT method'
     case(3)
        rt_coupling_adc=.true.
     case(4:)
        rt_coupling_dtrt=.true.
     end select
  !endif

  if(rt_err_grad_n .gt. 0. .or. rt_err_grad_xHII .gt. 0.                 &
       .or. rt_err_grad_xHI .gt. 0.) rt_refine=.true.

  rt_c_cgs = c_cgs * rt_c_fraction
  !call update_rt_c
  if(haardt_madau) rt_UV_hom=.true.                     ! UV in every cell
  call read_rt_pacs(nml_ok)
END SUBROUTINE read_rt_params

!*************************************************************************
SUBROUTINE read_rt_pacs(nml_ok)

! Read rt_pacs namelist
!-------------------------------------------------------------------------
  use amr_commons
  use rt_parameters
  use SED_module
  implicit none
  logical::nml_ok
  integer::i
  namelist/rt_pacs/pac_csn, pac_egy, spec2pac, pacL0, pacL1
!------------------------------------------------------------------------
  if(myid==1) then
     write(*,'(" Working with ",I2," photon packages and  " &
          ,I2, " ion species")')nPacs,nIons
     write(*,*) ''
  endif
   
  if(nPacs .le. 0) then
     rt = .false.
     return
  endif
  ! Use the ionization energies for HI, HeI, HeII as default package intervals
  pacL0(1:min(nPacs,3))=ionEvs(1:min(nPacs,3))  !    Lower interval bounds
  pacL1(1:min(nPacs,2))=ionEvs(2:min(nPacs+1,3))!    Upper interval bounds
  pacL1(min(nPacs,3))=0.                        !     Upper bound=infinity

  ! Default pacs are all blackbodies at E5 Kelvin
  pac_csn(1,:)=(/3.007d-18, 0d0, 0d0/)     ! avg photoion. c-section (cm2)
#if NPACS>1
  if(nPacs .ge. 2) pac_csn(2,:)=(/5.687d-19, 4.478d-18, 0d0/)
#endif
#if NPACS>2
  if(nPacs .ge. 3) pac_csn(3,:)=(/7.889d-20, 1.197d-18, 1.055d-18/)
#endif
  pac_egy(1,:)=(/17.439, 0., 0./)           ! cs weighted photon Energy (eV)
#if NPACS>1
  if(nPacs .ge. 2) pac_egy(2,:)=(/31.100, 32.359, 0./)
#endif
#if NPACS>2
  if(nPacs .ge. 3) pac_egy(3,:)=(/62.056, 62.696, 62.307/)
#endif
  do i=1,min(nIons,nPacs)
     spec2pac(i)=i                  ! Species contributions to packages
  end do
  
  ! Read namelist file
  rewind(1)
  read(1,NML=rt_pacs,END=101)
101 continue              ! no harm if no rt namelist

  call updateRTPac_CoolConstants
  call write_PacProps(.false.,6)
END SUBROUTINE read_rt_pacs

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
  integer::i,k,pac_ind
  real(dp)::vol,r,xn,yn,zn,en
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v, & 
            scale_np, scale_pf
!------------------------------------------------------------------------
  ! Initialize everything to zero
  !  uu=0.0d0
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  call rt_units(scale_np, scale_pf)
  dx_cgs=dx*scale_l
  dt_cgs=dt*scale_t
  ! Loop over RT regions
  do k=1,rt_nsource

     ! Find which photon package we should be contributing to
     if(rt_src_pac(k) .le. 0 .or. rt_src_pac(k) .gt. nPacs) cycle
     pac_ind = ipac(rt_src_pac(k))
     !pac_ind = 1+(ndim+1)*(rt_src_pac(k)-1)
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
              uu(i,pac_ind) = uu(i,pac_ind)+rt_n_source(k)/scale_Np
              ! The input flux is the fraction Fp/(c*Np) (Max 1 magnitude)
              uu(i,pac_ind+1) = uu(i,pac_ind+1)                           &
                        + rt_u_source(k) * rt_c * rt_n_source(k) / scale_Np
#if NDIM>1 
              uu(i,pac_ind+2) = uu(i,pac_ind+2)                           &
                        + rt_v_source(k) * rt_c * rt_n_source(k) / scale_Np
#endif
#if NDIM>2
              uu(i,pac_ind+3) = uu(i,pac_ind+3)                           &
                        + rt_w_source(k) * rt_c * rt_n_source(k) / scale_Np
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
              uu(i,pac_ind)=uu(i,pac_ind)                                &
                            + rt_n_source(k) / scale_Np * r / vol * dt_cgs
              ! The input flux is the fraction Fp/(c*Np) (Max 1 magnitude)
              uu(i,pac_ind+1)=uu(i,pac_ind+1) + rt_u_source(k) * rt_c    &
                            * rt_n_source(k) / scale_Np * r / vol * dt_cgs
#if NDIM>1
              uu(i,pac_ind+2)=uu(i,pac_ind+2) + rt_v_source(k) * rt_c    &
                            * rt_n_source(k) / scale_Np * r / vol * dt_cgs
#endif
#if NDIM>2
              uu(i,pac_ind+3)=uu(i,pac_ind+3) + rt_w_source(k) * rt_c    &
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
              uu(i,pac_ind)=rt_n_source(k) / scale_np
           endif
        end do
     end if
  end do

  return
END SUBROUTINE rt_sources_vsweep



