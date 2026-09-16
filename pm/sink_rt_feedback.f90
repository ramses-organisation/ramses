#ifdef RT
#if NDIM==3
!*************************************************************************
!*************************************************************************
!*************************************************************************
!*************************************************************************
SUBROUTINE update_sink_RT_feedback

! Turn on RT advection if needed.
! Update photon group properties from stellar populations.
!-------------------------------------------------------------------------
  use rt_parameters
  use sink_feedback_parameters
  implicit none

  if(nstellar>0)then
     rt_advect=.true.
  endif

END SUBROUTINE update_sink_RT_feedback
!*************************************************************************
!*************************************************************************
!*************************************************************************
!*************************************************************************
SUBROUTINE sink_RT_feedback(ilevel, dt)

! This routine adds photons from radiating sinks to appropriate cells in
! the  hydro grid. Emission is determined from massive stellar particles
! attached to the sinks.
! ilevel =>  grid level in which to perform the feedback
! dt     =>  real timestep length in code units
!-------------------------------------------------------------------------
  use pm_commons
  use amr_commons
  use rt_parameters
  use sink_feedback_parameters

  implicit none
  integer:: ilevel
  real(dp):: dt
  integer:: igrid, jgrid, ipart, jpart, next_part
  integer:: ig, ip, npart1, npart2, icpu
  integer,dimension(1:nvector),save:: ind_grid, ind_part, ind_grid_part
  !this array gathers the ionising flux by looping over stellar object
  !note that this array is local and therefore is not declared in pm_common
  real(dp),dimension(1:nsink,1:ngroups):: sink_ioni_flux
!-------------------------------------------------------------------------
  if(.not.rt_advect)RETURN
  if(nsink .le. 0 ) return
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Start by looping over the stellar objects and gather their fluxes
  call gather_ioni_flux(dt,sink_ioni_flux)

  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)
        npart2=0
        if(npart1 > 0)then
          ipart = headp(igrid)
           ! Loop over particles
           do jpart = 1, npart1
              next_part = nextp(ipart)
              ! only sink cloud particles
              if(idp(ipart) .lt. 0) then
                 npart2 = npart2+1
              endif
              ipart = next_part
           end do
        endif

        ! Gather sink and cloud particles within the grid
        if(npart2 > 0)then
           ig = ig+1
           ind_grid(ig) = igrid
           ipart = headp(igrid)
           ! Loop over particles
           do jpart = 1, npart1
              next_part = nextp(ipart)
              if(idp(ipart) .lt. 0) then
                 if(ig==0)then
                    ig=1
                    ind_grid(ig)=igrid
                 end if
                 ip = ip+1
                 ind_part(ip) = ipart
                 ind_grid_part(ip) = ig
              endif
              if(ip == nvector)then
                 call sink_RT_vsweep_stellar( &
                              ind_grid,ind_part,ind_grid_part,ig,ip,dt,ilevel,sink_ioni_flux)
                 ip = 0
                 ig = 0
              end if
              ipart = next_part  ! Go to next particle
           end do
           ! End loop over particles
        end if
        igrid = next(igrid)   ! Go to next grid
     end do
     ! End loop over grids
     if(ip > 0) then
         call sink_RT_vsweep_stellar( &
                     ind_grid,ind_part,ind_grid_part,ig,ip,dt,ilevel,sink_ioni_flux)
     endif
  end do
  ! End loop over cpus

111 format('   Entering sink_rt_feedback for level ',I2)

END SUBROUTINE sink_RT_feedback
!*************************************************************************
!*************************************************************************
!*************************************************************************
!*************************************************************************
SUBROUTINE gather_ioni_flux(dt,sink_ioni_flux)
! This routine is called by sink_RT_feedback is stellar objects are used
! It gathers the ionising flux on each sinks which is used to perform ionising radiation feedback

! sink_ioni_flux =>  the ionising flux of each sink
  use amr_commons
  use pm_commons
  use rt_parameters
  use sink_feedback_parameters

  implicit none

  real(dp),intent(in)::dt
  real(dp),dimension(1:nsink,1:ngroups),intent(out):: sink_ioni_flux !this arrays gathers the ionising flux by looping over stellar object
  integer:: istellar,isink,ig
  real(dp)::M_stellar,Flux_stellar
  real(dp),dimension(1:ngroups)::nphotons

  sink_ioni_flux = 0d0

  do istellar=1,nstellar
    ! find corresponding sink
     isink = 1
     do while ((isink.le.nsink) .and. (id_stellar(istellar) .ne. idsink(isink)))
       isink = isink + 1
     end do
     if (isink.gt.nsink) then
       write(*,*)"BUG: COULD NOT FIND SINK"
       call clean_stop
     endif
     M_stellar = mstellar(istellar)
     ! Reset the photon counter
     nphotons = 0d0
     Flux_stellar = 0
     ! Use fit to Vacca+ 1996
     !check whether the object is emitting
     if (t - tstellar(istellar) < hii_t) then
        !remember vaccafits is in code units because the corresponding parameters have been normalised in read_stellar_params (stf_K and stf_m0)
        call vaccafit(M_stellar,Flux_stellar)
        nphotons(feedback_photon_group) = Flux_stellar
     endif

     do ig=1,ngroups
        ! Remove negative photon counts
        nphotons(ig) = max(nphotons(ig),0d0)
        sink_ioni_flux(isink,ig) = sink_ioni_flux(isink,ig) + nphotons(ig)
     enddo

  enddo

END SUBROUTINE gather_ioni_flux
!*************************************************************************
!*************************************************************************
!*************************************************************************
!*************************************************************************
SUBROUTINE sink_RT_vsweep_stellar(ind_grid,ind_part,ind_grid_part,ng,np,dt,ilevel,sink_ioni_flux)
! This routine is called by subroutine sink_rt_feedback.
!
! Each sink and cloud particle dumps its share of the sink ionising flux into
! array rtunew, spread over the cells of its CIC stencil at level ilevel.
! Only the CIC volumes that are leaf cells at ilevel are deposited into
! (rt_upload_fine takes care of non-leaf cells).
!
! KNOWN ISSUE: A cloud particle is only ever visited at the level of
! the grid it is attached to, because sink_RT_feedback runs from rt_step, i.e.
! after kill_tree_fine and before merge_tree_fine. So any part of its CIC
! stencil reaching into a coarser neighbour would be dropped, with no other
! level's call to pick it up. This causes the total radiation to be lower
! than it should be. 
! This behaviour is different from the cloud mass in cic_amr, where rho_fine
! runs before kill_tree_fine and therefore deposits every particle at every
! level, so the coarse level holds the full mass.
! SOLUTION: Setting sink_refine=.true. ensures all CIC cells of every cloud
! particle sit at nlevelmax_sink. So nothing is ever in coarser grids.
!
! The ionising flux of each sink must be provided in sink_ioni_flux
!
! ind_grid       =>  grid indexes in amr_commons (1 to ng)
! ind_part       =>  sink indexes in pm_commons(1 to np)
! ind_grid_part  =>  points from star to grid (ind_grid) it resides in
! ng             =>  number of grids
! np             =>  number of sink and cloud particles
! dt             =>  timestep length in code units
! ilevel         =>  amr level at which we're adding radiation
! sink_ioni_flux =>  the ionising flux of each sink
!-------------------------------------------------------------------------
  use amr_commons
  use pm_commons
  use rt_hydro_commons
  use rt_parameters
  use sink_feedback_parameters
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  real(dp)::dt
  !-----------------------------------------------------------------------
  integer::j,ind,idim,nx_loc,isink,ig
  real(dp)::dx,dx_loc,scale,vol_cgs
  ! Particle based arrays
  real(dp),dimension(1:nvector,1:ndim),save::xpart
  real(dp),dimension(1:nvector,1:ndim,1:twotondim),save::xx
  real(dp),dimension(1:nvector,1:twotondim),save::vol
  integer ,dimension(1:nvector,1:twotondim),save::indp
  logical ,dimension(1:nvector,1:twotondim),save::ok
  real(dp),dimension(1:3)::skip_loc
  ! units and temporary quantities
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_Np,scale_Fp
  !this arrays gather the ionising flux by looping over stellar object
  real(dp),dimension(1:nsink,1:ngroups):: sink_ioni_flux

!-------------------------------------------------------------------------
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  call rt_units(scale_Np, scale_Fp)

  ! Mesh spacing in ilevel
  dx = 0.5D0**ilevel
  nx_loc = (icoarse_max - icoarse_min + 1)
  skip_loc = (/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1) = dble(icoarse_min)
  if(ndim>1)skip_loc(2) = dble(jcoarse_min)
  if(ndim>2)skip_loc(3) = dble(kcoarse_min)
  scale = boxlen/dble(nx_loc)
  dx_loc = dx*scale

  vol_cgs = (dx_loc*scale_l)**ndim

  ! Get cloud particle CIC weights
  do idim = 1, ndim
     do j = 1, np
        xpart(j,idim) = xp(ind_part(j),idim)
     end do
  end do
  call cic_get_cells(indp,xx,vol,ok,ind_grid,xpart,ind_grid_part,ng,np,ilevel)

  ! Increase photon density in cell due to accretion luminosity.
  ! Only the CIC volumes that are leaf cells at ilevel are used, as in cic_amr.
  ! Loop over eight CIC volumes
  do ind=1,twotondim
     do j=1,np
        if( ok(j,ind) ) then   ! ilevel cell
           ! Get sink index
           isink=-idp(ind_part(j))
           ! deposit the photons onto the grid
           do ig=1,ngroups
              rtunew(indp(j,ind),iGroups(ig))=rtunew(indp(j,ind),iGroups(ig)) + vol(j,ind) * &
                   sink_ioni_flux(isink,ig) * dt / dble(ncloud_sink) / vol_cgs / scale_Np
           end do
        endif
     end do
  end do

END SUBROUTINE sink_RT_vsweep_stellar
!*************************************************************************
!*************************************************************************
!*************************************************************************
!*************************************************************************
SUBROUTINE vaccafit(M,S)
  ! perform a fit of the Vacca et al. 96 ionising flux
  ! M - stellar mass / solar masses
  ! S - photon emission rate in / s

  use amr_parameters,only:dp
  use sink_feedback_parameters
  implicit none

  real(dp),intent(in)::M
  real(dp),intent(out)::S

  S = stf_K * (M / stf_m0)**stf_a / (1. + (M / stf_m0)**stf_b)**stf_c

END SUBROUTINE
!################################################################
!################################################################
!################################################################
!################################################################
#endif
#endif