!*************************************************************************
SUBROUTINE update_sink_RT_feedback(ilevel)
! Turn on RT advection if needed.
! Update photon group properties from stellar populations.
!-------------------------------------------------------------------------
  use amr_parameters
  use amr_commons
  use rt_parameters
  use pm_commons
  implicit none
  integer::ilevel
  logical,save::groupProps_init=.false.  
!-------------------------------------------------------------------------
  if (nsink > 0) then
     if (.not. rt_advect)print*,'turning RT advection on'
     rt_advect=.true.
  else
     rt_advect=.false.
  end if
END SUBROUTINE update_sink_RT_feedback

!*************************************************************************
SUBROUTINE sink_RT_feedback(ilevel, dt)

! This routine adds photons from radiating sinks to appropriate cells in 
! the  hydro grid.
! ilegel =>  grid level in which to perform the feedback
! ti     =>  initial time for the timestep (code units)
! dt     =>  real timestep length in code units
!-------------------------------------------------------------------------
  use pm_commons
  use amr_commons
  use rt_parameters
  implicit none
  integer:: ilevel
  real(dp):: dt
  integer:: igrid, jgrid, ipart, jpart, next_part
  integer:: i, ig, ip, npart1, npart2, icpu
  integer,dimension(1:nvector),save:: ind_grid, ind_part, ind_grid_part
!-------------------------------------------------------------------------
  if(.not.rt_advect)RETURN
  if(nsink .le. 0 ) return
  if(numbtot(1,ilevel)==0)return ! number of grids in the level
  if(verbose)write(*,111)ilevel
  ! Gather star particles only.
  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel) ! grid index
     ig=0                     ! index of grid with stars (in ind_grid)
     ip=0                     ! index of star particle   (in ind_part)
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)   ! Number of particles in the grid
        npart2=0              ! number of selected (i.e. sink and cloud) particles
        
        ! Count sink and cloud particles in the grid
        if(npart1 > 0)then
          ipart = headp(igrid)        ! particle index
           ! Loop over particles
           do jpart = 1, npart1
              ! Save next particle       <--- Very important !!!
              next_part = nextp(ipart)
              !if(idp(ipart) .gt. 0 .and. tp(ipart) .ne. 0.d0) then 
              ! changed .ne. to .lt. in order to get the cloud particles, as done in accrete_sink in sink_particle.f90
              if(idp(ipart) .lt. 0) then 
                 npart2 = npart2+1     ! only sink and cloud particles
              endif
              ipart = next_part        ! Go to next particle
           end do
        endif

        ! Gather sink and cloud particles within the grid
        if(npart2 > 0)then        
           ig = ig+1
           ind_grid(ig) = igrid
           ipart = headp(igrid)
           ! Loop over particles
           do jpart = 1, npart1
              ! Save next particle      <--- Very important !!!
              next_part = nextp(ipart)
              ! Select only cloud particles
              ! changed .ne. to .lt. in order to get the cloud particles, as done in accrete_sink in sink_particle.f90
              if(idp(ipart) .lt. 0) then
                 if(ig==0)then
                    ig=1      
                    ind_grid(ig)=igrid
                 end if
                 ip = ip+1
                 ind_part(ip) = ipart
                 ind_grid_part(ip) = ig ! points to grid a cloud particle is in  
              endif
              if(ip == nvector)then
                 call sink_RT_vsweep( &
                      ind_grid,ind_part,ind_grid_part,ig,ip,dt,ilevel)
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
        call sink_RT_vsweep( &
             ind_grid,ind_part,ind_grid_part,ig,ip,dt,ilevel)
     endif
  end do 
  ! End loop over cpus

111 format('   Entering star_rt_feedback for level ',I2)

END SUBROUTINE sink_RT_feedback

!*************************************************************************
SUBROUTINE sink_RT_vsweep(ind_grid,ind_part,ind_grid_part,ng,np,dt,ilevel)

! This routine is called by subroutine sink_rt_feedback.
! Each sink and cloud  particle dumps a number of photons into the nearest grid cell
! using array rtunew.
! Radiation is injected into cells at level ilevel, but it is important
! to know that ilevel-1 cells may also get some radiation. This is due
! to sink and cloud particles that have just crossed to a coarser level. 
!
! ind_grid      =>  grid indexes in amr_commons (1 to ng)
! ind_part      =>  sink indexes in pm_commons(1 to np)
! ind_grid_part =>  points from star to grid (ind_grid) it resides in
! ng            =>  number of grids
! np            =>  number of sink and cloud particles
! dt            =>  timestep length in code units
! ilevel        =>  amr level at which we're adding radiation
!-------------------------------------------------------------------------
  use amr_commons
  use pm_commons
  use rt_hydro_commons
  use rt_parameters
  use rt_cooling_module, only:iIR
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  real(dp)::dt
  !-----------------------------------------------------------------------
  integer::i,j,idim,nx_loc,ip,isink
  real(dp)::dx,dx_loc,scale,vol_loc
  logical::error
  ! Grid based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle based arrays
  integer,dimension(1:nvector),save::igrid_son,ind_son
  logical,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector,ngroups),save::part_NpInp
  real(dp),dimension(1:nvector,1:ndim),save::x
  integer ,dimension(1:nvector,3),save::id=0,igd=0,icd=0
  integer ,dimension(1:nvector),save::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc
  real(dp)::Ep2Np
  ! units and temporary quantities
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v, scale_Np  & 
            , scale_Fp, age, z, scale_Nphot                   &
            , dt_loc_Gyr, scale_msun
  real(dp),parameter::vol_factor=2**ndim   ! Vol factor for ilevel-1 cells
  ! changed:
  real(dp),dimension(1:nvector),save::dn
  real(dp)::hnu,scale_energy,ener
  real(dp),dimension(1:5)::frac_ener
  integer::g
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
  scale = boxlen/dble(nx_loc) ! usually scale == 1
  dx_loc = dx*scale
  vol_loc = dx_loc**ndim
  scale_nPhot = vol_loc * scale_np * scale_l**ndim / 1.d50
  scale_msun = scale_d * scale_l**ndim / m_sun    

  frac_ener = (/0.3,0.25,0.079,0.067,0.085/)


  ! Conversion factor from Energy density in code units to photon density in code units
  Ep2Np=(scale_d * scale_v**2)/( scale_Np * group_egy(1) * ev_to_erg)


  ! Lower left corners of 3x3x3 grid-cubes (with given grid in center)
  do idim = 1, ndim
     do i = 1, ng
        x0(i,idim) = xg(ind_grid(i),idim) - 3.0D0*dx
     end do
  end do

  ! Gather 27 neighboring father cells (should be present anytime !)
  do i=1,ng
     ind_cell(i) = father(ind_grid(i))
  end do
  call get3cubefather(&
          ind_cell, nbors_father_cells, nbors_father_grids, ng, ilevel)
  ! now nbors_father cells are a cube of 27 cells with ind_cell in the 
  ! middle and nbors_father_grids are the 8 grids that contain these 27
  ! cells (though only one of those is fully included in the cube)


  ! Rescale position of stars to positions within 3x3x3 cell supercube
  do idim = 1, ndim
     do j = 1, np
        x(j,idim) = xp(ind_part(j),idim)/scale + skip_loc(idim)
        x(j,idim) = x(j,idim) - x0(ind_grid_part(j),idim)
        x(j,idim) = x(j,idim)/dx 
        ! so 0<x<2 is bottom cell, ...., 4<x<6 is top cell 
     end do
  end do

   ! NGP at level ilevel
  do idim=1,ndim
     do j=1,np
        id(j,idim) = x(j,idim) ! So id=0-5 is the cell (in the 
     end do                    ! 3x3x3 supercube) containing the star
  end do

  ! Compute parent grids
  do idim = 1, ndim
     do j = 1, np
        igd(j,idim) = id(j,idim)/2 ! must be 0, 1 or 2
     end do
  end do
  do j = 1, np
     kg(j) = 1 + igd(j,1) + 3*igd(j,2) + 9*igd(j,3) ! 1 to 27
  end do
  do j = 1, np
     igrid(j) = son(nbors_father_cells(ind_grid_part(j),kg(j))) 
     ! grid (not cell) containing the cloud particle 
  end do

  ! Check if particles are entirely in level ilevel.
  ! This should always be the case in post-processing
  ok(1:np) = .true.
  do j = 1, np
     ok(j) = ok(j) .and. igrid(j) > 0
  end do ! if ok(j) is true then particle j's cell contains a grid. 
  ! Otherwise it is a leaf, and we need to fill it with radiation.

  ! Compute parent cell position within it's grid
  do idim = 1, ndim
     do j = 1, np
        if( ok(j) ) then
           icd(j,idim) = id(j,idim) - 2*igd(j,idim) ! 0 or 1
        end if
     end do
  end do
  do j = 1, np
     if( ok(j) ) then
        icell(j) = 1 + icd(j,1) + 2*icd(j,2) + 4*icd(j,3) ! 1 to 8
     end if
  end do

  ! Compute parent cell adress (and particle radiation contribution)
  do j = 1, np

     if( ok(j) )then
        indp(j) = ncoarse + (icell(j)-1)*ngridmax + igrid(j)
     else
        indp(j) = nbors_father_cells(ind_grid_part(j),kg(j))
     end if
  end do


  do g=1,nGroups
     iGroups(g)=1+(ndim+1)*(g-1)

     ! Increase photon density in cell due to accretion luminosity
     do j=1,np
        ! it's ok to respect only highest ilevel particles since sinks are refined to the max level
        if( ok(j) ) then 
          isink=-idp(ind_part(j))
          ! energy increase due to accretion luminosity 
          !ener=acc_lum(isink)*dt*frac_ener(g)/(group_egy(g)*1.6022e-12)
          ener=acc_lum(isink)*frac_ener(g)/(group_egy(g)*1.6022e-12)
          !increase in energy density per cloud particle (code units)
          dn(j)=ener/dble(ncloud_sink)/vol_loc

          !increase in photon number density (code units)
          dn(j)=dn(j)*Ep2Np

          ! depositing the photons onto the grid
          rtunew(indp(j),iGroups(g))=rtunew(indp(j),iGroups(g))+dn(j)
        endif
     end do
  end do
    
END SUBROUTINE sink_RT_vsweep
