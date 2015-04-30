!################################################################
!################################################################
!################################################################
!################################################################
subroutine create_sink
  use amr_commons
  use pm_commons
  use hydro_commons
  use clfind_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  !----------------------------------------------------------------------------
  ! sink creation routine
  ! -remove all cloud particles, keep only global sink arrays
  ! -call clumpfinder for potential relevant peaks
  ! -flag peaks which are eligible for sink formation (flag 2)
  ! -One true RAMSES particle is created 
  ! -Sink cloud particles are created
  ! -Cloud particles are scattered to grid
  ! -Accretion routine is called (with on_creation=.true.)
  !----------------------------------------------------------------------------

  integer::ilevel,ivar
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v

  if(verbose)write(*,*)' Entering create_sink'

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Merge all particles to level 1
  do ilevel=levelmin-1,1,-1
     call merge_tree_fine(ilevel)
  end do
  
  ! Remove all particle clouds around old sinks (including the central one)
  call kill_entire_cloud(1) 
  
  ! DO NOT MODIFY FLAG2 BETWEEN CLUMP_FINDER AND MAKE_SINK_FROM_CLUMP     
  if (create_sinks)then
     
     ! Run the clump finder,(produce no output, keep clump arrays allocated)
     call clump_finder(.false.,.true.)
     
     ! Trim clumps down to R_accretion ball around peaks 
     if(clump_core)call trim_clumps
     
     ! Compute simple additive quantities and means (1st moments)
     call compute_clump_properties(uold(1,1))
     
     ! Compute quantities relative to mean (2nd moments)
     call compute_clump_properties_round2(uold(1,1))
  
     ! Apply all checks and flag cells for sink formation
     call flag_formation_sites

     ! Create new sink particles if relevant
     do ilevel=levelmin,nlevelmax
        call make_sink_from_clump(ilevel)
     end do

     ! Deallocate clump finder arrays
     deallocate(npeaks_per_cpu)
     deallocate(clump_mass4)!!
     deallocate(ipeak_start)
     if (ntest>0)then
        deallocate(icellp)
        deallocate(levp)
        deallocate(testp_sort)
        deallocate(imaxp)
     endif
     call deallocate_all

  end if

  !merge sinks - for star formation runs
  if (merging_scheme == 'timescale')call merge_star_sink
  
  ! Merge sink for smbh runs 
  if (smbh)call merge_smbh_sink
  
  ! Create new cloud particles
  call create_cloud_from_sink

  ! Scatter particle to the grid
  ! Compute Bondi parameters
  do ilevel=1,nlevelmax
     call make_tree_fine(ilevel)
     call kill_tree_fine(ilevel)
     call virtual_tree_fine(ilevel)
     call collect_acczone_avg(ilevel)
  end do

  ! Perform first accretion with seed mass
  ! Gather particles to levelmin
  do ilevel=nlevelmax,levelmin,-1
     call grow_sink(ilevel,.true.)
     call merge_tree_fine(ilevel)
  end do

  ! Update hydro quantities for split cells
  if(hydro)then
     do ilevel=nlevelmax,levelmin,-1
        call upload_fine(ilevel)
        do ivar=1,nvar
           call make_virtual_fine_dp(uold(1,ivar),ilevel)
        end do
        ! Update boundaries 
        if(simple_boundary)call make_boundary_hydro(ilevel)
     end do
  end if

  ! Update the cloud particle properties at levelmin
  ! So they match their parent sink properties
  call update_cloud(levelmin)

  ! Effective accretion rate during last coarse step
  acc_rate(1:nsink)=acc_rate(1:nsink)/dtnew(levelmin)
  ! ir_eff and 5 are ratio of infalling energy which is radiated and protostellar radius
  if(ir_feedback)acc_lum(1:nsink)=ir_eff*acc_rate(1:nsink)*msink(1:nsink)/(5*6.955d10/scale_l)

  ! Compute and print accretion rates
  call compute_accretion_rate(.true.)
  
end subroutine create_sink
!################################################################
!################################################################
!################################################################
!################################################################
subroutine create_cloud_from_sink
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  !----------------------------------------------------------------------
  ! This routine creates the whole cloud of particles for one sink in one go, 
  ! including the central one. Particles are produced in the right MPI domain
  ! at level 1. Not very efficient...
  ! replaces: 
  !   - create_part_from_sink
  !   - create_cloud
  !   - mk_cloud
  !----------------------------------------------------------------------

  real(dp)::scale,dx_min,rr,rmax,rmass
  integer ::i,icpu,isink,indp,ii,jj,kk,nx_loc,idim
  integer ::ntot,ntot_all,info
  logical ::ok_free
  real(dp),dimension(1:ndim)::xrel
  real(dp),dimension(1:nvector,1:ndim)::xtest
  integer ,dimension(1:nvector)::ind_grid,ind_part,cc,ind_cloud
  logical ,dimension(1:nvector)::ok_true
  logical,dimension(1:ndim)::period

  ok_true=.true.

  if(numbtot(1,1)==0) return
  if(verbose)write(*,*)' Entering create_cloud_from_sink'

#if NDIM==3

  ! Level 1 linked list
  do icpu=1,ncpu
     if(numbl(icpu,1)>0)then
        ind_grid(1)=headl(icpu,1)
     endif
  end do

  period(1)=(nx==1)
  if(ndim>1)period(2)=(ny==1)
  if(ndim>2)period(3)=(nz==1)
  
  ! Mesh spacing in that level
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax/aexp

  rmax=dble(ir_cloud)*dx_min
  rmass=dble(ir_cloud_massive)*dx_min
  
  do kk=-2*ir_cloud,2*ir_cloud
     xrel(3)=dble(kk)*dx_min/2.0
     do jj=-2*ir_cloud,2*ir_cloud
        xrel(2)=dble(jj)*dx_min/2.0
        do ii=-2*ir_cloud,2*ir_cloud
           xrel(1)=dble(ii)*dx_min/2.0
           rr=sqrt(sum(xrel**2))
           if(rr<=rmax)then
              do isink=1,nsink
                 xtest(1,1:3)=xsink(isink,1:3)+xrel(1:3)
                 do idim=1,ndim
                    if (period(idim) .and. xtest(1,idim)>boxlen)xtest(1,idim)=xtest(1,idim)-boxlen
                    if (period(idim) .and. xtest(1,idim)<0.)xtest(1,idim)=xtest(1,idim)+boxlen
                 end do
                 call cmp_cpumap(xtest,cc,1)
                 if(cc(1).eq.myid)then                    
                    call remove_free(ind_cloud,1)
                    call add_list(ind_cloud,ind_grid,ok_true,1)
                    indp=ind_cloud(1)
                    idp(indp)=-isink
                    levelp(indp)=levelmin
                    if (rr<=rmass .and. msink(isink)<msink_direct)then
                       mp(indp)=msink(isink)/dble(ncloud_sink_massive)
                    else
                       mp(indp)=0.
                    end if
                    xp(indp,1:3)=xtest(1,1:3)
                    vp(indp,1:3)=vsink(isink,1:3)
                    tp(indp)=tsink(isink)     ! Birth epoch
                 end if
              end do
           end if
        end do
     end do
  end do
  
  do isink=1,nsink
     direct_force_sink(isink)=(msink(isink) .ge. msink_direct)
  end do

#endif
end subroutine create_cloud_from_sink
!################################################################
!################################################################
!################################################################
!################################################################
subroutine kill_entire_cloud(ilevel)
  use pm_commons
  use amr_commons
  implicit none
  integer::ilevel
  !------------------------------------------------------------------------
  ! This routine removes cloud particles (including the central one).
  !------------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,next_part
  integer::ig,ip,npart1,npart2,icpu
  integer,dimension(1:nvector)::ind_grid,ind_part,ind_grid_part
  logical,dimension(1:nvector)::ok=.true.

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel
  ! Gather sink and cloud particles.
  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0        
        ! Count sink and cloud particles
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if(idp(ipart).lt.0)then
                 npart2=npart2+1
              endif
              ipart=next_part  ! Go to next particle
           end do
        endif        
        ! Gather sink and cloud particles
        if(npart2>0)then        
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              ! Select only sink particles
              if(idp(ipart).lt.0)then
                 if(ig==0)then
                    ig=1
                    ind_grid(ig)=igrid
                 end if
                 ip=ip+1
                 ind_part(ip)=ipart
                 ind_grid_part(ip)=ig   
              endif
              if(ip==nvector)then
                 call remove_list(ind_part,ind_grid_part,ok,ip)
                 call add_free_cond(ind_part,ok,ip)
                 ip=0
                 ig=0
              end if
              ipart=next_part  ! Go to next particle
           end do
           ! End loop over particles
        end if
        
        igrid=next(igrid)   ! Go to next grid
     end do
     
     ! End loop over grids
     if(ip>0)then
        call remove_list(ind_part,ind_grid_part,ok,ip)
        call add_free_cond(ind_part,ok,ip)
     end if
  end do
111 format('   Entering kill_cloud for level ',I2)
end subroutine kill_entire_cloud
!################################################################
!################################################################
!################################################################
!################################################################
subroutine collect_acczone_avg(ilevel)
  use pm_commons
  use amr_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel

  !------------------------------------------------------------------------
  ! This routine is used to collect all relevant information to compute the 
  ! the accretion rates. The information is collected level-by-level when
  ! going down in the call tree (leafs at the bottom), while accretion is 
  ! performed on the way up.
  ! - first, the volume of each particle is computed. The volume is reduced if
  ! sinks are overlapping.
  ! - then a loop over all particles of ilevel (vectorized) is used to compute 
  ! the volume-weighted quantities.
  ! No gaussian accretion kernel is used anymore...
  !------------------------------------------------------------------------

  integer::igrid,jgrid,ipart,jpart,next_part,info,ind
  integer::ig,ip,npart1,npart2,icpu,nx_loc,isink
  integer,dimension(1:nvector)::ind_grid,ind_part,ind_grid_part
  real(dp)::dx_loc,dx_min,scale,factG
  character(LEN=15)::action

  if(ilevel<levelmin)return
  if(verbose)write(*,111)ilevel
  action='count'
  call count_clouds(ilevel,action)

  call make_virtual_reverse_dp(rho(1),ilevel)
  call make_virtual_fine_dp(rho(1),ilevel)

  action='weight'
  call count_clouds(ilevel,action)
  
  level_sink_new(1:nsinkmax,ilevel)=.false.
  
  ! Gravitational constant
  factG=1d0
  if(cosmo)factG=3d0/8d0/3.1415926*omega_m*aexp

  ! Mesh spacing in that level
  dx_loc=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax/aexp

  ! compute (volume weighted) averages over accretion zone
  wden=0d0; wvol=0d0; weth=0d0; wmom=0d0; wdiv=0d0
  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0
        
        ! Count sink and cloud particles
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if(idp(ipart).lt.0)then
                 npart2=npart2+1
              endif
              ipart=next_part  ! Go to next particle
           end do
        endif
        
        ! Gather sink and cloud particles
        if(npart2>0)then
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              ! Select only sink particles
              if(idp(ipart).lt.0)then
                 if(ig==0)then
                    ig=1
                    ind_grid(ig)=igrid
                 end if
                 ip=ip+1
                 ind_part(ip)=ipart
                 ind_grid_part(ip)=ig   
              endif
              if(ip==nvector)then
                 call collect_acczone_avg_np(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
                 ip=0
                 ig=0
              end if
              ipart=next_part  ! Go to next particle
           end do
           ! End loop over particles
        end if
        igrid=next(igrid)   ! Go to next grid
     end do

     ! End loop over grids
     if(ip>0)then
        call collect_acczone_avg_np(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
     end if
  end do
  ! End loop over cpus

  if(nsink>0)then
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(wden,wden_new,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(wvol,wvol_new,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(weth,weth_new,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(wmom,wmom_new,nsinkmax*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(wdiv,wdiv_new,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(level_sink_new,level_sink,nsinkmax*(nlevelmax-levelmin+1),MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,info)
#else
     wden_new=wden
     wvol_new=wvol
     weth_new=weth
     wmom_new=wmom
     wdiv_new=wdiv
     level_sink=level_sink_new
#endif
  endif
  do isink=1,nsink
     weighted_density(isink,ilevel)=wden_new(isink)
     weighted_volume(isink,ilevel)=wvol_new(isink)
     weighted_momentum(isink,ilevel,1:ndim)=wmom_new(isink,1:ndim)
     weighted_ethermal(isink,ilevel)=weth_new(isink)
     weighted_divergence(isink,ilevel)=wdiv_new(isink)
  end do

111 format('   Entering collect_acczone_avg for level ',I2)

end subroutine collect_acczone_avg
!################################################################
!################################################################
!################################################################
!################################################################
subroutine collect_acczone_avg_np(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid,ind_part,ind_grid_part

  !-----------------------------------------------------------------------
  ! inner loop of collect_acczone_avg
  ! weighted gas quantities and divergence for each particle are computed
  ! -> replaces bondi_average and divergence_sink
  ! no CIC averaging over quantities anymore as average over whole sink 
  ! accretion zone is computed
  !-----------------------------------------------------------------------

  integer::j,irad,nx_loc,isink,divdim,idim,ind
  real(dp)::d,u,v=0d0,w=0d0,e
  real(dp)::scale,weight,dx_min,one_over_dx_min
#ifdef SOLVERmhd
  real(dp)::bx1,bx2,by1,by2,bz1,bz2
#endif
  real(dp),dimension(1:nvector),save::egas,divpart
  real(dp),dimension(1:nvector,1:ndim),save::xpart
  real(dp) ,dimension(1:nvector,1:nvar),save::fluid_var_left,fluid_var_right,fluid_var
  integer ,dimension(1:nvector),save::cind,cind_right,cind_left
  
  do idim=1,ndim
     do j=1,np
        xpart(j,idim)=xp(ind_part(j),idim)
     end do
  end do

  call cic_get_vals(fluid_var,ind_grid,xpart,ind_grid_part,ng,np,ilevel,.true.)

  ! compute thermal energy
  do j=1,np
     d=fluid_var(j,1)
     if(d>0.)then
        u=fluid_var(j,2)/d
        v=fluid_var(j,3)/d
        w=fluid_var(j,4)/d
        e=fluid_var(j,5)
#ifdef SOLVERmhd
        bx1=fluid_var(j,6)
        by1=fluid_var(j,7)
        bz1=fluid_var(j,8)
        bx2=fluid_var(j,nvar+1)
        by2=fluid_var(j,nvar+2)
        bz2=fluid_var(j,nvar+3)
        e=e-0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)
#endif
        e=e-0.5*d*(u*u+v*v+w*w)
#if NENER>0
        do irad=1,nener
           e=e-fluid_var(j,5+irad)
        end do
#endif
        egas(j)=e
     else
        egas(j)=0.
     end if
  end do
    
  ! Mesh spacing in that level
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax/aexp
  one_over_dx_min=1./dx_min
  
  divpart(1:np)=0._dp

  do divdim=1,3

     !use min cell spacing to obtain right position, get right cell index
     do j=1,np
        xpart(j,divdim)=xpart(j,divdim)+0.5*dx_min
     end do
     call cic_get_vals(fluid_var_right,ind_grid,xpart,ind_grid_part,ng,np,ilevel,.false.)

     !use min cell spacing to obtain left position, get left cell index
     do j=1,np
        xpart(j,divdim)=xpart(j,divdim)-dx_min
     end do
     call cic_get_vals(fluid_var_left,ind_grid,xpart,ind_grid_part,ng,np,ilevel,.false.)

     !back to original position
     do j=1,np
        xpart(j,divdim)=xpart(j,divdim)+0.5*dx_min
     end do

     !compute divergence of (rho*v - rho*vsink) in one go
     do j=1,np
        isink=-idp(ind_part(j))
        divpart(j)=divpart(j)+(fluid_var_right(j,divdim+1)-fluid_var_right(j,1)*vsink(isink,divdim)-&
             fluid_var_left(j,divdim+1)+fluid_var_left(j,1)*vsink(isink,divdim))*one_over_dx_min
     end do
  end do

  do j=1,np
     isink=-idp(ind_part(j))
     weight=0._dp
     do ind=1,twotondim
        weight=weight+weightp(ind_part(j),ind)
     end do
     wvol(isink)=wvol(isink)+weight
     wden(isink)=wden(isink)+weight*fluid_var(j,1)
     wmom(isink,1:3)=wmom(isink,1:3)+weight*fluid_var(j,2:4)
     weth(isink)=weth(isink)+weight*egas(j)
     wdiv(isink)=wdiv(isink)+weight*divpart(j)
     if (weight>0.)level_sink_new(isink,ilevel)=.true.
  end do

end subroutine collect_acczone_avg_np
!################################################################
!################################################################
!################################################################
!################################################################
subroutine grow_sink(ilevel,on_creation)
  use pm_commons
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  logical::on_creation
  !------------------------------------------------------------------------
  ! This routine performs accretion onto the sink. It vectorizes the loop
  ! over all sink cloud particles and calls accrete_sink as soon as nvector 
  ! particles are collected
  ! -> replaces grow_bondi and grow_jeans
  !------------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,next_part,info
  integer::ig,ip,npart1,npart2,icpu,isink,lev,nx_loc
  integer,dimension(1:nvector)::ind_grid,ind_part,ind_grid_part
  real(dp)::scale,dx_min,vol_min
  real(dp),dimension(1:ndim)::old_loc,old_vel

  if(accretion_scheme=='none')return
  if(verbose)write(*,111)ilevel

  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_min=(0.5D0**nlevelmax)*scale
  vol_min=dx_min**ndim

#if NDIM==3

  ! Compute sink accretion rates
  if ((.not. on_creation) .and. (.not. threshold_accretion))call compute_accretion_rate(.false.)
  !if(agn) call average_AGN

  ! Reset new sink variables
  msink_new=0d0; xsink_new=0.d0; vsink_new=0d0; delta_mass_new=0d0; lsink_new=0d0

  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0
        ! Count sink and cloud particles
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if(idp(ipart).lt.0)then
                 npart2=npart2+1
              endif
              ipart=next_part  ! Go to next particle
           end do
        endif
        ! Gather sink and cloud particles
        if(npart2>0)then        
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              ! Select only sink particles
              if(idp(ipart).lt.0)then
                 if(ig==0)then
                    ig=1
                    ind_grid(ig)=igrid
                 end if
                 ip=ip+1
                 ind_part(ip)=ipart
                 ind_grid_part(ip)=ig   
              endif
              if(ip==nvector)then
                 call accrete_sink(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel,on_creation)
                 ip=0
                 ig=0
              end if
              ipart=next_part  ! Go to next particle
           end do
           ! End loop over particles
        end if
        igrid=next(igrid)   ! Go to next grid
     end do
     ! End loop over grids
     if(ip>0)call accrete_sink(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel,on_creation)
  end do
  ! End loop over cpus
  if(nsink>0)then
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(msink_new,msink_all,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(delta_mass_new,delta_mass_all,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(xsink_new,xsink_all,nsinkmax*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(vsink_new,vsink_all,nsinkmax*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(lsink_new,lsink_all,nsinkmax*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#else
     msink_all=msink_new
     delta_mass_all=delta_mass_new
     xsink_all=xsink_new
     vsink_all=vsink_new
     lsink_all=lsink_new
#endif
  endif
  do isink=1,nsink
     ! Reset jump in sink coordinates
     do lev=levelmin,nlevelmax
        sink_jump(isink,1:ndim,lev)=sink_jump(isink,1:ndim,lev)-xsink(isink,1:ndim)
     end do
     
     !save old velocity and location to compute deltas
     old_loc(1:ndim)=xsink(isink,1:ndim)
     old_vel(1:ndim)=vsink(isink,1:ndim)
     
     ! Change to conservative quantities
     xsink(isink,1:ndim)=xsink(isink,1:ndim)*msink(isink)
     vsink(isink,1:ndim)=vsink(isink,1:ndim)*msink(isink)
     

     ! Accrete to sink variables
     msink(isink)=msink(isink)+msink_all(isink)
     xsink(isink,1:ndim)=xsink(isink,1:ndim)+xsink_all(isink,1:ndim)
     vsink(isink,1:ndim)=vsink(isink,1:ndim)+vsink_all(isink,1:ndim)
     !compute lsink with reference point of old xsink
     lsink(isink,1:3)=lsink(isink,1:3)+lsink_all(isink,1:3)
     
     ! Change back
     xsink(isink,1:ndim)=xsink(isink,1:ndim)/msink(isink)
     vsink(isink,1:ndim)=vsink(isink,1:ndim)/msink(isink)

     !correct for new center of mass location/velocity
     lsink(isink,1:3)=lsink(isink,1:3)-msink(isink)*cross((xsink(isink,1:ndim)-old_loc(1:ndim)),(vsink(isink,1:ndim)-old_vel(1:ndim)))

     ! Store jump in sink coordinates
     do lev=levelmin,nlevelmax
        sink_jump(isink,1:ndim,lev)=sink_jump(isink,1:ndim,lev)+xsink(isink,1:ndim)
     end do

     ! Store accreted mass
     acc_rate(isink)=acc_rate(isink)+msink_all(isink)
     if (agn.and.ok_blast_agn(isink))then
       delta_mass(isink)=delta_mass_all(isink) ! only the most recent accretion
     else
       delta_mass(isink)=delta_mass(isink)+delta_mass_all(isink)
     end if

  end do
  
#endif
111 format('   Entering grow_sink for level ',I2)

end subroutine grow_sink
!################################################################
!################################################################
!################################################################
!################################################################
subroutine accrete_sink(ind_grid,ind_part,ind_grid_part,ng,np,ilevel,on_creation)
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  logical::on_creation

  !-----------------------------------------------------------------------
  ! This routine is called by subroutine grow_sink. It performs accretion
  ! for nvector particles. Routine is not very efficient. Optimize if taking too long...
  !-----------------------------------------------------------------------

  integer::i,j,nx_loc,isink,ivar,irad,idim,ind
  real(dp)::v2,d,e,d_floor,density,volume
#ifdef SOLVERmhd
  real(dp)::bx1,bx2,by1,by2,bz1,bz2
#endif
  real(dp),dimension(1:nvar)::z
  real(dp)::factG,scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::dx,dx_loc,dx_min,scale,vol_min,vol_loc,weight,acc_mass,temp,d_jeans
  ! Grid based arrays
  real(dp),dimension(1:nvector,1:ndim)::xpart
  real(dp),dimension(1:nvector,1:ndim,1:twotondim)::xx
  real(dp),dimension(1:nvector,1:twotondim)::vol
  integer ,dimension(1:nvector)::cell_lev
  ! Particle based arrays
  logical,dimension(1:nvector,1:twotondim)::ok
  integer ,dimension(1:nvector,1:twotondim)::indp
  real(dp),dimension(1:3)::vv

  real(dp),dimension(1:3)::r_rel,x_acc,p_acc,p_rel,p_rel_rad,p_rel_acc,p_rel_tan,delta_x,delta_p,drag
  real(dp)::r_abs,fbk_ener,T2_AGN
  logical,dimension(1:ndim)::period
  real(dp)::virt_acc_mass,delta_e_tot,Mred,Macc
  real(dp),dimension(1:nsinkmax)::delta_M

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  
  period(1)=(nx==1)
#if NDIM>1
  if(ndim>1)period(2)=(ny==1)
#endif
#if NDIM>2
  if(ndim>2)period(3)=(nz==1)
#endif
  
#if NDIM==3
  
  ! Gravitational constant
  factG=1d0
  if(cosmo)factG=3d0/8d0/3.1415926*omega_m*aexp
  
  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim
  dx_min=scale*0.5D0**nlevelmax/aexp
  vol_min=dx_min**ndim

  T2_AGN=0.15*1d12/scale_T2

  do idim=1,ndim
     do j=1,np
        xpart(j,idim)=xp(ind_part(j),idim)
     end do
  end do

  virt_acc_mass=0.d0; delta_M=0.d0

  call cic_get_cells(indp,xx,vol,ok,ind_grid,xpart,ind_grid_part,ng,np,ilevel)

  do ind=1,twotondim
     do j=1,np
        if(ok(j,ind))then
           
           ! Convert uold to primitive variables
           d=uold(indp(j,ind),1)
           vv(1)=uold(indp(j,ind),2)/d
           vv(2)=uold(indp(j,ind),3)/d
           vv(3)=uold(indp(j,ind),4)/d
           e=uold(indp(j,ind),5)

#ifdef SOLVERmhd
           bx1=uold(indp(j,ind),6)
           by1=uold(indp(j,ind),7)
           bz1=uold(indp(j,ind),8)
           bx2=uold(indp(j,ind),nvar+1)
           by2=uold(indp(j,ind),nvar+2)
           bz2=uold(indp(j,ind),nvar+3)
           e=e-0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)
#endif
           v2=(vv(1)**2+vv(2)**2+vv(3)**2)
           e=e-0.5d0*d*v2
#if NENER>0
           do irad=1,nener
              e=e-uold(indp(j,ind),5+irad)
           end do
#endif
           e=e/d
           do ivar=imetal,nvar
              z(ivar)=uold(indp(j,ind),ivar)/d
           end do

           ! Get sink index
           isink=-idp(ind_part(j))        
 
           ! Drag force based on virtual accretion
           if (sink_drag)then
              Macc=max((dMBHoverdt(isink)-dMsink_overdt(isink))*dtnew(ilevel),0.d0)
              Mred=msink(isink)*(rho_gas(isink)*volume_gas(isink))/(msink(isink)+(rho_gas(isink)*volume_gas(isink)))
              delta_M(isink)=Mred*Macc/(Mred+Macc)
           end if
           delta_M(isink)=min(delta_M(isink),msink(isink)*2*dx_loc/(sum((xsink(isink,1:ndim)-xx(j,1:ndim,ind))**2)**0.5))
           
           ! Compute sink average density
           weight=weightp(ind_part(j),ind)
           density=rho_gas(isink)
           volume=volume_gas(isink)
           if (volume==0. .or. density==0.)print*,'something might be going wrong here...',weight,volume,density,ilevel

           if (on_creation)then
              if (new_born(isink))then
                 ! on sink creation, new sinks
                 acc_mass=mseed(isink)*weight/volume*d/density
              else
                 ! on sink creation, preexisting sinks
                 acc_mass=0.         
              end if
           else
              if (flux_accretion .or. bondi_accretion)then              
                 acc_mass=dMsink_overdt(isink)*dtnew(ilevel)*weight/volume*d/density
                 virt_acc_mass=delta_M(isink)*weight/volume*d/density
                 fbk_ener=delta_mass(isink)*T2_AGN*weight/volume*d/density
              end if

              if (threshold_accretion)then
                 ! User defined density threshold
                 d_floor=d_sink           
                 ! Jeans length related density threshold  
                 if(d_sink<0.0)then
                    temp=max(e*(gamma-1.0),smallc**2)
                    d_jeans=temp*3.1415926/(4.0*dx_loc)**2/factG
                    d_floor=d_jeans
                 endif
                 acc_mass=c_acc*weight*(d-d_floor)
              end if

              ! No neg accretion
              acc_mass=max(acc_mass,0.0_dp)               

           end if

           ! reference frame relative to the sink
           r_rel(1:3)=xx(j,1:3,ind)-xsink(isink,1:3) 
           do idim=1,ndim
              if (period(idim) .and. r_rel(idim)>boxlen*0.5)xx(j,idim,ind)=xx(j,idim,ind)-boxlen
              if (period(idim) .and. r_rel(idim)<boxlen*(-0.5))xx(j,idim,ind)=xx(j,idim,ind)+boxlen
           end do
           r_rel(1:3)=xx(j,1:3,ind)-xsink(isink,1:3) 
           r_abs=sum(r_rel(1:3)**2)**0.5      

           ! momentum in relative motion
           p_rel=d*vol_loc*(vv(1:3)-vsink(isink,1:3))
           p_rel_rad=sum(r_rel(1:3)*p_rel(1:3))*r_rel(1:3)/(r_abs**2+tiny(0.d0))
           p_rel_tan=p_rel-p_rel_rad

           if(nol_accretion)then
              p_rel_acc=p_rel_rad*acc_mass/(d*vol_loc)
           else
              p_rel_acc=p_rel*acc_mass/(d*vol_loc)
           end if

           ! total accreted/deccreted momentum
           drag=(virt_acc_mass*vsink(isink,1:ndim)*density/d-virt_acc_mass*vv(1:ndim)) 
           delta_p=-drag
           p_acc(1:3)=p_rel_acc+vsink(isink,1:3)*acc_mass
          
           ! total accreted/deccreted energy
           delta_e_tot=(virt_acc_mass*eps_sink(isink)*density/d-virt_acc_mass*e)+sum(drag(1:ndim)*vv(1:ndim))

           ! total accreted center of mass
           x_acc(1:3)=acc_mass*xx(j,1:3,ind)
           delta_x(1:3)=virt_acc_mass*(xx(j,1:3,ind)-xsink(isink,1:3))
           
           ! Add accreted properties to sink variables
           msink_new(isink)=msink_new(isink)+acc_mass
           delta_mass_new(isink)=delta_mass_new(isink)+acc_mass
           xsink_new(isink,1:3)=xsink_new(isink,1:3)+x_acc(1:3)+delta_x(1:3)
           vsink_new(isink,1:3)=vsink_new(isink,1:3)+p_acc(1:3)+delta_p(1:3)
           lsink_new(isink,1:3)=lsink_new(isink,1:3)+cross(r_rel(1:3),p_rel_acc(1:3))

           ! Check for neg density inside sink accretion radius
           if (8*acc_mass/vol_loc>d .and. (.not. on_creation))then 
              write(*,*)'====================================================='
              write(*,*)'DANGER of neg density :-( at location'
              write(*,*)xx(j,1:3,ind)
              write(*,*)'due to',isink
              write(*,*)acc_mass/vol_loc,d/d_sink,density/d_sink
              write(*,*)indp(j,ind),myid
              write(*,*)'nol_accretion: ',nol_accretion
              do i=1,nsink
                 write(*,*)i,' distance ',sum((xx(j,1:3,ind)-xsink(i,1:3))**2)**0.5/dx_min
              end do
              write(*,*)'try to decrease c_acc in SINK_PARAMS'
              write(*,*)'====================================================='
              !           call clean_stop
           end if

           if (on_creation)then
              ! new born sinks accrete from uold
              d=d-acc_mass/vol_loc
              !new gas velocity
              vv(1:3)=(d*vol_loc*vv(1:3)-p_acc(1:3))/(d*vol_loc-acc_mass)
              !convert back to conservative variables
              v2=(vv(1)**2+vv(2)**2+vv(3)**2)
              e=e*d
#ifdef SOLVERmhd
              e=e+0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)
#endif
              e=e+0.5d0*d*v2
#if NENER>0
              do irad=1,nener
                 e=e+uold(indp(j,ind),5+irad)
              end do
#endif
              uold(indp(j,ind),1)=d
              uold(indp(j,ind),2)=d*vv(1)
              uold(indp(j,ind),3)=d*vv(2)
              uold(indp(j,ind),4)=d*vv(3)
              uold(indp(j,ind),5)=e
              do ivar=imetal,nvar
                 uold(indp(j,ind),ivar)=d*z(ivar)
              end do
           else
              ! regular accretion from unew
              unew(indp(j,ind),1)=unew(indp(j,ind),1)-acc_mass/vol_loc-virt_acc_mass/vol_loc+virt_acc_mass/vol_loc*density/d
              unew(indp(j,ind),2:4)=unew(indp(j,ind),2:4)-vv(1:3)*acc_mass/vol_loc-delta_p(1:3)/vol_loc
              ! fix the energy equation to account for the work of the truly accreted material 
              unew(indp(j,ind),5)=unew(indp(j,ind),5)-uold(indp(j,ind),5)*acc_mass/(d*vol_loc)+delta_e_tot/vol_loc
              do ivar=imetal,nvar
                 unew(indp(j,ind),ivar)=unew(indp(j,ind),ivar)-uold(indp(j,ind),ivar)*acc_mass/(d*vol_loc)
              end do

              ! put the tangential momentum back into the gas
              if(nol_accretion)then
                 unew(indp(j,ind),2:4)=unew(indp(j,ind),2:4)+acc_mass/(d*vol_loc)*p_rel_tan(1:3)/vol_loc
                 unew(indp(j,ind),5)=unew(indp(j,ind),5)+acc_mass/(d*vol_loc)*sum(p_rel_tan(1:3)*uold(indp(j,ind),2:4)/d)
              end if
              
              ! AGN feedback
              if(agn)then
                if(ok_blast_agn(isink).and.delta_mass(isink)>0.0)then
                  unew(indp(j,ind),5)=unew(indp(j,ind),5)+fbk_ener/vol_loc
                end if
              end if
           end if

        endif
     end do
  end do
#endif
end subroutine accrete_sink
!################################################################
!################################################################
!################################################################
!################################################################
subroutine compute_accretion_rate(write_sinks)
  use pm_commons
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  logical::write_sinks

  !------------------------------------------------------------------------
  ! This routine computes the accretion rate onto the sink particles based
  ! on the information collected in collect accretion 
  ! It also creates output for the sink particle positions
  !------------------------------------------------------------------------

  integer::i,nx_loc,isink
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m
  real(dp)::factG,d_star,boost,fa_fact
  real(dp)::r2,v2,c2,density,volume,ethermal,dx_min,scale,mgas,rho_inf,divergence,v_bondi
  real(dp),dimension(1:3)::velocity
  real(dp),dimension(1:nsinkmax)::dMEDoverdt
  real(dp)::T2_gas,delta_mass_min
  real(dp)::T2_min,T2_AGN

  dt_acc=huge(0._dp)

  ! Gravitational constant
  factG=1d0
  if(cosmo)factG=3d0/8d0/3.1415926*omega_m*aexp

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_m=scale_d*scale_l**3d0
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax/aexp
  d_star=n_star/scale_nH
  T2_min=1d7/scale_T2
  T2_AGN=0.15*1d12/scale_T2

  ! Compute sink particle accretion rate by averaging contributions from all levels
  do isink=1,nsink
     density=0.d0; volume=0.d0; velocity=0.d0; ethermal=0d0
     divergence=0.d0
     do i=levelmin,nlevelmax
        density=density+weighted_density(isink,i)
        ethermal=ethermal+weighted_ethermal(isink,i)
        velocity(1:3)=velocity(1:3)+weighted_momentum(isink,i,1:3)
        volume=volume+weighted_volume(isink,i)
        divergence=divergence+weighted_divergence(isink,i)
     end do
     mgas=density
     density=density/(volume+tiny(0.0_dp))
     if (density>0)then
        velocity(1:3)=velocity(1:3)/(density*volume+tiny(0.0_dp))
        ethermal=ethermal/(density*volume+tiny(0.0_dp))
        c2=MAX((gamma-1.0)*ethermal,smallc**2)
        c2sink(isink)=c2
        v2=SUM((velocity(1:3)-vsink(isink,1:3))**2)
        v_bondi=sqrt(c2+v2)
        
        ! Compute Bondi-Hoyle accretion rate in code units
        if (star.and.alpha_acc_boost.lt.0.0)then
           boost=max((density/d_star)**2,1.0_dp)
        else
           boost=abs(alpha_acc_boost)
        end if
        
        v_bondi=v_bondi*boost**(-1./3.)
        
        ! Bondi radius
        r2=(factG*(msink(isink)+mgas)/v_bondi**2)**2
        
        ! extrapolate to rho_inf
        rho_inf=density/(bondi_alpha(ir_cloud*0.5*dx_min/r2**0.5))

        ! use other values for lambda depending on the EOS, this is for isothermal
        dMBHoverdt(isink)=4.*3.1415926*rho_inf*r2*v_bondi
        dMEDoverdt(isink)=4.*3.1415926*6.67d-8*msink(isink)*1.66d-24/(0.1*6.652d-25*3d10)*scale_t
        
        ! set sink accretion rate to Bondi value
        if(smbh.or.bondi_accretion)dMsink_overdt(isink)=dMBHoverdt(isink)

        ! limit accretion to Eddington rate
        if(eddington_limit)then
             dMsink_overdt(isink)=min(dMBHoverdt(isink),dMEDoverdt(isink))
        endif

        ! accretion rate is based on mass flux onto the sink, bondirate is used for subsonic accretion
        if (flux_accretion)then

           ! average divergence over all cloud particles and multiply by cloud volume
           dMsink_overdt(isink)=-1.*divergence

           ! correct for some small factor to keep density close to threshold
           fa_fact=(log10(density)-log10(d_sink))*0.1+1.
           dMsink_overdt(isink)=dMsink_overdt(isink)*fa_fact

           ! If accretion is subsonic (sonic radius smaller than accretion radius), use Bondi-rate instead.
           if ((0.5*msink(isink)/c2) < (ir_cloud*dx_min))then
              dMsink_overdt(isink)=dMBHoverdt(isink)
           end if
        end if

        ! store average quantities for accretion and diagnostics
        eps_sink(isink)=ethermal
        rho_gas(isink)=density
        volume_gas(isink)=volume
        vel_gas(isink,1:ndim)=velocity(1:ndim)

        ! make sure, accretion rate is positive
        dMsink_overdt(isink)=max(0.d0,dMsink_overdt(isink))

        ! compute maximum timestep allowed by sink
        if (dMsink_overdt(isink)>0.0)then

           ! make sure that sink doesnt accrete more than the gas mass within one timestep
           if (bondi_accretion .or. flux_accretion)dt_acc(isink)=c_acc*mgas/dMsink_overdt(isink)
           ! make sure that sink doesnt accrete more than its own mass within one timestep
           dt_acc(isink)=min(dt_acc(isink),(c_acc*msink(isink)/dMsink_overdt(isink)))
           
           ! make sure that agn feedback doesn't increase the thermal energy of the gas too much
           if(agn)then

              ! check whether we should have AGN feedback
              ok_blast_agn(1:nsink)=.false.
              T2_gas=ethermal
              delta_mass_min = mgas*(T2_min-T2_gas)/(T2_AGN-T2_min)
              if((T2_gas.ge.T2_min).or.(delta_mass(isink).ge.mgas*(T2_min-T2_gas)/(T2_AGN-T2_min)))then
                ok_blast_agn(isink)=.true.
              end if
              if(myid==1.and.ok_blast_agn(isink).and.delta_mass(isink).gt.0.)then
                write(*,'("***BLAST***",I4,1X,3(1PE12.5,1X))')isink &
                    & ,msink(isink)*scale_d*scale_l**3/2d33 &  
                    & ,delta_mass(isink)*scale_d*scale_l**3/2d33 &
                    & ,((delta_mass(isink)*T2_AGN+mgas*T2_gas) &
                    & /(delta_mass(isink)+mgas))*scale_T2
              endif

              ! make sure that sink doesnt blast more than the gas internal energy within one timestep
              if((T2_gas.ge.T2_min).or.(delta_mass(isink).ge.delta_mass_min))then
                 dt_acc(isink)=min(dt_acc(isink),(c_acc*mgas*T2_gas/(dMsink_overdt(isink)*T2_AGN)))
              end if
           end if
        end if
     end if
  end do 
  
  if (write_sinks)then 
     call print_sink_properties(dMEDoverdt,rho_inf,r2,v_bondi)
  end if

contains
  ! Routine to return alpha, defined as rho/rho_inf, for a critical
  ! Bondi accretion solution. The argument is x = r / r_Bondi.
  ! This is from Krumholz et al. (AJC)
  REAL(dp) function bondi_alpha(x)
    implicit none
    REAL(dp) x
    REAL(dp), PARAMETER :: XMIN=0.01, xMAX=2.0
    INTEGER, PARAMETER :: NTABLE=51
    REAL(dp) lambda_c, xtable, xtablep1, alpha_exp
    integer idx
    !     Table of alpha values. These correspond to x values that run from
    !     0.01 to 2.0 with uniform logarithmic spacing. The reason for
    !     this choice of range is that the asymptotic expressions are 
    !     accurate to better than 2% outside this range.
    REAL(dp), PARAMETER, DIMENSION(NTABLE) :: alphatable = (/ &
         820.254, 701.882, 600.752, 514.341, 440.497, 377.381, 323.427, &
         277.295, 237.845, 204.1, 175.23, 150.524, 129.377, 111.27, 95.7613, &
         82.4745, 71.0869, 61.3237, 52.9498, 45.7644, 39.5963, 34.2989, &
         29.7471, 25.8338, 22.4676, 19.5705, 17.0755, 14.9254, 13.0714, &
         11.4717, 10.0903, 8.89675, 7.86467, 6.97159, 6.19825, 5.52812, &
         4.94699, 4.44279, 4.00497, 3.6246, 3.29395, 3.00637, 2.75612, &
         2.53827, 2.34854, 2.18322, 2.03912, 1.91344, 1.80378, 1.70804, &
         1.62439 /)
    !     Define a constant that appears in these formulae
    lambda_c    = 0.25 * exp(1.5)
    !     Deal with the off-the-table cases
    if (x .le. XMIN) then
       bondi_alpha = lambda_c / sqrt(2. * x**3)
    else if (x .ge. XMAX) then
       bondi_alpha = exp(1./x)
    else
       !     We are on the table
       idx = floor ((NTABLE-1) * log(x/XMIN) / log(XMAX/XMIN))
       xtable = exp(log(XMIN) + idx*log(XMAX/XMIN)/(NTABLE-1))
       xtablep1 = exp(log(XMIN) + (idx+1)*log(XMAX/XMIN)/(NTABLE-1))
       alpha_exp = log(x/xtable) / log(xtablep1/xtable)
       !     Note the extra +1s below because of fortran 1 offset arrays
       bondi_alpha = alphatable(idx+1) * (alphatable(idx+2)/alphatable(idx+1))**alpha_exp
    end if
  end function bondi_alpha


end subroutine compute_accretion_rate
!################################################################
!################################################################
!################################################################
!################################################################
subroutine print_sink_properties(dMEDoverdt,rho_inf,r2,v_bondi)
  use pm_commons
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  real(dp),dimension(1:nsinkmax)::dMEDoverdt  
  integer::i,isink
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m
  real(dp)::l_abs
  real(dp)::v_bondi,r2,rho_inf

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_m=scale_d*scale_l**3d0
  
  if (smbh) then
     if(myid==1.and.nsink>0)then
        xmsink(1:nsink)=msink(1:nsink)
        call quick_sort_dp(xmsink(1),idsink_sort(1),nsink)
        write(*,*)'Number of sink = ',nsink
        write(*,'(" ============================================================================================")')
        write(*,'(" Id     Mass(Msol) Bondi(Msol/yr)   Edd(Msol/yr)              x              y              z")')
        write(*,'(" ============================================================================================")')
        do i=nsink,max(nsink-10,1),-1
           isink=idsink_sort(i)
           write(*,'(I3,10(1X,1PE14.7))')idsink(isink),msink(isink)*scale_m/2d33 &
                & ,dMBHoverdt(isink)*scale_m/scale_t/(2d33/(365.*24.*3600.)) &
                & ,dMEDoverdt(isink)*scale_m/scale_t/(2d33/(365.*24.*3600.)) &
                & ,xsink(isink,1:ndim),delta_mass(isink)*scale_m/2d33
        end do
        write(*,'(" ============================================================================================")')
        if(smbh_verbose)then
          write(*,'(" Id     rho(H/cc)  rho_inf(H/cc) Mgas(Msol) cs(km/s) vgas(km/s)  vsink(km/s) rBondi(pc)")')
          write(*,'(" vgas(km/s):  x   y   z     vsink(km/s):  x   y   z")')
          write(*,'(" ============================================================================================")')
          do i=nsink,max(nsink-10,1),-1
            isink=idsink_sort(i)
            write(*,'(I3,12(1X,1PE14.7))')idsink(isink),rho_gas(isink)*scale_nH,rho_inf*scale_nH &
                & ,rho_gas(isink)*volume_gas(isink)*scale_m/2d33,sqrt(c2sink(isink))*scale_v/1e5 &
                & ,sqrt(r2)*scale_l/3.086e18
            write(*,'(6(1X,1PE14.7))')vel_gas(isink,1:ndim)*scale_v/1e5,vsink(isink,1:ndim)*scale_v/1e5
          end do
          write(*,'(" ============================================================================================")')
        end if
     endif
  end if
  if (.not. smbh)then    
     if(myid==1.and.nsink>0.and. mod(nstep_coarse,ncontrol)==0)then
        xmsink(1:nsink)=msink(1:nsink)
        call quick_sort_dp(xmsink(1),idsink_sort(1),nsink)
        write(*,*)'Number of sink = ',nsink
        write(*,*)'Total mass in sink = ',sum(msink(1:nsink))*scale_m/1.9891d33
        write(*,*)'simulation time = ',t
        write(*,'(" ================================================================================================================================================ ")')
        write(*,'("   Id     M[Msol]          x             y             z            vx            vy            vz    acc_rate[Msol/y] acc_lum[Lsol]     age      ")')
        write(*,'(" ================================================================================================================================================ ")')
        do i=nsink,1,-1
           isink=idsink_sort(i)
           l_abs=(lsink(isink,1)**2+lsink(isink,2)**2+lsink(isink,3)**2)**0.5+1.d10*tiny(0.d0)
           write(*,'(I5,10(2X,E12.5))')&
                idsink(isink),msink(isink)*scale_m/1.9891d33, &
                xsink(isink,1:ndim),vsink(isink,1:ndim),&
                acc_rate(isink)*scale_m/1.9891d33/(scale_t)*365.*24.*3600.,acc_lum(isink)/scale_t**2*scale_l**3*scale_d*scale_l**2/scale_t/3.9d33,&
                (t-tsink(isink))*scale_t/(3600*24*365.25)
        end do
        write(*,'(" ================================================================================================================================================ ")')
     endif
  endif
end subroutine print_sink_properties
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine quenching(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
  integer::ilevel

  !------------------------------------------------------------------------
  ! This routine selects regions which are eligible for SMBH formation.
  ! It is based on a stellar density threshold and on a stellar velocity
  ! dispersion threshold.
  ! On exit, flag2 array is set to 0 for AGN sites and to 1 otherwise.
  !------------------------------------------------------------------------

  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::dx,dx_loc,scale,vol_loc
  real(dp)::str_d,tot_m,ave_u,ave_v,ave_w,sig_u,sig_v,sig_w
  integer::igrid,ipart,jpart,next_part,ind_cell,iskip,ind
  integer::i,npart1,npart2,nx_loc
  real(dp),dimension(1:3)::skip_loc

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim

#if NDIM==3
  ! Gather star particles only.

  ! Loop over grids
  do i=1,active(ilevel)%ngrid
     igrid=active(ilevel)%igrid(i)
     ! Number of particles in the grid
     npart1=numbp(igrid)
     npart2=0
     
     ! Reset velocity moments
     str_d=0.0
     tot_m=0.0
     ave_u=0.0
     ave_v=0.0
     ave_w=0.0
     sig_u=0.0
     sig_v=0.0
     sig_w=0.0
     
     ! Count star particles
     if(npart1>0)then
        ipart=headp(igrid)
        ! Loop over particles
        do jpart=1,npart1
           ! Save next particle   <--- Very important !!!
           next_part=nextp(ipart)
           if(idp(ipart).gt.0.and.tp(ipart).ne.0)then
              npart2=npart2+1
              tot_m=tot_m+mp(ipart)
              ave_u=ave_u+mp(ipart)*vp(ipart,1)
              ave_v=ave_v+mp(ipart)*vp(ipart,2)
              ave_w=ave_w+mp(ipart)*vp(ipart,3)
              sig_u=sig_u+mp(ipart)*vp(ipart,1)**2
              sig_v=sig_v+mp(ipart)*vp(ipart,2)**2
              sig_w=sig_w+mp(ipart)*vp(ipart,3)**2
           endif
           ipart=next_part  ! Go to next particle
        end do
     endif
     
     ! Normalize velocity moments
     if(npart2.gt.0)then
        ave_u=ave_u/tot_m
        ave_v=ave_v/tot_m
        ave_w=ave_w/tot_m
        sig_u=sqrt(sig_u/tot_m-ave_u**2)*scale_v/1d5
        sig_v=sqrt(sig_v/tot_m-ave_v**2)*scale_v/1d5
        sig_w=sqrt(sig_w/tot_m-ave_w**2)*scale_v/1d5
        str_d=tot_m/(2**ndim*vol_loc)*scale_nH
     endif
     
     ! Loop over cells
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        ind_cell=iskip+igrid
        ! AGN formation sites
        ! if n_star>0.1 H/cc and v_disp>100 km/s
        if(str_d>0.1.and.MAX(sig_u,sig_v,sig_w)>100.)then
           flag2(ind_cell)=0
        else
           flag2(ind_cell)=1
        end if
     end do
  end do
  ! End loop over grids

#endif

111 format('   Entering quenching for level ',I2)

end subroutine quenching
!################################################################
!################################################################
!################################################################
!################################################################
!################################################################
!################################################################
subroutine make_sink_from_clump(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  use clfind_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel

  !----------------------------------------------------------------------
  ! This routine uses creates a sink in every cell which was flagged (flag2)
  ! The global sink variables are updated
  ! The true RAMSES particle is NOT produced here...
  !----------------------------------------------------------------------

  integer ::ncache,nnew,ivar,irad,ngrid,icpu,index_sink,index_sink_tot
  integer ::igrid,ix,iy,iz,ind,i,iskip,isink,nx_loc
  integer ::ntot,ntot_all,info
  integer ,dimension(1:nvector)::ind_grid,ind_cell
  integer ,dimension(1:nvector)::ind_grid_new,ind_cell_new
  integer ,dimension(1:ncpu)::ntot_sink_cpu,ntot_sink_all
  integer::ipart,peak_nr,grid
  logical ::ok_free
  real(dp)::d,u,v,w,e,factG,delta_d,v2
  real(dp)::birth_epoch
  real(dp)::dx,dx_loc,scale,vol_loc
  real(dp)::fourpi,threepi2,tff,tsal
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp)::d_ball,vol,mclump
  real(dp),dimension(1:nlevelmax)::volume
  real(dp),dimension(1:nvar)::z
  real(dp),dimension(1:3)::skip_loc,x,xcell,xpeak
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp)::T2_min,T2_AGN
#ifdef SOLVERmhd
  real(dp)::bx1,bx2,by1,by2,bz1,bz2
#endif
  
  if(verbose)write(*,*)'entering make_sink_from_clump for level ',ilevel

  ! Gravitational constant
  factG=1d0
  if(cosmo)factG=3d0/8d0/3.1415926*omega_m*aexp

  ! Mesh spacing in that level
  dx=0.5D0**ilevel 
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  T2_min=1d7/scale_T2
  T2_AGN=0.15*1d12/scale_T2

  ! Birth epoch as proper time
  if(use_proper_time)then
     birth_epoch=texp
  else
     birth_epoch=t
  endif

  ! Cells center position relative to grid center position
  do ind=1,twotondim  
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     xc(ind,1)=(dble(ix)-0.5D0)*dx
     xc(ind,2)=(dble(iy)-0.5D0)*dx
     xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do
  
  ! Set new sink variables to zero
  msink_new=0d0; mseed_new=0d0; tsink_new=0d0; delta_mass_new=0d0; xsink_new=0d0; vsink_new=0d0
  oksink_new=0d0; idsink_new=0; new_born_new=.false.

#if NDIM==3

  !------------------------------------------------
  ! and count number of new sinks (flagged cells)
  !------------------------------------------------
  ntot=0
  ntot_sink_cpu=0
  if(numbtot(1,ilevel)>0)then
     ncache=active(ilevel)%ngrid
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do
        do ind=1,twotondim  
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do
           do i=1,ngrid
              if(flag2(ind_cell(i))>0)then
                 ntot=ntot+1
              end if
           end do
        end do
     end do
  end if

  !---------------------------------
  ! Compute global sink statistics
  !---------------------------------
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(ntot,ntot_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  ntot_all=ntot
#endif
#ifndef WITHOUTMPI
  ntot_sink_cpu=0; ntot_sink_all=0
  ntot_sink_cpu(myid)=ntot
  call MPI_ALLREDUCE(ntot_sink_cpu,ntot_sink_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  ntot_sink_cpu(1)=ntot_sink_all(1)
  do icpu=2,ncpu
     ntot_sink_cpu(icpu)=ntot_sink_cpu(icpu-1)+ntot_sink_all(icpu)
  end do
#endif
  nsink=nsink+ntot_all  
  nindsink=nindsink+ntot_all
  if(myid==1)then
     if(ntot_all.gt.0)then
        write(*,'(" Level = ",I6," New sinks produced from clumps= ",I6," Total sinks =",I8)')&
             & ilevel,ntot_all,nsink
     endif
  end if

  !-------------------------------------------
  ! Check wether max number of sink is reached
  !------------------------------------------
  ok_free=(nsink+ntot_all<=nsinkmax)
  if(.not. ok_free)then
     if(myid==1)write(*,*)'global list of sink particles is too long'
     if(myid==1)write(*,*)'New sink particles',ntot_all
     if(myid==1)write(*,*)'Increase nsinkmax'
     call clean_stop
  end if
   
  !------------------------------
  ! Create new sink particles
  !------------------------------
  ! Starting identity number
  if(myid==1)then
     index_sink=nsink-ntot_all
     index_sink_tot=nindsink-ntot_all
  else
     index_sink=nsink-ntot_all+ntot_sink_cpu(myid-1)
     index_sink_tot=nindsink-ntot_all+ntot_sink_cpu(myid-1)
  end if
  
  ! Loop over grids
  if(numbtot(1,ilevel)>0)then
     ncache=active(ilevel)%ngrid
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do

        ! Loop over cells
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do

           ! Gather cells with a new sink
           nnew=0
           do i=1,ngrid
              if (flag2(ind_cell(i))>0)then
                 nnew=nnew+1
                 ind_grid_new(nnew)=ind_grid(i)
                 ind_cell_new(nnew)=ind_cell(i)
              end if
           end do

           ! Create new sink particles
           do i=1,nnew
              index_sink=index_sink+1
              index_sink_tot=index_sink_tot+1

              ! Convert uold to primitive variables
              d=uold(ind_cell_new(i),1)
              u=uold(ind_cell_new(i),2)/d
              v=uold(ind_cell_new(i),3)/d
              w=uold(ind_cell_new(i),4)/d
              e=uold(ind_cell_new(i),5)
#ifdef SOLVERmhd
              bx1=uold(ind_cell_new(i),6)
              by1=uold(ind_cell_new(i),7)
              bz1=uold(ind_cell_new(i),8)
              bx2=uold(ind_cell_new(i),nvar+1)
              by2=uold(ind_cell_new(i),nvar+2)
              bz2=uold(ind_cell_new(i),nvar+3)
              e=e-0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)
#endif
              v2=(u**2+v**2+w**2)
              e=e-0.5d0*d*v2
#if NENER>0
              do irad=1,nener
                 e=e-uold(ind_cell_new(i),5+irad)
              end do
#endif             
              e=e/d
              do ivar=imetal,nvar
                 z(ivar)=uold(ind_cell_new(i),ivar)/d
              end do
              
              ! Get density maximum by quadratic expansion around cell center
              x(1)=(xg(ind_grid_new(i),1)+xc(ind,1)-skip_loc(1))*scale
              x(2)=(xg(ind_grid_new(i),2)+xc(ind,2)-skip_loc(2))*scale
              x(3)=(xg(ind_grid_new(i),3)+xc(ind,3)-skip_loc(3))*scale              
              call true_max(x(1),x(2),x(3),nlevelmax)

              ! Mass of the new sink
              if(sink_seedmass>=0.0)then
                 mseed_new(index_sink)=sink_seedmass
              else
                 if(smbh)then
                    ! The SMBH/sink mass is the mass that will heat the gas to 10**7 K after creation
                    fourpi=4.0d0*ACOS(-1.0d0)
                    threepi2=3.0d0*ACOS(-1.0d0)**2
                    if(cosmo)fourpi=1.5d0*omega_m*aexp
                    tff=sqrt(threepi2/8./fourpi/max(d,smallr))
                    tsal=0.1*6.652d-25*3d10/4./3.1415926/6.67d-8/1.66d-24/scale_t
                    mclump=clump_mass4(flag2(ind_cell_new(i)))
                    mseed_new(index_sink)=min(1.d-5/0.15*mclump*tsal/tff,mclump/2.0)
                 end if
                 
                 if(smbh.and.agn)then
                    mclump=clump_mass4(flag2(ind_cell_new(i)))
                    mseed_new(index_sink)=T2_min/T2_AGN*mclump
                 end if
              endif

              ! Give a tiny bit of mass to the sink...
              delta_d=d*1.d-10
              msink_new(index_sink)=delta_d*vol_loc                    
              delta_mass_new(index_sink)=msink_new(index_sink)

              ! Global index of the new sink
              oksink_new(index_sink)=1d0
              idsink_new(index_sink)=index_sink_tot

              ! Store properties of the new sink
              tsink_new(index_sink)=birth_epoch
              xsink_new(index_sink,1:3)=x(1:3)
              vsink_new(index_sink,1)=u
              vsink_new(index_sink,2)=v
              vsink_new(index_sink,3)=w
              new_born_new(index_sink)=.true.

              ! Convert back to conservative variable                                             
              d=d-delta_d
              e=e*d
#ifdef SOLVERmhd
              e=e+0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)
#endif
              e=e+0.5d0*d*(u**2+v**2+w**2)
#if NENER>0
              do irad=1,nener
                 e=e+uold(ind_cell_new(i),5+irad)
              end do
#endif              
              uold(ind_cell_new(i),1)=d
              uold(ind_cell_new(i),2)=d*u
              uold(ind_cell_new(i),3)=d*v
              uold(ind_cell_new(i),4)=d*w
              uold(ind_cell_new(i),5)=e
              do ivar=imetal,nvar
                 uold(ind_cell_new(i),ivar)=d*z(ivar)
              end do
           end do
           ! End loop over new sink particle cells
        end do
        ! End loop over cells
     end do
     ! End loop over grids
  end if

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(oksink_new,oksink_all,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(idsink_new,idsink_all,nsinkmax,MPI_INTEGER         ,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(msink_new ,msink_all ,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(mseed_new ,mseed_all ,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(tsink_new ,tsink_all ,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(xsink_new ,xsink_all ,nsinkmax*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(vsink_new ,vsink_all ,nsinkmax*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(delta_mass_new,delta_mass_all,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(new_born_new,new_born_all,nsinkmax,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,info)
#else
  oksink_all=oksink_new
  idsink_all=idsink_new
  msink_all=msink_new
  mseed_all=mseed_new
  tsink_all=tsink_new
  xsink_all=xsink_new
  vsink_all=vsink_new
  delta_mass_all=delta_mass_new
  new_born_all=new_born_new
#endif
  do isink=1,nsink
     if(oksink_all(isink)==1)then
        idsink(isink)=idsink_all(isink)
        msink(isink)=msink_all(isink)
        mseed(isink)=mseed_all(isink)
        tsink(isink)=tsink_all(isink)
        xsink(isink,1:ndim)=xsink_all(isink,1:ndim)
        vsink(isink,1:ndim)=vsink_all(isink,1:ndim)
        delta_mass(isink)=delta_mass_all(isink)
        acc_rate(isink)=msink_all(isink)
        new_born(isink)=new_born_all(isink)
     endif
  end do
#endif

end subroutine make_sink_from_clump
!################################################################
!################################################################
!################################################################
!################################################################
subroutine true_max(x,y,z,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
  real(dp)::x,y,z
  integer::ilevel

  !----------------------------------------------------------------------------
  ! Description: This subroutine takes the cell of maximum density and computes
  ! the true maximum by expanding the density around the cell center to second order.
  !----------------------------------------------------------------------------

  integer::k,j,i,nx_loc
  integer,dimension(1:nvector)::cell_index,cell_lev
  real(dp)::det,dx,dx_loc,scale,disp_max
  real(dp),dimension(-1:1,-1:1,-1:1)::cube3
  real(dp),dimension(1:nvector,1:ndim)::xtest
  real(dp),dimension(1:ndim)::gradient,displacement
  real(dp),dimension(1:ndim,1:ndim)::hess,minor

#if NDIM==3

  dx=0.5D0**ilevel 
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale


  do i=-1,1
     do j=-1,1
        do k=-1,1

           xtest(1,1)=x+i*dx_loc
#if NDIM>1
           xtest(1,2)=y+j*dx_loc
#endif
#if NDIM>2
           xtest(1,3)=z+k*dx_loc
#endif
           call get_cell_index(cell_index,cell_lev,xtest,ilevel,1)
           cube3(i,j,k)=uold(cell_index(1),1)

        end do
     end do
  end do

! compute gradient
  gradient(1)=0.5*(cube3(1,0,0)-cube3(-1,0,0))/dx_loc
  gradient(2)=0.5*(cube3(0,1,0)-cube3(0,-1,0))/dx_loc
  gradient(3)=0.5*(cube3(0,0,1)-cube3(0,0,-1))/dx_loc

  if (maxval(abs(gradient(1:3)))==0.)return  

  ! compute hessian
  hess(1,1)=(cube3(1,0,0)+cube3(-1,0,0)-2*cube3(0,0,0))/dx_loc**2.
  hess(2,2)=(cube3(0,1,0)+cube3(0,-1,0)-2*cube3(0,0,0))/dx_loc**2.
  hess(3,3)=(cube3(0,0,1)+cube3(0,0,-1)-2*cube3(0,0,0))/dx_loc**2.
  
  hess(1,2)=0.25*(cube3(1,1,0)+cube3(-1,-1,0)-cube3(1,-1,0)-cube3(-1,1,0))/dx_loc**2.
  hess(2,1)=hess(1,2)
  hess(1,3)=0.25*(cube3(1,0,1)+cube3(-1,0,-1)-cube3(1,0,-1)-cube3(-1,0,1))/dx_loc**2.
  hess(3,1)=hess(1,3)
  hess(2,3)=0.25*(cube3(0,1,1)+cube3(0,-1,-1)-cube3(0,1,-1)-cube3(0,-1,1))/dx_loc**2.
  hess(3,2)=hess(2,3)

  !determinant
  det=hess(1,1)*hess(2,2)*hess(3,3)+hess(1,2)*hess(2,3)*hess(3,1)+hess(1,3)*hess(2,1)*hess(3,2) &
       -hess(1,1)*hess(2,3)*hess(3,2)-hess(1,2)*hess(2,1)*hess(3,3)-hess(1,3)*hess(2,2)*hess(3,1)

  !matrix of minors
  minor(1,1)=hess(2,2)*hess(3,3)-hess(2,3)*hess(3,2)
  minor(2,2)=hess(1,1)*hess(3,3)-hess(1,3)*hess(3,1)
  minor(3,3)=hess(1,1)*hess(2,2)-hess(1,2)*hess(2,1)

  minor(1,2)=-1.*(hess(2,1)*hess(3,3)-hess(2,3)*hess(3,1))
  minor(2,1)=minor(1,2)
  minor(1,3)=hess(2,1)*hess(3,2)-hess(2,2)*hess(3,1)
  minor(3,1)=minor(1,3)
  minor(2,3)=-1.*(hess(1,1)*hess(3,2)-hess(1,2)*hess(3,1))
  minor(3,2)=minor(2,3)


  !displacement of the true max from the cell center
  displacement=0.
  do i=1,3
     do j=1,3
        displacement(i)=displacement(i)-minor(i,j)/(det+1.d10*tiny(0.d0))*gradient(j)
     end do
  end do

  !clipping the displacement in order to keep max in the cell
  disp_max=maxval(abs(displacement(1:3)))
  if (disp_max > dx_loc*0.5)then
     displacement(1)=displacement(1)/disp_max*dx_loc*0.5
     displacement(2)=displacement(2)/disp_max*dx_loc*0.5
     displacement(3)=displacement(3)/disp_max*dx_loc*0.5
  end if

  x=x+displacement(1)
  y=y+displacement(2)
  z=z+displacement(3)

#endif
end subroutine true_max
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine update_sink(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This routine is called at the leafs of the tree structure (right after    
! update time). Here is where the global sink variables vsink and xsink are 
! updated by summing the conributions from all levels.                      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(dp)::dteff,vnorm_rel,mach,alpha,factor
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  integer::lev,isink
  integer::i,nx_loc
  real(dp),dimension(1:ndim)::vrel,vcom
  real(dp)::scale_m,boost,beta
  real(dp)::factG=1,fa_fact,v_bondi
  real(dp)::r2,v2,c2,density,volume,ethermal,dx_min,scale,mgas,rho_inf,divergence
  real(dp),dimension(1:3)::velocity


#if NDIM==3

  if(verbose)write(*,*)'Entering update_sink for level ',ilevel

  fsink=0.
  call f_sink_sink

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_m=scale_d*scale_l**3d0
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax/aexp

  vsold(1:nsink,1:ndim,ilevel)=vsnew(1:nsink,1:ndim,ilevel)
  vsnew(1:nsink,1:ndim,ilevel)=vsink(1:nsink,1:ndim)

  ! Loop over sinks
  do isink=1,nsink

     ! sum force contributions from all levels and gather 
     do lev=levelmin,nlevelmax
        fsink(isink,1:ndim)=fsink(isink,1:ndim)+fsink_partial(isink,1:ndim,lev)
     end do
     if (.not. direct_force_sink(isink))then
        fsink(isink,1:ndim)=fsink(isink,1:ndim)/dble(ncloud_sink)
     end if
     
     ! compute timestep for the synchronization
     if (new_born(isink))then
        ! no sync necessary for newly produced sink
        dteff=0d0 
     else
        if(sinkint_level>ilevel)then
           ! level at which sinks are integrated has increased
           ! newdt_fine at coarser level has not been computed 
           ! -> use dtnew from coarse level
           dteff=dtnew(sinkint_level)
        else 
           ! normal case: finish timestep using dtold from current level.
           dteff=dtold(sinkint_level)
        end if
     end if
     
     ! this is the kick-kick (half old half new timestep)
     vsink(isink,1:ndim)=0.5D0*(dtnew(ilevel)+dteff)*fsink(isink,1:ndim)+vsink(isink,1:ndim)

     ! save the velocity
     vsnew(isink,1:ndim,ilevel)=vsink(isink,1:ndim)

     ! this is the kick-kick (half old half new timestep)
     vsink(isink,1:ndim)=0.5D0*(dtnew(ilevel)+dteff)*fsink(isink,1:ndim)+vsink(isink,1:ndim)

     ! save the velocity
     vsnew(isink,1:ndim,ilevel)=vsink(isink,1:ndim)
     
     ! and this is the drift (only for the global sink variable)
     xsink(isink,1:ndim)=xsink(isink,1:ndim)+vsink(isink,1:ndim)*dtnew(ilevel)
     new_born(isink)=.false.

  end do
  ! End loop over sinks

  ! Store deepest level
  sinkint_level=ilevel  

#endif

end subroutine update_sink
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine update_cloud(ilevel)
  use amr_commons
  use pm_commons
  implicit none
  integer::ilevel

  !----------------------------------------------------------------------
  ! update sink cloud particle properties
  ! -the particles are moved whenever the level of the grid they sit in is updated
  ! -the amount of drift they get is according to their levelp
  ! -since this is happening on the way down, at level ilevel all particles with
  ! level >= ilevel will be moved. Therefore, the sink_jump for all levels >= ilevel
  ! is set ot zero on exit.
  !----------------------------------------------------------------------

  integer::igrid,jgrid,ipart,jpart,next_part,ig,ip,npart1,isink,nx_loc
  integer,dimension(1:nvector)::ind_grid,ind_part,ind_grid_part
  real(dp)::dx,dx_loc,scale

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Mesh spacing in that level
  dx=0.5D0**ilevel 
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  ! Update particles position and velocity
  ig=0
  ip=0
  ! Loop over grids 
  igrid=headl(myid,ilevel)
  do jgrid=1,numbl(myid,ilevel)
     npart1=numbp(igrid)  ! Number of particles in the grid 
     if(npart1>0)then
        ig=ig+1
        ind_grid(ig)=igrid
        ipart=headp(igrid)
        ! Loop over particles
        do jpart=1,npart1
           ! Save next particle  <---- Very important !!!
           next_part=nextp(ipart) !move only particles which do actually belong to that level
           if(ig==0)then
              ig=1
              ind_grid(ig)=igrid
           end if
           ip=ip+1
           ind_part(ip)=ipart
           ind_grid_part(ip)=ig
           if(ip==nvector)then
              call upd_cloud(ind_part,ip)
              ip=0
              ig=0
           end if
           ipart=next_part  ! Go to next particle
        end do
        ! End loop over particles
     end if
     igrid=next(igrid)   ! Go to next grid
  end do
  ! End loop over grids
  if(ip>0)call upd_cloud(ind_part,ip)
  
  if (myid==1.and.verbose)then
     write(*,*)'sink drift due to accretion relative to grid size at level ',ilevel
     do isink=1,nsink
        write(*,*)'#sink: ',isink,' drift: ',sink_jump(isink,1:ndim,ilevel)/dx_loc
     end do
  end if

  sink_jump(1:nsink,1:ndim,ilevel:nlevelmax)=0.d0

111 format('   Entering update_cloud for level ',I2)

end subroutine update_cloud
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine upd_cloud(ind_part,np)
  use amr_commons
  use pm_commons
  use poisson_commons
  implicit none
  integer::np
  integer,dimension(1:nvector)::ind_part

  !------------------------------------------------------------
  ! Vector loop called by update_cloud
  !------------------------------------------------------------

  integer::j,idim,isink,lev
  real(dp),dimension(1:nvector,1:ndim)::new_xp,new_vp
  integer,dimension(1:nvector)::level_p

  ! Overwrite cloud particle mass with sink mass
  do j=1,np
     isink=-idp(ind_part(j))
     if(isink>0 .and. mp(ind_part(j))>0.)then
        mp(ind_part(j))=msink(isink)/dble(ncloud_sink_massive)
     endif
  end do

  ! store velocity 
  do idim=1,ndim
     do j=1,np
        new_vp(j,idim)=vp(ind_part(j),idim)
     end do
  end do

  ! Overwrite cloud particle velocity with sink velocity  
  ! is going to be overwritten again before move
  do idim=1,ndim
     do j=1,np
        isink=-idp(ind_part(j))
        if(isink>0)then
              new_vp(j,idim)=vsink(isink,idim)
        end if
     end do
  end do

  ! write back velocity
  do idim=1,ndim
     do j=1,np
        vp(ind_part(j),idim)=new_vp(j,idim)
     end do
  end do

  ! read level
  do j=1,np
     level_p(j)=levelp(ind_part(j))
  end do
 
  ! Update position
  do idim=1,ndim
     do j=1,np
        new_xp(j,idim)=xp(ind_part(j),idim)
     end do
  end do
  do idim=1,ndim
     do j=1,np
        isink=-idp(ind_part(j))
        if(isink>0)then
           lev=level_p(j)
           new_xp(j,idim)=new_xp(j,idim)+sink_jump(isink,idim,lev)
        endif
     end do
  end do
 
 ! Write back postion
  do idim=1,ndim
     do j=1,np
        xp(ind_part(j),idim)=new_xp(j,idim)
     end do
  end do

end subroutine upd_cloud
!################################################################
!################################################################
!################################################################
!################################################################
subroutine merge_star_sink
  use pm_commons
  use amr_commons
  implicit none

  !------------------------------------------------------------------------
  ! This routine merges sink particles for the star formation case if 
  ! they are too close and one of them is younger than ~1000 years
  !------------------------------------------------------------------------

  integer::j,isink,jsink,i,nx_loc,mergers
  real(dp)::dx_loc,scale,dx_min,rr,rmax2,rmax,mnew,t_larson1
  logical::iyoung,jyoung,merge
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v

  if(nsink==0)return

  ! Mesh spacing in that level
  dx_loc=0.5D0**nlevelmax
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax/aexp
  rmax=dble(ir_cloud)*dx_min ! Linking length in physical units
  rmax2=rmax*rmax

  !lifetime of first larson core in code units 
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  t_larson1=merging_timescale*365.25*24*3600/scale_t

  mergers=0
  !loop over all possible pairs (n-square, problematic when there are zilions of sinks)
  do isink=1,nsink-1
     if (msink(isink)>-1.)then
        do jsink=isink+1,nsink
           iyoung=(t-tsink(isink)<t_larson1)
           jyoung=(t-tsink(jsink)<t_larson1)

           rr=(xsink(isink,1)-xsink(jsink,1))**2&
                +(xsink(isink,2)-xsink(jsink,2))**2&
                +(xsink(isink,3)-xsink(jsink,3))**2
           
           merge=(iyoung .or. jyoung).and.rr<rmax2
           merge=merge .or. (iyoung .and. jyoung .and. rr<4*rmax2)
           merge=merge .and. msink(jsink)>=0
           

           if (merge)then
              if (myid==1)write(*,*)'merged ', idsink(jsink),' to ',idsink(isink)
              mergers=mergers+1
              mnew=msink(isink)+msink(jsink)
              xsink(isink,1:3)=(xsink(isink,1:3)*msink(isink)+xsink(jsink,1:3)*msink(jsink))/mnew
              vsink(isink,1:3)=(vsink(isink,1:3)*msink(isink)+vsink(jsink,1:3)*msink(jsink))/mnew
              !wrong! Angular momentum is most likely dominated by last merger, change some day...
              lsink(isink,1:3)=lsink(isink,1:3)+lsink(jsink,1:3)
              msink(isink)=mnew

              acc_rate(isink)=acc_rate(isink)+msink(jsink)
              !acc_lum is recomputed anyway before its used
              if(ir_feedback)acc_lum(isink)=ir_eff*acc_rate(isink)/dtnew(levelmin)*msink(isink)/(5*6.955d10/scale_l)

              msink(jsink)=-10.
              tsink(isink)=min(tsink(isink),tsink(jsink))
              idsink(isink)=min(idsink(isink),idsink(jsink))
           endif
        end do
     end if
  end do

  if (myid==1 .and. mergers>0)write(*,*)'merged ',mergers,' sinks'

  i=1
  do while (mergers>0)

     if (msink(i)<-1.)then !if sink has been merged to another one        

        mergers=mergers-1
        nsink=nsink-1

        !let them all slide back one index
        do j=i,nsink
           xsink(j,1:3)=xsink(j+1,1:3)
           vsink(j,1:3)=vsink(j+1,1:3)
           lsink(j,1:3)=lsink(j+1,1:3)
           msink(j)=msink(j+1)
           tsink(j)=tsink(j+1)
           idsink(j)=idsink(j+1)
           acc_rate(j)=acc_rate(j+1)
           acc_lum(j)=acc_lum(j+1)
        end do

        !whipe last position in the sink list
        xsink(nsink+1,1:3)=0.
        vsink(nsink+1,1:3)=0.
        lsink(nsink+1,1:3)=0.
        msink(nsink+1)=0.
        tsink(nsink+1)=0.
        idsink(nsink+1)=0
        acc_rate(nsink+1)=0.
        acc_lum(nsink+1)=0.

     else
        i=i+1
     end if
  end do

end subroutine merge_star_sink
!################################################################
!################################################################
!################################################################
!################################################################
subroutine merge_smbh_sink
  use pm_commons
  use amr_commons
  use hydro_commons,ONLY:mass_sph
  implicit none

  !------------------------------------------------------------------------
  ! This routine merges sink particles for the smbh formation case if 
  ! they are too close
  !------------------------------------------------------------------------

  integer::j,isink,jsink,i,nx_loc,mergers
  real(dp)::dx_loc,scale,dx_min,rr,rmax2,rmax,mnew,v1_v2,factG
  logical::merge
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v

  if(nsink==0)return

  ! Mesh spacing in that level
  dx_loc=0.5D0**nlevelmax
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax/aexp
  rmax=dble(ir_cloud)*dx_min ! Linking length in physical units
  rmax2=rmax*rmax

  ! Gravitational constant
  factG=1d0
  if(cosmo)factG=3d0/8d0/3.1415926*omega_m*aexp

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  mergers=0
  !loop over all possible pairs (n-square, problematic when there are zilions of sinks)
  do isink=1,nsink-1
     if (msink(isink)>-1)then
        do jsink=isink+1,nsink
           
           ! spacing check
           rr=(xsink(isink,1)-xsink(jsink,1))**2&
                +(xsink(isink,2)-xsink(jsink,2))**2&
                +(xsink(isink,3)-xsink(jsink,3))**2
           
           merge=rr<4*ssoft**2
           merge=merge .and. msink(jsink)>0
           
           ! escape velocity check
           if (mass_vel_check>0)then
              if((msink(isink)+msink(jsink)).ge.mass_vel_check*mass_sph*m_refine(nlevelmax)) then
                 v1_v2=(vsink(isink,1)-vsink(jsink,1))**2+(vsink(isink,2)-vsink(jsink,2))**2+(vsink(isink,3)-vsink(jsink,3))**2
                 merge=merge .and. 2*factG*(msink(isink)+msink(jsink))/sqrt(rr)>v1_v2
                 !!if(myid==1)write(*,*)'DB',idsink(isink),idsink(jsink),msink(isink),msink(jsink),rr,v1_v2
                 !!if(merge.and.myid==1)write(*,*)'MASS VEL CHECK MERGER'
              end if
           end if


           if (merge)then
              if (myid==1)write(*,*)'merged ', idsink(jsink),' to ',idsink(isink)
              mergers=mergers+1
              mnew=msink(isink)+msink(jsink)
              xsink(isink,1:3)=(xsink(isink,1:3)*msink(isink)+xsink(jsink,1:3)*msink(jsink))/mnew
              vsink(isink,1:3)=(vsink(isink,1:3)*msink(isink)+vsink(jsink,1:3)*msink(jsink))/mnew
              !wrong! Angular momentum is most likely dominated by last merger, change some day...
              lsink(isink,1:3)=lsink(isink,1:3)+lsink(jsink,1:3)
              msink(isink)=mnew

              acc_rate(isink)=acc_rate(isink)+acc_rate(jsink)
              delta_mass(isink)=delta_mass(isink)+delta_mass(jsink)
              !acc_lum is recomputed anyway before its used
              if(ir_feedback)acc_lum(isink)=ir_eff*acc_rate(isink)/dtnew(levelmin)*msink(isink)/(5*6.955d10/scale_l)

              msink(jsink)=-10.
              tsink(isink)=min(tsink(isink),tsink(jsink))
              idsink(isink)=min(idsink(isink),idsink(jsink))
           endif
        end do
     end if
  end do

  if (myid==1 .and. mergers>0)write(*,*)'merged ',mergers,' sinks'

  i=1
  do while (mergers>0)

     if (msink(i)<-1.)then !if sink has been merged to another one        

        mergers=mergers-1
        nsink=nsink-1

        !let them all slide back one index
        do j=i,nsink
           xsink(j,1:3)=xsink(j+1,1:3)
           vsink(j,1:3)=vsink(j+1,1:3)
           lsink(j,1:3)=lsink(j+1,1:3)
           msink(j)=msink(j+1)
           tsink(j)=tsink(j+1)
           idsink(j)=idsink(j+1)
           acc_rate(j)=acc_rate(j+1)
           delta_mass(j)=delta_mass(j+1)
           acc_lum(j)=acc_lum(j+1)
        end do

        !whipe last position in the sink list
        xsink(nsink+1,1:3)=0.
        vsink(nsink+1,1:3)=0.
        lsink(nsink+1,1:3)=0.
        msink(nsink+1)=0.
        tsink(nsink+1)=0.
        idsink(nsink+1)=0
        acc_rate(nsink+1)=0.
        delta_mass(nsink+1)=0.
        acc_lum(nsink+1)=0.

     else
        i=i+1
     end if
  end do

end subroutine merge_smbh_sink
!################################################################
!################################################################
!################################################################
!################################################################
subroutine f_gas_sink(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif 
 integer::ilevel
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  ! In this subroutine the sink-gas force contributions are calculated.  
  ! A plummer-sphere with radius ssoft is used for softening   
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  integer::igrid,ngrid,ncache,i,ind,iskip,ix,iy,iz,isink
  integer::info,nx_loc,idim
  real(dp)::dx,dx_loc,scale,vol_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:3)::skip_loc
  
  logical ,dimension(1:nvector)::ok
  integer ,dimension(1:nvector)::ind_grid,ind_cell
  real(dp),dimension(1:nvector,1:ndim)::xx,ff
  real(dp),dimension(1:nvector)::d2,mcell,denom
  real(dp)::rho_tff,rho_tff_tot,d_min
  logical,dimension(1:ndim)::period

  !  Cell spacing at that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim

  ! Set position of cell centers relative to grid centre
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  rho_tff=0.

  fsink_new=0.

  period(1)=(nx==1)
#if NDIM>1
  if(ndim>1)period(2)=(ny==1)
#endif
#if NDIM>2
  if(ndim>2)period(3)=(nz==1)
#endif
  
  ! Loop over sinks 
  do isink=1,nsink
     if (direct_force_sink(isink))then

        d_min=boxlen

        ! Loop over myid grids by vector sweeps
        ncache=active(ilevel)%ngrid
        do igrid=1,ncache,nvector
           ngrid=MIN(nvector,ncache-igrid+1)
           do i=1,ngrid
              ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
           end do

           ! Loop over cells
           do ind=1,twotondim
              iskip=ncoarse+(ind-1)*ngridmax
              do i=1,ngrid
                 ind_cell(i)=iskip+ind_grid(i)
              end do

              ! Check if cell is refined
              do i=1,ngrid
                 ok(i)=son(ind_cell(i))==0
              end do

              ! Gas and dark matter mass in cell
              do i=1,ngrid
                 mcell(i)=rho(ind_cell(i))*vol_loc
              end do

              ! Cell center
              do idim=1,ndim
                 do i=1,ngrid
                    xx(i,idim)=(xg(ind_grid(i),idim)+xc(ind,idim)-skip_loc(idim))*scale
                 end do
              end do

              ! Relative position and distance
              d2=0.d0
              do idim=1,ndim    
                 if (period(idim))then
                    do i=1,ngrid    
                       ! zero order Ewald sum
                       ff(i,idim)=xsink(isink,idim)-xx(i,idim)
                       if(ff(i,idim)>0.5*boxlen)ff(i,idim)=ff(i,idim)-boxlen
                       if(ff(i,idim)<-0.5*boxlen)ff(i,idim)=ff(i,idim)+boxlen
                       d2(i)=d2(i)+ff(i,idim)**2
                    end do
                 else
                    do i=1,ngrid
                       ff(i,idim)=xsink(isink,idim)-xx(i,idim)
                       d2(i)=d2(i)+ff(i,idim)**2
                    end do
                 end if
              end do

              ! Store minimum distance of cell in current level to isink
              do i=1,ngrid
                 d_min=min(d_min,d2(i))
              end do

              ! Compute sqrt(1/(ssoft**2+d2(i))) to save time
              do i=1,ngrid
                 denom(i)=(ssoft**2+d2(i))**(-1.5)
              end do

              ! Compute gas acceleration due to sink
              do i=1,ngrid
                 ff(i,1:ndim)=denom(i)*ff(i,1:ndim)
              end do

              ! Add gas acceleration due to sink
              do i=1,ngrid
                 f(ind_cell(i),1:ndim)=f(ind_cell(i),1:ndim)+msink(isink)*ff(i,1:ndim)
              end do

              ! Add sink acceleration due to gas
              do i=1,ngrid
                 if(ok(i))then
                    fsink_new(isink,1:ndim)=fsink_new(isink,1:ndim)-mcell(i)*ff(i,1:ndim)
                 end if
              end do
           end do !end loop over cells
        end do !end loop over grids
        
        d_min=d_min**0.5
        d_min=max(ssoft,d_min)
        rho_tff=max(rho_tff,msink(isink)/(4./3.*4.13145*d_min**3))

     end if !end if direct force
  end do !end loop over sinks

  !collect sink acceleration from cpus
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(fsink_new,fsink_all,nsinkmax*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#else
     fsink_all=fsink_new
#endif
  do isink=1,nsink
     if (direct_force_sink(isink))then
        fsink_partial(isink,1:ndim,ilevel)=fsink_all(isink,1:ndim)
     end if
  end do  

  do idim=1,ndim
     call make_virtual_fine_dp(f(1,idim),ilevel)
  end do
  
  
  !collect rho due to sinks for current level - used for timestep computation
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(rho_tff,rho_tff_tot,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
#else
  rho_tff_tot=rho_tff
#endif
  rho_sink_tff(ilevel)=rho_tff_tot
  
  
  
  if (ilevel==nlevelmax)call make_virtual_fine_dp(phi(1),ilevel)
  
end subroutine f_gas_sink
   
!################################################################
!################################################################
!################################################################
!################################################################
subroutine f_sink_sink
  use amr_commons
  use pm_commons  
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! In this subroutine the sink-sink force contribution are calculated by direct
  ! n^2 - summation. A plummer-sphere with radius 4 cells is used for softening
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer::isink,idim,jsink
  real(dp),allocatable,dimension(:)::d2
  real(dp),allocatable,dimension(:,:)::ff
  logical,dimension(1:ndim)::period

  allocate(d2(1:nsink))
  allocate(ff(1:nsink,1:ndim))
  
!  fsink=0.

  period(1)=(nx==1)
#if NDIM>1
  if(ndim>1)period(2)=(ny==1)
#endif
#if NDIM>2
  if(ndim>2)period(3)=(nz==1)  
#endif
  
  do isink=1,nsink
     if (direct_force_sink(isink))then
        d2=0.d0
        ff=0.d0
        do idim=1,ndim
           !compute relative position and and distances
           if (period(idim))then
              do jsink=1,nsink
                 if (direct_force_sink(jsink))then           
                    ff(jsink,idim)=xsink(jsink,idim)-xsink(isink,idim)
                    if(ff(jsink,idim)>0.5*boxlen)ff(jsink,idim)=ff(jsink,idim)-boxlen
                    if(ff(jsink,idim)<-0.5*boxlen)ff(jsink,idim)=ff(jsink,idim)+boxlen
                    d2(jsink)=d2(jsink)+ff(jsink,idim)**2
                 end if
              end do
           else
              do jsink=1,nsink
                 if (direct_force_sink(jsink))then           
                    ff(jsink,idim)=xsink(jsink,idim)-xsink(isink,idim)
                    d2(jsink)=d2(jsink)+ff(jsink,idim)**2
                 end if
              end do
           end if
        end do
        !compute acceleration
        do jsink=1,nsink
           if (direct_force_sink(jsink))then
              ff(jsink,1:ndim)=msink(jsink)/(ssoft**2+d2(jsink))**1.5*ff(jsink,1:ndim)
           end if
        end do
        do jsink=1,nsink           
           fsink(isink,1:ndim)=fsink(isink,1:ndim)+ff(jsink,1:ndim)
           if(jsink<0)then
              print*,'This is just a stupid trick to prevent'
              print*,'the compiler from optimizing this loop!'
           end if
        end do             
     end if
  end do
end subroutine f_sink_sink
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine read_sink_params()
  use pm_commons
  use amr_commons
  implicit none

  !------------------------------------------------------------------------
  ! read sink related parameters and perform some 'sanity chekcs'
  !------------------------------------------------------------------------

  real(dp)::dx_min,scale,cty
  integer::nx_loc
  namelist/sink_params/n_sink,rho_sink,d_sink,accretion_scheme,nol_accretion,merging_scheme,merging_timescale,&
       ir_cloud_massive,sink_soft,msink_direct,ir_cloud,nsinkmax,c_acc,create_sinks,sink_seedmass,&
       eddington_limit,sink_drag,alpha_acc_boost,alpha_drag_boost,acc_threshold_creation,mass_vel_check,&
       clump_core,smbh_verbose
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)  

  nx_loc=(icoarse_max-icoarse_min+1)
  scale = boxlen/dble(nx_loc)

  ! Read namelist file 
  rewind(1)
  read(1,NML=sink_params,END=111)
  goto 112
111 if(myid==1)write(*,*)'You did not set up &SINK_PARAMS in the namelist file'
  if(myid==1)write(*,*)'Using default values '
112 rewind(1)

  if (sink .and. (ndim .ne. 3))then
     if(myid==1)write(*,*)'Sink particles are only implemented for 3d sims.'
     print*, ndim
     call clean_stop
  end if

  if (sink .and. (nlevelmax==levelmin))then
     if(myid==1)write(*,*)'sink particles do currently not work in a single-level cartesian grid'
     if(myid==1)write(*,*)'because they need level 1 to be activated.'
     call clean_stop
  end if
  
  if (create_sinks .and. accretion_scheme=='none')then
     if(myid==1)write(*,*)'formation of new sinks without subsequent accretion is pointless.' 
     if(myid==1)write(*,*)'Choose accretion_scheme=...!'
!     call clean_stop
  end if

  if ((create_sinks .or. accretion_scheme .ne. 'none') .and. (.not. hydro))then
     if(myid==1)write(*,*)'sink creation and accretion require hydro to be turned on'
     call clean_stop
  end if
  
  !check for accretion scheme
  if (accretion_scheme =='flux')flux_accretion=.true.
  if (accretion_scheme =='bondi')bondi_accretion=.true.
  if (accretion_scheme =='threshold')threshold_accretion=.true.

  ! for sink formation and accretion a threshold must be given
  if (create_sinks .or. (accretion_scheme .ne. 'none'))then

     ! check for threshold  
     if (n_sink<0. .and. cosmo)then
        if(myid==1)write(*,*)'specify n_sink for a cosmological simulation'
        call clean_stop
     end if

     if (rho_sink<0. .and. n_sink<0. .and. d_sink>0.) then
        if(myid==1)write(*,*)'Found d_sink! Assuming code units'
     else if (rho_sink>0. .and. n_sink<0. .and. d_sink<0.)then
        if(myid==1)write(*,*)'Found rho_sink! Assuming g/cc'
        d_sink=rho_sink/scale_d
     else if (rho_sink<0. .and. n_sink>0. .and. d_sink<0.)then
        if(myid==1)write(*,*)'Found n_sink! Assuming H/cc'
        d_sink=n_sink/scale_nH
     else if ((rho_sink>0. .and. n_sink>0.) .or. (rho_sink>0. .and. d_sink>0.) .or. (n_sink>0. .and. d_sink>0.))then
        if (myid==1)write(*,*)'Use n_sink [H/cc] OR rho_sink [g/cc] OR d_sink [code_units]'
        call clean_stop
     else
        if(myid==1)write(*,*)'Trying to setting sink threshold such that jeans length at '
        if(myid==1)write(*,*)'max resolution is resolved by 4 cells, assuming isothermal gas'
        if(T2_star==0.)then 
           if(myid==1)write(*,*)'No value for T2_star given. Do not know what to do...'
           call clean_stop
        else
           dx_min=0.5**nlevelmax*scale
           d_sink=T2_star/scale_T2 *3.14159/16./(dx_min**2)
           if(myid==1)write(*,*)'d_sink = ',d_sink
           if(myid==1)write(*,*)'rho_sink = ',d_sink*scale_d
           if(myid==1)write(*,*)'n_sink = ',d_sink*scale_nH
        end if
     end if
  end if  
  
  
  if (merging_scheme == 'timescale')then
     if (merging_timescale<0.)then
        if (myid==1)write(*,*)'You chose sink merging on a timescale but did not provide the timescale'
        if (myid==1)write(*,*)'choosing 1000y as lifetime...'
        merging_timescale=1000.
     end if
     cty=scale_t/(365.25*24.*3600.)
     cont_speed=-1./(merging_timescale/cty)
  end if

  ! a warging sign...
  if (merging_scheme == 'fof' .and. (.not. smbh))then
     write(*,*)'friend of friend merging has been developped for SMBH - experimental!!'
  end if
     
  ! for smbh runs use friends of friends merging
  if(smbh)merging_scheme='fof'

  ! nol_accretion requires a somewhat smaller timestep per default
  if(c_acc < 0.)then
     if (nol_accretion)then
        c_acc=0.25
     else
        c_acc=0.75
     end if
  end if

  !check for periodic boundary conditions
  if (nx==1 .or. ny==1 .or. nz==1)then
     if (msink_direct .ge. 0.)then
        if(myid==1)print*, 'periodic boundaries in combination with '
        if(myid==1)print*, 'direct force sinks are not treated accurately....'
     end if
  end if

  if(msink_direct<0.)then 
     msink_direct=huge(0._dp)
  else
     !convert msink_direct in code units
     msink_direct=msink_direct*1.9891d33/(scale_d*scale_l**3)
  end if

end subroutine read_sink_params
!################################################################
!################################################################
!################################################################
!################################################################
subroutine count_clouds(ilevel,action)
  use pm_commons
  use amr_commons
  use poisson_commons, only:rho
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  character(len=15)::action

  !------------------------------------------------------------------------
  ! loop over all sink cloud particles and pass nvector parts to count_clouds_np
  ! first pass (action 'count'): count the number of sink cloud parts per cell
  ! and store them in flag2
  ! second pass (action 'weight'): assign a weight corresponding to the physical
  ! volume to each cloud particle. reduce if sinks overlap (flag2 too big).
  !------------------------------------------------------------------------

  integer::igrid,jgrid,ipart,jpart,next_part
  integer::ig,ip,npart1,npart2,icpu,ind,cell_index
  integer,dimension(1:nvector)::ind_grid,ind_part,ind_grid_part
  
  if(numbtot(1,ilevel)==0)return


  if (action=='count')then
     ! Loop over cpus
     do icpu=1,ncpu
        igrid=headl(icpu,ilevel)
        ! Loop over grids
        do jgrid=1,numbl(icpu,ilevel)
           !set flag2 to zero for all cells in the grid
           do ind=1,twotondim
              cell_index=ncoarse+(ind-1)*ngridmax+igrid
              rho(cell_index)=0.
           end do
           igrid=next(igrid)   ! Go to next grid
        end do
     end do
  end if
  

  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0

        ! Count sink and cloud particles
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if(idp(ipart).lt.0)then
                 npart2=npart2+1
              endif
              ipart=next_part  ! Go to next particle
           end do
        endif
        
        ! Gather sink and cloud particles
        if(npart2>0)then        
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              ! Select only sink particles
              if(idp(ipart).lt.0)then
                 if(ig==0)then
                    ig=1
                    ind_grid(ig)=igrid
                 end if
                 ip=ip+1
                 ind_part(ip)=ipart
                 ind_grid_part(ip)=ig   
              endif
              if(ip==nvector)then
!                 call count_clouds_np(ind_part,ip,action,ilevel)
                 call count_clouds_np(ind_grid,ind_part,ind_grid_part,ig,ip,action,ilevel)
                 ip=0
                 ig=0
              end if
              ipart=next_part  ! Go to next particle
           end do
           ! End loop over particles
        end if
        igrid=next(igrid)   ! Go to next grid
     end do

     ! End loop over grids
     if(ip>0)then
!        call count_clouds_np(ind_part,ip,action,ilevel)
        call count_clouds_np(ind_grid,ind_part,ind_grid_part,ig,ip,action,ilevel)
    end if
  end do
  ! End loop over cpus

end subroutine count_clouds
!################################################################
!################################################################
!################################################################
!################################################################
subroutine count_clouds_np(ind_grid,ind_part,ind_grid_part,ng,np,action,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons, only:uold
  use poisson_commons, only:rho
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_part,ind_grid,ind_grid_part
  character(len=15)::action
  !--------------------------------------------------------------------------------
  ! inner loop of count_clouds
  !--------------------------------------------------------------------------------  
  integer::i,nx_loc,isink,idim,ind
  real(dp)::scale,dx_min,weight,parts_per_cell,vol_min
  integer,dimension(1:nvector),save::clevl
  integer,dimension(1:nvector,1:twotondim),save::cind_part
  real(dp),dimension(1:nvector,1:ndim),save::xpart
  real(dp),dimension(1:nvector,1:ndim,1:twotondim),save::xx
  real(dp),dimension(1:nvector,1:twotondim),save::vol
  logical,dimension(1:nvector,twotondim)::ok

  ! compute volume of smallest cell
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_min=(0.5D0**nlevelmax)*scale
  vol_min=dx_min**ndim

  do idim=1,ndim
     do i=1,np
        xpart(i,idim)=xp(ind_part(i),idim)
     end do
  end do

!  call get_cell_index(cind_part,clevl,xpart,nlevelmax,np)
!  call get_cell_index_for_particle(cind_part,xx,clevl,ind_grid,xpart,ind_grid_part,ng,np,ilevel,ok)
  call cic_get_cells(cind_part,xx,vol,ok,ind_grid,xpart,ind_grid_part,ng,np,ilevel)

!  do i=1,np
!     if (cind_part(i).ne.cind_part2(i))print*,'oups at count clouds'
!  end do




!  do i=1,np
!     if (clevl(i) .ne. ilevel)then 
!        write(*,*)'particle outside its level'
!        write(*,*)'this should not happen'
!     end if
!  end do

  !only particles "in their level" accrete
  do ind=1,twotondim
     do i=1,np
        if(.not. ok(i,ind))then
           vol(i,ind)=0.
        end if
     end do
  end do
  
  if (action=='count')then
     do ind=1,twotondim
        do i=1,np
           rho(cind_part(i,ind))=rho(cind_part(i,ind))+vol(i,ind)
        end do
     end do
  end if
  
  ! weight each particle with its actual volume (reduced in the case of overlapping sinks)  
  if (action=='weight')then     
     do ind=1,twotondim
        do i=1,np
           isink=-idp(ind_part(i))
           parts_per_cell=8*8.**(nlevelmax-ilevel)
           weight=vol_min*0.125*vol(i,ind)
           if (rho(cind_part(i,ind))>parts_per_cell)then
              weight=weight*parts_per_cell/rho(cind_part(i,ind))
           end if
           weightp(ind_part(i),ind)=weight
        end do
        
        if(threshold_accretion)then
           do i=1,np
              if(uold(cind_part(i,ind),1)<d_sink)then
                 weightp(ind_part(i),ind)=0.
              end if
           end do
        end if
     end do
  end if
  
end subroutine count_clouds_np
!################################################################
!################################################################
!################################################################
!################################################################
subroutine get_cell_index_for_particle(indp,xx,cell_lev,ind_grid,xpart,ind_grid_part,ng,np,ilevel,ok)
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid,indp,cell_lev,ind_grid_part
  real(dp),dimension(1:nvector,1:ndim)::xpart,xx
  logical,dimension(1:nvector)::ok


  !-----------------------------------------------------------------------
  ! This subroutine finds the leaf cell in which a particle sits 
  !-----------------------------------------------------------------------

  integer::i,j,idim,nx_loc,ind,ix,iy,iz
  real(dp)::dx,dx_loc,scale,one_over_dx,one_over_scale
  ! Grid based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle based arrays

  real(dp),dimension(1:nvector,1:ndim),save::x
  integer ,dimension(1:nvector,1:ndim),save::id,igd,icd,icd_fine
  integer ,dimension(1:nvector),save::igrid,icell,kg,icell_fine
  real(dp),dimension(1:3),save::skip_loc
  real(dp),dimension(1:twotondim,1:3),save::xc
 

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  one_over_dx=1./dx
  one_over_scale=1./scale

  ! Cells center position relative to grid center position
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     xc(ind,1)=(dble(ix)-0.5D0)*dx
     xc(ind,2)=(dble(iy)-0.5D0)*dx
     xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  ! Lower left corner of 3x3x3 grid-cube
  do idim=1,ndim
     do i=1,ng
        x0(i,idim)=xg(ind_grid(i),idim)-3.0D0*dx
     end do
  end do

  ! Gather 27 neighboring father cells (should be present anytime !)
  do i=1,ng
     ind_cell(i)=father(ind_grid(i))
  end do
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ng,ilevel)

  ! Rescale position at level ilevel
  do idim=1,ndim
     do j=1,np
        x(j,idim)=xpart(j,idim)*one_over_scale+skip_loc(idim)
     end do
  end do
  do idim=1,ndim
     do j=1,np
        x(j,idim)=x(j,idim)-x0(ind_grid_part(j),idim)
     end do
  end do
  do idim=1,ndim
     do j=1,np
        x(j,idim)=x(j,idim)*one_over_dx
     end do
  end do

  ! Check for illegal moves
  do idim=1,ndim
     do j=1,np
        if(x(j,idim)<=0.0D0.or.x(j,idim)>=6.0D0)then
           print*,'cpu ', myid, ' hosts an escaped particle'
           print*,'x: ',x(j,1:ndim)
           print*,'x0: ',x0(ind_grid_part(j),1:ndim)
           print*,'xg: ',xg(ind_grid(ind_grid_part(j)),1:ndim)
           print*,'xp: ',xpart(j,1:ndim)
           print*,'skip_loc: ',skip_loc(1:ndim)
           print*,'scale: ',scale
           call clean_stop
        end if
     end do
  end do
  
  ! NGP at level ilevel
  do idim=1,ndim
     do j=1,np
        id(j,idim)=int(x(j,idim))
     end do
  end do

   ! Compute parent grids
  do idim=1,ndim
     do j=1,np
        igd(j,idim)=id(j,idim)/2
     end do
  end do
#if NDIM==1
  do j=1,np
     kg(j)=1+igd(j,1)
  end do
#endif
#if NDIM==2
  do j=1,np
     kg(j)=1+igd(j,1)+3*igd(j,2)
  end do
#endif
#if NDIM==3
  do j=1,np
     kg(j)=1+igd(j,1)+3*igd(j,2)+9*igd(j,3)
  end do
#endif

  do j=1,np
     igrid(j)=son(nbors_father_cells(ind_grid_part(j),kg(j)))
  end do

  ! Check if particle has escaped to ilevel-1
  ok(1:np)=.true.
  do j=1,np
     if (igrid(j)==0)then
!        print*,'particle has escaped to ilevel-1'
        ok(j)=.false.
        indp(j)=nbors_father_cells(ind_grid_part(j),kg(j))
        cell_lev(j)=ilevel-1
        xx(j,1:ndim)=(xg(ind_grid(ind_grid_part(j)),1:ndim)+(igd(j,1:ndim)-1.)*2.*dx-skip_loc(1:ndim))*scale
        if (sum((xx(j,1:ndim)-xpart(j,1:ndim))**2)**0.5>1.000001*dx_loc*3**0.5)print*,'oups at ilevel-1'
     end if
  end do
  
  ! Compute parent cell position
  do idim=1,ndim
     do j=1,np
!        if(ok(j))then
           icd(j,idim)=id(j,idim)-2*igd(j,idim)
!        end if
     end do
  end do
! #if NDIM==1
!   do j=1,np
!      if(ok(j))then
!         icell(j)=1+icd(j,1)
!      end if
!   end do
! #endif
! #if NDIM==2
!   do j=1,np
!      if(ok(j))then
!         icell(j)=1+icd(j,1)+2*icd(j,2)
!      end if
!   end do
! #endif
!#if NDIM==3
!  do j=1,np
!     if(ok(j))then
!        icell(j)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
!     end if
!  end do
!#endif
        
  call geticell(icell,icd,np)
  
  ! Compute parent cell adress
  do j=1,np
     if(ok(j))then
        indp(j)=ncoarse+(icell(j)-1)*ngridmax+igrid(j)
     end if
  end do

  ! Check if particles have leaked into level ilevel+1
  do j=1,np
     if(ok(j))then
        if (son(indp(j))>0)then
!           print*,'particle has escaped to ilevel+1'
           ok(j)=.false.
           cell_lev(j)=ilevel+1
           do idim=1,ndim
              icd_fine(1,idim)=int(2*(x(j,idim)-int(x(j,idim))))
           end do
           call geticell(icell_fine,icd_fine,1)
           xx(j,1:ndim)=(xg(son(indp(j)),1:ndim)+xc(icell_fine(1),1:ndim)*0.5-skip_loc(1:ndim))*scale
           indp(j)=ncoarse+(icell_fine(1)-1)*ngridmax+son(indp(j))
           if (sum((xx(j,1:ndim)-xpart(j,1:ndim))**2)**0.5>1.0000001*0.25*dx_loc*3**0.5)then
              print*,icd_fine(1,1:ndim)
              print*,icell_fine(1)
              print*,'oups at ilevel+1'
           end if
        end if
     endif
  end do

  !cell center positions for particles which sit in the leve ilevel
  do j=1,np
     if (ok(j))then
        xx(j,1:ndim)=(xg(igrid(j),1:ndim)+xc(icell(j),1:ndim)-skip_loc(1:ndim))*scale
        cell_lev(j)=ilevel
        ! if (sum((xx(j,1:ndim)-xpart(j,1:ndim))**2)**0.5>1.000001*0.5*dx_loc*3**0.5)then
        !    print*,'oups at ilevel'        
        !    print*,xpart(j,1:ndim)
        !    print*,xx(j,1:ndim)
        !    print*,sum((xx(j,1:ndim)-xpart(j,1:ndim))**2)**0.5/dx
        ! end if
     end if
  end do

end subroutine get_cell_index_for_particle
!################################################################
!################################################################
!################################################################
!################################################################
subroutine geticell(icell,icd,np)
  use amr_parameters, only:nvector,ndim
  integer::np
  integer,dimension(1:nvector,1:ndim)::icd
  integer,dimension(1:nvector)::icell
  ! mini subroutine that gets the cell index (1 to 8)
  ! for certain coordinates (0,1 along each direction)
  ! put into a subroutine to make the code look less ugly
  integer::j
    
#if NDIM==1
  do j=1,np
     icell(j)=1+icd(j,1)
  end do
#endif
#if NDIM==2
  do j=1,np
     icell(j)=1+icd(j,1)+2*icd(j,2)
  end do
#endif
#if NDIM==3
  do j=1,np
     icell(j)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
  end do
#endif
  
end subroutine geticell
!################################################################
!################################################################
!################################################################
!################################################################
subroutine get_cell_center(xx,index,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
  real(dp),dimension(1:3)::xx
  integer::index,ilevel
  ! little helper routine that gets the cell center of a cell with a given index
  ! used only for debugging...
  integer::ind,ix,iy,iz,i,ind8,grid
  real(dp)::dx,scale
  real(dp),dimension(1:8,1:3)::xc
  real(dp),dimension(1:3)::skip_loc
  integer::nx_loc
  dx=0.5D0**ilevel 
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)


  ! Cells center position relative to grid center position
  do ind=1,8  
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     xc(ind,1)=(dble(ix)-0.5D0)*dx
     xc(ind,2)=(dble(iy)-0.5D0)*dx
     xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do
  
  ind8=(index-ncoarse)/ngridmax

  grid=mod((index-ncoarse),ngridmax)
  print*,'ind8',ind8
  print*,'grid',grid
  print*,'xg',xg(grid,1:3)
  print*,'xc',xc(ind8+1,1:3)


  do i=1,3
     xx(i)=(xg(grid,i)+xc(ind8+1,i)-skip_loc(i))*scale
  end do
  print*,'xx',xx(1:3)
  
end subroutine get_cell_center
!################################################################
!################################################################
!################################################################
!################################################################
subroutine count_parts
  use pm_commons
  use amr_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  ! ugly routine to count all particles in the simulation level by level
  ! used for debugging
  integer::ilevel
  integer::igrid,jgrid,i,ngrid,ncache
  integer::ig,ip,npart1,npart2,npart2_tot,icpu,info
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nlevelmax)::npts

  npts=0  
  do ilevel=1,nlevelmax
     npart2=0
     ! Loop over cpus
     do icpu=1,ncpu
        igrid=headl(icpu,ilevel)
        ig=0
        ip=0
        ! Loop over grids
        do jgrid=1,numbl(icpu,ilevel)
           npart1=numbp(igrid)  ! Number of particles in the grid
           npart2=npart2+npart1
           igrid=next(igrid)   ! Go to next grid                                                              
        
        end do
     end do
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(npart2,npart2_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#else
     npart2_tot=npart2
#endif          
     npts(ilevel)=npart2_tot
  end do
  if (myid==1)print*,'total'
  if (myid==1)print*,npts(1:nlevelmax)
  
  

   npts=0  
   do ilevel=1,nlevelmax
      npart2=0        
      
      ncache=active(ilevel)%ngrid
      do igrid=1,ncache,nvector
         ngrid=MIN(nvector,ncache-igrid+1)
         do i=1,ngrid
            ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
         end do
         do i=1,ngrid
            npart2=npart2+numbp(ind_grid(i))
         end do
        
      end do

#ifndef WITHOUTMPI
      call MPI_ALLREDUCE(npart2,npart2_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#else
      npart2_tot=npart2
#endif          
      npts(ilevel)=npart2_tot
   end do
   if (myid==1)print*,'active'
   if (myid==1)print*,npts(1:nlevelmax)        
   
   
end subroutine count_parts
!################################################################
!################################################################
!################################################################
!################################################################
subroutine cic_get_cells(indp,xx,vol,ok,ind_grid,xpart,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use pm_commons
  implicit none
  integer::ng,np,ilevel
  integer ,dimension(1:nvector)::ind_grid,ind_grid_part
  real(dp),dimension(1:nvector,1:ndim)::xpart
  real(dp),dimension(1:nvector,1:twotondim)::vol
  integer ,dimension(1:nvector,1:twotondim)::indp
  logical ,dimension(1:nvector,1:twotondim)::ok
  real(dp),dimension(1:nvector,1:ndim,twotondim)::xx
  !------------------------------------------------------------------
  ! This routine returns the CIC cells and volumes for np particles.
  !------------------------------------------------------------------
  logical::error
  integer::i,j,ind,idim,nx_loc,ix,iy,iz
  real(dp)::dx,dx_loc,scale,vol_loc
  ! Grid-based arrays
  integer ,dimension(1:nvector)::ind_cell
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle-based arrays
  real(dp),dimension(1:nvector),save::mmm
  real(dp),dimension(1:nvector),save::ttt=0d0
  real(dp),dimension(1:nvector),save::vol2
  real(dp),dimension(1:nvector,1:ndim),save::x,dd,dg
  integer ,dimension(1:nvector,1:ndim),save::ig,id,igg,igd,icg,icd
  integer ,dimension(1:nvector,1:twotondim),save::igrid,icell,kg
  real(dp),dimension(1:nvector,1:ndim),save::x0
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:ndim)::xc

  ! Mesh spacing in that level
  dx=0.5D0**ilevel 
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim

  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     xc(ind,1)=(dble(ix)-0.5D0)*dx
#if NDIM>1
     xc(ind,2)=(dble(iy)-0.5D0)*dx
#endif
#if NDIM>2
     xc(ind,3)=(dble(iz)-0.5D0)*dx
#endif
  end do

  ! Lower left corner of 3x3x3 grid-cube  do idim=1,ndim
  do idim=1,ndim
     do i=1,ng
        x0(i,idim)=xg(ind_grid(i),idim)-3.0D0*dx
     end do
  end do

  do i=1,ng
     ind_cell(i)=father(ind_grid(i))
  end do
  
  ! Gather neighboring father cells (should be present anytime !)
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ng,ilevel)

  ! Rescale particle position at level ilevel
  do idim=1,ndim
     do j=1,np
        x(j,idim)=xpart(j,idim)/scale+skip_loc(idim)
     end do
  end do
  do idim=1,ndim
     do j=1,np
        x(j,idim)=x(j,idim)-x0(ind_grid_part(j),idim)
     end do
  end do
  do idim=1,ndim
     do j=1,np
        x(j,idim)=x(j,idim)/dx
     end do
  end do

  do idim=1,ndim
     do j=1,np
        if (x(j,idim)<0.5 .or. x(j,idim)>5.5)then
           print*,'particle outside allowed boundary for cic_get_cell'
           print*,x(j,1:ndim)
           if (x(j,idim)<0. .or. x(j,idim)>6.)then
              print*,'particle outside allowed 3by3by3 grid cube'
           end if
           call clean_stop
        end if
     end do
  end do
  
  ! CIC at level ilevel (dd: right cloud boundary; dg: left cloud boundary)
  do idim=1,ndim
     do j=1,np
        dd(j,idim)=x(j,idim)+0.5D0
        id(j,idim)=dd(j,idim)
        dd(j,idim)=dd(j,idim)-id(j,idim)
        dg(j,idim)=1.0D0-dd(j,idim)
        ig(j,idim)=id(j,idim)-1
     end do
  end do

  ! Compute cloud volumes
#if NDIM==1
  do j=1,np
     vol(j,1)=dg(j,1)
     vol(j,2)=dd(j,1)
  end do
#endif
#if NDIM==2
  do j=1,np
     vol(j,1)=dg(j,1)*dg(j,2)
     vol(j,2)=dd(j,1)*dg(j,2)
     vol(j,3)=dg(j,1)*dd(j,2)
     vol(j,4)=dd(j,1)*dd(j,2)
  end do
#endif
#if NDIM==3
  do j=1,np
     vol(j,1)=dg(j,1)*dg(j,2)*dg(j,3)
     vol(j,2)=dd(j,1)*dg(j,2)*dg(j,3)
     vol(j,3)=dg(j,1)*dd(j,2)*dg(j,3)
     vol(j,4)=dd(j,1)*dd(j,2)*dg(j,3)
     vol(j,5)=dg(j,1)*dg(j,2)*dd(j,3)
     vol(j,6)=dd(j,1)*dg(j,2)*dd(j,3)
     vol(j,7)=dg(j,1)*dd(j,2)*dd(j,3)
     vol(j,8)=dd(j,1)*dd(j,2)*dd(j,3)
  end do
#endif
        
  ! Compute parent grids
  do idim=1,ndim
     do j=1,np
        igg(j,idim)=ig(j,idim)/2
        igd(j,idim)=id(j,idim)/2
     end do
  end do
#if NDIM==1
  do j=1,np
     kg(j,1)=1+igg(j,1)
     kg(j,2)=1+igd(j,1)
  end do
#endif
#if NDIM==2
  do j=1,np
     kg(j,1)=1+igg(j,1)+3*igg(j,2)
     kg(j,2)=1+igd(j,1)+3*igg(j,2)
     kg(j,3)=1+igg(j,1)+3*igd(j,2)
     kg(j,4)=1+igd(j,1)+3*igd(j,2)
  end do
#endif
#if NDIM==3
  do j=1,np
     kg(j,1)=1+igg(j,1)+3*igg(j,2)+9*igg(j,3)
     kg(j,2)=1+igd(j,1)+3*igg(j,2)+9*igg(j,3)
     kg(j,3)=1+igg(j,1)+3*igd(j,2)+9*igg(j,3)
     kg(j,4)=1+igd(j,1)+3*igd(j,2)+9*igg(j,3)
     kg(j,5)=1+igg(j,1)+3*igg(j,2)+9*igd(j,3)
     kg(j,6)=1+igd(j,1)+3*igg(j,2)+9*igd(j,3)
     kg(j,7)=1+igg(j,1)+3*igd(j,2)+9*igd(j,3)
     kg(j,8)=1+igd(j,1)+3*igd(j,2)+9*igd(j,3)
  end do
#endif
  do ind=1,twotondim
     do j=1,np
        igrid(j,ind)=son(nbors_father_cells(ind_grid_part(j),kg(j,ind)))
     end do
  end do
  
  ! Check if particle has escaped to ilevel-1
  ok(1:np,1:twotondim)=.true.
  do ind=1,twotondim
     do j=1,np
        if (igrid(j,ind)==0)then
           ! print*,'particle has partially escaped to ilevel-1' 
           ok(j,ind)=.false.
           indp(j,ind)=nbors_father_cells(ind_grid_part(j),kg(j,ind))
!           xx(j,1:ndim,ind)=(xg(ind_grid(ind_grid_part(j)),1:ndim)+ind3(kg(j,ind),1:ndim)*2.*dx-skip_loc(1:ndim))*scale
        end if
     end do
  end do


  ! Compute parent cell position
  do idim=1,ndim
     do j=1,np
        icg(j,idim)=ig(j,idim)-2*igg(j,idim)
        icd(j,idim)=id(j,idim)-2*igd(j,idim)
     end do
  end do
#if NDIM==1
  do j=1,np
     icell(j,1)=1+icg(j,1)
     icell(j,2)=1+icd(j,1)
  end do
#endif
#if NDIM==2
  do j=1,np
     icell(j,1)=1+icg(j,1)+2*icg(j,2)
     icell(j,2)=1+icd(j,1)+2*icg(j,2)
     icell(j,3)=1+icg(j,1)+2*icd(j,2)
     icell(j,4)=1+icd(j,1)+2*icd(j,2)
  end do
#endif
#if NDIM==3
  do j=1,np
     icell(j,1)=1+icg(j,1)+2*icg(j,2)+4*icg(j,3)
     icell(j,2)=1+icd(j,1)+2*icg(j,2)+4*icg(j,3)
     icell(j,3)=1+icg(j,1)+2*icd(j,2)+4*icg(j,3)
     icell(j,4)=1+icd(j,1)+2*icd(j,2)+4*icg(j,3)
     icell(j,5)=1+icg(j,1)+2*icg(j,2)+4*icd(j,3)
     icell(j,6)=1+icd(j,1)+2*icg(j,2)+4*icd(j,3)
     icell(j,7)=1+icg(j,1)+2*icd(j,2)+4*icd(j,3)
     icell(j,8)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
  end do
#endif


  do ind=1,twotondim
     do j=1,np
        if (ok(j,ind))then
           ! Compute parent cell adress for cells in ilevel or ilevel+1
           indp(j,ind)=ncoarse+(icell(j,ind)-1)*ngridmax+igrid(j,ind)
           ! Check if particles have leaked into level ilevel+1
           ! if so, set ok to false, but read values from split cell
           ok(j,ind)=(son(indp(j,ind))==0)
        end if
     end do
  end do


  do ind=1,twotondim
     do j=1,np
        if (ok(j,ind))then
           xx(j,1:ndim,ind)=(xg(igrid(j,ind),1:ndim)+xc(icell(j,ind),1:ndim)-skip_loc(1:ndim))*scale
        end if
     end do
  end do
  
end subroutine cic_get_cells



subroutine cic_get_vals(fluid_var,ind_grid,xpart,ind_grid_part,ng,np,ilevel,ilevel_only)
  use amr_commons
  use pm_commons
  use poisson_commons
  use hydro_commons, ONLY: nvar,uold
  implicit none
  integer::ng,np,ilevel
  logical::ilevel_only

  integer ,dimension(1:nvector)::ind_grid,ind_grid_part
#ifdef SOLVERmhd
  real(dp) ,dimension(1:nvector,1:nvar+3)::fluid_var
#else
  real(dp) ,dimension(1:nvector,1:nvar)::fluid_var
#endif
  real(dp) ,dimension(1:nvector,1:ndim)::xpart

  !------------------------------------------------------------------
  ! This routine returns the CIC cells and volumes for np particles.
  !------------------------------------------------------------------
  integer::i,j,ind,ivar

  ! Particle-based arrays

  integer ,dimension(1:nvector,1:twotondim),save::igrid,icell,kg
  real(dp),dimension(1:nvector,1:ndim),save::x0
  real(dp),dimension(1:nvector),save::vol_tot
  real(dp),dimension(1:nvector,1:ndim,1:twotondim),save::xx
  real(dp),dimension(1:nvector,1:twotondim),save::vol
  integer ,dimension(1:nvector,1:twotondim),save::indp
  logical ,dimension(1:nvector,twotondim),save::ok

  call cic_get_cells(indp,xx,vol,ok,ind_grid,xpart,ind_grid_part,ng,np,ilevel)

  fluid_var(1:np,1:nvar)=0._dp

  if (ilevel_only)then
     do ind=1,twotondim
        do j=1,np
           if(.not. ok(j,ind))vol(j,ind)=0.
        end do
     end do
     
     vol_tot(1:np)=0._dp
     do ind=1,twotondim
        do j=1,np
           vol_tot(j)=vol_tot(j)+vol(j,ind)
        end do
     end do
  end if

  do ivar=1,nvar
     do ind=1,twotondim
        do j=1,np
           
           fluid_var(j,ivar)=fluid_var(j,ivar)&
                +uold(indp(j,ind),ivar)*vol(j,ind)        
           
        end do
     end do
  end do

  if (ilevel_only)then
     do ivar=1,nvar
        do j=1,np
           if (vol_tot(j)>0.)fluid_var(j,ivar)=fluid_var(j,ivar)/vol_tot(j)
        end do
     end do
  end if
  

end subroutine cic_get_vals



