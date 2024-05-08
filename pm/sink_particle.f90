#if NDIM==3
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine create_sink
  use amr_commons
  use pm_commons
  use hydro_commons
  use clfind_commons
  use mpi_mod
  implicit none

  !----------------------------------------------------------------------------
  ! sink creation routine
  ! -remove all cloud particles, keep only global sink arrays
  ! -call clumpfinder for potential relevant peaks
  ! -flag peaks which are eligible for sink formation (flag 2)
  ! -new sink particles are created
  ! -new cloud particles are created
  ! -cloud particles are scattered to grid
  ! -accretion routine is called (with on_creation=.true.)
  !----------------------------------------------------------------------------

  integer::ilevel,ivar

  if(verbose)write(*,*)' Entering create_sink'

  ! Merge all particles to level 1
  do ilevel=levelmin-1,1,-1
     call merge_tree_fine(ilevel)
  end do

  ! Remove all particle clouds around old sinks (including the central one)
  call kill_entire_cloud(1)
#ifndef WITHOUTMPI
  call synchronize_sink_info
#endif

  ! DO NOT MODIFY FLAG2 BETWEEN CLUMP_FINDER AND MAKE_SINK_FROM_CLUMP
  if (create_sinks)then

     ! Run the clump finder,(produce no output, keep clump arrays allocated)
     call clump_finder(.false.,.true.)

     ! Trim clumps down to R_accretion ball around peaks
     if(clump_core)call trim_clumps

     ! Compute simple additive quantities and means (1st moments)
     call compute_clump_properties(uold(1,1))

     ! Compute stellar density field
     if(star.and.mass_star_AGN>0)call compute_rho_star

     ! Compute quantities relative to mean (2nd moments)
     call compute_clump_properties_round2

     ! Apply all checks and flag cells for sink formation
     call flag_formation_sites

     ! Create new sink particles if relevant
     do ilevel=levelmin,nlevelmax
        call make_sink_from_clump(ilevel)
     end do

     ! Deallocate clump finder arrays
     deallocate(npeaks_per_cpu)
     deallocate(ipeak_start)
     if (ntest>0)then
        deallocate(icellp)
        deallocate(levp)
        deallocate(testp_sort)
        deallocate(imaxp)
     endif
     call deallocate_all

  end if

  ! Remove merged sinks
  call clean_merged_sinks

  ! Create new cloud particles
  call create_cloud_from_sink

  ! Scatter particle to levelmax
  do ilevel=1,nlevelmax
     call make_tree_fine(ilevel)
     call kill_tree_fine(ilevel)
     call virtual_tree_fine(ilevel)
     ! Compute Bondi parameters
     if(hydro)call collect_acczone_avg(ilevel)
  end do

  ! Gather particles to levelmin
  do ilevel=nlevelmax,levelmin,-1
     if(hydro)call set_unew_sink(ilevel)
  end do
  do ilevel=nlevelmax,levelmin,-1
     ! Perform first accretion with seed mass
     if(hydro)call grow_sink(ilevel,.true.)
     call merge_tree_fine(ilevel)
  end do
  do ilevel=nlevelmax,levelmin,-1
     if(hydro)call set_uold_sink(ilevel)
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

  ! Compute and print accretion rates
  if(hydro)call compute_accretion_rate(.true.)

end subroutine create_sink
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine create_cloud_from_sink
  use amr_commons
  use pm_commons
  use hydro_commons
  use mpi_mod
  use constants, only: M_sun
  implicit none

  !----------------------------------------------------------------------------
  ! This routine creates the whole cloud of particles for each sink,
  ! Particles are produced in the right MPI domain and inserted in the
  ! linked list at level 1.
  ! The cloud radius is dble(ir_cloud)*dx_min, where dx_min is
  ! the cell size at levelmax_sink. For cosmo runs, the cloud radius is
  ! dx_min/aexp (therefore it is constant in *physical* units).
  !----------------------------------------------------------------------------

  real(dp)::scale,dx_min,rr,rmax,rmass
  integer ::icpu,isink,indp,ii,jj,kk,nx_loc,idim
  real(dp),dimension(1:ndim)::xrel
  real(dp),dimension(1:nvector,1:ndim)::xtest
  integer ,dimension(1:nvector)::ind_grid,cc,ind_cloud
  logical ,dimension(1:nvector)::ok_true
  logical ,dimension(1:ndim)::period
  logical ::in_box
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
#ifndef WITHOUTMPI
  integer ::info
  integer ::nsink_min,nsink_max
#endif

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ok_true=.true.

  if(numbtot(1,1)==0) return
  if(verbose)write(*,*)' Entering create_cloud_from_sink'


  if(nsink==0) return
#ifndef WITHOUTMPI
  ! Checking if all CPUs have the same number of sinks
  call MPI_ALLREDUCE(nsink,nsink_min,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(nsink,nsink_max,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
  if(nsink_min/=nsink_max)then
     if(myid==1)write(*,*)'Sink number mismatch!',nsink_min,nsink_max
     call clean_stop
  endif
#endif

  ! Level 1 linked list
  do icpu=1,ncpu
     if(numbl(icpu,1)>0)then
        ind_grid(1)=headl(icpu,1)
     endif
  end do

  period(1)=(nx==1)
  period(2)=(ny==1)
  period(3)=(nz==1)

  ! Mesh spacing in that level
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax_sink/aexp

  rmax=dble(ir_cloud)*dx_min
  rmass=dble(ir_cloud_massive)*dx_min

  do kk=-2*ir_cloud,2*ir_cloud
     xrel(3)=dble(kk)*dx_min/2
     do jj=-2*ir_cloud,2*ir_cloud
        xrel(2)=dble(jj)*dx_min/2
        do ii=-2*ir_cloud,2*ir_cloud
           xrel(1)=dble(ii)*dx_min/2
           rr=sqrt(sum(xrel**2))
           if(rr<=rmax)then
              do isink=1,nsink
                 xtest(1,1:ndim)=xsink(isink,1:ndim)+xrel(1:ndim)
                 in_box=.true.
                 do idim=1,ndim
                    if (period(idim) .and. xtest(1,idim)>boxlen)xtest(1,idim)=xtest(1,idim)-boxlen
                    if (period(idim) .and. xtest(1,idim)<0.)    xtest(1,idim)=xtest(1,idim)+boxlen
                    if (xtest(1,idim)<0.0 .or. xtest(1,idim)>boxlen)in_box=.false.
                 end do
                 cc(1)=0
                 if(in_box)call cmp_cpumap(xtest,cc,1)
                 if(cc(1).eq.myid)then
                    call remove_free(ind_cloud,1)
                    call add_list(ind_cloud,ind_grid,ok_true,1)
                    indp               = ind_cloud(1)
                    idp(indp)          = -isink
                    typep(indp)%family = FAM_CLOUD
                    typep(indp)%tag    = 0
                    levelp(indp)       = levelmin
                    if (rr<=rmass)then
                       ! check if direct_force is turned on
                       if(mass_sink_direct_force .ge. 0.0)then
                          if(msink(isink)<mass_sink_direct_force*M_sun/(scale_d*scale_l**ndim))then
                             mp(indp)=msink(isink)/dble(ncloud_sink_massive)
                          else
                             mp(indp)=0
                          endif
                       else
                          mp(indp)=msink(isink)/dble(ncloud_sink_massive)
                       endif
                    else
                       mp(indp)        = 0
                    end if
                    xp(indp,1:ndim)       = xtest(1,1:ndim)
                    vp(indp,1:ndim)       = vsink(isink,1:ndim)
                    tp(indp)              = tsink(isink)     ! Birth epoch
                 end if
              end do
           end if
        end do
     end do
  end do

  sink_jump(1:nsink,1:ndim,levelmin:nlevelmax)=0d0
  if(mass_sink_direct_force .ge. 0.0)then
     do isink=1,nsink
        direct_force_sink(isink)=(msink(isink) .ge. mass_sink_direct_force*M_sun/(scale_d*scale_l**ndim))
     end do
  else
     do isink=1,nsink
        direct_force_sink(isink)=.false.
     end do
  endif

end subroutine create_cloud_from_sink
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine kill_entire_cloud(ilevel)
  use pm_commons
  use amr_commons
  implicit none
  !----------------------------------------------------------------------------
  ! This routine removes cloud particles (including the central one)
  !----------------------------------------------------------------------------
  integer::ilevel
  integer::igrid,jgrid,ipart,jpart,next_part
  integer::ig,ip,npart1,npart2,icpu,ncache,istart
  integer,dimension(1:nvector)::ind_grid,ind_part,ind_grid_part
  logical,dimension(1:nvector)::ok=.true.

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel
  ! Gather sink and cloud particles.
  ! Loop over cpus
  do icpu=1,ncpu+nboundary
     if(icpu<=ncpu)then
        ncache=numbl(icpu,ilevel)
        istart=headl(icpu,ilevel)
     else
        ncache=numbb(icpu-ncpu,ilevel)
        istart=headb(icpu-ncpu,ilevel)
     end if
     igrid=istart
     ig=0
     ip=0
     ! Loop over grids
     do jgrid=1,ncache
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0
        ! Count sink and cloud particles
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if( is_cloud(typep(ipart)) ) then
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
              if( is_cloud(typep(ipart)) ) then
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
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine collect_acczone_avg(ilevel)
  use pm_commons
  use amr_commons
  use poisson_commons
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::info
#endif
  !----------------------------------------------------------------------------
  ! This routine is used to collect all relevant information to compute the
  ! the accretion rates. The information is collected level-by-level when
  ! going down in the call tree (leafs at the bottom), while accretion is
  ! performed on the way up.
  ! - first, the volume of each particle is computed. The volume is reduced if
  ! sinks are overlapping.
  ! - then a loop over all particles of ilevel (vectorized) is used to compute
  ! the volume-weighted quantities.
  !----------------------------------------------------------------------------

  integer::ilevel
  integer::igrid,jgrid,ipart,jpart,next_part
  integer::ig,ip,npart1,npart2,icpu,isink
  integer,dimension(1:nvector)::ind_grid,ind_part,ind_grid_part

  if(ilevel<levelmin)return
  if(verbose)write(*,111)ilevel

  ! Compute (volume weighted) averages over accretion zone
  wden=0d0; wvol=0d0; weth=0d0; wmom=0d0

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
              if( is_cloud(typep(ipart)) ) then
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
              if( is_cloud(typep(ipart)) ) then
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
#else
     wden_new=wden
     wvol_new=wvol
     weth_new=weth
     wmom_new=wmom
#endif
  endif

  do isink=1,nsink
     weighted_density(isink,ilevel)=wden_new(isink)
     weighted_volume(isink,ilevel)=wvol_new(isink)
     weighted_momentum(isink,ilevel,1:ndim)=wmom_new(isink,1:ndim)
     weighted_ethermal(isink,ilevel)=weth_new(isink)
  end do

111 format('   Entering collect_acczone_avg for level ',I2)

end subroutine collect_acczone_avg
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine collect_acczone_avg_np(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
  !----------------------------------------------------------------------------
  ! inner loop of collect_acczone_avg
  ! weighted gas quantities for each particle are computed
  ! no CIC averaging over quantities anymore as average over whole sink
  ! accretion zone is computed
  !----------------------------------------------------------------------------

  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid,ind_part,ind_grid_part
  integer::j,nx_loc,isink,idim,ind
#if NENER>0
  integer::irad
#endif
  real(dp)::d,e,v2
  real(dp)::scale,weight,dx_cloud,vol_cloud
  real(dp),dimension(1:ndim)::vv
#ifdef SOLVERmhd
  real(dp)::bx1,bx2,by1,by2,bz1,bz2
#endif
  real(dp),dimension(1:nvector,1:ndim),save::xpart
  integer ,dimension(1:nvector,1:twotondim),save::indp
  real(dp),dimension(1:nvector,1:ndim,1:twotondim)::xx
  real(dp),dimension(1:nvector,1:twotondim)::vol
  logical, dimension(1:nvector,1:twotondim)::ok

  ! Compute volume of each cloud particle
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_cloud=(0.5D0**nlevelmax_sink)*scale/aexp/2 ! factor of 2 hard-coded
  vol_cloud=dx_cloud**ndim

  ! Copy cloud particle coordinates
  do idim=1,ndim
     do j=1,np
        xpart(j,idim)=xp(ind_part(j),idim)
     end do
  end do

  ! Compute cloud particle CIC weights at the current level
  call cic_get_cells(indp,xx,vol,ok,ind_grid,xpart,ind_grid_part,ng,np,ilevel)

  do ind=1,twotondim
     do j=1,np
        if(ok(j,ind))then

           ! Convert uold to primitive variables
           d=max(uold(indp(j,ind),1),smallr)
           vv(1)=uold(indp(j,ind),2)/d
           vv(2)=uold(indp(j,ind),3)/d
           vv(3)=uold(indp(j,ind),4)/d

           ! Compute fluid total energy
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
#if NENER>0
           do irad=0,nener-1
              e=e-uold(indp(j,ind),inener+irad)
           end do
#endif
           e=e/d ! Specific energy
           v2=sum(vv**2)
           e=e-0.5d0*v2 ! Remove kinetic energy

           ! Get sink index
           isink=-idp(ind_part(j))

           ! Get cloud particle CIC weight
           weight=vol_cloud*vol(j,ind)

           ! Compute sink average quantities
           wvol(isink)=wvol(isink)+weight
           wden(isink)=wden(isink)+weight*d
           wmom(isink,1:ndim)=wmom(isink,1:ndim)+weight*d*vv(1:ndim)
           weth(isink)=weth(isink)+weight*d*e

        endif
     end do
  end do

end subroutine collect_acczone_avg_np
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine grow_sink(ilevel,on_creation)
  use pm_commons
  use amr_commons
  use hydro_commons
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::info
#endif
  !----------------------------------------------------------------------------
  ! This routine performs accretion onto the sink. It vectorizes the loop
  ! over all sink cloud particles and calls accrete_sink as soon as nvector
  ! particles are collected
  !----------------------------------------------------------------------------

  integer::ilevel
  logical::on_creation
  integer::igrid,jgrid,ipart,jpart,next_part
  integer::ig,ip,npart1,npart2,icpu,isink,lev
  integer,dimension(1:nvector)::ind_grid,ind_part,ind_grid_part

  if(accretion_scheme=='none'.and.(.not.on_creation))return
  if(verbose)write(*,111)ilevel

  ! Compute sink accretion rates
  call compute_accretion_rate(.false.)

  ! Reset new sink variables
  msink_new=0d0; msmbh_new=0d0; dmfsink_new=0d0
  xsink_new=0d0; vsink_new=0d0; lsink_new=0d0; delta_mass_new=0d0

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
              if (is_cloud(typep(ipart)) ) then
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
              if( is_cloud(typep(ipart)) ) then
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
     call MPI_ALLREDUCE(msmbh_new,msmbh_all,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(dmfsink_new,dmfsink_all,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(xsink_new,xsink_all,nsinkmax*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(vsink_new,vsink_all,nsinkmax*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(lsink_new,lsink_all,nsinkmax*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(delta_mass_new,delta_mass_all,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#else
     msink_all=msink_new
     msmbh_all=msmbh_new
     dmfsink_all=dmfsink_new
     xsink_all=xsink_new
     vsink_all=vsink_new
     lsink_all=lsink_new
     delta_mass_all=delta_mass_new
#endif
  endif

  do isink=1,nsink
     if (msink(isink)>0.)then ! if msink=0 then there was no accretion and thus no shift

        ! Update mass from accretion
        msink(isink)=msink(isink)+msink_all(isink)
        msmbh(isink)=msmbh(isink)+msmbh_all(isink)
        dmfsink(isink)=dmfsink(isink)+dmfsink_all(isink)

        ! Reset jump in old sink coordinates
        do lev=levelmin,nlevelmax
           sink_jump(isink,1:ndim,lev)=sink_jump(isink,1:ndim,lev)-xsink(isink,1:ndim)
        end do

        ! Accrete to sink variables in the frame of the initial sink position
        xsink(isink,1:ndim)=xsink(isink,1:ndim)+xsink_all(isink,1:ndim)/(msink(isink)+msmbh(isink))
        vsink(isink,1:ndim)=vsink(isink,1:ndim)+vsink_all(isink,1:ndim)/(msink(isink)+msmbh(isink))
        lsink(isink,1:ndim)=lsink(isink,1:ndim)+lsink_all(isink,1:ndim)-cross(xsink_all(isink,1:ndim),vsink_all(isink,1:ndim))/(msink(isink)+msmbh(isink))

        ! Store jump in new sink coordinates
        do lev=levelmin,nlevelmax
           sink_jump(isink,1:ndim,lev)=sink_jump(isink,1:ndim,lev)+xsink(isink,1:ndim)
        end do
     end if

     ! Store accreted mass
     if (agn.and.ok_blast_agn(isink))then
        delta_mass(isink)=delta_mass_all(isink) ! only the most recent accretion
     else
        delta_mass(isink)=delta_mass(isink)+delta_mass_all(isink)
     end if

  end do

111 format('   Entering grow_sink for level ',I2)

end subroutine grow_sink
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine accrete_sink(ind_grid,ind_part,ind_grid_part,ng,np,ilevel,on_creation)
  use amr_commons
  use pm_commons
  use hydro_commons
  use constants, only: pi, twopi, c_cgs, factG_in_cgs, M_sun, mH, sigma_t, ev2erg
#ifdef RT
  use rt_hydro_commons,only: rtunew
  use rt_parameters,only: nGroups, iGroups, group_egy, rt_AGN, group_egy_AGNfrac
#endif
  implicit none
  !----------------------------------------------------------------------------
  ! This routine is called by subroutine grow_sink. It performs accretion
  ! for nvector particles. Routine is not very efficient.
  ! Optimize if taking too long...
  !----------------------------------------------------------------------------

  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  logical::on_creation
  integer::j,nx_loc,isink,ivar,idim,ind
  real(dp)::d,e,density,volume
#ifdef SOLVERmhd
  real(dp)::bx1,bx2,by1,by2,bz1,bz2
#endif
#if NENER>0
  integer::irad
#endif
#ifdef RT
  integer::igroup
  real(dp)::scale_Np,scale_Fp
  real(dp)::scale_vol,scale_evtocode,dert,Np_inj
#endif
  real(dp)::factG,scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m
  real(dp)::dx,dx_loc,dx_min,dx_cloud,scale,vol_min,vol_loc,vol_cloud,weight,m_acc,m_acc_smbh
  ! Grid based arrays
  real(dp),dimension(1:nvector,1:ndim)::xpart
  real(dp),dimension(1:nvector,1:ndim,1:twotondim)::xx
  real(dp),dimension(1:nvector,1:twotondim)::vol
  ! Particle based arrays
  logical,dimension(1:nvector,1:twotondim)::ok
  integer ,dimension(1:nvector,1:twotondim)::indp
  real(dp),dimension(1:ndim)::vv

  real(dp),dimension(1:ndim)::r_rel,v_rel,x_acc,p_acc,l_acc
  real(dp)::fbk_ener_AGN,fbk_mom_AGN,r_len
  logical,dimension(1:ndim)::period

  real(dp)::tan_theta,cone_dist,orth_dist
  real(dp),dimension(1:ndim)::cone_dir
  real(dp)::acc_ratio,v_AGN

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_m=scale_d*scale_l**ndim

  period(1)=(nx==1)
  period(2)=(ny==1)
  period(3)=(nz==1)

  ! Gravitational constant
  factG=1d0
  if(cosmo)factG=3d0/4d0/twopi*omega_m*aexp

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim
  dx_min=scale*0.5D0**nlevelmax_sink/aexp
  vol_min=dx_min**ndim

  ! Compute volume of each cloud particle
  dx_cloud=dx_min/2 ! factor of 2 hard-coded
  vol_cloud=dx_cloud**ndim

#ifdef RT
  call rt_units(scale_Np, scale_Fp)
  scale_vol=scale_l**ndim
  scale_evtocode=eV2erg/(scale_d*scale_l**5/scale_t**2)
#endif

  ! Jet geometry safety net
  cone_opening = max(tiny(0.0d0),cone_opening)
  cone_opening = min(cone_opening, 180d0)
  tan_theta = tan(pi/180d0*cone_opening/2) ! tangent of half of the opening angle

  ! Get cloud particle CIC weights
  do idim=1,ndim
     do j=1,np
        xpart(j,idim)=xp(ind_part(j),idim)
     end do
  end do
  call cic_get_cells(indp,xx,vol,ok,ind_grid,xpart,ind_grid_part,ng,np,ilevel)

  ! Loop over eight CIC volumes
  do ind=1,twotondim
     do j=1,np
        if(ok(j,ind))then ! Only volumes at the current level

           ! Convert uold to primitive variables
           d=max(uold(indp(j,ind),1),smallr)
           vv(1)=uold(indp(j,ind),2)/d
           vv(2)=uold(indp(j,ind),3)/d
           vv(3)=uold(indp(j,ind),4)/d

           ! Compute the gas total specific energy (kinetic plus thermal)
           ! removing all the other energies, if any.
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
#if NENER>0
           do irad=0,nener-1
              e=e-uold(indp(j,ind),inener+irad)
           end do
#endif
           e=e/d ! Specific energy

           ! Get sink index
           isink=-idp(ind_part(j))

           ! Reference frame relative to the sink position
           r_rel(1:ndim)=xx(j,1:ndim,ind)-xsink(isink,1:ndim)
           do idim=1,ndim
              if (period(idim) .and. r_rel(idim)>boxlen*0.5d0)r_rel(idim)=r_rel(idim)-boxlen
              if (period(idim) .and. r_rel(idim)<boxlen*(-0.5d0))r_rel(idim)=r_rel(idim)+boxlen
           end do
           r_len=sqrt(sum(r_rel**2))

           ! Reference frame relative to the sink velocity
           v_rel(1:ndim)=vv(1:ndim)-vsink(isink,1:ndim)

           ! Cloud particle CIC weight
           weight=vol_cloud*vol(j,ind)

           ! Get sink average density
           density=rho_gas(isink)
           volume=volume_gas(isink)
           if (volume<=0. .or. density<=0.)then
              print*,'something might be going wrong here...',weight,volume,density,ilevel
           endif

           ! Compute accreted mass
           if (on_creation)then
              if (new_born(isink))then
                 ! on sink creation, new sinks
                 m_acc     =mass_sink_seed*M_sun/scale_m*weight/volume
                 m_acc_smbh=mass_smbh_seed*M_sun/scale_m*weight/volume
                 if (agn_acc_method=='mass') then
                    m_acc      = m_acc * (d/density)
                    m_acc_smbh = m_acc_smbh * (d/density)
                 endif
              else
                 ! on sink creation, preexisting sinks
                 m_acc     =0
                 m_acc_smbh=0
              end if
           else
              if (bondi_accretion)then
                 m_acc      = dMsink_overdt(isink)*dtnew(ilevel)*weight/volume
                 m_acc_smbh = dMsmbh_overdt(isink)*dtnew(ilevel)*weight/volume
              end if

              if (threshold_accretion.and.d_sink>0.0)then
                 m_acc=c_acc*weight*(d-d_sink)
              end if

              if (agn_acc_method=='mass') then
                 m_acc      = m_acc * (d/density)
                 m_acc_smbh = m_acc_smbh * (d/density)
              endif

              if(agn.and.msink(isink).gt.0)then
                 acc_ratio=dMsmbh_overdt(isink)/(4d0*pi*factG_in_cgs*msmbh(isink)*mH/(0.1d0*sigma_T*c_cgs)*scale_t)
                 if (AGN_fbk_mode_switch_threshold > 0.0) then
                    if (acc_ratio > AGN_fbk_mode_switch_threshold) then
                       ! Eddington ratio higher than AGN_fbk_mode_switch_threshold -> energy
                       AGN_fbk_frac_ener = 1
                       AGN_fbk_frac_mom = 0
                    else
                       ! Eddington ratio lower than AGN_fbk_mode_switch_threshold -> momentum
                       AGN_fbk_frac_ener = 0
                       AGN_fbk_frac_mom = 1
                    end if
                 end if
                 v_AGN = (2*0.1d0*epsilon_kin/kin_mass_loading)**0.5d0*c_cgs ! in cm/s
                 if (agn_inj_method=='mass') then
                    fbk_ener_AGN=AGN_fbk_frac_ener*min(delta_mass(isink)*T2_AGN/scale_T2*weight/volume*d/density,T2_max/scale_T2*weight*d)
                    fbk_mom_AGN=AGN_fbk_frac_mom*kin_mass_loading*delta_mass(isink)*v_AGN/scale_v*weight/volume*d/density/(1d0-cos(pi/180*cone_opening/2))
                 else if (agn_inj_method=='volume') then
                    fbk_ener_AGN=AGN_fbk_frac_ener*min(delta_mass(isink)*T2_AGN/scale_T2*weight/volume, T2_max/scale_T2*weight*d)
                    fbk_mom_AGN=AGN_fbk_frac_mom*kin_mass_loading*delta_mass(isink)*v_AGN/scale_v*weight/volume/(1d0-cos(pi/180*cone_opening/2))
                 endif
              end if
           end if

           m_acc     =max(m_acc,0.0_dp)
           m_acc_smbh=max(m_acc_smbh,0.0_dp)

           ! Accreted relative center of mass
           x_acc(1:ndim)=(m_acc+m_acc_smbh)*r_rel(1:ndim)

           ! Accreted relative momentum
           p_acc(1:ndim)=(m_acc+m_acc_smbh)*v_rel(1:ndim)

           ! Accreted relative angular momentum
           l_acc(1:ndim)=(m_acc+m_acc_smbh)*cross(r_rel(1:ndim),v_rel(1:ndim))

           ! Add accreted properties to sink variables
           msink_new(isink)=msink_new(isink)+m_acc
           msmbh_new(isink)=msmbh_new(isink)+m_acc_smbh
           dmfsink_new(isink)=dmfsink_new(isink)+m_acc
           xsink_new(isink,1:ndim)=xsink_new(isink,1:ndim)+x_acc(1:ndim)
           vsink_new(isink,1:ndim)=vsink_new(isink,1:ndim)+p_acc(1:ndim)
           lsink_new(isink,1:ndim)=lsink_new(isink,1:ndim)+l_acc(1:ndim)
           if(mass_smbh_seed>0.0)then
              delta_mass_new(isink)=delta_mass_new(isink)+m_acc_smbh
           else
              delta_mass_new(isink)=delta_mass_new(isink)+m_acc
           end if

           m_acc=m_acc+m_acc_smbh
           ! Accrete mass, momentum and gas total energy
           unew(indp(j,ind),1)=unew(indp(j,ind),1)-m_acc/vol_loc
           unew(indp(j,ind),2:ndim+1)=unew(indp(j,ind),2:ndim+1)-m_acc*vv(1:ndim)/vol_loc
           unew(indp(j,ind),ndim+2)=unew(indp(j,ind),ndim+2)-m_acc*e/vol_loc
           ! Note that we do not accrete magnetic fields and non-thermal energies.

           ! Accrete passive scalars
           do ivar=imetal,nvar
              unew(indp(j,ind),ivar)=unew(indp(j,ind),ivar)-m_acc*uold(indp(j,ind),ivar)/d/vol_loc
           end do

           ! AGN feedback
           ! Only for sinks that could accrete.
           ! Can deposit thermal dump or/and momentum kicks (and/or radiation if RT used).
           if( .not. on_creation)then
              if(agn)then
                 if(ok_blast_agn(isink).and.delta_mass(isink)>0.0)then
                    if(AGN_fbk_frac_ener.gt.0.0)then ! thermal AGN feedback
                       unew(indp(j,ind),ndim+2)=unew(indp(j,ind),ndim+2)+fbk_ener_AGN/vol_loc
                    end if

                    if(AGN_fbk_frac_mom.gt.0.0)then ! momentum AGN feedback
                       ! checking if particle is in cone
                       cone_dir(1:ndim)=lsink(isink,1:ndim)/sqrt(sum(lsink(isink,1:ndim)**2))
                       cone_dist=sum(r_rel(1:ndim)*cone_dir(1:ndim))
                       orth_dist=sqrt(sum((r_rel(1:ndim)-cone_dist*cone_dir(1:ndim))**2))
                       if (orth_dist.le.abs(cone_dist)*tan_theta)then
                          unew(indp(j,ind),2:ndim+1)=unew(indp(j,ind),2:ndim+1)+fbk_mom_AGN*r_rel(1:ndim)/(r_len)/vol_loc
                          unew(indp(j,ind),ndim+2)=unew(indp(j,ind),ndim+2)+sum(fbk_mom_AGN*r_rel(1:ndim)/(r_len)*vv(1:ndim))/vol_loc
                       end if
                    end if
                 end if
#ifdef RT
                 if(rt_AGN)then
                    dert=0.1d0*delta_mass(isink)*(c_cgs/scale_v)**2
                    do igroup=1,nGroups
                       Np_inj=dert * group_egy_AGNfrac(igroup) / (scale_evtocode*group_egy(igroup)) / (vol_loc*scale_vol) / scale_Np*weight/volume
                       rtunew(indp(j,ind),iGroups(igroup))=rtunew(indp(j,ind),iGroups(igroup))+Np_inj
                    enddo
                 end if
#endif
              end if
           end if
        endif
     end do
  end do
end subroutine accrete_sink
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine compute_accretion_rate(write_sinks)
  use pm_commons
  use amr_commons
  use hydro_commons
  use constants, only: pi, twopi, c_cgs, factG_in_cgs, M_sun, mH, sigma_T
  use mpi_mod
  implicit none
  logical::write_sinks

  !----------------------------------------------------------------------------
  ! This routine computes the accretion rate onto the sink particles based
  ! on the information collected in collect accretion.
  ! It also creates output for the sink particle positions
  !----------------------------------------------------------------------------

  integer::i,nx_loc,isink
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m
  real(dp)::factG,d_star,boost
  real(dp)::vrel2,c2,density,volume,ethermal,dx_min,scale,mgas,v_bondi
  real(dp)::r2_smbh,rho_inf_smbh
  real(dp),dimension(1:ndim)::velocity
  real(dp),dimension(1:nsinkmax)::dMEDoverdt,r2,rho_inf
  real(dp),dimension(1:nsinkmax)::dMEDoverdt_smbh
  real(dp)::T2_gas,delta_mass_min

  ! Gravitational constant
  factG=1d0
  if(cosmo)factG=3d0/4d0/twopi*omega_m*aexp

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_m=scale_d*scale_l**ndim
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax_sink/aexp
  d_star=n_star/scale_nH

  ! Compute sink particle accretion rate by averaging contributions from all levels
  do isink=1,nsink

     dMsink_overdt(isink)=0
     dMsmbh_overdt(isink)=0
     dMBHoverdt(isink)=0
     dMBHoverdt_smbh(isink)=0
     dMEDoverdt(isink)=0
     dMEDoverdt_smbh(isink)=0

     ! Compute sink sphere average quantities
     density=0d0; volume=0d0; velocity=0d0; ethermal=0d0
     do i=levelmin,nlevelmax
        density=density+weighted_density(isink,i)
        ethermal=ethermal+weighted_ethermal(isink,i)
        velocity(1:ndim)=velocity(1:ndim)+weighted_momentum(isink,i,1:ndim)
        volume=volume+weighted_volume(isink,i)
     end do
     mgas=density
     density=density/(volume+tiny(0.0_dp))
     if (volume<=0. .or. density<=0.)then
        if(myid==1)print*,'something might be going wrong here...',volume,density,isink,idsink(isink),xsink(isink,1:ndim)
        cycle
     endif
     ! Compute Bondi-Hoyle accretion rate in code units
     if (star.and.acc_sink_boost.lt.0.0)then
        boost=max((density/(boost_threshold_density/scale_nH))**2,1.0_dp)
     else
        boost=abs(acc_sink_boost)
     end if

     velocity(1:ndim)=velocity(1:ndim)/(density*volume+tiny(0.0_dp))
     ethermal=ethermal/(density*volume+tiny(0.0_dp))
     c2=MAX((gamma-1.0d0)*ethermal,smallc**2)*boost**(-2d0/3d0)
     c2sink(isink)=c2
     vrel2=SUM((velocity(1:ndim)-vsink(isink,1:ndim))**2)
     if(bondi_use_vrel)then
        v_bondi=sqrt(c2+vrel2)
     else
        v_bondi=sqrt(c2)
     endif

     ! Bondi radius
     r2(isink)=(factG*msink(isink)/v_bondi**2)**2

     ! Extrapolate to rho_inf
     rho_inf(isink)=density/(bondi_alpha(ir_cloud*0.5d0*dx_min/(r2(isink)+tiny(0.0_dp))**0.5d0))

     ! Compute Bondi-Hoyle accretion rate in code units
     dMBHoverdt(isink)=4*pi*rho_inf(isink)*r2(isink)*v_bondi

     ! Compute Eddington accretion rate in code units
     dMEDoverdt(isink)=4*pi*factG_in_cgs*msink(isink)*mH/(0.1d0*sigma_T*c_cgs)*scale_t

     ! Compute final sink accretion rate
     if(bondi_accretion)dMsink_overdt(isink)=dMBHoverdt(isink)
     if(eddington_limit)dMsink_overdt(isink)=min(dMBHoverdt(isink),dMEDoverdt(isink))

     if(smbh.and.mass_smbh_seed>0.0)then
        r2_smbh=(factG*msmbh(isink)/v_bondi**2)**2
        rho_inf_smbh=density/(bondi_alpha(ir_cloud*0.5d0*dx_min/(r2_smbh+tiny(0.0_dp))**0.5d0))
        dMBHoverdt_smbh(isink)=4*pi*rho_inf_smbh*r2_smbh*v_bondi
        dMEDoverdt_smbh(isink)=4*pi*factG_in_cgs*msmbh(isink)*mH/(0.1d0*sigma_T*c_cgs)*scale_t
        if(bondi_accretion)dMsmbh_overdt(isink)=dMBHoverdt_smbh(isink)
        if(eddington_limit)dMsmbh_overdt(isink)=min(dMBHoverdt(isink),dMEDoverdt_smbh(isink))
        dMsink_overdt(isink)=max(0d0,dMBHoverdt(isink)-dMsmbh_overdt(isink))
     end if

     ! Store average quantities for diagnostics
     eps_sink(isink)=ethermal
     rho_gas(isink)=density
     volume_gas(isink)=volume
     vel_gas(isink,1:ndim)=velocity(1:ndim)

     if (agn.and.dMsink_overdt(isink)>0.0)then
        ! Check whether we should have AGN feedback
        if(T2_min<=0.0)then ! If zero or less, we always deposit feedback
           ok_blast_agn(isink)=.true.
        else ! Checking temperature of surrounding gas
           ok_blast_agn(isink)=.false.
           T2_gas=ethermal*scale_T2 ! in Kelvin
           delta_mass_min = mgas*(T2_min-T2_gas)/(T2_AGN-T2_min)
           if((T2_gas.ge.T2_min).or.(delta_mass(isink).ge.mgas*(T2_min-T2_gas)/(T2_AGN-T2_min)))then
              ok_blast_agn(isink)=.true.
           end if
        end if
     end if

     if(msink(isink).ge.max_mass_nsc*M_sun/scale_m.and.mass_smbh_seed>0.0)dMsink_overdt(isink)=0.0

  end do

  if (write_sinks)then
     call print_sink_properties(dMEDoverdt,dMEDoverdt_smbh,rho_inf,r2)
  end if

contains
  ! Routine to return alpha, defined as rho/rho_inf, for a critical
  ! Bondi accretion solution. The argument is x = r / r_Bondi.
  ! This is from Krumholz et al. (AJC)
  REAL(dp) function bondi_alpha(x)
    implicit none
    REAL(dp) x
    REAL(dp), PARAMETER :: XMIN=0.01d0, xMAX=2.0d0
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
    lambda_c    = exp(1.5d0) / 4
    !     Deal with the off-the-table cases
    if (x .le. XMIN) then
       bondi_alpha = lambda_c / sqrt(2d0 * x**ndim)
    else if (x .ge. XMAX) then
       bondi_alpha = exp(1d0/x)
    else
       !     We are on the table
       idx = floor ((NTABLE-1) * log(x/XMIN) / log(XMAX/XMIN))
       xtable = exp(log(XMIN) + idx*log(XMAX/XMIN)/(NTABLE-1))
       xtablep1 = exp(log(XMIN) + (idx+1)*log(XMAX/XMIN)/(NTABLE-1d0))
       alpha_exp = log(x/xtable) / log(xtablep1/xtable)
       !     Note the extra +1s below because of fortran 1 offset arrays
       bondi_alpha = alphatable(idx+1) * (alphatable(idx+2)/alphatable(idx+1))**alpha_exp
    end if
  end function bondi_alpha

end subroutine compute_accretion_rate
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine print_sink_properties(dMEDoverdt,dMEDoverdt_smbh,rho_inf,r2)
  use pm_commons
  use amr_commons
  use hydro_commons
  use constants, only: twopi, M_sun, pc2cm, yr2sec
  use mpi_mod
  implicit none
  real(dp),dimension(1:nsinkmax)::dMEDoverdt,rho_inf,r2
  real(dp),dimension(1:nsinkmax)::dMEDoverdt_smbh
  integer::i,isink,nx_loc
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m
  real(dp)::l_abs,l_max,factG,scale,dx_min
  real(dp),dimension(1:ndim)::skip_loc

  ! Gravitational constant
  factG=1d0
  if(cosmo)factG=3d0/4d0/twopi*omega_m*aexp

  ! Mesh spacing in that level
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_min=0.5D0**nlevelmax_sink*scale/aexp

  ! Scaling factors
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_m=scale_d*scale_l**ndim

  if (smbh) then
     if(myid==1.and.nsink>0)then
        xmsink(1:nsink)=msink(1:nsink)
        call quick_sort_dp(xmsink(1),idsink_sort(1),nsink)
        write(*,*)'Number of sink = ',nsink
        write(*,'(" ============================================================================================")')
        if(mass_smbh_seed>0)then
           write(*,'(" Id  Mass(Msol)     Bondi(Msol/yr) Edd(Msol/yr)  BH: Mass(Msol)  Bondi(Msol/yr)  Edd(Msol/yr)    x              y              z      DeltaM(Msol)")')
        else
           write(*,'(" Id  Mass(Msol)     Bondi(Msol/yr) Edd(Msol/yr)   x              y              z              DeltaM(Msol)")')
        endif
        write(*,'(" ============================================================================================")')
        do i=nsink,max(nsink-30,1),-1
           isink=idsink_sort(i)
           if(mass_smbh_seed>0)then
              write(*,'(I3,12(1X,1PE14.7))')idsink(isink) &
                & ,msink(isink)*scale_m/M_sun &
                & ,dMBHoverdt(isink)*scale_m/scale_t/(M_sun/yr2sec) &
                & ,dMEDoverdt(isink)*scale_m/scale_t/(M_sun/yr2sec) &
                & ,msmbh(isink)*scale_m/M_sun &
                & ,dMBHoverdt_smbh(isink)*scale_m/scale_t/(M_sun/yr2sec) &
                & ,dMEDoverdt_smbh(isink)*scale_m/scale_t/(M_sun/yr2sec) &
                & ,xsink(isink,1:ndim),delta_mass(isink)*scale_m/M_sun
           else
              write(*,'(I3,12(1X,1PE14.7))')idsink(isink) &
                & ,msink(isink)*scale_m/M_sun &
                & ,dMBHoverdt(isink)*scale_m/scale_t/(M_sun/yr2sec) &
                & ,dMEDoverdt(isink)*scale_m/scale_t/(M_sun/yr2sec) &
                & ,xsink(isink,1:ndim),delta_mass(isink)*scale_m/M_sun
           endif
        end do
        write(*,'(" ============================================================================================")')
        if(verbose_AGN)then
           write(*,'(" Id  rho(H/cc)      rho_inf(H/cc)  Mgas(Msol)     cs(km/s)       rBondi(pc)     vgasx(km/s)    vgasy(km/s)    vgasz(km/s)    vsinkx(km/s)   vsinky(km/s)   vsinkz(km/s) ")')
           write(*,'(" ============================================================================================")')
           do i=nsink,max(nsink-30,1),-1
              isink=idsink_sort(i)
              write(*,'(I3,12(1X,1PE14.7))')idsink(isink),rho_gas(isink)*scale_nH,rho_inf(isink)*scale_nH &
                   & ,rho_gas(isink)*volume_gas(isink)*scale_m/M_sun,sqrt(c2sink(isink))*scale_v/1e5 &
                   & ,sqrt(r2(isink))*scale_l/pc2cm,vel_gas(isink,1:ndim)*scale_v/1e5,vsink(isink,1:ndim)*scale_v/1e5
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
        write(*,*)'Total mass in sink [Msol] = ',sum(msink(1:nsink))*scale_m/M_sun
        write(*,*)'simulation time [yr] = ',t*scale_t/yr2sec
        write(*,'(" =============================================================================================================================================")')
        write(*,'("   Id     M[Msol]          x             y             z         vx[km/s]      vy[km/s]      vz[km/s]     spin/spmax    Mdot[Msol/y]   age[yr]")')
        write(*,'(" =============================================================================================================================================")')
        do i=nsink,1,-1
           isink=idsink_sort(i)
           l_abs=(lsink(isink,1)**2+lsink(isink,2)**2+lsink(isink,3)**2)**0.5d0
           l_max=msink(isink)*sqrt(factG*msink(isink)/(dble(ir_cloud)*dx_min))*(dble(ir_cloud)*dx_min)
           write(*,'(I5,10(2X,1PE12.5))')&
                & idsink(isink),&
                & msink(isink)*scale_m/M_sun,&
                & xsink(isink,1:ndim),&
                & vsink(isink,1:ndim)*scale_v/1d5,&
                & l_abs/l_max,&
                & dMsink_overdt(isink)*scale_m/scale_t/(M_sun/yr2sec),&
                & (t-tsink(isink))*scale_t/yr2sec
        end do
        write(*,'(" =============================================================================================================================================")')
     endif
  endif
end subroutine print_sink_properties
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine make_sink_from_clump(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  use clfind_commons
  use constants, only:
  use mpi_mod
  implicit none

  !----------------------------------------------------------------------------
  ! This routine uses creates a sink in every cell which was flagged (flag2)
  ! The global sink variables are updated
  ! The true RAMSES particle is NOT produced here...
  !----------------------------------------------------------------------------

  integer ::ilevel
  integer ::ncache,nnew,ivar,ngrid,index_sink,index_sink_tot
  integer ::igrid,ix,iy,iz,ind,i,iskip,isink,nx_loc
  integer ::ntot,ntot_all
  integer ,dimension(1:nvector)::ind_grid,ind_cell
  integer ,dimension(1:nvector)::ind_grid_new,ind_cell_new
  integer ,dimension(1:ncpu)::ntot_sink_cpu
#ifndef WITHOUTMPI
  integer ::icpu, info
  integer ,dimension(1:ncpu)::ntot_sink_all
#endif
  logical ::ok_free
  real(dp)::d,u,v,w,e,delta_d,v2
  real(dp)::birth_epoch
  real(dp)::dx,dx_loc,scale,vol_loc
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp),dimension(1:nvar)::z
  real(dp),dimension(1:ndim)::skip_loc,x
  real(dp),dimension(1:twotondim,1:ndim)::xc
#ifdef SOLVERmhd
  real(dp)::bx1,bx2,by1,by2,bz1,bz2
#endif
#if NENER>0
  integer ::irad
#endif

#if NDIM==3
  if(verbose)write(*,*)'entering make_sink_from_clump for level ',ilevel

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

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
  msink_new=0d0; msmbh_new=0d0; dmfsink_new=0d0
  xsink_new=0d0; vsink_new=0d0; lsink_new=0d0; delta_mass_new=0d0
  tsink_new=0d0; oksink_new=0d0; idsink_new=0; new_born_new=.false.

  ! Count number of new sinks (flagged cells)
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

  ! Compute global sink statistics
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(ntot,ntot_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#else
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
        write(*,'(" Level = ",I6," Number of new sinks produced= ",I6," Total sinks =",I8)')&
             & ilevel,ntot_all,nsink
     endif
  end if

  ! Check wether max number of sink is reached
  ok_free=(nsink+ntot_all<=nsinkmax)
  if(.not. ok_free)then
     if(myid==1)write(*,*)'global list of sink particles is too long'
     if(myid==1)write(*,*)'New sink particles',ntot_all
     if(myid==1)write(*,*)'Increase nsinkmax'
     call clean_stop
  end if

  ! Create new sink particles
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
              d=max(uold(ind_cell_new(i),1),smallr)
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
              do irad=0,nener-1
                 e=e-uold(ind_cell_new(i),inener+irad)
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
              call true_max(x(1),x(2),x(3),nlevelmax_sink)

              ! Give a tiny bit of mass to the sink...
              delta_d=d*1d-10
              msink_new(index_sink)=delta_d*vol_loc
              msmbh_new(index_sink)=delta_d*vol_loc
              delta_mass_new(index_sink)=msmbh_new(index_sink)
              dmfsink_new(index_sink)=delta_d*vol_loc

              ! Global index of the new sink
              oksink_new(index_sink)=1d0
              idsink_new(index_sink)=index_sink_tot

              ! Store properties of the new sink
              tsink_new(index_sink)=birth_epoch
              xsink_new(index_sink,1:ndim)=x(1:ndim)
              vsink_new(index_sink,1)=u
              vsink_new(index_sink,2)=v
              vsink_new(index_sink,3)=w
              lsink_new(index_sink,1:ndim)=0
              new_born_new(index_sink)=.true.

              ! Convert back to conservative variable
              d=d-delta_d
              e=e*d
#ifdef SOLVERmhd
              e=e+0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)
#endif
              e=e+0.5d0*d*(u**2+v**2+w**2)
#if NENER>0
              do irad=0,nener-1
                 e=e+uold(ind_cell_new(i),inener+irad)
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
  call MPI_ALLREDUCE(msink_new ,msink_all ,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(msmbh_new ,msmbh_all ,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(dmfsink_new ,dmfsink_all ,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(xsink_new ,xsink_all ,nsinkmax*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(vsink_new ,vsink_all ,nsinkmax*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(lsink_new ,lsink_all ,nsinkmax*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(delta_mass_new,delta_mass_all,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(oksink_new,oksink_all,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(idsink_new,idsink_all,nsinkmax,MPI_INTEGER         ,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(tsink_new ,tsink_all ,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(new_born_new,new_born_all,nsinkmax,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,info)
#else
  msink_all=msink_new
  msmbh_all=msmbh_new
  dmfsink_all=dmfsink_new
  xsink_all=xsink_new
  vsink_all=vsink_new
  lsink_all=lsink_new
  delta_mass_all=delta_mass_new
  oksink_all=oksink_new
  idsink_all=idsink_new
  tsink_all=tsink_new
  new_born_all=new_born_new
#endif
  do isink=1,nsink
     if(oksink_all(isink)==1)then
        msink(isink)=msink_all(isink)
        msmbh(isink)=msmbh_all(isink)
        dmfsink(isink)=dmfsink_all(isink)
        xsink(isink,1:ndim)=xsink_all(isink,1:ndim)
        vsink(isink,1:ndim)=vsink_all(isink,1:ndim)
        lsink(isink,1:ndim)=lsink_all(isink,1:ndim)
        delta_mass(isink)=delta_mass_all(isink)
        idsink(isink)=idsink_all(isink)
        tsink(isink)=tsink_all(isink)
        new_born(isink)=new_born_all(isink)
        fsink(isink,1:ndim)=0.0
        fsink_partial(isink,1:ndim,levelmin:nlevelmax)=0.0
        vsold(isink,1:ndim,ilevel)=vsink_all(isink,1:ndim)
        vsnew(isink,1:ndim,ilevel)=vsink_all(isink,1:ndim)
     endif
  end do
#endif

end subroutine make_sink_from_clump
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine true_max(x,y,z,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use clfind_commons, only:ivar_clump
  use poisson_commons, only:rho
  implicit none
  real(dp)::x,y,z
  integer::ilevel

  !----------------------------------------------------------------------------
  ! Description: This subroutine takes the cell of maximum density and computes
  ! the true maximum by expanding the density around the cell center to second order.
  !----------------------------------------------------------------------------

  integer::k,j,i,nx_loc,counter, ioft, n
  integer,dimension(1:threetondim)::cell_index,cell_lev
  real(dp)::det,dx,dx_loc,scale,disp_max,numerator
  real(dp),dimension(-1:1,-1:1,-1:1)::cube3
  real(dp),dimension(1:threetondim,1:ndim)::xtest
  real(dp),dimension(1:nvector,1:ndim)::xtest_copy
  real(dp),dimension(1:ndim)::gradient,displacement
  real(dp),dimension(1:ndim,1:ndim)::hess,minor
  real(dp)::smallreal=1d-100

#if NDIM==3

  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  counter=0
  do i=-1,1
     do j=-1,1
        do k=-1,1
           counter=counter+1
           xtest(counter,1)=x+i*dx_loc
           xtest(counter,2)=y+j*dx_loc
           xtest(counter,3)=z+k*dx_loc
        end do
     end do
  end do

  do ioft=0, threetondim-1, nvector
     n=min(threetondim-ioft, nvector)
     do i=1,n
       xtest_copy(i,1:ndim) = xtest(ioft+i,1:ndim)
     enddo
     call get_cell_index(cell_index(ioft+1:ioft+n), cell_lev(ioft+1:ioft+n), &
          &  xtest_copy, ilevel, n)
  end do

  counter=0
  if(ivar_clump==0)then
     do i=-1,1
        do j=-1,1
           do k=-1,1
              counter=counter+1
              cube3(i,j,k)=rho(cell_index(counter))
           end do
        end do
     end do
  else if(hydro)then
     do i=-1,1
        do j=-1,1
           do k=-1,1
              counter=counter+1
              cube3(i,j,k)=uold(cell_index(counter),1)
           end do
        end do
     end do
  else
     return
  end if

  ! Compute gradient
  gradient(1)=0.5d0*(cube3(1,0,0)-cube3(-1,0,0))/dx_loc
  gradient(2)=0.5d0*(cube3(0,1,0)-cube3(0,-1,0))/dx_loc
  gradient(3)=0.5d0*(cube3(0,0,1)-cube3(0,0,-1))/dx_loc

  if (maxval(abs(gradient(1:ndim)))==0.)return

  ! Compute hessian
  hess(1,1)=(cube3(1,0,0)+cube3(-1,0,0)-2*cube3(0,0,0))/dx_loc**2
  hess(2,2)=(cube3(0,1,0)+cube3(0,-1,0)-2*cube3(0,0,0))/dx_loc**2
  hess(3,3)=(cube3(0,0,1)+cube3(0,0,-1)-2*cube3(0,0,0))/dx_loc**2

  hess(1,2)=0.25d0*(cube3(1,1,0)+cube3(-1,-1,0)-cube3(1,-1,0)-cube3(-1,1,0))/dx_loc**2
  hess(2,1)=hess(1,2)
  hess(1,3)=0.25d0*(cube3(1,0,1)+cube3(-1,0,-1)-cube3(1,0,-1)-cube3(-1,0,1))/dx_loc**2
  hess(3,1)=hess(1,3)
  hess(2,3)=0.25d0*(cube3(0,1,1)+cube3(0,-1,-1)-cube3(0,1,-1)-cube3(0,-1,1))/dx_loc**2
  hess(3,2)=hess(2,3)

  ! Determinant
  det=    hess(1,1)*hess(2,2)*hess(3,3)+hess(1,2)*hess(2,3)*hess(3,1)+hess(1,3)*hess(2,1)*hess(3,2) &
       & -hess(1,1)*hess(2,3)*hess(3,2)-hess(1,2)*hess(2,1)*hess(3,3)-hess(1,3)*hess(2,2)*hess(3,1)

  ! Matrix of minors
  minor(1,1)=hess(2,2)*hess(3,3)-hess(2,3)*hess(3,2)
  minor(2,2)=hess(1,1)*hess(3,3)-hess(1,3)*hess(3,1)
  minor(3,3)=hess(1,1)*hess(2,2)-hess(1,2)*hess(2,1)

  minor(1,2)=-1d0*(hess(2,1)*hess(3,3)-hess(2,3)*hess(3,1))
  minor(2,1)=minor(1,2)
  minor(1,3)=hess(2,1)*hess(3,2)-hess(2,2)*hess(3,1)
  minor(3,1)=minor(1,3)
  minor(2,3)=-1d0*(hess(1,1)*hess(3,2)-hess(1,2)*hess(3,1))
  minor(3,2)=minor(2,3)

  ! Displacement of the true max from the cell center
  displacement=0
  do i=1,ndim
     do j=1,ndim
        numerator = gradient(j)*minor(i,j)
        if(numerator>0) displacement(i)=displacement(i)-numerator/(det+smallreal*numerator)
     end do
  end do

  ! Clipping the displacement in order to keep max in the cell
  disp_max=maxval(abs(displacement(1:ndim)))
  if (disp_max > dx_loc*0.499999)then
     displacement(1)=displacement(1)/disp_max*dx_loc*0.499999d0
     displacement(2)=displacement(2)/disp_max*dx_loc*0.499999d0
     displacement(3)=displacement(3)/disp_max*dx_loc*0.499999d0
  end if

  x=x+displacement(1)
  y=y+displacement(2)
  z=z+displacement(3)

#endif
end subroutine true_max
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine update_sink(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use sink_feedback_parameters
  use constants, only: twopi, M_sun, yr2sec
  use mpi_mod
  implicit none
  integer::ilevel

  !----------------------------------------------------------------------------
  ! This routine is called at the leafs of the tree structure (right after
  ! update time). Here is where the global sink variables vsink and xsink are
  ! updated by summing the conributions from all levels.
  !----------------------------------------------------------------------------

  integer::lev,isink,jsink,nx_loc,idim,istellar
  logical::iyoung,jyoung,overlap,merge_flag
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::dteff,scale,dx_min
  real(dp)::t_larson1,rr,rmax,rmax2,factG,v1_v2,mcom,fsink_norm
  real(dp),dimension(1:ndim)::xcom,vcom,lcom,r_rel
  logical,dimension(1:ndim)::period
  real(dp),dimension(1:nsink,1:ndim)::xsinkold, fsinkold

#if NDIM==3

  if(verbose)write(*,*)'Entering update_sink for level ',ilevel

  period(1)=(nx==1)
  period(2)=(ny==1)
  period(3)=(nz==1)

  ! Mesh spacing in that level
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax_sink/aexp
  rmax=dble(ir_cloud)*dx_min ! Linking length in physical units
  rmax2=rmax*rmax

  ! Lifetime of first larson core in code units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  t_larson1=merging_timescale*yr2sec/scale_t

  ! Gravitational constant
  factG=1d0
  if(cosmo)factG=3d0/4d0/twopi*omega_m*aexp

  ! Set overlap mass to sink mass
  msum_overlap=msink

  ! Check for overlapping sinks
  do isink=1,nsink-1
     if (msink(isink)>0.)then
        do jsink=isink+1,nsink

           ! Compute relative distance
           r_rel(1:ndim)=xsink(jsink,1:ndim)-xsink(isink,1:ndim)
           do idim=1,ndim
              if (period(idim) .and. r_rel(idim)>boxlen*0.5d0)   r_rel(idim)=r_rel(idim)-boxlen
              if (period(idim) .and. r_rel(idim)<boxlen*(-0.5d0))r_rel(idim)=r_rel(idim)+boxlen
           end do
           rr=sum(r_rel**2)

           ! Check for overlap
           overlap=rr<4*rmax2 .and. msink(jsink)>0.

           if(overlap)then
              msum_overlap(isink)=msum_overlap(isink)+msink(jsink)
              msum_overlap(jsink)=msum_overlap(jsink)+msink(isink)

              ! Merging based on relative distance
              merge_flag=rr<4*dx_min**2 ! Sinks are within two cells from each other

              ! Merging based on relative velocity
              if(mass_merger_vel_check>0 .and. (msink(isink)+msink(jsink)).ge.mass_merger_vel_check*M_sun/(scale_d*scale_l**ndim)) then
                 v1_v2=(vsink(isink,1)-vsink(jsink,1))**2+(vsink(isink,2)-vsink(jsink,2))**2+(vsink(isink,3)-vsink(jsink,3))**2
                 merge_flag=merge_flag .and. 2*factG*(msink(isink)+msink(jsink))/sqrt(rr)>v1_v2
              end if

              ! Merging based on sink age
              if (merging_timescale>0d0)then
                 iyoung=(t-tsink(isink)<t_larson1)
                 jyoung=(t-tsink(jsink)<t_larson1)
                 merge_flag=merge_flag .and. (iyoung .or. jyoung)
              end if

              if (merge_flag.eqv..true.)then

                 if(myid==1)then
                    write(*,*)'> Merging sink ',idsink(jsink),' into sink ',idsink(isink)
                    if(verbose_AGN)then
                       write(*,*)'>> Sink #1: ',idsink(isink)
                       write(*,*)msink(isink)/M_sun*(scale_d*scale_l**ndim)
                       write(*,*)xsink(isink,1:ndim)
                       write(*,*)'>> Sink #2: ',idsink(jsink)
                       write(*,*)msink(jsink)/M_sun*(scale_d*scale_l**ndim)
                       write(*,*)xsink(jsink,1:ndim)
                    endif
                 endif

                 ! Set new values of remaining sink (keep one with larger index)
                 ! Compute centre of mass quantities
                 mcom     =(msink(isink)+msink(jsink))
                 xcom(1:ndim)=xsink(isink,1:ndim)+msink(jsink)*r_rel(1:ndim)/mcom
                 vcom(1:ndim)=(msink(isink)*vsink(isink,1:ndim)+msink(jsink)*vsink(jsink,1:ndim))/mcom
                 lcom(1:ndim)=msink(isink)*cross((xsink(isink,1:ndim)-xcom(1:ndim)),vsink(isink,1:ndim)-vcom(1:ndim))+ &
                      &    msink(jsink)*cross((xsink(jsink,1:ndim)-xcom(1:ndim)),vsink(jsink,1:ndim)-vcom(1:ndim))

                 ! Reset jump in old sink coordinates
                 do lev=levelmin,nlevelmax
                    sink_jump(isink,1:ndim,lev)=sink_jump(isink,1:ndim,lev)-xsink(isink,1:ndim)
                 end do

                 ! Compute merged quantities
                 msink(isink)        = mcom
                 msmbh(isink)        = msmbh(isink)+msmbh(jsink)
                 dmfsink(isink)      = dmfsink(isink)+dmfsink(jsink)
                 delta_mass(isink)   = delta_mass(isink)+delta_mass(jsink)
                 xsink(isink,1:ndim) = xcom(1:ndim)
                 vsink(isink,1:ndim) = vcom(1:3)
                 lsink(isink,1:ndim) = lcom(1:ndim)+lsink(isink,1:ndim)+lsink(jsink,1:ndim)
                 tsink(isink)        = min(tsink(isink),tsink(jsink))
                 idsink(isink)       = min(idsink(isink),idsink(jsink))

                 ! Store jump in new sink coordinates
                 do lev=levelmin,nlevelmax
                    sink_jump(isink,1:ndim,lev)=sink_jump(isink,1:ndim,lev)+xsink(isink,1:ndim)
                 end do

                 ! Zero mass of the sink that was merged in
                 msink(jsink)=0
                 msmbh(jsink)=0
                 dmfsink(jsink)=0
                 msum_overlap(jsink)=0
                 delta_mass(jsink)=0

                 ! check whether there are stellar particles attached to the merged in sink
                 if(stellar)then
                     do istellar = 1, nstellar
                         if(id_stellar(istellar).eq.idsink(jsink))then
                             id_stellar(istellar) = idsink(isink)
                         endif
                     end do
                 endif

              end if
           end if
        end do
     end if
  end do

  ! Store old xsink and fsink for the gradient descent timestep
  xsinkold=0.0
  fsinkold=0.0
  do isink=1,nsink
     if(msink(isink)>0.0)then
        fsinkold(isink,1:ndim)=fsink(isink,1:ndim)
        xsinkold(isink,1:ndim)=xsink(isink,1:ndim)
     endif
  enddo

  ! Updating sink positions

  fsink=0
  call f_sink_sink

  vsold(1:nsink,1:ndim,ilevel)=vsnew(1:nsink,1:ndim,ilevel)
  vsnew(1:nsink,1:ndim,ilevel)=vsink(1:nsink,1:ndim)

  ! Loop over sinks
  do isink=1,nsink
     if(msink(isink)>0.0)then
        ! Sum force contributions from all levels and gather
        do lev=levelmin,nlevelmax
           fsink(isink,1:ndim)=fsink(isink,1:ndim)+fsink_partial(isink,1:ndim,lev)
        end do
        if (.not. direct_force_sink(isink))then
           fsink(isink,1:ndim)=fsink(isink,1:ndim)/dble(ncloud_sink)
        end if

        ! Compute timestep for the synchronization
        if (new_born(isink))then
           ! No sync necessary for newly produced sink
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

        ! This is the kick-kick (half old half new timestep)
        vsink(isink,1:ndim)=0.5D0*(dtnew(ilevel)+dteff)*fsink(isink,1:ndim)+vsink(isink,1:ndim)

        ! Save the velocity
        vsnew(isink,1:ndim,ilevel)=vsink(isink,1:ndim)

        ! and this is the drift (only for the global sink variable)
        xsink(isink,1:ndim)=xsink(isink,1:ndim)+vsink(isink,1:ndim)*dtnew(ilevel)

        ! Gradient descent - Barzilai&Borwein-like
        if (sink_descent) then
           xsink_graddescent(1:nsink,1:ndim)=0.0
           fsink_norm=NORM2(fsink(isink,1:ndim))
           gamma_grad_descent = 0.0d0
           graddescent_over_dt = 0.0d0
           if (.not. new_born(isink))then
              do idim=1,ndim
                 gamma_grad_descent = gamma_grad_descent + (xsink(isink,idim)-xsinkold(isink,idim))*(fsink(isink,idim)-fsinkold(isink,idim))
              enddo
              if(gamma_grad_descent>0.0)then
                 gamma_grad_descent = fudge_graddescent*dtnew(ilevel)*SQRT(ABS(gamma_grad_descent)/(NORM2(fsink(isink,1:ndim)-fsinkold(isink,1:ndim)))**2)
                 ! Require thatthe sink cannot move more than half a grid
                 if(gamma_grad_descent*fsink_norm>dx_min/2.0) then
                    xsink_graddescent(isink,1:ndim) = fsink(isink,1:ndim) * dx_min/2.0/fsink_norm
                 else
                    xsink_graddescent(isink,1:ndim) = fsink(isink,1:ndim) * gamma_grad_descent
                 endif
                 ! Uopdate the sink position
                 xsink(isink,1:ndim)=xsink(isink,1:ndim)+ xsink_graddescent(isink,1:ndim)
                 ! Store the descent velocity for the time-stepping
                 graddescent_over_dt(isink) = NORM2(xsink_graddescent(isink,1:ndim))/dtnew(ilevel)
              endif
           endif
        endif

        new_born(isink)=.false.
     end if
  end do
  ! End loop over sinks

  ! Store deepest level
  sinkint_level=ilevel

#endif
end subroutine update_sink
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine update_cloud(ilevel)
  use amr_commons
  use pm_commons
  implicit none
  integer::ilevel

  !----------------------------------------------------------------------------
  ! Update sink cloud particle properties
  ! - the particles are moved whenever the level of the grid they sit in is updated
  ! - the amount of drift they get is according to their levelp
  ! - since this is happening on the way down, at level ilevel all particles with
  ! level >= ilevel will be moved. Therefore, the sink_jump for all levels >= ilevel
  ! is set to zero on exit.
  !----------------------------------------------------------------------------

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
     write(*,*)'Sink drift due to accretion relative to grid size at level ',ilevel
     do isink=1,nsink
        write(*,'("#sink: ",I6," drift: ",3(1PE12.5))')isink,sink_jump(isink,1:ndim,ilevel)/dx_loc
     end do
  end if

  sink_jump(1:nsink,1:ndim,ilevel:nlevelmax_sink)=0d0

111 format('   Entering update_cloud for level ',I2)

end subroutine update_cloud
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine upd_cloud(ind_part,np)
  use amr_commons
  use pm_commons
  use poisson_commons
  implicit none
  integer::np
  integer,dimension(1:nvector)::ind_part

  !----------------------------------------------------------------------------
  ! Vector loop called by update_cloud
  !----------------------------------------------------------------------------

  integer::j,idim,isink,lev
  real(dp),dimension(1:nvector,1:ndim)::new_xp,new_vp
  integer,dimension(1:nvector)::level_p

  ! Overwrite cloud particle mass with sink mass
  do j=1,np
     if ( is_cloud(typep(ind_part(j))) .and. mp(ind_part(j))>0.0d0 ) then
        isink = -idp(ind_part(j))
        mp(ind_part(j)) = msink(isink)/dble(ncloud_sink_massive)
     end if
  end do

  ! Store velocity
  do idim=1,ndim
     do j=1,np
        new_vp(j,idim) = vp(ind_part(j),idim)
     end do
  end do

  ! Overwrite cloud particle velocity with sink velocity
  ! is going to be overwritten again before move
  do idim=1,ndim
     do j=1,np
        if ( is_cloud(typep(ind_part(j))) ) then
           isink = -idp(ind_part(j))
           new_vp(j,idim) = vsink(isink,idim)
        end if
     end do
  end do

  ! Write back velocity
  do idim=1,ndim
     do j=1,np
        vp(ind_part(j),idim) = new_vp(j,idim)
     end do
  end do

  ! Read level
  do j=1,np
     level_p(j) = levelp(ind_part(j))
  end do

  ! Update position
  do idim=1,ndim
     do j=1,np
        new_xp(j,idim) = xp(ind_part(j),idim)
     end do
  end do
  do idim=1,ndim
     do j=1,np
        if ( is_cloud(typep(ind_part(j))) ) then
           isink = -idp(ind_part(j))
           lev = level_p(j)
           new_xp(j,idim) = new_xp(j,idim) + sink_jump(isink,idim,lev)
        endif
     end do
  end do

 ! Write back postion
  do idim=1,ndim
     do j=1,np
        xp(ind_part(j),idim) = new_xp(j,idim)
     end do
  end do

end subroutine upd_cloud
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine clean_merged_sinks
  use pm_commons
  use amr_commons
  implicit none

  !----------------------------------------------------------------------------
  ! This routine cleans up the list of merged sink particles.
  ! Merging is done in update_sink (mass and acc_rate strictly zero),
  ! accrete_sinks sets cloud masses to 0,
  ! here we just clean the list of sinks.
  !----------------------------------------------------------------------------
  integer::i,j,mergers

  if(nsink==0)return

  mergers=0
  do j=1,nsink
     if (msink(j)<=0d0)then ! if sink has been merged to another one
         mergers=mergers+1
         msink(j)=-10
     end if
  end do

  if (myid==1 .and. mergers>0)write(*,*)'Clean ',mergers,' merged sinks'

  ! Sort sink particle arrays to account for merged sinks that disappeared
  i=1
  do while (mergers>0)
     if (msink(i)==-10.)then ! if sink has been merged to another one
        if (myid==1 .and. mergers>0)write(*,*)'Cleaning ',idsink(i)

        mergers=mergers-1
        nsink=nsink-1

        ! Let them all slide back one index
        do j=i,nsink
           msink(j)=msink(j+1)
           msmbh(j)=msmbh(j+1)
           dmfsink(j)=dmfsink(j+1)
           xsink(j,1:ndim)=xsink(j+1,1:ndim)
           vsink(j,1:ndim)=vsink(j+1,1:ndim)
           lsink(j,1:ndim)=lsink(j+1,1:ndim)
           new_born(j)=new_born(j+1)
           tsink(j)=tsink(j+1)
           idsink(j)=idsink(j+1)
           msum_overlap(j)=msum_overlap(j+1)
           delta_mass(j)=delta_mass(j+1)
        end do

        ! Whipe last position in the sink list
        msink(nsink+1)=0d0
        msmbh(nsink+1)=0d0
        dmfsink(nsink+1)=0d0
        xsink(nsink+1,1:ndim)=0d0
        vsink(nsink+1,1:ndim)=0d0
        lsink(nsink+1,1:ndim)=0d0
        new_born(nsink+1)=.false.
        tsink(nsink+1)=0d0
        idsink(nsink+1)=0
        msum_overlap(nsink+1)=0d0
        delta_mass(nsink+1)=0d0
     else
        i=i+1
     end if
  end do

end subroutine clean_merged_sinks
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine f_gas_sink(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  use constants, only: pi, twopi
  use mpi_mod
  implicit none

  !----------------------------------------------------------------------------
  ! In this subroutine the sink-gas force contributions are calculated.
  ! A plummer-sphere with radius ssoft is used for softening
  !----------------------------------------------------------------------------
#ifndef WITHOUTMPI
  integer::info
#endif
  integer ::ilevel
  integer ::igrid,ngrid,ncache,i,ind,iskip,ix,iy,iz,isink
  integer ::nx_loc,idim
  real(dp)::dx,dx_loc,scale,vol_loc,dx_min,factG
  real(dp),dimension(1:twotondim,1:ndim)::xc
  real(dp),dimension(1:ndim)::skip_loc

  logical ,dimension(1:nvector)::ok
  integer ,dimension(1:nvector)::ind_grid,ind_cell
  real(dp),dimension(1:nvector,1:ndim)::xx,ff
  real(dp),dimension(1:nvector)::d2,mcell,denom
  real(dp)::rho_tff,rho_tff_tot,d_min
  logical ,dimension(1:ndim)::period

#if NDIM==3

  ! Gravitational constant
  factG=1d0
  if(cosmo)factG=3d0/4d0/twopi*omega_m*aexp

  !  Cell spacing at that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim
  dx_min=scale*0.5D0**nlevelmax_sink/aexp
  ssoft=sink_soft*dx_min

  ! Set position of cell centers relative to grid centre
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     xc(ind,1)=(dble(ix)-0.5D0)*dx
     xc(ind,2)=(dble(iy)-0.5D0)*dx
     xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  rho_tff=0

  fsink_all=0
  fsink_new=0

  period(1)=(nx==1)
  period(2)=(ny==1)
  period(3)=(nz==1)

  ! Loop over sinks
  do isink=1,nsink
     if (direct_force_sink(isink))then

        d_min=boxlen**2

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
              d2=0d0
              do idim=1,ndim
                 do i=1,ngrid
                    ! zero order Ewald sum
                    ff(i,idim)=xsink(isink,idim)-xx(i,idim)
                    if (period(idim))then
                       if(ff(i,idim)>0.5d0*boxlen)ff(i,idim)=ff(i,idim)-boxlen
                       if(ff(i,idim)<-0.5d0*boxlen)ff(i,idim)=ff(i,idim)+boxlen
                    end if
                    d2(i)=d2(i)+ff(i,idim)**2
                 end do
              end do

              ! Store minimum distance of cell in current level to isink
              do i=1,ngrid
                 d_min=min(d_min,d2(i))
              end do

              ! Compute sqrt(1/(ssoft**2+d2(i))) to save time
              do i=1,ngrid
                 denom(i)=(ssoft**2+d2(i))**(-1.5d0)
              end do

              ! Compute gas acceleration due to sink
              do i=1,ngrid
                 ff(i,1:ndim)=denom(i)*ff(i,1:ndim)
              end do

              ! Add gas acceleration due to sink
              do i=1,ngrid
                 f(ind_cell(i),1:ndim)=f(ind_cell(i),1:ndim)+factG*msink(isink)*ff(i,1:ndim)
              end do

              ! Add sink acceleration due to gas
              do i=1,ngrid
                 if(ok(i))then
                    fsink_new(isink,1:ndim)=fsink_new(isink,1:ndim)-factG*mcell(i)*ff(i,1:ndim)
                 end if
              end do
           end do !end loop over cells
        end do !end loop over grids

        d_min=d_min**0.5d0
        d_min=max(ssoft,d_min)
        rho_tff=max(rho_tff,max(msink(isink),msum_overlap(isink))/(4d0/3d0*pi*d_min**ndim))

     end if !end if direct force
  end do !end loop over sinks

  ! Collect sink acceleration from CPUs
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

  ! Collect rho due to sinks for current level - used for timestep computation
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(rho_tff,rho_tff_tot,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
#else
  rho_tff_tot=rho_tff
#endif
  rho_sink_tff(ilevel)=rho_tff_tot

  if (ilevel==nlevelmax_sink)call make_virtual_fine_dp(phi(1),ilevel)

#endif
end subroutine f_gas_sink
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine f_sink_sink
  use amr_commons
  use pm_commons
  use constants, only: twopi
  use mpi_mod
  implicit none

  !----------------------------------------------------------------------------
  ! In this subroutine the sink-sink force contribution is calculated by direct
  ! n^2 - summation. A plummer-sphere with radius 4 cells is used for softening
  !----------------------------------------------------------------------------
  integer::isink,idim,jsink
  real(dp),allocatable,dimension(:)::d2
  real(dp),allocatable,dimension(:,:)::ff
  logical,dimension(1:ndim)::period
  real(dp)::factG

#if NDIM==3

  allocate(d2(1:nsink))
  allocate(ff(1:nsink,1:ndim))

  ! Gravitational constant
  factG=1d0
  if(cosmo)factG=3d0/4d0/twopi*omega_m*aexp

  period(1)=(nx==1)
  period(2)=(ny==1)
  period(3)=(nz==1)

  do isink=1,nsink
     if(msink(isink)>0.0)then
        if (direct_force_sink(isink))then
           d2=0d0
           ff=0d0
           do idim=1,ndim
              ! Compute relative position and distances
              if (period(idim))then
                 do jsink=1,nsink
                    if (direct_force_sink(jsink))then
                       ff(jsink,idim)=xsink(jsink,idim)-xsink(isink,idim)
                       if(ff(jsink,idim)>0.5d0*boxlen)ff(jsink,idim)=ff(jsink,idim)-boxlen
                       if(ff(jsink,idim)<-0.5d0*boxlen)ff(jsink,idim)=ff(jsink,idim)+boxlen
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
           ! Compute acceleration
           do jsink=1,nsink
              if (direct_force_sink(jsink))then
                 ff(jsink,1:ndim)=factG*msink(jsink)/(ssoft**2+d2(jsink))**1.5d0*ff(jsink,1:ndim)
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
     end if
  end do

#endif
end subroutine f_sink_sink
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine read_sink_params()
  use pm_commons
  use amr_commons
  use constants, only: pi,yr2sec
  implicit none

  !----------------------------------------------------------------------------
  ! Read sink related parameters and perform some 'sanity chekcs'
  !----------------------------------------------------------------------------

  real(dp)::dx_min,scale,cty
  integer::nx_loc
  namelist/sink_params/n_sink,rho_sink,d_sink,accretion_scheme,merging_timescale,&
       ir_cloud_massive,sink_soft,mass_sink_direct_force,ir_cloud,nsinkmax,create_sinks,&
       check_energies,mass_sink_seed,mass_smbh_seed,c_acc,nlevelmax_sink,&
       eddington_limit,acc_sink_boost,mass_merger_vel_check,&
       clump_core,verbose_AGN,T2_AGN,T2_min,cone_opening,mass_halo_AGN,mass_clump_AGN,mass_star_AGN,&
       AGN_fbk_frac_ener,AGN_fbk_frac_mom,T2_max,boost_threshold_density,&
       epsilon_kin,AGN_fbk_mode_switch_threshold,kin_mass_loading,bondi_use_vrel,smbh,agn,max_mass_nsc,&
       agn_acc_method,agn_inj_method,sink_descent,gamma_grad_descent,fudge_graddescent
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v

  if(.not.cosmo) call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

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

  if (nlevelmax_sink .eq. 0) then
     nlevelmax_sink = nlevelmax
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

  if (mass_sink_seed <= 0)then
     if(myid==1)write(*,*)'Sink seed mass is not specified. Exiting.'
     call clean_stop
  end if


  ! Check for accretion scheme
  if (accretion_scheme=='bondi')bondi_accretion=.true.
  if (accretion_scheme=='threshold')threshold_accretion=.true.

  ! For sink formation and accretion a threshold must be given
  if (create_sinks .or. (accretion_scheme .ne. 'none'))then

     ! Check for threshold
     if (.not. cosmo)then

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
           dx_min=0.5d0**nlevelmax_sink*scale
           d_sink=T2_star/scale_T2 *pi/16/(dx_min**2)
           if(myid==1)write(*,*)'d_sink = ',d_sink
           if(myid==1)write(*,*)'rho_sink = ',d_sink*scale_d
           if(myid==1)write(*,*)'n_sink = ',d_sink*scale_nH
        end if
     end if

     endif
  end if

  if (merging_timescale > 0.)then
     cty=scale_t/yr2sec
     cont_speed=-1d0/(merging_timescale/cty)
  end if

  ! Check for periodic boundary conditions
  if (nx==1 .or. ny==1 .or. nz==1)then
     if (mass_sink_direct_force .ge. 0.)then
        if(myid==1)print*, 'periodic boundaries in combination with '
        if(myid==1)print*, 'direct force sinks are not treated accurately....'
     end if
  end if

end subroutine read_sink_params
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
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
  !----------------------------------------------------------------------------
  ! This routine returns the CIC cells and volumes for np particles
  !----------------------------------------------------------------------------
  integer::i,j,ind,idim,nx_loc,ix,iy,iz
  real(dp)::dx,dx_loc,scale,vol_loc
  ! Grid-based arrays
  integer ,dimension(1:nvector)::ind_cell
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle-based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x,dd,dg
  integer ,dimension(1:nvector,1:ndim),save::ig,id,igg,igd,icg,icd
  integer ,dimension(1:nvector,1:twotondim),save::igrid,icell,kg
  real(dp),dimension(1:nvector,1:ndim),save::x0
  real(dp),dimension(1:ndim)::skip_loc
  real(dp),dimension(1:twotondim,1:ndim)::xc

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim

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

  do i=1,ng
     ind_cell(i)=father(ind_grid(i))
  end do

  ! Gather neighboring father cells (should be present anytime!)
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
        if (x(j,idim)<0.5d0 .or. x(j,idim)>5.5)then
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
        id(j,idim)=int(dd(j,idim))
        dd(j,idim)=dd(j,idim)-id(j,idim)
        dg(j,idim)=1.0D0-dd(j,idim)
        ig(j,idim)=id(j,idim)-1
     end do
  end do

  ! Compute cloud volumes
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

  ! Compute parent grids
  do idim=1,ndim
     do j=1,np
        igg(j,idim)=ig(j,idim)/2
        igd(j,idim)=id(j,idim)/2
     end do
  end do
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
           ok(j,ind)=.false.
           indp(j,ind)=nbors_father_cells(ind_grid_part(j),kg(j,ind))
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
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
! not used
subroutine cic_get_vals(fluid_var,ind_grid,xpart,ind_grid_part,ng,np,ilevel,ilevel_only)
  use amr_commons
  use pm_commons
  use poisson_commons
  use hydro_commons, ONLY: nvar,uold
  implicit none

  !----------------------------------------------------------------------------
  ! This routine returns the CIC cells and volumes for np particles.
  !----------------------------------------------------------------------------
  integer::ng,np,ilevel
  integer::j,ind,ivar
  logical::ilevel_only

  integer ,dimension(1:nvector)::ind_grid,ind_grid_part
#ifdef SOLVERmhd
  real(dp) ,dimension(1:nvector,1:nvar+3)::fluid_var
#else
  real(dp) ,dimension(1:nvector,1:nvar)::fluid_var
#endif
  real(dp) ,dimension(1:nvector,1:ndim)::xpart

  ! Particle-based arrays
  real(dp),dimension(1:nvector),save::vol_tot
  real(dp),dimension(1:nvector,1:ndim,1:twotondim),save::xx
  real(dp),dimension(1:nvector,1:twotondim),save::vol
  integer ,dimension(1:nvector,1:twotondim),save::indp
  logical ,dimension(1:nvector,twotondim),save::ok

  call cic_get_cells(indp,xx,vol,ok,ind_grid,xpart,ind_grid_part,ng,np,ilevel)

  fluid_var(1:np,1:nvar)=0

  if (ilevel_only)then
     do ind=1,twotondim
        do j=1,np
           if(.not. ok(j,ind))vol(j,ind)=0
        end do
     end do

     vol_tot(1:np)=0
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
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine set_unew_sink(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !----------------------------------------------------------------------------
  ! This routine sets array unew to its initial value uold before creating
  ! new sinks. unew is set to zero in virtual boundaries.
  !----------------------------------------------------------------------------
  integer::i,ivar,ind,icpu,iskip

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Set unew to uold for myid cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
#ifdef SOLVERmhd
     do ivar=1,nvar+3
#else
     do ivar=1,nvar
#endif
        do i=1,active(ilevel)%ngrid
           unew(active(ilevel)%igrid(i)+iskip,ivar) = uold(active(ilevel)%igrid(i)+iskip,ivar)
        end do
     end do
  end do

  ! Set unew to 0 for virtual boundary cells
  do icpu=1,ncpu
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
#ifdef SOLVERmhd
     do ivar=1,nvar+3
#else
     do ivar=1,nvar
#endif
        do i=1,reception(icpu,ilevel)%ngrid
#ifdef LIGHT_MPI_COMM
           unew(reception(icpu,ilevel)%pcomm%igrid(i)+iskip,ivar)=0
#else
           unew(reception(icpu,ilevel)%igrid(i)+iskip,ivar)=0
#endif
        end do
     end do
  end do
  end do

111 format('   Entering set_unew_sink for level ',i2)

end subroutine set_unew_sink
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine set_uold_sink(ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none

  !----------------------------------------------------------------------------
  ! This routine sets array uold to its new value unew after the hydro step
  !----------------------------------------------------------------------------
  integer::ilevel
  integer::i,ivar,ind,iskip

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Reverse update boundaries
#ifdef SOLVERmhd
  do ivar=1,nvar+3
#else
  do ivar=1,nvar
#endif
     call make_virtual_reverse_dp(unew(1,ivar),ilevel)
#ifdef SOLVERmhd
  end do
#else
  end do
#endif

  ! Set uold to unew for myid cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
#ifdef SOLVERmhd
     do ivar=1,nvar+3
#else
     do ivar=1,nvar
#endif
        do i=1,active(ilevel)%ngrid
           uold(active(ilevel)%igrid(i)+iskip,ivar) = unew(active(ilevel)%igrid(i)+iskip,ivar)
        end do
     end do
  end do
111 format('   Entering set_uold_sink for level ',i2)

end subroutine set_uold_sink
#endif
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
#ifndef WITHOUTMPI
subroutine synchronize_sink_info
  use pm_commons
  use mpi_mod
  implicit none
  !----------------------------------------------------------------------------
  ! This routine syncronizes sink variables across all CPUs if MPI is used.
  ! This is done to prevent roundoff errors.
  !----------------------------------------------------------------------------
  integer::info

  call MPI_BCAST(msink,      nsinkmax, MPI_DOUBLE_PRECISION, 1, MPI_COMM_WORLD, info)
  call MPI_BCAST(msmbh,      nsinkmax, MPI_DOUBLE_PRECISION, 1, MPI_COMM_WORLD, info)
  call MPI_BCAST(dmfsink,    nsinkmax, MPI_DOUBLE_PRECISION, 1, MPI_COMM_WORLD, info)
  call MPI_BCAST(xsink,    3*nsinkmax, MPI_DOUBLE_PRECISION, 1, MPI_COMM_WORLD, info)
  call MPI_BCAST(vsink,    3*nsinkmax, MPI_DOUBLE_PRECISION, 1, MPI_COMM_WORLD, info)
  call MPI_BCAST(lsink,    3*nsinkmax, MPI_DOUBLE_PRECISION, 1, MPI_COMM_WORLD, info)
  call MPI_BCAST(delta_mass, nsinkmax, MPI_DOUBLE_PRECISION, 1, MPI_COMM_WORLD, info)
  call MPI_BCAST(idsink,     nsinkmax, MPI_INTEGER,          1, MPI_COMM_WORLD, info)
  call MPI_BCAST(tsink,      nsinkmax, MPI_DOUBLE_PRECISION, 1, MPI_COMM_WORLD, info)
  call MPI_BCAST(new_born,   nsinkmax, MPI_LOGICAL,          1, MPI_COMM_WORLD, info)

end subroutine synchronize_sink_info
#endif
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
