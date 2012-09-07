!################################################################
!################################################################
!################################################################
!################################################################
subroutine create_sink
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !----------------------------------------------------------------------------
  ! Description: This subroutine create sink particle in cells where a density 
  ! threshold has been crossed. It also removes from the gas the corresponding 
  ! particle mass. On exit, all fluid variables in the cell are modified.
  ! This routine is called only once per coarse step by routine amr_step.
  ! Romain Teyssier, October 7th, 2007 
  !
  ! The sink algorithm has been changed completely!
  !
  ! -If sink and clumpfind are both true, the new version is invoked.
  ! -Clumps which are good for sink production are flagged by the clumpfinder
  ! -The global sink variables are initialized, gas is accreted from the hostcell only.
  ! -One true RAMSES particle is created 
  ! -Sink cloud particles are created
  ! -Cloud particles are scattered to grid
  ! -Accretion routine is called
  !----------------------------------------------------------------------------
  ! local constants                                                                
  integer::ilevel,ivar,totparts
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m


  if(verbose)write(*,*)' Entering create_sink'

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_m=scale_d*scale_l**3


  !sink algorithm without clumpfinder
  if (clumpfind .eqv. .false.)then
     
     ! Remove particles to finer levels
     do ilevel=levelmin,nlevelmax
        call make_tree_fine(ilevel)
        call kill_tree_fine(ilevel)
        call virtual_tree_fine(ilevel)
     end do

     ! Create new sink particles
     ! and gather particle from the grid
     call make_sink(nlevelmax)
     do ilevel=nlevelmax-1,1,-1
        if(ilevel>=levelmin)call make_sink(ilevel)
        call merge_tree_fine(ilevel)
     end do
     
     ! Remove particle clouds around old sinks
     call kill_cloud(1)
     
     ! Merge sink using FOF 
     call merge_sink(1)
  end if

  !sink algorithm with clumpfinder
  if (clumpfind)then
     ! DO NOT USE FLAG2 BETWEEN CLUMP_FINDER AND MAKE_SINK_FROM_CLUMP
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     call clump_finder(.false.)                                  !!!!
     do ilevel=levelmin,nlevelmax                                !!!!
        call make_sink_from_clump(ilevel)                        !!!!
     end do                                                      !!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
     call create_part_from_sink
  end if

  ! Create new particle clouds
  call create_cloud(1)

! #ifndef WITHOUTMPI
!   call MPI_ALLREDUCE(npart,totparts,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
!   if(myid==1)print*,totparts
! #endif
  
  ! Scatter particle to the grid 
  do ilevel=1,nlevelmax
     call make_tree_fine(ilevel)
     call virtual_tree_fine(ilevel)
     call kill_tree_fine(ilevel)
     call virtual_tree_fine(ilevel)
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

  ! Compute Bondi parameters and gather particle
  do ilevel=nlevelmax,levelmin,-1
     if(bondi)then 
        call bondi_hoyle(ilevel)
     else 
        call grow_jeans(ilevel)
     end if
     call merge_tree_fine(ilevel)
  end do
  
  !update the cloud particle properties at levelmin
  call update_cloud(levelmin)

  call compute_accretion_rate(levelmin)
  
  if(agn)call agn_feedback
     
end subroutine create_sink
!################################################################
!################################################################
!################################################################
!################################################################
!################################################################
!################################################################
!################################################################
!################################################################
subroutine create_part_from_sink
  use amr_commons
  use pm_commons
  use hydro_commons
!  use cooling_module, ONLY: XH=>X, rhoc, mH 
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !----------------------------------------------------------------------
  ! Description: This subroutine create true RAMSES particles from the list
  ! of sink particles and stores them at level 1.
  !----------------------------------------------------------------------
  ! local constants
  integer ::icpu,index_sink
  integer ::i,isink
  integer ::indp
  integer ::ntot,ntot_all,info
  logical ::ok_free

  real(dp),dimension(1:nvector,1:ndim),save::xs
  integer ,dimension(1:nvector),save::ind_grid,cc
  integer ,dimension(1:nvector),save::ind_part
  logical ,dimension(1:nvector),save::ok_true=.true.
  integer ,dimension(1:ncpu)::ntot_sink_cpu,ntot_sink_all
 
 
  if(numbtot(1,1)==0) return
  if(.not. hydro)return
  if(ndim.ne.3)return

  if(verbose)write(*,*)' Entering create_part_from_sink'

#if NDIM==3

  ntot=0
  ! Loop over sinks
  do isink=1,nsink
     xs(1,1:ndim)=xsink(isink,1:ndim)
     call cmp_cpumap(xs,cc,1)
     if(cc(1).eq.myid)ntot=ntot+1
  end do

  !---------------------------------
  ! Check for free particle memory
  !--------------------------------
  ok_free=(numbp_free-ntot)>=0
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(numbp_free,numbp_free_tot,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  numbp_free_tot=numbp_free
#endif
   if(.not. ok_free)then
      write(*,*)'No more free memory for particles'
      write(*,*)'New sink particles',ntot
      write(*,*)'Increase npartmax'
#ifndef WITHOUTMPI
      call MPI_ABORT(MPI_COMM_WORLD,1,info)
#endif
#ifdef WITHOUTMPI
      stop
#endif
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
  if(myid==1)then
     if(ntot_all.gt.0)then
        write(*,'(" sinks to be created= ",I6," Tot =",I8)')&
             & ntot_all,nsink
     endif
  end if

  ! Starting identity number
  if(myid==1)then
     index_sink=nsink-ntot_all
  else
     index_sink=nsink-ntot_all+ntot_sink_cpu(myid-1)
  end if

  ! Level 1 linked list
  do icpu=1,ncpu
     if(numbl(icpu,1)>0)then
        ind_grid(1)=headl(icpu,1)
     endif
  end do

  !sort the sink according to mass
  if(nsink>0)then
     do i=1,nsink
        xmsink(i)=msink(i)
     end do
     call quick_sort(xmsink(1),idsink_sort(1),nsink)
  endif

  ! Loop over sinks
  do i=nsink,1,-1
     isink=idsink_sort(i)
     xs(1,1:ndim)=xsink(isink,1:ndim)
     call cmp_cpumap(xs,cc,1)

     ! Create new particles
     if(cc(1).eq.myid)then
           index_sink=index_sink+1
           ! Update linked list
           call remove_free(ind_part,1)
           call add_list(ind_part,ind_grid,ok_true,1)
           indp=ind_part(1)
           tp(indp)=tsink(isink)     ! Birth epoch
           mp(indp)=msink(isink)     ! Mass
           levelp(indp)=levelmin
           idp(indp)=-isink          ! Identity
           xp(indp,1)=xsink(isink,1) ! Position
           xp(indp,2)=xsink(isink,2)
           xp(indp,3)=xsink(isink,3)
           vp(indp,1)=vsink(isink,1) ! Velocity
           vp(indp,2)=vsink(isink,2)
           vp(indp,3)=vsink(isink,3)
        endif

  end do

#endif

end subroutine create_part_from_sink
!################################################################
!################################################################
!################################################################
!################################################################
subroutine create_cloud(ilevel)
  use pm_commons
  use amr_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !------------------------------------------------------------------------
  ! This routine creates a cloud of test particle around each sink particle.
  !------------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,next_part
  integer::ig,ip,npart1,npart2,icpu
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Gather sink particles only.

  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0
        
        ! Count sink particles
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
        
        ! Gather sink particles
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
                 call mk_cloud(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
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
     if(ip>0)call mk_cloud(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
  end do 
  ! End loop over cpus

111 format('   Entering create_cloud for level ',I2)

end subroutine create_cloud
!################################################################
!################################################################
!################################################################
!################################################################
subroutine mk_cloud(ind_part,ind_grid_part,np,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
  integer::np,ilevel
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !-----------------------------------------------------------------------
  ! This routine is called by subroutine create_cloud.
  !-----------------------------------------------------------------------
  integer::j,isink,ii,jj,kk,nx_loc,ncloud
  real(dp)::dx_loc,scale,dx_min,xx,yy,zz,rr,rmax
  ! Particle-based arrays
  integer ,dimension(1:nvector),save::ind_cloud
  logical ,dimension(1:nvector),save::ok_true=.true.

  ! Mesh spacing in that level
  dx_loc=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax/aexp

  rmax=dble(ir_cloud)*dx_min
  xx=0.0; yy=0.0;zz=0.0
  ncloud=0
  do kk=-2*ir_cloud,2*ir_cloud
     zz=dble(kk)*dx_min/2.0
     do jj=-2*ir_cloud,2*ir_cloud
        yy=dble(jj)*dx_min/2.0
        do ii=-2*ir_cloud,2*ir_cloud
           xx=dble(ii)*dx_min/2.0
           rr=sqrt(xx*xx+yy*yy+zz*zz)
           if(rr<=rmax)ncloud=ncloud+1
        end do
     end do
  end do
  ncloud_sink=ncloud

  do kk=-2*ir_cloud,2*ir_cloud
     zz=dble(kk)*dx_min/2.0
     do jj=-2*ir_cloud,2*ir_cloud
        yy=dble(jj)*dx_min/2.0
        do ii=-2*ir_cloud,2*ir_cloud
           xx=dble(ii)*dx_min/2.0
           rr=sqrt(xx*xx+yy*yy+zz*zz)
           if(rr>0.and.rr<=rmax)then
              call remove_free(ind_cloud,np)
              call add_list(ind_cloud,ind_grid_part,ok_true,np)
              do j=1,np
                 isink=-idp(ind_part(j))
                 idp(ind_cloud(j))=-isink
                 levelp(ind_cloud(j))=levelmin
                 mp(ind_cloud(j))=msink(isink)/dble(ncloud)
                 xp(ind_cloud(j),1)=xp(ind_part(j),1)+xx
                 vp(ind_cloud(j),1)=vsink(isink,1)
                 xp(ind_cloud(j),2)=xp(ind_part(j),2)+yy
                 vp(ind_cloud(j),2)=vsink(isink,2)
                 xp(ind_cloud(j),3)=xp(ind_part(j),3)+zz
                 vp(ind_cloud(j),3)=vsink(isink,3)
              end do
           end if
        end do
     end do
  end do

  ! Reduce sink particle mass
  do j=1,np
     isink=-idp(ind_part(j))
     mp(ind_part(j))=msink(isink)/dble(ncloud)
     vp(ind_part(j),1)=vsink(isink,1)
     vp(ind_part(j),2)=vsink(isink,2)
     vp(ind_part(j),3)=vsink(isink,3)
  end do

end subroutine mk_cloud
!################################################################
!################################################################
!################################################################
!################################################################
subroutine bondi_hoyle(ilevel)
  use pm_commons
  use amr_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !------------------------------------------------------------------------
  ! This routine computes the parameters of Bondi-Hoyle
  ! accretion for sink particles.
  ! It calls routine bondi_velocity and bondi_average.
  !------------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,next_part,idim,info
  integer::ig,ip,npart1,npart2,icpu,nx_loc,isink
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part
  real(dp)::r2,dx_loc,dx_min,scale,factG

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Gravitational constant
  factG=1d0
  if(cosmo)factG=3d0/8d0/3.1415926*omega_m*aexp

  ! Mesh spacing in that level
  dx_loc=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax/aexp

  ! Reset new sink variables
  v2sink_new=0d0; c2sink_new=0d0; oksink_new=0d0

  ! Gather sink particles only.

  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0

        ! Count only sink particles
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if(idp(ipart).lt.0)then
                 isink=-idp(ipart)
                 r2=0.0
                 do idim=1,ndim
                    r2=r2+(xp(ipart,idim)-xsink(isink,idim))**2
                 end do
                 if(r2==0.0)then
                    npart2=npart2+1
                 end if
              endif
              ipart=next_part  ! Go to next particle
           end do
        endif

        ! Gather only sink particles
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
                 isink=-idp(ipart)
                 r2=0.0
                 do idim=1,ndim
                    r2=r2+(xp(ipart,idim)-xsink(isink,idim))**2
                 end do
                 if(r2==0.0)then
                    if(ig==0)then
                       ig=1
                       ind_grid(ig)=igrid
                    end if
                    ip=ip+1
                    ind_part(ip)=ipart
                    ind_grid_part(ip)=ig
                 endif
              endif
              if(ip==nvector)then
                 call bondi_velocity(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
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
     if(ip>0)call bondi_velocity(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
  end do
  ! End loop over cpus

  if(nsink>0)then
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(oksink_new,oksink_all,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(c2sink_new,c2sink_all,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(v2sink_new,v2sink_all,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#else
     oksink_all=oksink_new
     c2sink_all=c2sink_new
     v2sink_all=v2sink_new
#endif
  endif

  do isink=1,nsink
     if(oksink_all(isink)==1d0)then
        c2sink(isink)=c2sink_all(isink)
        v2sink(isink)=v2sink_all(isink)
        ! Compute sink radius
        r2sink(isink)=(factG*msink(isink)/(v2sink(isink)+c2sink(isink)))**2
        r2k(isink)=min(max(r2sink(isink),(dx_min/4.0)**2),(2.*dx_min)**2)
     endif
  end do

  ! Gather sink and cloud particles.
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
                 call bondi_average(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
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
     if(ip>0)call bondi_average(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
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

111 format('   Entering bondi_hoyle for level ',I2)

end subroutine bondi_hoyle
!################################################################
!################################################################
!################################################################
!################################################################
subroutine bondi_velocity(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !-----------------------------------------------------------------------
  ! This routine is called by subroutine bondi_hoyle.
  ! It computes the gas velocity and soud speed in the cell
  ! each sink particle sits in.
  !-----------------------------------------------------------------------
  integer::i,j,idim,nx_loc,isink
  real(dp)::v2,c2,d,u,v,w,e,bx1,bx2,by1,by2,bz1,bz2
  real(dp)::dx,dx_loc,scale,vol_loc
  logical::error
  ! Grid based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle based arrays
  logical,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector,1:ndim),save::x
  integer ,dimension(1:nvector,1:ndim),save::id,igd,icd
  integer ,dimension(1:nvector),save::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc

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
        x(j,idim)=xp(ind_part(j),idim)/scale+skip_loc(idim)
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

  ! Check for illegal moves
  error=.false.
  do idim=1,ndim
     do j=1,np
        if(x(j,idim)<=0.0D0.or.x(j,idim)>=6.0D0)error=.true.
     end do
  end do
  if(error)then
     write(*,*)'problem in bondi_velocity'
     write(*,*)ilevel,ng,np
     stop
  end if

  ! NGP at level ilevel
  do idim=1,ndim
     do j=1,np
        id(j,idim)=x(j,idim)
     end do
  end do

   ! Compute parent grids
  do idim=1,ndim
     do j=1,np
        igd(j,idim)=id(j,idim)/2
     end do
  end do
  do j=1,np
     kg(j)=1+igd(j,1)+3*igd(j,2)+9*igd(j,3)
  end do
  do j=1,np
     igrid(j)=son(nbors_father_cells(ind_grid_part(j),kg(j)))
  end do

  ! Check if particles are entirely in level ilevel
  ok(1:np)=.true.
  do j=1,np
     ok(j)=ok(j).and.igrid(j)>0
  end do
  ! Compute parent cell position
  do idim=1,ndim
     do j=1,np
        if(ok(j))then
           icd(j,idim)=id(j,idim)-2*igd(j,idim)
        end if
     end do
  end do
  do j=1,np
     if(ok(j))then
        icell(j)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
     end if
  end do

  ! Compute parent cell adress
  do j=1,np
     if(ok(j))then
        indp(j)=ncoarse+(icell(j)-1)*ngridmax+igrid(j)
     end if
  end do

  ! Gather hydro variables
  do j=1,np
     if(ok(j))then
        d=uold(indp(j),1)
        u=uold(indp(j),2)/d
        v=uold(indp(j),3)/d
        w=uold(indp(j),4)/d
        e=uold(indp(j),5)/d
#ifdef SOLVERmhd
        bx1=uold(indp(j),6)
        by1=uold(indp(j),7)
        bz1=uold(indp(j),8)
        bx2=uold(indp(j),nvar+1)
        by2=uold(indp(j),nvar+2)
        bz2=uold(indp(j),nvar+3)
        e=e-0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)/d
#endif
        v2=(u**2+v**2+w**2)
        e=e-0.5d0*v2
        c2=MAX(gamma*(gamma-1.0)*e,smallc**2)
        isink=-idp(ind_part(j))
        v2=(u-vsink(isink,1))**2+(v-vsink(isink,2))**2+(w-vsink(isink,3))**2
        v2sink_new(isink)=v2
        c2sink_new(isink)=c2
        oksink_new(isink)=1d0
     endif
  end do

#endif

end subroutine bondi_velocity
!################################################################
!################################################################
!################################################################
!################################################################
subroutine bondi_average(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !-----------------------------------------------------------------------
  ! This routine is called by subroutine bondi_hoyle. Each cloud particle
  ! reads up the value of density, sound speed and velocity from its
  ! position in the grid.
  !-----------------------------------------------------------------------
  logical::error
  integer::i,j,ind,idim,nx_loc,isink
  real(dp)::d,u,v=0d0,w=0d0,e,bx1,bx2,by1,by2,bz1,bz2
  real(dp)::dx,scale,weight,r2
  ! Grid-based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle-based arrays
  logical ,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector),save::dgas,ugas,vgas,wgas,egas
  real(dp),dimension(1:nvector,1:ndim),save::x,dd,dg
  integer ,dimension(1:nvector,1:ndim),save::ig,id,igg,igd,icg,icd
  real(dp),dimension(1:nvector,1:twotondim),save::vol
  integer ,dimension(1:nvector,1:twotondim),save::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc

  ! Mesh spacing in that level
  dx=0.5D0**ilevel 
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)

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
        x(j,idim)=xp(ind_part(j),idim)/scale+skip_loc(idim)
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

  ! Check for illegal moves
  error=.false.
  do idim=1,ndim
     do j=1,np
        if(x(j,idim)<0.5D0.or.x(j,idim)>5.5D0)error=.true.
     end do
  end do
  if(error)then
     write(*,*)'problem in average_density'
     do idim=1,ndim
        do j=1,np
           if(x(j,idim)<0.5D0.or.x(j,idim)>5.5D0)then
              write(*,*)x(j,1:ndim)
           endif
        end do
     end do
     stop
  end if

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

  ! Check if particles are entirely in level ilevel
  ok(1:np)=.true.
  do ind=1,twotondim
     do j=1,np
        ok(j)=ok(j).and.igrid(j,ind)>0
     end do
  end do

  ! If not, rescale position at level ilevel-1
  do idim=1,ndim
     do j=1,np
        if(.not.ok(j))then
           x(j,idim)=x(j,idim)/2.0D0
        end if
     end do
  end do
  ! If not, redo CIC at level ilevel-1
  do idim=1,ndim
     do j=1,np
        if(.not.ok(j))then
           dd(j,idim)=x(j,idim)+0.5D0
           id(j,idim)=dd(j,idim)
           dd(j,idim)=dd(j,idim)-id(j,idim)
           dg(j,idim)=1.0D0-dd(j,idim)
           ig(j,idim)=id(j,idim)-1
        end if
     end do
  end do

 ! Compute parent cell position
  do idim=1,ndim
     do j=1,np
        if(ok(j))then
           icg(j,idim)=ig(j,idim)-2*igg(j,idim)
           icd(j,idim)=id(j,idim)-2*igd(j,idim)
        else
           icg(j,idim)=ig(j,idim)
           icd(j,idim)=id(j,idim)
        end if
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
     if(ok(j))then
        icell(j,1)=1+icg(j,1)+2*icg(j,2)
        icell(j,2)=1+icd(j,1)+2*icg(j,2)
        icell(j,3)=1+icg(j,1)+2*icd(j,2)
        icell(j,4)=1+icd(j,1)+2*icd(j,2)
     else
        icell(j,1)=1+icg(j,1)+3*icg(j,2)
        icell(j,2)=1+icd(j,1)+3*icg(j,2)
        icell(j,3)=1+icg(j,1)+3*icd(j,2)
        icell(j,4)=1+icd(j,1)+3*icd(j,2)
     end if
  end do
#endif
#if NDIM==3
  do j=1,np
     if(ok(j))then
        icell(j,1)=1+icg(j,1)+2*icg(j,2)+4*icg(j,3)
        icell(j,2)=1+icd(j,1)+2*icg(j,2)+4*icg(j,3)
        icell(j,3)=1+icg(j,1)+2*icd(j,2)+4*icg(j,3)
        icell(j,4)=1+icd(j,1)+2*icd(j,2)+4*icg(j,3)
        icell(j,5)=1+icg(j,1)+2*icg(j,2)+4*icd(j,3)
        icell(j,6)=1+icd(j,1)+2*icg(j,2)+4*icd(j,3)
        icell(j,7)=1+icg(j,1)+2*icd(j,2)+4*icd(j,3)
        icell(j,8)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
     else
        icell(j,1)=1+icg(j,1)+3*icg(j,2)+9*icg(j,3)
        icell(j,2)=1+icd(j,1)+3*icg(j,2)+9*icg(j,3)
        icell(j,3)=1+icg(j,1)+3*icd(j,2)+9*icg(j,3)
        icell(j,4)=1+icd(j,1)+3*icd(j,2)+9*icg(j,3)
        icell(j,5)=1+icg(j,1)+3*icg(j,2)+9*icd(j,3)
        icell(j,6)=1+icd(j,1)+3*icg(j,2)+9*icd(j,3)
        icell(j,7)=1+icg(j,1)+3*icd(j,2)+9*icd(j,3)
        icell(j,8)=1+icd(j,1)+3*icd(j,2)+9*icd(j,3)   
     end if
  end do
#endif
        
  ! Compute parent cell adresses
  do ind=1,twotondim
     do j=1,np
        if(ok(j))then
           indp(j,ind)=ncoarse+(icell(j,ind)-1)*ngridmax+igrid(j,ind)
        else
           indp(j,ind)=nbors_father_cells(ind_grid_part(j),icell(j,ind))
        end if
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

  dgas(1:np)=0.0D0  ! Gather gas density
  ugas(1:np)=0.0D0  ! Gather gas x-momentum (velocity ?)
  vgas(1:np)=0.0D0  ! Gather gas y-momentum (velocity ?)
  wgas(1:np)=0.0D0  ! Gather gas z-momentum (velocity ?)
  egas(1:np)=0.0D0  ! Gather gas thermal energy (specific ?)

  ! ROM to AJC: if you want, replace below vol(j,ind) by 1./twotondim
  do ind=1,twotondim
     do j=1,np
        d=uold(indp(j,ind),1)
        u=uold(indp(j,ind),2)/d
        v=uold(indp(j,ind),3)/d
        w=uold(indp(j,ind),4)/d
        e=uold(indp(j,ind),5)/d
#ifdef SOLVERmhd
        bx1=uold(indp(j,ind),6)
        by1=uold(indp(j,ind),7)
        bz1=uold(indp(j,ind),8)
        bx2=uold(indp(j,ind),nvar+1)
        by2=uold(indp(j,ind),nvar+2)
        bz2=uold(indp(j,ind),nvar+3)
        e=e-0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)/d
#endif
        e=e-0.5*(u*u+v*v+w*w)
        dgas(j)=dgas(j)+d*vol(j,ind)
        ugas(j)=ugas(j)+d*u*vol(j,ind)
        vgas(j)=vgas(j)+d*v*vol(j,ind)
        wgas(j)=wgas(j)+d*w*vol(j,ind)
        egas(j)=egas(j)+d*e*vol(j,ind)
     end do
  end do

  do j=1,np
     isink=-idp(ind_part(j))
     r2=0d0
     do idim=1,ndim
        r2=r2+(xp(ind_part(j),idim)-xsink(isink,idim))**2
     end do
     weight=exp(-r2/r2k(isink))
     wden(isink)=wden(isink)+weight*dgas(j)
     wmom(isink,1)=wmom(isink,1)+weight*ugas(j)
     wmom(isink,2)=wmom(isink,2)+weight*vgas(j)
     wmom(isink,3)=wmom(isink,3)+weight*wgas(j)
     weth(isink)=weth(isink)+weight*egas(j)
     wvol(isink)=wvol(isink)+weight
  end do

end subroutine bondi_average
!################################################################
!################################################################
!################################################################
!################################################################
subroutine grow_bondi(ilevel)
  use pm_commons
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !------------------------------------------------------------------------
  ! This routine performs Bondi-Hoyle accretion of the gas onto 
  ! sink particles. On exit, sink mass and velocity are modified.
  !------------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,next_part,info,iskip,ind
  integer::i,ig,ip,npart1,npart2,icpu,isink
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part


  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Compute sink accretion rates
  call compute_accretion_rate(0)

  ! Reset new sink variables                                                                        
  msink_new=0d0; vsink_new=0d0; delta_mass_new=0d0

  ! Store initial gas density in unew(:,1)
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,active(ilevel)%ngrid
        unew(active(ilevel)%igrid(i)+iskip,1) = uold(active(ilevel)%igrid(i)+iskip,1)
     enddo
  enddo
     
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
                 call accrete_bondi(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
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
     if(ip>0)call accrete_bondi(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
  end do 
  ! End loop over cpus

  if(nsink>0)then
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(msink_new,msink_all,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(delta_mass_new,delta_mass_all,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(vsink_new,vsink_all,nsinkmax*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#else
     msink_all=msink_new
     delta_mass_all=delta_mass_new
     vsink_all=vsink_new
#endif
  endif
  do isink=1,nsink
     vsink(isink,1:ndim)=vsink(isink,1:ndim)*msink(isink)+vsink_all(isink,1:ndim)
     msink(isink)=msink(isink)+msink_all(isink)
     delta_mass(isink)=delta_mass(isink)+delta_mass_all(isink)
     vsink(isink,1:ndim)=vsink(isink,1:ndim)/msink(isink)
  end do

111 format('   Entering grow_bondi for level ',I2)

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
end subroutine grow_bondi
!################################################################
!################################################################
!################################################################
!################################################################
subroutine accrete_bondi(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !-----------------------------------------------------------------------
  ! This routine is called by subroutine bondi_hoyle.
  !-----------------------------------------------------------------------
  integer::i,j,idim,nx_loc,isink,ivar
  real(dp)::r2,d,u,v,w,e,bx1,bx2,by1,by2,bz1,bz2,dini
  real(dp),dimension(1:nvar)::z
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::dx,dx_loc,scale,vol_loc,weight,acc_mass
  logical::error
  ! Grid based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle based arrays
  logical,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector,1:ndim),save::x
  integer ,dimension(1:nvector,1:ndim),save::id,igd,icd
  integer ,dimension(1:nvector),save::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc

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
        x(j,idim)=xp(ind_part(j),idim)/scale+skip_loc(idim)
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

  ! Check for illegal moves
  error=.false.
  do idim=1,ndim
     do j=1,np
        if(x(j,idim)<=0.0D0.or.x(j,idim)>=6.0D0)error=.true.
     end do
  end do
  if(error)then
     write(*,*)'problem in accrete_bondi'
     write(*,*)ilevel,ng,np
     stop
  end if

  ! NGP at level ilevel
  do idim=1,ndim
     do j=1,np
        id(j,idim)=x(j,idim)
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

  ! Check if particles are entirely in level ilevel
  ok(1:np)=.true.
  do j=1,np
     ok(j)=ok(j).and.igrid(j)>0
  end do

  ! Compute parent cell position
  do idim=1,ndim
     do j=1,np
        if(ok(j))then
           icd(j,idim)=id(j,idim)-2*igd(j,idim)
        end if
     end do
  end do
#if NDIM==1
  do j=1,np
     if(ok(j))then
        icell(j)=1+icd(j,1)
     end if
  end do
#endif
#if NDIM==2
  do j=1,np
     if(ok(j))then
        icell(j)=1+icd(j,1)+2*icd(j,2)
     end if
  end do
#endif
#if NDIM==3
  do j=1,np
     if(ok(j))then
        icell(j)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
     end if
  end do
#endif
        
  ! Compute parent cell adress
  do j=1,np
     if(ok(j))then
        indp(j)=ncoarse+(icell(j)-1)*ngridmax+igrid(j)
     end if
  end do

  ! Remove mass from hydro cells
  do j=1,np
     if(ok(j))then
        isink=-idp(ind_part(j))
        r2=0d0
        do idim=1,ndim
           r2=r2+(xp(ind_part(j),idim)-xsink(isink,idim))**2
        end do
        weight=exp(-r2/r2k(isink))
           
        d=uold(indp(j),1)
        dini=unew(indp(j),1) ! Initial density
        u=uold(indp(j),2)/d
        v=uold(indp(j),3)/d
        w=uold(indp(j),4)/d
        e=uold(indp(j),5)/d
#ifdef SOLVERmhd
        bx1=uold(indp(j),6)
        by1=uold(indp(j),7)
        bz1=uold(indp(j),8)
        bx2=uold(indp(j),nvar+1)
        by2=uold(indp(j),nvar+2)
        bz2=uold(indp(j),nvar+3)
        e=e-0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)/d
#endif
        do ivar=imetal,nvar
           z(ivar)=uold(indp(j),ivar)/d
        end do
        ! Compute accreted mass with cloud weighting
        acc_mass=dMBHoverdt(isink)*weight/total_volume(isink)*dtnew(ilevel)
        ! Cannot accrete more than 25% of initial gass mass in the cell
        acc_mass=max(min(acc_mass,(d-0.75*dini)*vol_loc),0.0_dp)
        msink_new(isink)=msink_new(isink)+acc_mass
        delta_mass_new(isink)=delta_mass_new(isink)+acc_mass
        vsink_new(isink,1)=vsink_new(isink,1)+acc_mass*u
        vsink_new(isink,2)=vsink_new(isink,2)+acc_mass*v
        vsink_new(isink,3)=vsink_new(isink,3)+acc_mass*w
        ! Remove accreted mass
        d=d-acc_mass/vol_loc
#ifdef SOLVERmhd
        e=e+0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)/d
#endif
        uold(indp(j),1)=d
        uold(indp(j),2)=d*u
        uold(indp(j),3)=d*v
        uold(indp(j),4)=d*w
        uold(indp(j),5)=d*e
        do ivar=imetal,nvar
           uold(indp(j),ivar)=d*z(ivar)
        end do
     endif

  end do

end subroutine accrete_bondi
!################################################################
!################################################################
!################################################################
!################################################################
subroutine compute_accretion_rate(ilevel)
  use pm_commons
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !------------------------------------------------------------------------
  ! This routine computes the accretion rate on the sink particles.
  !------------------------------------------------------------------------
  integer::i,nx_loc,isink
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m
  real(dp)::factG,d_star,boost,vel_max,l_abs,rot_period
  real(dp)::r2,v2,c2,density,volume,ethermal,dx_min,scale
  real(dp),dimension(1:3)::velocity

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

  ! Maximum relative velocity
  vel_max=10. ! in km/sec
  vel_max=vel_max*1d5/scale_v
  
  if(smbh)then

     ! Compute sink particle accretion rate
     do isink=1,nsink
        density=0d0
        volume=0d0
        velocity=0d0
        ethermal=0d0
        ! Loop over level: sink cloud can overlap several levels
        do i=levelmin,nlevelmax
           density=density+weighted_density(isink,i)
           ethermal=ethermal+weighted_ethermal(isink,i)
           velocity(1)=velocity(1)+weighted_momentum(isink,i,1)
           velocity(2)=velocity(2)+weighted_momentum(isink,i,2)
           velocity(3)=velocity(3)+weighted_momentum(isink,i,3)
           volume=volume+weighted_volume(isink,i)
        end do
        density=density/volume
        velocity(1:3)=velocity(1:3)/density/volume
        ethermal=ethermal/density/volume
        total_volume(isink)=volume
        c2=MAX(gamma*(gamma-1.0)*ethermal,smallc**2)
        v2=min(SUM((velocity(1:3)-vsink(isink,1:3))**2),vel_max**2)
        r2=(factG*msink(isink)/(c2+v2))**2
        
        ! Compute Bondi-Hoyle accretion rate in code units
        boost=1.0
        if(star)boost=max((density/d_star)**2,1.0_dp)
        !     dMBHoverdt(isink)=boost*4.*3.1415926*density*r2*sqrt(1.12**2*c2+v2)/bondi_alpha(1.2*dx_min/sqrt(r2))
        dMBHoverdt(isink)=boost*4.*3.1415926*density*r2*sqrt(c2+v2)
        
        ! Compute Eddington accretion rate in code units
        dMEDoverdt(isink)=4.*3.1415926*6.67d-8*msink(isink)*1.66d-24/(0.1*6.652d-25*3d10)*scale_t
        
        ! Compute BH radius of influence in kpc
        rBH(isink)=250.*(msink(isink)*scale_m/(5d9*2d33))**(0.5)
        
        ! Compute BH radius of influence in code units
        epsBH(isink)=max(4.*dx_min,rBH(isink)*3.08d21/scale_l)
        
     end do
     
     if(myid==1.and.ilevel==levelmin.and.nsink>0)then
        do i=1,nsink
           xmsink(i)=msink(i)
        end do
        call quick_sort(xmsink(1),idsink_sort(1),nsink)
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
     endif
     
     ! Take the minimum accretion rate
     do isink=1,nsink
        dMBHoverdt(isink)=min(dMBHoverdt(isink),dMEDoverdt(isink))
     end do
     
  else
     if (ilevel==levelmin)then    
        acc_rate(1:nsink)=acc_rate(1:nsink)/dtnew(levelmin)

        if(ir_feedback)then
           do i=1,nsink ! 0.75 and 5 are ratio of infalling energy which is radiated and protostellar radius
              acc_lum(i)=0.75*acc_rate(i)*msink(i)/(5*6.955d10/scale_l)
           end do
        end if

        if(myid==1.and.nsink>0.and. mod(nstep_coarse,ncontrol)==0)then
           do i=1,nsink
              xmsink(i)=msink(i)
           end do
           
           call quick_sort(xmsink(1),idsink_sort(1),nsink)
           write(*,*)'Number of sink = ',nsink
           write(*,'(" ====================================================================================================================================================== ")')
           write(*,'("  Id     M[Msol]    x           y           z           vx        vy        vz     rot_period[y] lx/|l|  ly/|l|  lz/|l| acc_rate[Msol/y] acc_lum[Lsol]  ")')
           write(*,'(" ====================================================================================================================================================== ")')
           do i=nsink,1,-1
              isink=idsink_sort(i)
              l_abs=(lsink(isink,1)**2+lsink(isink,2)**2+lsink(isink,3)**2)**0.5+1.d-99
              rot_period=32*3.1415*msink(isink)*(dx_min)**2/(5*l_abs)
              write(*,'(I5,2X,F9.5,3(2X,F10.8),3(2X,F7.4),2X,F13.5,3(2X,F6.3),2X,E11.3,4x,E11.3)')idsink(isink),msink(isink)*scale_m/2d33, &
                   xsink(isink,1:ndim),vsink(isink,1:ndim),&
                   rot_period*scale_t/(3600*24*365),lsink(isink,1)/l_abs,lsink(isink,2)/l_abs,lsink(isink,3)/l_abs,&
                   acc_rate(isink)*scale_m/2.d33/(scale_t)*365.*24.*3600.,acc_lum(isink)/scale_t**2*scale_l**3*scale_d*scale_l**2/scale_t/3.933d33
              
           end do
           write(*,'(" ====================================================================================================================================================== ")')
        endif
        !acc_rate=0. taken away because accretion rate must not be set to 0 before dump_all! now in amr_step just after dump_all
     end if
  endif

end subroutine compute_accretion_rate
!################################################################
!################################################################
!################################################################
!################################################################
subroutine grow_jeans(ilevel)
  use pm_commons
  use amr_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !------------------------------------------------------------------------
  ! This routine checks wether a cell hosting a cloud particle has crossed
  ! the accretion threshold and accretes mass, center of mass, momentum, angular momentum
  !------------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,next_part,info,lev
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  integer::ig,ip,npart1,npart2,icpu,isink
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  
  if(numbtot(1,ilevel)==0)return
  
  if(verbose)write(*,111)ilevel



  ! Reset new sink variables
  msink_new=0d0; vsink_new=0d0; lsink_new=0d0; xsink_new=0.d0; 

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
                 call accrete_jeans(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
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
     if(ip>0)call accrete_jeans(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
  end do 
  ! End loop over cpus

  if(nsink>0)then
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(msink_new,msink_all,nsinkmax     ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(xsink_new,xsink_all,nsinkmax*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(vsink_new,vsink_all,nsinkmax*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(lsink_new,lsink_all,nsinkmax*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#else
     xsink_all=xsink_new
     msink_all=msink_new
     vsink_all=vsink_new
     lsink_all=lsink_new
#endif
  endif
  
  do isink=1,nsink
     do lev=levelmin,nlevelmax
        sink_jump(isink,1:ndim,lev)=sink_jump(isink,1:ndim,lev)-xsink(isink,1:ndim)
     end do
     !change to conservative quantities
     xsink(isink,1:ndim)=xsink(isink,1:ndim)*msink(isink)
     vsink(isink,1:ndim)=vsink(isink,1:ndim)*msink(isink)
     !accrete to sink variables
     msink(isink)=msink(isink)+msink_all(isink)
     xsink(isink,1:ndim)=xsink(isink,1:ndim)+xsink_all(isink,1:ndim)
     vsink(isink,1:ndim)=vsink(isink,1:ndim)+vsink_all(isink,1:ndim)
     lsink(isink,1:3)=lsink(isink,1:3)+lsink_all(isink,1:3)
     !change back     
     xsink(isink,1:ndim)=xsink(isink,1:ndim)/msink(isink)
     vsink(isink,1:ndim)=vsink(isink,1:ndim)/msink(isink)
     do lev=levelmin,nlevelmax
        sink_jump(isink,1:ndim,lev)=sink_jump(isink,1:ndim,lev)+xsink(isink,1:ndim)
     end do
     acc_rate(isink)=acc_rate(isink)+msink_all(isink)
  end do



111 format('   Entering grow_jeans for level ',I2)

end subroutine grow_jeans
!################################################################
!################################################################
!################################################################
!################################################################
subroutine accrete_jeans(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !-----------------------------------------------------------------------
  ! This routine is called by grow_jeans
  !-----------------------------------------------------------------------
  integer::i,j,idim,nx_loc,isink,ivar,ind,ix,iy,iz
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,factG
  real(dp)::v2,d,u,v=0d0,w=0d0,e,bx1,bx2,by1,by2,bz1,bz2,delta_d
  real(dp),dimension(1:nvar)::z
  real(dp)::dx,dx_loc,scale,vol_loc,temp,d_jeans,acc_mass,d_sink,d_thres
  logical::error
  ! Grid based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle based arrays
  logical,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector,1:ndim),save::x
  integer ,dimension(1:nvector,1:ndim),save::id,igd,icd
  integer ,dimension(1:nvector),save::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc,xx
  real(dp),dimension(1:twotondim,1:3)::xc

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

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

  ! Cells center position relative to grid center position               
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     xc(ind,1)=(dble(ix)-0.5D0)*dx
     xc(ind,2)=(dble(iy)-0.5D0)*dx
     xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  ! Density threshold for sink particle formation
  d_sink=n_sink/scale_nH

#if NDIM==3
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
        x(j,idim)=xp(ind_part(j),idim)/scale+skip_loc(idim)
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

  ! Check for illegal moves
  error=.false.
  do idim=1,ndim
     do j=1,np
        if(x(j,idim)<=0.0D0.or.x(j,idim)>=6.0D0)error=.true.
     end do
  end do
  if(error)then
     write(*,*)'problem in accrete_jeans'
     write(*,*)ilevel,ng,np
     stop
  end if

  ! NGP at level ilevel
  do idim=1,ndim
     do j=1,np
        id(j,idim)=x(j,idim)
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

  ! Check if particles are entirely in level ilevel
  ok(1:np)=.true.
  do j=1,np
     ok(j)=ok(j).and.igrid(j)>0
  end do


  ! Compute parent cell position
  do idim=1,ndim
     do j=1,np
        if(ok(j))then
           icd(j,idim)=id(j,idim)-2*igd(j,idim)
        end if
     end do
  end do
#if NDIM==1
  do j=1,np
     if(ok(j))then
        icell(j)=1+icd(j,1)
     end if
  end do
#endif
#if NDIM==2
  do j=1,np
     if(ok(j))then
        icell(j)=1+icd(j,1)+2*icd(j,2)
     end if
  end do
#endif
#if NDIM==3
  do j=1,np
     if(ok(j))then
        icell(j)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
     end if
  end do
#endif
        
  ! Compute parent cell adress
  do j=1,np
     if(ok(j))then
        indp(j)=ncoarse+(icell(j)-1)*ngridmax+igrid(j)
     end if
  end do

  ! Check if particles are in a leaf cell
  do j=1,np
     ok(j)=ok(j).and.son(indp(j))==0
  end do


  ! Gather hydro variables
  do j=1,np
     if(ok(j))then

        !get cell center positions
        xx(1)=(x0(ind_grid_part(j),1)+3.0D0*dx+xc(icell(j),1)-skip_loc(1))*scale
        xx(2)=(x0(ind_grid_part(j),2)+3.0D0*dx+xc(icell(j),2)-skip_loc(2))*scale
        xx(3)=(x0(ind_grid_part(j),3)+3.0D0*dx+xc(icell(j),3)-skip_loc(3))*scale

        ! Convert uold to primitive variables
        d=uold(indp(j),1)
        u=uold(indp(j),2)/d
        v=uold(indp(j),3)/d
        w=uold(indp(j),4)/d
        e=uold(indp(j),5)/d
#ifdef SOLVERmhd
        bx1=uold(indp(j),6)
        by1=uold(indp(j),7)
        bz1=uold(indp(j),8)
        bx2=uold(indp(j),nvar+1)
        by2=uold(indp(j),nvar+2)
        bz2=uold(indp(j),nvar+3)
        e=e-0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)/d
#endif
        v2=(u**2+v**2+w**2)
        e=e-0.5d0*v2
        do ivar=imetal,nvar
           z(ivar)=uold(indp(j),ivar)/d
        end do


        ! User defined density threshold 
        d_thres=d_sink

        ! Jeans length related density threshold
        if(d_thres<0.0)then
           temp=max(e*(gamma-1.0),smallc**2)
           d_jeans=temp*3.1415926/(4.0*dx_loc)**2/factG
           d_thres=d_jeans
        endif

        if(d>d_thres)then
           isink=-idp(ind_part(j))
           delta_d=(d-d_thres)*0.0125
           acc_mass=delta_d*vol_loc
           msink_new(isink  )=msink_new(isink  )+acc_mass
           xsink_new(isink,1:ndim)=xsink_new(isink,1:ndim)+acc_mass*xx(1:ndim)
           vsink_new(isink,1)=vsink_new(isink,1)+acc_mass*u
           vsink_new(isink,2)=vsink_new(isink,2)+acc_mass*v
           vsink_new(isink,3)=vsink_new(isink,3)+acc_mass*w
           lsink_new(isink,1)=((xx(2)-xsink(isink,2))*(w-vsink(isink,3)/msink(isink))&
                -(xx(3)-xsink(isink,3))*(v-vsink(isink,2)/msink(isink)))*acc_mass
           lsink_new(isink,2)=((xx(3)-xsink(isink,3))*(u-vsink(isink,1)/msink(isink))&
                -(xx(1)-xsink(isink,1))*(w-vsink(isink,3)/msink(isink)))*acc_mass
           lsink_new(isink,3)=((xx(1)-xsink(isink,1))*(v-vsink(isink,2)/msink(isink))&
                -(xx(2)-xsink(isink,2))*(u-vsink(isink,1)/msink(isink)))*acc_mass


           d=d-delta_d
           ! Convert back to conservative variable                                             
#ifdef SOLVERmhd
           e=e+0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)/d
#endif
           e=e+0.5d0*(u**2+v**2+w**2)
           uold(indp(j),1)=d
           uold(indp(j),2)=d*u
           uold(indp(j),3)=d*v
           uold(indp(j),4)=d*w
           uold(indp(j),5)=d*e
           do ivar=imetal,nvar
              uold(indp(j),ivar)=d*z(ivar)
           end do


!            !remove the gas from unew
!            unew(indp(j),1)=unew(indp(j),1)-delta_d
!            unew(indp(j),2)=unew(indp(j),2)-delta_d*u
!            unew(indp(j),3)=unew(indp(j),3)-delta_d*v
!            unew(indp(j),4)=unew(indp(j),4)-delta_d*w
!            unew(indp(j),5)=unew(indp(j),5)-delta_d*(e+0.5d0*v2)
           
! #ifdef SOLVERmhd
!            unew(indp(j),6)=unew(indp(j),6)-delta_d*bx1
!            unew(indp(j),7)=unew(indp(j),7)-delta_d*by1
!            unew(indp(j),8)=unew(indp(j),8)-delta_d*bz1
!            unew(indp(j),nvar+1)=unew(indp(j),nvar+1)-delta_d*bx1
!            unew(indp(j),nvar+2)=unew(indp(j),nvar+2)-delta_d*by1
!            unew(indp(j),nvar+3)=unew(indp(j),nvar+3)-delta_d*bz1
! #endif
!            do ivar=imetal,nvar
!               unew(indp(j),ivar)=unew(indp(j),ivar)-delta_d*z(ivar)
!            end do
        endif
     endif
  end do
  
#endif
  
end subroutine accrete_jeans
!################################################################
!################################################################
!################################################################
!################################################################
subroutine make_sink_from_clump(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
!  use cooling_module, ONLY: XH=>X, rhoc, mH 
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
  real(dp),dimension(1:twotondim,1:3)::xc
  integer ::ncache,nnew,ivar,ngrid,icpu,index_sink,index_sink_tot
  integer ::igrid,ix,iy,iz,ind,i,iskip,isink,nx_loc
  integer ::ii,jj,kk,ncloud
  integer ::ntot,ntot_all,info
  logical ::ok_free

  real(dp),dimension(1:nvar)::z

  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m
  real(dp)::d,u,v,w,e,factG,delta_d,v2
  real(dp)::msink_max2,rsink_max2
  real(dp)::rmax,rmax2
  real(dp)::d_thres,d_sink
  real(dp)::birth_epoch,xx,yy,zz,rr
  real(dp),dimension(1:3)::skip_loc,x
  real(dp)::dx,dx_loc,scale,vol_loc,dx_min,vol_min
  real(dp)::bx1,bx2,by1,by2,bz1,bz2

  integer ,dimension(1:nvector),save::ind_grid,ind_cell
  integer ,dimension(1:nvector),save::ind_grid_new,ind_cell_new
  integer ,dimension(1:ncpu)::ntot_sink_cpu,ntot_sink_all
  
  if(.not. hydro)return
  if(ndim.ne.3)return

  if(verbose)write(*,*)'entering make_sink_from_clump for level ',ilevel

  ! Conversion factor from user units to cgs units                              
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_m=scale_d*scale_l**3

  ! Minimum radius to create a new sink from any other
  rsink_max=10.0 ! in kpc
  rsink_max2=(rsink_max*3.08d21/scale_l)**2

  ! Maximum value for the initial sink mass
  msink_max=1d5 ! in Msol
  msink_max2=msink_max*2d33/scale_m
  
  ! Gravitational constant
  factG=1d0
  if(cosmo)factG=3d0/8d0/3.1415926*omega_m*aexp

  ! Density threshold for sink particle creation
  d_sink=n_sink/scale_nH

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
  dx_min=(0.5D0**nlevelmax)*scale
  vol_min=dx_min**ndim

  rmax=dble(ir_cloud)*dx_min/aexp ! Linking length in physical units
  rmax2=rmax*rmax

  ! Birth epoch
  birth_epoch=t

  ! Cells center position relative to grid center position
  do ind=1,twotondim  
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     xc(ind,1)=(dble(ix)-0.5D0)*dx
     xc(ind,2)=(dble(iy)-0.5D0)*dx
     xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  xx=0.0; yy=0.0;zz=0.0
  ncloud=0
  do kk=-2*ir_cloud,2*ir_cloud
     zz=dble(kk)*0.5
     do jj=-2*ir_cloud,2*ir_cloud
        yy=dble(jj)*0.5
        do ii=-2*ir_cloud,2*ir_cloud
           xx=dble(ii)*0.5
           rr=sqrt(xx*xx+yy*yy+zz*zz)
           if(rr<=dble(ir_cloud))ncloud=ncloud+1
        end do
     end do
  end do
  ncloud_sink=ncloud

  ! Set new sink variables to zero
  msink_new=0d0; tsink_new=0d0; delta_mass_new=0d0; xsink_new=0d0; vsink_new=0d0
  oksink_new=0d0; idsink_new=0; new_born=0; level_sink_new=0

#if NDIM==3

  !------------------------------------------------
  ! Convert hydro variables to primitive variables
  ! and count number of new sinks
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
              ntot=ntot+flag2(ind_cell(i))
           end do
        end do
     end do
  end if

  !---------------------------------
  ! Check for free particle memory
  !--------------------------------
  ok_free=(numbp_free-ntot)>=0
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(numbp_free,numbp_free_tot,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  numbp_free_tot=numbp_free
#endif
   if(.not. ok_free)then
      write(*,*)'No more free memory for particles'
      write(*,*)'New sink particles',ntot
      write(*,*)'Increase npartmax'
#ifndef WITHOUTMPI
      call MPI_ABORT(MPI_COMM_WORLD,1,info)
#endif
#ifdef WITHOUTMPI
      stop
#endif
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
              e=uold(ind_cell_new(i),5)/d
#ifdef SOLVERmhd
              bx1=uold(ind_cell_new(i),6)
              by1=uold(ind_cell_new(i),7)
              bz1=uold(ind_cell_new(i),8)
              bx2=uold(ind_cell_new(i),nvar+1)
              by2=uold(ind_cell_new(i),nvar+2)
              bz2=uold(ind_cell_new(i),nvar+3)
              e=e-0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)/d
#endif
              v2=(u**2+v**2+w**2)
              e=e-0.5d0*v2
              do ivar=imetal,nvar
                 z(ivar)=uold(ind_cell_new(i),ivar)/d
              end do
              

              
              ! Get density maximum position
              x(1)=(xg(ind_grid_new(i),1)+xc(ind,1)-skip_loc(1))*scale
              x(2)=(xg(ind_grid_new(i),2)+xc(ind,2)-skip_loc(2))*scale
              x(3)=(xg(ind_grid_new(i),3)+xc(ind,3)-skip_loc(3))*scale              
              call true_max(x(1),x(2),x(3),ilevel)

              ! User defined density threshold
              d_thres=d_sink

              ! Mass of the new sink
              delta_d=(d-d_thres)*0.0125
              if (d>d_thres)then
                 msink_new(index_sink)=delta_d*vol_loc
              else
                 write(*,*)'sink production with negative mass'
                 call clean_stop
              endif

              delta_mass_new(index_sink)=0d0

              ! Global index of the new sink
              oksink_new(index_sink)=1d0
              idsink_new(index_sink)=index_sink_tot

              ! Store properties of the new sink
              tsink_new(index_sink)=birth_epoch
              xsink_new(index_sink,1:3)=x(1:3)
              vsink_new(index_sink,1)=u
              vsink_new(index_sink,2)=v
              vsink_new(index_sink,3)=w
              new_born(index_sink)=1
              level_sink_new(index_sink)=0 ! Newly created sinks are time-synchronized

!              unew(ind_cell_new(i),1)=unew(ind_cell_new(i),1)-msink_new(index_sink)/vol_loc

              d=d-delta_d
              ! Convert back to conservative variable                                             
#ifdef SOLVERmhd
              e=e+0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)/d
#endif
              e=e+0.5d0*(u**2+v**2+w**2)
              uold(ind_cell_new(i),1)=d
              uold(ind_cell_new(i),2)=d*u
              uold(ind_cell_new(i),3)=d*v
              uold(ind_cell_new(i),4)=d*w
              uold(ind_cell_new(i),5)=d*e
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
  call MPI_ALLREDUCE(tsink_new ,tsink_all ,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(xsink_new ,xsink_all ,nsinkmax*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(vsink_new ,vsink_all ,nsinkmax*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(delta_mass_new,delta_mass_all,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(new_born,new_born_all,nsinkmax,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(level_sink_new,level_sink_all,nsinkmax,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#else
  oksink_all=oksink_new
  idsink_all=idsink_new
  msink_all=msink_new
  tsink_all=tsink_new
  xsink_all=xsink_new
  vsink_all=vsink_new
  delta_mass_all=delta_mass_new
  new_born_all=new_born
  level_sink_all=level_sink_new
#endif
  do isink=1,nsink
     if(oksink_all(isink)==1)then
        idsink(isink)=idsink_all(isink)
        msink(isink)=msink_all(isink)
        tsink(isink)=tsink_all(isink)
        xsink(isink,1:ndim)=xsink_all(isink,1:ndim)
        vsink(isink,1:ndim)=vsink_all(isink,1:ndim)
        delta_mass(isink)=delta_mass_all(isink)
        level_sink(isink)=level_sink_all(isink)
        acc_rate(isink)=msink_all(isink)
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
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  real(dp)::x,y,z
  integer::ilevel
  !----------------------------------------------------------------------------
  ! Description: This subroutine takes the cell of maximum density and computes
  ! the maximum by expanding the density around the cell center to second order.
  !----------------------------------------------------------------------------
  ! local constants                                                                
  integer::k,j,i,nx_loc
  real(dp)::det,dx,dx_loc,scale,disp_length
  real(dp),dimension(-1:1,-1:1,-1:1)::cube3
  real(dp),dimension(1:nvector,1:ndim)::xtest
  real(dp),dimension(1:ndim)::gradient,displacement
  integer,dimension(1:nvector)::cell_index,cell_lev
  real(dp),dimension(1:ndim,1:ndim)::hess,minor
  real(dp),dimension(1:3)::skip_loc

#if NDIM==3

  dx=0.5D0**ilevel 
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
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
        displacement(i)=displacement(i)-minor(i,j)/det*gradient(j)
     end do
  end do

  disp_length=(displacement(1)**2+displacement(2)**2+displacement(3)**2)**0.5

  if (disp_length > 0.86*dx_loc)then
     write(*,*)'displacement to true maximum is too big and therefore shortened to dx.'
     displacement(1)=displacement(1)/disp_length*dx_loc*0.86
     displacement(2)=displacement(2)/disp_length*dx_loc*0.86
     displacement(3)=displacement(3)/disp_length*dx_loc*0.86
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
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This routine is called at the leafs of the tree structure (right after    
!update time). Here is where the global sink variables vsink and xsink are 
!updated by summing the conributions from all levels.                      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real(kind=8)::dteff
  integer::lev,isink

  if(verbose)write(*,*)'Entering update_sink for level ',ilevel

  vsold(1:nsink,1:ndim,ilevel)=vsnew(1:nsink,1:ndim,ilevel)
  vsnew(1:nsink,1:ndim,ilevel)=vsink(1:nsink,1:ndim)

  do isink=1,nsink


     ! sum force contributions from all levels
     fsink(isink,1:ndim)=0.
     do lev=levelmin,nlevelmax
        fsink(isink,1:ndim)=fsink(isink,1:ndim)+fsink_partial(isink,1:ndim,lev)
     end do
     fsink(isink,1:ndim)=fsink(isink,1:ndim)/dble(ncloud_sink)


     ! compute timestep for the synchronization
     if(level_sink(isink)>ilevel)then
        dteff=dtnew(level_sink(isink))
     else if(level_sink(isink)>0)then
        dteff=dtold(level_sink(isink))
     else
        dteff=0d0 ! timestep must be zero for newly produced sink
     end if

     ! this is the kick-kick (half old half new timestep)
     ! old timestep might be the one of a different level
     vsink(isink,1:ndim)=0.5D0*(dtnew(ilevel)+dteff)*fsink(isink,1:ndim)+vsink(isink,1:ndim)

     ! safe the velocity
     vsnew(isink,1:ndim,ilevel)=vsink(isink,1:ndim)

     ! and this is the drift(only for the global sink variable)
     xsink(isink,1:ndim)=xsink(isink,1:ndim)+vsink(isink,1:ndim)*dtnew(ilevel)
     level_sink(isink)=ilevel
  end do

end subroutine update_sink

!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine update_cloud(ilevel)
  use amr_commons
  use pm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
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
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part
  real(dp)::dx,dx_loc,scale,vol_loc
  real(dp),dimension(1:3)::skip_loc

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

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
        write(*,*),'#sink: ',isink,' drift: ',sink_jump(isink,1:ndim,ilevel)/dx_loc
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
  use hydro_commons, ONLY: uold
  implicit none
  integer::np
  integer,dimension(1:nvector)::ind_part
  !------------------------------------------------------------
  !vector loop called by update_cloud
  !------------------------------------------------------------
  integer::j,idim,isink
  real(dp),dimension(1:nvector,1:ndim),save::new_xp,new_vp
  integer,dimension(1:nvector)::level_p

  ! Overwrite cloud particle mass with sink mass
  do j=1,np
     isink=-idp(ind_part(j))
     if(isink>0)then
        mp(ind_part(j))=msink(isink)/dble(ncloud_sink)
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
           new_xp(j,idim)=new_xp(j,idim)+sink_jump(isink,idim,level_p(j))
        endif
     end do
  end do
  !write back postion
  do idim=1,ndim
     do j=1,np
        xp(ind_part(j),idim)=new_xp(j,idim)
     end do
  end do

end subroutine upd_cloud
