subroutine clump_finder
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !----------------------------------------------------------------------------
  ! Description: 
  ! Andreas Bleuler
  !----------------------------------------------------------------------------
  ! local constants
  integer::istep,ilevel,ivar,info,icpu,igrid,npartbound,isink,nmove,nmove_all,npeaks_tot,i
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m
  character(LEN=5)::nchar
  character(LEN=80)::filename

  integer,dimension(1:ncpu)::npeaks_per_cpu_tot

#if NDIM==3

  if(verbose)write(*,*)' Entering clump_finder'

  ! Create test particle
  nstar_tot=0
  do ilevel=levelmin,nlevelmax
     call create_test_particle(ilevel)
  end do
  do ilevel=nlevelmax,levelmin,-1
     call merge_tree_fine(ilevel)
  end do

  ! Move particle to the nearest peak
  nmove=nstar_tot
  istep=0
  do while(nmove>0)
     
     if(myid==1)then
        write(*,*)"istep=",istep
        write(*,*)"nmove=",nmove
     endif

     ! Move particle across oct and processor boundaries
     do ilevel=levelmin,nlevelmax
        call make_tree_fine(ilevel)
        call kill_tree_fine(ilevel)
        call virtual_tree_fine(ilevel)
     end do
     
     ! Proceed one step to densest neighbor
     nmove=0; nmove_all=0
     istep=istep+1
     do ilevel=nlevelmax,levelmin,-1
        call move_test(nmove,ilevel)
        call merge_tree_fine(ilevel)
     end do
#ifndef WITHOUTMPI     
     call MPI_ALLREDUCE(nmove,nmove_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
     nmove=nmove_all
#endif
     call title(istep,nchar)
     filename='clump/part_'//TRIM(nchar)//'.out'
     call backup_part(filename)

  end do
  print*,"escaped the loop"

  call assign_part_to_peak(npeaks_per_cpu_tot,npeaks_tot)

  do ilevel=levelmin-1,1,-1
     call merge_tree_fine(ilevel)
  end do

  ! Re-assign xp to the old position stored in vp (argument is the level in which
  ! the particles sit)
  
  
  call move_back_to_origin(1)

  do ilevel=1,nlevelmax
     call make_tree_fine(ilevel)
     call kill_tree_fine(ilevel)
     call virtual_tree_fine(ilevel)
  end do
  
  
  do ilevel=nlevelmax-1,1,-1       !remove this part and the do loop
     call merge_tree_fine(ilevel)
  end do
  
  do ilevel=1,nlevelmax
     call make_tree_fine(ilevel)
     call kill_tree_fine(ilevel)
     call virtual_tree_fine(ilevel)
  end do

  ! Get cell index, assign each cell to a global peak index
  ! Compute peak mass etc...

  call create_peak_array(npeaks_per_cpu_tot,npeaks_tot) !loop on all levels, myid only

  ! Remove cleanly test particle
  do ilevel=nlevelmax,levelmin,-1
     call merge_tree_fine(ilevel)
  end do
  call remove_test_particle(levelmin)

#endif

end subroutine clump_finder
!################################################################
!################################################################
!################################################################
!################################################################
subroutine create_test_particle(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use cooling_module, ONLY: XH=>X, rhoc, mH 
  use random
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !----------------------------------------------------------------------
  ! Description: This subroutine spawns star-particle of constant mass
  ! using a Poisson probability law if some gas condition are fulfilled. 
  ! It modifies hydrodynamic variables according to mass conservation 
  ! and assumes an isothermal transformation... 
  ! On exit, the gas velocity and sound speed are unchanged.
  ! New star particles are synchronized with other collisionless particles.
  ! Array flag2 is used as temporary work space.
  ! Yann Rasera  10/2002-01/2003
  !----------------------------------------------------------------------
  ! local constants
  real(dp)::t0,d0,e0,mgas,mcell
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:twotondim,1:3)::xc
  ! other variables
  integer ::ncache,nnew,ivar,ngrid,icpu,index_star,ndebris_tot
  integer ::igrid,ix,iy,iz,ind,i,j,n,iskip,istar,inew,nx_loc
  integer ::ntot,ntot_all,info,nstar_corrected,ideb,ndeb,nparts
#ifdef SOLVERhydro
  integer ::imetal=6
#endif
#ifdef SOLVERmhd
  integer ::imetal=9
#endif
  logical ::ok_free,ok_all
  real(dp)::d,x,y,z,u,v,w,e,zg,vdisp,dgas
  real(dp)::mstar,dstar,tstar,nISM,nCOM
  real(dp)::velc,uc,vc,wc,mass_load
  real(dp)::vxgauss,vygauss,vzgauss,birth_epoch
  real(kind=8)::mlost,mtot,mlost_all,mtot_all
  real(kind=8)::RandNum,GaussNum,PoissMean   
  real(dp)::vsn,costheta,sintheta,phi,cosphi,sinphi
  real(dp),dimension(1:3)::skip_loc
  real(dp)::dx,dx_loc,scale,vol_loc,dx_min,vol_min
  real(dp)::mdebris,vdebris,zdebris,rdebris
  real(dp)::bx1,bx2,by1,by2,bz1,bz2
  integer ,dimension(1:ncpu,1:IRandNumSize)::allseed
  integer ,dimension(1:nvector)::ind_grid,ind_cell,nstar
  integer ,dimension(1:nvector)::ind_grid_new,ind_cell_new,ind_part
  integer ,dimension(1:nvector)::list_debris,ind_debris
  logical ,dimension(1:nvector)::ok,ok_new=.true.,ok_true=.true.
  integer ,dimension(1:ncpu)::ntot_star_cpu,ntot_star_all
  
  if(numbtot(1,ilevel)==0) return
  if(.not. hydro)return
  if(ndim.ne.3)return

  if(verbose)write(*,*)' Entering star_formation'
  
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
  dx_min=(0.5D0**nlevelmax)*scale
  vol_min=dx_min**ndim

  ! Clump density threshold from H/cc to code units
  d0   = 5.d9/scale_nH

  ! Cells center position relative to grid center position
  do ind=1,twotondim  
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     xc(ind,1)=(dble(ix)-0.5D0)*dx
     xc(ind,2)=(dble(iy)-0.5D0)*dx
     xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

#if NDIM==3
  !------------------------------------------------
  ! Compute number of new stars in each cell
  !------------------------------------------------
  ntot=0 
  ! Loop over grids
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     ! Test particle formation ---> logical array ok(i)
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        ! Flag leaf cells
        do i=1,ngrid
           ok(i)=son(ind_cell(i))==0
        end do
        ! Density criterion
        do i=1,ngrid
           d=uold(ind_cell(i),1)
           if(d<=d0)ok(i)=.false. 
        end do
        ! Compute test particle map
        do i=1,ngrid
           flag2(ind_cell(i))=0
           if(ok(i))then
              flag2(ind_cell(i))=1
              ntot=ntot+1
           endif
        end do
     end do
  end do

  !---------------------------------
  ! Check for free particle memory
  !---------------------------------
  ok_free=(numbp_free-ntot)>=0
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(numbp_free,numbp_free_tot,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  numbp_free_tot=numbp_free
#endif
  if(.not. ok_free)then
     write(*,*)'No more free memory for particles'
     write(*,*)'Increase npartmax'
#ifndef WITHOUTMPI
    call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
    stop
#endif
  end if

  !---------------------------------
  ! Compute test particle statistics
  !---------------------------------
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(ntot,ntot_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  ntot_all=ntot
#endif
  ntot_star_cpu=0; ntot_star_all=0
  ntot_star_cpu(myid)=ntot
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(ntot_star_cpu,ntot_star_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  ntot_star_cpu(1)=ntot_star_all(1)
#endif
  do icpu=2,ncpu
     ntot_star_cpu(icpu)=ntot_star_cpu(icpu-1)+ntot_star_all(icpu)
  end do
  nstar_tot=nstar_tot+ntot_all
  if(myid==1)then
     if(ntot_all.gt.0)then
        write(*,'(" Level=",I6," New test particle=",I6," Tot=",I10)')&
             & ilevel,ntot_all,nstar_tot
     endif
  end if

  !------------------------------
  ! Create new test particles
  !------------------------------
  ! Starting identity number
  if(myid==1)then
     index_star=nstar_tot-ntot_all
  else
     index_star=nstar_tot-ntot_all+ntot_star_cpu(myid-1)
  end if

  ! Loop over grids
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

        ! Flag cells with test particle
        do i=1,ngrid
           ok(i)=flag2(ind_cell(i))>0
        end do

        ! Gather new test particle arrays
        nnew=0
        do i=1,ngrid
           if (ok(i))then
              nnew=nnew+1
              ind_grid_new(nnew)=ind_grid(i)
              ind_cell_new(nnew)=ind_cell(i)
           end if
        end do

        ! Update linked list for test particles
        call remove_free(ind_part,nnew)
        call add_list(ind_part,ind_grid_new,ok_new,nnew)

        ! Calculate new test particle positions
        do i=1,nnew
           index_star=index_star+1

           ! Get cell centre positions
           x=(xg(ind_grid_new(i),1)+xc(ind,1)-skip_loc(1))*scale
           y=(xg(ind_grid_new(i),2)+xc(ind,2)-skip_loc(2))*scale
           z=(xg(ind_grid_new(i),3)+xc(ind,3)-skip_loc(3))*scale

           ! Set test particle variables
           levelp(ind_part(i))=ilevel   ! Level
           idp(ind_part(i))=index_star  ! Star identity
           xp(ind_part(i),1)=x
           xp(ind_part(i),2)=y
           xp(ind_part(i),3)=z
           vp(ind_part(i),1)=x
           vp(ind_part(i),2)=y
           vp(ind_part(i),3)=z

        end do
        ! End loop over new test particles

     end do
     ! End loop over cells
  end do
  ! End loop over grids
  
#endif

end subroutine create_test_particle
!################################################################
!################################################################
!################################################################
!################################################################
subroutine remove_test_particle(ilevel)
  use amr_commons
  use pm_commons
  implicit none
  integer::ilevel
  !------------------------------------------------------------------
  !------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,idim,icpu
  integer::i,ig,ip,npart1
  real(dp)::dx

  integer,dimension(1:nvector)::ind_grid,ind_cell
  integer,dimension(1:nvector)::ind_part,ind_grid_part
  real(dp),dimension(1:nvector,1:ndim)::x0
    
  ! Loop over cpus
  do icpu=1,ncpu
     ! Loop over grids
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0   
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        if(npart1>0)then        
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           
           ! Loop over particles
           do jpart=1,npart1
              if(ig==0)then
                 ig=1
                 ind_grid(ig)=igrid
              end if
              ip=ip+1
              ind_part(ip)=ipart
              ind_grid_part(ip)=ig
              if(ip==nvector)then
                 call rm_test(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
                 ip=0
                 ig=0
              end if
              ipart=nextp(ipart)  ! Go to next particle
           end do
           ! End loop over particles
           
        end if

        igrid=next(igrid)   ! Go to next grid
     end do
     ! End loop over grids

     if(ip>0)then
        call rm_test(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
     end if

  end do
  ! End loop over cpus

end subroutine remove_test_particle
!################################################################
!################################################################
!################################################################
!################################################################
subroutine rm_test(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  logical::error
  integer::j,isink,ii,jj,kk,ind,idim,nx_loc
  real(dp)::dx_loc,scale,dx_min,xx,yy,zz,rr,r2,r2_eps
  ! Particle-based arrays
  logical,dimension(1:nvector)::ok

  do j=1,np
     ok(j)=.true.
  end do
  
  ! Remove particles from parent linked list
  call remove_list(ind_part,ind_grid_part,ok,np)
  call add_free_cond(ind_part,ok,np)

end subroutine rm_test
!################################################################
!################################################################
!################################################################
!################################################################




subroutine move_test(nmove,ilevel)
  use amr_commons
  use pm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h' 
#endif
  integer::nmove,ilevel
  !----------------------------------------------------------------------
  ! Update particle position and time-centred velocity at level ilevel. 
  ! If particle sits entirely in level ilevel, then use fine grid force
  ! for CIC interpolation. Otherwise, use coarse grid (ilevel-1) force.
  !----------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,next_part,ig,ip,npart1,info,isink
  integer,dimension(1:nvector)::ind_grid,ind_part,ind_grid_part

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

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
           next_part=nextp(ipart)
           if(ig==0)then
              ig=1
              ind_grid(ig)=igrid
           end if
           ip=ip+1
           ind_part(ip)=ipart
           ind_grid_part(ip)=ig   
           if(ip==nvector)then
              call movet(ind_grid,ind_part,ind_grid_part,ig,ip,nmove,ilevel)
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
  if(ip>0)call movet(ind_grid,ind_part,ind_grid_part,ig,ip,nmove,ilevel)

111 format('   Entering move_test for level ',I2)

end subroutine move_test
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine movet(ind_grid,ind_part,ind_grid_part,ng,np,nm,ilevel)
  use amr_commons
  use pm_commons
  use poisson_commons
  use hydro_commons, ONLY: uold
  implicit none
  integer::ng,np,nm,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !------------------------------------------------------------
  ! This routine computes the force on each particle by
  ! inverse CIC and computes new positions for all particles.
  ! If particle sits entirely in fine level, then CIC is performed
  ! at level ilevel. Otherwise, it is performed at level ilevel-1.
  ! This routine is called by move_fine.
  !------------------------------------------------------------
  logical::error
  integer::i,j,ind,idim,nx_loc,isink
  integer::i1,j1,k1,i2,j2,k2
  real(dp)::dx,length,dx_loc,scale,vol_loc,r2
  integer::i1min,i1max,j1min,j1max,k1min,k1max
  integer::i2min,i2max,j2min,j2max,k2min,k2max
  integer::i3min,i3max,j3min,j3max,k3min,k3max
  ! Grid-based arrays
  integer ,dimension(1:nvector)::father_cell
  real(dp),dimension(1:nvector,1:ndim)::x0
  integer ,dimension(1:nvector,1:threetondim)::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim)::nbors_father_grids
  ! Particle-based arrays
  logical ,dimension(1:nvector)::ok
  real(dp),dimension(1:nvector,1:ndim)::x,ff,new_xp,xtest,xmax,dd,dg
  real(dp),dimension(1:nvector)::density_max,rr
  integer ,dimension(1:nvector,1:ndim)::ig,id,igg,igd,icg,icd
  integer ,dimension(1:nvector)::cell_index,cell_levl,ind_max
  real(dp),dimension(1:nvector,1:twotondim)::vol
  real(dp),dimension(1:nvector,1:ndim,1:twotondim)::xpart
  integer ,dimension(1:nvector,1:twotondim)::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc

  ! Meshspacing in that level
  dx=0.5D0**ilevel 
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**3

  ! Integer constants
  i1min=0; i1max=0; i2min=0; i2max=0; i3min=1; i3max=1
  j1min=0; j1max=0; j2min=0; j2max=0; j3min=1; j3max=1
  k1min=0; k1max=0; k2min=0; k2max=0; k3min=1; k3max=1
  if(ndim>0)then
     i1max=2; i2max=3; i3max=2
  end if
  if(ndim>1)then
     j1max=2; j2max=3; j3max=2
  end if
  if(ndim>2)then
     k1max=2; k2max=3; k3max=2
  end if

  !====================================================
  ! Get particle density and cell
  !====================================================
  do j=1,np
     xtest(j,1:ndim)=xp(ind_part(j),1:ndim)
  end do
  call get_cell_index(cell_index,cell_levl,xtest,ilevel,np)
  do j=1,np
     density_max(j)=uold(cell_index(j),1)
     ind_max(j)=cell_index(j)
     xmax(j,1:ndim)=xp(ind_part(j),1:ndim)
  end do

  !====================================================
  ! Check for potential new positions at level ilevel-1
  !====================================================
  if(ilevel>levelmin)then

  ! Lower left corner of 3x3x3 grid-cube
  do idim=1,ndim
     do i=1,ng
        x0(i,idim)=xg(ind_grid(i),idim)-3.0D0*dx
     end do
  end do
  
  ! Rescale particle position at level ilevel-1
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
  do idim=1,ndim
     do j=1,np
        x(j,idim)=x(j,idim)/2.0D0
     end do
  end do

  ! Check for illegal moves
  error=.false.
  do idim=1,ndim
     do j=1,np
        if(x(j,idim)<0.5D0.or.x(j,idim)>2.5D0)error=.true.
     end do
  end do
  if(error)then
     write(*,*)'problem in move'
     do idim=1,ndim
        do j=1,np
           if(x(j,idim)<0.5D0.or.x(j,idim)>2.5D0)then
              write(*,*)x(j,1:ndim)
           endif
        end do
     end do
     stop
  end if

  ! Do CIC at level ilevel-1
  do idim=1,ndim
     do j=1,np
        dd(j,idim)=x(j,idim)+0.5D0
        id(j,idim)=dd(j,idim)
        dd(j,idim)=dd(j,idim)-id(j,idim)
        dg(j,idim)=1.0D0-dd(j,idim)
        ig(j,idim)=id(j,idim)-1
     end do
  end do

 ! Compute parent cell position
#if NDIM==1
  do j=1,np
     xpart(j,1,1)=0.5+ig(j,1)
     xpart(j,1,2)=0.5+id(j,1)
  end do
#endif
#if NDIM==2
  do j=1,np
     ! Particle 1
     xpart(j,1,1)=0.5+ig(j,1)
     xpart(j,2,1)=0.5+ig(j,2)
     ! Particle 2
     xpart(j,1,2)=0.5+id(j,1)
     xpart(j,2,2)=0.5+ig(j,2)
     ! Particle 3
     xpart(j,1,3)=0.5+ig(j,1)
     xpart(j,2,3)=0.5+id(j,2)
     ! Particle 4
     xpart(j,1,4)=0.5+id(j,1)
     xpart(j,2,4)=0.5+id(j,2)
  end do
#endif
#if NDIM==3
  do j=1,np
     ! Particle 1
     xpart(j,1,1)=0.5+ig(j,1)
     xpart(j,2,1)=0.5+ig(j,2)
     xpart(j,3,1)=0.5+ig(j,3)
     ! Particle 2
     xpart(j,1,2)=0.5+id(j,1)
     xpart(j,2,2)=0.5+ig(j,2)
     xpart(j,3,2)=0.5+ig(j,3)
     ! Particle 3
     xpart(j,1,3)=0.5+ig(j,1)
     xpart(j,2,3)=0.5+id(j,2)
     xpart(j,3,3)=0.5+ig(j,3)
     ! Particle 4
     xpart(j,1,4)=0.5+id(j,1)
     xpart(j,2,4)=0.5+id(j,2)
     xpart(j,3,4)=0.5+ig(j,3)
     ! Particle 5
     xpart(j,1,5)=0.5+ig(j,1)
     xpart(j,2,5)=0.5+ig(j,2)
     xpart(j,3,5)=0.5+id(j,3)
     ! Particle 6
     xpart(j,1,6)=0.5+id(j,1)
     xpart(j,2,6)=0.5+ig(j,2)
     xpart(j,3,6)=0.5+id(j,3)
     ! Particle 7
     xpart(j,1,7)=0.5+ig(j,1)
     xpart(j,2,7)=0.5+id(j,2)
     xpart(j,3,7)=0.5+id(j,3)
     ! Particle 8
     xpart(j,1,8)=0.5+id(j,1)
     xpart(j,2,8)=0.5+id(j,2)
     xpart(j,3,8)=0.5+id(j,3)
  end do
#endif
        
  ! Test those particles
  do ind=1,twotondim
     do idim=1,ndim
        do j=1,np
           xtest(j,idim)=xpart(j,idim,ind)*2.*dx+x0(ind_grid_part(j),idim)
        end do
        do j=1,np
           xtest(j,idim)=(xtest(j,idim)-skip_loc(idim))*scale
        end do
     end do
     call get_cell_index(cell_index,cell_levl,xtest,ilevel-1,np)
     do j=1,np
        if(son(cell_index(j))==0)then
           if(uold(cell_index(j),1)>density_max(j))then
              density_max(j)=uold(cell_index(j),1)
              ind_max(j)=cell_index(j)
              xmax(j,1:ndim)=xtest(j,1:ndim)
           endif
        endif
     end do
  end do

  endif

  !====================================================
  ! Check for potential new positions at level ilevel
  !====================================================
  ! Generate 3x3x3 neighboring cells at level ilevel
  do k1=k1min,k1max
  do j1=j1min,j1max
  do i1=i1min,i1max

     do j=1,np
        xtest(j,1)=xp(ind_part(j),1)+(i1-1)*dx_loc
#if NDIM>1
        xtest(j,2)=xp(ind_part(j),2)+(j1-1)*dx_loc
#endif     
#if NDIM>2
        xtest(j,3)=xp(ind_part(j),3)+(k1-1)*dx_loc
#endif     
     end do

     call get_cell_index(cell_index,cell_levl,xtest,ilevel,np)

     do j=1,np
        if(son(cell_index(j))==0.and.cell_levl(j)==ilevel)then
           if(uold(cell_index(j),1)>density_max(j))then
              !if(myid==1.and.ind_part(j)==959)then
                 !write(*,*)'Particle ',j,ind_part(j)
                 !write(*,*)uold(cell_index(j),1),density_max(j)
                 !write(*,*)xp(ind_part(j),1:ndim)
                 !write(*,*)xtest(j,1:ndim)
              !endif
              density_max(j)=uold(cell_index(j),1)
              ind_max(j)=cell_index(j)
              xmax(j,1:ndim)=xtest(j,1:ndim)
           endif
        endif
     end do

  end do
  end do
  end do

  !====================================================
  ! Check for potential new positions at level ilevel+1
  !====================================================
  if(ilevel<nlevelmax)then

  ! Generate 4x4x4 neighboring cells at level ilevel+1
  do k2=k2min,k2max
  do j2=j2min,j2max
  do i2=i2min,i2max

     do j=1,np
        xtest(j,1)=xp(ind_part(j),1)+(i2-1.5)*dx_loc/2.0
#if NDIM>1
        xtest(j,2)=xp(ind_part(j),2)+(j2-1.5)*dx_loc/2.0
#endif     
#if NDIM>2
        xtest(j,3)=xp(ind_part(j),3)+(k2-1.5)*dx_loc/2.0
#endif     
     end do
     call get_cell_index(cell_index,cell_levl,xtest,ilevel+1,np)
     do j=1,np
        if(son(cell_index(j))==0.and.cell_levl(j)==(ilevel+1))then
           if(uold(cell_index(j),1)>density_max(j))then
              density_max(j)=uold(cell_index(j),1)
              ind_max(j)=cell_index(j)
              xmax(j,1:ndim)=xtest(j,1:ndim)
           endif
        endif
     end do
  end do
  end do
  end do

  endif

  !====================================================
  ! Update position
  !====================================================
  rr(1:np)=0.0
  do idim=1,ndim
     do j=1,np
        rr(j)=rr(j)+(xp(ind_part(j),idim)-xmax(j,idim))**2
     end do
  end do
  do j=1,np
     xp(ind_part(j),1:ndim)=xmax(j,1:ndim)
     idp(ind_part(j))=ind_max(j)
  end do
  do j=1,np
     if(rr(j)>1d-3*dx_loc**2)then
        nm=nm+1
     endif
  end do

end subroutine movet

!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine assign_part_to_peak(npeaks_per_cpu_tot,npeaks_tot)
  use amr_commons
  use hydro_commons
  use pm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::igrid,jgrid,ipart,jpart,next_part,ig,ip,ilevel,npart1,info,nparts,nparts_tot
  integer::npeaks,npeaks_tot,jj,peak_nr,icpu
  integer*8,dimension(1)::n_cls

  integer,dimension(1:nvector)::ind_grid,ind_cell,init_ind_cell,init_cell_lev,cell_lev
  integer,dimension(1:nvector)::ind_part,ind_grid_part
  real(dp),dimension(1:nvector,1:ndim)::pos

  integer,dimension(1:ncpu)::npeaks_per_cpu,npeaks_per_cpu_tot



  nparts=0
  flag2=0   !use flag2 as temporary array 

  !loop over all particles (on levelmin) and write -1 into the flag2 of each cell which 
  !contains a particle
  ilevel=levelmin
  ig=0
  ip=0
  ! Loop over grids 
  do icpu=1,ncpu !loop cpus
     igrid=headl(icpu,ilevel) 
     do jgrid=1,numbl(icpu,ilevel) ! Number of grids in the level ilevel on process myid?
        npart1=numbp(igrid)  ! Number of particles in the grid
        nparts=npart1+nparts
        if(npart1>0)then
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle  <---- Very important !!!
              next_part=nextp(ipart)
              if (xp(ipart,1)>0.75.or.xp(ipart,1)<0.25)print*,'alert_x',ipart,xp(ipart,1)
              if (xp(ipart,2)>0.75.or.xp(ipart,2)<0.25)print*,'alert_y',ipart,xp(ipart,2)
              if (xp(ipart,3)>0.75.or.xp(ipart,3)<0.25)print*,'alert_z',ipart,xp(ipart,3)
              if (vp(ipart,1)>0.75.or.vp(ipart,1)<0.25)print*,'alert_vx',ipart,vp(ipart,1)
              if (vp(ipart,2)>0.75.or.vp(ipart,2)<0.25)print*,'alert_vy',ipart,vp(ipart,2)
              if (vp(ipart,3)>0.75.or.vp(ipart,3)<0.25)print*,'alert_vz',ipart,vp(ipart,3)
              if(ig==0)then
                 ig=1
                 ind_grid(ig)=igrid
              end if
              ip=ip+1
              ind_part(ip)=ipart
              ind_grid_part(ip)=ig   
              if(ip==nvector)then 
                 call get_cell_indices(init_ind_cell,init_cell_lev,ind_cell,cell_lev,ind_part,pos,ip,nlevelmax)
                 do jj=1,nvector
                    flag2(ind_cell(jj))=-1
                 end do
                 ip=0
                 ig=0
              end if
              ipart=next_part  ! Go to next particle
           end do
           ! End loop over particles
        end if
        igrid=next(igrid)   ! Go to next grid
     end do
     if(ip>0)then 
        call get_cell_indices(init_ind_cell,init_cell_lev,ind_cell,cell_lev,ind_part,pos,ip,nlevelmax)
        do jj=1,ip
           flag2(ind_cell(jj))=-1
        end do
     end if
  end do

#ifndef WITHOUTMPI     
  call MPI_ALLREDUCE(nparts,nparts_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#endif


  !end loop over all particles
  print*,'parts on myid ',myid,' = ',nparts,nparts_tot


  !calculate peaks per cpu array
  npeaks=0  
  n_cls=shape(flag2)
  do jj=0,n_cls(1)
     npeaks=npeaks-flag2(jj)
  end do
  print*,'n_peaks on processor number',myid,'= ',npeaks

  npeaks_per_cpu=0; npeaks_per_cpu_tot=0
  npeaks_per_cpu(myid)=npeaks

#ifndef WITHOUTMPI     
  call MPI_ALLREDUCE(npeaks_per_cpu,npeaks_per_cpu_tot,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI     
  npeaks_per_cpu_tot=npeaks_per_cpu
#endif
#ifndef WITHOUTMPI     
  call MPI_ALLREDUCE(npeaks,npeaks_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI     
  npeaks_tot=npeaks
#endif

  !determine writing positions for cpu
  peak_nr=1
  do jj=1,myid-1
     peak_nr=peak_nr+npeaks_per_cpu_tot(jj)
  end do

  !write the peak_number into each cell above the threshold
  n_cls=shape(flag2)
  do jj=0,n_cls(1)
     if(flag2(jj)==-1)then
        flag2(jj)=peak_nr
        peak_nr=peak_nr+1
     end if
  end do

  !loop over all particles (on levelmin) and write the peak_nr into the mass variable of each   
  !particle
  ilevel=levelmin
  ig=0
  ip=0
  do icpu=1,ncpu !loop cpus
     ! Loop over grids
     igrid=headl(icpu,ilevel) 
     do jgrid=1,numbl(icpu,ilevel) ! Number of grids in the level ilevel on process myid?
        npart1=numbp(igrid)  ! Number of particles in the grid
        if(npart1>0)then
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle  <---- Very important !!!
              next_part=nextp(ipart)
              if(ig==0)then
                 ig=1
                 ind_grid(ig)=igrid
              end if
              ip=ip+1
              ind_part(ip)=ipart
              ind_grid_part(ip)=ig   
              if(ip==nvector)then 
                 call get_cell_indices(init_ind_cell,init_cell_lev,ind_cell,cell_lev,ind_part,pos,ip,nlevelmax)
                 do jj=1,nvector
                    mp(ind_part(jj))=1.*flag2(ind_cell(jj))
                 end do
                 ip=0
                 ig=0
              end if
              ipart=next_part  ! Go to next particle
           end do
           ! End loop over particles
        end if
        igrid=next(igrid)   ! Go to next grid
     end do
     if(ip>0)then 
        call get_cell_indices(init_ind_cell,init_cell_lev,ind_cell,cell_lev,ind_part,pos,ip,nlevelmax)
        do jj=1,ip
           mp(ind_part(jj))=1.*flag2(ind_cell(jj))  
        end do
     end if
  end do
  !end loop over all particles



end subroutine assign_part_to_peak


!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine create_peak_array(npeaks_per_cpu_tot,npeaks_tot)
  use amr_commons
  use hydro_commons
  use pm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::igrid,jgrid,ipart,jpart,next_part,ig,ip,npart1,info,ilevel,nx_loc
  integer::nparts,nparts_tot,npeaks,npeaks_tot,jj,kk,peak_map_pos,i,peak_nr
  integer*8,dimension(1)::n_cls

  integer,dimension(1:nvector)::ind_grid,ind_cell,init_ind_cell,init_cell_lev,cell_lev
  integer,dimension(1:nvector)::ind_part,ind_grid_part
  real(dp),dimension(1:nvector,1:ndim)::pos
  real(dp),dimension(1:ndim)::position


  integer,dimension(1:ncpu)::npeaks_per_cpu,npeaks_per_cpu_tot

  !allocatable arrays with the size of number of particles on myid, one for the particle positions
  !and one for the map linking the initial cell index and level to the peak number
  integer,allocatable,dimension(:,:)::peak_map
  real(dp),allocatable,dimension(:,:)::part_pos


  !peak_properties
  real(dp),dimension(1:npeaks_tot,1:ndim)::peak_pos,peak_pos_tot,clump_size,clump_size_tot
  real(dp),dimension(1:npeaks_tot)::min_dens,min_dens_tot,max_dens,max_dens_tot
  real(dp),dimension(1:npeaks_tot)::clump_mass,clump_mass_tot,clump_vol,clump_vol_tot
  real(dp)::tot_mass

  character(LEN=5)::myidstring



  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::dx,dx_loc,scale,vol_loc,dx_min,vol_min
  real(dp),dimension(1:nlevelmax)::volume
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  
#if NDIM==3

  do ilevel=1,nlevelmax
     ! Mesh spacing in that level
     dx=0.5D0**ilevel 
     nx_loc=(icoarse_max-icoarse_min+1)
     scale=boxlen/dble(nx_loc)
     dx_loc=dx*scale
     vol_loc=dx_loc**ndim
     volume(ilevel)=vol_loc
  end do

  !find number of particles on cpu
  nparts=0
  !loop over all particles
  do ilevel=nlevelmax,levelmin,-1
     ! Loop over grids
     igrid=headl(myid,ilevel) 
     do jgrid=1,numbl(myid,ilevel) 
        nparts=nparts+numbp(igrid)
        igrid=next(igrid)   ! Go to next grid
     end do
  end do
  !end loop over all particles


#ifndef WITHOUTMPI     
  call MPI_ALLREDUCE(nparts,nparts_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI     
  nparts_tot=nparts
#endif
  print*,'number of particles on ',myid,' after moving back = ',nparts,' of ',nparts_tot


  allocate(part_pos(nparts,ndim))
  allocate(peak_map(nparts,3))


  peak_map=0; part_pos=0.

  peak_map_pos=1
  !loop over all particles
  do ilevel=nlevelmax,levelmin,-1
     ig=0
     ip=0
     ! Loop over grids
     igrid=headl(myid,ilevel) 
     do jgrid=1,numbl(myid,ilevel) ! Number of grids in the level ilevel on process myid?
        npart1=numbp(igrid)  ! Number of particles in the grid
        if(npart1>0)then
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle  <---- Very important !!!
              if (xp(ipart,1)>0.75.or.xp(ipart,1)<0.25)print*,'alert_x',ipart,xp(ipart,1)
              if (xp(ipart,2)>0.75.or.xp(ipart,2)<0.25)print*,'alert_y',ipart,xp(ipart,2)
              if (xp(ipart,3)>0.75.or.xp(ipart,3)<0.25)print*,'alert_z',ipart,xp(ipart,3)
              if (vp(ipart,1)>0.75.or.vp(ipart,1)<0.25)print*,'alert_vx',ipart,vp(ipart,1)
              if (vp(ipart,2)>0.75.or.vp(ipart,2)<0.25)print*,'alert_vy',ipart,vp(ipart,2)
              if (vp(ipart,3)>0.75.or.vp(ipart,3)<0.25)print*,'alert_vz',ipart,vp(ipart,3)
              next_part=nextp(ipart)
              if(ig==0)then
                 ig=1
                 ind_grid(ig)=igrid
              end if
              ip=ip+1
              ind_part(ip)=ipart
              ind_grid_part(ip)=ig   
              if(ip==nvector)then 
                 call get_cell_indices(init_ind_cell,init_cell_lev,ind_cell,cell_lev,ind_part,pos,ip,nlevelmax)
                 do jj=1,nvector
                    peak_map(peak_map_pos,1)=ind_cell(jj)
                    peak_map(peak_map_pos,2)=cell_lev(jj)
                    peak_map(peak_map_pos,3)=int(mp(ind_part(jj)))
                    do i=1,ndim
                       part_pos(peak_map_pos,i)=pos(jj,i)
                    end do
                    peak_map_pos=peak_map_pos+1
                 end do
                 ip=0
                 ig=0
              end if
              ipart=next_part  ! Go to next particle
           end do
           ! End loop over particles
        end if
        igrid=next(igrid)   ! Go to next grid
     end do
     if(ip>0)then 
        call get_cell_indices(init_ind_cell,init_cell_lev,ind_cell,cell_lev,ind_part,pos,ip,nlevelmax)
        do jj=1,ip
           peak_map(peak_map_pos,1)=ind_cell(jj)
           peak_map(peak_map_pos,2)=cell_lev(jj)
           peak_map(peak_map_pos,3)=int(mp(ind_part(jj)))
           do i=1,ndim
              part_pos(peak_map_pos,i)=pos(jj,i)
           end do
           peak_map_pos=peak_map_pos+1
        end do
     end if
  end do
  !end loop over all particles



  !write peak_map to file
  call title(myid,myidstring)
  open(unit=20,file=TRIM('clump/peak_map.txt'//myidstring),form='formatted')
  write(20,*)nparts
  do kk=1,ndim
     do jj=1,nparts
        write(20,'(F8.6)')part_pos(jj,kk)
     end do
  end do
  do jj=1,nparts
     write(20,*)peak_map(jj,3)
  end do
  close(20)




  !allocate all peak based arrays
  peak_pos=0.; peak_pos_tot=0.
  clump_size=0.; clump_size_tot=0.
  min_dens=1.d99; min_dens_tot=1.d99;
  max_dens=0.; max_dens_tot=0.
  clump_mass=0.;clump_mass_tot=0.
  clump_vol=0.;clump_vol_tot=0.

  n_cls=shape(flag2)
  do jj=0,n_cls(1)
     if(flag2(jj)>0)then
        call get_cell_center(jj,position,ilevel)
        peak_pos(flag2(jj),1:ndim)=position(1:ndim)
     end if
  end do

  !collect information related to the clumps
  do jj=1,nparts
     peak_nr=peak_map(jj,3) !find peak_number
     !find min density
     if(uold(peak_map(jj,1),1)<=min_dens(peak_nr))then 
        min_dens(peak_nr)=uold(peak_map(jj,1),1)
     end if    
     !find max density
     if(uold(peak_map(jj,1),1)>=max_dens(peak_nr))then 
        max_dens(peak_nr)=uold(peak_map(jj,1),1)
     end if
     !find clump mass
     clump_mass(peak_nr)=clump_mass(peak_nr)+ &
          volume(peak_map(jj,2))*uold(peak_map(jj,1),1)
     !find clump size (take center of mass instead of peak as reference point)
     do i=1,ndim
        clump_size(peak_nr,i)=clump_size(peak_nr,i)+ & 
             (part_pos(jj,i)-peak_pos(peak_nr,i))**2 * volume(peak_map(jj,2))
     end do
     !clump volume
     clump_vol(peak_nr)=clump_vol(peak_nr)+volume(peak_map(jj,2))
  end do

#ifndef WITHOUTMPI     
  call MPI_ALLREDUCE(min_dens,min_dens_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI     
  min_dens_tot=min_dens
#endif
#ifndef WITHOUTMPI     
  call MPI_ALLREDUCE(max_dens,max_dens_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI     
  max_dens_tot=max_dens
#endif
#ifndef WITHOUTMPI     
  call MPI_ALLREDUCE(clump_mass,clump_mass_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI     
  clump_mass_tot=clump_mass
#endif
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(clump_vol,clump_vol_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,\
info)
#endif
#ifdef WITHOUTMPI
  clump_vol_tot=clump_vol
#endif
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(clump_size,clump_size_tot,3*npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,\
info)
#endif
#ifdef WITHOUTMPI
  clump_size_tot=clump_size
#endif
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(peak_pos,peak_pos_tot,3*npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,\
info)
#endif
#ifdef WITHOUTMPI
  peak_pos_tot=peak_pos
#endif


  if(myid==1) then 
     tot_mass=0.
     print*,'Cl_N  peak_x [uu] peak_y [uu] peak_z [uu] size_x [AU]'&
          ,' size_y [AU] size_z [AU] rho- [g/cc] rho+ [g/cc] rho_av [g/cc] M_cl [M_sol] V_cl [AU^3]'
     do jj=1,npeaks_tot
        write(*,'(I5,3(X,F11.5),3(X,F11.2),3(XE11.2E2),X,F13.5,XE11.2E2)'),jj&
             ,peak_pos_tot(jj,1),peak_pos_tot(jj,2),peak_pos_tot(jj,3)&
             ,(5.*clump_size_tot(jj,1)/clump_vol_tot(jj))**0.5*scale_l/1.496d13&
             ,(5.*clump_size_tot(jj,2)/clump_vol_tot(jj))**0.5*scale_l/1.496d13&
             ,(5.*clump_size_tot(jj,3)/clump_vol_tot(jj))**0.5*scale_l/1.496d13&
             ,min_dens_tot(jj)*scale_d,max_dens_tot(jj)*scale_d&
             ,clump_mass_tot(jj)/clump_vol_tot(jj)*scale_d&
             ,clump_mass_tot(jj)*scale_d*scale_l**3/1.98d33&
             ,clump_vol_tot(jj)*(scale_l/1.496d13)**3
        tot_mass=tot_mass+clump_mass_tot(jj)*scale_d*scale_l**3/1.98d33
     end do
     print*,'total mass in clumps =',tot_mass
  end if


  
#ifndef WITHOUTMPI     
  call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif

#endif

end subroutine create_peak_array





!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine get_cell_indices(init_cell_index,init_cell_lev,cell_index,cell_lev,ind_part,init_pos,np,ilevel)
  use amr_commons
  use hydro_commons
  use pm_commons
  implicit none
  real(dp),dimension(1:nvector,1:ndim)::init_pos
  real(dp),dimension(1:nvector,1:ndim)::xtest
  integer::ilevel,j,np
  integer,dimension(1:nvector)::ind_part,init_cell_index,cell_index,cell_lev,init_cell_lev

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!gets the index of the initial cell and the index of the final cell (peak!
!cell) and sets the temporary array flag2 1 at each peak                 !   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do j=1,np
     init_pos(j,1)=vp(ind_part(j),1)
     xtest(j,1)=xp(ind_part(j),1)
#if NDIM>1
     init_pos(j,2)=vp(ind_part(j),2)
     xtest(j,2)=xp(ind_part(j),2)
#endif     
#if NDIM>2
     init_pos(j,3)=vp(ind_part(j),3)  
     xtest(j,3)=xp(ind_part(j),3)
#endif     
  end do

  call get_cell_index(cell_index,cell_lev,xtest,ilevel,np)
  call get_cell_index(init_cell_index,init_cell_lev,init_pos,ilevel,np)

  
end subroutine get_cell_indices
  
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine get_cell_index(cell_index,cell_levl,xpart,ilevel,np)
  use amr_commons
  implicit none
  integer::np,ilevel
  integer,dimension(1:nvector)::cell_index,cell_levl
  real(dp),dimension(1:nvector,1:3)::xpart
  ! This function returns the index of the cell, at maximum level
  ! ilevel, in which the input particle sits
  real(dp)::xx,yy,zz
  integer::i,j,ii,jj,kk,ind,iskip,igrid,ind_cell,igrid0

  if ((nx.eq.1).and.(ny.eq.1).and.(nz.eq.1)) then
  else if ((nx.eq.3).and.(ny.eq.3).and.(nz.eq.3)) then
  else
     write(*,*)"nx=ny=nz != 1,3 is not supported."
     call clean_stop
  end if

  igrid0=son(1+icoarse_min+jcoarse_min*nx+kcoarse_min*nx*ny)
  do i=1,np
     xx = xpart(i,1) + (nx-1)/2.0
     yy = xpart(i,2) + (ny-1)/2.0
     zz = xpart(i,3) + (nz-1)/2.0
     igrid=igrid0
     do j=1,ilevel !!!!!!!!!!!!!!!!changed here 1 to levelmin (doesn't work)
        ii=1; jj=1; kk=1
        if(xx<xg(igrid,1))ii=0
        if(yy<xg(igrid,2))jj=0
        if(zz<xg(igrid,3))kk=0
        ind=1+ii+2*jj+4*kk
        iskip=ncoarse+(ind-1)*ngridmax
        ind_cell=iskip+igrid
        igrid=son(ind_cell)
        if(igrid==0.or.j==ilevel)exit
     end do
     cell_index(i)=ind_cell
     cell_levl(i)=j
  end do
end subroutine get_cell_index



!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine get_cell_center(ind_cell,cell_pos,ilevel)
  use amr_commons
  use hydro_commons
  use pm_commons
  implicit none
  real(dp),dimension(1:ndim)::cell_pos
  integer::ilevel,jj,nc
  integer::ind_cell,ind_ind,ind_grid


  real(dp)::t0,d0,e0,mgas,mcell
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:twotondim,1:3)::xc
  integer ::igrid,ix,iy,iz,ind,i,j,n,iskip,istar,inew,nx_loc
  real(dp)::dx,dx_loc,scale,vol_loc,dx_min,vol_min
  real(dp),dimension(1:3)::skip_loc
  logical::leafe

#if NDIM==3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!gets the position of the cell_center and the level of the leave-cell as !
!a function of cell_id                                                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  leafe=.false.
  ilevel=levelmin
  do while(leafe .EQV. .false.)
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
     
     
     ! Cells center position relative to grid center position
     do ind=1,twotondim  
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        xc(ind,1)=(dble(ix)-0.5D0)*dx
        xc(ind,2)=(dble(iy)-0.5D0)*dx
        xc(ind,3)=(dble(iz)-0.5D0)*dx
     end do
     
     ind_grid=mod((ind_cell-ncoarse),ngridmax)
     ind_ind=1+(ind_cell-ncoarse)/ngridmax
     cell_pos(1)=(xg(ind_grid,1)+xc(ind,1)-skip_loc(1))*scale
     cell_pos(2)=(xg(ind_grid,2)+xc(ind,2)-skip_loc(2))*scale
     cell_pos(3)=(xg(ind_grid,3)+xc(ind,3)-skip_loc(3))*scale
     
     !check weather cell is refined
     if (son(ind_cell)==0)leafe=.true. 
     ilevel=ilevel+1
     
  end do

#endif
end subroutine get_cell_center






!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine move_back_to_origin(ilevel)
  use amr_commons
  use hydro_commons
  use pm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::igrid,jgrid,ipart,jpart,next_part,ilevel,npart1,info,icpu,nparts

  nparts=0

  !loop over all particles
  do icpu=1,ncpu
     ! Loop over grids
     igrid=headl(icpu,ilevel) 
     do jgrid=1,numbl(icpu,ilevel) ! Number of grids 
        npart1=numbp(igrid)  ! Number of particles in the grid
        nparts=nparts+npart1
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle  <---- Very important !!!
              next_part=nextp(ipart)
              xp(ipart,1:ndim)=vp(ipart,1:ndim)!
              ipart=next_part  ! Go to next particle
           end do
           ! End loop over particles
        end if
        igrid=next(igrid)   ! Go to next grid
     end do
  end do
  !end loop over all particles
  
  print*,'number of particles moved on ',myid,' = ',nparts


end subroutine move_back_to_origin
