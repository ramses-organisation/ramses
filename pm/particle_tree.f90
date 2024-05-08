!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_tree
  use pm_commons
  use amr_commons
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::info
#endif
  !------------------------------------------------------
  ! This subroutine build the particle linked list at the
  ! coarse level for ALL the particles in the box.
  ! This routine should be used only as initial set up for
  ! the particle tree.
  !------------------------------------------------------
  integer::ipart,idim,i,nxny,ilevel
  integer::npart1,icpu,nx_loc
  logical::error
  real(dp),dimension(1:3)::xbound
  integer,dimension(1:nvector),save::ix,iy,iz
  integer,dimension(1:nvector),save::ind_grid,ind_part
  logical,dimension(1:nvector),save::ok=.true.
  real(dp),dimension(1:3)::skip_loc
  real(dp)::scale

  if(verbose)write(*,*)'  Entering init_tree'

  ! Local constants
  nxny=nx*ny
  xbound(1:3)=(/dble(nx),dble(ny),dble(nz)/)
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)

  !----------------------------------
  ! Initialize particle linked list
  !----------------------------------
  prevp(1)=0; nextp(1)=2
  do ipart=2,npartmax-1
     prevp(ipart)=ipart-1
     nextp(ipart)=ipart+1
  end do
  prevp(npartmax)=npartmax-1; nextp(npartmax)=0
  ! Free memory linked list
  headp_free=npart+1
  tailp_free=npartmax
  numbp_free=tailp_free-headp_free+1
  if(numbp_free>0)then
     prevp(headp_free)=0
  end if
  nextp(tailp_free)=0
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(numbp_free,numbp_free_tot,1,MPI_INTEGER,MPI_MIN,&
       & MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  numbp_free_tot=numbp_free
#endif

  !--------------
  ! Periodic box
  !--------------
  do idim=1,ndim
     do ipart=1,npart
        if(xp(ipart,idim)/scale+skip_loc(idim)<0.0d0) &
             & xp(ipart,idim)=xp(ipart,idim)+(xbound(idim)-skip_loc(idim))*scale
        if(xp(ipart,idim)/scale+skip_loc(idim)>=xbound(idim)) &
             & xp(ipart,idim)=xp(ipart,idim)-(xbound(idim)-skip_loc(idim))*scale
     end do
     if(sink)then
        do ipart=1,nsink
           if(xsink(ipart,idim)/scale+skip_loc(idim)<0.0d0) &
                & xsink(ipart,idim)=xsink(ipart,idim)+(xbound(idim)-skip_loc(idim))*scale
           if(xsink(ipart,idim)/scale+skip_loc(idim)>=xbound(idim)) &
                & xsink(ipart,idim)=xsink(ipart,idim)-(xbound(idim)-skip_loc(idim))*scale
        end do
     endif
  end do

  !----------------------------------
  ! Reset all linked lists at level 1
  !----------------------------------
  do i=1,active(1)%ngrid
     headp(active(1)%igrid(i))=0
     tailp(active(1)%igrid(i))=0
     numbp(active(1)%igrid(i))=0
  end do
  do icpu=1,ncpu
     do i=1,reception(icpu,1)%ngrid
#ifdef LIGHT_MPI_COMM
        headp(reception(icpu,1)%pcomm%igrid(i))=0
        tailp(reception(icpu,1)%pcomm%igrid(i))=0
        numbp(reception(icpu,1)%pcomm%igrid(i))=0
#else
        headp(reception(icpu,1)%igrid(i))=0
        tailp(reception(icpu,1)%igrid(i))=0
        numbp(reception(icpu,1)%igrid(i))=0
#endif
     end do
  end do

  !------------------------------------------------
  ! Build linked list at level 1 by vector sweeps
  !------------------------------------------------
  do ipart=1,npart,nvector
     npart1=min(nvector,npart-ipart+1)
     ! Gather particles
     do i=1,npart1
        ind_part(i)=ipart+i-1
     end do
     ! Compute coarse cell
#if NDIM>0
     do i=1,npart1
        ix(i)=int(xp(ind_part(i),1)/scale+skip_loc(1))
     end do
#endif
#if NDIM>1
     do i=1,npart1
        iy(i)=int(xp(ind_part(i),2)/scale+skip_loc(2))
     end do
#endif
#if NDIM>2
     do i=1,npart1
        iz(i)=int(xp(ind_part(i),3)/scale+skip_loc(3))
     end do
#endif
     ! Compute level 1 grid index
     error=.false.
     do i=1,npart1
        ind_grid(i)=son(1+ix(i)+nx*iy(i)+nxny*iz(i))
        if(ind_grid(i)==0)error=.true.
     end do
     if(error)then
        write(*,*)'Error in init_tree'
        write(*,*)'Particles appear in unrefined regions'
        call clean_stop
     end if
     ! Add particle to level 1 linked list
     call add_list(ind_part,ind_grid,ok,npart1)
  end do

  ! destroy and recreate cloud particles to account for changes in sink
  ! radius, newly added sinks, etc
  do ilevel=levelmin-1,1,-1
     call merge_tree_fine(ilevel)
  end do

#if NDIM==3
  if(sink)then
     call kill_entire_cloud(1)
     call create_cloud_from_sink
  endif
#endif

  ! Sort particles down to levelmin
  do ilevel=1,levelmin-1
     call make_tree_fine(ilevel)
     call kill_tree_fine(ilevel)
     ! Update boundary conditions for remaining particles
     call virtual_tree_fine(ilevel)
  end do

end subroutine init_tree
!################################################################
!################################################################
!################################################################
!################################################################
subroutine make_tree_fine(ilevel)
  use pm_commons
  use amr_commons
  implicit none
  integer::ilevel
  !-----------------------------------------------------------------------
  ! This subroutine checks if particles have moved from their parent grid
  ! to one of the 3**ndim neighboring sister grids. The particle is then
  ! disconnected from the parent grid linked list, and connected to the
  ! corresponding sister grid linked list. If the sister grid does
  ! not exist, the particle is left to its original parent grid.
  ! Particles must not move to a distance greater than direct neighbors
  ! boundaries. Otherwise an error message is issued and the code stops.
  !-----------------------------------------------------------------------
  integer::idim,nx_loc
  real(dp)::dx,scale
  real(dp),dimension(1:3)::xbound
  real(dp),dimension(1:3)::skip_loc
  integer::igrid,jgrid,ipart,jpart,next_part
  integer::ig,ip,npart1,icpu
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  xbound(1:3)=(/dble(nx),dble(ny),dble(nz)/)
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)

  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        if(npart1>0)then
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle  <--- Very important !!!
              next_part=nextp(ipart)
              if(ig==0)then
                 ig=1
                 ind_grid(ig)=igrid
              end if
              ip=ip+1
              ind_part(ip)=ipart
              ind_grid_part(ip)=ig
              ! Gather nvector particles
              if(ip==nvector)then
                 call check_tree(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
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
     if(ip>0)call check_tree(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
  end do
  ! End loop over cpus

  ! Periodic boundaries
  if(sink)then
     do idim=1,ndim
        do ipart=1,nsink
           if(xsink(ipart,idim)/scale+skip_loc(idim)<0.0d0) &
                & xsink(ipart,idim)=xsink(ipart,idim)+(xbound(idim)-skip_loc(idim))*scale
           if(xsink(ipart,idim)/scale+skip_loc(idim)>=xbound(idim)) &
                & xsink(ipart,idim)=xsink(ipart,idim)-(xbound(idim)-skip_loc(idim))*scale
        end do
     end do
  endif

111 format('   Entering make_tree_fine for level ',I2)

end subroutine make_tree_fine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine check_tree(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use pm_commons
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !-----------------------------------------------------------------------
  ! This routine is called by make_tree_fine.
  !-----------------------------------------------------------------------
  logical::error
  integer::i,j,idim,nx_loc
  real(dp)::dx,xxx,scale
  real(dp),dimension(1:3)::xbound
  ! Grid-based arrays
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector),save::ind_father
  ! Particle-based arrays
  integer,dimension(1:nvector),save::ind_son,igrid_son
  integer,dimension(1:nvector),save::list1,list2
  logical,dimension(1:nvector),save::ok
  real(dp),dimension(1:3)::skip_loc

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  xbound(1:3)=(/dble(nx),dble(ny),dble(nz)/)
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
     ind_father(i)=father(ind_grid(i))
  end do
  call get3cubefather(ind_father,nbors_father_cells,nbors_father_grids,ng,ilevel)

  ! Compute particle position in 3-cube
  error=.false.
  ind_son(1:np)=1
  ok(1:np)=.false.
  do idim=1,ndim
     do j=1,np
        i=floor((xp(ind_part(j),idim)/scale+skip_loc(idim)-x0(ind_grid_part(j),idim))/dx/2.0D0)
        if(i<0.or.i>2)error=.true.
        i=MAX(i,0)
        i=MIN(i,2)
        ind_son(j)=ind_son(j)+i*3**(idim-1)
        ! Check if particle has escaped from its parent grid
        ok(j)=ok(j).or.i.ne.1
     end do
  end do
  if(error)then
     write(*,*)'Problem in check_tree at level ',ilevel
     write(*,*)'A particle has moved outside allowed boundaries'
     do idim=1,ndim
        do j=1,np
           i=floor((xp(ind_part(j),idim)/scale+skip_loc(idim)-x0(ind_grid_part(j),idim))/dx/2.0D0)
           if(i<0.or.i>2)then
              write(*,*)xp(ind_part(j),1:ndim)
              write(*,*)x0(ind_grid_part(j),1:ndim)*scale
           endif
        end do
     end do
     stop
  end if

  ! Compute neighboring grid index
  do j=1,np
     igrid_son(j)=son(nbors_father_cells(ind_grid_part(j),ind_son(j)))
  end do

  ! If escaped particle sits in unrefined cell, leave it to its parent grid.
  ! For ilevel=levelmin, this should never happen.
  do j=1,np
     if(igrid_son(j)==0)ok(j)=.false.
  end do

  ! Periodic box
  do idim=1,ndim
     do j=1,np
        if(ok(j))then
           xxx=xp(ind_part(j),idim)/scale+skip_loc(idim)-xg(igrid_son(j),idim)
           if(xxx> xbound(idim)/2.0)then
              xp(ind_part(j),idim)=xp(ind_part(j),idim)-(xbound(idim)-skip_loc(idim))*scale
           endif
           if(xxx<-xbound(idim)/2.0)then
              xp(ind_part(j),idim)=xp(ind_part(j),idim)+(xbound(idim)-skip_loc(idim))*scale
           endif
        endif
     enddo
  enddo

  ! Switch particles linked list
  do j=1,np
     if(ok(j))then
        list1(j)=ind_grid(ind_grid_part(j))
        list2(j)=igrid_son(j)
     end if
  end do
  call remove_list(ind_part,list1,ok,np)
  call add_list(ind_part,list2,ok,np)

end subroutine check_tree
!################################################################
!################################################################
!################################################################
!################################################################
subroutine kill_tree_fine(ilevel)
  use pm_commons
  use amr_commons
  implicit none
  integer::ilevel
  !------------------------------------------------------------------------
  ! This routine sorts particle between ilevel grids and their
  ! ilevel+1 children grids. Particles are disconnected from their parent
  ! grid linked list and connected to their corresponding child grid linked
  ! list. If the  child grid does not exist, the particle is left to its
  ! original parent grid.
  !------------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,next_part
  integer::i,ig,ip,npart1,icpu
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part

  if(numbtot(1,ilevel)==0)return
  if(ilevel==nlevelmax)return
  if(numbtot(1,ilevel+1)==0)return
  if(verbose)write(*,111)ilevel

  ! Reset all linked lists at level ilevel+1
  do i=1,active(ilevel+1)%ngrid
     headp(active(ilevel+1)%igrid(i))=0
     tailp(active(ilevel+1)%igrid(i))=0
     numbp(active(ilevel+1)%igrid(i))=0
  end do
  do icpu=1,ncpu
     do i=1,reception(icpu,ilevel+1)%ngrid
#ifdef LIGHT_MPI_COMM
        headp(reception(icpu,ilevel+1)%pcomm%igrid(i))=0
        tailp(reception(icpu,ilevel+1)%pcomm%igrid(i))=0
        numbp(reception(icpu,ilevel+1)%pcomm%igrid(i))=0
#else
        headp(reception(icpu,ilevel+1)%igrid(i))=0
        tailp(reception(icpu,ilevel+1)%igrid(i))=0
        numbp(reception(icpu,ilevel+1)%igrid(i))=0
#endif
     end do
  end do

  ! Sort particles between ilevel and ilevel+1

  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        if(npart1>0)then
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if(ig==0)then
                 ig=1
                 ind_grid(ig)=igrid
              end if
              ip=ip+1
              ind_part(ip)=ipart
              ind_grid_part(ip)=ig
              if(ip==nvector)then
                 call kill_tree(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
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
     if(ip>0)call kill_tree(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
  end do
  ! End loop over cpus

111 format('   Entering kill_tree_fine for level ',I2)

end subroutine kill_tree_fine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine kill_tree(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use pm_commons
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !-----------------------------------------------------------------------
  ! This routine is called by subroutine kill_tree_fine.
  !-----------------------------------------------------------------------
  integer::i,j,idim,nx_loc
  real(dp)::dx,xxx,scale
  ! Grid based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
  ! Particle based arrays
  integer,dimension(1:nvector),save::igrid_son,ind_son
  integer,dimension(1:nvector),save::list1,list2
  logical,dimension(1:nvector),save::ok
  real(dp),dimension(1:3)::skip_loc

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)

  ! Compute lower left corner of grid
  do idim=1,ndim
     do i=1,ng
        x0(i,idim)=xg(ind_grid(i),idim)-dx
     end do
  end do

  ! Select only particles within grid boundaries
  ok(1:np)=.true.
  do idim=1,ndim
     do j=1,np
        xxx=(xp(ind_part(j),idim)/scale+skip_loc(idim)-x0(ind_grid_part(j),idim))/dx
        ok(j)=ok(j) .and. (xxx >= 0d0 .and. xxx < 2.0d0)
     end do
  end do

  ! Determines in which son particles sit
  ind_son(1:np)=0
  do idim=1,ndim
     do j=1,np
        i=int((xp(ind_part(j),idim)/scale+skip_loc(idim)-x0(ind_grid_part(j),idim))/dx)
        ind_son(j)=ind_son(j)+i*2**(idim-1)
     end do
  end do
  do j=1,np
     ind_son(j)=ncoarse+ind_son(j)*ngridmax+ind_grid(ind_grid_part(j))
  end do

  ! Determine which son cell is refined
  igrid_son(1:np)=0
  do j=1,np
     if(ok(j))igrid_son(j)=son(ind_son(j))
  end do
  do j=1,np
     ok(j)=igrid_son(j)>0
  end do

  ! Compute particle linked list
  do j=1,np
     if(ok(j))then
        list1(j)=ind_grid(ind_grid_part(j))
        list2(j)=igrid_son(j)
     end if
  end do

  ! Remove particles from their original linked lists
  call remove_list(ind_part,list1,ok,np)
  ! Add particles to their new linked lists
  call add_list(ind_part,list2,ok,np)

end subroutine kill_tree
!################################################################
!################################################################
!################################################################
!################################################################
subroutine merge_tree_fine(ilevel)
  use pm_commons
  use amr_commons
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------
  ! This routine disconnects all particles contained in children grids
  ! and connects them to their parent grid linked list.
  !---------------------------------------------------------------
  integer::igrid,iskip,icpu
  integer::i,ind,ncache,ngrid
  integer,dimension(1:nvector),save::ind_grid,ind_cell,ind_grid_son
  logical,dimension(1:nvector),save::ok

  if(numbtot(1,ilevel)==0)return
  if(ilevel==nlevelmax)return
  if(verbose)write(*,111)ilevel

  ! Loop over cpus
  do icpu=1,ncpu
     if(icpu==myid)then
        ncache=active(ilevel)%ngrid
     else
        ncache=reception(icpu,ilevel)%ngrid
     end if
     ! Loop over grids by vector sweeps
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        if(icpu==myid)then
           do i=1,ngrid
              ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
           end do
        else
           do i=1,ngrid
#ifdef LIGHT_MPI_COMM
              ind_grid(i)=reception(icpu,ilevel)%pcomm%igrid(igrid+i-1)
#else
              ind_grid(i)=reception(icpu,ilevel)%igrid(igrid+i-1)
#endif
           end do
        end if
        ! Loop over children grids
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do
           do i=1,ngrid
              ind_grid_son(i)=son(ind_cell(i))
           end do
           do i=1,ngrid
              ok(i)=ind_grid_son(i)>0
           end do
           do i=1,ngrid
           if(ok(i))then
           if(numbp(ind_grid_son(i))>0)then
              if(numbp(ind_grid(i))>0)then
                 ! Connect son linked list at the tail of father linked list
                 nextp(tailp(ind_grid(i)))=headp(ind_grid_son(i))
                 prevp(headp(ind_grid_son(i)))=tailp(ind_grid(i))
                 numbp(ind_grid(i))=numbp(ind_grid(i))+numbp(ind_grid_son(i))
                 tailp(ind_grid(i))=tailp(ind_grid_son(i))
              else
                 ! Initialize father linked list
                 headp(ind_grid(i))=headp(ind_grid_son(i))
                 tailp(ind_grid(i))=tailp(ind_grid_son(i))
                 numbp(ind_grid(i))=numbp(ind_grid_son(i))
              end if

           end if
           end if
           end do
        end do
        ! End loop over children
     end do
     ! End loop over grids
  end do
  ! End loop over cpus

111 format('   Entering merge_tree_fine for level ',I2)

end subroutine merge_tree_fine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine virtual_tree_fine(ilevel)
  use pm_commons
  use amr_commons
  use mpi_mod
  implicit none
  integer::ilevel
  !-----------------------------------------------------------------------
  ! This subroutine move particles across processors boundaries.
  !-----------------------------------------------------------------------
#ifndef WITHOUTMPI
#ifdef LIGHT_MPI_COMM
  integer:: offset_np, iactive
#endif
  integer::ip,ipcom,npart1,next_part,ncache,ncache_tot
  integer::icpu,igrid,ipart,jpart
  integer::info,buf_count,tagf=102,tagu=102
  integer::countsend,countrecv
  integer,dimension(MPI_STATUS_SIZE,2*ncpu)::statuses
  integer,dimension(2*ncpu)::reqsend,reqrecv
  integer,dimension(ncpu)::sendbuf,recvbuf
  logical::ok_free
  integer::particle_data_width, particle_data_width_int
  integer,dimension(1:nvector),save::ind_part,ind_list,ind_com
  ! MC tracer
  real(dp) :: dx, d2min, d2, x1(1:ndim), x2(1:ndim)
  integer :: ipart2, jpart2
#endif

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

#ifndef WITHOUTMPI
  dx=0.5D0**ilevel
  ! Count particle sitting in virtual boundaries
  do icpu=1,ncpu
     reception(icpu,ilevel)%npart=0
     do igrid=1,reception(icpu,ilevel)%ngrid
#ifdef LIGHT_MPI_COMM
        reception(icpu,ilevel)%npart=reception(icpu,ilevel)%npart+&
             & numbp(reception(icpu,ilevel)%pcomm%igrid(igrid))
#else
        reception(icpu,ilevel)%npart=reception(icpu,ilevel)%npart+&
             & numbp(reception(icpu,ilevel)%igrid(igrid))
#endif
     end do
     sendbuf(icpu)=reception(icpu,ilevel)%npart
  end do

  ! Calculate how many particle properties are being transferred
  ! igrid, level, id, families
  particle_data_width_int = 4
  if (MC_tracer) then
     ! Also send partp
     particle_data_width_int = particle_data_width_int + 1
  end if
  particle_data_width = twondim+1
  if(star.or.sink) then
     if(metal) then
        particle_data_width=twondim+3
     else
        particle_data_width=twondim+2
     endif
  endif

#ifdef OUTPUT_PARTICLE_POTENTIAL
  particle_data_width=particle_data_width+1
#endif

  ! Allocate communication buffer in emission
  do icpu=1,ncpu
     ncache=reception(icpu,ilevel)%npart
     if(ncache>0)then
        ! Allocate reception buffer
#ifdef LIGHT_MPI_COMM
        allocate(reception(icpu,ilevel)%pcomm%f8(1:particle_data_width_int,1:ncache))
        allocate(reception(icpu,ilevel)%pcomm%u(1:particle_data_width,1:ncache))
#else
        allocate(reception(icpu,ilevel)%fp(1:ncache,1:particle_data_width_int))
        allocate(reception(icpu,ilevel)%up(1:ncache,1:particle_data_width))
#endif
     end if
  end do
  if (MC_tracer) then
     ! Use itmpp to store the index within communicator
     ! Note: itmpp is also used in `sink_particle_tracer` for
     ! `gas_tracers`, so there is no interference here.
     do icpu=1,ncpu
        if(reception(icpu,ilevel)%npart>0)then
           ! Gather particles by vector sweeps
           ipcom=0
           do igrid=1,reception(icpu,ilevel)%ngrid
#ifdef LIGHT_MPI_COMM
              npart1=numbp(reception(icpu,ilevel)%pcomm%igrid(igrid))
              ipart =headp(reception(icpu,ilevel)%pcomm%igrid(igrid))
#else
              npart1=numbp(reception(icpu,ilevel)%igrid(igrid))
              ipart =headp(reception(icpu,ilevel)%igrid(igrid))
#endif
              ! Store index within communicator for stars
              do jpart = 1, npart1
                 ipcom = ipcom+1
                 if (is_star(typep(ipart))) then
                    itmpp(ipart) = ipcom
                 end if

                 ipart = nextp(ipart)
              end do
           end do
        end if
     end do
  end if

  ! Gather particle in communication buffer
  do icpu=1,ncpu
     if(reception(icpu,ilevel)%npart>0)then
     ! Gather particles by vector sweeps
     ipcom=0
     ip=0
     do igrid=1,reception(icpu,ilevel)%ngrid
#ifdef LIGHT_MPI_COMM
        npart1=numbp(reception(icpu,ilevel)%pcomm%igrid(igrid))
        ipart =headp(reception(icpu,ilevel)%pcomm%igrid(igrid))
#else
        npart1=numbp(reception(icpu,ilevel)%igrid(igrid))
        ipart =headp(reception(icpu,ilevel)%igrid(igrid))
#endif
        ! Loop over particles
        do jpart=1,npart1
           ! Save next particle  <--- Very important !!!
           next_part=nextp(ipart)
           ip=ip+1
           ipcom=ipcom+1
           ind_com (ip)=ipcom
           ind_part(ip)=ipart
#ifdef LIGHT_MPI_COMM
           ind_list(ip)=reception(icpu,ilevel)%pcomm%igrid(igrid)
           reception(icpu,ilevel)%pcomm%f8(1,ipcom)=igrid
#else
           ind_list(ip)=reception(icpu,ilevel)%igrid(igrid)
           reception(icpu,ilevel)%fp(ipcom,1)=igrid
#endif
           if(ip==nvector)then
              call fill_comm(ind_part,ind_com,ind_list,ip,ilevel,icpu)
              ip=0
           end if
           ipart=next_part  ! Go to next particle
        end do
     end do
     if(ip>0)call fill_comm(ind_part,ind_com,ind_list,ip,ilevel,icpu)
     end if
  end do

  ! Communicate virtual particle number to parent cpu
  call MPI_ALLTOALL(sendbuf,1,MPI_INTEGER,recvbuf,1,MPI_INTEGER,MPI_COMM_WORLD,info)

#ifdef LIGHT_MPI_COMM
  ! Count particles
  emission_part(ilevel)%nactive=0
  emission_part(ilevel)%nparts_tot=0
  do icpu=1,ncpu
     ncache=recvbuf(icpu)
     ! Cumulated counter of particles to send
     if(ncache>0) then
        emission_part(ilevel)%nactive=emission_part(ilevel)%nactive+1
        emission_part(ilevel)%nparts_tot=emission_part(ilevel)%nparts_tot+ncache
     end if
  end do
#endif

  ! Allocate communication buffer in reception
#ifdef LIGHT_MPI_COMM
  if (emission_part(ilevel)%nactive > 0) then
    ! Allocate reception buffer
    allocate(emission_part(ilevel)%cpuid(emission_part(ilevel)%nactive))
    allocate(emission_part(ilevel)%nparts(emission_part(ilevel)%nactive))
    allocate(emission_part(ilevel)%f8(particle_data_width_int,emission_part(ilevel)%nparts_tot))
    allocate(emission_part(ilevel)%u(particle_data_width,emission_part(ilevel)%nparts_tot))
  endif

  iactive=1
  do icpu=1,ncpu
    ncache=recvbuf(icpu)
    if (ncache>0) then
      emission_part(ilevel)%nparts(iactive) = ncache
      emission_part(ilevel)%cpuid(iactive) = icpu
      iactive=iactive+1
    endif
  end do
#else
  do icpu=1,ncpu
      emission(icpu,ilevel)%npart=recvbuf(icpu)
     ncache=emission(icpu,ilevel)%npart
     if(ncache>0)then
        ! Allocate reception buffer
        allocate(emission(icpu,ilevel)%fp(1:ncache,1:particle_data_width_int))
        allocate(emission(icpu,ilevel)%up(1:ncache,1:particle_data_width))
     end if
  end do
#endif

#ifdef LIGHT_MPI_COMM
  ! Receive particles
  countrecv=0
  offset_np=1
  do iactive=1,emission_part(ilevel)%nactive
     ncache=emission_part(ilevel)%nparts(iactive)
     buf_count=ncache*particle_data_width_int
     countrecv=countrecv+1
#ifndef LONGINT
     call MPI_IRECV(emission_part(ilevel)%f8(1,offset_np),buf_count, &
          & MPI_INTEGER,emission_part(ilevel)%cpuid(iactive)-1, &
          & tagf,MPI_COMM_WORLD,reqrecv(countrecv),info)
#else
     call MPI_IRECV(emission_part(ilevel)%f8(1,offset_np),buf_count, &
          & MPI_INTEGER8,emission_part(ilevel)%cpuid(iactive)-1, &
          & tagf,MPI_COMM_WORLD,reqrecv(countrecv),info)
#endif
     buf_count=ncache*particle_data_width
     countrecv=countrecv+1
     call MPI_IRECV(emission_part(ilevel)%u(1,offset_np),buf_count, &
          & MPI_DOUBLE_PRECISION,emission_part(ilevel)%cpuid(iactive)-1, &
          & tagu,MPI_COMM_WORLD,reqrecv(countrecv),info)
     offset_np=offset_np+ncache
  end do
#else
  ! Receive particles
  countrecv=0
  do icpu=1,ncpu
     ncache=emission(icpu,ilevel)%npart
     if(ncache>0)then
        buf_count=ncache*particle_data_width_int
        countrecv=countrecv+1
#ifndef LONGINT
        call MPI_IRECV(emission(icpu,ilevel)%fp,buf_count, &
             & MPI_INTEGER,icpu-1,&
             & tagf,MPI_COMM_WORLD,reqrecv(countrecv),info)
#else
        call MPI_IRECV(emission(icpu,ilevel)%fp,buf_count, &
             & MPI_INTEGER8,icpu-1,&
             & tagf,MPI_COMM_WORLD,reqrecv(countrecv),info)
#endif
        buf_count=ncache*particle_data_width
        countrecv=countrecv+1
        call MPI_IRECV(emission(icpu,ilevel)%up,buf_count, &
             & MPI_DOUBLE_PRECISION,icpu-1,&
             & tagu,MPI_COMM_WORLD,reqrecv(countrecv),info)
     end if
  end do
#endif

  ! Send particles
  countsend=0
  do icpu=1,ncpu
     ncache=reception(icpu,ilevel)%npart
     if(ncache>0)then
        buf_count=ncache*particle_data_width_int
        countsend=countsend+1
#ifdef LIGHT_MPI_COMM
#ifndef LONGINT
        call MPI_ISEND(reception(icpu,ilevel)%pcomm%f8,buf_count, &
             & MPI_INTEGER,icpu-1,&
             & tagf,MPI_COMM_WORLD,reqsend(countsend),info)
#else
        call MPI_ISEND(reception(icpu,ilevel)%pcomm%f8,buf_count, &
             & MPI_INTEGER8,icpu-1,&
             & tagf,MPI_COMM_WORLD,reqsend(countsend),info)
#endif
        buf_count=ncache*particle_data_width
        countsend=countsend+1
        call MPI_ISEND(reception(icpu,ilevel)%pcomm%u,buf_count, &
             & MPI_DOUBLE_PRECISION,icpu-1,&
             & tagu,MPI_COMM_WORLD,reqsend(countsend),info)
#else
#ifndef LONGINT
        call MPI_ISEND(reception(icpu,ilevel)%fp,buf_count, &
             & MPI_INTEGER,icpu-1,&
             & tagf,MPI_COMM_WORLD,reqsend(countsend),info)
#else
        call MPI_ISEND(reception(icpu,ilevel)%fp,buf_count, &
             & MPI_INTEGER8,icpu-1,&
             & tagf,MPI_COMM_WORLD,reqsend(countsend),info)
#endif
        buf_count=ncache*particle_data_width
        countsend=countsend+1
        call MPI_ISEND(reception(icpu,ilevel)%up,buf_count, &
             & MPI_DOUBLE_PRECISION,icpu-1,&
             & tagu,MPI_COMM_WORLD,reqsend(countsend),info)
#endif
     end if
  end do

  ! Wait for full completion of receives
  call MPI_WAITALL(countrecv,reqrecv,statuses,info)

  ! Compute total number of newly created particles
#ifdef LIGHT_MPI_COMM
  ncache_tot = emission_part(ilevel)%nparts_tot
#else
  ncache_tot=0
  do icpu=1,ncpu
     ncache_tot=ncache_tot+emission(icpu,ilevel)%npart
  end do
#endif

  ! Wait for full completion of sends
  call MPI_WAITALL(countsend,reqsend,statuses,info)

  call MPI_ALLREDUCE(numbp_free,numbp_free_tot,1,MPI_INTEGER,MPI_MIN,&
       & MPI_COMM_WORLD,info)
  ok_free=(numbp_free-ncache_tot)>=0
  if(.not. ok_free)then
     write(*,*)'No more free memory for particles'
     write(*,*)'Increase npartmax'
     write(*,*)numbp_free,ncache_tot
     write(*,*)myid
#ifdef LIGHT_MPI_COMM
     ! Write list of communicating CPU ids + list of number of exchanged particles
     write(*,*)emission_part(ilevel)%nparts(:)
     write(*,*)emission_part(ilevel)%cpuid(:)
#else
     write(*,*)emission(1:ncpu,ilevel)%npart
#endif
     write(*,*)'============================'
     write(*,*)reception(1:ncpu,ilevel)%npart
     call MPI_ABORT(MPI_COMM_WORLD,1,info)
  end if

  ! Scatter new particles from communication buffer
#ifdef LIGHT_MPI_COMM
  offset_np=1
  do iactive=1,emission_part(ilevel)%nactive
     ! Loop over particles by vector sweeps
     icpu=emission_part(ilevel)%cpuid(iactive)
     ncache=emission_part(ilevel)%nparts(iactive)
#else
  do icpu=1,ncpu
     ! Loop over particles by vector sweeps
     ncache=emission(icpu,ilevel)%npart
#endif
     do ipart=1,ncache,nvector
        npart1=min(nvector,ncache-ipart+1)
        do ip=1,npart1
           ind_com(ip)=ipart+ip-1
        end do
#ifdef LIGHT_MPI_COMM
        call empty_comm(ind_com,npart1,ilevel,iactive,offset_np,particle_data_width, particle_data_width_int)
#else
        call empty_comm(ind_com,npart1,ilevel,icpu)
#endif
     end do
     ! Loop on star tracers in the communicator
     if (MC_tracer) then
        do ipart = 1, ncache
#ifdef LIGHT_MPI_COMM
           jpart = emission_part(ilevel)%f8(1,offset_np+ipart-1)
#else
           jpart = emission(icpu,ilevel)%fp(ipart,1)
#endif
           ! Get index of star within current CPU
           if (is_star_tracer(typep(jpart))) then
              ! Note: the partp array should store the index of the
              ! star within the communicator. However, sometimes
              ! (why?) this index is out of bounds (either 0 or
              ! greater than size of communicator). In this case, we
              ! find the star at the position of the tracer.
#ifdef LIGHT_MPI_COMM
              if ( (partp(jpart) > 0) .and. &
                   (partp(jpart) <= ncache)) then
                 partp(jpart) = emission_part(ilevel)%f8(1,offset_np+partp(jpart)-1)
#else
              if ( (partp(jpart) > 0) .and. &
                   (partp(jpart) <= size(emission(icpu,ilevel)%fp(:, 1))) ) then
                 partp(jpart) = emission(icpu,ilevel)%fp(partp(jpart), 1)
#endif
              else
                 d2min = (2*dx)**2
                 ! Try to find the star in the emission buffer
                 partp(jpart) = 0
                 x1(:) = xp(jpart, :)

                 do ipart2 = 1, ncache
#ifdef LIGHT_MPI_COMM
                    jpart2 = emission_part(ilevel)%f8(1,offset_np+ipart2-1)
#else
                    jpart2 = emission(icpu,ilevel)%fp(ipart2, 1)
#endif
                    if (is_star(typep(jpart2))) then
                       ! Check there is a star closer than dx. If
                       ! there are multiple, take the closest.
                       x2(:) = xp(jpart2, :)
                       if (all(abs(x2(:) - x1(:)) <= dx)) then
                          d2 = sum((x2(:) - x1(:))**2)
                          if (d2 < d2min) then
                             partp(jpart) = jpart2
                             d2min = d2
                          end if
                       end if
                    end if
                 end do
                 if (partp(jpart) == 0) then
                    write(*, *) 'An error occurred in virtual_tree_fine while converting star ids'
                    write(*, *) myid, '<-', icpu, '>< converting back', jpart, partp(jpart), xp(jpart, :)
                    ! stop
                    typep(jpart)%family = FAM_TRACER_GAS
                 end if
              end if
           end if
        end do
     end if
#ifdef LIGHT_MPI_COMM
     offset_np=offset_np+ncache
#endif
  end do

  ! Deallocate temporary communication buffers
#ifdef LIGHT_MPI_COMM
  if(emission_part(ilevel)%nactive>0)then
    deallocate(emission_part(ilevel)%cpuid)
    deallocate(emission_part(ilevel)%nparts)
    deallocate(emission_part(ilevel)%f8)
    deallocate(emission_part(ilevel)%u)
  endif
  do icpu=1,ncpu
     ncache=reception(icpu,ilevel)%npart
     if(ncache>0)then
       deallocate(reception(icpu,ilevel)%pcomm%f8)
       deallocate(reception(icpu,ilevel)%pcomm%u)
     end if
  end do
#else
  do icpu=1,ncpu
     ncache=emission(icpu,ilevel)%npart
     if(ncache>0)then
        deallocate(emission(icpu,ilevel)%fp)
        deallocate(emission(icpu,ilevel)%up)
     end if
     ncache=reception(icpu,ilevel)%npart
     if(ncache>0)then
        deallocate(reception(icpu,ilevel)%fp)
        deallocate(reception(icpu,ilevel)%up)
     end if
  end do
#endif

#endif

111 format('   Entering virtual_tree_fine for level ',I2)
end subroutine virtual_tree_fine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine fill_comm(ind_part,ind_com,ind_list,np,ilevel,icpu)
  use pm_commons
  use amr_commons
  implicit none
  integer::np,ilevel,icpu
  integer,dimension(1:nvector)::ind_part,ind_com,ind_list
  integer::current_property
  integer::i,idim
  logical,dimension(1:nvector),save::ok=.true.

  ! Gather particle level and identity
  do i=1,np
#ifdef LIGHT_MPI_COMM
     reception(icpu,ilevel)%pcomm%f8(2,ind_com(i))=levelp(ind_part(i))
     reception(icpu,ilevel)%pcomm%f8(3,ind_com(i))=idp   (ind_part(i))
     reception(icpu,ilevel)%pcomm%f8(4,ind_com(i))=part2int(typep(ind_part(i)))
#else
     reception(icpu,ilevel)%fp(ind_com(i),2)=levelp(ind_part(i))
     reception(icpu,ilevel)%fp(ind_com(i),3)=idp   (ind_part(i))
     reception(icpu,ilevel)%fp(ind_com(i),4)=part2int(typep(ind_part(i)))
#endif
  end do

  ! Gather particle position and velocity
  do idim=1,ndim
     do i=1,np
#ifdef LIGHT_MPI_COMM
        reception(icpu,ilevel)%pcomm%u(idim,ind_com(i)     )=xp(ind_part(i),idim)
        reception(icpu,ilevel)%pcomm%u(idim+ndim,ind_com(i))=vp(ind_part(i),idim)
#else
        reception(icpu,ilevel)%up(ind_com(i),idim     )=xp(ind_part(i),idim)
        reception(icpu,ilevel)%up(ind_com(i),idim+ndim)=vp(ind_part(i),idim)
#endif
     end do
  end do

  current_property = twondim+1
  ! Gather particle mass
  do i=1,np
#ifdef LIGHT_MPI_COMM
     reception(icpu,ilevel)%pcomm%u(current_property,ind_com(i))=mp(ind_part(i))
#else
     reception(icpu,ilevel)%up(ind_com(i),current_property)=mp(ind_part(i))
#endif
  end do
  current_property = current_property+1

#ifdef OUTPUT_PARTICLE_POTENTIAL
  ! Gather particle potential
  do i=1,np
#ifdef LIGHT_MPI_COMM
     reception(icpu,ilevel)%pcomm%u(current_property,ind_com(i))=ptcl_phi(ind_part(i))
#else
     reception(icpu,ilevel)%up(ind_com(i),current_property)=ptcl_phi(ind_part(i))
#endif
  end do
  current_property = current_property+1
#endif

  ! Gather particle birth epoch
  if(star.or.sink)then
     do i=1,np
#ifdef LIGHT_MPI_COMM
        reception(icpu,ilevel)%pcomm%u(current_property,ind_com(i))=tp(ind_part(i))
#else
        reception(icpu,ilevel)%up(ind_com(i),current_property)=tp(ind_part(i))
#endif
     end do
     current_property = current_property+1
     if(metal)then
        do i=1,np
#ifdef LIGHT_MPI_COMM
           reception(icpu,ilevel)%pcomm%u(current_property,ind_com(i))=zp(ind_part(i))
#else
           reception(icpu,ilevel)%up(ind_com(i),current_property)=zp(ind_part(i))
#endif
        end do
        current_property = current_property+1
     end if
  end if
  ! MC Tracer
  if (MC_tracer) then
     do i=1,np
        if (is_star_tracer(typep(ind_part(i)))) then
           ! Store index of the star *within* communicator
#ifdef LIGHT_MPI_COMM
           reception(icpu, ilevel)%pcomm%f8(5, ind_com(i)) = &
                itmpp(partp(ind_part(i)))
        else
           reception(icpu, ilevel)%pcomm%f8(5, ind_com(i)) = &
                partp(ind_part(i))
#else
           reception(icpu, ilevel)%fp(ind_com(i), 5) = &
                itmpp(partp(ind_part(i)))
        else
           reception(icpu, ilevel)%fp(ind_com(i), 5) = &
                partp(ind_part(i))
#endif
        end if
     end do
  end if

  ! Remove particles from parent linked list
  call remove_list(ind_part,ind_list,ok,np)
  call add_free(ind_part,np)

end subroutine fill_comm
!################################################################
!################################################################
!################################################################
!################################################################
#ifdef LIGHT_MPI_COMM
subroutine empty_comm(ind_com,np,ilevel,iactive,offset_np,particle_data_width,particle_data_width_int)
#else
subroutine empty_comm(ind_com,np,ilevel,icpu)
#endif
  use pm_commons
  use amr_commons
  implicit none
  integer::np,ilevel
#ifdef LIGHT_MPI_COMM
  integer::iactive,offset_np,offset_ig,found_cpu,particle_data_width,particle_data_width_int,j,nparts
#else
  integer::icpu
#endif
  integer,dimension(1:nvector)::ind_com

  integer::i,idim,igrid
  integer,dimension(1:nvector),save::ind_list,ind_part
  logical,dimension(1:nvector),save::ok=.true.
  integer::current_property

#ifdef LIGHT_MPI_COMM
  offset_ig=0
  found_cpu=-1
  do j=1,emission(ilevel)%nactive
    if(emission(ilevel)%cpuid(j)==emission_part(ilevel)%cpuid(iactive)) then
      found_cpu=1
      exit
    end if
    offset_ig=offset_ig+emission(ilevel)%ngrids(j)
  end do
  if(np>0.and.found_cpu==-1) then
    write(*, *) '[CPU #',myid,'] An error occurred in particle_tree->empty_comm while '
    write(*, *) 'searching for CPU #',emission_part(ilevel)%cpuid(iactive),' grid index.'
  endif
  nparts=emission_part(ilevel)%nparts(iactive)
#endif

  ! Compute parent grid index
  do i=1,np
#ifdef LIGHT_MPI_COMM
     igrid=int(emission_part(ilevel)%f8(1, offset_np+ind_com(i)-1), 4)
     ind_list(i)=emission(ilevel)%igrid(offset_ig+igrid)
#else
     igrid=int(emission(icpu,ilevel)%fp(ind_com(i),1), 4)
     ind_list(i)=emission(icpu,ilevel)%igrid(igrid)
#endif
  end do

  ! Add particle to parent linked list
  call remove_free(ind_part,np)
  call add_list(ind_part,ind_list,ok,np)

  ! Scatter particle level and identity
  do i=1,np
#ifdef LIGHT_MPI_COMM
     levelp(ind_part(i))=int(emission_part(ilevel)%f8(2, offset_np+ind_com(i)-1), 4)
     idp   (ind_part(i))=int(emission_part(ilevel)%f8(3, offset_np+ind_com(i)-1))
     typep(ind_part(i)) =int2part(int(emission_part(ilevel)%f8(4, offset_np+ind_com(i)-1), 4))
#else
     levelp(ind_part(i))=int(emission(icpu,ilevel)%fp(ind_com(i),2), 4)
     idp   (ind_part(i))=int(emission(icpu,ilevel)%fp(ind_com(i),3))
     typep(ind_part(i)) =int2part(int(emission(icpu,ilevel)%fp(ind_com(i),4), 4))
#endif
  end do

  ! Scatter particle position and velocity
  do idim=1,ndim
    do i=1,np
#ifdef LIGHT_MPI_COMM
     xp(ind_part(i),idim)=emission_part(ilevel)%u(idim, offset_np+ind_com(i)-1)
     vp(ind_part(i),idim)=emission_part(ilevel)%u(ndim+idim,offset_np+ind_com(i)-1)
#else
     xp(ind_part(i),idim)=emission(icpu,ilevel)%up(ind_com(i),idim     )
     vp(ind_part(i),idim)=emission(icpu,ilevel)%up(ind_com(i),idim+ndim)
#endif
    end do
  end do

  current_property = twondim+1

  ! Scatter particle mass
  do i=1,np
#ifdef LIGHT_MPI_COMM
     mp(ind_part(i))=emission_part(ilevel)%u(current_property, offset_np+ind_com(i)-1)
#else
     mp(ind_part(i))=emission(icpu,ilevel)%up(ind_com(i),current_property)
#endif
  end do
  current_property = current_property+1

#ifdef OUTPUT_PARTICLE_POTENTIAL
  ! Scatter particle phi
  do i=1,np
#ifdef LIGHT_MPI_COMM
     ptcl_phi(ind_part(i))=emission_part(ilevel)%u(current_property, offset_np+ind_com(i)-1)
#else
     ptcl_phi(ind_part(i))=emission(icpu,ilevel)%up(ind_com(i),current_property)
#endif
  end do
  current_property = current_property+1
#endif

  ! Scatter particle birth epoch
  if(star.or.sink)then
     do i=1,np
#ifdef LIGHT_MPI_COMM
        tp(ind_part(i))=emission_part(ilevel)%u(current_property, offset_np+ind_com(i)-1)
#else
        tp(ind_part(i))=emission(icpu,ilevel)%up(ind_com(i),current_property)
#endif
     end do
     current_property = current_property+1
     if(metal)then
        do i=1,np
#ifdef LIGHT_MPI_COMM
           zp(ind_part(i))=emission_part(ilevel)%u(current_property, offset_np+ind_com(i)-1)
#else
           zp(ind_part(i))=emission(icpu,ilevel)%up(ind_com(i),current_property)
#endif
        end do
        current_property = current_property+1
     end if
  end if

  ! MC Tracer
  if (MC_tracer) then
     do i=1,np
        ! Store the target
        ! NB: this 'partp' contains for star tracers: the adress in
        ! the communicator of the star particle
#ifdef LIGHT_MPI_COMM
        partp(ind_part(i)) = emission_part(ilevel)%f8(5 ,offset_np+ind_com(i)-1)
#else
        partp(ind_part(i)) = emission(icpu,ilevel)%fp(ind_com(i), 5)
#endif

        ! Use the communicator as a tmp array mapping index in comm to index in array
        ! of all particles
#ifdef LIGHT_MPI_COMM
        emission_part(ilevel)%f8(1, offset_np+ind_com(i)-1) = ind_part(i)
#else
        emission(icpu,ilevel)%fp(ind_com(i), 1) = ind_part(i)
#endif
     end do
  end if

end subroutine empty_comm
!################################################################
!################################################################
!################################################################
!################################################################
