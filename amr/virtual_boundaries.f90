!################################################################
!################################################################
!################################################################
!################################################################
subroutine authorize_coarse
  use amr_commons
  implicit none
  !----------------------------------------------------------------------
  ! This routine authorizes all base cells for refinement.
  ! This duplicates all base grids over all cpu's.
  !----------------------------------------------------------------------
  integer::nxny,i,j,k,ind

  if(verbose)write(*,*)'  Entering authorize_coarse'
  ! Constants
  nxny=nx*ny
  ! Initialize flag2(0) to zero
  flag2(0)=0
  ! Duplicate full domain over cpus
  do k=0,nz-1
  do j=0,ny-1
  do i=0,nx-1
     ind=1+i+j*nx+k*nxny
     flag2(ind)=1
  end do
  end do
  end do

end subroutine authorize_coarse
!################################################################
!################################################################
!################################################################
!################################################################
subroutine authorize_fine(ilevel)
  use amr_commons
  use bisection
  implicit none
  integer::ilevel
  ! -------------------------------------------------------------------
  ! This routine computes the authorization map (flag2) for level ilevel.
  ! All myid cells are first marked for authorization.
  ! All virtual cells that intersect the local ordering domain are
  ! also marked for authorization. Finally, the routine performs
  ! a dilatation of the authorization map of one cell width.
  ! Array flag1 for virtual cells is used as temporary work space.
  ! -------------------------------------------------------------------
  integer::ismooth,ibound,ngrid,i,ncache,iskip,igrid,ind,icpu
  integer::ix,iy,iz,idim,nx_loc,isub
  integer,dimension(1:3)::n_nbor
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  integer,dimension(1:nvector,0:twondim),save::igridn
  real(dp)::dx,dx_loc,scale
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(qdp),dimension(1:nvector),save::order_min,order_max
  logical::test
  real(dp),dimension(1:ndim)::xmin,xmax

  if(ilevel==nlevelmax)return
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

  ! Scaling factor
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  ! Authorize all myid grids (needed for uploads)
  ncache=active(ilevel)%ngrid
  ! Loop over grids by vector sweeps
  do igrid=1,ncache,nvector
     ! Gather nvector grids
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
        do i=1,ngrid
           flag2(ind_cell(i))=1
        end do
     end do
     ! End loop over cells
  end do
  ! End loop over grids

  ! Authorize virtual cells that contains myid children cells
  do icpu=1,ncpu
     ncache=reception(icpu,ilevel)%ngrid
     ! Loop over grids by vector sweeps
     do igrid=1,ncache,nvector
        ! Gather nvector grids
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
#ifdef LIGHT_MPI_COMM
           ind_grid(i)=reception(icpu,ilevel)%pcomm%igrid(igrid+i-1)
#else
           ind_grid(i)=reception(icpu,ilevel)%igrid(igrid+i-1)
#endif
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
           ! Reset flag2
           do i=1,ngrid
              flag2(ind_cell(i))=0
           end do

           if (ordering /= 'bisection') then
              ! Compute minimum and maximum ordering key
              call cmp_minmaxorder(xx,order_min,order_max,dx_loc,ngrid)
              ! Determine if cell is authorized
              do isub=1,overload
                 do i=1,ngrid
                    if(    order_max(i)>bound_key(myid-1+(isub-1)*ncpu).and.&
                         & order_min(i)<bound_key(myid  +(isub-1)*ncpu) )then
                       flag2(ind_cell(i))=1
                    endif
                 end do
              end do
           else ! recursive bisection method
               do i=1,ngrid
                  ! Test if cell overlaps the cpu
                  test=.true.
                  xmin=xx(i,:)-0.5d0*dx_loc
                  xmax=xx(i,:)+0.5d0*dx_loc
                  do idim=1,ndim
                     ! This needs to be a >=, not a >, to precisely match the
                     ! ordering/=case for refinement flagging
                     test=test .and. (bisec_cpubox_max(myid,idim).ge.xmin(idim) &
                                          .and. bisec_cpubox_min(myid,idim).le.xmax(idim))
                  end do
                  if(test) flag2(ind_cell(i))=1
               end do
           endif

           ! For load balancing operations
           if(balance)then
              if(ordering/='bisection') then
                 do isub=1,overload
                    do i=1,ngrid
                       if(    order_max(i)>bound_key2(myid-1+(isub-1)*ncpu).and.&
                            & order_min(i)<bound_key2(myid  +(isub-1)*ncpu) )then
                          flag2(ind_cell(i))=1
                       endif
                    end do
                 end do
              else
                 do i=1,ngrid
                    ! Test if cell overlaps the cpu with new cpu map
                    test=.true.
                    xmin=xx(i,:)-0.5d0*dx_loc
                    xmax=xx(i,:)+0.5d0*dx_loc
                    do idim=1,ndim
                       ! This needs to be a >=, not a >, to precisely match the
                       ! ordering/=case for refinement flagging
                       test=test .and. (bisec_cpubox_max2(myid,idim).ge.xmin(idim) &
                            .and. bisec_cpubox_min2(myid,idim).le.xmax(idim))
                    end do
                    if(test) flag2(ind_cell(i))=1
                 end do
              end if
              do i=1,ngrid
                 if(cpu_map2(father(ind_grid(i)))==myid)then
                    flag2(ind_cell(i))=1
                 endif
              end do
           end if
        end do
        ! End loop over cells
     end do
     ! End loop over grids
  end do
  ! End loop over cpus

  ! Apply dilatation operator over flag2 cells on virtual cells only

  flag2(0)=0
  ! Set flag2 to 0 for physical boundary grids
  do ibound=1,nboundary
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,boundary(ibound,ilevel)%ngrid
        flag2(boundary(ibound,ilevel)%igrid(i)+iskip)=0
     end do
  end do
  end do

  ! Loop over steps
  do ibound=1,nexpand_bound
  n_nbor(1:3)=(/1,2,3/)
  do ismooth=1,ndim
     ! Initialize flag1 to 0 in virtual cells
     do icpu=1,ncpu
        ncache=reception(icpu,ilevel)%ngrid
        do igrid=1,ncache,nvector
           ngrid=MIN(nvector,ncache-igrid+1)
           do i=1,ngrid
#ifdef LIGHT_MPI_COMM
              ind_grid(i)=reception(icpu,ilevel)%pcomm%igrid(igrid+i-1)
#else
              ind_grid(i)=reception(icpu,ilevel)%igrid(igrid+i-1)
#endif
           end do
           do ind=1,twotondim
              iskip=ncoarse+(ind-1)*ngridmax
              do i=1,ngrid
                 ind_cell(i)=iskip+ind_grid(i)
              end do
              do i=1,ngrid
                 flag1(ind_cell(i))=0
              end do
           end do
        end do
     end do

     ! Count neighbors and set flag2 accordingly
     do icpu=1,ncpu
        ncache=reception(icpu,ilevel)%ngrid
        do igrid=1,ncache,nvector
           ngrid=MIN(nvector,ncache-igrid+1)
           do i=1,ngrid
#ifdef LIGHT_MPI_COMM
              ind_grid(i)=reception(icpu,ilevel)%pcomm%igrid(igrid+i-1)
#else
              ind_grid(i)=reception(icpu,ilevel)%igrid(igrid+i-1)
#endif
           end do
           call getnborgrids(ind_grid,igridn,ngrid)
           do ind=1,twotondim
              call count_nbors2(igridn,ind,n_nbor(ismooth),ngrid)
           end do
        end do
     end do

     ! Set flag2=1 for cells with flag1=1
     do icpu=1,ncpu
        ncache=reception(icpu,ilevel)%ngrid
        do igrid=1,ncache,nvector
           ngrid=MIN(nvector,ncache-igrid+1)
           do i=1,ngrid
#ifdef LIGHT_MPI_COMM
              ind_grid(i)=reception(icpu,ilevel)%pcomm%igrid(igrid+i-1)
#else
              ind_grid(i)=reception(icpu,ilevel)%igrid(igrid+i-1)
#endif
           end do
           do ind=1,twotondim
              iskip=ncoarse+(ind-1)*ngridmax
              do i=1,ngrid
                 ind_cell(i)=iskip+ind_grid(i)
              end do
              do i=1,ngrid
                 if(flag1(ind_cell(i))==1)flag2(ind_cell(i))=1
              end do
           end do
        end do
     end do

  end do
  ! End loop over steps
  end do

  ! Compute authorization map for physical boundaries
  if(simple_boundary)call init_boundary_fine(ilevel)

  ! Restore boundaries for flag1
  call make_virtual_fine_int(flag1(1),ilevel)
  if(simple_boundary)call make_boundary_flag(ilevel)

111 format('   Entering authorize_fine for level ',I2)

end subroutine authorize_fine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine make_virtual_coarse_int(xx)
  use amr_commons
  use mpi_mod
  implicit none
  integer,dimension(1:ncoarse+ngridmax*twotondim)::xx
  !-----------------------------------------------------------
  ! This routine communicates virtual boundary conditions
  ! at the coarse level for integer arrays.
  !-----------------------------------------------------------
  integer::nxny,ncell

#ifndef WITHOUTMPI
  integer::i,j,k
  integer::icell,info
  integer,dimension(:),allocatable::ind_cell,fff,ffg
#endif

  ! Constants
  nxny=nx*ny
  ncell=  (icoarse_max-icoarse_min+1) &
       & *(jcoarse_max-jcoarse_min+1) &
       & *(kcoarse_max-kcoarse_min+1)

#ifndef WITHOUTMPI
  ! Allocate local arrays
  allocate(ind_cell(1:ncell),fff(1:ncell),ffg(1:ncell))

  ! Compute cell indices
  icell=0
  do k=kcoarse_min,kcoarse_max
  do j=jcoarse_min,jcoarse_max
  do i=icoarse_min,icoarse_max
     icell=icell+1
     ind_cell(icell)=1+i+j*nx+k*nxny
  end do
  end do
  end do

  ! Communications
  fff=0; ffg=0
  do icell=1,ncell
     if(cpu_map(ind_cell(icell))==myid)fff(icell)=xx(ind_cell(icell))
  end do
  call MPI_ALLREDUCE(fff,ffg,ncell,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  do icell=1,ncell
     xx(ind_cell(icell))=ffg(icell)
  end do

  ! Dealocate local arrays
  deallocate(ind_cell,fff,ffg)
#endif

end subroutine make_virtual_coarse_int
!################################################################
!################################################################
!################################################################
!################################################################
subroutine make_virtual_fine_dp(xx,ilevel)
  use amr_commons
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer,dimension(MPI_STATUS_SIZE,ncpu)::statuses
#endif
  integer::ilevel
  real(dp),dimension(1:ncoarse+ngridmax*twotondim)::xx
  ! -------------------------------------------------------------------
  ! This routine communicates virtual boundaries among all cpu's.
  ! at level ilevel for any double precision array in the AMR grid.
  ! -------------------------------------------------------------------
#ifndef WITHOUTMPI

  integer::icpu,i,j,ncache,iskip,step
#ifdef LIGHT_MPI_COMM
  integer::ind,idx,offset
#endif
  integer::countsend,countrecv
  integer::info,tag=101
#ifdef LIGHT_MPI_COMM
  integer,dimension(ncpu)::reqrecv
  integer,allocatable,dimension(:)::reqsend
#else
   integer,dimension(ncpu)::reqsend,reqrecv
#endif

#endif

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

#ifndef WITHOUTMPI

#ifdef LIGHT_MPI_COMM
  allocate(reqsend(emission(ilevel)%nactive))

  ! Allocate temporary communication buffers
  do icpu=1,ncpu
     ncache=reception(icpu,ilevel)%ngrid
     if(ncache>0) then
        if (.not. associated(reception(icpu,ilevel)%pcomm)) then
           allocate(reception(icpu,ilevel)%pcomm)
        end if
        allocate(reception(icpu,ilevel)%pcomm%u(1:ncache*twotondim, 1:1))
     end if
  end do
  allocate(emission(ilevel)%u(1:emission(ilevel)%ngrids_tot*twotondim, 1:1))
#endif

  ! Receive all messages
  countrecv=0
  do icpu=1,ncpu
     ncache=reception(icpu,ilevel)%ngrid
     if(ncache>0) then
       countrecv=countrecv+1
#ifdef LIGHT_MPI_COMM
       call MPI_IRECV(reception(icpu,ilevel)%pcomm%u,ncache*twotondim, &
            & MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
#else
       call MPI_IRECV(reception(icpu,ilevel)%u,ncache*twotondim, &
            & MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
#endif
     end if
  end do

  ! Gather emission array
#ifdef LIGHT_MPI_COMM
  offset=0
  do idx=1,emission(ilevel)%nactive
    do j=1,twotondim
      step=(j-1)*emission(ilevel)%ngrids(idx)
      iskip=ncoarse+(j-1)*ngridmax
      do i=1,emission(ilevel)%ngrids(idx)
        emission(ilevel)%u(offset*twotondim+i+step,1)=xx(emission(ilevel)%igrid(offset+i)+iskip)
      end do
    end do
    offset=offset+emission(ilevel)%ngrids(idx)
  end do
#else
  do icpu=1,ncpu
     if (emission(icpu,ilevel)%ngrid>0) then
       do j=1,twotondim
         step=(j-1)*emission(icpu,ilevel)%ngrid
         iskip=ncoarse+(j-1)*ngridmax
         do i=1,emission(icpu,ilevel)%ngrid
           emission(icpu,ilevel)%u(i+step,1)=xx(emission(icpu,ilevel)%igrid(i)+iskip)
         end do
       end do
     end if
  end do
#endif

  ! Send all messages
#ifdef LIGHT_MPI_COMM
  offset=1
  do idx=1,emission(ilevel)%nactive
    ncache=emission(ilevel)%ngrids(idx)*twotondim
    call MPI_ISEND(emission(ilevel)%u(offset,1),ncache,MPI_DOUBLE_PRECISION,&
                   emission(ilevel)%cpuid(idx)-1,tag,MPI_COMM_WORLD,reqsend(idx),info)
    offset=offset+ncache
  end do
#else
  countsend=0
  do icpu=1,ncpu
     ncache=emission(icpu,ilevel)%ngrid
     if(ncache>0) then
       countsend=countsend+1
       call MPI_ISEND(emission(icpu,ilevel)%u,ncache*twotondim, &
            & MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
     end if
  end do
#endif

  ! Wait for full completion of receives
  call MPI_WAITALL(countrecv,reqrecv,statuses,info)

  ! Scatter reception array
  do icpu=1,ncpu
    if (reception(icpu,ilevel)%ngrid>0) then
      do j=1,twotondim
        step=(j-1)*reception(icpu,ilevel)%ngrid
        iskip=ncoarse+(j-1)*ngridmax
        do i=1,reception(icpu,ilevel)%ngrid
#ifdef LIGHT_MPI_COMM
          xx(reception(icpu,ilevel)%pcomm%igrid(i)+iskip)=reception(icpu,ilevel)%pcomm%u(i+step,1)
#else
         xx(reception(icpu,ilevel)%igrid(i)+iskip)=reception(icpu,ilevel)%u(i+step,1)
#endif
        end do
      end do
    end if
  end do

  ! Wait for full completion of sends
#ifdef LIGHT_MPI_COMM
  call MPI_WAITALL(emission(ilevel)%nactive,reqsend,statuses,info)

  ! Deallocate temporary communication buffers
  do icpu=1,ncpu
     ncache=reception(icpu,ilevel)%ngrid
     if(ncache>0)deallocate(reception(icpu,ilevel)%pcomm%u)
  end do
  deallocate(emission(ilevel)%u)

  deallocate(reqsend)
#else
  call MPI_WAITALL(countsend,reqsend,statuses,info)
#endif

#endif

111 format('   Entering make_virtual_fine for level ',I2)

end subroutine make_virtual_fine_dp
!################################################################
!################################################################
!################################################################
!################################################################
subroutine make_virtual_fine_int(xx,ilevel)
  use amr_commons
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer,dimension(MPI_STATUS_SIZE,ncpu)::statuses
#endif
  integer::ilevel
  integer,dimension(1:ncoarse+ngridmax*twotondim)::xx
  ! -------------------------------------------------------------------
  ! This routine communicates virtual boundaries among all cpu's.
  ! at level ilevel for any integer array in the AMR grid.
  ! -------------------------------------------------------------------
#ifndef WITHOUTMPI

  integer::icpu,i,j,ncache,iskip,step
#ifdef LIGHT_MPI_COMM
  integer::ind,idx,offset
#endif
  integer::countsend,countrecv
  integer::info,tag=101
#ifdef LIGHT_MPI_COMM
  integer,dimension(ncpu)::reqrecv
  integer,allocatable,dimension(:)::reqsend
#else
  integer,dimension(ncpu)::reqsend,reqrecv
#endif

#endif

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

#ifndef WITHOUTMPI

#ifdef LIGHT_MPI_COMM
  allocate(reqsend(emission(ilevel)%nactive))

  ! Allocate temporary communication buffers
  do icpu=1,ncpu
     ncache=reception(icpu,ilevel)%ngrid
     if(ncache>0) then
        if (.not. associated(reception(icpu,ilevel)%pcomm)) then
           allocate(reception(icpu,ilevel)%pcomm)
        end if
        allocate(reception(icpu,ilevel)%pcomm%f(1:ncache*twotondim,1:1))
     end if
  end do
  allocate(emission(ilevel)%f(1:emission(ilevel)%ngrids_tot*twotondim, 1:1))
#endif

  ! Receive all messages
  countrecv=0
  do icpu=1,ncpu
     ncache=reception(icpu,ilevel)%ngrid
     if(ncache>0) then
       countrecv=countrecv+1
#ifdef LIGHT_MPI_COMM
       call MPI_IRECV(reception(icpu,ilevel)%pcomm%f,ncache*twotondim, &
            & MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
#else
       call MPI_IRECV(reception(icpu,ilevel)%f,ncache*twotondim, &
            & MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
#endif
     end if
  end do

  ! Gather emission array
#ifdef LIGHT_MPI_COMM
  offset=0
  do idx=1,emission(ilevel)%nactive
    do j=1,twotondim
      step=(j-1)*emission(ilevel)%ngrids(idx)
      iskip=ncoarse+(j-1)*ngridmax
      do i=1,emission(ilevel)%ngrids(idx)
        emission(ilevel)%f(offset*twotondim+i+step,1)=xx(emission(ilevel)%igrid(offset+i)+iskip)
      end do
    end do
    offset=offset+emission(ilevel)%ngrids(idx)
  end do
#else
  do icpu=1,ncpu
     if (emission(icpu,ilevel)%ngrid>0) then
       do j=1,twotondim
         step=(j-1)*emission(icpu,ilevel)%ngrid
         iskip=ncoarse+(j-1)*ngridmax
         do i=1,emission(icpu,ilevel)%ngrid
           emission(icpu,ilevel)%f(i+step,1)=xx(emission(icpu,ilevel)%igrid(i)+iskip)
         end do
       end do
     end if
  end do
#endif

  ! Send all messages
#ifdef LIGHT_MPI_COMM
  offset=1
  do idx=1,emission(ilevel)%nactive
    ncache=emission(ilevel)%ngrids(idx)*twotondim
    call MPI_ISEND(emission(ilevel)%f(offset,1),ncache,MPI_INTEGER,&
                   emission(ilevel)%cpuid(idx)-1,tag,MPI_COMM_WORLD,reqsend(idx),info)
    offset=offset+ncache
  end do
#else
  countsend=0
  do icpu=1,ncpu
     ncache=emission(icpu,ilevel)%ngrid
     if(ncache>0) then
       countsend=countsend+1
       call MPI_ISEND(emission(icpu,ilevel)%f,ncache*twotondim, &
            & MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
     end if
  end do
#endif

  ! Wait for full completion of receives
  call MPI_WAITALL(countrecv,reqrecv,statuses,info)

  ! Scatter reception array
  do icpu=1,ncpu
    if (reception(icpu,ilevel)%ngrid>0) then
      do j=1,twotondim
        step=(j-1)*reception(icpu,ilevel)%ngrid
        iskip=ncoarse+(j-1)*ngridmax
        do i=1,reception(icpu,ilevel)%ngrid
#ifdef LIGHT_MPI_COMM
          xx(reception(icpu,ilevel)%pcomm%igrid(i)+iskip)=reception(icpu,ilevel)%pcomm%f(i+step,1)
#else
          xx(reception(icpu,ilevel)%igrid(i)+iskip)=reception(icpu,ilevel)%f(i+step,1)
#endif
        end do
      end do
    end if
  end do

  ! Wait for full completion of sends
#ifdef LIGHT_MPI_COMM
  call MPI_WAITALL(emission(ilevel)%nactive,reqsend,statuses,info)

  ! Deallocate temporary communication buffers
  do icpu=1,ncpu
     ncache=reception(icpu,ilevel)%ngrid
     if(ncache>0)deallocate(reception(icpu,ilevel)%pcomm%f)
  end do
  deallocate(emission(ilevel)%f)

  deallocate(reqsend)
#else
  call MPI_WAITALL(countsend,reqsend,statuses,info)
#endif

#endif

111 format('   Entering make_virtual_fine for level ',I2)

end subroutine make_virtual_fine_int
!################################################################
!################################################################
!################################################################
!################################################################
subroutine make_virtual_reverse_dp(xx,ilevel)
  use amr_commons
  use mpi_mod
  implicit none
  integer::ilevel
  real(dp),dimension(1:ncoarse+ngridmax*twotondim)::xx
  ! -------------------------------------------------------------------
  ! This routine communicates virtual boundaries among all cpu's.
  ! at level ilevel in a reverse way for double precision arrays.
  ! -------------------------------------------------------------------
#ifndef WITHOUTMPI
  integer::icpu,i,j,ncache,iskip,step,icell,ibuf
#ifdef LIGHT_MPI_COMM
  integer::ind,idx,offset
#endif
  integer::countsend,countrecv
  integer::info,tag=101
#ifdef LIGHT_MPI_COMM
  integer,dimension(ncpu)::reqsend
  integer,allocatable,dimension(:)::reqrecv
#else
  integer,dimension(ncpu)::reqsend,reqrecv
#endif
  integer,dimension(MPI_STATUS_SIZE,ncpu)::statuses
  integer::switchlevel=3
#endif

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

#ifndef WITHOUTMPI

#ifdef LIGHT_MPI_COMM
  ! Allocate temporary communication buffers
  do icpu=1,ncpu
     ncache=reception(icpu,ilevel)%ngrid
     if(ncache>0) then
        if (.not. associated(reception(icpu,ilevel)%pcomm)) then
           allocate(reception(icpu,ilevel)%pcomm)
        end if
        allocate(reception(icpu,ilevel)%pcomm%u(1:ncache*twotondim, 1:1))
     end if
  end do
  allocate(emission(ilevel)%u(1:emission(ilevel)%ngrids_tot*twotondim, 1:1))
#endif

  if(ilevel.LE.switchlevel)then

 ! Gather emission array
  do icpu=1,ncpu
    if (reception(icpu,ilevel)%ngrid>0) then
      do j=1,twotondim
        step=(j-1)*reception(icpu,ilevel)%ngrid
        iskip=ncoarse+(j-1)*ngridmax
        do i=1,reception(icpu,ilevel)%ngrid
           ibuf=i+step
#ifdef LIGHT_MPI_COMM
           icell=reception(icpu,ilevel)%pcomm%igrid(i)+iskip
           reception(icpu,ilevel)%pcomm%u(ibuf,1)=xx(icell)
#else
           icell=reception(icpu,ilevel)%igrid(i)+iskip
           reception(icpu,ilevel)%u(ibuf,1)=xx(icell)
#endif
        end do
      end do
   end if
  end do

  ! Receive all messages
#ifdef LIGHT_MPI_COMM
  countrecv=0
  offset=1
  do idx=1,emission(ilevel)%nactive
     if (emission(ilevel)%cpuid(idx) > myid) exit
     ncache=emission(ilevel)%ngrids(idx)*twotondim
     countrecv=countrecv+1
     ! request to send
     call MPI_SEND(countrecv,0, MPI_INTEGER, emission(ilevel)%cpuid(idx)-1,101,MPI_COMM_WORLD,info)
     call MPI_RECV(emission(ilevel)%u(offset,1),ncache,MPI_DOUBLE_PRECISION,&
          & emission(ilevel)%cpuid(idx)-1,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,info)
     offset=offset+ncache
  end do
#else
  countrecv=0
  do icpu=1,myid-1
     ncache=emission(icpu,ilevel)%ngrid
     if(ncache>0) then
        countrecv=countrecv+1
        ! request to send
        call MPI_SEND(countrecv,0, MPI_INTEGER, icpu-1,101,MPI_COMM_WORLD,info)
        call MPI_RECV(emission(icpu,ilevel)%u,ncache*twotondim, &
             & MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,info)
     end if
  end do
#endif

  ! Send all messages
  countsend=0
  do icpu=1,ncpu
     ncache=reception(icpu,ilevel)%ngrid
     if(ncache>0) then
        countsend=countsend+1
        ! wait for request to send
        call MPI_RECV(countrecv,0, MPI_INTEGER, icpu-1,101,MPI_COMM_WORLD, &
             & MPI_STATUS_IGNORE, info)
#ifdef LIGHT_MPI_COMM
        call MPI_SEND(reception(icpu,ilevel)%pcomm%u,ncache*twotondim, &
             & MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,info)
#else
        call MPI_SEND(reception(icpu,ilevel)%u,ncache*twotondim, &
             & MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,info)
#endif
     end if
  end do

  ! Receive all messages
#ifdef LIGHT_MPI_COMM
  countrecv=0
  offset=1
  do idx=1,emission(ilevel)%nactive
     ncache=emission(ilevel)%ngrids(idx)*twotondim
     if (emission(ilevel)%cpuid(idx) > myid) then
         countrecv=countrecv+1
         ! request to send
         call MPI_SEND(countrecv,0, MPI_INTEGER, emission(ilevel)%cpuid(idx)-1,101,MPI_COMM_WORLD,info)
         call MPI_RECV(emission(ilevel)%u(offset,1),ncache,MPI_DOUBLE_PRECISION,&
                       & emission(ilevel)%cpuid(idx)-1,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,info)
     end if
     offset=offset+ncache
  end do
#else
  countrecv=0
  do icpu=myid+1,ncpu
     ncache=emission(icpu,ilevel)%ngrid
     if(ncache>0) then
        countrecv=countrecv+1
        ! request to send
        call MPI_SEND(countrecv,0, MPI_INTEGER, icpu-1,101,MPI_COMM_WORLD,info)
        call MPI_RECV(emission(icpu,ilevel)%u,ncache*twotondim, &
             & MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,info)
     end if
  end do
#endif

  ! Scatter reception array
#ifdef LIGHT_MPI_COMM
  offset=0
  do idx=1,emission(ilevel)%nactive
    do j=1,twotondim
      step=(j-1)*emission(ilevel)%ngrids(idx)
      iskip=ncoarse+(j-1)*ngridmax
      do i=1,emission(ilevel)%ngrids(idx)
        xx(emission(ilevel)%igrid(offset+i)+iskip)= xx(emission(ilevel)%igrid(offset+i)+iskip) &
                                                    & + emission(ilevel)%u(offset*twotondim+i+step,1)
      end do
    end do
    offset=offset+emission(ilevel)%ngrids(idx)
  end do
#else
  do icpu=1,ncpu
     if (emission(icpu,ilevel)%ngrid>0) then
        do j=1,twotondim
           step=(j-1)*emission(icpu,ilevel)%ngrid
           iskip=ncoarse+(j-1)*ngridmax
           do i=1,emission(icpu,ilevel)%ngrid
              xx(emission(icpu,ilevel)%igrid(i)+iskip)= &
                   & xx(emission(icpu,ilevel)%igrid(i)+iskip) + emission(icpu,ilevel)%u(i+step,1)
           end do
        end do
     end if
  end do
#endif

  else

#ifdef LIGHT_MPI_COMM
  allocate(reqrecv(emission(ilevel)%nactive))

  ! Receive all messages
  offset=1
  do idx=1,emission(ilevel)%nactive
    ncache=emission(ilevel)%ngrids(idx)*twotondim
    call MPI_IRECV(emission(ilevel)%u(offset,1),ncache,MPI_DOUBLE_PRECISION,&
                   emission(ilevel)%cpuid(idx)-1,tag,MPI_COMM_WORLD,reqrecv(idx),info)
    offset=offset+ncache
  end do
#else
  ! Receive all messages
  countrecv=0
  do icpu=1,ncpu
     ncache=emission(icpu,ilevel)%ngrid
     if(ncache>0) then
        countrecv=countrecv+1
        call MPI_IRECV(emission(icpu,ilevel)%u,ncache*twotondim, &
             & MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
     end if
  end do
#endif

  ! Gather emission array
  do icpu=1,ncpu
    if (reception(icpu,ilevel)%ngrid>0) then
      do j=1,twotondim
        step=(j-1)*reception(icpu,ilevel)%ngrid
        iskip=ncoarse+(j-1)*ngridmax
        do i=1,reception(icpu,ilevel)%ngrid
#ifdef LIGHT_MPI_COMM
          reception(icpu,ilevel)%pcomm%u(i+step,1)=xx(reception(icpu,ilevel)%pcomm%igrid(i)+iskip)
#else
          reception(icpu,ilevel)%u(i+step,1)=xx(reception(icpu,ilevel)%igrid(i)+iskip)
#endif
        end do
      end do
   end if
  end do

  ! Send all messages
  countsend=0
  do icpu=1,ncpu
     ncache=reception(icpu,ilevel)%ngrid
     if(ncache>0) then
       countsend=countsend+1
#ifdef LIGHT_MPI_COMM
        call MPI_ISEND(reception(icpu,ilevel)%pcomm%u,ncache*twotondim, &
             & MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
#else
        call MPI_ISEND(reception(icpu,ilevel)%u,ncache*twotondim, &
             & MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
#endif
     end if
  end do

  ! Wait for full completion of receives
#ifdef LIGHT_MPI_COMM
  call MPI_WAITALL(emission(ilevel)%nactive,reqrecv,statuses,info)
#else
  call MPI_WAITALL(countrecv,reqrecv,statuses,info)
#endif

  ! Scatter reception array
#ifdef LIGHT_MPI_COMM
  offset=0
  do idx=1,emission(ilevel)%nactive
    do j=1,twotondim
      step=(j-1)*emission(ilevel)%ngrids(idx)
      iskip=ncoarse+(j-1)*ngridmax
      do i=1,emission(ilevel)%ngrids(idx)
        xx(emission(ilevel)%igrid(offset+i)+iskip)= xx(emission(ilevel)%igrid(offset+i)+iskip) &
                                                    +emission(ilevel)%u(offset*twotondim+i+step,1)
      end do
    end do
    offset=offset+emission(ilevel)%ngrids(idx)
  end do
#else
  do icpu=1,ncpu
     if (emission(icpu,ilevel)%ngrid>0) then
       do j=1,twotondim
         step=(j-1)*emission(icpu,ilevel)%ngrid
         iskip=ncoarse+(j-1)*ngridmax
         do i=1,emission(icpu,ilevel)%ngrid
            xx(emission(icpu,ilevel)%igrid(i)+iskip)= &
               & xx(emission(icpu,ilevel)%igrid(i)+iskip) + emission(icpu,ilevel)%u(i+step,1)
         end do
       end do
     end if
   end do
#endif

  ! Wait for full completion of sends
  call MPI_WAITALL(countsend,reqsend,statuses,info)

#ifdef LIGHT_MPI_COMM
  deallocate(reqrecv)
#endif

  endif

#ifdef LIGHT_MPI_COMM
  ! Deallocate temporary communication buffers
  do icpu=1,ncpu
     ncache=reception(icpu,ilevel)%ngrid
     if(ncache>0)deallocate(reception(icpu,ilevel)%pcomm%u)
  end do
  deallocate(emission(ilevel)%u)
#endif

#endif

111 format('   Entering make_virtual_reverse for level ',I2)

end subroutine make_virtual_reverse_dp
!################################################################
!################################################################
!################################################################
!################################################################
subroutine make_virtual_reverse_int(xx,ilevel)
  use amr_commons
  use mpi_mod
  implicit none
  integer::ilevel
  integer,dimension(1:ncoarse+ngridmax*twotondim)::xx
  ! -------------------------------------------------------------------
  ! This routine communicates virtual boundaries among all cpu's.
  ! at level ilevel in a reverse way for integer arrays.
  ! -------------------------------------------------------------------
#ifndef WITHOUTMPI

  integer::icpu,i,j,ncache,iskip,step,icell,ibuf
#ifdef LIGHT_MPI_COMM
  integer::ind,idx,offset
#endif
  integer::countsend,countrecv
  integer::info,tag=101
#ifdef LIGHT_MPI_COMM
  integer,dimension(ncpu)::reqsend
  integer,allocatable,dimension(:)::reqrecv
#else
   integer,dimension(ncpu)::reqsend,reqrecv
#endif
  integer,dimension(MPI_STATUS_SIZE,ncpu)::statuses
  integer::switchlevel=3

#endif

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

#ifndef WITHOUTMPI

#ifdef LIGHT_MPI_COMM
  ! Allocate temporary communication buffers
  do icpu=1,ncpu
     ncache=reception(icpu,ilevel)%ngrid
     if(ncache>0) then
        if (.not. associated(reception(icpu,ilevel)%pcomm)) then
           allocate(reception(icpu,ilevel)%pcomm)
        end if
        allocate(reception(icpu,ilevel)%pcomm%f(1:ncache*twotondim, 1:1))
     end if
  end do
  allocate(emission(ilevel)%f(1:emission(ilevel)%ngrids_tot*twotondim, 1:1))
#endif

  if(ilevel.le.switchlevel) then

  ! Gather emission array
  do icpu=1,ncpu
     if (reception(icpu,ilevel)%ngrid>0) then
        do j=1,twotondim
           step=(j-1)*reception(icpu,ilevel)%ngrid
           iskip=ncoarse+(j-1)*ngridmax
           do i=1,reception(icpu,ilevel)%ngrid
              ibuf=i+step
#ifdef LIGHT_MPI_COMM
              icell=reception(icpu,ilevel)%pcomm%igrid(i)+iskip
              reception(icpu,ilevel)%pcomm%f(ibuf,1)=xx(icell)
#else
              icell=reception(icpu,ilevel)%igrid(i)+iskip
              reception(icpu,ilevel)%f(ibuf,1)=xx(icell)
#endif
           end do
        end do
     end if
  end do

  ! Receive all messages
#ifdef LIGHT_MPI_COMM
  countrecv=0
  offset=1
  do idx=1,emission(ilevel)%nactive
     if (emission(ilevel)%cpuid(idx) > myid) exit
     ncache=emission(ilevel)%ngrids(idx)*twotondim
     countrecv=countrecv+1
     ! request to send
     call MPI_SEND(countrecv,0, MPI_INTEGER, emission(ilevel)%cpuid(idx)-1,101,MPI_COMM_WORLD,info)
     call MPI_RECV(emission(ilevel)%f(offset,1),ncache, MPI_INTEGER, &
          & emission(ilevel)%cpuid(idx)-1,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,info)
     offset=offset+ncache
  end do
#else
  countrecv=0
  do icpu=1,myid-1
     ncache=emission(icpu,ilevel)%ngrid
     if(ncache>0) then
        countrecv=countrecv+1
        ! request to send
        call MPI_SEND(countrecv,0, MPI_INTEGER, icpu-1,101,MPI_COMM_WORLD,info)
        call MPI_RECV(emission(icpu,ilevel)%f,ncache*twotondim, &
             & MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,info)
     end if
  end do
#endif

  ! Send all messages
  countsend=0
  do icpu=1,ncpu
     ncache=reception(icpu,ilevel)%ngrid
     if(ncache>0) then
        countsend=countsend+1
        ! wait for request to send
        call MPI_RECV(countrecv,0, MPI_INTEGER, icpu-1,101,MPI_COMM_WORLD, &
             & MPI_STATUS_IGNORE, info)
#ifdef LIGHT_MPI_COMM
        call MPI_SEND(reception(icpu,ilevel)%pcomm%f,ncache*twotondim, &
                      & MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,info)
#else
        call MPI_SEND(reception(icpu,ilevel)%f,ncache*twotondim, &
                      & MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,info)
#endif
     end if
  end do

  ! Receive all messages
#ifdef LIGHT_MPI_COMM
  countrecv=0
  offset = 1
  do idx=1,emission(ilevel)%nactive
     ncache=emission(ilevel)%ngrids(idx)*twotondim
     if (emission(ilevel)%cpuid(idx) > myid) then
        countrecv=countrecv+1
        ! request to send
        call MPI_SEND(countrecv,0, MPI_INTEGER, emission(ilevel)%cpuid(idx)-1,101,MPI_COMM_WORLD,info)
        call MPI_RECV(emission(ilevel)%f(offset,1),ncache, MPI_INTEGER, &
                      & emission(ilevel)%cpuid(idx)-1,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,info)
     end if
     offset=offset+ncache
  end do
#else
  countrecv=0
  do icpu=myid+1,ncpu
     ncache=emission(icpu,ilevel)%ngrid
     if(ncache>0) then
        countrecv=countrecv+1
        ! request to send
        call MPI_SEND(countrecv,0, MPI_INTEGER, icpu-1,101,MPI_COMM_WORLD,info)
        call MPI_RECV(emission(icpu,ilevel)%f,ncache*twotondim, &
                      & MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,info)
     end if
  end do
#endif

  ! Scatter reception array
#ifdef LIGHT_MPI_COMM
  offset=0
  do idx=1,emission(ilevel)%nactive
     do j=1,twotondim
        step=(j-1)*emission(ilevel)%ngrids(idx)
        iskip=ncoarse+(j-1)*ngridmax
        do i=1,emission(ilevel)%ngrids(idx)
           xx(emission(ilevel)%igrid(offset+i)+iskip)= xx(emission(ilevel)%igrid(offset+i)+iskip) &
                                                       & + emission(ilevel)%f(offset*twotondim+i+step,1)
        end do
     end do
     offset=offset+emission(ilevel)%ngrids(idx)
  end do
#else
  do icpu=1,ncpu
     if (emission(icpu,ilevel)%ngrid>0) then
        do j=1,twotondim
           step=(j-1)*emission(icpu,ilevel)%ngrid
           iskip=ncoarse+(j-1)*ngridmax
           do i=1,emission(icpu,ilevel)%ngrid
              xx(emission(icpu,ilevel)%igrid(i)+iskip)= xx(emission(icpu,ilevel)%igrid(i)+iskip) &
                                                        & + emission(icpu,ilevel)%f(i+step,1)
           end do
        end do
     end if
  end do
#endif

  else

#ifdef LIGHT_MPI_COMM
  allocate(reqrecv(emission(ilevel)%nactive))

  ! Receive all messages
  offset=1
  do idx=1,emission(ilevel)%nactive
    ncache=emission(ilevel)%ngrids(idx)*twotondim
    call MPI_IRECV(emission(ilevel)%f(offset,1),ncache,MPI_INTEGER,&
                   emission(ilevel)%cpuid(idx)-1,tag,MPI_COMM_WORLD,reqrecv(idx),info)
    offset=offset+ncache
  end do
#else
  ! Receive all messages
  countrecv=0
  do icpu=1,ncpu
     ncache=emission(icpu,ilevel)%ngrid
     if(ncache>0) then
        countrecv=countrecv+1
        call MPI_IRECV(emission(icpu,ilevel)%f,ncache*twotondim, &
             & MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
     end if
  end do
#endif

  ! Gather emission array
  do icpu=1,ncpu
     if (reception(icpu,ilevel)%ngrid>0) then
        do j=1,twotondim
           step=(j-1)*reception(icpu,ilevel)%ngrid
           iskip=ncoarse+(j-1)*ngridmax
           do i=1,reception(icpu,ilevel)%ngrid
#ifdef LIGHT_MPI_COMM
              reception(icpu,ilevel)%pcomm%f(i+step,1)=xx(reception(icpu,ilevel)%pcomm%igrid(i)+iskip)
#else
              reception(icpu,ilevel)%f(i+step,1)=xx(reception(icpu,ilevel)%igrid(i)+iskip)
#endif
           end do
        end do
     end if
  end do

  ! Send all messages
  countsend=0
  do icpu=1,ncpu
     ncache=reception(icpu,ilevel)%ngrid
     if(ncache>0) then
       countsend=countsend+1
#ifdef LIGHT_MPI_COMM
        call MPI_ISEND(reception(icpu,ilevel)%pcomm%f,ncache*twotondim, &
             & MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
#else
        call MPI_ISEND(reception(icpu,ilevel)%f,ncache*twotondim, &
             & MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
#endif
     end if
  end do

  ! Wait for full completion of receives
#ifdef LIGHT_MPI_COMM
  call MPI_WAITALL(emission(ilevel)%nactive,reqrecv,statuses,info)
#else
  call MPI_WAITALL(countrecv,reqrecv,statuses,info)
#endif

  ! Scatter reception array
#ifdef LIGHT_MPI_COMM
  offset=0
  do idx=1,emission(ilevel)%nactive
    do j=1,twotondim
      step=(j-1)*emission(ilevel)%ngrids(idx)
      iskip=ncoarse+(j-1)*ngridmax
      do i=1,emission(ilevel)%ngrids(idx)
        xx(emission(ilevel)%igrid(offset+i)+iskip)= xx(emission(ilevel)%igrid(offset+i)+iskip) &
                                                    & + emission(ilevel)%f(offset*twotondim+i+step,1)
      end do
    end do
    offset=offset+emission(ilevel)%ngrids(idx)
  end do
#else
  do icpu=1,ncpu
     if (emission(icpu,ilevel)%ngrid>0) then
        do j=1,twotondim
           step=(j-1)*emission(icpu,ilevel)%ngrid
           iskip=ncoarse+(j-1)*ngridmax
           do i=1,emission(icpu,ilevel)%ngrid
              xx(emission(icpu,ilevel)%igrid(i)+iskip)= xx(emission(icpu,ilevel)%igrid(i)+iskip) &
                                                        & + emission(icpu,ilevel)%f(i+step,1)
           end do
        end do
     end if
  end do
#endif

  ! Wait for full completion of sends
  call MPI_WAITALL(countsend,reqsend,statuses,info)

#ifdef LIGHT_MPI_COMM
  deallocate(reqrecv)
#endif

  endif

#ifdef LIGHT_MPI_COMM
  ! Deallocate temporary communication buffers
  do icpu=1,ncpu
     ncache=reception(icpu,ilevel)%ngrid
     if(ncache>0)deallocate(reception(icpu,ilevel)%pcomm%f)
  end do
  deallocate(emission(ilevel)%f)
#endif


#endif

111 format('   Entering make_virtual_reverse for level ',I2)

end subroutine make_virtual_reverse_int
!################################################################
!################################################################
!################################################################
!################################################################
subroutine build_comm(ilevel)
  use amr_commons
  use poisson_commons, only: lookup_mg
  use mpi_mod
  implicit none
  integer::ilevel
  ! -------------------------------------------------------------------
  ! This routine builds the communication structure for level ilevel.
  ! Array flag2 is used as temporary work space.
  ! -------------------------------------------------------------------
  integer::icpu,ibound
#ifdef LIGHT_MPI_COMM
  integer::idx,offset, ngrids
#endif
  integer::ncache,ind,iskip
  integer::i,j,k,nxny
  integer::igrid,jgrid,ngrid
#ifndef WITHOUTMPI
  integer::info,tag=101
#ifdef LIGHT_MPI_COMM
  integer,dimension(ncpu)::reqsend
  integer,allocatable,dimension(:)::reqrecv
#else
   integer,dimension(ncpu)::reqsend,reqrecv
#endif
  integer,dimension(ncpu)::sendbuf,recvbuf
  integer,dimension(MPI_STATUS_SIZE,ncpu)::statuses
  integer::countsend,countrecv
#endif
  integer,dimension(1:nvector),save::ind_grid,ind_cell

  if(verbose)write(*,111)ilevel
  nxny=nx*ny

  !----------------------------------------------------------------
  ! Compute grids global adress using flag2 array at level ilevel-1
  !----------------------------------------------------------------
  if(ilevel==1)then
     do k=kcoarse_min,kcoarse_max
     do j=jcoarse_min,jcoarse_max
     do i=icoarse_min,icoarse_max
        ind=1+i+j*nx+k*nxny
        if(cpu_map(ind)==myid)then
           flag2(ind)=son(ind)
        else
           flag2(ind)=0
        end if
     end do
     end do
     end do
     call make_virtual_coarse_int(flag2(1))
  else
     ! Initialize flag2 to local adress for cpu map = myid cells
     ncache=active(ilevel-1)%ngrid
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel-1)%igrid(igrid+i-1)
        end do
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do
           do i=1,ngrid
              if(cpu_map(ind_cell(i))==myid)then
                 flag2(ind_cell(i))=son(ind_cell(i))
              else
                 flag2(ind_cell(i))=0
              end if
           end do
        end do
     end do
     do icpu=1,ncpu
        ncache=reception(icpu,ilevel-1)%ngrid
        do igrid=1,ncache,nvector
           ngrid=MIN(nvector,ncache-igrid+1)
           do i=1,ngrid
#ifdef LIGHT_MPI_COMM
              ind_grid(i)=reception(icpu,ilevel-1)%pcomm%igrid(igrid+i-1)
#else
              ind_grid(i)=reception(icpu,ilevel-1)%igrid(igrid+i-1)
#endif
           end do
           do ind=1,twotondim
              iskip=ncoarse+(ind-1)*ngridmax
              do i=1,ngrid
                 ind_cell(i)=iskip+ind_grid(i)
              end do
              do i=1,ngrid
                 if(cpu_map(ind_cell(i))==myid)then
                    flag2(ind_cell(i))=son(ind_cell(i))
                 else
                    flag2(ind_cell(i))=0
                 end if
              end do
           end do
        end do
     end do
     call make_virtual_reverse_int(flag2(1),ilevel-1)
     call make_virtual_fine_int   (flag2(1),ilevel-1)
  end if

  !--------------------------------------------------------
  ! Compute number and index of active grid at level ilevel
  !--------------------------------------------------------
  ncache=numbl(myid,ilevel)
  ! Reset old communicator
  if(active(ilevel)%ngrid>0)then
     active(ilevel)%ngrid=0
     deallocate(active(ilevel)%igrid)
  end if
  if(ncache>0)then
     ! Allocate grid index to new communicator
     active(ilevel)%ngrid=ncache
     allocate(active(ilevel)%igrid(1:ncache))
     ! Gather all grids
     igrid=headl(myid,ilevel)
     do jgrid=1,numbl(myid,ilevel)
        active(ilevel)%igrid(jgrid)=igrid
        igrid=next(igrid)
     end do
  end if
  ! Fill up lookup_mg for active
  if(poisson)then
     igrid=headl(myid,ilevel)
     do jgrid=1,numbl(myid,ilevel)
        lookup_mg(igrid)=0
        igrid=next(igrid)
     end do
  end if

  !----------------------------------------------------
  ! Compute number and index of physical boundary grids
  !----------------------------------------------------
  do ibound=1,nboundary
     ncache=numbb(ibound,ilevel)
     ! Reset old communicator
     if(boundary(ibound,ilevel)%ngrid>0)then
        boundary(ibound,ilevel)%ngrid=0
        deallocate(boundary(ibound,ilevel)%igrid)
     end if
     if(ncache>0)then
        ! Allocate grid index to new communicator
        boundary(ibound,ilevel)%ngrid=ncache
        allocate(boundary(ibound,ilevel)%igrid(1:ncache))
        ! Gather all grids
        igrid=headb(ibound,ilevel)
        do jgrid=1,numbb(ibound,ilevel)
           boundary(ibound,ilevel)%igrid(jgrid)=igrid
           igrid=next(igrid)
        end do
     end if
  end do

  !----------------------------------------------------
  ! Compute number and index of virtual boundary grids
  !----------------------------------------------------
#ifndef WITHOUTMPI

#ifdef LIGHT_MPI_COMM
   if (emission(ilevel)%nactive>0) then
     emission(ilevel)%nactive=0
     emission(ilevel)%ngrids_tot=0
     deallocate(emission(ilevel)%cpuid)
     deallocate(emission(ilevel)%ngrids)
     deallocate(emission(ilevel)%igrid)
   end if
#endif
   do icpu=1,ncpu
      ncache=0
      if(icpu.ne.myid)ncache=numbl(icpu,ilevel)
      ! Reset old communicators
#ifndef LIGHT_MPI_COMM
      if(emission(icpu,ilevel)%ngrid>0)then
         emission(icpu,ilevel)%ngrid=0
         deallocate(emission(icpu,ilevel)%igrid)
         deallocate(emission(icpu,ilevel)%u)
         deallocate(emission(icpu,ilevel)%f)
      end if
#endif
      if(reception(icpu,ilevel)%ngrid>0)then
         reception(icpu,ilevel)%ngrid=0
#ifdef LIGHT_MPI_COMM
         deallocate(reception(icpu,ilevel)%pcomm%igrid)
#else
         deallocate(reception(icpu,ilevel)%igrid)
         deallocate(reception(icpu,ilevel)%u)
         deallocate(reception(icpu,ilevel)%f)
#endif
      end if
      if(ncache>0)then
         ! Allocate grid index to new communicator
         reception(icpu,ilevel)%ngrid=ncache
#ifdef LIGHT_MPI_COMM
         if (.not. associated(reception(icpu,ilevel)%pcomm)) then
            allocate(reception(icpu,ilevel)%pcomm)
         end if
         allocate(reception(icpu,ilevel)%pcomm%igrid(1:ncache))
#else
         allocate(reception(icpu,ilevel)%igrid(1:ncache))
#endif
         ! Gather all grids
         igrid=headl(icpu,ilevel)
         do jgrid=1,numbl(icpu,ilevel)
#ifdef LIGHT_MPI_COMM
            reception(icpu,ilevel)%pcomm%igrid(jgrid)=igrid
#else
            reception(icpu,ilevel)%igrid(jgrid)=igrid
#endif
            igrid=next(igrid)
         end do
         ! Allocate temporary communication buffer
#ifdef LIGHT_MPI_COMM
         allocate(reception(icpu,ilevel)%pcomm%f(1:ncache, 1:1))
         do i=1,ncache
            reception(icpu,ilevel)%pcomm%f(i,1) = &
            & flag2(father(reception(icpu,ilevel)%pcomm%igrid(i)))
         end do
#else
         allocate(reception(icpu,ilevel)%f(1:ncache,1:1))
         do i=1,ncache
            reception(icpu,ilevel)%f(i,1) = &
            & flag2(father(reception(icpu,ilevel)%igrid(i)))
         end do
#endif

         ! Fill up lookup_mg for reception
         if(poisson)then
            do i=1,ncache
#ifdef LIGHT_MPI_COMM
               lookup_mg(reception(icpu,ilevel)%pcomm%igrid(i))= -reception(icpu,ilevel)%pcomm%f(i,1)
#else
               lookup_mg(reception(icpu,ilevel)%igrid(i))= -reception(icpu,ilevel)%f(i,1)
#endif
            end do
         end if
      end if
      sendbuf(icpu)=reception(icpu,ilevel)%ngrid
   end do

  !--------------------------------------------------------
  ! Communicate virtual grid number and index to parent cpu
  !--------------------------------------------------------
  call MPI_ALLTOALL(sendbuf,1,MPI_INTEGER,recvbuf,1,MPI_INTEGER,MPI_COMM_WORLD,info)

  ! Allocate grid index
#ifdef LIGHT_MPI_COMM
  ! Count entries
  emission(ilevel)%nactive=0
  emission(ilevel)%ngrids_tot=0
  do icpu=1,ncpu
    if(icpu==myid) cycle
    ngrids=recvbuf(icpu)
    if(ngrids>0) then
      emission(ilevel)%nactive=emission(ilevel)%nactive+1
      emission(ilevel)%ngrids_tot=emission(ilevel)%ngrids_tot+ngrids
    end if
  end do
  ! Allocate structures
  if(emission(ilevel)%nactive>0)then
    allocate(emission(ilevel)%cpuid (emission(ilevel)%nactive))
    allocate(emission(ilevel)%ngrids(emission(ilevel)%nactive))
    allocate(emission(ilevel)%igrid (emission(ilevel)%ngrids_tot))
  end if
  allocate(reqrecv(emission(ilevel)%nactive))
#else
  ! Allocate grid index
  do icpu=1,ncpu
   emission(icpu,ilevel)%ngrid=recvbuf(icpu)
   ncache=emission(icpu,ilevel)%ngrid
   if(ncache>0)allocate(emission(icpu,ilevel)%igrid(1:ncache))
end do
#endif

  ! Receive grid list
#ifdef LIGHT_MPI_COMM
  idx=1
  offset=1
  do icpu=1,ncpu
    if(icpu==myid) cycle
    ngrids=recvbuf(icpu)
    if(ngrids>0) then
      emission(ilevel)%cpuid(idx)=icpu
      emission(ilevel)%ngrids(idx)=ngrids
      call MPI_IRECV(emission(ilevel)%igrid(offset),ngrids, &
           & MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqrecv(idx),info)
      idx=idx+1
      offset=offset+ngrids
    end if
    if(idx>emission(ilevel)%nactive) exit
  end do
#else
  countrecv=0
  do icpu=1,ncpu
     ncache=emission(icpu,ilevel)%ngrid
     if(ncache>0) then
        countrecv=countrecv+1
        call MPI_IRECV(emission(icpu,ilevel)%igrid,ncache, &
             & MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
     end if
  end do
#endif

  ! Send global index
  countsend=0
  do icpu=1,ncpu
     ncache=reception(icpu,ilevel)%ngrid
     if(ncache>0) then
        countsend=countsend+1
#ifdef LIGHT_MPI_COMM
        call MPI_ISEND(reception(icpu,ilevel)%pcomm%f,ncache, &
             & MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
#else
        call MPI_ISEND(reception(icpu,ilevel)%f,ncache, &
             & MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
#endif
     end if
  end do

  ! Wait for full completion of sends
  call MPI_WAITALL(countsend,reqsend,statuses,info)

  ! Wait for full completion of receives
#ifdef LIGHT_MPI_COMM
  call MPI_WAITALL(emission(ilevel)%nactive,reqrecv,statuses,info)
#else
  call MPI_WAITALL(countrecv,reqrecv,statuses,info)
#endif

  ! Deallocate temporary communication buffers
  do icpu=1,ncpu
     ncache=reception(icpu,ilevel)%ngrid
#ifdef LIGHT_MPI_COMM
     if(ncache>0)deallocate(reception(icpu,ilevel)%pcomm%f)
#else
     if(ncache>0)deallocate(reception(icpu,ilevel)%f)
#endif
  end do

#ifdef LIGHT_MPI_COMM
  deallocate(reqrecv)
#else
  ! Allocate temporary communication buffers
  do icpu=1,ncpu
    ncache=emission(icpu,ilevel)%ngrid
    if(ncache>0)then
      allocate(emission(icpu,ilevel)%u(1:ncache*twotondim,1:1))
      allocate(emission(icpu,ilevel)%f(1:ncache*twotondim,1:1))
    endif
    ncache=reception(icpu,ilevel)%ngrid
    if(ncache>0)then
      allocate(reception(icpu,ilevel)%u(1:ncache*twotondim,1:1))
      allocate(reception(icpu,ilevel)%f(1:ncache*twotondim,1:1))
    endif
  end do
#endif

#endif

111 format('   Entering build_comm for level ',I2)

end subroutine build_comm
!################################################################
!################################################################
!################################################################
!################################################################
