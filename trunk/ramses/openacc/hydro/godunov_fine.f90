!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine godunov_fine_acc(ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  use acc_commons
#if defined(_OPENACC) 
!|| defined(USE_ACC_VERSION)
  use device_memory
#endif
  implicit none
#ifndef WITHOUTMPI
#ifndef INCLUDEOK
#define INCLUDEOK
  include 'mpif.h'
#endif
#endif  
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine is a wrapper to the second order Godunov solver.
  ! Patches of size (8x8x8), (4x4x4) and (2x2x2) are gathered from level ilevel 
  ! and sent to the hydro solver.
  ! On entry, hydro variables are gathered from array uold.
  ! On exit, unew has been updated. 
  !--------------------------------------------------------------------------
  integer::i,j,ivar,nx_ok,nx_cell,ngrid_ok,ilev
  integer :: ncache,ngrid,levelup,iskip,tot_cubes
  integer,dimension(1:nlevelmax)::ngroup
  integer,allocatable,dimension(:)::isort
  integer ::nrefinedu
  integer ::nrefinedg
  integer ::npatches
  integer::info
  real(dp) :: dx, dx_loc
#ifdef MPROFILING
#ifndef INCLUDEOK
#define INCLUDEOK
  include 'mpif.h'
#endif
  integer::info
  real(kind=8)::tt1,tt2
  tt1=MPI_WTIME(info)
#endif
#undef INCLUDEOK
    
#if NDIM==1
  print *,"NDIM=1 is not supported in the OpenACC configuration."
#ifndef WITHOUTMPI
  call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
  stop
#endif
#endif
#if NDIM==2
  print *,"NDIM=2 is not supported in the OpenACC configuration."
#ifndef WITHOUTMPI
  call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
  stop
#endif
#endif 
  
  
  if(numbtot(1,ilevel)==0)return
  if(static)return
  if(verbose)write(*,111)ilevel
  
  ncache=active(ilevel)%ngrid
  if(ncache==0) return
  
  levelup=max(ilevel-3,1)
  allocate(isort(ncache)) 
  !$acc data create(isort)
  
  ! Sort grids in large patches.
  call sort_group_grid(isort,ilevel,levelup,ngroup,ncache)
    
  ! Calculate the maximal number of cells and allocate enought memory to store the patches.
  nrefinedu = 0
  nrefinedg = 0
  do ilev=levelup,ilevel-1
     nx_ok=2**(ilevel-1-ilev)
     ngrid_ok=nx_ok**ndim
     nx_cell=2*nx_ok 
     npatches = ngroup(ilev)/ngrid_ok

     nrefinedu = max(nrefinedu, ncube*(nx_cell+4)**3*nvar)
     nrefinedg = max(nrefinedg, ncube*(nx_cell+4)**3*ndim)
  enddo

  ! Allocate auxiliary arrays
  allocate(uloc_tot(nrefinedu))
  allocate(gloc_tot(nrefinedg)) 

  ! Mesh spacing at that level
  dx=0.5D0**ilevel
  dx_loc=boxlen*dx/dble(icoarse_max-icoarse_min+1)
  dx_acc=dx_loc
  dtdx = dtnew(ilevel)/dx_acc
  dtdy = dtnew(ilevel)/dx_acc
  dtdz = dtnew(ilevel)/dx_acc
      
!$acc data create(uloc_tot, gloc_tot) &
!$acc copyin(dtdx,dtdy,dtdz)
  
#ifdef _OPENACC
  ! Print available gpu memory. Here we call a cuda function.
  if(myid==1) call print_dev_mem
#endif
  
#ifndef _OPENACC
  ! This is nedded only when compiled without OpenACC.
  if(pressure_fix) then
     do j=1,size(divu_acc)
        divu_acc(j) = 0.0d0
        enew_acc(j) = 0.0d0
     end do
  end if
#endif
  
  ! Builds and solves the patches. Before (8x8x8), then (4x4x4) and (2x2x2)
  iskip=1
  do ilev=levelup,ilevel-1
     nx_ok=2**(ilevel-1-ilev)
     ngrid_ok=nx_ok**ndim
     nx_cell=2*nx_ok          ! Number of cells in the patch (8, 4 and 2)
     if(ngroup(ilev)>0)then
        tot_cubes = ngroup(ilev)/ngrid_ok  ! Number of patches of the current ilev
        call solve_level(ilevel,ilev,tot_cubes,iskip,ngroup,ngrid_ok,isort,nx_cell,dx,dx_loc)
     endif
  end do
    
!$acc end data
!$acc end data

#if defined(_OPENACC) 
   call update_divu_and_enew(ilevel)
#else
  if(pressure_fix) then
     ! If we are running on the cpu we dont need to update unew but just divu and enew.
     do j=1,size(divu)
        divu(j) = divu(j) + divu_acc(j)
        enew(j) = enew(j) + enew_acc(j)
     end do
  end if
#endif

  ! Deallocate local arrays
  deallocate(isort,uloc_tot,gloc_tot)

#ifdef MPROFILING
  tt2=MPI_WTIME(info)
  acc_t_godunov_fine = acc_t_godunov_fine + (tt2-tt1)
#endif
  
111 format('   Entering godunov_fine for level ',i2)
999 format(' Level',I3,' found ',I6,' groups of size ',I2,'')
end subroutine godunov_fine_acc
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine solve_level(ilevel,ilev,tot_cubes,iskip,ngroup,ngrid_ok,isort,nxp,dx,dx_loc)
   use acc_commons
   use amr_commons
   use hydro_parameters
   implicit none
   integer::ilevel,ilev,tot_cubes,iskip,ngrid_ok,nxp
   integer,dimension(active(ilevel)%ngrid)::isort
   integer,dimension(1:nlevelmax)::ngroup
   real(dp)::dx,dx_loc
   !--------------------------------------------------------------------------
   ! This routine builds tot_cubes patches of size nxp**3 and solves them in 
   ! chunks of ncube patches.
   ! The ncube parameter has a relevant impact on the gpu performance.
   !--------------------------------------------------------------------------
   real(dp)::dx_box
   integer::nxp4,nc1,nc2,i,j
   logical, dimension(-1:nxp+2,-1:nxp+2,-1:nxp+2,1:ncube)::ok
   real(dp),dimension(1:ndim,1:ncube)::box_xmin
   real(dp),dimension(1:nxp+1,1:nxp+1,1:nxp+1,1:nvar,1:ndim,1:ncube)::flux
   real(dp),dimension(1:nxp+1,1:nxp+1,1:nxp+1,1:2   ,1:ndim,1:ncube)::tmp
   
   ! Only for fill_hydro_grids_interpol_2 and 3
   integer, dimension(1:(nxp/2+2)**3,1:ncube)::index_list,grid_list,bindex
   integer, dimension(1:ncube)::z
   integer::maxz
   
   allocate(qin(-1:nxp+2,-1:nxp+2,-1:nxp+2,1:nvar,1:ncube), &
          & cin(-1:nxp+2,-1:nxp+2,-1:nxp+2,1:ncube), &
          & dq(-1:nxp+2,-1:nxp+2,-1:nxp+2,1:nvar,1:ndim,1:ncube), &
          & qm(-1:nxp+2,-1:nxp+2,-1:nxp+2,1:nvar,1:ndim,1:ncube), &
          & qp(-1:nxp+2,-1:nxp+2,-1:nxp+2,1:nvar,1:ndim,1:ncube) )
               
   nxp4 = nxp+4   ! number of cells + boundary
   dx_box = nxp*dx
   
   !$acc data create(ok,box_xmin) &
   !$acc copyin(dtnew(ilevel),nxp,ilevel,dx,dx_box,dx_loc,nxp4) &
   !$acc present(uloc_tot, gloc_tot) &
   !$acc create(index_list,grid_list,bindex,z,maxz) &
   !$acc create(flux,tmp,qin,cin,dq,qm,qp)
   
   ! This loop builds and solves chunks of ncube patches.
   ! ncube is the equivalent of nvector in godunov_fine.
   ! Every patch has a number nc
   nc2=0
   do i=iskip,iskip+ngroup(ilev)-1,ngrid_ok*ncube
     nc1=nc2+1 ! number of the first patch in the chunk
     nc2 = nc2 + min(ncube,(iskip+ngroup(ilev)-i)/ngrid_ok)  ! number of the last patch in the chunk
         
#if defined(_OPENACC) || defined(USE_ACC_VERSION) 
     ! Three different routines are available to build the patches
     
     !call fill_hydro_grids_interpol_1(nxp,ilevel,dx,dx_box,nxp4,nc1,nc2,i,ngrid_ok,isort,ok(:,:,:,1:nc2-nc1+1),box_xmin(:,1:nc2-nc1+1))
     call fill_hydro_grids_interpol_2(nxp,ilevel,dx,dx_box,nxp4,nc1,nc2,i,ngrid_ok,isort,ok(:,:,:,1:nc2-nc1+1),box_xmin(:,1:nc2-nc1+1),index_list,grid_list,bindex,z,maxz)
     !call fill_hydro_grids_interpol_3(nxp,ilevel,dx,dx_box,nxp4,nc1,nc2,i,ngrid_ok,isort,ok(:,:,:,1:nc2-nc1+1),box_xmin(:,1:nc2-nc1+1),index_list,grid_list,bindex,z,maxz)
     
     ! One routine that solves the patches and updates unew at ilevel and ilevel-1
     call solve_hydro_grids_upd_flux(nxp,ilevel,dx,dx_box,dx_loc,nxp4,nc1,nc2,ok(:,:,:,1:nc2-nc1+1), &
                                    & box_xmin(:,1:nc2-nc1+1),flux(:,:,:,:,:,1:nc2-nc1+1), &
                                    & tmp(:,:,:,:,:,1:nc2-nc1+1))
#endif
   end do
   iskip=iskip+ngroup(ilev)
   
   !$acc end data
  
   deallocate(qin,cin,dq,qm,qp)
   
end subroutine solve_level
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine sort_group_grid(isort_fin,ilevel,levelup,ngroup,ncache)
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel,levelup,ncache
  integer,dimension(1:ncache)::isort_fin
  integer,dimension(1:nlevelmax)::ngroup
  !--------------------------------------------------------------------------
  ! This routine uses hilbert keys in order to build the patches. 
  ! First it tries to build patches 8x8x8 large. In order to do that it  
  ! parses the cells at level levelup=ilevel-3 and counts how many octs at 
  ! level ilevel are in the cell. Then it parses the cells at level ilevel-2 
  ! and builds patches 4x4x4. Similar at ilevel-1.
  !--------------------------------------------------------------------------
  integer::i,j,k,bit_length,ncode,nxny,nx_loc
  integer::ifirst,iskip,ngrid_ok,ilevelup,info
  integer,dimension(1:ncache)::ix,iy,iz

  integer,dimension(1:ncache)::hkeys,inext

  real(kind=8)::bscale
  real(dp)::scale
  
  integer :: nb,ntosort,hprev
  integer, dimension(1:ncache) :: itosort

  ! Local constants
  nxny=nx*ny
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)

  ! Initialisation
  ifirst=1
  do j=1,size(ngroup)
   ngroup(j)=0
  end do
  ntosort = ncache
  itosort = active(ilevel)%igrid(1:ncache)

  
  !$acc data present(isort_fin) &
  !$acc create(hkeys,ix,iy,iz) copyin(itosort)
  
  do ilevelup=levelup,ilevel-1

     bscale=2**ilevelup
     ncode=nx_loc*bscale
     do bit_length=1,32
        ncode=ncode/2
        if(ncode<=1) exit
     end do
     if(bit_length==32) then
        write(*,*)'Error in cmp_minmaxorder'
#ifndef WITHOUTMPI
        call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
        stop
#endif
     end if
          
     ! Coordinates of the octs
     !$acc parallel loop gang vector present(xg)
     do nb=1,ntosort

            ix(nb)=int((xg(itosort(nb),1)-dble(icoarse_min))*bscale)
#if NDIM>1
            iy(nb)=int((xg(itosort(nb),2)-dble(jcoarse_min))*bscale)
#endif
#if NDIM>2
            iz(nb)=int((xg(itosort(nb),3)-dble(kcoarse_min))*bscale)
#endif
     end do
     !$acc end parallel loop
     
     ! Computing the hilbert key of each oct. Octs in the same cell will
     ! have the same key.
     call hilbert3d_acc(ix(1:ntosort),iy(1:ntosort),iz(1:ntosort),&
                      & hkeys(1:ntosort),bit_length,ntosort)
     
     ! The sorting is much faster on the host.
     !$acc update host(hkeys(1:ntosort))
     call quick_sort_int(hkeys(1:ntosort),itosort(1:ntosort),ntosort)  
     ! after the sorting the octs in the same cell (so same key) will be 
     ! close one to the other in the list
     
     ! Counting how many octs have the same key
     j=1
     inext(j)=1
     hprev=hkeys(j)
     do i=2,ntosort
        if(hkeys(i)==hprev)then
           inext(j)=inext(j)+1
        else
           hprev=hkeys(i)
           j=j+1
           inext(j)=1
        endif
     end do
     ! at this point j contains the number of cells at ilevelup which contains
     ! at least one oct, and inext contains the number of octs in this cell.
          
     ! Here we check if we have enough octs in order to fulfill a cell
     ngrid_ok=(2**(ilevel-1-ilevelup))**ndim ! number of expected octs
     iskip=1
     ntosort=0
     do i=1,j
        if(inext(i)==ngrid_ok)then  ! the cell is full of octs, let's store their index
           do k=1,ngrid_ok
              isort_fin(ifirst+k-1)=itosort(iskip+k-1)
           end do
           ifirst=ifirst+ngrid_ok
           ngroup(ilevelup) = ngroup(ilevelup)+ngrid_ok
        else if(inext(i)>0) then ! the number of octs is not enough to fulfill a cell
           do k=1,inext(i)       ! at this level. So store their index for the next iteration
              itosort(ntosort+k) = itosort(iskip+k-1)
           end do
           ntosort = ntosort+inext(i)
        endif
        iskip=iskip+inext(i)
     end do
          
     if(ntosort==0) exit; ! there's no more octs to sort
     
     ! This is done only if still some octs to sort.
     ! Update the gpu for the parallel part of the next iteration.
     !$acc update device(itosort(1:ntosort))

  end do
  
  !$acc end data
  
  ! isort_fin has been updated on the host, so now we update the device
  !$acc update device(isort_fin)
  
end subroutine sort_group_grid
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine update_divu_and_enew(ilevel)
   use amr_commons
   use hydro_commons
   use acc_commons
#if defined(_OPENACC) 
  use device_memory
#endif
   implicit none
   integer::ilevel
   integer::ngrid,ncell,igrid,i,j,iskip,cell,icpu,ivar
   integer,dimension(:),allocatable::grid
   
   ngrid = active(ilevel)%ngrid + sum(reception(:,ilevel)%ngrid)
   ncell = twotondim*ngrid
   
   if(ngrid==0) return
   
   allocate(grid(1:ngrid))

   !$acc data create(grid) present(divu_acc,enew_acc,divu,enew)

   !$acc kernels
   grid(1:active(ilevel)%ngrid) = active(ilevel)%igrid(1:active(ilevel)%ngrid)
   !$acc end kernels 

   !$acc parallel loop gang
   do icpu=1,ncpu
      !$acc loop vector
      do igrid=1,reception(icpu,ilevel)%ngrid
         grid(active(ilevel)%ngrid+sum(reception(1:icpu,ilevel)%ngrid)-reception(icpu,ilevel)%ngrid+igrid) = reception(icpu,ilevel)%igrid(igrid)
      end do
   end do
   
   if(pressure_fix) then
      !$acc parallel loop collapse(2) gang vector
      do i=1,twotondim
         do igrid=1,ngrid
            iskip = ncoarse + (i-1)*ngridmax          
            cell = iskip + grid(igrid)
            divu(cell) = divu(cell) + divu_acc(cell)
            enew(cell) = enew(cell) + enew_acc(cell)
         end do
      end do
   end if

   !$acc end data
  
   deallocate(grid)

end subroutine update_divu_and_enew
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine copy_hydro_to_host(ilevel)
   use amr_commons
   use hydro_commons
   use acc_commons
#if defined(_OPENACC) 
!|| defined(USE_ACC_VERSION)
  use device_memory
#endif
   implicit none
   integer::ilevel
   !----------------------------------------------------------------------------------
   ! This routine copies the grid of uold at level ilevel in the host memory. 
   ! Moreover divu_acc and enew_acc are set to zero on the device so that they
   ! are ready to use in the next iteration.
   !----------------------------------------------------------------------------------
   integer::ngrid,ncell,igrid,ivar,i,j,k,iskip,cell,icpu
   integer,dimension(:),allocatable::grid
   
   ngrid = active(ilevel)%ngrid + sum(reception(:,ilevel)%ngrid)
   ncell = twotondim*ngrid
   
   if(ngrid==0) return
   
   allocate(grid(1:ngrid))
   
   !$acc data create(grid) present(active,reception)
   !$acc data if(pressure_fix) present(divu_acc,enew_acc)

   !$acc kernels
   grid(1:active(ilevel)%ngrid) = active(ilevel)%igrid(1:active(ilevel)%ngrid)
   !$acc end kernels 

   !$acc parallel loop gang independent
   do icpu=1,ncpu
      !$acc loop vector
      do igrid=1,reception(icpu,ilevel)%ngrid
         grid(active(ilevel)%ngrid+sum(reception(1:icpu,ilevel)%ngrid)-reception(icpu,ilevel)%ngrid+igrid) = reception(icpu,ilevel)%igrid(igrid)
      end do
   end do

   !$acc update host(grid)
   
   if(pressure_fix) then
      !$acc parallel loop collapse(2) gang vector
      do i=1,twotondim
         do igrid=1,ngrid
            iskip = ncoarse + (i-1)*ngridmax          
            cell = iskip + grid(igrid)
            divu_acc(cell) = 0.0d0
            enew_acc(cell) = 0.0d0
         end do
      end do
      !$acc end parallel loop
   end if

   do ivar=1,nvar
      call update_globalvar_dp_to_host(uold(1,ivar),ilevel)
   end do
   
   !$acc end data
   !$acc end data
  
   deallocate(grid)
   
end subroutine copy_hydro_to_host
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine godunov_fine(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine is a wrapper to the second order Godunov solver.
  ! Small grids (2x2x2) are gathered from level ilevel and sent to the
  ! hydro solver. On entry, hydro variables are gathered from array uold.
  ! On exit, unew has been updated. 
  !--------------------------------------------------------------------------
  integer::i,ivar,igrid,ncache,ngrid
  integer,dimension(1:nvector),save::ind_grid

  if(numbtot(1,ilevel)==0)return
  if(static)return
  if(verbose)write(*,111)ilevel

  ! Loop over active grids by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     call godfine1(ind_grid,ngrid,ilevel)
  end do

111 format('   Entering godunov_fine for level ',i2)

end subroutine godunov_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine set_unew(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine sets array unew to its initial value uold before calling
  ! the hydro scheme. unew is set to zero in virtual boundaries.
  !--------------------------------------------------------------------------
  integer::i,ivar,irad,ind,icpu,iskip
  real(dp)::d,u,v,w,e
  integer::ngrid,ncell,igrid,ivar,i,j,k,iskip,cell,icpu
  real(dp),dimension(:,:),allocatable::tmp_uold
  integer,dimension(:),allocatable::grid
#ifdef MPROFILING
  include 'mpif.h'
  integer::info
  real(kind=8)::tt1,tt2
  tt1=MPI_WTIME(info)
#endif

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  !$acc data present(unew,uold,divu,enew,active,reception)

  ! Set unew to uold for myid cells
  do ind=1,twotondim
  !$acc parallel async(ind)
     iskip=ncoarse+(ind-1)*ngridmax
     !$acc loop collapse(2)
     do ivar=1,nvar
        do i=1,active(ilevel)%ngrid
           unew(active(ilevel)%igrid(i)+iskip,ivar) = uold(active(ilevel)%igrid(i)+iskip,ivar)
        end do
     end do
     if(pressure_fix)then
        !$acc loop
        do i=1,active(ilevel)%ngrid
           divu(active(ilevel)%igrid(i)+iskip) = 0.0
        end do
        !$acc loop
        do i=1,active(ilevel)%ngrid
           d=max(uold(active(ilevel)%igrid(i)+iskip,1),smallr)
           u=0.0; v=0.0; w=0.0
           if(ndim>0)u=uold(active(ilevel)%igrid(i)+iskip,2)/d
           if(ndim>1)v=uold(active(ilevel)%igrid(i)+iskip,3)/d
           if(ndim>2)w=uold(active(ilevel)%igrid(i)+iskip,4)/d
           e=uold(active(ilevel)%igrid(i)+iskip,ndim+2)-0.5*d*(u**2+v**2+w**2)
#if NENER>0
           do irad=1,nener
              e=e-uold(active(ilevel)%igrid(i)+iskip,ndim+2+irad)
           end do
#endif  
           enew(active(ilevel)%igrid(i)+iskip)=e
        end do
     end if
  !$acc end parallel
  end do

  ! Set unew to 0 for virtual boundary cells
  do icpu=1,ncpu
  !$acc parallel async(icpu)
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     !$acc loop collapse(2)
     do ivar=1,nvar
        do i=1,reception(icpu,ilevel)%ngrid
           unew(reception(icpu,ilevel)%igrid(i)+iskip,ivar)=0.0
        end do
     end do
     if(pressure_fix)then
        !$acc loop
        do i=1,reception(icpu,ilevel)%ngrid
           divu(reception(icpu,ilevel)%igrid(i)+iskip) = 0.0
           enew(reception(icpu,ilevel)%igrid(i)+iskip) = 0.0
        end do
     end if
  end do
  !$acc end parallel
  end do
  !$acc wait

  !$acc end data

#ifdef MPROFILING
  tt2=MPI_WTIME(info)
  acc_t_set_unew = acc_t_set_unew + (tt2-tt1)
#endif  

111 format('   Entering set_unew for level ',i2)

end subroutine set_unew
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine set_uold(ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ilevel
  !---------------------------------------------------------
  ! This routine sets array uold to its new value unew 
  ! after the hydro step.
  !---------------------------------------------------------
  integer::i,ivar,irad,ind,iskip,nx_loc,ind_cell,icpu
  real(dp)::scale,d,u,v,w
  real(dp)::e_kin,e_cons,e_prim,e_trunc,div,dx,fact,d_old
#ifdef MPROFILING
  include 'mpif.h'
  integer::info
  real(kind=8)::tt1,tt2
  tt1=MPI_WTIME(info)
#endif

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5d0**ilevel*scale

  ! Add gravity source term at time t with half time step
  if(poisson)then
     call add_gravity_source_terms(ilevel)
  end if

  ! Add non conservative pdV terms to unew 
  ! for thermal and/or non-thermal energies
  if(pressure_fix.OR.nener>0)then
     call add_pdv_source_terms(ilevel)
  endif

  !$acc data present(unew,uold,divu,enew,active) copyin(dtnew)

  ! Set uold to unew for myid cells
  !$acc parallel loop collapse(2) gang independent
  do ind=1,twotondim
     do ivar=1,nvar 
        !$acc loop vector
        do i=1,active(ilevel)%ngrid
           iskip=ncoarse+(ind-1)*ngridmax
           uold(active(ilevel)%igrid(i)+iskip,ivar) = unew(active(ilevel)%igrid(i)+iskip,ivar)
        end do
     end do
  end do

  if(pressure_fix)then
  !$acc parallel loop independent
  do ind=1,twotondim
     !$acc loop vector
     do i=1,active(ilevel)%ngrid
        iskip=ncoarse+(ind-1)*ngridmax
           ind_cell=active(ilevel)%igrid(i)+iskip
           d=max(uold(ind_cell,1),smallr)
           u=0.0; v=0.0; w=0.0
           if(ndim>0)u=uold(ind_cell,2)/d
           if(ndim>1)v=uold(ind_cell,3)/d
           if(ndim>2)w=uold(ind_cell,4)/d
           e_kin=0.5*d*(u**2+v**2+w**2)
#if NENER>0
           do irad=1,nener
              e_kin=e_kin+uold(ind_cell,ndim+2+irad)
           end do
#endif
           e_cons=uold(ind_cell,ndim+2)-e_kin
           e_prim=enew(ind_cell)
           ! Note: here divu=-div.u*dt
           div=abs(divu(ind_cell))*dx/dtnew(ilevel)
           e_trunc=beta_fix*d*max(div,3.0*hexp*dx)**2
           if(e_cons<e_trunc)then
              uold(ind_cell,ndim+2)=e_prim+e_kin
           end if
      end do
  end do
  end if

  !$acc end data
  
#ifdef MPROFILING
  tt2=MPI_WTIME(info)
  acc_t_set_uold = acc_t_set_uold + (tt2-tt1)
#endif  

111 format('   Entering set_uold for level ',i2)

end subroutine set_uold
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine godfine1(ind_grid,ncache,ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ilevel,ncache
  integer,dimension(1:nvector)::ind_grid
  !-------------------------------------------------------------------
  ! This routine gathers first hydro variables from neighboring grids
  ! to set initial conditions in a 6x6x6 grid. It interpolate from
  ! coarser level missing grid variables. It then calls the
  ! Godunov solver that computes fluxes. These fluxes are zeroed at 
  ! coarse-fine boundaries, since contribution from finer levels has
  ! already been taken into account. Conservative variables are updated 
  ! and stored in array unew(:), both at the current level and at the 
  ! coarser level if necessary.
  !-------------------------------------------------------------------
  integer ,dimension(1:nvector,1:threetondim     ),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim       ),save::nbors_father_grids
  integer ,dimension(1:nvector,0:twondim         ),save::ibuffer_father
  real(dp),dimension(1:nvector,0:twondim  ,1:nvar),save::u1
  real(dp),dimension(1:nvector,1:twotondim,1:nvar),save::u2
  real(dp),dimension(1:nvector,0:twondim  ,1:ndim),save::g1=0.0d0
  real(dp),dimension(1:nvector,1:twotondim,1:ndim),save::g2=0.0d0

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar),save::uloc
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim),save::gloc=0.0d0
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:nvar,1:ndim),save::flux
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:2,1:ndim),save::tmp
  logical ,dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::ok

  integer,dimension(1:nvector),save::igrid_nbor,ind_cell,ind_buffer,ind_exist,ind_nexist

  integer::i,j,ivar,idim,ind_son,ind_father,iskip,nbuffer,ibuffer
  integer::i0,j0,k0,i1,j1,k1,i2,j2,k2,i3,j3,k3,nx_loc,nb_noneigh,nexist
  integer::i1min,i1max,j1min,j1max,k1min,k1max
  integer::i2min,i2max,j2min,j2max,k2min,k2max
  integer::i3min,i3max,j3min,j3max,k3min,k3max
  real(dp)::dx,scale,oneontwotondim

  oneontwotondim = 1.d0/dble(twotondim)

  ! Mesh spacing in that level
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5D0**ilevel*scale
  
  ! Integer constants
  i1min=0; i1max=0; i2min=0; i2max=0; i3min=1; i3max=1
  j1min=0; j1max=0; j2min=0; j2max=0; j3min=1; j3max=1
  k1min=0; k1max=0; k2min=0; k2max=0; k3min=1; k3max=1
  if(ndim>0)then
     i1max=2; i2max=1; i3max=2
  end if
  if(ndim>1)then
     j1max=2; j2max=1; j3max=2
  end if
  if(ndim>2)then
     k1max=2; k2max=1; k3max=2
  end if

  !------------------------------------------
  ! Gather 3^ndim neighboring father cells
  !------------------------------------------
  do i=1,ncache
     ind_cell(i)=father(ind_grid(i))
  end do
  ! This was a call introduced with OpenACC
  ! call get3cubefather_godfine(ind_cell,nbors_father_cells,nbors_father_grids,ncache,ilevel)
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ncache,ilevel)
  
  !---------------------------
  ! Gather 6x6x6 cells stencil
  !---------------------------
  ! Loop over 3x3x3 neighboring father cells
  do k1=k1min,k1max
  do j1=j1min,j1max
  do i1=i1min,i1max
     
     ! Check if neighboring grid exists
     nbuffer=0
     nexist=0
     ind_father=1+i1+3*j1+9*k1
     do i=1,ncache
        igrid_nbor(i)=son(nbors_father_cells(i,ind_father))
        if(igrid_nbor(i)>0) then
           nexist=nexist+1
           ind_exist(nexist)=i
        else
          nbuffer=nbuffer+1
          ind_nexist(nbuffer)=i
          ind_buffer(nbuffer)=nbors_father_cells(i,ind_father)
        end if
     end do
     
     ! If not, interpolate hydro variables from parent cells
     if(nbuffer>0)then
        call getnborfather(ind_buffer,ibuffer_father,nbuffer,ilevel)
        do j=0,twondim
           do ivar=1,nvar
              do i=1,nbuffer
                 u1(i,j,ivar)=uold(ibuffer_father(i,j),ivar)
              end do
           end do
        end do
        call interpol_hydro(u1,u2,nbuffer)
     endif

     ! Loop over 2x2x2 cells
     do k2=k2min,k2max
     do j2=j2min,j2max
     do i2=i2min,i2max

        ind_son=1+i2+2*j2+4*k2
        iskip=ncoarse+(ind_son-1)*ngridmax
        do i=1,nexist
           ind_cell(i)=iskip+igrid_nbor(ind_exist(i))
        end do
        
        i3=1; j3=1; k3=1
        if(ndim>0)i3=1+2*(i1-1)+i2
        if(ndim>1)j3=1+2*(j1-1)+j2
        if(ndim>2)k3=1+2*(k1-1)+k2
        
        ! Gather hydro variables
        do ivar=1,nvar
           do i=1,nexist
              uloc(ind_exist(i),i3,j3,k3,ivar)=uold(ind_cell(i),ivar)
           end do
           do i=1,nbuffer
              uloc(ind_nexist(i),i3,j3,k3,ivar)=u2(i,ind_son,ivar)
           end do
        end do
        
        ! Gather gravitational acceleration
        if(poisson)then
           do idim=1,ndim
              do i=1,nexist
                 gloc(ind_exist(i),i3,j3,k3,idim)=f(ind_cell(i),idim)
              end do
              ! Use straight injection for buffer cells
              do i=1,nbuffer
                 gloc(ind_nexist(i),i3,j3,k3,idim)=f(ibuffer_father(i,0),idim)
              end do
           end do
        end if
        
        ! Gather refinement flag
        do i=1,nexist
           ok(ind_exist(i),i3,j3,k3)=son(ind_cell(i))>0
        end do
        do i=1,nbuffer
           ok(ind_nexist(i),i3,j3,k3)=.false.
        end do
        
     end do
     end do
     end do
     ! End loop over cells

  end do
  end do
  end do
  ! End loop over neighboring grids

  !-----------------------------------------------
  ! Compute flux using second-order Godunov method
  !-----------------------------------------------
  call unsplit(uloc,gloc,flux,tmp,dx,dx,dx,dtnew(ilevel),ncache)

  !------------------------------------------------
  ! Reset flux along direction at refined interface    
  !------------------------------------------------
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     do k3=k3min,k3max+k0
     do j3=j3min,j3max+j0
     do i3=i3min,i3max+i0
        do ivar=1,nvar
           do i=1,ncache
              if(ok(i,i3-i0,j3-j0,k3-k0) .or. ok(i,i3,j3,k3))then
                 flux(i,i3,j3,k3,ivar,idim)=0.0d0
              end if
           end do
        end do
        if(pressure_fix)then
        do ivar=1,2
           do i=1,ncache
              if(ok(i,i3-i0,j3-j0,k3-k0) .or. ok(i,i3,j3,k3))then
                 tmp (i,i3,j3,k3,ivar,idim)=0.0d0
              end if
           end do
        end do
        end if
     end do
     end do
     end do
  end do
  !--------------------------------------
  ! Conservative update at level ilevel
  !--------------------------------------
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     do k2=k2min,k2max
     do j2=j2min,j2max
     do i2=i2min,i2max
        ind_son=1+i2+2*j2+4*k2
        iskip=ncoarse+(ind_son-1)*ngridmax
        do i=1,ncache
           ind_cell(i)=iskip+ind_grid(i)
        end do
        i3=1+i2
        j3=1+j2
        k3=1+k2
        ! Update conservative variables new state vector
        do ivar=1,nvar
           do i=1,ncache
              unew(ind_cell(i),ivar)=unew(ind_cell(i),ivar)+ &
                   & (flux(i,i3   ,j3   ,k3   ,ivar,idim) &
                   & -flux(i,i3+i0,j3+j0,k3+k0,ivar,idim))
           end do
        end do
        if(pressure_fix)then
        ! Update velocity divergence
        do i=1,ncache
           divu(ind_cell(i))=divu(ind_cell(i))+ &
                & (tmp(i,i3   ,j3   ,k3   ,1,idim) &
                & -tmp(i,i3+i0,j3+j0,k3+k0,1,idim))
        end do
        ! Update internal energy
        do i=1,ncache
           enew(ind_cell(i))=enew(ind_cell(i))+ &
                & (tmp(i,i3   ,j3   ,k3   ,2,idim) &
                & -tmp(i,i3+i0,j3+j0,k3+k0,2,idim))
        end do
        end if
     end do
     end do
     end do
  end do

  !--------------------------------------
  ! Conservative update at level ilevel-1
  !--------------------------------------
  ! Loop over dimensions
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     
     !----------------------
     ! Left flux at boundary
     !----------------------     
     ! Check if grids sits near left boundary
     ! and gather neighbor father cells index
     nb_noneigh=0
     do i=1,ncache
        if (son(nbor(ind_grid(i),2*idim-1))==0) then
           nb_noneigh = nb_noneigh + 1
           ind_buffer(nb_noneigh) = nbor(ind_grid(i),2*idim-1)
           ind_cell(nb_noneigh) = i
        end if
     end do
     ! Conservative update of new state variables
     do ivar=1,nvar
        ! Loop over boundary cells
        do k3=k3min,k3max-k0
        do j3=j3min,j3max-j0
        do i3=i3min,i3max-i0
           do i=1,nb_noneigh
              unew(ind_buffer(i),ivar)=unew(ind_buffer(i),ivar) &
                   & -flux(ind_cell(i),i3,j3,k3,ivar,idim)*oneontwotondim
           end do
        end do
        end do
        end do
     end do
     if(pressure_fix)then
     ! Update velocity divergence
     do k3=k3min,k3max-k0
     do j3=j3min,j3max-j0
     do i3=i3min,i3max-i0
        do i=1,nb_noneigh
           divu(ind_buffer(i))=divu(ind_buffer(i)) &
                & -tmp(ind_cell(i),i3,j3,k3,1,idim)*oneontwotondim
        end do
     end do
     end do
     end do
     ! Update internal energy
     do k3=k3min,k3max-k0
     do j3=j3min,j3max-j0
     do i3=i3min,i3max-i0
        do i=1,nb_noneigh
           enew(ind_buffer(i))=enew(ind_buffer(i)) &
                & -tmp(ind_cell(i),i3,j3,k3,2,idim)*oneontwotondim
        end do
     end do
     end do
     end do
     end if
     
     !-----------------------
     ! Right flux at boundary
     !-----------------------     
     ! Check if grids sits near right boundary
     ! and gather neighbor father cells index
     nb_noneigh=0
     do i=1,ncache
        if (son(nbor(ind_grid(i),2*idim))==0) then
           nb_noneigh = nb_noneigh + 1
           ind_buffer(nb_noneigh) = nbor(ind_grid(i),2*idim)
           ind_cell(nb_noneigh) = i
        end if
     end do
     ! Conservative update of new state variables
     do ivar=1,nvar
        ! Loop over boundary cells
        do k3=k3min+k0,k3max
        do j3=j3min+j0,j3max
        do i3=i3min+i0,i3max
           do i=1,nb_noneigh
              unew(ind_buffer(i),ivar)=unew(ind_buffer(i),ivar) &
                   & +flux(ind_cell(i),i3+i0,j3+j0,k3+k0,ivar,idim)*oneontwotondim
           end do
        end do
        end do
        end do
     end do
     if(pressure_fix)then
     ! Update velocity divergence
     do k3=k3min+k0,k3max
     do j3=j3min+j0,j3max
     do i3=i3min+i0,i3max
        do i=1,nb_noneigh
           divu(ind_buffer(i))=divu(ind_buffer(i)) &
                & +tmp(ind_cell(i),i3+i0,j3+j0,k3+k0,1,idim)*oneontwotondim
        end do
     end do
     end do
     end do
     ! Update internal energy
     do k3=k3min+k0,k3max
     do j3=j3min+j0,j3max
     do i3=i3min+i0,i3max
        do i=1,nb_noneigh
           enew(ind_buffer(i))=enew(ind_buffer(i)) &
                & +tmp(ind_cell(i),i3+i0,j3+j0,k3+k0,2,idim)*oneontwotondim
        end do
     end do
     end do
     end do
     end if

  end do
  ! End loop over dimensions

end subroutine godfine1

!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine add_gravity_source_terms(ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine adds to unew the gravity source terms
  ! with only half a time step. Only the momentum and the
  ! total energy are modified in array unew.
  !--------------------------------------------------------------------------
  integer::i,ivar,ind,iskip,nx_loc,ind_cell
  real(dp)::d,u,v,w,e_kin,e_prim,d_old,fact
#ifdef MPROFILING
  include 'mpif.h'
  integer::info
  real(kind=8)::tt1,tt2
  tt1=MPI_WTIME(info)
#endif

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  !$acc data present(unew,uold,active) copyin(dtnew)

  ! Add gravity source term at time t with half time step
  do ind=1,twotondim
  !$acc parallel async(ind)
     iskip=ncoarse+(ind-1)*ngridmax
     !$acc loop
     do i=1,active(ilevel)%ngrid
        ind_cell=active(ilevel)%igrid(i)+iskip
        d=max(unew(ind_cell,1),smallr)
        u=0.0; v=0.0; w=0.0
        if(ndim>0)u=unew(ind_cell,2)/d
        if(ndim>1)v=unew(ind_cell,3)/d
        if(ndim>2)w=unew(ind_cell,4)/d
        e_kin=0.5*d*(u**2+v**2+w**2)
        e_prim=unew(ind_cell,ndim+2)-e_kin
        d_old=max(uold(ind_cell,1),smallr)
        fact=d_old/d*0.5*dtnew(ilevel)
        if(ndim>0)then
           u=u+f(ind_cell,1)*fact
           unew(ind_cell,2)=d*u
        endif
        if(ndim>1)then
           v=v+f(ind_cell,2)*fact
           unew(ind_cell,3)=d*v
        end if
        if(ndim>2)then
           w=w+f(ind_cell,3)*fact
           unew(ind_cell,4)=d*w
        endif
        e_kin=0.5*d*(u**2+v**2+w**2)
        unew(ind_cell,ndim+2)=e_prim+e_kin
     end do
  !$acc end parallel
  end do
  !$acc wait

  !$acc end data
  
#ifdef MPROFILING
  tt2=MPI_WTIME(info)
  acc_t_add_gravity_source_terms = acc_t_add_gravity_source_terms + (tt2-tt1)
#endif  

111 format('   Entering add_gravity_source_terms for level ',i2)

end subroutine add_gravity_source_terms
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine add_pdv_source_terms(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !---------------------------------------------------------
  ! This routine adds the pdV source term to the internal
  ! energy equation and to the non-thermal energy equations.
  !---------------------------------------------------------
  integer::i,ivar,irad,ind,iskip,nx_loc,ind_cell1
  integer::ncache,igrid,ngrid,idim,id1,ig1,ih1,id2,ig2,ih2
  integer,dimension(1:3,1:2,1:8)::iii,jjj
  real(dp)::scale,dx,dx_loc,d,u,v,w,eold

  integer ,dimension(1:nvector),save::ind_grid,ind_cell
  integer ,dimension(1:nvector,0:twondim),save::igridn
  integer ,dimension(1:nvector,1:ndim),save::ind_left,ind_right
  real(dp),dimension(1:nvector,1:ndim,1:ndim),save::velg,veld
  real(dp),dimension(1:nvector,1:ndim),save::dx_g,dx_d
  real(dp),dimension(1:nvector),save::divu_loc
#ifdef MPROFILING
  include 'mpif.h'
  integer::info
  real(kind=8)::tt1,tt2
  tt1=MPI_WTIME(info)
#endif

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5d0**ilevel
  dx_loc=dx*scale

  iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

#if NENER>0
    do irad=1,nener
       call update_globalvar_dp_to_host(unew(1,ndim+2+irad),ilevel)
    end do
#endif
  call update_globalvar_dp_to_host(enew,ilevel)

  ! Loop over myid grids by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
   
     ! Gather nvector grids
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     
     ! Gather neighboring grids
     do i=1,ngrid
        igridn(i,0)=ind_grid(i)
     end do
     do idim=1,ndim
        do i=1,ngrid
           ind_left (i,idim)=nbor(ind_grid(i),2*idim-1)
           ind_right(i,idim)=nbor(ind_grid(i),2*idim  )
           igridn(i,2*idim-1)=son(ind_left (i,idim))
           igridn(i,2*idim  )=son(ind_right(i,idim))
        end do
     end do
     
     ! Loop over cells
     do ind=1,twotondim
        
        ! Compute central cell index
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        
        ! Gather all neighboring velocities
        do idim=1,ndim
           id1=jjj(idim,1,ind); ig1=iii(idim,1,ind)
           ih1=ncoarse+(id1-1)*ngridmax
           do i=1,ngrid
              if(igridn(i,ig1)>0)then
                 velg(i,idim,1:ndim) = uold(igridn(i,ig1)+ih1,2:ndim+1)/max(uold(igridn(i,ig1)+ih1,1),smallr)
                 dx_g(i,idim) = dx_loc
              else
                 velg(i,idim,1:ndim) = uold(ind_left(i,idim),2:ndim+1)/max(uold(ind_left(i,idim),1),smallr)
                 dx_g(i,idim) = dx_loc*1.5_dp
              end if
           enddo
           id2=jjj(idim,2,ind); ig2=iii(idim,2,ind)
           ih2=ncoarse+(id2-1)*ngridmax
           do i=1,ngrid
              if(igridn(i,ig2)>0)then
                 veld(i,idim,1:ndim)= uold(igridn(i,ig2)+ih2,2:ndim+1)/max(uold(igridn(i,ig2)+ih2,1),smallr)
                 dx_d(i,idim)=dx_loc
              else 
                 veld(i,idim,1:ndim)= uold(ind_right(i,idim),2:ndim+1)/max(uold(ind_right(i,idim),1),smallr)
                 dx_d(i,idim)=dx_loc*1.5_dp
              end if
           enddo
        end do
        ! End loop over dimensions
  
        ! Compute divu = Trace G
        divu_loc(1:ngrid)=0.0d0
        do i=1,ngrid
           do idim=1,ndim
              divu_loc(i) = divu_loc(i) + (veld(i,idim,idim)-velg(i,idim,idim)) &
                   &                    / (dx_g(i,idim)     +dx_d(i,idim))
           enddo
        end do

        ! Update thermal internal energy 
        if(pressure_fix)then
           do i=1,ngrid
              ! Compute old thermal energy
              d=max(uold(ind_cell(i),1),smallr)
              u=0.0; v=0.0; w=0.0
              if(ndim>0)u=uold(ind_cell(i),2)/d
              if(ndim>1)v=uold(ind_cell(i),3)/d
              if(ndim>2)w=uold(ind_cell(i),4)/d
              eold=uold(ind_cell(i),ndim+2)-0.5*d*(u**2+v**2+w**2)
#if NENER>0
              do irad=1,nener
                 eold=eold-uold(ind_cell(i),ndim+2+irad)
              end do
#endif
              ! Add -pdV term
              enew(ind_cell(i))=enew(ind_cell(i)) &
                   & -(gamma-1.0d0)*eold*divu_loc(i)*dtnew(ilevel)
           end do
        end if

#if NENER>0
        do irad=1,nener
           do i=1,ngrid
              ! Add -pdV term
              unew(ind_cell(i),ndim+2+irad)=unew(ind_cell(i),ndim+2+irad) &
                & -(gamma_rad(irad)-1.0d0)*uold(ind_cell(i),ndim+2+irad)*divu_loc(i)*dtnew(ilevel)
           end do
        end do
#endif

     enddo
     ! End loop over cells
  end do
  ! End loop over grids

#if NENER>0
    do irad=1,nener
       call update_globalvar_dp_to_device(unew(1,ndim+2+irad),ilevel)
    end do
#endif
  call update_globalvar_dp_to_device(enew,ilevel)

#ifdef MPROFILING
  tt2=MPI_WTIME(info)
  acc_t_add_pdv_source_terms = acc_t_add_pdv_source_terms + (tt2-tt1)
#endif 

111 format('   Entering add_pdv_source_terms for level ',i2)

end subroutine add_pdv_source_terms

