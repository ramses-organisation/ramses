!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine fill_hydro_grids_interpol_1(nxp,ilevel,dx,dx_box,nxp4,nc1,nc2,icube,ngrid_ok,isort,ok,box_xmin)
  use amr_commons
  use hydro_commons
  use poisson_commons
  use acc_commons
  implicit none
  integer::ilevel,nxp,nxp4,nc1,nc2,ngrid_ok,icube
  integer,dimension(active(ilevel)%ngrid) :: isort
  real(dp)::dx,dx_box
  real(dp),dimension(1:ndim,nc1:nc2):: box_xmin
  logical, dimension(-1:nxp+2,-1:nxp+2,-1:nxp+2,nc1:nc2)::ok
  !--------------------------------------------------------------------------
  ! This routine fills up a Cartesian grid from AMR data.
  ! It iterates over the patches and the cells. For each cell in the patch it 
  ! computes its coordinates and checks if this cell exists in the AMR grid.
  ! If yes then it simply takes the values from uold. Otherwise it searches
  ! for the neighbouring cells and does and interpolation.
  !--------------------------------------------------------------------------
  ! fill_hydro_grids_interpol_1 is faster than the two other filling routines
  ! when ncube = O(50). In a test with ncube=56 it tourned out that it is 1.2
  ! times faster. 
  !--------------------------------------------------------------------------
  integer::i,j,k,l,idim,ivar,nc,ic,jc,kc
  integer::cell_index,cell_levl,cell_grid
  integer::uindex, gindex, bindex
  integer, dimension(0:twondim)::ibuffer_father
  real(dp),dimension(0:twondim,1:nvar)::u1
  real(dp),dimension(1:twotondim,1:nvar)::u2
  real(dp),dimension(1:ndim):: xx_dp
  
  ! Private variables for getnborfather_acc
  integer,dimension(1:3)::ibound,iskip1,iskip2
  integer,dimension(1:3)::ix
  integer,dimension(0:twondim)::igridn
  integer,dimension(1:twondim)::icelln_ok 

  
  !$acc data present(igrid_acc,ucount0,gcount0) &
  !$acc present(dx,dx_box,ilevel,nxp4) &
  !$acc present(uold,uloc_tot,gloc_tot,f) &
  !$acc present(isort,xg,son,nbor,ok,box_xmin)
  
  !$acc parallel loop
  do nc=nc1,nc2
     ! Index of one of the octs in the patch.
     igrid_acc(nc-nc1+1)=isort(icube+ngrid_ok*(nc-nc1))
     ! Computation of the starintg and last index of the patch
     ucount0(nc-nc1+1) = 1 + (nc-nc1)*nvar*nxp4**3
     gcount0(nc-nc1+1) = 1 + (nc-nc1)*ndim*nxp4**3
  end do
  !$acc end parallel loop
  
  ! Compute box coordinates in normalized unites
  !$acc parallel loop collapse(2)
  do nc=nc1,nc2
   do idim=1,ndim
      box_xmin(idim,nc)=int(xg(igrid_acc(nc-nc1+1),idim)/dx_box)*dx_box ! corner of the box
   end do
  end do
  !$acc end parallel loop
    
  !$acc parallel loop collapse(4) gang vector present(ggg_acc,hhh_acc) &
  !$acc private(xx_dp,u1,u2,ibuffer_father,ibound,iskip1,iskip2,ix,igridn,icelln_ok)
  do nc=nc1,nc2
     
   do k=-1,nxp+2
#if NDIM>1
      do j=-1,nxp+2
#endif
#if NDIM>2
         do i=-1,nxp+2
#endif         
            ! Coordinate of the cell
            xx_dp(1) = box_xmin(1,nc) + (dble(i)-0.5)*dx
            xx_dp(2) = box_xmin(2,nc) + (dble(j)-0.5)*dx
            xx_dp(3) = box_xmin(3,nc) + (dble(k)-0.5)*dx
           
            ! Compute cell index
            call hydro_get_cell_index_scalar(cell_index,cell_grid,cell_levl,xx_dp,ilevel)
            
            ! Index of the cell in the patch
            bindex = i+2 + nxp4*(j+1) + nxp4*nxp4*(k+1)
            
            if(cell_levl==ilevel)then   ! the cell exists in the AMR grid        
               ok(i,j,k,nc) = son(cell_index)>0  ! checks if the cell is refined
               do ivar=1,nvar
                  uindex = (ucount0(nc-nc1+1)-1) + bindex + (ivar-1)*nxp4**3 ! index of ivar in uloc_tot
                  uloc_tot(uindex) = uold(cell_index,ivar) ! Store hydro variable in local Cartesian grid
               end do
            else ! the does not exists in the AMR grid ...
               ok(i,j,k,nc)=.false. ! ... and so it has no sons
               
               ! Searching for father and its neighbours
               call getnborfather_acc(cell_index,ibuffer_father,ilevel,ibound,iskip1,iskip2,ix,igridn,icelln_ok)
               
               ! Storing the values of father and uncles in u1
               do idim=0,twondim
                  do ivar=1,nvar
                    u1(idim,ivar)=uold(ibuffer_father(idim),ivar)
                   end do
               end do 
               
               ! This interpolation will compute a full oct but just one cell is needed.
               call interpol_hydro_one(u1,u2)
               
               ! index of the cell in u2
               ic=0;jc=0;kc=0;
               if(mod(i,2)==0) ic=1;
               if(mod(j,2)==0) jc=1;
               if(mod(k,2)==0) kc=1;
               ! updating uloc_tot with the interpolated value
               do ivar=1,nvar
                  uindex = (ucount0(nc-nc1+1)-1) + bindex + (ivar-1)*nxp4**3
                  uloc_tot(uindex) = u2(1+ic+2*jc+4*kc,ivar)
               end do
            end if
            
            ! f is not interpolated. If the cell does not exists in the amr grid then the value
            ! of the father is taken.
            do idim=1,ndim
               gindex = (gcount0(nc-nc1+1)-1) + bindex + (idim-1)*nxp4**3
               if(poisson) then
                 gloc_tot(gindex) = f(cell_index,idim)
               else
                 gloc_tot(gindex) = 0.0d0
               end if
            end do

#if NDIM>2
         end do
#endif
#if NDIM>1
      end do
#endif
   end do
  
  end do
  !$acc end parallel loop
    
  !$acc end data
  
end subroutine fill_hydro_grids_interpol_1
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine fill_hydro_grids_interpol_2(nxp,ilevel,dx,dx_box,nxp4,nc1,nc2,icube,ngrid_ok,isort,ok,box_xmin,index_list,grid_list,bindex,z,maxz)
  use amr_commons
  use hydro_commons
  use poisson_commons
  use acc_commons
  implicit none
  integer::ilevel,nxp,nxp4,nc1,nc2,ngrid_ok,icube,maxz
  integer,dimension(active(ilevel)%ngrid) :: isort
  real(dp)::dx,dx_box
  real(dp),dimension(1:ndim,nc1:nc2):: box_xmin
  integer, dimension(1:(nxp/2+2)**3,1:ncube)::index_list,grid_list,bindex
  integer, dimension(1:ncube)::z
  logical, dimension(-1:nxp+2,-1:nxp+2,-1:nxp+2,nc1:nc2)::ok
  !--------------------------------------------------------------------------
  ! This routine fills up a Cartesian grid from AMR data.
  ! It iterates over the patches and the father cells. For each father cell 
  ! in the patch it computes its coordinates, checks if it has sons and
  ! stores their index in a list. Then the list is reordered and the father
  ! cells which need interpolation are put at the beginning of the list. At
  ! this point the interpolation occurs. Finally the patch is built using
  ! uold and the interpolated values.
  !--------------------------------------------------------------------------
  ! fill_hydro_grids_interpol_2 is much faster than  
  ! fill_hydro_grids_interpol_1 when ncube = O(500). In a test with ncube=500
  ! it was 2.4 times faster.
  !--------------------------------------------------------------------------
  integer::i,j,k,l,idim,ivar,nc,pos,ic,jc,kc
  integer::cell_index,cell_levl,cell_grid
  integer :: uindex, gindex
  integer, dimension(0:twondim)::ibuffer_father
  real(dp),dimension(:,:,:),allocatable::u2
  real(dp),dimension(1:ndim):: xx_dp
  
  ! Private variables for getnborfather_acc
  integer,dimension(1:3)::ibound,iskip1,iskip2
  integer,dimension(1:3)::ix
  integer,dimension(0:twondim)::igridn
  integer,dimension(1:twondim)::icelln_ok 

  !$acc data present(igrid_acc,isort,xg,son,ucount0,gcount0,dx,dx_box,ilevel,nxp4,uold,uloc_tot,gloc_tot,f,ok,box_xmin) &
  !$acc present(index_list,grid_list,bindex,z,maxz)

  !$acc parallel loop
  do nc=nc1,nc2
     ! Index of one of the octs in the patch.
     igrid_acc(nc-nc1+1)=isort(icube+ngrid_ok*(nc-nc1))
     ! Computation of the starintg index of the patch
     ucount0(nc-nc1+1) = 1 + (nc-nc1)*nvar*nxp4**3
     gcount0(nc-nc1+1) = 1 + (nc-nc1)*ndim*nxp4**3
  end do
  !$acc end parallel loop
   
  ! Compute box coordinates in normalized unites
  !$acc parallel loop collapse(2)
  do nc=nc1,nc2
   do idim=1,ndim
      box_xmin(idim,nc)=int(xg(igrid_acc(nc-nc1+1),idim)/dx_box)*dx_box ! corner of the box
   end do
  end do
  !$acc end parallel loop
  
  
  !$acc parallel loop collapse(4) gang vector private(xx_dp)
  do nc=nc1,nc2
   do k=0,nxp/2+1
#if NDIM>1
      do j=0,nxp/2+1
#endif
#if NDIM>2
         do i=0,nxp/2+1
#endif
            ! Computing the coordinate of the father cell
            xx_dp(1) = box_xmin(1,nc) + (dble(i)-0.5)*2.0d0*dx
            xx_dp(2) = box_xmin(2,nc) + (dble(j)-0.5)*2.0d0*dx
            xx_dp(3) = box_xmin(3,nc) + (dble(k)-0.5)*2.0d0*dx         
              
            call hydro_get_cell_index_scalar(cell_index,cell_grid,cell_levl,xx_dp,ilevel-1)
  
            pos = i+1 + j*(nxp/2+2) + k*(nxp/2+2)**2 ! index of the father cell in the following lists
            index_list(pos,nc-nc1+1) = cell_index    ! index of the father cell
            grid_list(pos,nc-nc1+1) = cell_grid      ! index of the grid
            bindex(pos,nc-nc1+1) =  1 + 2*i + 2*nxp4*j + 2*nxp4*nxp4*k ! index in the patch of the first cell of the grid
                                                                    
            if(cell_grid==0) then ! the father cell has no son...
               do kc=0,1
               do jc=0,1
               do ic=0,1
                  ok(2*i+ic-1,2*j+jc-1,2*k+kc-1,nc)=.false. !...so the cells neither
               end do
               end do
               end do
            else  ! the father cell has a son
               do kc=0,1
               do jc=0,1
               do ic=0,1
                  cell_index = ncoarse+(ic+2*jc+4*kc)*ngridmax+cell_grid ! index of the cells
                  ok(2*i+ic-1,2*j+jc-1,2*k+kc-1,nc) = son(cell_index)>0  ! check if they are refined
               end do
               end do
               end do
            end if 
            
#if NDIM>2
         end do
#endif
#if NDIM>1
      end do
#endif
   end do
  
  end do  
  !$acc end parallel loop
  
  !$acc kernels
  maxz = 0.0d0
  !$acc end kernels
    
  ! Putting the father cells which need interpolation at the beginning of the list

  !$acc parallel loop gang vector reduction(max:maxz)
  do nc=nc1,nc2     
     ! putting the zeros in grid_list(:,nc) at the beginning of the list. index_list and bindex will be reordered
     ! z will contain the number of zeros found, so the number of cells to be interpolated
     call sort_zeros(grid_list(:,nc-nc1+1), index_list(:,nc-nc1+1), bindex(:,nc-nc1+1), (nxp/2+2)**3, z(nc-nc1+1))
     maxz = max(maxz,z(nc-nc1+1))
  end do
  !$acc end parallel loop
  
  !$acc update host(maxz)
  allocate(u2(1:maxz*(nc2-nc1+1),1:twotondim,1:nvar))
  
  !$acc data create(u2)
  
  ! Storing the variables needed for the interpolation
  !$acc parallel loop collapse(2) gang vector &
  !$acc present(ggg_acc,hhh_acc,nbor) private(ibuffer_father,ibound,iskip1,iskip2,ix,igridn,icelln_ok)
  do nc=nc1,nc2
     do l=1,maxz
        
        if(l<=z(nc-nc1+1)) then
         call getnborfather_acc(index_list(l,nc-nc1+1),ibuffer_father,ilevel,&
                               &ibound,iskip1,iskip2,ix,igridn,icelln_ok)

         ! Here an empty chunk of uloc_tot is used in order to store the variables for the interpolation.
         do j=0,twondim
            do ivar=1,nvar
               uloc_tot(l+(nc-nc1)*maxz + j*(nc2-nc1+1)*maxz + (ivar-1)*(twondim+1)*(nc2-nc1+1)*maxz)=uold(ibuffer_father(j),ivar)
            end do
         end do
        end if
        
     end do
  end do
  !$acc end parallel loop
  
  ! Parallel interpolation. The new values will be in u2
  call interpol_hydro_acc(uloc_tot(1:nvar*(twondim+1)*(nc2-nc1+1)*maxz),u2,maxz*(nc2-nc1+1),nc1,nc2,maxz,z)
  
  ! Filling uloc_tot and gloc_tot
  !$acc parallel loop collapse(5) gang vector
  do nc=nc1,nc2 ! loop over patches
     do l=1,(nxp/2+2)**3 ! loop over father cells
        do k=0,1 ! loops over cells
        do j=0,1
        do i=0,1
           pos = bindex(l,nc-nc1+1) + i + j*nxp4 + k*nxp4*nxp4 ! index in the patch of the cell
           cell_index = ncoarse + (i+2*j+4*k)*ngridmax + grid_list(l,nc-nc1+1)
           do ivar=1,nvar
               uindex = (ucount0(nc-nc1+1)-1) + pos + (ivar-1)*nxp4**3 ! index in uloc_tot of ivar
               if(l<=z(nc-nc1+1)) then ! use interpolated value
                   uloc_tot(uindex) = u2(l+(nc-nc1)*maxz,i+1+2*j+4*k,ivar)
               else ! use value in uold
                   uloc_tot(uindex) = uold(cell_index,ivar)
               end if
           end do  
           do idim=1,ndim
               gindex = (gcount0(nc-nc1+1)-1) + pos + (idim-1)*nxp4**3
               if(poisson) then
                 if(l<=z(nc-nc1+1)) then ! use the index of the father
                     gloc_tot(gindex) = f(index_list(l,nc-nc1+1),idim) 
                 else   ! use the index of the cell
                     gloc_tot(gindex) = f(cell_index,idim)
                 end if
               else
                 gloc_tot(gindex) = 0.0d0
               end if
            end do
        end do
        end do
        end do
     end do
  end do
  !$acc end parallel loop
  
  !$acc end data
  !$acc end data
  
  deallocate(u2)
  
end subroutine fill_hydro_grids_interpol_2
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine fill_hydro_grids_interpol_3(nxp,ilevel,dx,dx_box,nxp4,nc1,nc2,icube,ngrid_ok,isort,ok,box_xmin,index_list,grid_list,bindex,z,maxz)
  use amr_commons
  use hydro_commons
  use poisson_commons
  use acc_commons
  implicit none
  integer::ilevel,nxp,nxp4,nc1,nc2,ngrid_ok,icube,maxz
  integer,dimension(active(ilevel)%ngrid) :: isort
  real(dp)::dx,dx_box
  
  real(dp),dimension(1:ndim,nc1:nc2):: box_xmin
  integer, dimension(1:(nxp/2+2)**3,1:ncube)::index_list,grid_list,bindex
  integer, dimension(1:ncube)::z
  logical, dimension(-1:nxp+2,-1:nxp+2,-1:nxp+2,nc1:nc2)::ok
  !--------------------------------------------------------------------------
  ! This routine fills up a Cartesian grid from AMR data.
  ! It is almost identical to fill_hydro_grids_interpol_2. The only
  ! difference is that here the interpolated values are not stored in the 
  ! temporary variable u2 but directly in uloc_tot. In this way less memory
  ! is used but there's a small overhead in the interpolation routine due to 
  ! the computation of the indeces of uloc_tot and the access of non 
  ! contiguous memory.
  !--------------------------------------------------------------------------
  ! In a test with ncube=500 it turned out that fill_hydro_grids_interpol_3
  ! is 1.025 times slower than fill_hydro_grids_interpol_2. 
  !--------------------------------------------------------------------------
  integer::i,j,k,l,idim,ivar,nc,pos,ic,jc,kc
  integer::cell_index,cell_levl,cell_grid
  integer :: uindex, gindex
  integer, dimension(0:twondim)::ibuffer_father
  real(dp),dimension(1:ndim):: xx_dp
  
  ! Private variables for getnborfather_acc
  integer,dimension(1:3)::ibound,iskip1,iskip2
  integer,dimension(1:3)::ix
  integer,dimension(0:twondim)::igridn
  integer,dimension(1:twondim)::icelln_ok 
  
  !$acc data present(igrid_acc,isort,xg,son,ucount0,gcount0,dx,dx_box,ilevel,nxp4,uold,uloc_tot,gloc_tot,f,ok,box_xmin) &
  !$acc present(index_list,grid_list,bindex,z,maxz)
 
  !$acc parallel loop
  do nc=nc1,nc2
     ! Index of one of the octs in the patch.
     igrid_acc(nc-nc1+1)=isort(icube+ngrid_ok*(nc-nc1))
     ! Computation of the starintg index of the patch
     ucount0(nc-nc1+1) = 1 + (nc-nc1)*nvar*nxp4**3
     gcount0(nc-nc1+1) = 1 + (nc-nc1)*ndim*nxp4**3
  end do
  !$acc end parallel loop
   
  ! Compute box coordinates in normalized unites
  !$acc parallel loop collapse(2)
  do nc=nc1,nc2
   do idim=1,ndim
      box_xmin(idim,nc)=int(xg(igrid_acc(nc-nc1+1),idim)/dx_box)*dx_box ! corner of the box
   end do
  end do
  !$acc end parallel loop
  
  
  !$acc parallel loop collapse(4) gang vector private(xx_dp)
  do nc=nc1,nc2
   do k=0,nxp/2+1
#if NDIM>1
      do j=0,nxp/2+1
#endif
#if NDIM>2
         do i=0,nxp/2+1
#endif
            ! Computing the coordinate of the father cell
            xx_dp(1) = box_xmin(1,nc) + (dble(i)-0.5)*2.0d0*dx
            xx_dp(2) = box_xmin(2,nc) + (dble(j)-0.5)*2.0d0*dx
            xx_dp(3) = box_xmin(3,nc) + (dble(k)-0.5)*2.0d0*dx         
              
            call hydro_get_cell_index_scalar(cell_index,cell_grid,cell_levl,xx_dp,ilevel-1)
  
            pos = i+1 + j*(nxp/2+2) + k*(nxp/2+2)**2 ! index of the father cell in the following lists
            index_list(pos,nc-nc1+1) = cell_index    ! index of the father cell
            grid_list(pos,nc-nc1+1) = cell_grid      ! index of the grid
            bindex(pos,nc-nc1+1) =  1 + 2*i + 2*nxp4*j + 2*nxp4*nxp4*k ! index in the patch of the first cell of the grid
            
            if(cell_grid==0) then ! the father cell has no son...
               do kc=0,1
               do jc=0,1
               do ic=0,1
                  ok(2*i+ic-1,2*j+jc-1,2*k+kc-1,nc)=.false. !...so the cells neither
               end do
               end do
               end do
            else ! the father cell has a son
               do kc=0,1
               do jc=0,1
               do ic=0,1
                  cell_index = ncoarse+(ic+2*jc+4*kc)*ngridmax+cell_grid ! index of the cells
                  ok(2*i+ic-1,2*j+jc-1,2*k+kc-1,nc) = son(cell_index)>0  ! check if they are refined
               end do
               end do
               end do
            end if 
            
#if NDIM>2
         end do
#endif
#if NDIM>1
      end do
#endif
   end do
  
  end do  
  !$acc end parallel loop
  
  !$acc kernels
  maxz = 0.0d0
  !$acc end kernels
    
  ! Putting the father cells which need interpolation at the beginning of the list
  !$acc parallel loop gang vector reduction(max:maxz)
  do nc=nc1,nc2     
     ! putting the zeros in grid_list(:,nc) at the beginning of the list. index_list and bindex will be reordered
     ! z will contain the number of zeros found, so the number of cells to be interpolated
     call sort_zeros(grid_list(:,nc-nc1+1), index_list(:,nc-nc1+1), bindex(:,nc-nc1+1), (nxp/2+2)**3, z(nc-nc1+1))
     maxz = max(maxz,z(nc-nc1+1))
  end do
  !$acc end parallel loop
  
  !$acc update host(maxz)  
  
  ! Storing the variables needed for the interpolation
  !$acc parallel loop collapse(2) gang vector &
  !$acc present(ggg_acc,hhh_acc,nbor) private(ibuffer_father,ibound,iskip1,iskip2,ix,igridn,icelln_ok)
  do nc=nc1,nc2
     do l=1,maxz
        
        if(l<=z(nc-nc1+1)) then
         call getnborfather_acc(index_list(l,nc-nc1+1),ibuffer_father,ilevel,&
                               &ibound,iskip1,iskip2,ix,igridn,icelln_ok)

         do ivar=1,nvar
            do j=0,twondim
               uloc_tot(l+(nc-nc1)*maxz + j*(nc2-nc1+1)*maxz + (ivar-1)*(twondim+1)*(nc2-nc1+1)*maxz)=uold(ibuffer_father(j),ivar)
            end do
         end do
        end if
        
     end do
  end do
  !$acc end parallel loop
  
  ! Parallel interpolation. The new values will be in uloc_tot at the right place.
  call interpol_hydro_acc_bis(uloc_tot(1:+nvar*(twondim+1)*(nc2-nc1+1)*maxz),maxz*(nc2-nc1+1),nc1,nc2,maxz,z,bindex,nxp4)
  
  ! Filling uloc_tot and gloc_tot
  !$acc parallel loop collapse(5) gang vector
  do nc=nc1,nc2 ! loop over patches
     do l=1,(nxp/2+2)**3  ! loop over father cells
        do k=0,1  ! loops over cells
        do j=0,1
        do i=0,1
           pos = bindex(l,nc-nc1+1) + i + j*nxp4 + k*nxp4*nxp4 ! index in the patch of the cell
           cell_index = ncoarse + (i+2*j+4*k)*ngridmax + grid_list(l,nc-nc1+1)
           do ivar=1,nvar
               uindex = (ucount0(nc-nc1+1)-1) + pos + (ivar-1)*nxp4**3 ! index in uloc_tot of ivar
               if(l>z(nc-nc1+1)) then ! use value in uold
                   uloc_tot(uindex) = uold(cell_index,ivar)
               end if
           end do  
           do idim=1,ndim
               gindex = (gcount0(nc-nc1+1)-1) + pos + (idim-1)*nxp4**3
               if(poisson) then
                 if(l<=z(nc-nc1+1)) then ! use the index of the father
                     gloc_tot(gindex) = f(index_list(l,nc-nc1+1),idim)
                 else ! use the index of the cell
                     gloc_tot(gindex) = f(cell_index,idim)
                 end if
               else
                 gloc_tot(gindex) = 0.0d0
               end if
            end do
        end do
        end do
        end do
     end do
  end do
  !$acc end parallel loop
  
  !$acc end data
  
end subroutine fill_hydro_grids_interpol_3
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine solve_hydro_grids_upd_flux(nxp,ilevel,dx,dx_box,dx_loc,nxp4,nc1,nc2,ok,box_xmin,flux,tmp)
  use amr_commons
  use hydro_commons
  use poisson_commons
  use acc_commons
  implicit none
  integer::ilevel,nxp,nxp4,nc1,nc2
  real(dp)::dx,dx_loc,dx_box
  logical, dimension(-1:nxp+2,-1:nxp+2,-1:nxp+2,nc1:nc2)::ok
  real(dp),dimension(1:nxp+1,1:nxp+1,1:nxp+1,1:nvar,1:ndim,nc1:nc2)::flux
  real(dp),dimension(1:nxp+1,1:nxp+1,1:nxp+1,1:2   ,1:ndim,nc1:nc2)::tmp
  real(dp),dimension(1:ndim,nc1:nc2)::box_xmin
  !--------------------------------------------------------------------------
  ! This routine computes the fluxes across the cells of the patches by
  ! calling the unsplit_gpu routine. Then it set to zero the fluxes at
  ! refined interface. Finally it updates unew at level ilevel and ilevel-1.
  !--------------------------------------------------------------------------
  real(dp),dimension(1:ndim) :: xx_dp
  integer :: i,j,k,l,idim,ivar,uindex,cell_index,cell_grid,cell_levl,i0,j0,k0,nc,ic,jc,kc,cell,m,p,q,pc,qc
    
!$acc data present(tmp,flux) &
!$acc present(unew,uloc_tot,gloc_tot,xg,divu_acc,enew_acc) &
!$acc present(ucount0,ok,box_xmin,son,nbor,ilevel,nxp,nxp4,oneontwotondim)
  
  !-----------------------------------------------
  ! Compute flux using second-order Godunov method
  !-----------------------------------------------
  call unsplit_gpu(uloc_tot,gloc_tot, &
                    & dx_loc,nxp,dtnew(ilevel),nc2-nc1+1,ncube,flux,tmp) 
                     
  !------------------------------------------------
  ! Reset flux along direction at refined interface    
  !------------------------------------------------
  !$acc parallel loop collapse(5) gang vector
  do nc=nc1,nc2
  do idim=1,ndim
  do k=1,nxp+1
  do j=1,nxp+1
  do i=1,nxp+1
     i0=0; j0=0; k0=0;
     if(idim==1) i0=1
     if(idim==2) j0=1
     if(idim==3) k0=1
     if(ok(i-i0,j-j0,k-k0,nc) .or. ok(i,j,k,nc)) then
        flux(i,j,k,:,idim,nc)=0.0d0
        if(pressure_fix) tmp(i,j,k,:,idim,nc)=0.0d0
     end if
  end do
  end do
  end do
  end do
  end do
  !$acc end parallel loop
  
  
  !--------------------------------------
  ! Conservative update at level ilevel
  !--------------------------------------
  !$acc parallel loop collapse(4) gang vector private(xx_dp)
  do nc=nc1,nc2
   do k=1,nxp
#if NDIM>1
      do j=1,nxp
#endif
#if NDIM>2
         do i=1,nxp           
#endif
            xx_dp(1) = box_xmin(1,nc) + (dble(i)-0.5)*dx
            xx_dp(2) = box_xmin(2,nc) + (dble(j)-0.5)*dx
            xx_dp(3) = box_xmin(3,nc) + (dble(k)-0.5)*dx
            call hydro_get_cell_index_scalar(cell_index,cell_grid,cell_levl,xx_dp,ilevel)
            
            do ivar=1,nvar
#if NDIM==1
            uindex = (ucount0(nc-nc1+1)-1) + (i+2) + (ivar-1)*nxp4
#endif
#if NDIM==2
            uindex = (ucount0(nc-nc1+1)-1) + (i+2) + (j+1)*nxp4 + (ivar-1)*nxp4**2
#endif
#if NDIM==3
            uindex = (ucount0(nc-nc1+1)-1) + ((i+2) + (j+1)*nxp4 + (k+1)*nxp4*nxp4) + (ivar-1)*nxp4**3
#endif
            do idim=1,ndim
               i0=0; j0=0; k0=0
               if(idim==1)i0=1
               if(idim==2)j0=1
               if(idim==3)k0=1
               
               ! Update conservative variables new state vector
               unew(cell_index,ivar)  = unew(cell_index,ivar) + (flux(i,j,k,ivar,idim,nc)-flux(i+i0,j+j0,k+k0,ivar,idim,nc))
            end do

            end do
            
           ! Update velocity divergence and internal energy
           if(pressure_fix)then
              do idim=1,ndim
                  i0=0; j0=0; k0=0
                  if(idim==1)i0=1
                  if(idim==2)j0=1
                  if(idim==3)k0=1
#if NDIM==1
                  divu_acc(cell_index) = divu_acc(cell_index)+(tmp(i,1,idim,nc)-tmp(i+i0,1,idim,nc))
                  enew_acc(cell_index) = enew_acc(cell_index)+(tmp(i,2,idim,nc)-tmp(i+i0,2,idim,nc))
#endif
#if NDIM==2
                  divu_acc(cell_index) = divu_acc(cell_index)+(tmp(i,j,1,idim,nc)-tmp(i+i0,j+j0,1,idim,nc))
                  enew_acc(cell_index) = enew_acc(cell_index)+(tmp(i,j,2,idim,nc)-tmp(i+i0,j+j0,2,idim,nc))
#endif
#if NDIM==3
                  divu_acc(cell_index) = divu_acc(cell_index)+(tmp(i,j,k,1,idim,nc)-tmp(i+i0,j+j0,k+k0,1,idim,nc))
                  enew_acc(cell_index) = enew_acc(cell_index)+(tmp(i,j,k,2,idim,nc)-tmp(i+i0,j+j0,k+k0,2,idim,nc))
#endif
              end do
            end if
            
#if NDIM>2
         end do
#endif
#if NDIM>1
      end do
#endif
   end do
   
  end do
  !$acc end parallel loop
  
  
  !--------------------------------------
  ! Conservative update at level ilevel-1
  !--------------------------------------
  ! ATTENTION: The parallel loop can't be on the outside loop due to race condition. 
  do l=1,twondim ! Loop over faces of the patch
  !$acc parallel loop collapse(3) gang vector private(xx_dp,idim)
  do nc=nc1,nc2  ! Loop over patches
     do p=1,nxp/2! Loops over grids on the face
     do q=1,nxp/2
        
        i0=0; j0=0; k0=0
        if(l==1 .or. l==2) then
           idim=1
           i0=1
        endif
        if(l==3 .or. l==4) then
           idim=2
           j0=1
        end if
        if(l==5 .or. l==6) then
           idim=3
           k0=1
        end if
        m = 1-mod(l,2)
        
        xx_dp(1) = box_xmin(1,nc) + ((dble(p)-0.5)*j0 + (dble(q)-0.5)*k0 + ((1-m)*(0.0-0.5)+m*(nxp/2+1-0.5))*i0)*2.0d0*dx
        xx_dp(2) = box_xmin(2,nc) + ((dble(p)-0.5)*k0 + (dble(q)-0.5)*i0 + ((1-m)*(0.0-0.5)+m*(nxp/2+1-0.5))*j0)*2.0d0*dx
        xx_dp(3) = box_xmin(3,nc) + ((dble(p)-0.5)*i0 + (dble(q)-0.5)*j0 + ((1-m)*(0.0-0.5)+m*(nxp/2+1-0.5))*k0)*2.0d0*dx
        
        call hydro_get_cell_index_scalar(cell_index,cell_grid,cell_levl,xx_dp,ilevel-1)
  
         if(son(cell_index)==0) then
            do ivar=1,nvar
               !$acc loop seq
               do pc=0,1
               do qc=0,1
                  ic = ((1-m)+m*(nxp+1))*i0 + (2*p-pc)*j0 + (2*q-qc)*k0
                  jc = ((1-m)+m*(nxp+1))*j0 + (2*p-pc)*k0 + (2*q-qc)*i0
                  kc = ((1-m)+m*(nxp+1))*k0 + (2*p-pc)*i0 + (2*q-qc)*j0
                  unew(cell_index,ivar)=unew(cell_index,ivar) &
                                     & +(2*m-1)*flux(ic,jc,kc,ivar,idim,nc)*oneontwotondim
               end do
               end do
               !$acc end loop
            end do
            if(pressure_fix) then
               !$acc loop seq
               do pc=0,1
               do qc=0,1
                  ic = ((1-m)+m*(nxp+1))*i0 + (2*p-pc)*j0 + (2*q-qc)*k0
                  jc = ((1-m)+m*(nxp+1))*j0 + (2*p-pc)*k0 + (2*q-qc)*i0
                  kc = ((1-m)+m*(nxp+1))*k0 + (2*p-pc)*i0 + (2*q-qc)*j0
                  divu_acc(cell_index) = divu_acc(cell_index) &
                             & +(2*m-1)*tmp(ic,jc,kc,1,idim,nc)*oneontwotondim
                  enew_acc(cell_index) = enew_acc(cell_index) &
                             & +(2*m-1)*tmp(ic,jc,kc,2,idim,nc)*oneontwotondim
               end do
               end do
               !$acc end loop
            end if
         end if
         
     end do
     end do
  end do
  !$acc end parallel loop
  end do

!$acc end data
  
end subroutine solve_hydro_grids_upd_flux
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine hydro_get_cell_index_scalar(cell_index,cell_grid,cell_levl,xpart,ilevel)
  use amr_commons
  implicit none
  integer::ilevel
  integer::cell_index,cell_levl,cell_grid
  real(dp),dimension(1:ndim)::xpart
  !----------------------------------------------------------------------------
  ! This routine returns the index of the cell, at maximum level
  ! ilevel, in which the input particle sits.
  ! Warning: coordinates are supposed to be in normalized units
  ! for the inner computational box so always between 0 and 1.
  !----------------------------------------------------------------------------
  real(dp)::xx,yy,zz
  integer::j,ii,jj,kk,ind,iskip,igrid,ind_cell,igrid0
  
     ind_cell=0
     ii=0; jj=0; kk=0
     xx=xpart(1)
     if(xx<0d0)xx=xx+dble(nx)
     if(xx>=dble(nx))xx=xx-dble(nx)
     ii=int(xx)
#if NDIM>1
     yy=xpart(2)
     if(yy<0d0)yy=yy+dble(ny)
     if(yy>=dble(ny))yy=yy-dble(ny)
     jj=int(yy)
#endif
#if NDIM>2
     zz=xpart(3)
     if(zz<0d0)zz=zz+dble(nz)
     if(zz>=dble(nz))zz=zz-dble(nz)
     kk=int(zz)
#endif
     igrid=son(1+ii+jj*nx+kk*nx*ny)
     do j=1,ilevel
        ii=0; jj=0; kk=0
        if(xx>xg(igrid,1))ii=1
#if NDIM>1
        if(yy>xg(igrid,2))jj=1
#endif
#if NDIM>2
        if(zz>xg(igrid,3))kk=1
#endif
        ind=1+ii+2*jj+4*kk
        iskip=ncoarse+(ind-1)*ngridmax
        ind_cell=iskip+igrid
        igrid=son(ind_cell)
        if(igrid==0.or.j==ilevel)exit
     end do
     cell_index=ind_cell
     cell_grid=igrid
     cell_levl=j

end subroutine hydro_get_cell_index_scalar
