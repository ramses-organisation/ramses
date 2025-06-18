! This file contains subroutines to help access neighboring cells or grids.
!
!   GETTING 3^NDIM SURROUNDING NEIGHBORS, i.e. 3x3x3 cube incl. diagonal
!
!   * get3cubefather: IN: cell, OUT: 3^ndim neighbor cells
!                     Obtain the neighboring cells of the cell.
!                     Calls get3cubepos.
!
!   * get3cubefather_grids: IN: cell, OUT: 2^ndim neighbor grids
!                           Obtain the grids that contain the neighboring cells.
!                           Calls get_grids_of_nbor_cells.

!   (*) get3cubepos: IN: grid & ind, OUT: 3^ndim neighbor cells
!                    Obtain the neighboring cells of a cell sitting at position
!                    ind in the input grid.
!                    Helper for get3cubefather. Calls get_grids_of_nbor_cells.

!   (*) get_grids_of_nbor_cells: IN: grid & ind, OUT: 2^ndim neighbor grids
!                                Obtain the grids that contain the neighbor cells
!                                of a cell sitting at position ind in the input grid.
!                                Helper for get3cubepos and get3cubefather_grids.
!
!   GETTING 2xNDIM DIRECT NEIGHBORS
!
!   * getnborfather: IN: cell, OUT: 2*ndim neighbor cells
!                    Obtain the direct neighboring cells of the cell.
!                    Calls getnborgrids and getnborcells.
!
!   * getnborgrids: IN: grid, OUT: 2*ndim neighbor grids
!                   Obtain the direct neighboring grids of a grid.
!                   The result includes the grid itself at
!                   position 0 in the array.
!                   Often called before getnborcells.
!
!   * getnborcells: IN: neighbor grids & ind, OUT: 2*ndim neighbor cells
!                   Obtain the 2*ndim direct neighboring cells
!                   of a cell at position ind in its grid, of which the neighbor
!                   grids are provided.
!                   The result contains only the neighbors, not
!                   the cell itself.
!                   The input for this function is the output of
!                   getnborgrids and ind.
!
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine get3cubefather(ind_cell,nbor_cells,ncell,ilevel)
  use amr_parameters, only:nx,ny,nz,ndim,ngridmax,nvector,threetondim
  use amr_commons, only:ncoarse
  implicit none
  integer,intent(in)::ncell,ilevel
  integer,dimension(1:nvector),intent(in)::ind_cell
  integer,dimension(1:nvector,1:threetondim),intent(out)::nbor_cells
  !------------------------------------------------------------------
  ! This subroutine determines the 3^ndim neighboring cells
  ! of the input cell ind_cell. This includes diagonal neighbors and
  ! the center input cell itself.
  ! For example in 3D, as output, nbor_cells contains the 3x3x3 cube
  ! of 27 cells around, and including, ind_cell.
  ! According to the refinement rule, they should be present anytime.
  !------------------------------------------------------------------
  integer::i,j,nxny,i1,j1,k1,ind,iok,iskip
  integer::i1min,i1max,j1min,j1max,k1min,k1max,ind_father
  integer,dimension(1:nvector),save::ix,iy,iz,iix,iiy,iiz
  integer,dimension(1:nvector),save::pos,ind_grid_father

  nxny=nx*ny

  if(ilevel==1)then  ! Easy...
     ! ncell is always 1 in this case
     ! (because currently only cubic domains are supported, meaning
     ! there is only 1 root cell)
     do i=1,ncell
        iz(i)=(ind_cell(i)-1)/nxny
        iy(i)=(ind_cell(i)-1-iz(i)*nxny)/nx
        ix(i)=(ind_cell(i)-1-iy(i)*nx-iz(i)*nxny)
     end do

     i1min=0; i1max=2
     j1min=0; j1max=0
     if(ndim > 1)j1max=2
     k1min=0; k1max=0
     if(ndim > 2)k1max=2

     ! Loop over 3^ndim neighboring cells by looking
     ! left, center and right of the root cell
     do k1=k1min,k1max
        iiz=iz
        if(ndim > 2)then
           do i=1,ncell
              iiz(i)=iz(i)+k1-1
              if(iiz(i) < 0   )iiz(i)=nz-1
              if(iiz(i) > nz-1)iiz(i)=0
           end do
        end if
        do j1=j1min,j1max
           iiy=iy
           if(ndim > 1)then
              do i=1,ncell
                 iiy(i)=iy(i)+j1-1
                 if(iiy(i) < 0   )iiy(i)=ny-1
                 if(iiy(i) > ny-1)iiy(i)=0
              end do
           end if
           do i1=i1min,i1max
              iix=ix
              if(ndim > 0)then
                 do i=1,ncell
                    iix(i)=ix(i)+i1-1
                    if(iix(i) < 0   )iix(i)=nx-1
                    if(iix(i) > nx-1)iix(i)=0
                 end do
              end if
              ind_father=1+i1+3*j1+9*k1
              do i=1,ncell
                 nbor_cells(i,ind_father)=1 &
                      & +iix(i) &
                      & +iiy(i)*nx &
                      & +iiz(i)*nxny
              end do
           end do
        end do
     end do

  else    ! else, more complicated...

     ! Get the cell's position in its grid, that is the
     ! index ind, between 1 and twotondim.
     do i=1,ncell
        pos(i)=(ind_cell(i)-ncoarse-1)/ngridmax+1  !integer devision
     end do

     ! Convert the cell's index to the index of the grid to which the cell belongs
     do i=1,ncell
        iskip=ncoarse+(pos(i)-1)*ngridmax
        ind_grid_father(i)=ind_cell(i)-iskip
     end do

     ! Using the grid index and cell position, get the neighbor cells
     call get3cubepos(ind_grid_father,pos,nbor_cells,ncell)

  end if

end subroutine get3cubefather
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine get3cubefather_grids(ind_cell,nbors_grids,ncell,ilevel)
  use amr_commons
  implicit none
  integer,intent(in)::ncell,ilevel
  integer,dimension(1:nvector),intent(in)::ind_cell
  integer,dimension(1:nvector,1:twotondim),intent(out)::nbors_grids
  !------------------------------------------------------------------
  ! This subroutine determines the two^ndim grids that contain the
  ! 3x3x3 neighboring cells of the cell ind_cell.
  ! According to the refinement rule,they should be present anytime.
  ! Example 1D: see get_grids_of_nbor_cells.
  !------------------------------------------------------------------
  integer::i,j,nxny,i1,j1,k1,ind,iok,iskip
  integer::i1min,i1max,j1min,j1max,k1min,k1max,ind_father
  integer,dimension(1:nvector),save::ix,iy,iz,iix,iiy,iiz
  integer,dimension(1:nvector),save::pos,ind_grid_father
  integer,dimension(1:nvector,1:twotondim),save::nbors_grids_ok

  nxny=nx*ny

  if(ilevel==1)then  ! Easy...

     do i=1,ncell
        iz(i)=(ind_cell(i)-1)/nxny
        iy(i)=(ind_cell(i)-1-iz(i)*nxny)/nx
        ix(i)=(ind_cell(i)-1-iy(i)*nx-iz(i)*nxny)
     end do

     i1min=0; i1max=2
     j1min=0; j1max=0
     if(ndim > 1)j1max=2
     k1min=0; k1max=0
     if(ndim > 2)k1max=2

     ! Loop over 2^ndim neighboring father grids
     do k1=k1min,k1max
        iiz=iz
        if(ndim > 2)then
           do i=1,ncell
              iiz(i)=iz(i)+2*k1-1
              if(iiz(i) < 0   )iiz(i)=nz-1
              if(iiz(i) > nz-1)iiz(i)=0
           end do
        end if
        do j1=j1min,j1max
           iiy=iy
           if(ndim > 1)then
              do i=1,ncell
                 iiy(i)=iy(i)+2*j1-1
                 if(iiy(i) < 0   )iiy(i)=ny-1
                 if(iiy(i) > ny-1)iiy(i)=0
              end do
           end if
           do i1=i1min,i1max
              iix=ix
              if(ndim > 0)then
                 do i=1,ncell
                    iix(i)=ix(i)+2*i1-1
                    if(iix(i) < 0   )iix(i)=nx-1
                    if(iix(i) > nx-1)iix(i)=0
                 end do
              end if
              ind_father=1+i1+2*j1+4*k1
              do i=1,ncell
                 nbors_grids(i,ind_father)=1 &
                      & +(iix(i)/2) &
                      & +(iiy(i)/2)*(nx/2) &
                      & +(iiz(i)/2)*(nxny/4)
              end do
           end do
        end do
     end do

  else    ! else, more complicated...

     ! Get the cell's position in its grid, that is the
     ! index ind, between 1 and twotondim.
     do i=1,ncell
        pos(i)=(ind_cell(i)-ncoarse-1)/ngridmax+1  !integer devision
     end do

     ! Convert the cell's index to the index of the grid to which the cell belongs
     do i=1,ncell
        iskip=ncoarse+(pos(i)-1)*ngridmax
        ind_grid_father(i)=ind_cell(i)-iskip
     end do

     ! Using the grid index and cell position to get neighbor grids
     call get_grids_of_nbor_cells(ind_grid_father,pos,nbors_grids,ncell)

  end if

end subroutine get3cubefather_grids
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine get3cubepos(ind_grid,ind,nbors_cells,ng)
  use amr_commons
  use amr_constants, only:lll,mmm
  implicit none
  integer,intent(in)::ng
  integer,dimension(1:nvector),intent(in)::ind_grid,ind
  integer,dimension(1:nvector,1:threetondim)::nbors_cells
  !--------------------------------------------------------------------
  ! This subroutine determines the 3^ndim neighboring cells
  ! of the input cell at position ind in grid ind_grid.
  ! According to the refinements rules and since the input cell is refined,
  ! they should be present anytime.
  !
  ! example 1D:
  !
  !   input:            || IND_GRID ||
  !                       containing
  !              ||  cell=ind  | another cell ||
  !
  !   The grids of that contain the neighbors are determined by
  !   the subroutine get_grids_of_nbor_cells:
  !
  !     ||   NBOR GRID   ||    IND_GRID   ||
  !     || cell |  nbor  ||  IND  |  nbor ||
  !
  !   Then, to extract the cells from the grids, we use the indices from
  !   lll and mmm
  !
  !--------------------------------------------------------------------
  integer::i,j,iskip
  integer::icell,igrid
  integer,dimension(1:nvector,1:twotondim),save::nbors_grids
  integer,dimension(1:twotondim),save::lll_loc,mmm_loc

  ! Get the grids that contain all neighbor cells
  call get_grids_of_nbor_cells(ind_grid,ind,nbors_grids,ng)

  ! Extract each of the 3x3x3 neighbor cells out of the grids
  do j=1,threetondim

     ! Fetch magic indices
     lll_loc = lll(j,:)
     mmm_loc = mmm(j,:)

     do i=1,ng
        ! index in the list of nbor_grids of the grid in which the j-th neighbor cell should be located
        igrid = lll_loc(ind(i))
        ! position in which the neighbor cell should be located in its grid
        icell = mmm_loc(ind(i))

        iskip=ncoarse+(icell-1)*ngridmax

        ! If the grid for neighbor cell j exists,
        ! determine its index based on the index of the grid in which is sits
        ! otherwise set 0.
        nbors_cells(i,j) = merge(iskip+nbors_grids(i,igrid), & ! if-part
                              &                           0, & ! else-part
                              & (nbors_grids(i,igrid)>0))      ! condition
     end do
  end do

end subroutine get3cubepos
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine get_grids_of_nbor_cells(ind_grid,ind,nbors_grids,ng)
  use amr_commons
  implicit none
  integer,intent(in)::ng
  integer,dimension(1:nvector),intent(in)::ind_grid,ind
  integer,dimension(1:nvector,1:twotondim),intent(out)::nbors_grids
  !--------------------------------------------------------------------
  ! This subroutine determines the two^ndim grids that contain the
  ! 3x3x3 neighboring cells of a cell at position ind in the grid ind_grid.
  ! According to the refinements rules and since the input cell is refined,
  ! they should be present anytime.
  !
  ! example in 1D:
  !
  !   input:     ||       IND_GRID        ||
  !              || cell=IND | other cell ||
  !
  !   for IND_GRID we have access to the father cell (through father(grid))
  !   and the neighboring cells of the father cell (through nbor(grid,:)):
  !
  !      || another cell | neighbor f. cell || father cell | neighbor f. cell ||
  !                                         ||  IND_GRID  ||
  !
  !   For those neighbor cells, we can get the child grid (through son(neighbor_cell)):
  !
  !      || another cell | neighbor f. cell || father cell | neighbor f. cell ||
  !                     || son grid of nbor ||  IND_GRID  || son grid of nbor ||
  !
  !   We need only the grids which contain the neighbor cells. This will be ind_grid
  !   itself and either its left or right neighbor:
  !
  !                     || SON GRID OF NBOR ||  IND_GRID  || son grid of nbor ||
  !                     ||  cell   | nbor   || ind  | nbor||  cell  |   cell  ||
  !
  !   In this example, the requested grids are the left and middle grid (in capital letters).
  !--------------------------------------------------------------------
  integer::i,inbor
  integer::ii,iimin,iimax
  integer::jj,jjmin,jjmax
  integer::kk,kkmin,kkmax
  integer,dimension(1:8)::iii=(/1,2,1,2,1,2,1,2/)
  integer,dimension(1:8)::jjj=(/3,3,4,4,3,3,4,4/)
  integer,dimension(1:8)::kkk=(/5,5,5,5,6,6,6,6/)
  integer,dimension(1:nvector),save::ind_grid1,ind_grid2,ind_grid3

  iimin=0; iimax=1
  jjmin=0; jjmax=0
  if(ndim>1)jjmax=1
  kkmin=0; kkmax=0
  if(ndim>2)kkmax=1

  do kk=kkmin,kkmax
     do i=1,ng
        ind_grid1(i)=ind_grid(i)
     end do
     if(kk>0)then
        do i=1,ng
           if(ind_grid(i)>0)then
              ind_grid1(i)=son(nbor(ind_grid(i),kkk(ind(i))))
           endif
        end do
     end if

     do jj=jjmin,jjmax
        do i=1,ng
           ind_grid2(i)=ind_grid1(i)
        end do
        if(jj>0)then
           do i=1,ng
              if(ind_grid1(i)>0)then
                 ind_grid2(i)=son(nbor(ind_grid1(i),jjj(ind(i))))
              endif
           end do
        end if

        do ii=iimin,iimax
           do i=1,ng
              ind_grid3(i)=ind_grid2(i)
           end do
           if(ii>0)then
              do i=1,ng
                 if(ind_grid2(i)>0)then
                    ind_grid3(i)=son(nbor(ind_grid2(i),iii(ind(i))))
                 endif
              end do
           end if

           inbor=1+ii+2*jj+4*kk
           do i=1,ng
              nbors_grids(i,inbor)=ind_grid3(i)
           end do

        end do
     end do
  end do

end subroutine get_grids_of_nbor_cells
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine getnborcells(igridn,ind,icelln,ng)
  use amr_commons
  use amr_constants, only:iii,jjj
  implicit none
  integer,intent(in)::ng,ind
  integer,dimension(1:nvector,0:twondim),intent(in)::igridn
  integer,dimension(1:nvector,1:twondim),intent(out)::icelln
  !--------------------------------------------------------------
  ! This routine computes the index of 6-neighboring cells (non-diagonal).
  ! The user must provide nbors_grids = index of the 6 neighboring
  ! grids with the cell's grid at position 0 in the list,
  ! as returned by the routine getnborgrids.
  ! ind is the position index of the center cell in its grid.
  !--------------------------------------------------------------
  integer::i,in,igrid,icell,iskip,inbor,idim

  ! Reset indices
  icelln(1:ng,1:twondim)=0

  ! Extract, out of the input grids, the left and right
  ! neighbor cell in each dimension.
  do inbor=1,2
     do idim=1,ndim
        ! index for neighbor order
        in=(idim-1)*2 + inbor
        ! index in the input list nbors_grids, of the grid in which the j-th neighbor cell should be located
        igrid=iii(idim,inbor,ind)
        ! position in which the neighbor cell should be located in its grid
        icell=jjj(idim,inbor,ind)

        iskip=ncoarse+(icell-1)*ngridmax
        do i=1,ng
           if(igridn(i,igrid)>0)then
              icelln(i,in)=iskip+igridn(i,igrid)
           end if
        end do
     enddo
  end do

end subroutine getnborcells
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine getnborfather(ind_cell,ind_father,ncell,ilevel)
  use amr_commons
  implicit none
  integer::ncell,ilevel
  integer,dimension(1:nvector)::ind_cell
  integer,dimension(1:nvector,0:twondim)::ind_father
  !-----------------------------------------------------------------
  ! This subroutine determines the 2*ndim neighboring cells
  ! of the input cell (ind_cell).
  ! If for some reasons they don't exist, the routine returns
  ! the neighboring father cells of the input cell.
  !-----------------------------------------------------------------
  integer::nxny,i,idim,j,iok,ind,iskip
  integer,dimension(1:3)::ibound,iskip1,iskip2
  integer,dimension(1:nvector,1:3),save::ix
  integer,dimension(1:nvector),save::ind_grid_father,pos
  integer,dimension(1:nvector,0:twondim),save::igridn,igridn_ok
  integer,dimension(1:nvector,1:twondim),save::icelln_ok

  nxny=nx*ny

  if(ilevel==1)then

     ibound(1)=nx-1
     iskip1(1)=1
     iskip2(1)=nx-1
     ibound(2)=ny-1
     iskip1(2)=nx
     iskip2(2)=(ny-1)*nx
     ibound(3)=nz-1
     iskip1(3)=nxny
     iskip2(3)=(nz-1)*nxny

     ! Get father cell
     do i=1,ncell
        ind_father(i,0)=ind_cell(i)
     end do

     do i=1,ncell
        ix(i,3)=(ind_father(i,0)-1)/nxny
        ix(i,2)=(ind_father(i,0)-1-ix(i,3)*nxny)/nx
        ix(i,1)=(ind_father(i,0)-1-ix(i,2)*nx-ix(i,3)*nxny)
     end do

     do idim=1,ndim
        do i=1,ncell
           if(ix(i,idim)>0)then
              ind_father(i,2*idim-1)=ind_father(i,0)-iskip1(idim)
           else
              ind_father(i,2*idim-1)=ind_father(i,0)+iskip2(idim)
           end if
        end do
        do i=1,ncell
           if(ix(i,idim)<ibound(idim))then
              ind_father(i,2*idim)=ind_father(i,0)+iskip1(idim)
           else
              ind_father(i,2*idim)=ind_father(i,0)-iskip2(idim)
           end if
        end do
     end do

  else

     ! Get father cell
     do i=1,ncell
        ind_father(i,0)=ind_cell(i)
     end do

     ! Get the cell's position in its grid, that is the
     ! index ind, between 1 and twotondim.
     do i=1,ncell
        pos(i)=(ind_cell(i)-ncoarse-1)/ngridmax+1  !integer devision
     end do

     ! Convert the cell's index to the index of the grid to which the cell belongs
     do i=1,ncell
        iskip=ncoarse+(pos(i)-1)*ngridmax
        ind_grid_father(i)=ind_cell(i)-iskip
     end do

     ! Get the 2**ndim grids that contain the neighboring cells
     call getnborgrids(ind_grid_father,igridn,ncell)

     ! Loop over position
     do ind=1,twotondim

        ! Select father cells that sit at position ind
        do j=0,twondim
           iok=0
           do i=1,ncell
              if(pos(i)==ind)then
                 iok=iok+1
                 igridn_ok(iok,j)=igridn(i,j)
              end if
           end do
        end do

        ! Get neighboring cells for selected cells
        if(iok>0)call getnborcells(igridn_ok,ind,icelln_ok,iok)

        ! Update neighboring father cells for selected cells
        do j=1,twondim
           iok=0
           do i=1,ncell
              if(pos(i)==ind)then
                 iok=iok+1
                 if(icelln_ok(iok,j)>0)then
                    ind_father(i,j)=icelln_ok(iok,j)
                 else
                    ind_father(i,j)=nbor(ind_grid_father(i),j)
                 end if
              end if
           end do
        end do

     end do

  end if

end subroutine getnborfather
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine getnborgrids(igrid,igridn,ngrid)
  use amr_commons
  implicit none
  integer,intent(in)::ngrid
  integer,dimension(1:nvector),intent(in)::igrid
  integer,dimension(1:nvector,0:twondim),intent(out)::igridn
  !---------------------------------------------------------
  ! This routine computes the index of the 2 x ndim neighboring
  ! grids for grid igrid(:). The index for the central
  ! grid is stored in igridn(:,0).  We need it for when we want
  ! to later derive the neighbor cells from these grids (see
  ! getnborcells).
  ! If for some reasons the neighboring grids don't exist,
  ! then igridn(:,j) = 0 (because son is 0, i.e. it's a leaf cell).
  !---------------------------------------------------------
  integer::i,j

  ! Store central grid at position 0
  do i=1,ngrid
     igridn(i,0)=igrid(i)
  end do
  ! Store neighboring grids
  do j=1,twondim
     do i=1,ngrid
        ! nbor(igrid(i),j) is the jth neighbor cell of the father cell
        ! of the grid with index igrid
        ! Son points to the child grid of that cell.
        igridn(i,j)=son(nbor(igrid(i),j))
     end do
  end do

end subroutine getnborgrids
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine getnborgrids_check(igrid,igridn,ngrid)
  use amr_commons
  implicit none
  integer::ngrid
  integer,dimension(1:nvector)::igrid
  integer,dimension(1:nvector,0:twondim)::igridn
  !---------------------------------------------------------
  ! This routine does EXACTLY the same as getnborgrids
  ! but it checks if the neighbor of igrid exists in order
  ! to avoid son(0) leading to crash.
  !---------------------------------------------------------
  integer::i,j

  ! Store central grid
  do i=1,ngrid
     igridn(i,0)=igrid(i)
  end do
  ! Store neighboring grids
  do j=1,twondim
     do i=1,ngrid
        if (nbor(igrid(i),j)>0)igridn(i,j)=son(nbor(igrid(i),j))
     end do
  end do

end subroutine getnborgrids_check
