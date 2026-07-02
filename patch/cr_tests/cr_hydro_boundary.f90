!################################################################
!################################################################
!################################################################
!################################################################
subroutine cr_make_boundary_hydro(ilevel)
  use amr_commons
  use cr_hydro_commons
  use cr_parameters
  implicit none
  integer::ilevel
  ! -------------------------------------------------------------------
  ! Boundary conditions for the separated CR arrays cruold(:,1:ncrvar). Gas
  ! rho/v/B/MHD bookkeeping is handled by the gas make_boundary_hydro, not here.
  ! Reflexive walls flip the sign of each CR flux component normal to the
  ! boundary (energy unchanged). Periodic boundaries never reach this routine
  ! (served by the cruold virtual-cell exchange).
  ! -------------------------------------------------------------------
  integer::ibound,boundary_dir,idim,inbor
  integer::i,ncache,ivar,igrid,ngrid,ind
  integer::iskip,iskip_ref,nx_loc,ix,iy,iz
  integer::crType
  integer,dimension(1:8)::ind_ref
  integer,dimension(1:nvector),save::ind_grid,ind_grid_ref
  integer,dimension(1:nvector),save::ind_cell,ind_cell_ref

  real(dp)::switch,dx,dx_loc,scale
  real(dp),dimension(1:3)::gs,skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp),dimension(1:nvector,1:ncrvar),save::uu

  if(.not. simple_boundary)return
  if(.not. cr_advect)return
  if(verbose)write(*,111)ilevel

  ! Mesh size at level ilevel
  dx=0.5D0**ilevel

  ! Rescaling factors
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  ! Loop over boundaries
  do ibound=1,nboundary

     call set_boundary_references(ibound,ind_ref,boundary_dir,inbor)

     ! Velocity sign switch for reflexive boundary conditions (applied to
     ! the CR flux component normal to the boundary).
     gs=(/1,1,1/)
     if(boundary_type(ibound)==1.or.boundary_type(ibound)==2)gs(1)=-1
     if(boundary_type(ibound)==3.or.boundary_type(ibound)==4)gs(2)=-1
     if(boundary_type(ibound)==5.or.boundary_type(ibound)==6)gs(3)=-1

     ! Loop over grids by vector sweeps
     ncache=boundary(ibound,ilevel)%ngrid
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=boundary(ibound,ilevel)%igrid(igrid+i-1)
        end do

        ! Gather neighboring reference grid
        do i=1,ngrid
           ind_grid_ref(i)=son(nbor(ind_grid(i),inbor))
        end do

        ! Loop over cells
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do

           ! Gather neighboring reference cell
           iskip_ref=ncoarse+(ind_ref(ind)-1)*ngridmax
           do i=1,ngrid
              ind_cell_ref(i)=iskip_ref+ind_grid_ref(i)
           end do

           ! Wall or reflexive boundary conditions
           if((boundary_type(ibound)/10).eq.0)then

              ! Gather reference CR variables
              do ivar=1,ncrvar
                 do i=1,ngrid
                    uu(i,ivar)=cruold(ind_cell_ref(i),ivar)
                 end do
              end do

              ! Scatter to boundary region, flipping CR flux normal to wall.
              ! crType=0 -> CR energy (no flip); 1..ndim -> flux component.
              do ivar=1,ncrvar
                 switch=1
                 crType=mod(ivar-Ecr_idx(1),ndim+1)
                 if(crType.ne.0)switch=gs(crType)
                 do i=1,ngrid
                    cruold(ind_cell(i),ivar)=uu(i,ivar)*switch
                 end do
              end do

              ! CR-test reflexive-wall energy pin (Jiang & Oh 413/424): pin
              ! group-1 CR energy while the wall still reflects the flux.
              if(trim(jiang_test)=='413' .or. trim(jiang_test)=='424')then
                 do i=1,ngrid
                    cruold(ind_cell(i),Ecr_idx(1))=3.0d0
                 end do
              end if

           ! Free or outflowing or zero gradient boundary conditions
           else if((boundary_type(ibound)/10).eq.1)then

              ! Gather reference CR variables
              do ivar=1,ncrvar
                 do i=1,ngrid
                    uu(i,ivar)=cruold(ind_cell_ref(i),ivar)
                 end do
              end do

              ! Scatter to boundary region (zero gradient: plain copy)
              do ivar=1,ncrvar
                 do i=1,ngrid
                    cruold(ind_cell(i),ivar)=uu(i,ivar)
                 end do
              end do

           ! Imposed boundary conditions
           else

              ! Compute cell center in code units
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

              call cr_boundana(xx,uu,dx_loc,ibound,ngrid)

              ! Scatter variables
              do ivar=1,ncrvar
                 do i=1,ngrid
                    cruold(ind_cell(i),ivar)=uu(i,ivar)
                 end do
              end do

           end if

        end do
        ! End loop over cells

     end do
     ! End loop over grids

  end do
  ! End loop over boundaries

111 format('   Entering cr_make_boundary_hydro for level ',I2)

end subroutine cr_make_boundary_hydro
