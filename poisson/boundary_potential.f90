!################################################################
!################################################################
!################################################################
!################################################################
subroutine make_boundary_force(ilevel)
  use amr_commons
  use poisson_commons
  implicit none
  integer::ilevel
  ! -------------------------------------------------------------------
  ! This routine set up boundary conditions for fine levels.
  ! -------------------------------------------------------------------
  integer::ibound,boundary_dir,idim,inbor=1
  integer::i,ncache,ivar,igrid,ngrid,ind
  integer::iskip,iskip_ref,nx_loc,ix,iy,iz
  integer,dimension(1:8)::ind_ref
  integer,dimension(1:nvector),save::ind_grid,ind_grid_ref
  integer,dimension(1:nvector),save::ind_cell,ind_cell_ref

  real(dp)::switch,dx,dx_loc,scale
  real(dp),dimension(1:3)::gs,skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp),dimension(1:nvector,1:ndim),save::ff

  if(.not. simple_boundary)return
  if(verbose)write(*,111)ilevel

  ! Mesh size at level ilevel
  dx=0.5D0**ilevel

  ! Rescaling factors
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=[0, 0, 0]
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

     ! Vector sign switch for reflexive boundary conditions
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

           ! Wall and free boundary conditions
           if((boundary_type(ibound)/10).ne.2)then

              ! Gather reference hydro variables
              do ivar=1,ndim
                 do i=1,ngrid
                    ff(i,ivar)=f(ind_cell_ref(i),ivar)
                 end do
              end do
              ! Scatter to boundary region
              do ivar=1,ndim
                 switch=gs(ivar)
                 do i=1,ngrid
                    f(ind_cell(i),ivar)=ff(i,ivar)*switch
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

              call gravana(xx,ff,dx_loc,ngrid)

              ! Scatter variables
              do ivar=1,ndim
                 do i=1,ngrid
                    f(ind_cell(i),ivar)=ff(i,ivar)
                 end do
              end do

           end if

        end do
        ! End loop over cells

     end do
     ! End loop over grids

  end do
  ! End loop over boundaries

111 format('   Entering make_boundary_force for level ',I2)

end subroutine make_boundary_force
!##########################################################################
!##########################################################################
!##########################################################################
!##########################################################################
subroutine make_boundary_phi(ilevel)
  use amr_commons
  use poisson_commons
  use constants, only: twopi
  implicit none
  integer::ilevel
  ! -------------------------------------------------------------------
  ! This routine set up boundary conditions for fine levels.
  ! -------------------------------------------------------------------
  integer::ibound,idim
  integer::i,ncache,igrid,ngrid,ind
  integer::iskip,nx_loc,ix,iy,iz
  integer,dimension(1:nvector),save::ind_grid,ind_cell

  real(dp)::dx,dx_loc,scale,fourpi,boxlen2
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector),save::rr,pp
  real(dp),dimension(1:nvector,1:ndim),save::xx

  if(.not. simple_boundary)return
  if(verbose)write(*,111)ilevel

  ! Mesh size at level ilevel
  dx=0.5D0**ilevel
  fourpi = 2*twopi

  ! Rescaling factors
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=[0, 0, 0]
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  boxlen2=boxlen**2

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

     ! Loop over grids by vector sweeps
     ncache=boundary(ibound,ilevel)%ngrid
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=boundary(ibound,ilevel)%igrid(igrid+i-1)
        end do

        ! Loop over cells
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do

           ! Compute cell center in code units
           do idim=1,ndim
              do i=1,ngrid
                 xx(i,idim)=xg(ind_grid(i),idim)+xc(ind,idim)
              end do
           end do

           ! Rescale position from code units to user units
           rr(1:ngrid)=0
           do idim=1,ndim
              do i=1,ngrid
                 xx(i,idim)=(xx(i,idim)-skip_loc(idim))*scale
                 rr(i)=rr(i)+(xx(i,idim)-multipole(idim+1)/multipole(1))**2
              end do
           end do

           do i=1,ngrid
               rr(i)=sqrt(rr(i))
           end do

           if(ngrid>0) call phi_ana(rr,pp,ngrid)

           ! Scatter variables
           do i=1,ngrid
              phi(ind_cell(i))=pp(i)/scale
           end do

        end do
        ! End loop over cells

     end do
     ! End loop over grids

  end do
  ! End loop over boundaries

111 format('   Entering make_boundary_phi for level ',I2)

end subroutine make_boundary_phi
!##########################################################################
!##########################################################################
!##########################################################################
!##########################################################################
subroutine make_boundary_mask(ilevel)
  use amr_commons
  use poisson_commons
  implicit none
  integer::ilevel
  ! -------------------------------------------------------------------
  ! This routine set up boundary conditions for fine levels.
  ! -------------------------------------------------------------------
  integer::ibound
  integer::i,ncache,igrid,ngrid,ind
  integer::iskip
  integer,dimension(1:nvector),save::ind_grid,ind_cell

  if(.not. simple_boundary)return
  if(verbose)write(*,111)ilevel

  ! Loop over boundaries
  do ibound=1,nboundary

     ! Loop over grids by vector sweeps
     ncache=boundary(ibound,ilevel)%ngrid
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=boundary(ibound,ilevel)%igrid(igrid+i-1)
        end do

        ! Loop over cells
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do

           ! Set mask to -1d0
           do i=1,ngrid
              f(ind_cell(i),3)=-1
           end do

        end do
        ! End loop over cells

     end do
     ! End loop over grids

  end do
  ! End loop over boundaries

111 format('   Entering make_boundary_mask for level ',I2)

end subroutine make_boundary_mask
!##########################################################################
!##########################################################################
!##########################################################################
!##########################################################################
