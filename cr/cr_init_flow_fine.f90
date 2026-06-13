!################################################################
!################################################################
!################################################################
!################################################################
subroutine cr_init_flow
  ! Initialise the separated cosmic-ray state arrays cruold(:,1:ncrvars)
  ! on refinement, mirroring rt/rt_init_flow_fine.f90 (rt_init_flow) but
  ! writing CR variables from cruold/crunew with iCRu=1 instead of the RT
  ! rtuold columns. Inert unless cr_advect.
  use amr_commons
  use cr_parameters
  use cr_hydro_commons, ONLY: cruold
  implicit none

  integer::ilevel,ivar

  if(.not.cr_advect)return
  if(verbose)write(*,*)'Entering cr_init_flow'
  do ilevel=nlevelmax,1,-1
     if(ilevel>=levelmin)call cr_init_flow_fine(ilevel)
     ! Restrict the SEPARATE CR array cruold(:,1:ncrvars) fine->coarse onto
     ! coarse split cells. The mirrored rt_init_flow calls rt_upload_fine here;
     ! the generic upload_fine restricts only the gas array uold (it 'use
     ! hydro_commons' and never touches cruold), so cr_upload_fine is the
     ! faithful CR analog (see also cr_interpol_hydro for prolongation).
     call cr_upload_fine(ilevel)
     do ivar=1,ncrvars
        call make_virtual_fine_dp(cruold(1,ivar),ilevel)
     end do
     if(simple_boundary)call cr_make_boundary_hydro(ilevel)
  end do
  if(verbose)write(*,*)'Complete cr_init_flow'

end subroutine cr_init_flow
!################################################################
!################################################################
!################################################################
!################################################################
subroutine cr_init_flow_fine(ilevel)
  ! Per-level CR initialisation: gather new cell centre positions, call
  ! cr_condinit, and scatter the result into cruold(:,1:ncrvars). Modeled
  ! on rt_init_flow_fine's condinit path (CR has no ic_cr_* IC-file reader
  ! in phase 1, so only the built-in cr_condinit branch is used).
  use amr_commons
  use cr_parameters
  use cr_hydro_commons
  implicit none
  integer::ilevel

  integer::i,igrid,ncache,iskip,ngrid
  integer::ind,idim,ivar,ix,iy,iz,nx_loc
  integer,dimension(1:nvector),save::ind_grid,ind_cell

  real(dp)::dx,dx_loc,scale
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp),dimension(1:nvector,1:ncrvars),save::uu

  if(.not.cr_advect)return
  if(numbtot(1,ilevel)==0)return
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

  ! Local constants
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  ncache=active(ilevel)%ngrid

  !-------------------------------------------------------
  ! Compute CR initial conditions from subroutine cr_condinit
  !-------------------------------------------------------
  ! Loop over grids by vector sweeps
  do igrid=1,ncache,nvector
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

        ! Call CR initial condition routine
        call cr_condinit(xx,uu,dx_loc,ngrid,ilevel)
        ! Scatter variables
        do ivar=1,ncrvars
           do i=1,ngrid
              cruold(ind_cell(i),ivar)=uu(i,ivar)
           end do
        end do
     end do
     ! End loop over cells
  end do
  ! End loop over grids

111 format('   Entering cr_init_flow_fine for level ',I2)

end subroutine cr_init_flow_fine
!################################################################
!################################################################
!################################################################
!################################################################
SUBROUTINE cr_upload_fine(ilevel)
! Restriction operator (averaging down) for the SEPARATE CR variables only.
! Faithful analog of rt/rt_interpol_hydro.f90's rt_upload_fine (the
! separate-array sibling of the generic upload_fine, which restricts the gas
! array uold alone). Central transformation vs RT: CR state lives in
! cruold(:,1:ncrvars) with iCRu=1, so the variable loop runs over 1:ncrvars
! and the RT photon-tracer (NTRACEGROUPS) branches have no CR analog.
!------------------------------------------------------------------------
  use amr_commons
  use cr_parameters
  use cr_hydro_commons
  implicit none
  integer::ilevel
  integer::i,ncache,igrid,ngrid,ind,iskip,nsplit,icell
  integer,dimension(1:nvector),save::ind_grid,ind_cell,ind_split
  logical,dimension(1:nvector),save::ok
!------------------------------------------------------------------------
  if(.not.cr_advect)return
  if(ilevel==nlevelmax)return
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Loop over active grids by vector sweeps
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

        ! Gather split cells
        do i=1,ngrid
           ok(i)=son(ind_cell(i))>0
        end do

        ! Count split cells
        nsplit=0
        do i=1,ngrid
           if(ok(i))nsplit=nsplit+1
        end do

        ! Upload for selected cells
        if(nsplit>0)then
           icell=0
           do i=1,ngrid
              if(ok(i))then
                 icell=icell+1
                 ind_split(icell)=ind_cell(i)
              end if
           end do
           call cr_upl(ind_split,nsplit)
        end if

     end do
     ! End loop over cells

  end do
  ! End loop over grids

111 format('   Entering cr_upload_fine for level',i2)

END SUBROUTINE cr_upload_fine
!################################################################
!################################################################
!################################################################
!################################################################
SUBROUTINE cr_upl(ind_cell,ncell)
! Restriction operation (averaging down) for the CR variables only:
! average each coarse split cell's 2^ndim children into the parent over
! 1:ncrvars. Faithful analog of rt_upl (rt/rt_interpol_hydro.f90:66-130)
! with cruold in place of rtuold and no photon-tracer normalisation.
!------------------------------------------------------------------------
  use amr_commons
  use cr_parameters
  use cr_hydro_commons
  implicit none
  integer::ncell
  integer,dimension(1:nvector)::ind_cell
  integer ::ivar,i,ind_son,iskip_son
  integer ,dimension(1:nvector),save::igrid_son,ind_cell_son
  real(dp),dimension(1:nvector),save::getx
!------------------------------------------------------------------------
  ! Get child oct index
  do i=1,ncell
     igrid_son(i)=son(ind_cell(i))
  end do

  ! Loop over CR variables
  do ivar=1,ncrvars
     getx(1:ncell)=0.0d0
     do ind_son=1,twotondim
        iskip_son=ncoarse+(ind_son-1)*ngridmax
        do i=1,ncell
           ind_cell_son(i)=iskip_son+igrid_son(i)
        end do
        do i=1,ncell
           getx(i)=getx(i)+cruold(ind_cell_son(i),ivar)
        end do
     end do

     ! Scatter result to cells
     do i=1,ncell
        cruold(ind_cell(i),ivar)=getx(i)/dble(twotondim)
     end do

  end do
  ! End loop over variables

END SUBROUTINE cr_upl
!################################################################
!################################################################
!################################################################
!################################################################
