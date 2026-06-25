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

!************************************************************************
SUBROUTINE cr_region_condinit(x,u,dx,nn,ilevel)
  ! Fill the separated CR buffer u(1:nn,1:ncrvars) from the namelist
  ! crmom_region per CR-OWNED region geometry (cr_nregion / cr_region_type /
  ! cr_reg_*_center / cr_reg_length_* / cr_exp_region in &cr_params), the
  ! analogue of rt_region_condinit (rt/rt_init_flow_fine.f90) for the RT
  ! module. Faithful port of the CR-handling part of ramses_cral
  ! mhd/init_flow_fine.f90 region_condinit (the square/point region loop),
  ! keeping ONLY the CR slots: cral wrote q(i,icru+ivar-1) with icru=nvar+4;
  ! here we write u(i,ivar) with the separated layout (iCRu=1), ivar=1..ncrvars
  ! maps energy + ndim fluxes per group identically. The region-membership
  ! arithmetic is verbatim from region_condinit; only the geometry source is
  ! the CR-owned cr_reg_* set (not the gas nregion/x_center/length_x/...). The
  ! cr_reg_* in every test namelist equal the gas &init_params geometry, so
  ! this is bit-identical to the previous gas-geometry version. Lives in core
  ! (mirroring rt_region_condinit) and is called from the patch cr_condinit;
  ! ilevel/dx/vol are unused in the body (the square branch ignores them, the
  ! point branch is a CR no-op) but ilevel is in the signature to mirror
  ! rt_region_condinit. Default (no region covers a cell) leaves the smallcr
  ! floor / zero flux that cr_condinit set before calling this routine; cral's
  ! default is 0, but every 422 cell is covered by a region so the difference
  ! never shows.
  use amr_parameters
  use cr_parameters, only: cr_nregion,cr_region_type,cr_reg_x_center &
       & ,cr_reg_y_center,cr_reg_z_center,cr_reg_length_x,cr_reg_length_y &
       & ,cr_reg_length_z,cr_exp_region,cr_reg_group,crmom_region,ncrvars &
       & ,iCRu,smallcr
  implicit none
  integer ::nn
  integer::ilevel
  real(dp)::dx
  real(dp),dimension(1:nvector,1:ncrvars)::u
  real(dp),dimension(1:nvector,1:ndim  )::x
  integer::i,k,ivar
  real(dp)::vol,r,xn,yn,zn,en

  ! Loop over CR initial-conditions regions
  do k=1,cr_nregion

     ! For "square" regions only:
     if(cr_region_type(k) .eq. 'square')then
        en=cr_exp_region(k)
        do i=1,nn
           xn=0.0d0; yn=0.0d0; zn=0.0d0
           xn=2.0d0*abs(x(i,1)-cr_reg_x_center(k))/cr_reg_length_x(k)
#if NDIM>1
           yn=2.0d0*abs(x(i,2)-cr_reg_y_center(k))/cr_reg_length_y(k)
#endif
#if NDIM>2
           zn=2.0d0*abs(x(i,3)-cr_reg_z_center(k))/cr_reg_length_z(k)
#endif
           if(cr_exp_region(k)<10)then
              r=(xn**en+yn**en+zn**en)**(1.0/en)
           else
              r=max(xn,yn,zn)
           end if
           ! If cell lies within region, REPLACE CR variables by region values
           if(r<1.0)then
              do ivar=1,ncrvars
                 u(i,ivar)=crmom_region(k,ivar)
              end do
           end if
        end do
     end if

     ! For "point" regions only: cral's region_condinit point branch updates
     ! ONLY the gas primitives (density/velocity/pressure/NENER) and never
     ! touches the CR slots (mhd/init_flow_fine.f90:612-645). Faithful port =>
     ! the CR buffer is left untouched here, keeping whatever the enclosing
     ! 'square' region (or the smallcr floor) set. 423's point region has
     ! crmom_region(2,1)=0 so this was already a numerical no-op, but matching
     ! cral exactly avoids any CR injection if a future test uses a point region
     ! with nonzero crmom_region.
  end do

END SUBROUTINE cr_region_condinit
!################################################################
!################################################################
!################################################################
!################################################################
