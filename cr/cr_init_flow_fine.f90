!################################################################
!################################################################
!################################################################
!################################################################
subroutine cr_init_flow
  ! Initialise the separated cosmic-ray state arrays cruold(:,1:ncrvar)
  ! on refinement. Inert unless cr_advect.
  use amr_commons
  use cr_parameters
  use cr_hydro_commons, ONLY: cruold
  implicit none

  integer::ilevel,ivar

  if(.not.cr_advect)return
  if(verbose)write(*,*)'Entering cr_init_flow'
  do ilevel=nlevelmax,1,-1
     if(ilevel>=levelmin)call cr_init_flow_fine(ilevel)
     ! Restrict the separate CR array cruold(:,1:ncrvar) fine->coarse onto
     ! coarse split cells (generic upload_fine restricts only the gas uold).
     call cr_upload_fine(ilevel)
     do ivar=1,ncrvar
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
  ! cr_condinit, and scatter the result into cruold(:,1:ncrvar).
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
  real(dp),dimension(1:nvector,1:ncrvar),save::uu

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
        do ivar=1,ncrvar
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
  ! Fill the separated CR buffer u(1:nn,1:ncrvar) from namelist crmom_region
  ! per CR-owned region geometry (cr_nregion / cr_region_type / cr_reg_* in &cr_params).
  use amr_parameters
  use cr_parameters, only: cr_nregion,cr_region_type,cr_reg_x_center &
       & ,cr_reg_y_center,cr_reg_z_center,cr_reg_length_x,cr_reg_length_y &
       & ,cr_reg_length_z,cr_exp_region,cr_reg_group,crmom_region,ncrvar &
       & ,iCRu,smallcr
  implicit none
  integer ::nn
  integer::ilevel
  real(dp)::dx
  real(dp),dimension(1:nvector,1:ncrvar)::u
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
              do ivar=1,ncrvar
                 u(i,ivar)=crmom_region(k,ivar)
              end do
           end if
        end do
     end if

     ! Point regions are a CR no-op: they update only the gas primitives and
     ! never touch the CR slots, so the CR buffer is left untouched here.
  end do

END SUBROUTINE cr_region_condinit
!################################################################
!################################################################
!################################################################
!################################################################
