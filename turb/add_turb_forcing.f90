#if USE_TURB==1
!################################################################
!################################################################
!################################################################
!################################################################
subroutine calc_turb_forcing(ilevel)
  use amr_commons
  use hydro_commons
  use turb_commons
  implicit none
  integer::ilevel
  !-------------------------------------------------------------------
  ! Calculate forcing from turbulent accelerations
  !-------------------------------------------------------------------
  integer::ncache,ngrid,i,igrid,iskip,ind
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  real(kind=dp) :: x_cell(1:ndim,1:nvector)     ! Cell positions
  real(kind=dp) :: rho(1:nvector)               ! Cell densities
  real(kind=dp) :: aturb(1:ndim,1:nvector)      ! Turbulent acceleration
  real(dp),dimension(1:3,1:twotondim)::xc
  real(dp),dimension(1:3)::skip_loc
  real(dp)::dx,dx_loc,scale
  integer::ix,iy,iz,idim
  integer::nx_loc

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Mesh size at level ilevel in coarse cell units
  dx=0.5D0**ilevel
  
  ! Rescaling factors
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
!   scale=boxlen/dble(nx_loc)
  scale=turb_gs_real/dble(nx_loc)
  dx_loc=dx*scale

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(1,ind)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(2,ind)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(3,ind)=(dble(iz)-0.5D0)*dx
  end do

  ! Loop over active grids by vector sweeps
  ncache=active(ilevel)%ngrid
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
        do i=1,ngrid
           do idim=1,ndim
              x_cell(idim,i)=xg(ind_grid(i),idim)+xc(idim,ind)
           end do
        end do
        ! Rescale position from code units to 0->turb_gs_real units
        do i=1,ngrid
           do idim=1,ndim
              x_cell(idim,i)=(x_cell(idim,i)-skip_loc(idim))*scale
           end do
        end do
        
        ! Gather cell densities
        do i=1,ngrid
           rho(i) = uold(ind_cell(i), 1)
        end do
        
        call turb_force_calc(ngrid, x_cell, rho, aturb)
        
        do i=1,ngrid
           do idim=1,ndim
              fturb(ind_cell(i), idim) = aturb(idim, i)
           end do
        end do
     
     end do
     ! End loop over cells

  end do
  ! End loop over grids

111 format('   Entering calc_turb_forcing for level',i2)

end subroutine calc_turb_forcing
!################################################################
!################################################################
!################################################################
!################################################################
subroutine add_turb_forcing(ilevel, dteff)
  use amr_commons
  use hydro_commons
  use turb_commons
  implicit none
  integer::ilevel
  real(dp)::dteff
  !-------------------------------------------------------------------
  ! Update velocity from turbulent accelerations
  !-------------------------------------------------------------------
  integer::ncache,ngrid,i,igrid,iskip,ind,nx_loc
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  
  !if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,112)ilevel
  
  ! Mesh spacing at that level
  nx_loc=icoarse_max-icoarse_min+1
  
  ! Loop over active grids by vector sweeps
  ncache=active(ilevel)%ngrid
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
        
        call turb_forcing1(ind_cell,ngrid,dteff)
     end do
     ! End loop over cells

  end do
  ! End loop over grids


112 format('   Entering add_turb_forcing for level',i2)
  
end subroutine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine turb_forcing1(ind_cell,ncell,dteff)
  use amr_commons
  use hydro_commons
  use turb_commons
  implicit none
  integer::ncell
  real(dp)::dteff
  integer,dimension(1:nvector)::ind_cell
  !-------------------------------------------------------------------
  ! Turbulence update for hydro variables
  !-------------------------------------------------------------------
  integer::i,idim,neul=5,nndim=3
  real(dp),dimension(1:nvector),save::pp

  ! Compute internal + magnetic + radiative energy
  ! TC: correct index?
  do i=1,ncell
     pp(i)=uold(ind_cell(i),neul)
  end do
  do idim=1,nndim
     do i=1,ncell
        pp(i)=pp(i)-0.5*uold(ind_cell(i),idim+1)**2/max(uold(ind_cell(i),1),smallr)
     end do
  end do
  do i=1,ncell
     uold(ind_cell(i),neul)=pp(i)
  end do
  
  ! Update momentum
  do idim=1,ndim
     do i=1,ncell
        pp(i)=uold(ind_cell(i),idim+1)+ &
             & max(uold(ind_cell(i),1),smallr)*fturb(ind_cell(i),idim)*dteff
     end do
     do i=1,ncell
        uold(ind_cell(i),idim+1)=pp(i)
     end do
  end do
  
  ! Update total energy
  do i=1,ncell
     pp(i)=uold(ind_cell(i),neul)
  end do
  do idim=1,nndim
     do i=1,ncell
        pp(i)=pp(i)+0.5*uold(ind_cell(i),idim+1)**2/max(uold(ind_cell(i),1),smallr)
     end do
  end do
  do i=1,ncell
     uold(ind_cell(i),neul)=pp(i)
  end do
  
end subroutine turb_forcing1
!################################################################
!################################################################
!################################################################
!################################################################

#endif
