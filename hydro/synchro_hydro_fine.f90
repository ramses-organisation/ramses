!################################################################
!################################################################
!################################################################
!################################################################
subroutine synchro_hydro_fine(ilevel,dteff,which_force)
  use amr_commons
  use hydro_commons
#if USE_TURB==1
  use turb_commons
#endif
  implicit none
  integer::ilevel
  real(dp)::dteff
  integer::which_force !gravity=1, turbulence=2
  !-------------------------------------------------------------------
  ! Update velocity  from gravitational acceleration
  !-------------------------------------------------------------------
  integer::ncache,ngrid,i,igrid,iskip,ind
  integer,dimension(1:nvector),save::ind_grid,ind_cell

#if USE_TURB==1
  if(.not. (poisson.or.turb))return
#else
  if(.not. poisson)return
  if(numbtot(1,ilevel)==0)return
#endif
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
           ind_cell(i)=ind_grid(i)+iskip
        end do
        call synchydrofine1(ind_cell,ngrid,dteff,which_force)
     end do
     ! End loop over cells

  end do
  ! End loop over grids

111 format('   Entering synchro_hydro_fine for level',i2)

end subroutine synchro_hydro_fine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine synchydrofine1(ind_cell,ncell,dteff,which_force)
  use amr_commons
  use hydro_commons
  use poisson_commons
#if USE_TURB==1
  use turb_commons
#endif
  implicit none
  integer::ncell
  integer::which_force !gravity=1, turbulence=2
  real(dp)::dteff
  integer,dimension(1:nvector)::ind_cell
  !-------------------------------------------------------------------
  ! Force (gravity or turbulence) update for hydro variables
  !-------------------------------------------------------------------
  integer::i,idim
#ifdef SOLVERmhd
  integer::neul=5,nndim=3
#else
  integer::neul=ndim+2,nndim=ndim
#endif
  real(dp),dimension(1:nvector),save::pp

  ! Compute internal + magnetic + radiative energy
  do i=1,ncell
     pp(i)=uold(ind_cell(i),neul)
  end do
  do idim=1,nndim
     do i=1,ncell
        pp(i)=pp(i)-0.5d0*uold(ind_cell(i),idim+1)**2/max(uold(ind_cell(i),1),smallr)
     end do
  end do
  do i=1,ncell
     uold(ind_cell(i),neul)=pp(i)
  end do

  ! Update momentum
  do idim=1,ndim
#if USE_TURB==1
     if (which_force==2) then
        !turbulence
        do i=1,ncell
           pp(i)=uold(ind_cell(i),idim+1)+ &
             & max(uold(ind_cell(i),1),smallr)*fturb(ind_cell(i),idim)*dteff
        end do
     else
#endif
        !gravity
        do i=1,ncell
           pp(i)=uold(ind_cell(i),idim+1)+max(uold(ind_cell(i),1),smallr)*f(ind_cell(i),idim)*dteff
        end do
#if SOLVERhydro
        if(strict_equilibrium>0)then
           do i=1,ncell
              pp(i)=pp(i)-rho_eq(ind_cell(i))*f(ind_cell(i),idim)*dteff
           end do
        endif
#endif
#if USE_TURB==1
     endif
#endif
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
        pp(i)=pp(i)+0.5d0*uold(ind_cell(i),idim+1)**2/max(uold(ind_cell(i),1),smallr)
     end do
  end do
  do i=1,ncell
     uold(ind_cell(i),neul)=pp(i)
  end do

end subroutine synchydrofine1
!################################################################
!################################################################
!################################################################
!################################################################
