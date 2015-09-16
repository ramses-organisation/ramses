!################################################################
!################################################################
!################################################################
!################################################################
subroutine synchro_hydro_fine(ilevel,dteff)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ilevel
  real(dp)::dteff
  !-------------------------------------------------------------------
  ! Update velocity  from gravitational acceleration
  !-------------------------------------------------------------------
  integer::ncache,ngrid,i,igrid,iskip,ind
  ! Variables because of not calling synchydrofine1
  integer::i,idim,neul=ndim+2,nndim=ndim
#ifdef MPROFILING
  include 'mpif.h'
  integer::info
  real(kind=8)::tt1,tt2
  tt1=MPI_WTIME(info)
#endif

  if(.not. poisson)return
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  !$acc data present(uold,f,active,pp_shfv,ind_cell_shfv,ind_grid_shfv)

  ! Loop over active grids by vector sweeps
  ncache=active(ilevel)%ngrid
  !$acc parallel loop gang worker
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     ! Loop over cells
     do ind=1,twotondim
!**************************************************************************
!       call synchydrofine1(ind_cell,ngrid,dteff) - Explicitly in the code
!**************************************************************************
     !$acc loop vector
     do i=1,ngrid
        do idim=1,nndim
        ind_grid_shfv(i)=active(ilevel)%igrid(igrid+i-1)
        iskip=ncoarse+(ind-1)*ngridmax
        ind_cell_shfv(i)=ind_grid_shfv(i)+iskip
     ! Compute internal + magnetic + radiative energy
        pp_shfv(i)=uold(ind_cell_shfv(i),neul)
        pp_shfv(i)=pp_shfv(i)-0.5*uold(ind_cell_shfv(i),idim+1)**2/uold(ind_cell_shfv(i),1)
        uold(ind_cell_shfv(i),neul)=pp_shfv(i)

     ! Update momentum
        pp_shfv(i)=uold(ind_cell_shfv(i),idim+1)+ &
             & uold(ind_cell_shfv(i),1)*f(ind_cell_shfv(i),idim)*dteff
        uold(ind_cell_shfv(i),idim+1)=pp_shfv(i)

     ! Update total energy
       pp_shfv(i)=uold(ind_cell_shfv(i),neul)
       pp_shfv(i)=pp_shfv(i)+0.5*uold(ind_cell_shfv(i),idim+1)**2/uold(ind_cell_shfv(i),1)

       uold(ind_cell_shfv(i),neul)=pp_shfv(i)
       end do
    end do
    end do
    ! End loop over cells
  end do
  ! End loop over grids

  !$acc end data
  
#ifdef MPROFILING
  tt2=MPI_WTIME(info)
  acc_t_synchro_hydro_fine = acc_t_synchro_hydro_fine + (tt2-tt1)
#endif  

111 format('   Entering synchro_hydro_fine for level',i2)

end subroutine synchro_hydro_fine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine synchydrofine1(ind_cell,ncell,dteff)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ncell
  real(dp)::dteff
  integer,dimension(1:nvector)::ind_cell
  !-------------------------------------------------------------------
  ! Gravity update for hydro variables
  !-------------------------------------------------------------------
  integer::i,idim,neul=ndim+2,nndim=ndim
  real(dp),dimension(1:nvector),save::pp

  ! Compute internal + magnetic + radiative energy
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
             & max(uold(ind_cell(i),1),smallr)*f(ind_cell(i),idim)*dteff
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

end subroutine synchydrofine1
!################################################################
!################################################################
!################################################################
!################################################################

