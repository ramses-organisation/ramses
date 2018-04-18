subroutine courant_fine(ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  use mpi_mod
  implicit none
  integer::ilevel
  !----------------------------------------------------------------------
  ! Using the Courant-Friedrich-Levy stability condition,               !
  ! this routine computes the maximum allowed time-step.                !
  !----------------------------------------------------------------------
  integer::i,ivar,idim,ind,ncache,igrid,iskip
  integer::info,nleaf,ngrid,nx_loc
  integer,dimension(1:nvector),save::ind_grid,ind_cell,ind_leaf

  real(dp)::dt_lev,dx,vol,scale
  real(kind=8)::mass_loc,e_loc,momx_loc,momy_loc,momz_loc,dt_loc,gam_loc=1d0,mom_all
  real(kind=8)::mass_all,dt_all,gam_all=1d0,e_all,momx_all,momy_all,momz_all
  real(kind=8),dimension(3)::comm_buffin,comm_buffout
  real(dp),dimension(1:nvector,1:nvar),save::uu
  real(dp),dimension(1:nvector,1:ndim),save::gg
  real(dp),dimension(1:nvector,1:ndim),save::q
  real(dp)::lor,v
!  real(dp)::lor_max_loc,lor_max_all

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  gam_loc=0d0
  mass_all=0.0d0; mass_loc=0.0d0
  e_all=0d0; e_loc=0d0
  mom_all=0d0;momx_all=0d0;momy_all=0d0;momz_all=0d0
  momx_loc=0d0 ;momy_loc=0d0; momz_loc=0d0
  dt_all=dtnew(ilevel); dt_loc=dt_all

  ! Mesh spacing at that level
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5D0**ilevel*scale
  vol=dx**ndim

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

        ! Gather leaf cells
        nleaf=0
        do i=1,ngrid
           if(son(ind_cell(i))==0)then
              nleaf=nleaf+1
              ind_leaf(nleaf)=ind_cell(i)
           end if
        end do


        ! Gather hydro variables
        do ivar=1,nvar
           do i=1,nleaf
              uu(i,ivar)=uold(ind_leaf(i),ivar)
           end do
        end do

        call ctoprimbis(uu,nleaf,q)
        ! Compute total mass,energy,momentum
        do i=1,nleaf
           mass_loc=mass_loc+uu(i,1)*vol
           v=sqrt(q(i,2)**2+q(i,3)**2+q(i,4)**2)
           lor=(1d0-v**2)**(-1./2.)
           gam_loc=max(gam_loc,lor)
           momx_loc=momx_loc+uu(i,2)*vol
           momy_loc=momy_loc+uu(i,3)*vol
           momz_loc=momz_loc+uu(i,4)*vol
           e_loc=e_loc+uu(i,5)*vol
       end do
        ! Compute CFL time-step
        if(nleaf>0)then
           call cmpdt(uu,gg,dx,dt_lev,nleaf)
           dt_loc=min(dt_loc,dt_lev)
        end if

     end do
     ! End loop over cells

  end do
  ! End loop over grids
! gam_loc saves the max lorentz factor of each processor

  ! Compute global quantities
#ifndef WITHOUTMPI
  comm_buffin(1)=mass_loc
  comm_buffin(2)=momx_loc
  comm_buffin(3)=momy_loc
  comm_buffin(4)=momz_loc
  call MPI_ALLREDUCE(comm_buffin,comm_buffout,4,MPI_DOUBLE_PRECISION,MPI_SUM,&
       &MPI_COMM_WORLD,info)
  mass_all=comm_buffout(1)
  momx_all=comm_buffout(2)
  momy_all=comm_buffout(3)
  momz_all=comm_buffout(4)


 comm_buffin(1)=e_loc
  call MPI_ALLREDUCE(comm_buffin,comm_buffout,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
       &MPI_COMM_WORLD,info)
  e_all=comm_buffout(1)


  mom_all=sqrt(momx_all**2+momy_all**2+momz_all**2)

  call MPI_ALLREDUCE(dt_loc,dt_all,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
       &MPI_COMM_WORLD,info)

  comm_buffin(1)=gam_loc
  call MPI_ALLREDUCE(comm_buffin,comm_buffout,1,MPI_DOUBLE_PRECISION,MPI_MAX,&
       &MPI_COMM_WORLD,info)
  gam_all=comm_buffout(1)

#endif

#ifdef WITHOUTMPI
  mass_all=mass_loc
  gam_all=gam_loc
  dt_all=dt_loc
  e_all=e_loc
  mom_all=sqrt(momx_loc**2+momy_loc**2+momz_loc**2)
#endif
  lor_max=gam_all
  mass_tot=mass_tot+mass_all


  e_tot=e_tot+e_all
  mom_tot= mom_tot+mom_all
  dtnew(ilevel)=MIN(dtnew(ilevel),dt_all)


111 format('   Entering courant_fine for level ',I2)

end subroutine courant_fine
