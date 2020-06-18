subroutine courant_fine(ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::info
  real(kind=8),dimension(3)::comm_buffin,comm_buffout
#endif
  integer::ilevel
  !----------------------------------------------------------------------
  ! Using the Courant-Friedrich-Levy stability condition,               !
  ! this routine computes the maximum allowed time-step.                !
  !----------------------------------------------------------------------
  integer::i,ivar,idim,ind,ncache,igrid,iskip
  integer::nleaf,ngrid,nx_loc
  integer,dimension(1:nvector),save::ind_grid,ind_cell,ind_leaf

  real(dp)::dt_lev,dx,vol,scale
  real(kind=8)::mass_loc,ekin_loc,eint_loc,dt_loc
  real(kind=8)::mass_all,ekin_all,eint_all,dt_all
  real(dp),dimension(1:nvector,1:nvar),save::uu
  real(dp),dimension(1:nvector,1:ndim),save::gg

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  mass_all=0.0d0; mass_loc=0.0d0
  ekin_all=0.0d0; ekin_loc=0.0d0
  eint_all=0.0d0; eint_loc=0.0d0
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

      ! Gather gravitational acceleration
      gg=0.0d0
      if(poisson)then
          do idim=1,ndim
            do i=1,nleaf
                gg(i,idim)=f(ind_leaf(i),idim)
            end do
          end do
      end if

      ! Compute total mass
      do i=1,nleaf
          mass_loc=mass_loc+uu(i,1)*vol
      end do

      ! Compute total energy
      do i=1,nleaf
          ekin_loc=ekin_loc+uu(i,ndim+2)*vol
      end do

      ! Compute total internal energy
      do i=1,nleaf
          eint_loc=eint_loc+uu(i,ndim+2)*vol
      end do
      do ivar=1,ndim
          do i=1,nleaf
            eint_loc=eint_loc-0.5d0*uu(i,1+ivar)**2/max(uu(i,1),smallr)*vol
          end do
      end do
#if NENER>0
      do ivar=1,nener
          do i=1,nleaf
            eint_loc=eint_loc-uu(i,ndim+2+ivar)*vol
          end do
      end do
#endif

      ! Compute CFL time-step
      if(nleaf>0)then
          call cmpdt(uu,gg,dx,dt_lev,nleaf)
          dt_loc=min(dt_loc,dt_lev)
      end if

    end do
    ! End loop over cells

  end do
  ! End loop over grids

  ! Compute global quantities
#ifndef WITHOUTMPI
  comm_buffin(1)=mass_loc
  comm_buffin(2)=ekin_loc
  comm_buffin(3)=eint_loc
  call MPI_ALLREDUCE(comm_buffin,comm_buffout,3,MPI_DOUBLE_PRECISION,MPI_SUM,&
       &MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(dt_loc,dt_all,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
       &MPI_COMM_WORLD,info)
  mass_all=comm_buffout(1)
  ekin_all=comm_buffout(2)
  eint_all=comm_buffout(3)
#endif
#ifdef WITHOUTMPI
  mass_all=mass_loc
  ekin_all=ekin_loc
  eint_all=eint_loc
  dt_all=dt_loc
#endif
  mass_tot=mass_tot+mass_all
  ekin_tot=ekin_tot+ekin_all
  eint_tot=eint_tot+eint_all
  dtnew(ilevel)=MIN(dtnew(ilevel),dt_all)

111 format('   Entering courant_fine for level ',I2)

end subroutine courant_fine





!! Begin subroutine Energy_fine !!

subroutine energy_fine(ilevel)
  use amr_commons
  use hydro_commons

  implicit none

  integer :: ilevel
  !----------------------------------------------------------
  ! This routine add heating at the bottom of the 
  ! convection zone
  !----------------------------------------------------------
  integer::igrid,ngrid,ncache,i,ind,iskip,ix,iy,iz
  integer::nx_loc,idim,ivar
#ifdef SOLVERmhd
  integer::neul=5
#else
  integer::neul=ndim+2
#endif
  real(dp)::dx,dx_loc,scale,d,u,v,w,A,B,C
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:3)::skip_loc

  integer ,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp),dimension(1:nvector)::ee
  real(dp),dimension(1:nvector,1:nvar)::uut ! JRCC
  real(dp),dimension(1:nvector)::req,peq
  real(kind=8)::uud,ekin,uuv

  if(numbtot(1,ilevel)==0)return

  ! Mesh size at level ilevel in coarse cell units
  dx=0.5D0**ilevel

  ! Rescaling factors
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=dble(nx_loc)/boxlen
  dx_loc=dx/scale

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  !-------------------------------------
  ! Add analytical heating/cooling
  !-------------------------------------
  ncache=active(ilevel)%ngrid

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
              xx(i,idim)=(xx(i,idim)-skip_loc(idim))/scale
           end do
        end do
        

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! TODO
        ! 1. Compute kinetic energy
        ! 2. Remove kinetic energy from total energy (uold(neul))
        ! 3. Call eneana and save heating/cooling and damping (in vv)
        ! 4. Add heating cooling to total energy
        ! 5. Compute new momentum with damping values vv
        ! 6. Add back kinetic energy with new momentum variable

        ! Remove kinetic energy
        do i=1,ngrid
          ekin = 0.0d0
          uud = max(uold(ind_cell(i),1),smallr)
          do idim=1,ndim
             uuv = uold(ind_cell(i),idim+1)/uud
             ekin = ekin+0.5d0*uud*uuv**2
          end do
          uold(ind_cell(i),neul) = uold(ind_cell(i),neul)-ekin
        end do
        
        ! JRCC
        ! Scatter variables
        ! See init_flow_fine.f90, same procedure for condinit.f90 ?
        do ivar=1,nvar
          do i=1,ngrid
             uut(i,ivar)=uold(ind_cell(i),ivar)
          end do
        end do

        do i=1,ngrid
          req(i) = rho_eq(ind_cell(i))
          peq(i) = p_eq(ind_cell(i))
        end do 

        ! Spongy layers that damp the profiles towards the equilibrium ones
        ! In: rho, rho*u, rho*v, (rho*w), e_int - e_kin
        ! Out: new variables "closer" to equilibrium 
        call spongelayers(xx,uut,req,peq,t,ngrid)
        
        ! rescatter variables
        do ivar=1,nvar
          do i=1,ngrid
            ! if (abs((uold(ind_cell(i),ivar)-uut(i,ivar))/uold(ind_cell(i),ivar))>10.0**(-8.0)) then
            !     write(*,*)'HERE var=',ivar,', x=',xx(i,1), ', y=',xx(i,2)
            !     write(*,*)'  uold=',uold(ind_cell(i),ivar),', uu=',uut(i,ivar)
            ! end if 
            uold(ind_cell(i),ivar) = uut(i,ivar)
          end do
        end do 

        ! Impose analytical energy field and damping
        call eneana(xx,ee,dx_loc,t,ngrid)
        ! Update total energy
        do i=1,ngrid
          uold(ind_cell(i),neul)=uold(ind_cell(i),neul)+ee(i)*dtnew(ilevel)
#if NVAR>NDIM+2
          uold(ind_cell(i),ndim+3)=uold(ind_cell(i),ndim+3)+(gamma-1.0)/uut(i,1)**(gamma-1.0)*ee(i)*dtnew(ilevel)
#endif
        end do

        ! Add back kinetic energy
        do i=1,ngrid
          ekin = 0.0d0
          uud = max(uold(ind_cell(i),1),smallr)
          do idim=1,ndim
             uuv = uold(ind_cell(i),idim+1)/uud
             ekin = ekin+0.5d0*uud*uuv**2
          end do
          uold(ind_cell(i),neul) = uold(ind_cell(i),neul)+ekin
        end do
        
    

     end do
     ! End loop over cells

  end do
  ! End loop over grids

end subroutine energy_fine 


subroutine check_cons(ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::info
  real(kind=8),dimension(3)::comm_buffin,comm_buffout
#endif
  integer::ilevel
  !----------------------------------------------------------------------
  ! Check mass and energy conservation
  !----------------------------------------------------------------------
  integer::i,ivar,ind,ncache,igrid,iskip
  integer::nleaf,ngrid
  integer,dimension(1:nvector),save::ind_grid,ind_cell,ind_leaf

  real(dp)::dx,vol
  real(kind=8)::mass_loc,ekin_loc,eint_loc
  real(kind=8)::mass_all,ekin_all,eint_all
  real(dp),dimension(1:nvector,1:nvar),save::uu

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  mass_all=0.0d0; mass_loc=0.0d0
  ekin_all=0.0d0; ekin_loc=0.0d0
  eint_all=0.0d0; eint_loc=0.0d0

  ! Mesh spacing at that level
  dx=0.5D0**ilevel*boxlen
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

      ! Compute total mass
      do i=1,nleaf
          mass_loc=mass_loc+uu(i,1)*vol
      end do

      ! Compute total energy
      do i=1,nleaf
          ekin_loc=ekin_loc+uu(i,ndim+2)*vol
      end do

      ! Compute total internal energy
      do i=1,nleaf
          eint_loc=eint_loc+uu(i,ndim+2)*vol
      end do
      do ivar=1,ndim
          do i=1,nleaf
            eint_loc=eint_loc-0.5d0*uu(i,1+ivar)**2/max(uu(i,1),smallr)*vol
          end do
      end do
#if NENER>0
      do ivar=1,nener
          do i=1,nleaf
            eint_loc=eint_loc-uu(i,ndim+2+ivar)*vol
          end do
      end do
#endif
    end do
    ! End loop over cells
  end do
  ! End loop over grids

  ! Compute global quantities
#ifndef WITHOUTMPI
  comm_buffin(1)=mass_loc
  comm_buffin(2)=ekin_loc
  comm_buffin(3)=eint_loc
  call MPI_ALLREDUCE(comm_buffin,comm_buffout,3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  mass_all=comm_buffout(1)
  ekin_all=comm_buffout(2)
  eint_all=comm_buffout(3)
#endif
#ifdef WITHOUTMPI
  mass_all=mass_loc
  ekin_all=ekin_loc
  eint_all=eint_loc
#endif
  mass_tot=mass_tot+mass_all
  ekin_tot=ekin_tot+ekin_all
  eint_tot=eint_tot+eint_all

111 format('   Entering check_cons for level ',I2)

end subroutine check_cons
