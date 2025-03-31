subroutine courant_fine(ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::info
#endif
  integer::ilevel
  !----------------------------------------------------------------------
  ! Using the Courant-Friedrich-Levy stability condition,               !
  ! this routine computes the maximum allowed time-step.                !
  !----------------------------------------------------------------------
  integer::i,ivar,idim,imat,ind,ncache,igrid,iskip
  integer::nleaf,ngrid,nx_loc
  integer,dimension(1:nvector),save::ind_grid,ind_cell,ind_leaf

  real(dp)::dt_lev,dx,vol,skip_loc,scale,dx_loc,fourpi,twopi
  real(kind=8)::mass_loc,ekin_loc,eint_loc,dt_loc
  real(kind=8)::mass_all,ekin_all,eint_all,dt_all
  real(dp),dimension(1:nvector,1:nvar),save::uu
  real(dp),dimension(1:nvector,1:ndim),save::gg
  real(dp),dimension(1:nvector),save::dtot

  real(dp),dimension(1:nvector),save::rloc
  real(dp),dimension(1:8)::xc
  integer ::ix,iy,iz

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  mass_all=0.0d0; mass_loc=0.0d0
  ekin_all=0.0d0; ekin_loc=0.0d0
  eint_all=0.0d0; eint_loc=0.0d0
  dt_all=dtnew(ilevel); dt_loc=dt_all

  ! Mesh spacing at that level
  twopi=2.0d0*ACOS(-1.0d0)
  fourpi=4.0d0*ACOS(-1.0d0)
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol=dx_loc**ndim
  skip_loc=dble(icoarse_min)

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     xc(ind)=(dble(ix)-0.5D0)*dx
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
              rloc(nleaf)=(xg(ind_grid(i),1)+xc(ind)-skip_loc)*scale
           end if
        end do

        ! Gather hydro variables
        do ivar=1,nvar
           do i=1,nleaf
              uu(i,ivar)=uold(ind_leaf(i),ivar)
           end do
        end do

        ! Calculate total density
        dtot(1:nleaf) = 0.0
        do imat=1,nmat
           do i=1,nleaf
              dtot(i) = dtot(i) + uu(i,nmat+imat)
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
           mass_loc = mass_loc + dtot(i)*vol
        end do

        ! Compute total energy
        do imat=1,nmat
           do i=1,nleaf
              ekin_loc = ekin_loc + uu(i,2*nmat+ndim+imat)*vol
           end do
        end do

        ! Compute CFL time-step
        if(nleaf>0)then
           call cmpdt(uu,gg,rloc,dx_loc,dt_lev,nleaf)
           dt_loc=min(dt_loc,dt_lev)
        end if

     end do
     ! End loop over cells

  end do
  ! End loop over grids

  ! Compute global quantities
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(mass_loc,mass_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
       & MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(ekin_loc,ekin_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
       & MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(eint_loc,eint_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
       & MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(  dt_loc,  dt_all,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
       & MPI_COMM_WORLD,info)
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
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmpdt(uu,grav,rr,dx,dt,ncell)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none
  integer::ncell
  logical::inv
  real(dp)::dx,dt
  real(dp),dimension(1:nvector,1:nvar)::uu
  real(dp),dimension(1:nvector,1:ndim)::grav
  real(dp),dimension(1:nvector)::rr,dtot

  real(dp),dimension(1:nvector,1:npri),save::qq
  real(dp),dimension(1:nvector,1:nmat),save::ff,gg
  real(dp),dimension(1:nvector),save::gg_mat,ee_mat,pp_mat,cc_mat
  real(dp),dimension(1:nvector),save::ekin,cc,st
  real(dp)::dtcell,eps
  integer::k,idim,imat

  ! Convert to primitive variable

  ! Volume fraction and fluid density
  do imat = 1,nmat
     do k = 1,ncell
        ff(k,imat) = uu(k,imat)
        gg(k,imat) = uu(k,nmat+imat)/ff(k,imat)
     end do
  end do

  ! Compute total density
  dtot(1:ncell) = 0.0
  do imat=1,nmat
    do k = 1,ncell
       dtot(k)  = dtot(k) + uu(k,nmat+imat)
    end do
  end do

  ! Compute velocity and specific kinetic energy
  ekin(1:ncell)    = 0.0
  do idim = 1,ndim
     do k = 1,ncell
        qq(k,idim) = uu(k,2*nmat+idim)/dtot(k)
        ekin(k)    = ekin(k) + half*qq(k,idim)**2
     end do
  end do

  ! Compute partial internal energies
  do imat=1,nmat
    do k = 1,ncell
       qq(k,ndim+nmat+imat) = uu(k,2*nmat+ndim+imat)/ff(k,imat) - gg(k,imat)*ekin(k)
    end do
  end do

  ! Calculate the total speed of sound
  cc(1:ncell)=0
  inv=.false.
  do imat=1,nmat
    do k=1,ncell
      gg_mat(k) = gg(k,imat)
      ee_mat(k) = qq(k,ndim+nmat+imat)
    end do
    ! Call eos routine
    call eos(gg_mat,ee_mat,pp_mat,cc_mat,imat,inv,ncell)
    do k=1,ncell
      cc(k) = cc(k) + ff(k,imat)*gg(k,imat) * cc_mat(k)**2
    end do
  end do
  ! Convert rho c^2 to c
  cc(1:ncell)=sqrt(cc(1:ncell)/dtot(1:ncell))

  ! Compute wave speed
  do k = 1,ncell
     cc(k) = abs(qq(k,1))+cc(k)
  end do
  do idim = 2,ndim
     do k = 1,ncell
        cc(k) = cc(k) + abs(qq(k,idim))+cc(k)
     end do
  end do

  ! Compute gravity strength ratio
  do k = 1,ncell
     st(k) = zero
  end do
  do idim = 1,ndim
     do k = 1,ncell
        st(k) = st(k) + abs(grav(k,idim))
     end do
  end do
  do k = 1,ncell
     st(k) = st(k)*dx/cc(k)**2
     st(k) = MAX(st(k),0.0001_dp)
  end do

  ! Compute maximum time step for each authorized cell
  dt = courant_factor*dx/smallc
  do k = 1,ncell
     dtcell = dx/cc(k)*(sqrt(one+two*courant_factor*st(k))-one)/st(k)
     dt = min(dt,dtcell)
  end do

end subroutine cmpdt
