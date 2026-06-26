subroutine courant_fine(ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  use mpi_mod
#if USE_TURB==1
  use turb_commons
#endif
#ifdef NIMHD
  use nimhd_parameters
#endif
  implicit none
#ifndef WITHOUTMPI
  integer::info
  real(kind=8),dimension(4)::comm_buffin,comm_buffout
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
  real(kind=8)::mass_loc,ekin_loc,eint_loc,emag_loc,dt_loc
  real(kind=8)::mass_all,ekin_all,eint_all,emag_all,dt_all
  real(dp),dimension(1:nvector,1:nvar_all),save::uu
  real(dp),dimension(1:nvector,1:ndim),save::gg
#ifdef NIMHD
  ! maybe to see if ambipolar dt at a certain level restricts the global dt?
  ! can probably remove
  real(dp)::dtwad_loc,dtwad_all
  real(dp)::dtambdiff_loc,dtambdiff_lev,dtambdiff_all
  real(dp)::dtmagdiff_loc,dtmagdiff_lev,dtmagdiff_all
  real(dp)::tmag1,tmag2
#endif

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  mass_all=0.0d0; mass_loc=0.0d0
  ekin_all=0.0d0; ekin_loc=0.0d0
  emag_all=0.0d0; emag_loc=0.0d0
  eint_all=0.0d0; eint_loc=0.0d0
  dt_all=dtnew(ilevel); dt_loc=dt_all
#ifdef NIMHD
  ! maybe to see if ambipolar dt at a certain level restricts the global dt?
  ! can probably remove
  dtambdiff_all=dtambdiff(ilevel); dtambdiff_loc=dtambdiff_all
  dtmagdiff_all=dtmagdiff(ilevel); dtmagdiff_loc=dtmagdiff_all
  dtwad_all=dtwad(ilevel); dtwad_loc=dtwad_all
#endif

  ! Mesh spacing at that level
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5D0**ilevel*scale
  vol=dx**ndim

  if (ischeme .eq. 1) then
     CALL velocity_fine(ilevel)
     do ivar=1,nvar_all
        call make_virtual_fine_dp(uold(1,ivar),ilevel)
     end do
     if(simple_boundary)call make_boundary_hydro(ilevel)
  endif

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
        do ivar=1,nvar_all
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

#if USE_TURB==1
        if (turb .AND. turb_type/=3) then
           do idim=1,ndim
              do i=1,nleaf
                 gg(i,idim)=gg(i,idim)+fturb(ind_leaf(i),idim)
              end do
           end do
        end if
#endif

        ! Compute total mass
        do i=1,nleaf
           mass_loc=mass_loc+uu(i,1)*vol
        end do

        ! Compute total energy
        do i=1,nleaf
           ekin_loc=ekin_loc+uu(i,neul)*vol
        end do

        ! Compute total magnetic energy
        do ivar=1,3
           do i=1,nleaf
              emag_loc=emag_loc+0.125d0*(uu(i,5+ivar)+uu(i,nvar+ivar))**2*vol
           end do
        end do

        ! Compute total internal energy
        do i=1,nleaf
           eint_loc=eint_loc+uu(i,neul)*vol
        end do
        do ivar=1,3
           do i=1,nleaf
              eint_loc=eint_loc-0.5d0*uu(i,1+ivar)**2/uu(i,1)*vol &
                   & -0.125d0*(uu(i,neul+ivar)+uu(i,nvar+ivar))**2*vol
           end do
        end do
#if NENER>0
        do ivar=1,nener
           do i=1,nleaf
              eint_loc=eint_loc-uu(i,nhydro+ivar)*vol
           end do
        end do
#endif

        ! Compute CFL time-step
        if(nleaf>0)then
#ifdef NIMHD
           ! Warning: cmpdt_nimhd needs to be done before cmpdt because there uu is altered
           call cmpdt_nimhd(uu,dx,nleaf,dtambdiff_lev,dtmagdiff_lev)
           dtambdiff_loc=min(dtambdiff_loc,dtambdiff_lev)
           dtmagdiff_loc=min(dtmagdiff_loc,dtmagdiff_lev)
#endif
           call cmpdt(uu,gg,dx,dt_lev,nleaf)
           dt_loc=min(dt_loc,dt_lev)
#ifdef NIMHD
           ! store the ideal MHD timestep (without effect of gravity, etc.)
           dtwad_loc=min(dtwad_loc,dt_lev)
#endif
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
  comm_buffin(4)=emag_loc
  call MPI_ALLREDUCE(comm_buffin,comm_buffout,4,MPI_DOUBLE_PRECISION,MPI_SUM,&
       &MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(dt_loc     ,dt_all      ,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
       &MPI_COMM_WORLD,info)
  mass_all=comm_buffout(1)
  ekin_all=comm_buffout(2)
  eint_all=comm_buffout(3)
  emag_all=comm_buffout(4)
#endif
#ifdef WITHOUTMPI
  mass_all=mass_loc
  ekin_all=ekin_loc
  eint_all=eint_loc
  emag_all=emag_loc
  dt_all=dt_loc
#endif

  mass_tot=mass_tot+mass_all
  ekin_tot=ekin_tot+ekin_all
  eint_tot=eint_tot+eint_all
  emag_tot=emag_tot+emag_all
  dtnew(ilevel)=MIN(dtnew(ilevel),dt_all)

#ifdef NIMHD

  ! Compute global quantities
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(dtambdiff_loc,dtambdiff_all,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
       &MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(dtmagdiff_loc,dtmagdiff_all,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
       &MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(dtwad_loc    ,dtwad_all    ,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
       &MPI_COMM_WORLD,info)
#else
  dtambdiff_all=dtambdiff_loc
  dtmagdiff_all=dtmagdiff_loc
  dtwad_all=dtwad_loc
#endif

  dtambdiff(ilevel)=MIN(dtambdiff(ilevel), dtambdiff_all)
  dtmagdiff(ilevel)=MIN(dtmagdiff(ilevel), dtmagdiff_all)
  dtwad(ilevel)=MIN(dtwad(ilevel), dtwad_all) 
  
  ! timestep reduction due to ambipolar diffusion
  if(nambipolar) then
     ! WARNING this should not be done for tests
     if (nminitimestep) then
        ! alfven time alone maybe not correct
        ! comparison with global time step
        tmag1=max(dtambdiff(ilevel),dtwad(ilevel)*coefalfven)
     else
        tmag1=dtambdiff(ilevel)
     endif
     dtnew(ilevel)=MIN(dtnew(ilevel),tmag1)
   endif

  ! timestep reduction due to Ohmic diffusion
  if(nmagdiffu) then
     if (nminitimestep) then
        ! alfven time alone maybe not correct
        ! comparison with global time step
        tmag2=max(dtmagdiff(ilevel),dtwad(ilevel)*coefdtohm)
     else
        tmag2=dtmagdiff(ilevel)
     endif
     dtnew(ilevel)=MIN(dtnew(ilevel),tmag2)
   end if
#endif
  
111 format('   Entering courant_fine for level ',I2)

end subroutine courant_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmpdt(uu,gg,dx,dt,ncell)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none
  integer::ncell
  real(dp)::dx,dt
  real(dp),dimension(1:nvector,1:nvar+3)::uu
  real(dp),dimension(1:nvector,1:ndim)::gg
  real(dp),dimension(1:nvector),save::a2,B2,rho,ctot

  real(dp)::dtcell,smallp,cf,cc,bc,bn
  integer::k,idim
#if NENER>0
  integer::irad
#endif

  smallp = smallr*smallc**2/gamma

  ! Convert to primitive variables
  do k = 1,ncell
     uu(k,1)=max(uu(k,1),smallr)
     rho(k)=uu(k,1)
  end do
  do idim = 1,3
     do k = 1, ncell
        uu(k,idim+1) = uu(k,idim+1)/rho(k)
     end do
  end do

  do k = 1,ncell
     B2(k)=zero
  end do
  do idim = 1,3
     do k = 1, ncell
        Bc = half*(uu(k,neul+idim)+uu(k,nvar+idim))
        B2(k)=B2(k)+Bc**2
        uu(k,neul) = uu(k,neul)-half*uu(k,1)*uu(k,idim+1)**2-half*Bc**2
     end do
  end do
#if NENER>0
  do irad = 1,nener
     do k = 1, ncell
        uu(k,neul) = uu(k,neul)-uu(k,nhydro+irad)
     end do
  end do
#endif

  ! Compute thermal sound speed
  do k = 1, ncell
     uu(k,neul) = max((gamma-one)*uu(k,neul),smallp)
     a2(k)=gamma*uu(k,neul)/uu(k,1)
  end do
#if NENER>0
  do irad = 1,nener
     do k = 1, ncell
        a2(k) = a2(k) + gamma_rad(irad)*(gamma_rad(irad)-1)*uu(k,nhydro+irad)/uu(k,1)
     end do
  end do
#endif

  ! Compute maximum wave speed (fast magnetosonic)
  do k = 1, ncell
     ctot(k)=zero
  end do
  if(ischeme.eq.1)then
     do idim = 1,ndim   ! WARNING: ndim instead of 3
        do k = 1, ncell
           ctot(k)=ctot(k)+abs(uu(k,idim+1))
        end do
     end do
  else
     do idim = 1,ndim   ! WARNING: ndim instead of 3
        do k = 1, ncell
           cc=half*(B2(k)/rho(k)+a2(k))
           BN=half*(uu(k,5+idim)+uu(k,nvar+idim))
           cf=sqrt(cc+sqrt(cc**2-a2(k)*BN**2/rho(k)))
           ctot(k)=ctot(k)+abs(uu(k,idim+1))+cf
        end do
     end do
  endif

  ! Compute gravity strength ratio
  do k = 1, ncell
     rho(k)=zero
  end do
  do idim = 1,ndim
     do k = 1, ncell
        rho(k)=rho(k)+abs(gg(k,idim))
     end do
  end do
  do k = 1, ncell
     rho(k)=rho(k)*dx/ctot(k)**2
     rho(k)=MAX(rho(k),0.0001_dp)
  end do

  ! Compute maximum time step for each authorized cell
  dt = courant_factor*dx/smallc
  do k = 1,ncell
     dtcell=dx/ctot(k)*(sqrt(one+two*courant_factor*rho(k))-one)/rho(k)
     dt = min(dt,dtcell)
  end do

end subroutine cmpdt
!###########################################################
!###########################################################
!###########################################################
!###########################################################
! Non-ideal MHD timesteps
#ifdef NIMHD
subroutine cmpdt_nimhd(uu,dx,ncell,dtambdiff,dtohmdiss)
   use amr_parameters
   use hydro_parameters
   use nimhd_parameters

   implicit none
   integer::ncell
   real(dp)::dx
   real(dp)::dtambdiff,dtohmdiss    ! ambipolar and Ohmic diffusion times
   real(dp),dimension(1:nvector,1:nvar+3)::uu

   real(dp)::xx,betaad,etaohmdiss
   real(dp),dimension(1:nvector),save::B2,rho,tcell
   integer::k,idim
   integer :: ht
  
   do k = 1,ncell
      rho(k)=max(uu(k,1),smallr)
      call temperature_eos(rho(k), uu(k,nvar), tcell(k),ht)
   end do
 
   do k = 1,ncell
      B2(k)=0
   end do
   do idim = 1,3
      do k = 1, ncell
         B2(k)=B2(k) + (0.5d0*(uu(k,5+idim)+uu(k,nvar+idim)))**2
      end do
   end do
  
   ! Ohmic dissipation
   dtohmdiss=1d35
   if (nmagdiffu) then
      do k = 1,ncell
         xx=etaohmdiss(rho(k),B2(k),tcell(k),0d0,0d0,.false.)
         if(xx.gt.0d0) then
            dtohmdiss=min(dtohmdiss, coefohm*dx*dx/xx)
         endif
      end do
   endif
   
   ! ambipolar diffusion
   dtambdiff=1d36 !TC: can we put HUGE or something here?
   if (nambipolar) then
      do k = 1,ncell
         xx=B2(k)*betaad(rho(k),rho(k),0.0d0,B2(k),B2(k),0,tcell(k),.false.)
         if (xx.gt.0d0) then
            dtambdiff=min(dtambdiff, coefad*dx*dx/xx)
         endif
      end do
   endif
 
end subroutine cmpdt_nimhd
#endif
!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine velocity_fine(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !----------------------------------------------------------
  ! This routine computes the gravitational acceleration,
  ! the maximum density rho_max, and the potential energy
  !----------------------------------------------------------
  integer::igrid,ngrid,ncache,i,ind,iskip,ix,iy,iz
  integer::nx_loc,idim
  real(dp)::dx,dx_loc,scale,d,u,v,w,A,B,C
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:3)::skip_loc

  integer ,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp),dimension(1:nvector,1:3),save::vv

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
  ! Compute analytical velocity field
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

        ! Impose analytical velocity field
        select case (condinit_kind)
            case('ponomarenko')
               call velana_ponomarenko(xx,vv,dx_loc,t,ngrid)
            case('nimhd_diffusion')
               call velana_nimhd_diffusion(xx,vv,dx_loc,t,ngrid)
            case default
               call velana(xx,vv,dx_loc,t,ngrid)
         end select

        ! Impose induction variables
        do i=1,ngrid
           uold(ind_cell(i),1)=1.0
        end do
        do idim=1,3
           do i=1,ngrid
              uold(ind_cell(i),idim+1)=vv(i,idim)
           end do
        end do
        ! Update total energy
        do i=1,ngrid
           d=uold(ind_cell(i),1)
           u=uold(ind_cell(i),2)/d
           v=uold(ind_cell(i),3)/d
           w=uold(ind_cell(i),4)/d
           A=0.5*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
           B=0.5*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
           C=0.5*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))
           uold(ind_cell(i),neul)=1.0/(gamma-1.)+0.5*d*(u**2+v**2+w**2)+0.5*(A**2+B**2+C**2)
        end do

     end do
     ! End loop over cells

  end do
  ! End loop over grids

end subroutine velocity_fine
!#########################################################
!#########################################################
!#########################################################
!#########################################################
