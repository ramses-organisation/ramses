!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine godunov_fine(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine is a wrapper to the second order Godunov solver.
  ! Small grids (2x2x2) are gathered from level ilevel and sent to the
  ! hydro solver. On entry, hydro variables are gathered from array uold.
  ! On exit, unew has been updated.
  !--------------------------------------------------------------------------
  integer::i,igrid,ncache,ngrid
  integer,dimension(1:nvector),save::ind_grid

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Loop over active grids by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     call godfine1(ind_grid,ngrid,ilevel)
  end do
  call make_virtual_reverse_dp(divu(1),ilevel)

111 format('   Entering godunov_fine for level ',i2)

end subroutine godunov_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine set_unew(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine sets array unew to its initial value uold before calling
  ! the hydro scheme. unew is set to zero in virtual boundaries.
  !--------------------------------------------------------------------------
  integer::i,ivar,ind,icpu,iskip
  real(dp)::d,u,v,w,e

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Set unew to uold for myid cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do ivar=1,nvar
        do i=1,active(ilevel)%ngrid
           unew(active(ilevel)%igrid(i)+iskip,ivar) = uold(active(ilevel)%igrid(i)+iskip,ivar)
        end do
     end do
     do i=1,active(ilevel)%ngrid
        divu(active(ilevel)%igrid(i)+iskip) = 0
     end do
  end do

  ! Set unew to 0 for virtual boundary cells
  do icpu=1,ncpu
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do ivar=1,nvar
        do i=1,reception(icpu,ilevel)%ngrid
           unew(reception(icpu,ilevel)%igrid(i)+iskip,ivar)=0
        end do
     end do
     do i=1,reception(icpu,ilevel)%ngrid
        divu(reception(icpu,ilevel)%igrid(i)+iskip) = 0
     end do
  end do
  end do

111 format('   Entering set_unew for level ',i2)

end subroutine set_unew
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine set_uold(ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ilevel
  !---------------------------------------------------------
  ! This routine sets array uold to its new value unew
  ! after the hydro step.
  !---------------------------------------------------------
  integer::i,ivar,ind,iskip

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Update unew using non-conservative source terms
  call add_pdv_source_terms(ilevel)

  ! Add gravity source terms to unew
  if(poisson)then
     call add_gravity_source_terms(ilevel)
  end if

  ! Set uold to unew for myid cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     if(static)then
     do ivar=1,nvar
        do i=1,active(ilevel)%ngrid
        ! Embedded body is material #1
        if(.not. (uold(active(ilevel)%igrid(i)+iskip,npri+1)>0.01))then
           uold(active(ilevel)%igrid(i)+iskip,ivar) = unew(active(ilevel)%igrid(i)+iskip,ivar)
        endif
        end do
     end do
     else
     do ivar=1,nvar
        do i=1,active(ilevel)%ngrid
           uold(active(ilevel)%igrid(i)+iskip,ivar) = unew(active(ilevel)%igrid(i)+iskip,ivar)
        end do
     end do
     endif
  end do

111 format('   Entering set_uold for level ',i2)

end subroutine set_uold
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine add_pdv_source_terms(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !-------------------------------------------------------------------
  ! Update volume fraction using compressibility source terms
  !-------------------------------------------------------------------
  integer::i,ivar,imat,idim,ind,iskip,ncache,igrid,ngrid
  logical,dimension(1:nvector),save::body
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp),dimension(1:nvector),save::pp,cc,ekin,kappa_hat,ffmax
  real(dp),dimension(1:nvector,1:npri),save::qq
  real(dp),dimension(1:nvector,1:nmat),save::ff,gg,fg,kappa_mat
  real(dp)::g0,p0,a0,b0,df_over_f,rloc,skip_loc,dx,eps,scale,dx_loc
  real(dp)::kappa_max
  real(dp)::one=1.0_dp, half=0.5_dp, zero=0.0_dp
  real(dp),dimension(1:8)::xc
  integer::ix,iy,iz,nx_loc
  logical::error

  dx=0.5d0**ilevel
  skip_loc=dble(icoarse_min)
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

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

     ! Gather nvector grids
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do

     ! Loop over cells
     do ind=1,twotondim

        ! Compute cell index
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=ind_grid(i)+iskip
        end do

        ! Source terms for volume fraction (Godunov-like advection)
        do imat=1,nmat
           ivar=npri+imat
           do i=1,ngrid
              unew(ind_cell(i),ivar)=unew(ind_cell(i),ivar)-uold(ind_cell(i),ivar)*divu(ind_cell(i))
           end do
        end do

        ! Volume fraction and fluid density
        do imat=1,nmat
           do i=1,ngrid
              ff(i,imat)=uold(ind_cell(i),imat+npri)
              gg(i,imat)=uold(ind_cell(i),imat+npri+nmat)
           end do
        end do
        ! Total density
        do i=1,ngrid
           qq(i,1)=uold(ind_cell(i),1)
        end do
        ! Specific kinetic energy
        ekin(1:ngrid)=zero
        do idim=1,ndim
           do i=1,ngrid
              qq(i,idim+1)=uold(ind_cell(i),idim+1)/qq(i,1)
              ekin(i)=ekin(i)+half*qq(i,idim+1)**2
           end do
        end do
        ! Total internal energy
        do i=1,ngrid
           qq(i,npri)=uold(ind_cell(i),npri)-qq(i,1)*ekin(i)
        end do

        ! Pressure from eos
        call eos(ff,gg,qq,pp,cc,kappa_mat,kappa_hat,ngrid)

        ! Source terms for fluid density (Godunov-like advection)
        do imat=1,nmat
           ivar=npri+nmat+imat
           do i=1,ngrid
              ! Material compressibility
              df_over_f =  one !kappa_hat(i)/kappa_mat(i,imat)

!!$           ! Explicit time integration
!!$           unew(ind_cell(i),ivar)=unew(ind_cell(i),ivar)* &
!!$                & (one+divu(ind_cell(i))*(df_over_f-one))

              ! Implicit time integration
              unew(ind_cell(i),ivar)=unew(ind_cell(i),ivar)/ &
                   & (one-divu(ind_cell(i))*(df_over_f-one))

!!$           ! Exponential time integration
!!$           unew(ind_cell(i),ivar)=unew(ind_cell(i),ivar)* &
!!$                & exp(+divu(ind_cell(i))*(df_over_f-one))

           end do
        end do

        ! Source terms for volume fraction (Godunov-like advection)
        do imat=1,nmat
           ivar=npri+imat
           do i=1,ngrid
              ! No compressibility in volume fraction
              df_over_f  = one !kappa_hat(i)/kappa_mat(i,imat)

              ! Explicit time integration
              unew(ind_cell(i),ivar)=unew(ind_cell(i),ivar)* &
                   & (one-divu(ind_cell(i))*(df_over_f-one))

!!$           ! Implicit time integration
!!$           unew(ind_cell(i),ivar)=unew(ind_cell(i),ivar)/ &
!!$                & (one+divu(ind_cell(i))*(df_over_f-one))

!!$           ! Exponential time integration
!!$           unew(ind_cell(i),ivar)=unew(ind_cell(i),ivar)* &
!!$                & exp(-divu(ind_cell(i))*(df_over_f-one))

           end do
        end do

     end do
     ! End loop over cells
  end do
  ! End loop over grids

end subroutine add_pdv_source_terms
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine add_gravity_source_terms(ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine adds to unew the gravity source terms
  ! with only half a time step. Only the momentum and the
  ! total energy are modified in array unew.
  !--------------------------------------------------------------------------
  integer::i,ind,iskip,ind_cell
  real(dp)::d,u,v,w,e_kin,e_prim,d_old,fact

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Add gravity source term at time t with half time step
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,active(ilevel)%ngrid
        ind_cell=active(ilevel)%igrid(i)+iskip
        d=max(unew(ind_cell,1),smallr)
        u=0; v=0; w=0
        if(ndim>0)u=unew(ind_cell,2)/d
        if(ndim>1)v=unew(ind_cell,3)/d
        if(ndim>2)w=unew(ind_cell,4)/d
        e_kin=0.5d0*d*(u**2+v**2+w**2)
        e_prim=unew(ind_cell,ndim+2)-e_kin
        d_old=max(uold(ind_cell,1),smallr)
        fact=d_old/d*0.5d0*dtnew(ilevel)
        if(ndim>0)then
           u=u+f(ind_cell,1)*fact
           unew(ind_cell,2)=d*u
        endif
        if(ndim>1)then
           v=v+f(ind_cell,2)*fact
           unew(ind_cell,3)=d*v
        end if
        if(ndim>2)then
           w=w+f(ind_cell,3)*fact
           unew(ind_cell,4)=d*w
        endif
        e_kin=0.5d0*d*(u**2+v**2+w**2)
        unew(ind_cell,ndim+2)=e_prim+e_kin
     end do
  end do

111 format('   Entering add_gravity_source_terms for level ',i2)

end subroutine add_gravity_source_terms
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine godfine1(ind_grid,ncache,ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ilevel,ncache
  integer,dimension(1:nvector)::ind_grid
  !-------------------------------------------------------------------
  ! This routine gathers first hydro variables from neighboring grids
  ! to set initial conditions in a 6x6x6 grid. It interpolate from
  ! coarser level missing grid variables. It then calls the
  ! Godunov solver that computes fluxes. These fluxes are zeroed at
  ! coarse-fine boundaries, since contribution from finer levels has
  ! already been taken into account. Conservative variables are updated
  ! and stored in array unew(:), both at the current level and at the
  ! coarser level if necessary.
  !-------------------------------------------------------------------
  integer ,dimension(1:nvector,1:threetondim     ),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim       ),save::nbors_father_grids
  integer ,dimension(1:nvector,0:twondim         ),save::ibuffer_father
  real(dp),dimension(1:nvector,0:twondim  ,1:nvar),save::u1
  real(dp),dimension(1:nvector,0:twondim  ,1:ndim),save::g1
  real(dp),dimension(1:nvector,1:twotondim,1:nvar),save::u2
  real(dp),dimension(1:nvector,1:twotondim,1:ndim),save::g2

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar),save::uloc
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim),save::gloc
  real(dp),dimension(1:nvector,iu1:iu2),save::rloc,wloc_left,wloc_right
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:nvar,1:ndim),save::flux
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:2,1:ndim),save::tmp
  logical ,dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::ok

  integer,dimension(1:nvector),save::igrid_nbor,ind_cell,ind_buffer
  logical,dimension(1:nvector),save::exist_nbor

  integer::i,j,ivar,idim,ind_son,ind_father,iskip,nbuffer,ibuffer,nx_loc
  integer::i0,j0,k0,i1,j1,k1,i2,j2,k2,i3,j3,k3
  integer::i1min,i1max,j1min,j1max,k1min,k1max
  integer::i2min,i2max,j2min,j2max,k2min,k2max
  integer::i3min,i3max,j3min,j3max,k3min,k3max
  real(dp)::dx,eps,scale,skip_loc,dx_loc

  ! Rescaling factors
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=dble(icoarse_min)
  scale=boxlen/dble(nx_loc)

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  dx_loc=dx*scale

  ! Integer constants
  i1min=0; i1max=0; i2min=0; i2max=0; i3min=1; i3max=1
  j1min=0; j1max=0; j2min=0; j2max=0; j3min=1; j3max=1
  k1min=0; k1max=0; k2min=0; k2max=0; k3min=1; k3max=1
  if(ndim>0)then
     i1max=2; i2max=1; i3max=2
  end if
  if(ndim>1)then
     j1max=2; j2max=1; j3max=2
  end if
  if(ndim>2)then
     k1max=2; k2max=1; k3max=2
  end if

  !------------------------------------------
  ! Gather 3^ndim neighboring father cells
  !------------------------------------------
  do i=1,ncache
     ind_cell(i)=father(ind_grid(i))
  end do
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ncache,ilevel)

  !---------------------------
  ! Gather 6x6x6 cells stencil
  !---------------------------
  gloc=0.0d0
  ! Loop over 3x3x3 neighboring father cells
  do k1=k1min,k1max
  do j1=j1min,j1max
  do i1=i1min,i1max

     ! Check if neighboring grid exists
     ind_father=1+i1+3*j1+9*k1
     do i=1,ncache
        igrid_nbor(i)=son(nbors_father_cells(i,ind_father))
        exist_nbor(i)=igrid_nbor(i)>0
     end do

     ! If not, interpolate variables from parent cells
     nbuffer=0
     do i=1,ncache
        if(.not. exist_nbor(i))then
           nbuffer=nbuffer+1
           ind_buffer(nbuffer)=nbors_father_cells(i,ind_father)
        end if
     end do
     call getnborfather(ind_buffer,ibuffer_father,nbuffer,ilevel)
     do j=0,twondim
        do ivar=1,nvar
           do i=1,nbuffer
              u1(i,j,ivar)=uold(ibuffer_father(i,j),ivar)
           end do
        end do
     end do
     call interpol_hydro(u1,u2,nbuffer)

     ! Loop over 2x2x2 cells
     do k2=k2min,k2max
     do j2=j2min,j2max
     do i2=i2min,i2max

        ind_son=1+i2+2*j2+4*k2
        iskip=ncoarse+(ind_son-1)*ngridmax
        do i=1,ncache
           ind_cell(i)=iskip+igrid_nbor(i)
        end do

        i3=1; j3=1; k3=1
        if(ndim>0)i3=1+2*(i1-1)+i2
        if(ndim>1)j3=1+2*(j1-1)+j2
        if(ndim>2)k3=1+2*(k1-1)+k2

        ! Gather hydro variables
        do ivar=1,nvar
           ibuffer=0
           do i=1,ncache
              if(exist_nbor(i))then
                 uloc(i,i3,j3,k3,ivar)=uold(ind_cell(i),ivar)
              else
                 ibuffer=ibuffer+1
                 uloc(i,i3,j3,k3,ivar)=u2(ibuffer,ind_son,ivar)
              end if
           end do
        end do

        ! Gather gravitational acceleration
        if(poisson)then
           do idim=1,ndim
              ibuffer=0
              do i=1,ncache
                 if(exist_nbor(i))then
                    gloc(i,i3,j3,k3,idim)=f(ind_cell(i),idim)
                 else
                    ibuffer=ibuffer+1
                    gloc(i,i3,j3,k3,idim)=f(ibuffer_father(ibuffer,0),idim)
                 end if
              end do
           end do
        end if

        ! Gather refinement flag
        do i=1,ncache
           if(exist_nbor(i))then
              ok(i,i3,j3,k3)=son(ind_cell(i))>0
           else
              ok(i,i3,j3,k3)=.false.
           end if
        end do

     end do
     end do
     end do
     ! End loop over cells

  end do
  end do
  end do
  ! End loop over neighboring grids

  !-----------------------
  ! Compute cell radii
  !-----------------------
  do i3=iu1,iu2
     do i=1,ncache
        rloc(i,i3)=(xg(ind_grid(i),1)-skip_loc+(dble(i3+1)-2.5d0)*dx)*scale
     end do
  end do

  !---------------------------------------------------------
  ! Compute left and right interface geometrical weighting
  !---------------------------------------------------------
  do i3=iu1,iu2
     if(geom==3)then
        do i=1,ncache
           eps=dx_loc/rloc(i,i3)/2.0
           wloc_left (i,i3)=(1.0-eps)**2/(1.0+eps**2/3.0)
           wloc_right(i,i3)=(1.0+eps)**2/(1.0+eps**2/3.0)
        end do
     else if (geom==2) then
        do i=1,ncache
           eps=dx_loc/rloc(i,i3)/2.0
           wloc_left (i,i3)=1.0-eps
           wloc_right(i,i3)=1.0+eps
        end do
     else
        do i=1,ncache
           wloc_left (i,i3)=1.0
           wloc_right(i,i3)=1.0
        end do
     endif
  end do

  !-----------------------------------------------
  ! Compute flux using second-order Godunov method
  !-----------------------------------------------
  flux=0.0d0; tmp=0.0d0
  call unsplit(uloc,gloc,rloc,flux,tmp,dx_loc,dx_loc,dx_loc,dtnew(ilevel),ncache)

  !------------------------------------------------
  ! Reset flux along direction at refined interface
  !------------------------------------------------
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     do k3=k3min,k3max+k0
     do j3=j3min,j3max+j0
     do i3=i3min,i3max+i0
        do ivar=1,nvar
           do i=1,ncache
              if(ok(i,i3-i0,j3-j0,k3-k0) .or. ok(i,i3,j3,k3))then
                 flux(i,i3,j3,k3,ivar,idim)=0.0d0
              end if
           end do
        end do
        do ivar=1,2
           do i=1,ncache
              if(ok(i,i3-i0,j3-j0,k3-k0) .or. ok(i,i3,j3,k3))then
                 tmp (i,i3,j3,k3,ivar,idim)=0.0d0
              end if
           end do
        end do
     end do
     end do
     end do
  end do
  !--------------------------------------
  ! Conservative update at level ilevel
  !--------------------------------------
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     do k2=k2min,k2max
     do j2=j2min,j2max
     do i2=i2min,i2max
        ind_son=1+i2+2*j2+4*k2
        iskip=ncoarse+(ind_son-1)*ngridmax
        do i=1,ncache
           ind_cell(i)=iskip+ind_grid(i)
        end do
        i3=1+i2
        j3=1+j2
        k3=1+k2
        ! Update conservative variables new state vector and velocity divergence
        if(geom>1.and.idim==1)then
           ! Use geometrical weighting
           do ivar=1,nvar
              do i=1,ncache
                 unew(ind_cell(i),ivar) = unew(ind_cell(i),ivar) + &
                      & (flux(i,i3   ,j3   ,k3   ,ivar,idim)*wloc_left (i,i3) &
                      & -flux(i,i3+i0,j3+j0,k3+k0,ivar,idim)*wloc_right(i,i3))
              end do
           end do
           do i=1,ncache
              divu(ind_cell(i))=divu(ind_cell(i))+ &
                   & (tmp(i,i3   ,j3   ,k3   ,1,idim)*wloc_left (i,i3) &
                   & -tmp(i,i3+i0,j3+j0,k3+k0,1,idim)*wloc_right(i,i3))
           end do
           ! Add -grad(P) term for radial momentum
           do i=1,ncache
              unew(ind_cell(i),2) = unew(ind_cell(i),2) + &
                   & (tmp (i,i3   ,j3   ,k3   ,2,idim) &
                   & -tmp (i,i3+i0,j3+j0,k3+k0,2,idim))
           end do
        else
           do ivar=1,nvar
              do i=1,ncache
                 unew(ind_cell(i),ivar)=unew(ind_cell(i),ivar)+ &
                      & (flux(i,i3   ,j3   ,k3   ,ivar,idim) &
                      & -flux(i,i3+i0,j3+j0,k3+k0,ivar,idim))
              end do
           end do
           ! Update velocity divergence
           do i=1,ncache
              divu(ind_cell(i))=divu(ind_cell(i))+ &
                   & (tmp(i,i3   ,j3   ,k3   ,1,idim) &
                   & -tmp(i,i3+i0,j3+j0,k3+k0,1,idim))
           end do
        endif
     end do
     end do
     end do
  end do

  !--------------------------------------
  ! Conservative update at level ilevel-1
  !--------------------------------------
  ! Loop over dimensions
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1

     !----------------------
     ! Left flux at boundary
     !----------------------
     ! Check if grids sits near left boundary
     do i=1,ncache
        exist_nbor(i)=son(nbor(ind_grid(i),2*idim-1))>0
     end do
     ! Gather neighbor father cells index
     do i=1,ncache
        ind_buffer(i)=nbor(ind_grid(i),2*idim-1)
     end do
     ! Conservative update of new state variables
     do ivar=1,nvar
        ! Loop over boundary cells
        do k3=k3min,k3max-k0
        do j3=j3min,j3max-j0
        do i3=i3min,i3max-i0
           do i=1,ncache
              if(.not.exist_nbor(i))then
                 unew(ind_buffer(i),ivar)=unew(ind_buffer(i),ivar) &
                      & -flux(i,i3,j3,k3,ivar,idim)/dble(twotondim)
              end if
           end do
        end do
        end do
        end do
     end do
     ! Update velocity divergence
     do k3=k3min,k3max-k0
     do j3=j3min,j3max-j0
     do i3=i3min,i3max-i0
        do i=1,ncache
           if(.not.exist_nbor(i))then
              divu(ind_buffer(i))=divu(ind_buffer(i)) &
                   & -tmp(i,i3,j3,k3,1,idim)/dble(twotondim)
           end if
        end do
     end do
     end do
     end do

     !-----------------------
     ! Right flux at boundary
     !-----------------------
     ! Check if grids sits near right boundary
     do i=1,ncache
        exist_nbor(i)=son(nbor(ind_grid(i),2*idim))>0
     end do
     ! Gather buffer indices
     do i=1,ncache
        ind_buffer(i)=nbor(ind_grid(i),2*idim)
     end do
     ! Conservative update of new state variables
     do ivar=1,nvar
        ! Loop over boundary cells
        do k3=k3min+k0,k3max
        do j3=j3min+j0,j3max
        do i3=i3min+i0,i3max
           do i=1,ncache
              if(.not.exist_nbor(i))then
                 unew(ind_buffer(i),ivar)=unew(ind_buffer(i),ivar) &
                      & +flux(i,i3+i0,j3+j0,k3+k0,ivar,idim)/dble(twotondim)
              end if
           end do
        end do
        end do
        end do
     end do
     ! Update velocity divergence
     do k3=k3min+k0,k3max
     do j3=j3min+j0,j3max
     do i3=i3min+i0,i3max
        do i=1,ncache
           if(.not.exist_nbor(i))then
              divu(ind_buffer(i))=divu(ind_buffer(i)) &
                   & +tmp(i,i3+i0,j3+j0,k3+k0,1,idim)/dble(twotondim)
           end if
        end do
     end do
     end do
     end do

  end do
  ! End loop over dimensions

end subroutine godfine1
