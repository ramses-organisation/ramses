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
  if(static)return
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
#if NENER>0
  integer::irad
#endif

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Set unew to uold for myid cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do ivar=1,nvar_all
        do i=1,active(ilevel)%ngrid
           unew(active(ilevel)%igrid(i)+iskip,ivar) = uold(active(ilevel)%igrid(i)+iskip,ivar)
        end do
     end do
     if(momentum_feedback>0)then
        do i=1,active(ilevel)%ngrid
           pstarnew(active(ilevel)%igrid(i)+iskip) = 0
        end do
     endif
     if(pressure_fix)then
        do i=1,active(ilevel)%ngrid
           divu(active(ilevel)%igrid(i)+iskip) = 0
        end do
        do i=1,active(ilevel)%ngrid
           d=max(uold(active(ilevel)%igrid(i)+iskip,1),smallr)
           u=0; v=0; w=0
           if(ndim>0)u=uold(active(ilevel)%igrid(i)+iskip,2)/d
           if(ndim>1)v=uold(active(ilevel)%igrid(i)+iskip,3)/d
           if(ndim>2)w=uold(active(ilevel)%igrid(i)+iskip,4)/d
           e=uold(active(ilevel)%igrid(i)+iskip,neul)-0.5d0*d*(u**2+v**2+w**2)
#if NENER>0
           do irad=1,nener
              e=e-uold(active(ilevel)%igrid(i)+iskip,nhydro+irad)
           end do
#endif
           enew(active(ilevel)%igrid(i)+iskip)=e
        end do
     end if
  end do

  ! Set unew to 0 for virtual boundary cells
  do icpu=1,ncpu
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do ivar=1,nvar_all
        do i=1,reception(icpu,ilevel)%ngrid
#ifdef LIGHT_MPI_COMM
           unew(reception(icpu,ilevel)%pcomm%igrid(i)+iskip,ivar)=0
#else
           unew(reception(icpu,ilevel)%igrid(i)+iskip,ivar)=0
#endif
        end do
     end do
     if(momentum_feedback>0)then
        do i=1,reception(icpu,ilevel)%ngrid
#ifdef LIGHT_MPI_COMM
           pstarnew(reception(icpu,ilevel)%pcomm%igrid(i)+iskip) = 0
#else
           pstarnew(reception(icpu,ilevel)%igrid(i)+iskip) = 0
#endif
        end do
     endif
     if(pressure_fix)then
        do i=1,reception(icpu,ilevel)%ngrid
#ifdef LIGHT_MPI_COMM
           divu(reception(icpu,ilevel)%pcomm%igrid(i)+iskip) = 0
           enew(reception(icpu,ilevel)%pcomm%igrid(i)+iskip) = 0
#else
           divu(reception(icpu,ilevel)%igrid(i)+iskip) = 0
           enew(reception(icpu,ilevel)%igrid(i)+iskip) = 0
#endif
        end do
     end if
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
  integer::i,ivar,ind,iskip,nx_loc,ind_cell
  real(dp)::scale,d,u,v,w
  real(dp)::e_kin,e_cons,e_prim,e_trunc,div,dx
#if NENER>0
  integer::irad
#endif

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5d0**ilevel*scale

  ! Set uold to unew for myid cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax

     ! -------------------------------------------------------------------------------------------------------------------------------------------------------------
     ! L. Romano 13.06.2023 -- Catch advection errors due to smallr
#if NVAR > NHYDRO+NENER
     do i=1,active(ilevel)%ngrid
        if(uold(active(ilevel)%igrid(i)+iskip,1).lt.smallr.and.unew(active(ilevel)%igrid(i)+iskip,1).gt.uold(active(ilevel)%igrid(i)+iskip,1))then
           ! inflow into previously floored cell: fix concentrations
           do ivar = nhydro+1+nener, nvar
              unew(active(ilevel)%igrid(i)+iskip,ivar) = uold(active(ilevel)%igrid(i)+iskip,ivar) * max(unew(active(ilevel)%igrid(i)+iskip, 1), smallr) / smallr
           end do
        else if(unew(active(ilevel)%igrid(i)+iskip,1).lt.smallr.and.uold(active(ilevel)%igrid(i)+iskip,1).gt.unew(active(ilevel)%igrid(i)+iskip,1))then
           ! outflow leading to density below floor: apply density floor to scalar density
           do ivar = nhydro+1+nener, nvar
              unew(active(ilevel)%igrid(i)+iskip,ivar) = uold(active(ilevel)%igrid(i)+iskip,ivar) * smallr / max(uold(active(ilevel)%igrid(i)+iskip, 1), smallr)
           end do
        end if
     end do
#endif
     ! -------------------------------------------------------------------------------------------------------------------------------------------------------------

     do ivar=1,nvar_all
        do i=1,active(ilevel)%ngrid
           uold(active(ilevel)%igrid(i)+iskip,ivar) = unew(active(ilevel)%igrid(i)+iskip,ivar)
        end do
     end do
     if(momentum_feedback>0)then
        do i=1,active(ilevel)%ngrid
           pstarold(active(ilevel)%igrid(i)+iskip) = pstarnew(active(ilevel)%igrid(i)+iskip)
        end do
     endif
     if(pressure_fix)then
        ! Correct total energy if internal energy is too small
        do i=1,active(ilevel)%ngrid
           ind_cell=active(ilevel)%igrid(i)+iskip
           d=max(uold(ind_cell,1),smallr)
           u=0; v=0; w=0
           if(ndim>0)u=uold(ind_cell,2)/d
           if(ndim>1)v=uold(ind_cell,3)/d
           if(ndim>2)w=uold(ind_cell,4)/d
           e_kin=0.5d0*d*(u**2+v**2+w**2)
#if NENER>0
           do irad=1,nener
              e_kin=e_kin+uold(ind_cell,nhydro+irad)
           end do
#endif
           e_cons=uold(ind_cell,neul)-e_kin
           e_prim=enew(ind_cell)
           ! Note: here divu=-div.u*dt
           div=abs(divu(ind_cell))*dx/dtnew(ilevel)
           e_trunc=beta_fix*d*max(div,3.0d0*hexp*dx)**2
           if(e_cons<e_trunc)then
              uold(ind_cell,neul)=e_prim+e_kin
           end if
        end do
     end if
  end do

111 format('   Entering set_uold for level ',i2)

end subroutine set_uold
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
  real(dp)::req=0_dp

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
        e_prim=unew(ind_cell,neul)-e_kin
        d_old=max(uold(ind_cell,1),smallr)
        if(strict_equilibrium>0)req=rho_eq(ind_cell)
        fact=(d_old-req)/d*0.5d0*dtnew(ilevel)
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
        unew(ind_cell,neul)=e_prim+e_kin
     end do
  end do

111 format('   Entering add_gravity_source_terms for level ',i2)

end subroutine add_gravity_source_terms
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine add_pdv_source_terms(ilevel)
  use amr_commons
  use hydro_commons
  use amr_constants, only:iii,jjj
  implicit none
  integer::ilevel
  !---------------------------------------------------------
  ! This routine adds the pdV source term to the internal
  ! energy equation and to the non-thermal energy equations.
  !---------------------------------------------------------
  integer::i,ind,iskip,nx_loc,ind_cell1
  integer::ncache,igrid,ngrid,idim,id1,ig1,ih1,id2,ig2,ih2
  real(dp)::scale,dx,dx_loc,d,u,v,w,eold

  integer ,dimension(1:nvector),save::ind_grid,ind_cell
  integer ,dimension(1:nvector,0:twondim),save::igridn
  integer ,dimension(1:nvector,1:ndim),save::ind_left,ind_right
  real(dp),dimension(1:nvector,1:ndim,1:ndim),save::velg,veld
  real(dp),dimension(1:nvector,1:ndim),save::dx_g,dx_d
  real(dp),dimension(1:nvector),save::divu_loc
#if NENER>0
  integer::irad
#endif

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5d0**ilevel
  dx_loc=dx*scale

  ! Loop over myid grids by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector

     ! Gather nvector grids
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do

     ! Gather neighboring grids
     do i=1,ngrid
        igridn(i,0)=ind_grid(i)
     end do
     do idim=1,ndim
        do i=1,ngrid
           ind_left (i,idim)=nbor(ind_grid(i),2*idim-1)
           ind_right(i,idim)=nbor(ind_grid(i),2*idim  )
           igridn(i,2*idim-1)=son(ind_left (i,idim))
           igridn(i,2*idim  )=son(ind_right(i,idim))
        end do
     end do

     ! Loop over cells
     do ind=1,twotondim

        ! Compute central cell index
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do

        ! Gather all neighboring velocities
        do idim=1,ndim
           id1=jjj(idim,1,ind); ig1=iii(idim,1,ind)
           ih1=ncoarse+(id1-1)*ngridmax
           do i=1,ngrid
              if(igridn(i,ig1)>0)then
                 velg(i,idim,1:ndim) = uold(igridn(i,ig1)+ih1,2:ndim+1)/max(uold(igridn(i,ig1)+ih1,1),smallr)
                 dx_g(i,idim) = dx_loc
              else
                 velg(i,idim,1:ndim) = uold(ind_left(i,idim),2:ndim+1)/max(uold(ind_left(i,idim),1),smallr)
                 dx_g(i,idim) = dx_loc*1.5_dp
              end if
           enddo
           id2=jjj(idim,2,ind); ig2=iii(idim,2,ind)
           ih2=ncoarse+(id2-1)*ngridmax
           do i=1,ngrid
              if(igridn(i,ig2)>0)then
                 veld(i,idim,1:ndim)= uold(igridn(i,ig2)+ih2,2:ndim+1)/max(uold(igridn(i,ig2)+ih2,1),smallr)
                 dx_d(i,idim)=dx_loc
              else
                 veld(i,idim,1:ndim)= uold(ind_right(i,idim),2:ndim+1)/max(uold(ind_right(i,idim),1),smallr)
                 dx_d(i,idim)=dx_loc*1.5_dp
              end if
           enddo
        end do
        ! End loop over dimensions

        ! Compute divu = Trace G
        divu_loc(1:ngrid)=0.0d0
        do i=1,ngrid
           do idim=1,ndim
              divu_loc(i) = divu_loc(i) + (veld(i,idim,idim)-velg(i,idim,idim)) &
                   &                    / (dx_g(i,idim)     +dx_d(i,idim))
           enddo
        end do

        ! Update thermal internal energy
        if(pressure_fix)then
           do i=1,ngrid
              ! Compute old thermal energy
              d=max(uold(ind_cell(i),1),smallr)
              u=0; v=0; w=0
              if(ndim>0)u=uold(ind_cell(i),2)/d
              if(ndim>1)v=uold(ind_cell(i),3)/d
              if(ndim>2)w=uold(ind_cell(i),4)/d
              eold=uold(ind_cell(i),neul)-0.5d0*d*(u**2+v**2+w**2)
#if NENER>0
              do irad=1,nener
                 eold=eold-uold(ind_cell(i),nhydro+irad)
              end do
#endif
              ! Add -pdV term
              enew(ind_cell(i))=enew(ind_cell(i)) &
                   & -(gamma-1.0d0)*eold*divu_loc(i)*dtnew(ilevel)
           end do
        end if

#if NENER>0
        do irad=1,nener
           do i=1,ngrid
              ! Add -pdV term
              unew(ind_cell(i),nhydro+irad)=unew(ind_cell(i),nhydro+irad) &
                & -(gamma_rad(irad)-1.0d0)*uold(ind_cell(i),nhydro+irad)*divu_loc(i)*dtnew(ilevel)
           end do
        end do
#endif

     enddo
     ! End loop over cells
  end do
  ! End loop over grids

  return

  ! This is the old technique based on the "pressure fix" option.

  ! Update thermal internal energy
  if(pressure_fix)then
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,active(ilevel)%ngrid
           ind_cell1=active(ilevel)%igrid(i)+iskip
           ! Compute old thermal energy
           d=max(uold(ind_cell1,1),smallr)
           u=0; v=0; w=0
           if(ndim>0)u=uold(ind_cell1,2)/d
           if(ndim>1)v=uold(ind_cell1,3)/d
           if(ndim>2)w=uold(ind_cell1,4)/d
           eold=uold(ind_cell1,neul)-0.5d0*d*(u**2+v**2+w**2)
#if NENER>0
           do irad=1,nener
              eold=eold-uold(ind_cell1,nhydro+irad)
           end do
#endif
           ! Add pdV term
           enew(ind_cell1)=enew(ind_cell1) &
                & +(gamma-1.0d0)*eold*divu(ind_cell1) ! Note: here divu=-div.u*dt
        end do
     end do
  end if

#if NENER>0
  do irad=1,nener
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,active(ilevel)%ngrid
           ind_cell1=active(ilevel)%igrid(i)+iskip
           unew(ind_cell1,nhydro+irad)=unew(ind_cell1,nhydro+irad) &
                & +(gamma_rad(irad)-1.0d0)*uold(ind_cell1,nhydro+irad)*divu(ind_cell1) ! Note: here divu=-div.u*dt
        end do
     end do
  end do
#endif

111 format('   Entering add_pdv_source_terms for level ',i2)

end subroutine add_pdv_source_terms
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine add_viscosity_source_terms(ilevel)
   use amr_commons
   use hydro_commons
   use poisson_commons
   use pm_commons


   implicit none
   integer::ilevel
   !--------------------------------------------------------------------------
   ! This routine adds to unew the viscosity terms
   ! Only the momentum and the
   ! total energy are modified in array unew.
   !--------------------------------------------------------------------------
   integer::i,ind,iskip,nx_loc,ix,iy,iz, j, iskip_son
   integer::ncache,igrid,ngrid,idim,id1,ig1,ih1,id2,ig2,ih2,jdim
   integer,dimension(1:3,1:2,1:8)::iii,jjj
   integer,dimension(1:8)::i3l,j3l, k3l, i3r, j3r, k3r, i3c, j3c, k3c
   real(kind=8)::scale,dx,dx_loc
   real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2

   integer ,dimension(1:nvector),save::ind_grid,ind_cell
   integer ,dimension(1:nvector,0:twondim),save::igridn
   integer ,dimension(1:nvector,1:ndim),save::ind_left,ind_right
   real(dp),dimension(1:nvector,1:ndim,1:ndim),save::vel_left,vel_right
   real(dp),dimension(1:nvector,1:2,1:2,1:ndim),save::vel_diag
   real(dp),dimension(1:nvector,1:2,1:ndim),save::dx_diag
   real(dp) :: vproml1, vproml2, vpromr1, vpromr2
   integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
   integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
   real(dp),dimension(1:nvector,1:ndim),save::den_left,den_right
   real(dp),dimension(1:nvector,1:ndim),save::dx_left,dx_right
   real(dp),dimension(1:nvector,1:ndim),save::laplacian_du_loc
   real(dp) :: mu_viscosity
   real(dp), dimension(1:ndim) :: vel
   real(dp) :: den
   real(dp) :: den_dvel_left, den_dvel_right ! density times velocity derivative left and right
   real(dp) :: dxf
   real(dp) :: d=0,u=0,v=0,w=0,du=0,dv=0,dw=0,e_kin=0,e_nokin=0
#ifdef SOLVERMHD
   real(dp) :: A,B,C
#endif

   real(dp),dimension(1:nvector,1:ndim),save::x
   real(dp),dimension(1:twotondim,1:3)::xc
   real(dp),dimension(1:3)::skip_loc


   real(dp),dimension(1:2):: mu_viscosity_left, mu_viscosity_right

   integer ,dimension(1:nvector,0:twondim         ),save::ibuffer_father
   real(dp),dimension(1:nvector,0:twondim  ,1:nvar),save::u1
   real(dp),dimension(1:nvector,1:twotondim,1:nvar),save::u2
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar),save::uloc
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim),save::x_diag

   integer::i1,j1,k1,i2,j2,k2,i3,j3,k3,nexist,nbuffer,ind_son,ind_father,ivar
   integer,dimension(1:nvector),save::igrid_nbor,ind_buffer,ind_exist,ind_nexist
   integer::i1min,i1max,j1min,j1max,k1min,k1max
   integer::i2min,i2max,j2min,j2max,k2min,k2max
   integer::i3min,i3max,j3min,j3max,k3min,k3max

   if(numbtot(1,ilevel)==0)return
   if(verbose)write(*,111)ilevel

   nx_loc=icoarse_max-icoarse_min+1
   skip_loc = (/0.0d0, 0.0d0, 0.0d0/)
   if (ndim > 0) skip_loc(1) = dble(icoarse_min)
   if (ndim > 1) skip_loc(2) = dble(jcoarse_min)
   if (ndim > 2) skip_loc(3) = dble(kcoarse_min)
   scale=boxlen/dble(nx_loc)
   dx=0.5D0**ilevel
   dx_loc=dx*scale


   call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

   iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
   iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
   iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
   iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
   iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
   iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)


! These arrays are neccesary to obtain the correct neighbors according to each cell
   i3l(1:8) = (/-1,0,-1,0,-1,0,-1,0/); j3l(1:8) = (/-1,-1,0,0,-1,-1,0,0/);  k3l(1:8) = (/-1,-1,-1,-1,0,0,0,0/)
   i3c(1:8) = (/1,2,1,2,1,2,1,2/);     j3c(1:8) = (/1,1,2,2,1,1,2,2/);      k3c(1:8) = (/1,1,1,1,2,2,2,2/)
   i3r(1:8) = (/3,4,3,4,3,4,3,4/);     j3r(1:8) = (/3,3,4,4,3,3,4,4/);      k3r(1:8) = (/3,3,3,3,4,4,4,4/)


   ! Loop over myid grids by vector sweeps
   ncache=active(ilevel)%ngrid
   do igrid=1,ncache,nvector

      ! Gather nvector grids
      ngrid=MIN(nvector,ncache-igrid+1)
      do i=1,ngrid
         ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
      end do

      ! Gather neighboring grids
      do i=1,ngrid
         igridn(i,0)=ind_grid(i)
      end do
      do idim=1,ndim
         do i=1,ngrid
            ind_left (i,idim)=nbor(ind_grid(i),2*idim-1)
            ind_right(i,idim)=nbor(ind_grid(i),2*idim  )
            igridn(i,2*idim-1)=son(ind_left (i,idim))
            igridn(i,2*idim  )=son(ind_right(i,idim))
         end do
      end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Here begins the diagonal neighbors block of code !!!!!!!!!!!!!!!!!!!!!!!!!!!
!This was taken from the godfine1 subroutine in this same file

     ! Compute father cell index
      do i=1,ngrid
         ind_cell(i)=father(ind_grid(i))
      end do

      call get3cubefather(ind_cell,nbors_father_cells,ngrid,ilevel)

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


      !---------------------------
      ! Gather 6x6x6 cells stencil
      !---------------------------
      ! Loop over 3x3x3 neighboring father cells
      do k1=k1min,k1max
      do j1=j1min,j1max
      do i1=i1min,i1max

         ! Check if neighboring grid exists
         nbuffer=0
         nexist=0
         ind_father=1+i1+3*j1+9*k1
         do i=1,ngrid
            igrid_nbor(i)=son(nbors_father_cells(i,ind_father))
            if(igrid_nbor(i)>0) then
               nexist=nexist+1
               ind_exist(nexist)=i
            else
              nbuffer=nbuffer+1
              ind_nexist(nbuffer)=i
              ind_buffer(nbuffer)=nbors_father_cells(i,ind_father)
            end if
         end do

         ! If not, interpolate hydro variables from parent cells
         if(nbuffer>0)then
            call getnborfather(ind_buffer,ibuffer_father,nbuffer,ilevel)
            do j=0,twondim
               do ivar=1,nvar
                  do i=1,nbuffer
                     u1(i,j,ivar)=uold(ibuffer_father(i,j),ivar)
                  end do
               end do
            end do
            call interpol_hydro(u1,u2,nbuffer)
         endif

         ! Loop over 2x2x2 cells
         do k2=k2min,k2max
         do j2=j2min,j2max
         do i2=i2min,i2max

            ind_son=1+i2+2*j2+4*k2
            iskip=ncoarse+(ind_son-1)*ngridmax
            do i=1,nexist
               ind_cell(i)=iskip+igrid_nbor(ind_exist(i))
            end do


            i3=1; j3=1; k3=1
            if(ndim>0)then
               i3=1+2*(i1-1)+i2
            end if
            if(ndim>1)then
               j3=1+2*(j1-1)+j2
            end if
            if(ndim>2)then
               k3=1+2*(k1-1)+k2
            end if



            ! Gather hydro variables
            do ivar=1,nvar
               do i=1,nexist
                  uloc(ind_exist(i),i3,j3,k3,ivar)=uold(ind_cell(i),ivar)
               end do
               do i=1,nbuffer
                  uloc(ind_nexist(i),i3,j3,k3,ivar)=u2(i,ind_son,ivar)
               end do
            end do
            ! Gather neighboring positions for the dx and neighboring viscosities for the laplacian
            do idim=1,ndim
               !do i=1,nexist
               do i=1,ngrid
                  x_diag(i,i3,j3,k3,idim) = (xg(igrid_nbor(i),idim)+xc(ind_son,idim)-skip_loc(idim))*scale
               end do
            end do

         end do
         end do
         end do
         ! End loop over cells

      end do
      end do
      end do
      ! End loop over neighboring grids
      ! Loop over grid cells
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!End of diagonal indices block !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Loop over cells
      do ind=1,twotondim

         iskip_son=ncoarse+(ind-1)*ngridmax

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Here the left right neighbors are calculated !!!!!!!!!!!!!!!!!!

          ! Compute central cell index
         iskip=ncoarse+(ind-1)*ngridmax
         do i=1,ngrid
            ind_cell(i)=iskip+ind_grid(i)
         end do


         ! Gather all neighboring velocities
         do idim=1,ndim
            id1=jjj(idim,1,ind); ig1=iii(idim,1,ind)
            ih1=ncoarse+(id1-1)*ngridmax
            do i=1,ngrid
               if(igridn(i,ig1)>0)then
                  den_left(i,idim) = max(uold(igridn(i,ig1)+ih1,1),smallr)
                  vel_left(i,idim,1:ndim) = uold(igridn(i,ig1)+ih1,2:ndim+1)/max(uold(igridn(i,ig1)+ih1,1),smallr)
                  dx_left(i,idim) = dx_loc
               else
                  den_left(i,idim) = max(uold(ind_left(i,idim),1),smallr)
                  vel_left(i,idim,1:ndim) = uold(ind_left(i,idim),2:ndim+1)/max(uold(ind_left(i,idim),1),smallr)
                  dx_left(i,idim) = dx_loc*1.5_dp
               end if
            enddo
            id2=jjj(idim,2,ind); ig2=iii(idim,2,ind)
            ih2=ncoarse+(id2-1)*ngridmax
            do i=1,ngrid
               if(igridn(i,ig2)>0)then
                  den_right(i,idim)= max(uold(igridn(i,ig2)+ih2,1),smallr)
                  vel_right(i,idim,1:ndim)= uold(igridn(i,ig2)+ih2,2:ndim+1)/max(uold(igridn(i,ig2)+ih2,1),smallr)
                  dx_right(i,idim)=dx_loc
               else
                  den_right(i,idim)= max(uold(ind_right(i,idim),1),smallr)
                  vel_right(i,idim,1:ndim)= uold(ind_right(i,idim),2:ndim+1)/max(uold(ind_right(i,idim),1),smallr)
                  dx_right(i,idim)=dx_loc*1.5_dp
               end if
            end do
         end do
         ! End loop over dimensions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Here are the diagonal neighbors !!!!!!!!!!!!!!!!!!!!!!!
        ! Gather cell centre positions
         do idim=1,ndim
            do i=1,ngrid
               x(i,idim)=xg(ind_grid(i),idim)+xc(ind,idim)
            end do
         end do

         ! Rescale position from code units to user units
         do idim = 1, ndim
            do i = 1, ngrid
               x(i, idim) = (x(i, idim) - skip_loc(idim))*scale
            end do
         end do

         do i = 1,ngrid

            den_left(i,1)= max(uloc(i,i3l(ind),j3c(ind),1,1),smallr)
            den_left(i,2)= max(uloc(i,i3c(ind),j3l(ind),1,1),smallr)
            vel_left(i,1,1:ndim)= uloc(i,i3l(ind),j3c(ind),1,2:ndim+1)/max(uloc(i,i3l(ind),j3c(ind),1,1),smallr)
            vel_left(i,2,1:ndim)= uloc(i,i3c(ind),j3l(ind),1,2:ndim+1)/max(uloc(i,i3c(ind),j3l(ind),1,1),smallr)

            dx_left(i,1) = dx_loc*2.0_dp
            dx_left(i,2) = dx_loc*2.0_dp



            den_right(i,1)= max(uloc(i,i3r(ind),j3c(ind),1,1),smallr)
            den_right(i,2)= max(uloc(i,i3c(ind),j3r(ind),1,1),smallr)
            vel_right(i,1,1:ndim)= uloc(i,i3r(ind),j3c(ind),1,2:ndim+1)/max(uloc(i,i3r(ind),j3c(ind),1,1),smallr)
            vel_right(i,2,1:ndim)= uloc(i,i3c(ind),j3r(ind),1,2:ndim+1)/max(uloc(i,i3c(ind),j3r(ind),1,1),smallr)

            dx_right(i,1) = dx_loc*2.0_dp
            dx_right(i,2) = dx_loc*2.0_dp


            vel_diag(i,1,1,1:ndim) = uloc(i,i3l(ind),j3l(ind),1,2:ndim+1)/max(uloc(i,i3l(ind),j3l(ind),1,1),smallr)
            vel_diag(i,1,2,1:ndim) = uloc(i,i3l(ind),j3r(ind),1,2:ndim+1)/max(uloc(i,i3l(ind),j3r(ind),1,1),smallr)
            vel_diag(i,2,1,1:ndim) = uloc(i,i3r(ind),j3l(ind),1,2:ndim+1)/max(uloc(i,i3r(ind),j3l(ind),1,1),smallr)
            vel_diag(i,2,2,1:ndim)= uloc(i,i3r(ind),j3r(ind),1,2:ndim+1)/max(uloc(i,i3r(ind),j3r(ind),1,1),smallr)

            dx_diag(i,1,1) = dx_loc*2.0_dp
            dx_diag(i,2,1) = dx_loc*2.0_dp
            dx_diag(i,1,2) = dx_loc*2.0_dp
            dx_diag(i,2,2) = dx_loc*2.0_dp

         end do

! Compute the laplacian (I haven't programed the structure for averaging viscosities with neighbors, I hope to do that once the laplacian is working properly)
         laplacian_du_loc(1:ngrid,1:ndim)=0.0d0
         do i = 1,ngrid

            den = max(uold(ind_cell(i),1), smallr) ! density of the cell
            vel(1:ndim) = uold(ind_cell(i), 2:ndim+1) / max(uold(ind_cell(i),1), smallr) ! velocity of the cell

           ! Viscosity of the cell
            select case (viscosity_kind)
               case('constant_uniform')

                  mu_viscosity = mu_viscosity_constant
                  mu_viscosity_left(1) = mu_viscosity_constant
                  mu_viscosity_left(2) = mu_viscosity_constant
                  mu_viscosity_right(1) = mu_viscosity_constant
                  mu_viscosity_right(2) = mu_viscosity_constant

            end select


            ! Add non crossed terms (d_i (sigma nu d_i v_j ))
            do jdim=1,ndim ! component of the laplacian and the velocity
               do idim=1,ndim ! direction for derivatives

                  den_dvel_left  =  ((den+den_left(i,idim))/2.0 ) * (mu_viscosity+mu_viscosity_left(idim))/2.0 * ((vel(jdim) - vel_left(i,idim,jdim) ) / dx_left(i,idim) )!*(1/2.0)
                  den_dvel_right =  ((den+den_right(i,idim))/2.0) * (mu_viscosity+mu_viscosity_right(idim))/2.0 * ((vel_right(i,idim,jdim)-vel(jdim) ) /  dx_right(i,idim))!*(1/2.0)

                  dxf = (dx_left(i,idim)+dx_right(i,idim))/2.0

                  laplacian_du_loc(i,jdim) = laplacian_du_loc(i,jdim) + (den_dvel_right - den_dvel_left) /  dxf
               end do
            end do

         ! Add crossed terms to the laplacian

            den = max(uold(ind_cell(i),1), smallr) ! density of the cell
            vel(1:ndim) = uold(ind_cell(i), 2:ndim+1) / max(uold(ind_cell(i),1), smallr) ! velocity of the cell

            ! Crossed terms d_y (sigma nu d_x v_i )

            do idim=1,ndim

               vproml1 = (vel(idim)+vel_left(i,2,idim)+vel_left(i,1,idim)+vel_diag(i,1,1,idim))/4.0
               vproml2 = (vel(idim)+vel_left(i,2,idim)+vel_right(i,1,idim)+vel_diag(i,2,1,idim))/4.0

               vpromr1 = (vel(idim)+vel_right(i,2,idim)+vel_left(i,1,idim)+vel_diag(i,1,2,idim))/4.0
               vpromr2 = (vel(idim)+vel_right(i,2,idim)+vel_right(i,1,idim)+vel_diag(i,2,2,idim))/4.0

               den_dvel_left  =  ( (den+den_left(i,2))/2.0 ) * (mu_viscosity+mu_viscosity_left(2))/2.0 * ( ( vproml2 - vproml1 )/dx_diag(i,1,1) )! *(1/2.0)
               den_dvel_right =  ( (den+den_right(i,2))/2.0 ) * (mu_viscosity+mu_viscosity_right(2))/2.0 * ( (vpromr2 - vpromr1 ) /dx_diag(i,2,1) )!*(1/2.0)

               dxf = (dx_left(i,2)+dx_right(i,2))/2.0

               laplacian_du_loc(i,3-idim) = laplacian_du_loc(i,3-idim) + ((-1.0)**(idim))*(den_dvel_right - den_dvel_left) / dxf

            end do

            ! Crossed terms d_x (sigma nu d_y v_i )
            do idim=1,ndim

               vproml1 = (vel(idim)+vel_left(i,1,idim)+vel_left(i,2,idim)+vel_diag(i,1,1,idim))/4.0
               vproml2 = (vel(idim)+vel_left(i,1,idim)+vel_right(i,2,idim)+vel_diag(i,1,2,idim))/4.0

               vpromr1 = (vel(idim)+vel_right(i,1,idim)+vel_left(i,2,idim)+vel_diag(i,2,1,idim))/4.0
               vpromr2 = (vel(idim)+vel_right(i,1,idim)+vel_right(i,2,idim)+vel_diag(i,2,2,idim))/4.0


               den_dvel_left  = ( (den+den_left(i,1))/2.0 ) * (mu_viscosity+mu_viscosity_left(1))/2.0 * ( ( vproml2 - vproml1 )/dx_diag(i,1,2) )! *(1/2.0)
               den_dvel_right = ( (den+den_right(i,1))/2.0 ) * (mu_viscosity+mu_viscosity_right(1))/2.0 * ( (vpromr2 - vpromr1 ) / dx_diag(i,2,2) )! *(1/2.0)


               dxf = (dx_left(i,1)+dx_right(i,1))/2.0

               laplacian_du_loc(i,3-idim) = laplacian_du_loc(i,3-idim) + ((-1.0)**(idim-1))*(den_dvel_right - den_dvel_left) /  dxf

            end do
         end do

         do i = 1,ngrid
         ! Add viscosity term at time t
            d = max(unew(ind_cell(i),1),smallr)
            u=0.0d0; v=0.0d0; w=0.0d0
            u=unew(ind_cell(i),2)/d
            v=unew(ind_cell(i),3)/d
#if NDIM > 2
            w=unew(ind_cell(i),4)/d
#endif
            e_kin = 0.5d0*d*(u**2 + v**2 + w**2)
            e_nokin = unew(ind_cell(i),ndim+2) - e_kin

            du=0.0d0; dv=0.0d0; dw=0.0d0
            du=unew(ind_cell(i),2)
            dv=unew(ind_cell(i),3)
#if NDIM > 2
            dw=unew(ind_cell(i),4)
#endif
            du = du + laplacian_du_loc(i, 1)*dtnew(ilevel)
            unew(ind_cell(i), 2) = du
            dv = dv + laplacian_du_loc(i, 2)*dtnew(ilevel)
            unew(ind_cell(i), 3) = dv
#if NDIM > 2
            dw = dw + laplacian_du_loc(i, 3)*dtnew(ilevel)
            unew(ind_cell(i), 4) = dw
#endif

            u = du/d
            v = dv/d
            w = dw/d

            e_kin = 0.5d0*d*(u**2 + v**2 + w**2)
            unew(ind_cell(i), ndim + 2) = e_nokin + e_kin
         end do

      end do
      ! End loop over cells
   end do
   ! End loop over sweeps

 111 format('    Entering add_viscosity_terms for level ',i2)

 end subroutine add_viscosity_source_terms




!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine godfine1(ind_grid,ncache,ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  use amr_constants, only:i1min,i1max,j1min,j1max,k1min,k1max, &
                       &  i2min,i2max,j2min,j2max,k2min,k2max, &
                       &  i3min,i3max,j3min,j3max,k3min,k3max
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

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar),save::uloc
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim),save::gloc
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::ploc=0.0d0
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::req_loc=0.0d0
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::peq_loc=0.0d0
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:nvar,1:ndim),save::flux
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:2,1:ndim),save::tmp
  logical ,dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::ok

  integer,dimension(1:nvector),save::ind_cell,ind_buffer

  integer::i,j,ivar,idim,ind_son,iskip
  integer::i0,j0,k0,i1,j1,k1,i2,j2,k2,i3,j3,k3,nx_loc,nb_noneigh
  real(dp)::dx,scale,oneontwotondim,d

  oneontwotondim = 1d0/dble(twotondim)

  ! Mesh spacing in that level
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5d0**ilevel*scale

  gloc=0

  !------------------------------------------
  ! Gather 3^ndim neighboring father cells
  !------------------------------------------
  do i=1,ncache
     ind_cell(i)=father(ind_grid(i))
  end do
  call get3cubefather(ind_cell,nbors_father_cells,ncache,ilevel)

  !---------------------------
  ! Gather 6x6x6 cells stencil
  !---------------------------
  if(levelmin==nlevelmax)then
     call gather_stencil_unigrid(nbors_father_cells,uloc,gloc,req_loc,peq_loc,ok,ncache,ilevel)
  else
     call gather_stencil_amr(nbors_father_cells,uloc,gloc,req_loc,peq_loc,ok,ncache,ilevel)
  end if

  !-----------------------------------------------
  ! Compute flux using second-order Godunov method
  !-----------------------------------------------
  call unsplit(uloc,gloc,ploc,flux,tmp,dx,dx,dx,dtnew(ilevel),ncache)
  !--------------------------------------
  ! Store the fluxes for later use
  !--------------------------------------
  if (MC_tracer) then
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
                 do i=1,ncache
                    d = max(uold(ind_cell(i),1), smallr)
                    ! Copy left flux
                    fluxes(ind_cell(i),(idim-1)*2+1)= flux(i,i3   ,j3   ,k3,   1,idim)&
                         / d
                    ! Copy right flux
                    fluxes(ind_cell(i),(idim-1)*2+2)=-flux(i,i3+i0,j3+j0,k3+k0,1,idim)&
                         / d
                 end do
              end do
           end do
        end do
     end do
  end if

  !------------------------------------------------
  ! Reset flux along direction at refined interface
  !------------------------------------------------
  if(levelmin.lt.nlevelmax)then
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
        if(pressure_fix)then
        do ivar=1,2
           do i=1,ncache
              if(ok(i,i3-i0,j3-j0,k3-k0) .or. ok(i,i3,j3,k3))then
                 tmp (i,i3,j3,k3,ivar,idim)=0.0d0
              end if
           end do
        end do
        end if
     end do
     end do
     end do
  end do
  end if

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
        ! Update conservative variables new state vector
        do ivar=1,nvar
           do i=1,ncache
              unew(ind_cell(i),ivar)=unew(ind_cell(i),ivar)+ &
                   & (flux(i,i3   ,j3   ,k3   ,ivar,idim) &
                   & -flux(i,i3+i0,j3+j0,k3+k0,ivar,idim))
           end do
        end do
        if(pressure_fix)then
        ! Update velocity divergence
        do i=1,ncache
           divu(ind_cell(i))=divu(ind_cell(i))+ &
                & (tmp(i,i3   ,j3   ,k3   ,1,idim) &
                & -tmp(i,i3+i0,j3+j0,k3+k0,1,idim))
        end do
        ! Update internal energy
        do i=1,ncache
           enew(ind_cell(i))=enew(ind_cell(i))+ &
                & (tmp(i,i3   ,j3   ,k3   ,2,idim) &
                & -tmp(i,i3+i0,j3+j0,k3+k0,2,idim))
        end do
        end if
     end do
     end do
     end do
  end do

  !--------------------------------------
  ! Conservative update at level ilevel-1
  !--------------------------------------
  if(levelmin.ne.nlevelmax)then
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
     ! and gather neighbor father cells index
     nb_noneigh=0
     do i=1,ncache
        if (son(nbor(ind_grid(i),2*idim-1))==0) then
           nb_noneigh = nb_noneigh + 1
           ind_buffer(nb_noneigh) = nbor(ind_grid(i),2*idim-1)
           ind_cell(nb_noneigh) = i
        end if
     end do
     ! Conservative update of new state variables
     do ivar=1,nvar
        ! Loop over boundary cells
        do k3=k3min,k3max-k0
        do j3=j3min,j3max-j0
        do i3=i3min,i3max-i0
           do i=1,nb_noneigh
              unew(ind_buffer(i),ivar)=unew(ind_buffer(i),ivar) &
                   & -flux(ind_cell(i),i3,j3,k3,ivar,idim)*oneontwotondim
           end do
        end do
        end do
        end do
     end do
     if(pressure_fix)then
     ! Update velocity divergence
     do k3=k3min,k3max-k0
     do j3=j3min,j3max-j0
     do i3=i3min,i3max-i0
        do i=1,nb_noneigh
           divu(ind_buffer(i))=divu(ind_buffer(i)) &
                & -tmp(ind_cell(i),i3,j3,k3,1,idim)*oneontwotondim
        end do
     end do
     end do
     end do
     ! Update internal energy
     do k3=k3min,k3max-k0
     do j3=j3min,j3max-j0
     do i3=i3min,i3max-i0
        do i=1,nb_noneigh
           enew(ind_buffer(i))=enew(ind_buffer(i)) &
                & -tmp(ind_cell(i),i3,j3,k3,2,idim)*oneontwotondim
        end do
     end do
     end do
     end do
     end if

     !-----------------------
     ! Right flux at boundary
     !-----------------------
     ! Check if grids sits near right boundary
     ! and gather neighbor father cells index
     nb_noneigh=0
     do i=1,ncache
        if (son(nbor(ind_grid(i),2*idim))==0) then
           nb_noneigh = nb_noneigh + 1
           ind_buffer(nb_noneigh) = nbor(ind_grid(i),2*idim)
           ind_cell(nb_noneigh) = i
        end if
     end do
     ! Conservative update of new state variables
     do ivar=1,nvar
        ! Loop over boundary cells
        do k3=k3min+k0,k3max
        do j3=j3min+j0,j3max
        do i3=i3min+i0,i3max
           do i=1,nb_noneigh
              unew(ind_buffer(i),ivar)=unew(ind_buffer(i),ivar) &
                   & +flux(ind_cell(i),i3+i0,j3+j0,k3+k0,ivar,idim)*oneontwotondim
           end do
        end do
        end do
        end do
     end do
     if(pressure_fix)then
     ! Update velocity divergence
     do k3=k3min+k0,k3max
     do j3=j3min+j0,j3max
     do i3=i3min+i0,i3max
        do i=1,nb_noneigh
           divu(ind_buffer(i))=divu(ind_buffer(i)) &
                & +tmp(ind_cell(i),i3+i0,j3+j0,k3+k0,1,idim)*oneontwotondim
        end do
     end do
     end do
     end do
     ! Update internal energy
     do k3=k3min+k0,k3max
     do j3=j3min+j0,j3max
     do i3=i3min+i0,i3max
        do i=1,nb_noneigh
           enew(ind_buffer(i))=enew(ind_buffer(i)) &
                & +tmp(ind_cell(i),i3+i0,j3+j0,k3+k0,2,idim)*oneontwotondim
        end do
     end do
     end do
     end do
     end if

  end do
  ! End loop over dimensions
  endif

end subroutine godfine1
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine gather_stencil_unigrid(nbors_father_cells,uloc,gloc,req_loc,peq_loc,ok,ncache,ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  use amr_constants, only:i1min,i1max,j1min,j1max,k1min,k1max, &
                       &  i2min,i2max,j2min,j2max,k2min,k2max, &
                       &  i3min,i3max,j3min,j3max,k3min,k3max
  implicit none
  integer ,dimension(1:nvector,1:threetondim     ),intent(in)::nbors_father_cells
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar),intent(inout)::uloc
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim),intent(inout)::gloc
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),intent(inout)::req_loc
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),intent(inout)::peq_loc
  logical ,dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),intent(inout)::ok
  integer,intent(in)::ilevel,ncache
  !-------------------------------------------------------------------------
  ! Gather 6x6x6 cells stencil, in the case of a uniform grid, meaning the
  ! neighboring grids on this level are guaranteed to exist.
  !-------------------------------------------------------------------------
  integer,dimension(1:nvector),save::igrid_nbor,ind_cell
  integer::i,ivar,idim,iskip
  integer::i1,j1,k1,i2,j2,k2,i3,j3,k3,ind_son,ind_father

  ! Loop over 3x3x3 neighboring father cells
  do k1=k1min,k1max
  do j1=j1min,j1max
  do i1=i1min,i1max

     ! Get neighbor grid index
     ind_father=1+i1+3*j1+9*k1
     do i=1,ncache
        igrid_nbor(i)=son(nbors_father_cells(i,ind_father))
     end do

     ! Loop over 2x2x2 cells
     do k2=k2min,k2max
     do j2=j2min,j2max
     do i2=i2min,i2max

        ! Get cell index
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
           do i=1,ncache
              uloc(i,i3,j3,k3,ivar)=uold(ind_cell(i),ivar)
           end do
        end do

        ! Gather equilibrium model
        if(strict_equilibrium>0)then
           do idim=1,ndim
              do i=1,ncache
                 req_loc(i,i3,j3,k3)=rho_eq(ind_cell(i))
                 peq_loc(i,i3,j3,k3)=p_eq(ind_cell(i))
              end do
           end do
        end if

        ! Gather gravitational acceleration
        if(poisson)then
           do idim=1,ndim
              do i=1,ncache
                 gloc(i,i3,j3,k3,idim)=f(ind_cell(i),idim)
              end do
           end do
        end if

        ! Gather refinement flag
        do i=1,ncache
           ok(i,i3,j3,k3)=son(ind_cell(i))>0
        end do

     end do
     end do
     end do
     ! End loop over cells

  end do
  end do
  end do
  ! End loop over neighboring grids

end subroutine gather_stencil_unigrid
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine gather_stencil_amr(nbors_father_cells,uloc,gloc,req_loc,peq_loc,ok,ncache,ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  use amr_constants, only:i1min,i1max,j1min,j1max,k1min,k1max, &
                       &  i2min,i2max,j2min,j2max,k2min,k2max, &
                       &  i3min,i3max,j3min,j3max,k3min,k3max
  implicit none
  integer ,dimension(1:nvector,1:threetondim     ),intent(in)::nbors_father_cells
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar),intent(inout)::uloc
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim),intent(inout)::gloc
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),intent(inout)::req_loc
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),intent(inout)::peq_loc
  logical ,dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),intent(inout)::ok
  integer,intent(in)::ilevel,ncache
  !-------------------------------------------------------------------------
  ! Gather 6x6x6 cells stencil, in the general case of AMR. Interpolate from
  ! coarser level missing grid variables.
  !-------------------------------------------------------------------------
  integer ,dimension(1:nvector,0:twondim         ),save::ibuffer_father
  real(dp),dimension(1:nvector,0:twondim  ,1:nvar),save::u1
  real(dp),dimension(1:nvector,1:twotondim,1:nvar),save::u2
  real(dp),dimension(1:nvector,1:twotondim       ),save::req2=0.0d0
  real(dp),dimension(1:nvector,1:twotondim       ),save::peq2=0.0d0

  integer,dimension(1:nvector),save::igrid_nbor,ind_cell,ind_buffer
  integer,dimension(1:nvector),save::ind_exist,ind_nexist
  integer::nexist,nbuffer
  integer::i,j,ivar,idim,iskip
  integer::i1,j1,k1,i2,j2,k2,i3,j3,k3,ind_son,ind_father

  ! Loop over 3x3x3 neighboring father cells
  do k1=k1min,k1max
  do j1=j1min,j1max
  do i1=i1min,i1max

     ! Check if neighboring grid exists
     nbuffer=0
     nexist=0
     ind_father=1+i1+3*j1+9*k1
     do i=1,ncache
        igrid_nbor(i)=son(nbors_father_cells(i,ind_father))
        if(igrid_nbor(i)>0) then
           nexist=nexist+1
           ind_exist(nexist)=i
        else
          nbuffer=nbuffer+1
          ind_nexist(nbuffer)=i
          ind_buffer(nbuffer)=nbors_father_cells(i,ind_father)
        end if
     end do

     ! If not, interpolate hydro variables from parent cells
     if(nbuffer>0)then
        call getnborfather(ind_buffer,ibuffer_father,nbuffer,ilevel)
        do j=0,twondim
           do ivar=1,nvar
              do i=1,nbuffer
                 u1(i,j,ivar)=uold(ibuffer_father(i,j),ivar)
              end do
           end do
        end do
        call interpol_hydro(u1,u2,nbuffer)
     endif

     ! Loop over 2x2x2 cells
     do k2=k2min,k2max
     do j2=j2min,j2max
     do i2=i2min,i2max

        ind_son=1+i2+2*j2+4*k2
        iskip=ncoarse+(ind_son-1)*ngridmax
        do i=1,nexist
           ind_cell(i)=iskip+igrid_nbor(ind_exist(i))
        end do

        i3=1; j3=1; k3=1
        if(ndim>0)i3=1+2*(i1-1)+i2
        if(ndim>1)j3=1+2*(j1-1)+j2
        if(ndim>2)k3=1+2*(k1-1)+k2

        ! Gather hydro variables
        do ivar=1,nvar
           do i=1,nexist
              uloc(ind_exist(i),i3,j3,k3,ivar)=uold(ind_cell(i),ivar)
           end do
           do i=1,nbuffer
              uloc(ind_nexist(i),i3,j3,k3,ivar)=u2(i,ind_son,ivar)
           end do
        end do

        ! Gather equilibrium model
        if(strict_equilibrium>0)then
           do idim=1,ndim
              do i=1,nexist
                 req_loc(ind_exist(i),i3,j3,k3)=rho_eq(ind_cell(i))
                 peq_loc(ind_exist(i),i3,j3,k3)=p_eq(ind_cell(i))
              end do
              ! Use straight injection for buffer cells
              do i=1,nbuffer
                 req_loc(ind_nexist(i),i3,j3,k3)=req2(i,ind_son)
                 peq_loc(ind_nexist(i),i3,j3,k3)=peq2(i,ind_son)
              end do
           end do
        end if

        ! Gather gravitational acceleration
        if(poisson)then
           do idim=1,ndim
              do i=1,nexist
                 gloc(ind_exist(i),i3,j3,k3,idim)=f(ind_cell(i),idim)
              end do
              ! Use straight injection for buffer cells
              do i=1,nbuffer
                 gloc(ind_nexist(i),i3,j3,k3,idim)=f(ibuffer_father(i,0),idim)
              end do
           end do
        end if

        ! Gather refinement flag
        do i=1,nexist
           ok(ind_exist(i),i3,j3,k3)=son(ind_cell(i))>0
        end do
        do i=1,nbuffer
           ok(ind_nexist(i),i3,j3,k3)=.false.
        end do

     end do
     end do
     end do
     ! End loop over cells

  end do
  end do
  end do
  ! End loop over neighboring grids

end subroutine gather_stencil_amr
