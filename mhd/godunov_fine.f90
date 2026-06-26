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
#if USE_FLD==1
  use radiation_parameters,ONLY:eray_min
#endif
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine sets array unew to its initial value uold before calling
  ! the hydro scheme. unew is set to zero in virtual boundaries.
  !--------------------------------------------------------------------------
  integer::i,ivar,ind,icpu,iskip
  real(dp)::d,u,v,w,e,A,B,C
#if NENER>0
  integer::irad
#endif

#if USE_FLD==1
  real(dp)::scale_nH,scale_T2,scale_t,scale_v,scale_d,scale_l
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
#endif
  
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Set unew to uold for myid cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do ivar=1,nvar_all
        do i=1,active(ilevel)%ngrid
           unew(active(ilevel)%igrid(i)+iskip,ivar) = uold(active(ilevel)%igrid(i)+iskip,ivar)
#if USE_FLD==1
           if(ngrp .gt. 0 .and. ivar .gt. firstindex_er .and.ivar .le. firstindex_er+ngrp)then
              unew(active(ilevel)%igrid(i)+iskip,ivar)=max(unew(active(ilevel)%igrid(i)+iskip,ivar),eray_min/(scale_d*scale_v**2))
           end if
#endif
        end do
     end do
     if(pressure_fix)then
        do i=1,active(ilevel)%ngrid
           divu(active(ilevel)%igrid(i)+iskip) = 0
        end do
        do i=1,active(ilevel)%ngrid
           d=max(uold(active(ilevel)%igrid(i)+iskip,1),smallr)
           u=uold(active(ilevel)%igrid(i)+iskip,2)/d
           v=uold(active(ilevel)%igrid(i)+iskip,3)/d
           w=uold(active(ilevel)%igrid(i)+iskip,4)/d
           A=0.5*(uold(active(ilevel)%igrid(i)+iskip,6)+uold(active(ilevel)%igrid(i)+iskip,nvar+1))
           B=0.5*(uold(active(ilevel)%igrid(i)+iskip,7)+uold(active(ilevel)%igrid(i)+iskip,nvar+2))
           C=0.5*(uold(active(ilevel)%igrid(i)+iskip,8)+uold(active(ilevel)%igrid(i)+iskip,nvar+3))
           e=uold(active(ilevel)%igrid(i)+iskip,neul)-0.5*d*(u**2+v**2+w**2)-0.5*(A**2+B**2+C**2)
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
subroutine scale_cosmomag(ind_cell,exp_scale)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ind_cell
  !--------------------------------------------------------------------------
  ! This routine updates magnetic field to scale with cosmic expansion
  !--------------------------------------------------------------------------
  real(dp)::A,B,C,exp_scale,e_mag

  ! Compute old e_mag
  A=0.5*(unew(ind_cell,6)+unew(ind_cell,nvar+1))
  B=0.5*(unew(ind_cell,7)+unew(ind_cell,nvar+2))
  C=0.5*(unew(ind_cell,8)+unew(ind_cell,nvar+3))
  e_mag=0.5*(A**2+B**2+C**2)

  ! Remove from internal energy
  unew(ind_cell,neul) = unew(ind_cell,neul) - e_mag

  ! Rescale B
  unew(ind_cell,6:8) = unew(ind_cell,6:8) * exp_scale
  unew(ind_cell,nvar+1:nvar+3) = unew(ind_cell,nvar+1:nvar+3) * exp_scale

  ! Compute new e_mag
  A=0.5*(unew(ind_cell,6)+unew(ind_cell,nvar+1))
  B=0.5*(unew(ind_cell,7)+unew(ind_cell,nvar+2))
  C=0.5*(unew(ind_cell,8)+unew(ind_cell,nvar+3))
  e_mag=0.5*(A**2+B**2+C**2)

  ! Add back to internal energy
  unew(ind_cell,neul) = unew(ind_cell,neul) + e_mag
end subroutine scale_cosmomag
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine update_cosmomag(ilevel,exp_scale)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine updates magnetic field to scale with cosmic expansion
  !--------------------------------------------------------------------------
  integer::i,ind,iskip,ind_cell,icpu
  real(dp)::exp_scale

  do ind=1,twotondim
    iskip=ncoarse+(ind-1)*ngridmax

    ! Update the active cells
    do i=1,active(ilevel)%ngrid
      ind_cell = active(ilevel)%igrid(i)+iskip
      call scale_cosmomag(ind_cell,exp_scale)
    end do

    ! Do the same for reception cells
    do icpu=1,ncpu
      do i=1,reception(icpu,ilevel)%ngrid
#ifdef LIGHT_MPI_COMM
        ind_cell = reception(icpu,ilevel)%pcomm%igrid(i)+iskip
#else
        ind_cell = reception(icpu,ilevel)%igrid(i)+iskip
#endif
        call scale_cosmomag(ind_cell,exp_scale)
      end do
    end do
  end do
end subroutine update_cosmomag
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine set_uold(ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
#if USE_FLD==1
  use cooling_module,ONLY:kB,mH
  use radiation_parameters,ONLY:Tr_floor
#endif
  implicit none
  integer::ilevel
  !---------------------------------------------------------
  ! This routine sets array uold to its new value unew
  ! after the hydro step.
  !---------------------------------------------------------
  integer::i,ivar,ind,iskip,nx_loc,ind_cell
  real(dp)::scale,d,u,v,w,A,B,C
  real(dp)::e_mag,e_kin,e_cons,e_prim,e_trunc,div,dx
#if NENER>0
  integer::irad
#endif

#if USE_FLD==1
  real(dp)::cv,dd,ee,TP_loc
  real(dp)::scale_nH,scale_T2,scale_t,scale_v,scale_d,scale_l
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
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
     ! L. Romano 14.06.2023 -- Catch advection errors due to smallr
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
        do i=1,active(ilevel)%ngrid
           ind_cell=active(ilevel)%igrid(i)+iskip
           d=max(uold(ind_cell,1),smallr)
           u=uold(ind_cell,2)/d
           v=uold(ind_cell,3)/d
           w=uold(ind_cell,4)/d
           A=0.5*(uold(ind_cell,6)+uold(ind_cell,nvar+1))
           B=0.5*(uold(ind_cell,7)+uold(ind_cell,nvar+2))
           C=0.5*(uold(ind_cell,8)+uold(ind_cell,nvar+3))
           e_kin=0.5*d*(u**2+v**2+w**2)
#if NENER>0
           do irad=1,nener
              e_kin=e_kin+uold(ind_cell,nhydro+irad)
           end do
#endif
           e_mag=0.5*(A**2+B**2+C**2)
           e_cons=uold(ind_cell,neul)-e_kin-e_mag
           e_prim=enew(ind_cell)
           ! Note: here divu=-div.u*dt
           div=abs(divu(ind_cell))*dx/dtnew(ilevel)
           ! Estimate of the local truncation errors
           e_trunc=beta_fix*d*max(div,3.0*hexp*dx)**2
           if(e_cons<e_trunc)then
              uold(ind_cell,neul)=e_prim+e_kin+e_mag
           end if

#if USE_FLD==1
           e_prim = uold(ind_cell,5)-e_kin-e_mag ! uncomment this for radiative shock
           if(energy_fix)then
              e_prim = uold(ind_cell,nvar)
           end if

           ! Compute temperature for perfect gas to prevent crash in the interpolation routine of the EOS
           
           ! Compute gas temperature in cgs
           dd=d *scale_d
           ee=e_prim *scale_d*scale_v**2

           Cv= dd*kB/(mu_gas*mH*(gamma-1.0d0))!(cgs)
           Tp_loc =  ee/(Cv)

           !if(Tp_loc .lt. 0.3d0*Tr_floor .and. ntestDADM.eq.0)then
           if(Tp_loc .lt. 0.3d0*Tr_floor)then
              e_prim = uold(ind_cell,nvar)
              uold(ind_cell,5)=e_prim+e_kin+e_mag
              ee=e_prim *scale_d*scale_v**2
              Cv= dd*kB/(mu_gas*mH*(gamma-1.0d0))!(cgs)
              Tp_loc =  ee/Cv
           endif

           !if(Tp_loc .lt. 0.3d0*Tr_floor.and. ntestDADM.eq.0)then
           if(Tp_loc .lt. 0.3d0*Tr_floor)then
              ! Prevent articifial cooling due to negligible internal energy (wto the total energy)
              call enerint_eos(d,0.3d0*Tr_floor,e_prim)
              !if((.not.(barotrop)).and. (ntestDADM.eq.0))print*,'WARNING TP_LOC',tp_loc,e_trunc,e_prim,e_cons,dd
              uold(ind_cell,5)=e_prim+e_kin+e_mag
              uold(ind_cell,nvar)=e_prim
           end if

           if(energy_fix)then
              uold(ind_cell,5)=e_prim+e_kin+e_mag
              uold(ind_cell,nvar)=e_prim
           end if
#endif
        end do
     end if
  end do

  ! Overwrite state if using induction scheme

  if(ischeme == 1) call velocity_fine(ilevel)

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
#if USE_FLD==1
  use cooling_module!,ONLY:clight
  use radiation_parameters!,ONLY:eray_min,nu_min_hz,nu_max_hz,stellar_photon
  use constants, only:clight
  use units_commons, only:scale_kappa
#endif
  implicit none
  integer::ilevel
  !---------------------------------------------------------
  ! This routine adds the pdV source term to the internal
  ! energy equation and to the non-thermal energy equations.
  !---------------------------------------------------------
  integer::i,ind,iskip,nx_loc,ind_cell1
  integer::ncache,igrid,ngrid,idim,id1,ig1,ih1,id2,ig2,ih2
  real(dp)::scale,dx,dx_loc,d,u,v,w,eold,A,B,C

  integer ,dimension(1:nvector),save::ind_grid,ind_cell
  integer ,dimension(1:nvector,0:twondim),save::igridn
  integer ,dimension(1:nvector,1:ndim),save::ind_left,ind_right
  real(dp),dimension(1:nvector,1:ndim,1:ndim),save::velg,veld
  real(dp),dimension(1:nvector,1:ndim),save::dx_g,dx_d
#if USE_FLD==0
  real(dp),dimension(1:nvector),save::divu_loc
#else
  real(dp),dimension(1:nvector,1:ndim,1:ndim),save::divu_loc
#endif
#if NENER>0
  integer::irad
#endif

#if USE_FLD==1 || USE_M_1==1
  integer::j,ivar
  real(dp)::usquare,emag,ekin,erad_loc,eps,cv,pp_eos

  real(dp),dimension(1:nvector,1:ndim,1:ngrp),save::Erg,Erd
  real(dp) ,dimension(1:nvector,1:ndim,1:ngrp)::gradEr

  integer::k,igroup
  real(dp)::e_mag,e_kin,e_cons,e_prim,e_trunc,div,fact,e_r
  real(dp)::rosseland_ana,radiation_source
  real(dp)::Pgdivu,u_square,d_loc,Tp_loc,Tr_loc,cal_Teg
  real(dp)::kappa_R,gradEr_norm,gradEr_norm2,R,lambda,lambda_fld,chi,PgmErdivu,gradEru
  real(dp) ,dimension(1:ndim,1:ndim,1:ngrp)::Pg
  real(dp) ,dimension(1:ndim       )::u_loc
#if USE_M_1==1
  real(dp), dimension(3,3         ) :: Dedd,Dedd_dE
  real(dp), dimension(3,3,3       ) :: Dedd_dF
  real(dp), dimension(    1:3     ) :: Fr_temp
  real(dp), dimension(3,3,3       ) :: Hr
  real(dp), dimension(3,3,3,1:ngrp) :: Qg
  real(dp), dimension(1:ndim      ) :: nuQrDivu
  real(dp) :: nuQr,nuQl,Qr_nu
#endif
  real(dp) :: nuPrDivu,nuPr,nuPl,Pr_nu

  !  EOS
  real(dp) :: dd,ee,cmp_Cv_eos
  integer  :: ht
#endif

#if USE_FLD==1
  real(dp)::scale_nH,scale_T2,scale_t,scale_v,scale_d,scale_l
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
#endif

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5d0**ilevel
  dx_loc=dx*scale

  velg=0.0; veld=0.0d0

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
#if USE_FLD==1
                 do igroup=1,ngrp
                    Erg(i,idim,igroup) = max(uold(igridn(i,ig1)+ih1,firstindex_er+igroup),eray_min/(scale_d*scale_v**2))
                 end do
#endif                
              else
                 velg(i,idim,1:ndim) = uold(ind_left(i,idim),2:ndim+1)/max(uold(ind_left(i,idim),1),smallr)
                 dx_g(i,idim) = dx_loc*1.5_dp
#if USE_FLD==1
                 do igroup=1,ngrp
                    Erg(i,idim,igroup) = max(uold(ind_left(i,idim),firstindex_er+igroup),eray_min/(scale_d*scale_v**2))
                 end do
#endif
              end if
           enddo
           id2=jjj(idim,2,ind); ig2=iii(idim,2,ind)
           ih2=ncoarse+(id2-1)*ngridmax
           do i=1,ngrid
              if(igridn(i,ig2)>0)then
                 veld(i,idim,1:ndim)= uold(igridn(i,ig2)+ih2,2:ndim+1)/max(uold(igridn(i,ig2)+ih2,1),smallr)
                 dx_d(i,idim)=dx_loc
#if USE_FLD==1
                 do igroup=1,ngrp
                    Erd(i,idim,igroup) = max(uold(igridn(i,ig2)+ih2,firstindex_er+igroup),eray_min/(scale_d*scale_v**2))
                 end do
#endif
              else
                 veld(i,idim,1:ndim)= uold(ind_right(i,idim),2:ndim+1)/max(uold(ind_right(i,idim),1),smallr)
                 dx_d(i,idim)=dx_loc*1.5_dp
#if USE_FLD==1
                 do igroup=1,ngrp
                    Erd(i,idim,igroup) = max(uold(ind_right(i,idim),firstindex_er+igroup),eray_min/(scale_d*scale_v**2))
                 end do
#endif
              end if
           enddo
        end do
        ! End loop over dimensions

#if USE_FLD==0
        ! Compute divu = Trace G
        divu_loc(1:ngrid)=0.0d0
        do i=1,ngrid
           do idim=1,ndim
              divu_loc(i) = divu_loc(i) + (veld(i,idim,idim)-velg(i,idim,idim)) &
                   &                    / (dx_g(i,idim)     +dx_d(i,idim))
           enddo
        end do
#else
        divu_loc(1:ngrid,1:ndim,1:ndim)=0.0d0
        do i=1,ngrid
           do j=1,ndim
              do idim=1,ndim
                 divu_loc(i,j,idim) = (veld(i,j,idim)-velg(i,j,idim)) &
                      &                    / (dx_g(i,j)     +dx_d(i,j))
              end do
              do igroup=1,ngrp
                 gradEr(i,j,igroup) = (Erd(i,j,igroup)-Erg(i,j,igroup))/(dx_g(i,j)+dx_d(i,j))
              enddo
           enddo
        end do
#endif

        ! Update thermal internal energy
        if(pressure_fix)then
           do i=1,ngrid
              ! Compute old thermal energy
              d=max(uold(ind_cell(i),1),smallr)
              u=0; v=0; w=0
              if(ndim>0)u=uold(ind_cell(i),2)/d
              if(ndim>1)v=uold(ind_cell(i),3)/d
              if(ndim>2)w=uold(ind_cell(i),4)/d
              A=0.5*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
              B=0.5*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
              C=0.5*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))
              eold=uold(ind_cell(i),neul)-0.5d0*d*(u**2+v**2+w**2)-0.5*(A**2+B**2+C**2)
#if NENER>0
              do irad=1,nener
                 eold=eold-uold(ind_cell(i),nhydro+irad)
              end do
#endif
              ! Add -pdV term
#if USE_FLD==0
              enew(ind_cell(i))=enew(ind_cell(i)) &
                   & -(gamma-1.0d0)*eold*divu_loc(i)*dtnew(ilevel)
#else
              do idim=1,ndim
                 enew(ind_cell(i))=enew(ind_cell(i)) &
                      & -(gamma-1.0d0)*eold*divu_loc(i,idim,idim)*dtnew(ilevel)
              end do
#endif
           end do
        end if

#if USE_FLD==0
#if NENER>0
        do irad=1,nener
           do i=1,ngrid
              ! Add -pdV term
              unew(ind_cell(i),nhydro+irad)=unew(ind_cell(i),nhydro+irad) &
                & -(gamma_rad(irad)-1.0d0)*uold(ind_cell(i),nhydro+irad)*divu_loc(i)*dtnew(ilevel)
           end do
        end do
#endif
#else
!***********
#if USE_FLD==1 || USE_M_1==1
        do i=1,ngrid
           d_loc = uold(ind_cell(i),1)*scale_d
           u_loc(1:ndim) = uold(ind_cell(i),2:ndim+1)/uold(ind_cell(i),1)
           
           usquare=0.0
           do idim=1,ndim
              usquare=usquare+(uold(ind_cell(i),idim+1)/uold(ind_cell(i),1))**2
           end do
           
           ! Compute total magnetic energy
           emag = 0.0d0
           do ivar=1,3
              emag = emag + 0.125d0*(uold(ind_cell(i),5+ivar) &
                   &  +uold(ind_cell(i),nvar+ivar))**2
           end do
           erad_loc=0.0D0
#if NENER>0
           do igroup=1,nener
              erad_loc = erad_loc + uold(ind_cell(i),8+igroup)
           end do
#endif
           d     = uold(ind_cell(i),1)
           ekin  = d*usquare/2.0
           ! Compute gas temperature in cgs
           eps   = uold(ind_cell(i),5)-ekin-emag-erad_loc
           if(energy_fix)eps   = uold(ind_cell(i),nvar)
  
           call temperature_eos(d,eps,Tp_loc,ht)

           ! Compute radiative pressure in all groups
           do igroup=1,ngrp
#if USE_FLD==1
              ! Compute radiative pressure
              Tr_loc = cal_Teg(uold(ind_cell(i),firstindex_er+igroup)*scale_d*scale_v**2,igroup)
              kappa_R = rosseland_ana(d_loc,Tp_loc,Tr_loc,igroup,in_sink(ind_cell(i)))/scale_kappa
              gradEr_norm2 = 0.0d0
              do idim=1,ndim
                 gradEr_norm2 = gradEr_norm2 + gradEr(i,idim,igroup)**2
              end do
              gradEr_norm  = (gradEr_norm2)**0.5
              R =  max(1.d-10,gradEr_norm/(max(uold(ind_cell(i),firstindex_er+igroup),eray_min/(scale_d*scale_v**2))*kappa_R))
              lambda = lambda_fld(R)
              chi = lambda + (lambda*R)**2
              do j=1,ndim
                 do k=1,ndim
                    Pg(j,k,igroup)=0.0d0
                    if(j .eq. k) Pg(j,k,igroup) = (1.0d0-chi)/2.0d0
                    if(R .gt. 1.d-8)Pg(j,k,igroup) = Pg(j,k,igroup) + (3.0d0*chi-1.0d0)/2.0d0*gradEr(i,j,igroup)*gradEr(i,k,igroup)/gradEr_norm2
                 enddo
              enddo
              Pg(1:ndim,1:ndim,igroup)=Pg(1:ndim,1:ndim,igroup)*uold(ind_cell(i),firstindex_er+igroup)
#endif
#if USE_M_1==1
              Fr_temp=zero
              do j=1,ndim
                 Fr_temp(j)=uold(ind_cell(i),firstindex_er+j*ngrp+igroup)
              enddo
              ! NOTE: due to normalisation at the initialisation, we have to multiply F by (scale_v/clight)
              !       due to the definition of F in cal_Dedd, we also have to divide F by an additional clight
!!$              call cal_Dedd(uold(ind_cell(i),firstindex_er+igroup),Fr_temp*scale_v/(clight*clight),Dedd,Dedd_dE,Dedd_dF)
              call cal_Dedd(uold(ind_cell(i),firstindex_er+igroup),Fr_temp*scale_v/clight,Dedd,Dedd_dE,Dedd_dF)
              Pg(1:ndim,1:ndim,igroup)=Dedd(1:ndim,1:ndim)*uold(ind_cell(i),firstindex_er+igroup)

              !WARNING WITH NORMALIZATION !!!
!              print*,'normalization to be assessed cal_hr called in godunov_fine and qr_nu...'
!              read(*,*)
              call cal_Hr(uold(ind_cell(i),firstindex_er+igroup),Fr_temp*scale_v,Hr)
              Qg(1:ndim,1:ndim,1:ndim,igroup)=Hr(1:ndim,1:ndim,1:ndim)*uold(ind_cell(i),firstindex_er+igroup)*clight/scale_v
#endif
           enddo
              
           do igroup=1,ngrp
#if USE_FLD==1
              ! Compute radiative pressure
              Tr_loc = cal_Teg(uold(ind_cell(i),firstindex_er+igroup)*scale_d*scale_v**2,igroup)
              kappa_R = rosseland_ana(d_loc,Tp_loc,Tr_loc,igroup,in_sink(ind_cell(i)))/scale_kappa
              gradEr_norm2 = 0.0d0
              do idim=1,ndim
                 gradEr_norm2 = gradEr_norm2 + gradEr(i,idim,igroup)**2
              end do
              gradEr_norm  = (gradEr_norm2)**0.5
              R =   max(1.d-10,gradEr_norm/(max(uold(ind_cell(i),firstindex_er+igroup),eray_min/(scale_d*scale_v**2))*kappa_R))
              lambda = lambda_fld(R)
#endif

              Pgdivu     = 0.0d0
              PgmErdivu  = 0.0d0
              gradEru    = 0.0d0
              nuPrDivu   = 0.0d0 ! multig: this is the Doppler shift term
#if USE_M_1==1
              nuQrDivu   = zero ! multig: this is the Doppler shift term
#endif

              do j=1,ndim
                 do k=1,ndim
                 
                    ! compure Pr:divU term
                    Pgdivu    = Pgdivu    + Pg(j,k,igroup)*divu_loc(i,j,k)
                    PgmErdivu = PgmErdivu + Pg(j,k,igroup)*divu_loc(i,j,k)
                    
                    ! compute -[nu P]*Div(u) term
                    if(divu_loc(i,j,k) > 0.0d0)then

                       if(igroup == ngrp)then
                          nuPr = 0.0d0
                       else
                          nuPr = nu_max_hz(igroup)*Pr_nu(Pg(j,k,igroup+1),uold(ind_cell(i),firstindex_er+igroup+1),igroup+1,nu_max_hz(igroup),scale_d,scale_v)
                       endif
                       if(igroup == 1 .or. (stellar_photon .and. igroup==2))then
                          nuPl = 0.0d0
                       else
                          nuPl = nu_min_hz(igroup)*Pr_nu(Pg(j,k,igroup),uold(ind_cell(i),firstindex_er+igroup),igroup,nu_min_hz(igroup),scale_d,scale_v)
                       endif

                    else

                       if(igroup == ngrp)then
                          nuPr = 0.0d0
                       else
                          nuPr = nu_max_hz(igroup)*Pr_nu(Pg(j,k,igroup),uold(ind_cell(i),firstindex_er+igroup),igroup,nu_max_hz(igroup),scale_d,scale_v)
                       endif
                       if(igroup == 1 .or. (stellar_photon .and. igroup==2))then
                          nuPl = 0.0d0
                       else
                          nuPl = nu_min_hz(igroup)*Pr_nu(Pg(j,k,igroup-1),uold(ind_cell(i),firstindex_er+igroup-1),igroup-1,nu_min_hz(igroup),scale_d,scale_v)
                       endif

                    endif
                    
                    nuPrDivu = nuPrDivu + (nuPr - nuPl)*divu_loc(i,j,k)

#if USE_M_1==1
                    do l=1,ndim

                       ! compute [nu Q]*Div(u) term
                       if(divu_loc(i,j,k) > zero)then

                          if(igroup == ngrp)then
                             nuQr = zero
                          else
                             nuQr = nu_max_hz(igroup)*Qr_nu(Qg(j,k,l,igroup+1),uold(ind_cell(i),firstindex_er+igroup+1),igroup+1,nu_max_hz(igroup),scale_d,scale_v)
                          endif
                          if(igroup == 1 .or. (stellar_photon .and. igroup==2))then
                             nuQl = zero
                          else
                             nuQl = nu_min_hz(igroup)*Qr_nu(Qg(j,k,l,igroup),uold(ind_cell(i),firstindex_er+igroup),igroup,nu_min_hz(igroup),scale_d,scale_v)
                          endif

                       else

                          if(igroup == ngrp)then
                             nuQr = zero
                          else
                             nuQr = nu_max_hz(igroup)*Qr_nu(Qg(j,k,l,igroup),uold(ind_cell(i),firstindex_er+igroup),igroup,nu_max_hz(igroup),scale_d,scale_v)
                          endif
                          if(igroup == 1 .or. (stellar_photon .and. igroup==2))then
                             nuQl = zero
                          else
                             nuQl = nu_min_hz(igroup)*Qr_nu(Qg(j,k,l,igroup-1),uold(ind_cell(i),firstindex_er+igroup-1),igroup-1,nu_min_hz(igroup),scale_d,scale_v)
                          endif

                       endif

                       nuQrDivu(l) = nuQrDivu(l) + (nuQr - nuQl)*divu_loc(i,j,k)

                    enddo
#endif

                 enddo
                 PgmErdivu = PgmErdivu - uold(ind_cell(i),firstindex_er+igroup)/3.0d0*divu_loc(i,j,j)
                 gradEru = gradEru + gradEr(i,j,igroup)*u_loc(j)
              enddo

#if USE_FLD==1
              if(stellar_photon .and. igroup==1)nuPrdivu=0.0
              unew(ind_cell(i),firstindex_er+igroup) = unew(ind_cell(i),firstindex_er+igroup)  - (Pgdivu - nuPrDivu)*dtnew(ilevel)
              unew(ind_cell(i),firstindex_er+igroup)=max(unew(ind_cell(i),firstindex_er+igroup),eray_min/(scale_d*scale_v**2))

              unew(ind_cell(i),5) = unew(ind_cell(i),5) - dtnew(ilevel)*( &
                   & (lambda-1.0d0/3.0d0)*gradEru + &
                   & PgmErdivu - nuPrDivu)

              do idim=1,ndim
                 unew(ind_cell(i),idim+1) = unew(ind_cell(i),1+idim) - dtnew(ilevel)*( &
                      & (lambda -1.0d0/3.0d0)*gradEr(i,idim,igroup))
              enddo    
#endif
#if USE_M_1==1
              unew(ind_cell(i),5) = unew(ind_cell(i),5) - dtnew(ilevel)*(Pgdivu-nuPrDivu)
              unew(ind_cell(i),firstindex_er+igroup) = unew(ind_cell(i),firstindex_er+igroup) - dtnew(ilevel)*(Pgdivu-nuPrDivu)

              do idim=1,ndim
                 print*,'Warning in godunov_fine with M1, update/check Tr_loc for rosseland_ana'
                 stop

                 unew(ind_cell(i),1+idim) = unew(ind_cell(i),1+idim) - dtnew(ilevel)*rosseland_ana(d_loc,Tp_loc,Tr_loc,igroup,in_sink(ind_cell(i)))/scale_kappa*uold(ind_cell(i),firstindex_er+idim*ngrp+igroup)*scale_v/clight
                 unew(ind_cell(i),5) = unew(ind_cell(i),5) - dtnew(ilevel)*rosseland_ana(d_loc,Tp_loc,Tr_loc,igroup,in_sink(ind_cell(i)))/scale_kappa*uold(ind_cell(i),firstindex_er+idim*ngrp+igroup)*scale_v/clight*uold(ind_cell(i),1+idim)

                 do j = 1,ndim
!                    unew(ind_cell(i),firstindex_er+idim*ngrp+igroup) = unew(ind_cell(i),firstindex_er+idim*ngrp+igroup) - dtnew(ilevel)*Fr_temp(j)*divu_loc(i,idim,j)
                    unew(ind_cell(i),firstindex_er+idim*ngrp+igroup) = unew(ind_cell(i),firstindex_er+idim*ngrp+igroup) - dtnew(ilevel)*(Fr_temp(j)*divu_loc(i,j,idim)-nuQrDivu(j))
                 enddo
              enddo
#endif
           enddo !end loop over rad groups
        end do ! end loop over ngrid
#endif

#if NENER>NGRP
        do irad=1,nent
           do i=1,ngrid
              ! Add -pdV term
              do idim=1,ndim
                 unew(ind_cell(i),8+irad)=unew(ind_cell(i),8+irad) &
                      & -(gamma_rad(irad)-1.0d0)*uold(ind_cell(i),8+irad)*divu_loc(i,idim,idim)*dtnew(ilevel)
              end do
           end do
        end do
#endif

        !update internal energy in unew(nvar)
        do i=1,ngrid

           usquare=0.0
           do idim=1,ndim
              usquare=usquare+(uold(ind_cell(i),idim+1)/uold(ind_cell(i),1))**2
           end do
           
           ! Compute total magnetic energy
           emag = 0.0d0
           do ivar=1,3
              emag = emag + 0.125d0*(uold(ind_cell(i),5+ivar) &
                   &  +uold(ind_cell(i),nvar+ivar))**2
           end do
           erad_loc=0.0D0
#if NENER>0
           do igroup=1,nener
              erad_loc = erad_loc + uold(ind_cell(i),8+igroup)
           end do
#endif

           d     = uold(ind_cell(i),1)
           ekin  = d*usquare/2.0
           ! Compute gas pressure in cgs
           eps   = uold(ind_cell(i),5)-ekin-emag-erad_loc
           if(energy_fix)eps   = uold(ind_cell(i),nvar)
           
           call pressure_eos(d,eps,pp_eos)
           do idim=1,ndim
              unew(ind_cell(i),nvar) = unew(ind_cell(i),nvar) &
                   & + (-pp_eos*divu_loc(i,idim,idim))*dtnew(ilevel)
           end do
        end do
!***********
#endif
        
     enddo
     ! End loop over cells
  end do
  ! End loop over grids

  return

  print*,'**********stop******'
  stop
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
  integer ,dimension(1:nvector,0:twondim         ),save::ibuffer_father
  integer ,dimension(1:nvector,0:twondim         ),save::ind1
  real(dp),dimension(1:nvector,0:twondim  ,1:nvar+3),save::u1
  real(dp),dimension(1:nvector,1:twotondim,1:nvar+3),save::u2

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3),save::uloc
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim),save::gloc=0.0d0
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:nvar,1:ndim),save::flux
  real(dp),dimension(1:nvector,1:3,1:3,1:3),save::emfx=0.0d0,emfy=0.0d0,emfz=0.0d0
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:2,1:ndim),save::tmp
  logical ,dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::ok

  integer,dimension(1:nvector),save::igrid_nbor,ind_cell,ind_buffer,ind_exist,ind_nexist

  integer::ind_buffer1,ind_buffer2,ind_buffer3
  integer::ind_father1,ind_father2,ind_father3
  integer::i,j,ivar,idim,ind_son,ind_father,iskip,nbuffer
  integer::i0,j0,k0,i1,j1,k1,i2,j2,k2,i3,j3,k3,nx_loc,nb_noneigh,nexist
  real(dp)::dflux_x,dflux_y,dflux_z
  real(dp)::dx,scale,oneontwotondim,d
  real(dp)::dflux,weight

  oneontwotondim = 1d0/dble(twotondim)

  ! Mesh spacing in that level
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5d0**ilevel*scale

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

     ! If not, interpolate variables from parent cells
     if(nbuffer>0)then
        call getnborfather(ind_buffer,ibuffer_father,nbuffer,ilevel)
        do j=0,twondim
           do ivar=1,nvar+3
              do i=1,nbuffer
                 u1(i,j,ivar)=uold(ibuffer_father(i,j),ivar)
              end do
           end do
           do i=1,nbuffer
              ind1(i,j)=son(ibuffer_father(i,j))
           end do
        end do
        call interpol_hydro(u1,ind1,u2,nbuffer)
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
        do ivar=1,nvar+3
           do i=1,nexist
              uloc(ind_exist(i),i3,j3,k3,ivar)=uold(ind_cell(i),ivar)
           end do
           do i=1,nbuffer
              uloc(ind_nexist(i),i3,j3,k3,ivar)=u2(i,ind_son,ivar)
           end do
        end do

        ! Gather gravitational acceleration
        if(poisson)then
           do idim=1,ndim
              do i=1,nexist
                 gloc(ind_exist(i),i3,j3,k3,idim)=f(ind_cell(i),idim)
              end do
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

  !-----------------------------------------------
  ! Compute flux using second-order Godunov method
  !-----------------------------------------------
  call mag_unsplit(uloc,gloc,flux,emfx,emfy,emfz,tmp,dx,dx,dx,dtnew(ilevel),ncache)
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

  if(ischeme.eq.1)then
  !---------------------------------
  ! Reset all Euler variables fluxes
  !---------------------------------
  do idim=1,ndim
  do ivar=1,neul
     do k3=k3min,k3max+1
     do j3=j3min,j3max+1
     do i3=i3min,i3max+1
        do i=1,ncache
           flux(i,i3,j3,k3,ivar,idim)=0.0d0
        end do
     end do
     end do
     end do
  end do
  end do

  endif

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
  !---------------------------------------------------------------
  ! Reset Euler fluxes for Bx
  !---------------------------------------------------------------
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     do k3=k3min,k3max+k0
     do j3=j3min,j3max+j0
     do i3=i3min,i3max+i0
        do i=1,ncache
           flux(i,i3,j3,k3,6,idim)=0.0d0
        end do
     end do
     end do
     end do
  end do
#if NDIM>1
  !---------------------------------------------------------------
  ! Reset electromotive force along direction z at refined edges
  !---------------------------------------------------------------
  do k3=k3min,k3max
  do j3=1,3
  do i3=1,3
     do i=1,ncache
        if(ok(i,i3  ,j3  ,k3) .or. ok(i,i3  ,j3-1,k3) .or.  &
         & ok(i,i3-1,j3  ,k3) .or. ok(i,i3-1,j3-1,k3))then
           emfz(i,i3,j3,k3)=0.0d0
#if NDIM==2
           emfz(i,i3,j3,k3+1)=0.0d0
#endif
        end if
     end do
  end do
  end do
  end do
  !---------------------------------------------------------------
  ! Reset Euler fluxes for By
  !---------------------------------------------------------------
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     do k3=k3min,k3max+k0
     do j3=j3min,j3max+j0
     do i3=i3min,i3max+i0
        do i=1,ncache
           flux(i,i3,j3,k3,7,idim)=0.0d0
        end do
     end do
     end do
     end do
  end do
#endif
#if NDIM>2
  !---------------------------------------------------------------
  ! Reset electromotive force along direction y at refined edges
  !---------------------------------------------------------------
  do k3=1,3
  do j3=1,2
  do i3=1,3
     do i=1,ncache
        if(ok(i,i3  ,j3,k3  ) .or. ok(i,i3  ,j3,k3-1) .or.  &
         & ok(i,i3-1,j3,k3  ) .or. ok(i,i3-1,j3,k3-1))then
           emfy(i,i3,j3,k3)=0.0d0
        end if
     end do
  end do
  end do
  end do
  !---------------------------------------------------------------
  ! Reset electromotive force along direction x at refined edges
  !---------------------------------------------------------------
  do k3=1,3
  do j3=1,3
  do i3=1,2
     do i=1,ncache
        if(ok(i,i3,j3  ,k3  ) .or. ok(i,i3,j3  ,k3-1) .or.  &
         & ok(i,i3,j3-1,k3  ) .or. ok(i,i3,j3-1,k3-1))then
           emfx(i,i3,j3,k3)=0.0d0
        end if
     end do
  end do
  end do
  end do
  !---------------------------------------------------------------
  ! Reset Euler fluxes for Bz
  !---------------------------------------------------------------
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     do k3=k3min,k3max+k0
     do j3=j3min,j3max+j0
     do i3=i3min,i3max+i0
        do i=1,ncache
           flux(i,i3,j3,k3,8,idim)=0.0d0
        end do
     end do
     end do
     end do
  end do
#endif

  !-----------------------------------------------------
  ! Conservative update at level ilevel for Euler system
  !-----------------------------------------------------
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
        do ivar=1,3
           do i=1,ncache
              unew(ind_cell(i),nvar+ivar)=unew(ind_cell(i),nvar+ivar)+ &
                   & (flux(i,i3   ,j3   ,k3   ,neul+ivar,idim) &
                   & -flux(i,i3+i0,j3+j0,k3+k0,neul+ivar,idim))
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

  !---------------------------------------------------------
  ! Conservative update at level ilevel for induction system
  !---------------------------------------------------------
  do k3=k3min,k3max
  do j3=j3min,j3max
  do i3=1,2
     ind_son=i3+2*(j3-1)+4*(k3-1)
     iskip=ncoarse+(ind_son-1)*ngridmax
     do i=1,ncache
        ind_cell(i)=iskip+ind_grid(i)
     end do
     ! Update Bx using constraint transport
     do i=1,ncache
        dflux_x=( emfy(i,i3,j3,k3)-emfy(i,i3,j3,k3+1) ) &
          &    -( emfz(i,i3,j3,k3)-emfz(i,i3,j3+1,k3) )
        unew(ind_cell(i),neul+1)=unew(ind_cell(i),neul+1)+dflux_x
        dflux_x=( emfy(i,i3+1,j3,k3)-emfy(i,i3+1,j3,k3+1) ) &
          &    -( emfz(i,i3+1,j3,k3)-emfz(i,i3+1,j3+1,k3) )
        unew(ind_cell(i),nvar+1)=unew(ind_cell(i),nvar+1)+dflux_x
     end do
  end do
  end do
  end do
#if NDIM>1
  do k3=k3min,k3max
  do j3=1,2
  do i3=1,2
     ind_son=i3+2*(j3-1)+4*(k3-1)
     iskip=ncoarse+(ind_son-1)*ngridmax
     do i=1,ncache
        ind_cell(i)=iskip+ind_grid(i)
     end do
     ! Update By using constraint transport
     do i=1,ncache
        dflux_y=( emfz(i,i3,j3,k3)-emfz(i,i3+1,j3,k3) ) &
          &    -( emfx(i,i3,j3,k3)-emfx(i,i3,j3,k3+1) )
        unew(ind_cell(i),neul+2)=unew(ind_cell(i),neul+2)+dflux_y
        dflux_y=( emfz(i,i3,j3+1,k3)-emfz(i,i3+1,j3+1,k3) ) &
          &    -( emfx(i,i3,j3+1,k3)-emfx(i,i3,j3+1,k3+1) )
        unew(ind_cell(i),nvar+2)=unew(ind_cell(i),nvar+2)+dflux_y
     end do
  end do
  end do
  end do
#endif
#if NDIM>2
  do k3=1,2
  do j3=1,2
  do i3=1,2
     ind_son=i3+2*(j3-1)+4*(k3-1)
     iskip=ncoarse+(ind_son-1)*ngridmax
     do i=1,ncache
        ind_cell(i)=iskip+ind_grid(i)
     end do
     ! Update Bz using constraint transport
     do i=1,ncache
        dflux_z=( emfx(i,i3,j3,k3)-emfx(i,i3,j3+1,k3) ) &
          &    -( emfy(i,i3,j3,k3)-emfy(i,i3+1,j3,k3) )
        unew(ind_cell(i),neul+3)=unew(ind_cell(i),neul+3)+dflux_z
        dflux_z=( emfx(i,i3,j3,k3+1)-emfx(i,i3,j3+1,k3+1) ) &
          &    -( emfy(i,i3,j3,k3+1)-emfy(i,i3+1,j3,k3+1) )
        unew(ind_cell(i),nvar+3)=unew(ind_cell(i),nvar+3)+dflux_z
     end do
  end do
  end do
  end do
#endif

  if(ilevel>levelmin)then

  !-----------------------------------------------------------
  ! Conservative update at level ilevel-1 for the Euler system
  !-----------------------------------------------------------
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
     do ivar=1,3
        ! Loop over boundary cells
        do k3=k3min,k3max-k0
        do j3=j3min,j3max-j0
        do i3=i3min,i3max-i0
           do i=1,nb_noneigh
              unew(ind_buffer(i),nvar+ivar)=unew(ind_buffer(i),nvar+ivar) &
                   & -flux(ind_cell(i),i3,j3,k3,neul+ivar,idim)*oneontwotondim
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
     do ivar=1,3
        ! Loop over boundary cells
        do k3=k3min+k0,k3max
        do j3=j3min+j0,j3max
        do i3=i3min+i0,i3max
           do i=1,nb_noneigh
              unew(ind_buffer(i),nvar+ivar)=unew(ind_buffer(i),nvar+ivar) &
                   & +flux(ind_cell(i),i3+i0,j3+j0,k3+k0,neul+ivar,idim)*oneontwotondim
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

#if NDIM>1
  !---------------------------------------------------------------
  ! Conservative update at level ilevel-1 for the induction system
  !---------------------------------------------------------------
  i1=1; j1=0; k1=0
  if(ndim>1)j1=1
  if(ndim>2)k1=1

  !--------------------------------------
  ! Deal with 4 EMFz edges
  !--------------------------------------

  ! Update coarse Bx and By using fine EMFz on X=0 and Y=0 grid edge
  ind_father1=1+(i1  )+3*(j1-1)+9*(k1  )
  ind_father2=1+(i1-1)+3*(j1-1)+9*(k1  )
  ind_father3=1+(i1-1)+3*(j1  )+9*(k1  )
  do i=1,ncache
     ind_buffer1=nbors_father_cells(i,ind_father1)
     ind_buffer2=nbors_father_cells(i,ind_father2)
     ind_buffer3=nbors_father_cells(i,ind_father3)
     weight=1.0
     if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
     if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
     dflux=(emfz(i,1,1,1)+emfz(i,1,1,2))*0.25*weight
     unew(ind_buffer1,1+neul)=unew(ind_buffer1,1+neul)+dflux
     unew(ind_buffer2,1+nvar)=unew(ind_buffer2,1+nvar)+dflux
     unew(ind_buffer2,2+nvar)=unew(ind_buffer2,2+nvar)-dflux
     unew(ind_buffer3,2+neul)=unew(ind_buffer3,2+neul)-dflux
     if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
        unew(ind_buffer3,1+nvar)=unew(ind_buffer3,1+nvar)-dflux*0.5
        unew(ind_buffer1,2+nvar)=unew(ind_buffer1,2+nvar)+dflux*0.5
     endif
  end do

  ! Update coarse Bx and By using fine EMFz on X=0 and Y=1 grid edge
  ind_father1=1+(i1-1)+3*(j1  )+9*(k1  )
  ind_father2=1+(i1-1)+3*(j1+1)+9*(k1  )
  ind_father3=1+(i1  )+3*(j1+1)+9*(k1  )
  do i=1,ncache
     ind_buffer1=nbors_father_cells(i,ind_father1)
     ind_buffer2=nbors_father_cells(i,ind_father2)
     ind_buffer3=nbors_father_cells(i,ind_father3)
     weight=1.0
     if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
     if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
     dflux=(emfz(i,1,3,1)+emfz(i,1,3,2))*0.25*weight
     unew(ind_buffer1,2+nvar)=unew(ind_buffer1,2+nvar)-dflux
     unew(ind_buffer2,2+neul)=unew(ind_buffer2,2+neul)-dflux
     unew(ind_buffer2,1+nvar)=unew(ind_buffer2,1+nvar)-dflux
     unew(ind_buffer3,1+neul)=unew(ind_buffer3,1+neul)-dflux
     if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
        unew(ind_buffer3,2+neul)=unew(ind_buffer3,2+neul)+dflux*0.5
        unew(ind_buffer1,1+nvar)=unew(ind_buffer1,1+nvar)+dflux*0.5
     endif
  end do

  ! Update coarse Bx and By using fine EMFz on X=1 and Y=1 grid edge
  ind_father1=1+(i1  )+3*(j1+1)+9*(k1  )
  ind_father2=1+(i1+1)+3*(j1+1)+9*(k1  )
  ind_father3=1+(i1+1)+3*(j1  )+9*(k1  )
  do i=1,ncache
     ind_buffer1=nbors_father_cells(i,ind_father1)
     ind_buffer2=nbors_father_cells(i,ind_father2)
     ind_buffer3=nbors_father_cells(i,ind_father3)
     weight=1.0
     if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
     if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
     dflux=(emfz(i,3,3,1)+emfz(i,3,3,2))*0.25*weight
     unew(ind_buffer1,1+nvar)=unew(ind_buffer1,1+nvar)-dflux
     unew(ind_buffer2,1+neul)=unew(ind_buffer2,1+neul)-dflux
     unew(ind_buffer2,2+neul)=unew(ind_buffer2,2+neul)+dflux
     unew(ind_buffer3,2+nvar)=unew(ind_buffer3,2+nvar)+dflux
     if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
        unew(ind_buffer3,1+neul)=unew(ind_buffer3,1+neul)+dflux*0.5
        unew(ind_buffer1,2+neul)=unew(ind_buffer1,2+neul)-dflux*0.5
     endif
  end do

  ! Update coarse Bx and By using fine EMFz on X=1 and Y=0 grid edge
  ind_father1=1+(i1+1)+3*(j1  )+9*(k1  )
  ind_father2=1+(i1+1)+3*(j1-1)+9*(k1  )
  ind_father3=1+(i1  )+3*(j1-1)+9*(k1  )
  do i=1,ncache
     ind_buffer1=nbors_father_cells(i,ind_father1)
     ind_buffer2=nbors_father_cells(i,ind_father2)
     ind_buffer3=nbors_father_cells(i,ind_father3)
     weight=1.0
     if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
     if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
     dflux=(emfz(i,3,1,1)+emfz(i,3,1,2))*0.25*weight
     unew(ind_buffer1,2+neul)=unew(ind_buffer1,2+neul)+dflux
     unew(ind_buffer2,2+nvar)=unew(ind_buffer2,2+nvar)+dflux
     unew(ind_buffer2,1+neul)=unew(ind_buffer2,1+neul)+dflux
     unew(ind_buffer3,1+nvar)=unew(ind_buffer3,1+nvar)+dflux
     if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
        unew(ind_buffer3,2+nvar)=unew(ind_buffer3,2+nvar)-dflux*0.5
        unew(ind_buffer1,1+neul)=unew(ind_buffer1,1+neul)-dflux*0.5
     endif
  end do
#endif
  !--------------------------------------
  ! Deal with 4 EMFx edges
  !--------------------------------------
#if NDIM>2
  ! Update coarse By and Bz using fine EMFx on Y=0 and Z=0 grid edge
  ind_father1=1+(i1  )+3*(j1  )+9*(k1-1)
  ind_father2=1+(i1  )+3*(j1-1)+9*(k1-1)
  ind_father3=1+(i1  )+3*(j1-1)+9*(k1  )
  do i=1,ncache
     ind_buffer1=nbors_father_cells(i,ind_father1)
     ind_buffer2=nbors_father_cells(i,ind_father2)
     ind_buffer3=nbors_father_cells(i,ind_father3)
     weight=1.0
     if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
     if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
     dflux=(emfx(i,1,1,1)+emfx(i,2,1,1))*0.25*weight
     unew(ind_buffer1,2+neul)=unew(ind_buffer1,2+neul)+dflux
     unew(ind_buffer2,2+nvar)=unew(ind_buffer2,2+nvar)+dflux
     unew(ind_buffer2,3+nvar)=unew(ind_buffer2,3+nvar)-dflux
     unew(ind_buffer3,3+neul)=unew(ind_buffer3,3+neul)-dflux
     if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
        unew(ind_buffer1,3+nvar)=unew(ind_buffer1,3+nvar)+dflux*0.5
        unew(ind_buffer3,2+nvar)=unew(ind_buffer3,2+nvar)-dflux*0.5
     endif
  end do

  ! Update coarse By and Bz using fine EMFx on Y=0 and Z=1 grid edge
  ind_father1=1+(i1  )+3*(j1-1)+9*(k1  )
  ind_father2=1+(i1  )+3*(j1-1)+9*(k1+1)
  ind_father3=1+(i1  )+3*(j1  )+9*(k1+1)
  do i=1,ncache
     ind_buffer1=nbors_father_cells(i,ind_father1)
     ind_buffer2=nbors_father_cells(i,ind_father2)
     ind_buffer3=nbors_father_cells(i,ind_father3)
     weight=1.0
     if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
     if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
     dflux=(emfx(i,1,1,3)+emfx(i,2,1,3))*0.25*weight
     unew(ind_buffer1,3+nvar)=unew(ind_buffer1,3+nvar)-dflux
     unew(ind_buffer2,3+neul)=unew(ind_buffer2,3+neul)-dflux
     unew(ind_buffer2,2+nvar)=unew(ind_buffer2,2+nvar)-dflux
     unew(ind_buffer3,2+neul)=unew(ind_buffer3,2+neul)-dflux
     if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
        unew(ind_buffer1,2+nvar)=unew(ind_buffer1,2+nvar)+dflux*0.5
        unew(ind_buffer3,3+neul)=unew(ind_buffer3,3+neul)+dflux*0.5
     endif
  end do

  ! Update coarse By and Bz using fine EMFx on Y=1 and Z=1 grid edge
  ind_father1=1+(i1  )+3*(j1  )+9*(k1+1)
  ind_father2=1+(i1  )+3*(j1+1)+9*(k1+1)
  ind_father3=1+(i1  )+3*(j1+1)+9*(k1  )
  do i=1,ncache
     ind_buffer1=nbors_father_cells(i,ind_father1)
     ind_buffer2=nbors_father_cells(i,ind_father2)
     ind_buffer3=nbors_father_cells(i,ind_father3)
     weight=1.0
     if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
     if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
     dflux=(emfx(i,1,3,3)+emfx(i,2,3,3))*0.25*weight
     unew(ind_buffer1,2+nvar)=unew(ind_buffer1,2+nvar)-dflux
     unew(ind_buffer2,2+neul)=unew(ind_buffer2,2+neul)-dflux
     unew(ind_buffer2,3+neul)=unew(ind_buffer2,3+neul)+dflux
     unew(ind_buffer3,3+nvar)=unew(ind_buffer3,3+nvar)+dflux
     if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
        unew(ind_buffer3,2+neul)=unew(ind_buffer3,2+neul)+dflux*0.5
        unew(ind_buffer1,3+neul)=unew(ind_buffer1,3+neul)-dflux*0.5
     endif
  end do

  ! Update coarse By and Bz using fine EMFx on Y=1 and Z=0 grid edge
  ind_father1=1+(i1  )+3*(j1+1)+9*(k1  )
  ind_father2=1+(i1  )+3*(j1+1)+9*(k1-1)
  ind_father3=1+(i1  )+3*(j1  )+9*(k1-1)
  do i=1,ncache
     ind_buffer1=nbors_father_cells(i,ind_father1)
     ind_buffer2=nbors_father_cells(i,ind_father2)
     ind_buffer3=nbors_father_cells(i,ind_father3)
     weight=1.0
     if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
     if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
     dflux=(emfx(i,1,3,1)+emfx(i,2,3,1))*0.25*weight
     unew(ind_buffer1,3+neul)=unew(ind_buffer1,3+neul)+dflux
     unew(ind_buffer2,3+nvar)=unew(ind_buffer2,3+nvar)+dflux
     unew(ind_buffer2,2+neul)=unew(ind_buffer2,2+neul)+dflux
     unew(ind_buffer3,2+nvar)=unew(ind_buffer3,2+nvar)+dflux
     if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
        unew(ind_buffer3,3+nvar)=unew(ind_buffer3,3+nvar)-dflux*0.5
        unew(ind_buffer1,2+neul)=unew(ind_buffer1,2+neul)-dflux*0.5
     endif
  end do

  !--------------------------------------
  ! Deal with 4 EMFy edges
  !--------------------------------------

  ! Update coarse Bx and Bz using fine EMFy on X=0 and Z=0 grid edge
  ind_father1=1+(i1  )+3*(j1  )+9*(k1-1)
  ind_father2=1+(i1-1)+3*(j1  )+9*(k1-1)
  ind_father3=1+(i1-1)+3*(j1  )+9*(k1  )
  do i=1,ncache
     ind_buffer1=nbors_father_cells(i,ind_father1)
     ind_buffer2=nbors_father_cells(i,ind_father2)
     ind_buffer3=nbors_father_cells(i,ind_father3)
     weight=1.0
     if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
     if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
     dflux=(emfy(i,1,1,1)+emfy(i,1,2,1))*0.25*weight
     unew(ind_buffer1,1+neul)=unew(ind_buffer1,1+neul)-dflux
     unew(ind_buffer2,1+nvar)=unew(ind_buffer2,1+nvar)-dflux
     unew(ind_buffer2,3+nvar)=unew(ind_buffer2,3+nvar)+dflux
     unew(ind_buffer3,3+neul)=unew(ind_buffer3,3+neul)+dflux
     if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
        unew(ind_buffer3,1+nvar)=unew(ind_buffer3,1+nvar)+dflux*0.5
        unew(ind_buffer1,3+nvar)=unew(ind_buffer1,3+nvar)-dflux*0.5
     endif
  end do

  ! Update coarse Bx and Bz using fine EMFy on X=0 and Z=1 grid edge
  ind_father1=1+(i1-1)+3*(j1  )+9*(k1  )
  ind_father2=1+(i1-1)+3*(j1  )+9*(k1+1)
  ind_father3=1+(i1  )+3*(j1  )+9*(k1+1)
  do i=1,ncache
     ind_buffer1=nbors_father_cells(i,ind_father1)
     ind_buffer2=nbors_father_cells(i,ind_father2)
     ind_buffer3=nbors_father_cells(i,ind_father3)
     weight=1.0
     if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
     if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
     dflux=(emfy(i,1,1,3)+emfy(i,1,2,3))*0.25*weight
     unew(ind_buffer1,3+nvar)=unew(ind_buffer1,3+nvar)+dflux
     unew(ind_buffer2,3+neul)=unew(ind_buffer2,3+neul)+dflux
     unew(ind_buffer2,1+nvar)=unew(ind_buffer2,1+nvar)+dflux
     unew(ind_buffer3,1+neul)=unew(ind_buffer3,1+neul)+dflux
     if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
        unew(ind_buffer3,3+neul)=unew(ind_buffer3,3+neul)-dflux*0.5
        unew(ind_buffer1,1+nvar)=unew(ind_buffer1,1+nvar)-dflux*0.5
     endif
  end do

  ! Update coarse Bx and Bz using fine EMFy on X=1 and Z=1 grid edge
  ind_father1=1+(i1  )+3*(j1  )+9*(k1+1)
  ind_father2=1+(i1+1)+3*(j1  )+9*(k1+1)
  ind_father3=1+(i1+1)+3*(j1  )+9*(k1  )
  do i=1,ncache
     ind_buffer1=nbors_father_cells(i,ind_father1)
     ind_buffer2=nbors_father_cells(i,ind_father2)
     ind_buffer3=nbors_father_cells(i,ind_father3)
     weight=1.0
     if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
     if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
     dflux=(emfy(i,3,1,3)+emfy(i,3,2,3))*0.25*weight
     unew(ind_buffer1,1+nvar)=unew(ind_buffer1,1+nvar)+dflux
     unew(ind_buffer2,1+neul)=unew(ind_buffer2,1+neul)+dflux
     unew(ind_buffer2,3+neul)=unew(ind_buffer2,3+neul)-dflux
     unew(ind_buffer3,3+nvar)=unew(ind_buffer3,3+nvar)-dflux
     if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
        unew(ind_buffer3,1+neul)=unew(ind_buffer3,1+neul)-dflux*0.5
        unew(ind_buffer1,3+neul)=unew(ind_buffer1,3+neul)+dflux*0.5
     endif
  end do

  ! Update coarse Bx and Bz using fine EMFy on X=1 and Z=0 grid edge
  ind_father1=1+(i1+1)+3*(j1  )+9*(k1  )
  ind_father2=1+(i1+1)+3*(j1  )+9*(k1-1)
  ind_father3=1+(i1  )+3*(j1  )+9*(k1-1)
  do i=1,ncache
     ind_buffer1=nbors_father_cells(i,ind_father1)
     ind_buffer2=nbors_father_cells(i,ind_father2)
     ind_buffer3=nbors_father_cells(i,ind_father3)
     weight=1.0
     if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
     if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
     dflux=(emfy(i,3,1,1)+emfy(i,3,2,1))*0.25*weight
     unew(ind_buffer1,3+neul)=unew(ind_buffer1,3+neul)-dflux
     unew(ind_buffer2,3+nvar)=unew(ind_buffer2,3+nvar)-dflux
     unew(ind_buffer2,1+neul)=unew(ind_buffer2,1+neul)-dflux
     unew(ind_buffer3,1+nvar)=unew(ind_buffer3,1+nvar)-dflux
     if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
        unew(ind_buffer3,3+nvar)=unew(ind_buffer3,3+nvar)+dflux*0.5
        unew(ind_buffer1,1+neul)=unew(ind_buffer1,1+neul)+dflux*0.5
     endif
  end do
#endif

  end if

end subroutine godfine1
!###########################################################
!###########################################################
!###########################################################
!###########################################################
#if USE_FLD==1
subroutine rad_force_fine(ilevel)
  use amr_commons
  use hydro_commons
  use cooling_module,ONLY:kB,mH
  use constants, only : c_cgs
  use radiation_parameters,ONLY:Tr_floor,eray_min,nu_min_hz,nu_max_hz,frad
  use const
  use units_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine sets array uold to its new value unew after the
  ! hydro step.
  !--------------------------------------------------------------------------
  integer::i,j,k,ivar,ind,iskip,nx_loc,info
  real(dp)::scale,d,u,v,w,A,B,C,d_old
  real(dp)::e_mag,e_kin,e_cons,e_prim,e_trunc,div,dx,fact,e_r

  integer ,dimension(1:nvector),save::ind_grid,ind_cell
  integer ,dimension(1:nvector,0:twondim),save::igridn
  integer ,dimension(1:nvector,1:ndim),save::ind_left,ind_right
  real(dp),dimension(1:nvector,1:ndim,1:ngrp),save::Erg,Erd
  real(dp),dimension(1:nvector,1:ndim,1:ndim),save::velg,veld
  real(dp),dimension(1:nvector,1:ndim)::dx_g,dx_d
  real(dp)::rosseland_ana
  real(dp)::Pgdivu,u_square,d_loc,Tp_loc,Tr_loc,cal_Teg

  integer::ncache,igrid,ngrid,idim,id1,ig1,ih1,id2,ig2,ih2,igroup
  integer  ,dimension(1:3,1:2,1:8)::iii,jjj
  real(dp)::dx_loc,surf_loc,vol_loc,usquare,emag,erad_loc,ekin,eps,cv,pp_eos
  real(dp)::kappa_R,gradEr_norm,gradEr_norm2,R,lambda,lambda_fld,chi,PgmErdivu,gradEru
  real(dp) ,dimension(1:3)::skip_loc
  real(dp) ,dimension(1:ndim,1:ngrp)::gradEr
  real(dp) ,dimension(1:ndim,1:ndim)::divu_loc
  real(dp) ,dimension(1:ndim,1:ndim,1:ngrp)::Pg
  real(dp) ,dimension(1:ndim       )::u_loc
  real(dp) :: nuPrDivu,nuPr,nuPl,Pr_nu
  real(dp), dimension(1:5) :: Pr_temp

  !  EOS
  real(dp) :: dd,ee,cmp_Cv_eos
  integer  :: ht 

#if USE_FLD==1
  real(dp)::scale_nH,scale_T2,scale_t,scale_v,scale_d,scale_l
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
#endif

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5d0**ilevel

  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale ! Warning: scale factor already done in dx
  vol_loc=dx_loc**ndim
  surf_loc=dx_loc**(ndim-1)

  iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

if(fld)then    
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
        
#if USE_FLD==1        
        ! Gather neighboring temperature
        do idim=1,ndim
           id1=jjj(idim,1,ind); ig1=iii(idim,1,ind)
           ih1=ncoarse+(id1-1)*ngridmax
           do i=1,ngrid
              if(igridn(i,ig1)>0)then
                 do igroup=1,ngrp
                    Erg(i,idim,igroup) = max(uold(igridn(i,ig1)+ih1,firstindex_er+igroup),eray_min/(scale_d*scale_v**2))
                 end do
                 velg(i,idim,1:ndim) = uold(igridn(i,ig1)+ih1,2:ndim+1)/uold(igridn(i,ig1)+ih1,1)
                 dx_g(i,idim) = dx_loc
              else
                 do igroup=1,ngrp
                    Erg(i,idim,igroup) = max(uold(ind_left(i,idim),firstindex_er+igroup),eray_min/(scale_d*scale_v**2))
                 end do
                 velg(i,idim,1:ndim) = uold(ind_left(i,idim),2:ndim+1)/uold(ind_left(i,idim),1)
                 dx_g(i,idim) = dx_loc*1.5_dp
              end if
           enddo
           id2=jjj(idim,2,ind); ig2=iii(idim,2,ind)
           ih2=ncoarse+(id2-1)*ngridmax
           do i=1,ngrid
              if(igridn(i,ig2)>0)then
                 do igroup=1,ngrp
                    Erd(i,idim,igroup) = max(uold(igridn(i,ig2)+ih2,firstindex_er+igroup),eray_min/(scale_d*scale_v**2))
                 end do
                 veld(i,idim,1:ndim)= uold(igridn(i,ig2)+ih2,2:ndim+1)/uold(igridn(i,ig2)+ih2,1)
                 dx_d(i,idim)=dx_loc
              else 
                 do igroup=1,ngrp
                    Erd(i,idim,igroup) = max(uold(ind_right(i,idim),firstindex_er+igroup),eray_min/(scale_d*scale_v**2))
                 end do
                 veld(i,idim,1:ndim)= uold(ind_right(i,idim),2:ndim+1)/uold(ind_right(i,idim),1)
                 dx_d(i,idim)=dx_loc*1.5_dp
              end if
           enddo
        end do
       ! End loop over dimensions
  
        do i=1,ngrid
           !compute divu
           do j=1,ndim
              do k=1,ndim
                 divu_loc(j,k) = (veld(i,j,k)-velg(i,j,k))/(dx_g(i,j)+dx_d(i,j))
              enddo
              do igroup=1,ngrp
                 gradEr(j,igroup) = (Erd(i,j,igroup)-Erg(i,j,igroup))/(dx_g(i,j)+dx_d(i,j))
              enddo
           enddo

           d_loc = uold(ind_cell(i),1)*scale_d
           u_loc(1:ndim) = uold(ind_cell(i),2:ndim+1)/uold(ind_cell(i),1)
           
           usquare=0.0
           do idim=1,ndim
              usquare=usquare+(uold(ind_cell(i),idim+1)/uold(ind_cell(i),1))**2
           end do
           
           ! Compute total magnetic energy
           emag = 0.0d0
           do ivar=1,3
              emag = emag + 0.125d0*(uold(ind_cell(i),5+ivar) &
                   &  +uold(ind_cell(i),nvar+ivar))**2
           end do
           erad_loc=0.0D0
#if NENER>0
           do igroup=1,nener
              erad_loc = erad_loc + uold(ind_cell(i),8+igroup)
           end do
#endif
           d     = uold(ind_cell(i),1)
           ekin  = d*usquare/2.0
           ! Compute gas temperature in cgs
           eps   = uold(ind_cell(i),5)-ekin-emag-erad_loc
           if(energy_fix)eps   = uold(ind_cell(i),nvar) ! comment this for radiative shock
           ! Compute gas temperature in cgs
           call temperature_eos(d,eps,Tp_loc,ht)

           frad(ind_cell(i),1:ndim)=0.0d0
           
           ! Compute radiative pressure in all groups
           do igroup=1,ngrp
              
              ! Compute radiative pressure
              Tr_loc = cal_Teg(uold(ind_cell(i),firstindex_er+igroup)*scale_d*scale_v**2,igroup)              
              kappa_R = rosseland_ana(d_loc,Tp_loc,Tr_loc,igroup,in_sink(ind_cell(i)))/scale_kappa
              gradEr_norm2 = (sum(gradEr(1:ndim,igroup)**2))
              gradEr_norm  = (gradEr_norm2)**0.5
              R =   max(1.d-10,gradEr_norm/(max(uold(ind_cell(i),firstindex_er+igroup),eray_min/(scale_d*scale_v**2))*kappa_R))
              lambda = lambda_fld(R)
              chi = lambda + (lambda*R)**2
              
              frad(ind_cell(i),1:ndim) =  frad(ind_cell(i),1:ndim) + lambda*gradEr(1:ndim,igroup)/d
           enddo !end loop over rad groups

        end do
#endif
#if USE_M_1==1
        do i=1,ngrid
           ! Compute density and temperature for opacity
           d_loc = uold(ind_cell(i),1)*scale_d
           
           usquare=zero
           do idim=1,ndim
              usquare=usquare+(uold(ind_cell(i),idim+1)/uold(ind_cell(i),1))**2
           end do

           emag = zero
           do ivar=1,3
              emag = emag + 0.125d0*(uold(ind_cell(i),5+ivar) &
                   &  +uold(ind_cell(i),nvar+ivar))**2
           end do
           erad_loc=zero
#if NENER>0
           do igroup=1,nener
              erad_loc = erad_loc + uold(ind_cell(i),8+igroup)
           end do
#endif
           d     = uold(ind_cell(i),1)
           ekin  = d*usquare/2.0
           ! Compute gas temperature in cgs
           eps   = uold(ind_cell(i),5)-ekin-emag-erad_loc
           call temperature_eos(d,eps,Tp_loc,ht)

           frad(ind_cell(i),1:ndim)=zero

           do igroup=1,ngrp

              Tr_loc = cal_Teg(uold(ind_cell(i),firstindex_er+igroup)*scale_d*scale_v**2,igroup)
              kappa_R = rosseland_ana(d_loc,Tp_loc,Tr_loc,igroup,in_sink(ind_cell(i)))/scale_kappa

              ! divide by d because equation over u and not d*u
              frad(ind_cell(i),1) =  frad(ind_cell(i),1) + kappa_R*uold(ind_cell(i),firstindex_fr+igroup)/(c_cgs/scale_v)/d
              frad(ind_cell(i),2) =  frad(ind_cell(i),2) + kappa_R*uold(ind_cell(i),firstindex_fr+igroup+ngrp)/(c_cgs/scale_v)/d
              frad(ind_cell(i),3) =  frad(ind_cell(i),3) + kappa_R*uold(ind_cell(i),firstindex_fr+igroup+2*ngrp)/(c_cgs/scale_v)/d

           enddo !end loop over rad groups
        enddo
#endif

     enddo
     ! End loop over cells
  end do
  ! End loop over grids
endif
  
111 format('   Entering rad_force_fine for level ',i2)

end subroutine rad_force_fine
#endif
