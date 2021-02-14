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
  integer::i,igrid,ncache,ngrid,imat
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
  do imat=1,nmat
     call make_virtual_reverse_dp(dive(1,imat),ilevel)
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
  integer::i,ivar,ind,icpu,iskip,imat
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
     do imat=1,nmat
        do i=1,active(ilevel)%ngrid
           dive(active(ilevel)%igrid(i)+iskip,imat) = 0
        end do
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
     do imat=1,nmat
        do i=1,reception(icpu,ilevel)%ngrid
           dive(reception(icpu,ilevel)%igrid(i)+iskip,imat) = 0
        end do
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
     do ivar=1,nvar
        do i=1,active(ilevel)%ngrid
           uold(active(ilevel)%igrid(i)+iskip,ivar) = unew(active(ilevel)%igrid(i)+iskip,ivar)
        end do
     end do
  end do

  if(nmat>1)call pressure_relaxation2(ilevel)
  
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
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp),dimension(1:nvector),save::ekin,dtot,ptot,divpv
  real(dp),dimension(1:nvector),save::gg_mat,ee_mat,pp_mat,cc_mat
  real(dp),dimension(1:nvector,1:nmat),save::ff,gg,ee,pp,yy
  real(dp),dimension(1:nvector,1:ndim),save::vv
  real(dp)::skip_loc,dx,eps,scale,dx_loc
  real(dp)::one=1.0_dp, half=0.5_dp, zero=0.0_dp
  real(dp),dimension(1:8)::xc
  integer::ix,iy,iz,nx_loc
  logical::error,inv

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

        ! Non-conservative source term for volume fraction
        do imat=1,nmat
           ivar=imat
           do i=1,ngrid
              unew(ind_cell(i),ivar) = unew(ind_cell(i),ivar) - uold(ind_cell(i),ivar)*divu(ind_cell(i))
           end do
        end do
      
        ! Source terms for partial energies

        ! Compute divergence of total pressure
        divpv(1:ngrid)=0.0d0
        do imat=1,nmat
           do i=1,ngrid
              divpv(i) = divpv(i) + dive(ind_cell(i),imat)
           end do
        end do

        ! Compute volume fraction and fluid density
        dtot(1:ngrid)=0
        do imat = 1,nmat
           do i = 1,ngrid
              ! Volume fractions
              ff(i,imat) = uold(ind_cell(i),imat)
              ! True densities
              gg(i,imat) = uold(ind_cell(i),nmat+imat)/ff(i,imat)
              ! Total density
              dtot(i) = dtot(i) + uold(ind_cell(i),nmat+imat)
           end do
        end do

        ! Compute mass fraction
        do imat = 1,nmat
           do i = 1,ngrid
              yy(i,imat) = ff(i,imat)*gg(i,imat)/dtot(i)
           end do
        end do

        ! Compute velocity
        do idim = 1,ndim
           do i = 1,ngrid
              vv(i,idim) = uold(ind_cell(i),2*nmat+idim)/dtot(i)
           end do
        end do
        
        ! Compute specific kinetic energy
        ekin(1:ngrid)=0.0
        do idim = 1,ndim
           do i = 1,ngrid
              ekin(i) = ekin(i) + half*vv(i,idim)**2 ! This is 0.5*u^2
           end do
        end do
        
        ! Compute individual internal energies
        do imat=1,nmat
           do i=1,ngrid
              ee(i,imat) = uold(ind_cell(i),2*nmat+ndim+imat)/ff(i,imat) - gg(i,imat)*ekin(i)
           end do
        end do

        ! Call eos routine
        inv=.false.
        ptot(1:ngrid)=0.0
        do imat=1,nmat
           do i=1,ngrid
              gg_mat(i) = gg(i,imat)
              ee_mat(i) = ee(i,imat)
           end do
           call eos(gg_mat,ee_mat,pp_mat,cc_mat,imat,inv,ngrid)
           do i=1,ngrid
              ! Individual pressures
              pp(i,imat) = pp_mat(i)
              ptot(i) = ptot(i) + ff(i,imat)*pp(i,imat)
           end do
        end do

        ! Update new partial total energies using non-conservative source terms
        do imat=1,nmat
           ivar=2*nmat+ndim+imat
           do i=1,ngrid
              unew(ind_cell(i),ivar)=unew(ind_cell(i),ivar) &
                   & - dive(ind_cell(i),imat) &
                   & + yy(i,imat)*divpv(i) &
                   & + ff(i,imat)*pp(i,imat)*divu(ind_cell(i)) &
                   & - yy(i,imat)*ptot(i)*divu(ind_cell(i))
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
  integer::i,ind,iskip,ind_cell,imat
  real(dp)::d,d_mat,u,v,w,d_old,fact
  real(dp),dimension(1:nmat)::e_kin,e_prim
  
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel
  
  ! Add gravity source term at time t with half time step
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,active(ilevel)%ngrid
        ind_cell = active(ilevel)%igrid(i) + iskip
        
        d = 0
        do imat=1,nmat
          d = d + max(unew(ind_cell,nmat+imat),smallr)
        end do

        u=0; v=0; w=0
        if(ndim>0) u = unew(ind_cell,2*nmat+1)/max(d,smallr)
        if(ndim>1) v = unew(ind_cell,2*nmat+2)/max(d,smallr)
        if(ndim>2) w = unew(ind_cell,2*nmat+3)/max(d,smallr)
        
        do imat=1,nmat
          d_mat         = unew(ind_cell,nmat+imat)/max(unew(ind_cell,imat),smallf)
          e_kin(imat)   = 0.5d0*d_mat*(u**2+v**2+w**2)
          e_prim(imat)  = unew(ind_cell,2*nmat+ndim+imat)/max(unew(ind_cell,imat),smallf) - e_kin(imat)
        end do
        
        d_old = 0
        do imat=1,nmat
          d_old = d_old + max(uold(ind_cell,nmat+imat),smallr)
        end do
        fact    = d_old/d*0.5d0*dtnew(ilevel)
        
        if(ndim>0)then
           u = u + f(ind_cell,1)*fact
           unew(ind_cell,2*nmat+1) = d*u
        endif
        if(ndim>1)then
           v = v + f(ind_cell,2)*fact
           unew(ind_cell,2*nmat+2) = d*v
        end if
        if(ndim>2)then
           w = w + f(ind_cell,3)*fact
           unew(ind_cell,2*nmat+3) = d*w
        endif
        
        do imat=1,nmat
          d_mat       = unew(ind_cell,nmat+imat)/max(unew(ind_cell,imat),smallf)
          e_kin(imat) = 0.5d0*d_mat*(u**2+v**2+w**2)
          unew(ind_cell,2*nmat+ndim+imat) = unew(ind_cell,imat)*(e_prim(imat) + e_kin(imat))
        end do
     end do
  end do
  
111 format('   Entering add_gravity_source_terms for level ',i2)

end subroutine add_gravity_source_terms
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine pressure_relaxation(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !-------------------------------------------------------------------
  ! Update volume fraction using compressibility source terms
  !-------------------------------------------------------------------
  integer::i,ivar,imat,idim,ind,iskip,ncache,igrid,ngrid
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp),dimension(1:nvector),save::ekin,dtot,etot,ptot
  real(dp),dimension(1:nvector),save::gg_mat,ee_mat,pp_mat,cc_mat
  real(dp),dimension(1:nvector),save::ec_hat,pc_hat,alpha_hat
  real(dp),dimension(1:nvector,1:nmat),save::ff,gg,ee,pp,yy
  real(dp),dimension(1:nvector,1:ndim),save::vv
  real(dp)::smallgamma,biggamma,p_0,e_c,p_c,a0,rho_0,eta
  real(dp)::skip_loc,dx,eps,scale,dx_loc
  real(dp)::one=1.0_dp, half=0.5_dp, zero=0.0_dp
  real(dp),dimension(1:8)::xc
  integer::ix,iy,iz,nx_loc
  logical::error,inv

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

        ! Compute volume fraction and fluid density
        dtot(1:ngrid)=0
        do imat = 1,nmat
           do i = 1,ngrid
              ! Volume fractions
              ff(i,imat) = uold(ind_cell(i),imat)
              ! True densities
              gg(i,imat) = uold(ind_cell(i),nmat+imat)/ff(i,imat)
              ! Total density
              dtot(i) = dtot(i) + uold(ind_cell(i),nmat+imat)
           end do
        end do

        ! Compute velocity
        do idim = 1,ndim
           do i = 1,ngrid
              vv(i,idim) = uold(ind_cell(i),2*nmat+idim)/dtot(i)
           end do
        end do
        
        ! Compute specific kinetic energy
        ekin(1:ngrid)=0.0
        do idim = 1,ndim
           do i = 1,ngrid
              ekin(i) = ekin(i) + half*vv(i,idim)**2 ! This is 0.5*u^2
           end do
        end do
        
        ! Compute total internal energy
        etot(1:ngrid)=0.0
        do imat=1,nmat
           do i=1,ngrid
              etot(i) = etot(i) + uold(ind_cell(i),2*nmat+ndim+imat)
           end do
        end do
        etot(1:ngrid)=etot(1:ngrid)-dtot(1:ngrid)*ekin(1:ngrid)
        
        ! Mie-Gruneisen
        alpha_hat(1:ngrid)=zero
        pc_hat (1:ngrid)=zero
        ec_hat(1:ngrid)=zero
        do imat = 1,nmat
           ! Get Mie-Grueneisen EOS parameters
           smallgamma=eos_params(imat,1);biggamma=eos_params(imat,2);p_0=eos_params(imat,3);rho_0=eos_params(imat,4)
           ! P - P_c = (gamma - one) * (e - e_c) ; e = P/(gamma-1) + (e_c-P_c/(gamma-1))
           a0 = one / (smallgamma-one)
           do i = 1,ngrid
              ! Update Mie-Gruneisen terms for each material
              eta = max(gg(i,imat),smallr)/rho_0
              p_c = p_0 * eta**biggamma
              e_c = p_c / (biggamma-one)
              ! Update total values
              alpha_hat(i) = alpha_hat(i) + ff(i,imat) * a0
              ec_hat(i) = ec_hat(i) + ff(i,imat) * e_c
              pc_hat(i) = pc_hat(i) + ff(i,imat) * p_c * a0
           end do
        end do
        pc_hat(1:ngrid) = pc_hat(1:ngrid)/alpha_hat(1:ngrid)
        
        ! Calculate pressure for given internal energy
        do i = 1,ngrid
           ptot(i) = (etot(i) - ec_hat(i)) / alpha_hat(i) + pc_hat(i)
        end do
        
        ! Call inverse eos routine
        inv=.true.
        do imat=1,nmat
           do i=1,ngrid
              gg_mat(i) = gg(i,imat)
              pp_mat(i) = ptot(i)
           end do
           call eos(gg_mat,ee_mat,pp_mat,cc_mat,imat,inv,ngrid)
           do i=1,ngrid
              ! Individual internal energies
              ee(i,imat) = ee_mat(i)
           end do
        end do

        ! Update new partial total energies
        do imat=1,nmat
           ivar=2*nmat+ndim+imat
           do i=1,ngrid
              if(son(ind_cell(i))==0)then
                 uold(ind_cell(i),ivar)=ff(i,imat)*ee(i,imat)+ff(i,imat)*gg(i,imat)*ekin(i)
              endif
           end do
        end do
      
     end do
     ! End loop over cells
  end do
  ! End loop over grids
     
end subroutine pressure_relaxation
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine pressure_relaxation2(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !-------------------------------------------------------------------
  ! Update volume fraction using compressibility source terms
  !-------------------------------------------------------------------
  integer::i,ivar,imat,idim,ind,iskip,ncache,igrid,ngrid
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp)::ekin,dtot,etot,ptot
  real(dp),dimension(1:nmat),save::ff,gg,ee,pp,rc2
  real(dp),dimension(1:nmat),save::ff_new,ee_new,w
  real(dp),dimension(1:ndim),save::vv
  real(dp)::smallgamma,biggamma,p_0,e_c,p_c,delpc,a0,rho_0,eta
  real(dp)::skip_loc,dx,eps,scale,dx_loc,dd,ddd
  real(dp)::t12,t13,t14,t23,t24,t34
  real(dp)::one=1.0_dp, half=0.5_dp, zero=0.0_dp
  real(dp),dimension(1:8)::xc
  integer::ix,iy,iz,nx_loc,iter
  logical::error,inv

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

        do i=1,ngrid
           if(son(ind_cell(i))==0)then
              
              ! Compute total mass density
              dtot=0
              do imat = 1,nmat
                 dtot = dtot + uold(ind_cell(i),nmat+imat)
              end do

              ! Compute velocity
              do idim = 1,ndim
                 vv(idim) = uold(ind_cell(i),2*nmat+idim)/dtot
              end do
        
              ! Compute specific kinetic energy
              ekin=0.0
              do idim = 1,ndim
                 ekin = ekin + half*vv(idim)**2
              end do

              do iter=1,10
              
              ! Compute volume fraction and true density
              do imat = 1,nmat
                 ff(imat) = uold(ind_cell(i),imat)
                 gg(imat) = uold(ind_cell(i),nmat+imat)/ff(imat)
              end do

              ! Compute specific internal energy
              do imat = 1,nmat
                 ee(imat) = uold(ind_cell(i),2*nmat+ndim+imat)/ff(imat)/gg(imat)-ekin
              end do
        
              ! Mie-Gruneisen
              ptot = 0.0
              do imat = 1,nmat
                 smallgamma=eos_params(imat,1);biggamma=eos_params(imat,2);p_0=eos_params(imat,3);rho_0=eos_params(imat,4)
                 eta = gg(imat)/rho_0
                 p_c = p_0 * eta**biggamma
                 e_c = p_c / (biggamma-one)
                 delpc = biggamma * p_c 
                 pp(imat) = (smallgamma-1)*(gg(imat)*ee(imat)-e_c) + p_c
                 rc2(imat) = delpc + smallgamma * (pp(imat)-p_c)
                 ptot = ptot + ff(imat)*pp(imat)
              end do

              do imat = 1,nmat
                 smallgamma = eos_params(imat,1)
                 rc2(imat) = rc2(imat) + (smallgamma-1)*(ptot-pp(imat))
              end do

              ! Compute weights
              do imat = 1,nmat
                 w(imat) = rc2(imat)/ff(imat)
              end do

              ! Compute new volume fraction
#if NMAT==2
              dd=w(1)+w(2)
              
              ff_new(1) = ff(1) + (pp(1)-pp(2))/dd
              ff_new(2) = ff(2) + (pp(2)-pp(1))/dd
#endif
#if NMAT==3
              dd=3*w(1)*w(2)+3*w(1)*w(3)+3*w(2)*w(3)
              
              ff_new(1) = ff(1) + ((w(2)+2*w(3))*(pp(1)-pp(2))+(2*w(2)+w(3))*(pp(1)-pp(3))+(w(2)-w(3))*(pp(2)-pp(3)))/dd
              ff_new(2) = ff(2) + ((w(3)+2*w(1))*(pp(2)-pp(3))+(2*w(3)+w(1))*(pp(2)-pp(1))+(w(3)-w(1))*(pp(3)-pp(1)))/dd
              ff_new(3) = ff(3) + ((w(1)+2*w(2))*(pp(3)-pp(1))+(2*w(1)+w(2))*(pp(3)-pp(2))+(w(1)-w(2))*(pp(1)-pp(2)))/dd
#endif
#if NMAT==4
              t12=w(1)*w(2)
              t13=w(1)*w(3)
              t14=w(1)*w(4)
              t23=w(2)*w(3)
              t24=w(2)*w(4)
              t34=w(3)*w(4)
              dd=4*w(1)*w(2)*w(3)+4*w(1)*w(2)*w(4)+4*w(1)*w(3)*w(4)+4*w(2)*w(3)*w(4)
              
              ff_new(1) = ff(1) + ((t23+t24+2*t34)*(pp(1)-pp(2))+(t23+2*t24+t34)*(pp(1)-pp(3))+(2*t23+t24+t34)*(pp(1)-pp(4))+ &
                   &                     (t24-t34)*(pp(2)-pp(3))+      (t23-t34)*(pp(2)-pp(4))+      (t23-t24)*(pp(3)-pp(4)))/dd
              ff_new(2) = ff(2) + ((t34+t13+2*t14)*(pp(2)-pp(3))+(t34+2*t13+t14)*(pp(2)-pp(4))+(2*t34+t13+t14)*(pp(2)-pp(1))+ &
                   &                     (t13-t14)*(pp(3)-pp(4))+      (t34-t14)*(pp(3)-pp(1))+      (t34-t13)*(pp(4)-pp(1)))/dd
              ff_new(3) = ff(3) + ((t14+t24+2*t12)*(pp(3)-pp(4))+(t14+2*t24+t12)*(pp(3)-pp(1))+(2*t14+t24+t12)*(pp(3)-pp(2))+ &
                   &                     (t24-t12)*(pp(4)-pp(1))+      (t14-t12)*(pp(4)-pp(2))+      (t14-t24)*(pp(1)-pp(2)))/dd
              ff_new(4) = ff(4) + ((t12+t13+2*t23)*(pp(4)-pp(1))+(t12+2*t13+t23)*(pp(4)-pp(2))+(2*t12+t13+t23)*(pp(4)-pp(3))+ &
                   &                     (t13-t23)*(pp(1)-pp(2))+      (t12-t23)*(pp(1)-pp(3))+      (t12-t13)*(pp(2)-pp(3)))/dd
#endif
              ! Compute new specific internal energy
              do imat = 1,nmat
                 ee_new(imat) = ee(imat) - ptot * (ff_new(imat) - ff(imat)) / (ff(imat)*gg(imat))
              end do

              ! Store new volume fraction
              do imat = 1,nmat
                 uold(ind_cell(i),imat) = ff_new(imat)
              end do

              ! Store new partial total energy
              do imat=1,nmat
                 uold(ind_cell(i),2*nmat+ndim+imat) = ff(imat)*gg(imat)*(ee_new(imat) + ekin)
              end do

!              write(*,*)iter,ff(1),ff(2),pp(1),pp(2)

              end do
              
           endif
        end do
      
     end do
     ! End loop over cells
  end do
  ! End loop over grids
     
end subroutine pressure_relaxation2
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
  ! to set initial conditions in a 6x6x6 grid. It interpolates from
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
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:nvar,1:ndim),save::flux
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:nmat+1,1:ndim),save::tmp
  logical ,dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::ok

  integer,dimension(1:nvector),save::igrid_nbor,ind_cell,ind_buffer
  logical,dimension(1:nvector),save::exist_nbor

  integer::i,j,ivar,imat,idim,ind_son,ind_father,iskip,nbuffer,ibuffer,nx_loc
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

  !-----------------------------------------------
  ! Compute flux using second-order Godunov method
  !-----------------------------------------------
  flux=0.0d0; tmp=0.0d0
  call unsplit(uloc,gloc,flux,tmp,dx_loc,dx_loc,dx_loc,dtnew(ilevel),ncache)

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
        do ivar=1,nmat+1
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
        ! Update energy divergence
        do imat=1,nmat
           do i=1,ncache
              dive(ind_cell(i),imat)=dive(ind_cell(i),imat)+ &
                   & (tmp(i,i3   ,j3   ,k3   ,1+imat,idim) &
                   & -tmp(i,i3+i0,j3+j0,k3+k0,1+imat,idim))
           end do
        end do
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
     ! Update energy divergence
     do k3=k3min,k3max-k0
     do j3=j3min,j3max-j0
     do i3=i3min,i3max-i0
        do imat=1,nmat
           do i=1,ncache
              if(.not.exist_nbor(i))then
                 dive(ind_buffer(i),imat)=dive(ind_buffer(i),imat) &
                      & -tmp(i,i3,j3,k3,1+imat,idim)/dble(twotondim)
              end if
           end do
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
     ! Update energy divergence
     do k3=k3min+k0,k3max
     do j3=j3min+j0,j3max
     do i3=i3min+i0,i3max
        do imat=1,nmat
           do i=1,ncache
              if(.not.exist_nbor(i))then
                 dive(ind_buffer(i),imat)=dive(ind_buffer(i),imat) &
                      & +tmp(i,i3+i0,j3+j0,k3+k0,1+imat,idim)/dble(twotondim)
              end if
           end do
        end do
     end do
     end do
     end do

  end do
  ! End loop over dimensions

end subroutine godfine1
