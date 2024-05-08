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

  if(nmat>1)then
     call pressure_relaxation2(ilevel)
     call phase_transition(ilevel)
  end if

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
  real(dp)::d,d_mat,u,v,w,e_kin,e_tot,d_old,fact
  real(dp)::one=1.0_dp,half=0.5_dp,zero=0.0_dp
  real(dp),dimension(1:nmat)::e_prim
#if NVAR > NDIM + 3*NMAT
  integer::ipscal,npscal
  real(dp)::s_entry,e_th,e_cold
  real(dp),dimension(1:nvector),save::g_mat,e_mat,s_mat
#endif
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Add gravity source term at time t with half time step
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,active(ilevel)%ngrid
        ind_cell = active(ilevel)%igrid(i) + iskip

        d = 0
        do imat = 1,nmat
          d = d + unew(ind_cell,nmat+imat)
        end do

        u=0; v=0; w=0; e_kin=0.0d0
        if(ndim>0) u = unew(ind_cell,2*nmat+1)/d
        e_kin = e_kin + half*u**2
        if(ndim>1) v = unew(ind_cell,2*nmat+2)/d
        e_kin = e_kin + half*v**2
        if(ndim>2) w = unew(ind_cell,2*nmat+3)/d
        e_kin = e_kin + half*w**2

#if NVAR > NDIM + 3*NMAT
        do imat=1,nmat

           ! Check that entropy is correct on entry
           s_entry = unew(ind_cell,3*nmat+ndim+imat)/unew(ind_cell,nmat+imat)
           if(s_entry<0)then
              write(*,*)'gravity: negative entropy on entry',imat,s_entry
           endif

           ! Compute thermal energy using the total energy
           e_tot = unew(ind_cell,2*nmat+ndim+imat)/unew(ind_cell,nmat+imat)
           e_prim(imat) = e_tot - e_kin
           g_mat(1) = unew(ind_cell,nmat+imat)/unew(ind_cell,imat)
           e_mat(1) = e_prim(imat)*g_mat(1)
           call eos_s(g_mat,e_mat,s_mat,imat,.false.,1)
           e_th = s_mat(1)*g_mat(1)**eos_params(imat,1)
           e_cold=eos_params(imat,3)*(g_mat(1)/eos_params(imat,4))**eos_params(imat,2)
           if(e_th < 0.001*e_kin .or. e_th < 0.001*e_cold)then
!           if(.true.)then
              ! If thermal energy is too small then use entropy instead
              s_mat(1) = s_entry
              call eos_s(g_mat,e_mat,s_mat,imat,.true.,1)
              e_prim(imat) = e_mat(1)/g_mat(1)
           else
              ! Otherwise keep thermal energy but update entropy accordingly
              unew(ind_cell,3*nmat+ndim+imat) = unew(ind_cell,nmat+imat)*s_mat(1)
           endif
        end do
#endif

        d_old = 0
        do imat = 1,nmat
          d_old = d_old + uold(ind_cell,nmat+imat)
        end do
        fact = d_old/d*0.5d0*dtnew(ilevel)

        e_kin=0d0
        if(ndim>0)then
           u = u + f(ind_cell,1)*fact
           unew(ind_cell,2*nmat+1) = d*u
           e_kin = e_kin + half*u**2
        endif
        if(ndim>1)then
           v = v + f(ind_cell,2)*fact
           unew(ind_cell,2*nmat+2) = d*v
           e_kin = e_kin + half*v**2
        end if
        if(ndim>2)then
           w = w + f(ind_cell,3)*fact
           unew(ind_cell,2*nmat+3) = d*w
           e_kin = e_kin + half*w**2
        endif

        do imat=1,nmat
           e_tot = e_prim(imat) + e_kin
           unew(ind_cell,2*nmat+ndim+imat) = unew(ind_cell,nmat+imat)*e_tot
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
  real(dp)::one=1.0_dp,half=0.5_dp,zero=0.0_dp
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
  real(dp)::dtot,ptot
  real(dp),dimension(1:nmat),save::ff_old,ff,gg,ee,pp,rc2
  real(dp),dimension(1:nmat),save::ff_new,ee_new,w,pp_new
  real(dp),dimension(1:ndim),save::vv
  real(dp)::error
  real(dp)::smallgamma,biggamma,p_0,e_c,p_c,delpc,a0,rho_0,eta,q
  real(dp)::E_1,E_2,A_1,A_2,C_v,T_0,E_0,p_c_1,p_c_2
  real(dp)::skip_loc,dx,eps,scale,dx_loc,dd,ddd
  real(dp)::t12,t13,t14,t23,t24,t34
  real(dp)::t123,t124,t125,t134,t135,t145
  real(dp)::t234,t235,t245,t345
  real(dp)::one=1.0_dp,half=0.5_dp, zero=0.0_dp
  real(dp),dimension(1:8)::xc
  integer::ix,iy,iz,nx_loc,iter
  integer::iter_max=100,iter_mean
  logical::inv
  real(dp),dimension(1:nvector),save::g_mat,e_mat
  real(dp)::e_kin,e_tot,ee_mat,cc_mat
#if NVAR > NDIM + 3*NMAT
  real(dp),dimension(1:nvector),save::s_mat
  real(dp)::s_entry,e_th,e_cold
#endif

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

              ! Compute old volume fraction
              do imat = 1,nmat
                 ff_old(imat) = uold(ind_cell(i),imat)
              end do

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
              e_kin=0.0
              do idim = 1,ndim
                 e_kin = e_kin + half*vv(idim)**2
              end do

              iter = 0
              error = 1
              do while(error.GT.1d-6.AND.iter<iter_max)

                 ! Compute volume fraction and true density
                 do imat = 1,nmat
                    ff(imat) = uold(ind_cell(i),imat)
                    gg(imat) = uold(ind_cell(i),nmat+imat)/ff(imat)
                 end do

                 ! Compute specific internal energy
                 do imat = 1,nmat
                    e_tot = uold(ind_cell(i),2*nmat+ndim+imat)/uold(ind_cell(i),nmat+imat)
                    ee(imat) = e_tot-e_kin
                    if(ee(imat)<0)then
                       ! write(*,*)'relaxation: negative energy',iter,ff(1),ff(2),gg(1),gg(2),ee(imat),e_tot,e_kin
                    endif
                 end do

                 ptot = 0.0
                 do imat = 1,nmat
                    if(eos_name =='stiffened gas')then
                      ! Get Stiffened Gas EOS parameters
                      smallgamma=eos_params(imat,1);q=eos_params(imat,2);p_0=eos_params(imat,3)
                      pp(imat)  = (smallgamma-1)*gg(imat)*(ee(imat)-q) - smallgamma*p_0
                      rc2(imat) = smallgamma * (pp(imat)+p_0)
                      ptot = ptot + ff(imat)*pp(imat)
                    else
                      ! Mie-Gruneisen
                      if(eos_name == 'mie-grueneisen')then
                         smallgamma=eos_params(imat,1);biggamma=eos_params(imat,2);p_0=eos_params(imat,3);rho_0=eos_params(imat,4)
                         eta = gg(imat)/rho_0
                         p_c = p_0 * eta**biggamma
                         e_c = p_c / (biggamma-one)
                         delpc = biggamma * p_c

                      ! Cochran-Chan
                      else if(eos_name == 'cochran-chan')then
                         smallgamma=eos_params(imat,1);rho_0=eos_params(imat,2)
                         E_1=eos_params(imat,3);E_2=eos_params(imat,4)
                         A_1=eos_params(imat,5);A_2=eos_params(imat,6)
                         C_v=eos_params(imat,7);T_0=eos_params(imat,8)

                         ! Define the Cochran-Chan constant term
                         E_0 = A_1 / (E_1-one) - A_2 / (E_2-one) + rho_0 * C_v * T_0

                         ! Update Mie-Gruneisen terms for each material
                         eta   = gg(imat)/rho_0
                         p_c_1 = A_1 * eta**E_1
                         p_c_2 = A_2 * eta**E_2
                         p_c   = p_c_1 - p_c_2
                         e_c   = p_c_1 / (E_1-1.0) - p_c_2 / (E_2-1.0) - eta * E_0
                         delpc = p_c_1 * E_1 - p_c_2 * E_2
                      end if
                      pp(imat) = (smallgamma-1)*(gg(imat)*ee(imat)-e_c) + p_c
                      rc2(imat) = delpc + smallgamma * (pp(imat)-p_c)
                      ptot = ptot + ff(imat)*pp(imat)
                      if(rc2(imat)<0.0)then
                        ! write(*,*) "Sound speed",imat,ff(imat),gg(imat),ee(imat),e_c,pp(imat),rc2(imat)
                      end if
                    end if
                 end do

                 if(eos_name .ne. 'stiffened gas')then
                   do imat = 1,nmat
                      smallgamma = eos_params(imat,1)
                      rc2(imat) = rc2(imat) + (smallgamma-1)*(ptot-pp(imat))
                      if(rc2(imat)<0.0)then
                        ! write(*,*) "Correction sound speed", rc2(imat), ptot,pp(imat)
                      end if
                   end do
                 end if

                 ! Compute weights
                 do imat = 1,nmat
                    w(imat) = rc2(imat)/ff(imat)
                 end do

                 ! Compute new volume fraction
#if NMAT==2
                 dd=w(1)+w(2)

                 ff_new(1) = ff(1) + (pp(1)-pp(2))/dd
                 ff_new(2) = ff(2) + (pp(2)-pp(1))/dd

                 if(ff_new(1)<0.0)then
!                    write(*,*) "relaxation: negative volume fraction",ff_new(1),ff_new(2), ff(1), ff(2)
                    ff_new(1) = ff(1) - ff(1)/2.0 !1d-8 !ff(1) + (pp(1)-pp(2))/dd/2.0
                    ff_new(2) = ff(2) + ff(1)/2.0 !1d0 - 1d-8 !ff(2) + (pp(2)-pp(1))/dd/2.0
                 end if

                 if(ff_new(2)<0.0)then
!                    write(*,*) "relaxation: negative volume fraction",ff_new(1),ff_new(2), ff(1), ff(2)
                    ff_new(1) = ff(1) + ff(2)/2.0 !1d0 - 1d-8 !ff(1) + (pp(1)-pp(2))/dd/2.0
                    ff_new(2) = ff(2) - ff(2)/2.0 !1d-8 !ff(2) + (pp(2)-pp(1))/dd/2.0
                 end if
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
#if NMAT==5
                 t123=w(1)*w(2)*w(3)
                 t124=w(1)*w(2)*w(4)
                 t125=w(1)*w(2)*w(5)
                 t134=w(1)*w(3)*w(4)
                 t135=w(1)*w(3)*w(5)
                 t145=w(1)*w(4)*w(5)
                 t234=w(2)*w(3)*w(4)
                 t235=w(2)*w(3)*w(5)
                 t245=w(2)*w(4)*w(5)
                 t345=w(3)*w(4)*w(5)
                 dd=5*w(1)*w(2)*w(3)*w(4)+5*w(1)*w(2)*w(3)*w(5)+5*w(1)*w(2)*w(4)*w(5)+5*w(1)*w(3)*w(4)*w(5)+5*w(2)*w(3)*w(4)*w(5)

                 ff_new(1) = ff(1) + ((t234+t235+t245+2*t345)*(pp(1)-pp(2))+(t234+t235+2*t245+t345)*(pp(1)-pp(3))+(t234+2*t235+t245+t345)*(pp(1)-pp(4))+(2*t234+t235+t245+t345)*(pp(1)-pp(5))+ &
                      &                           (t245-t345)*(pp(2)-pp(3))+            (t235-t345)*(pp(2)-pp(4))+            (t234-t345)*(pp(2)-pp(5))+ &
                      &                           (t235-t245)*(pp(3)-pp(4))+            (t234-t245)*(pp(3)-pp(5))+            (t234-t235)*(pp(4)-pp(5)))/dd
                 ff_new(2) = ff(2) + ((t345+t134+t135+2*t145)*(pp(2)-pp(3))+(t345+t134+2*t135+t145)*(pp(2)-pp(4))+(t345+2*t134+t135+t145)*(pp(2)-pp(5))+(2*t345+t134+t135+t145)*(pp(2)-pp(1))+ &
                      &                           (t135-t145)*(pp(3)-pp(4))+            (t134-t145)*(pp(3)-pp(5))+            (t345-t145)*(pp(3)-pp(1))+ &
                      &                           (t134-t135)*(pp(4)-pp(5))+            (t345-t135)*(pp(4)-pp(1))+            (t345-t134)*(pp(5)-pp(1)))/dd
                 ff_new(3) = ff(3) + ((t145+t245+t124+2*t125)*(pp(3)-pp(4))+(t145+t245+2*t124+t125)*(pp(3)-pp(5))+(t145+2*t245+t124+t125)*(pp(3)-pp(1))+(2*t145+t245+t124+t125)*(pp(3)-pp(2))+ &
                      &                           (t124-t125)*(pp(4)-pp(5))+            (t245-t125)*(pp(4)-pp(1))+            (t145-t125)*(pp(4)-pp(2))+ &
                      &                           (t245-t124)*(pp(5)-pp(1))+            (t145-t124)*(pp(5)-pp(2))+            (t145-t245)*(pp(1)-pp(2)))/dd
                 ff_new(4) = ff(4) + ((t125+t135+t235+2*t123)*(pp(4)-pp(5))+(t125+t135+2*t235+t123)*(pp(4)-pp(1))+(t125+2*t135+t235+t123)*(pp(4)-pp(2))+(2*t125+t135+t235+t123)*(pp(4)-pp(3))+ &
                      &                           (t235-t123)*(pp(5)-pp(1))+            (t135-t123)*(pp(5)-pp(2))+            (t125-t123)*(pp(5)-pp(3))+ &
                      &                           (t135-t235)*(pp(1)-pp(2))+            (t125-t235)*(pp(1)-pp(3))+            (t125-t135)*(pp(2)-pp(3)))/dd
                 ff_new(5) = ff(5) + ((t123+t124+t134+2*t234)*(pp(5)-pp(1))+(t123+t124+2*t134+t234)*(pp(5)-pp(2))+(t123+2*t124+t134+t234)*(pp(5)-pp(3))+(2*t123+t124+t134+t234)*(pp(5)-pp(4))+ &
                      &                           (t134-t234)*(pp(1)-pp(5))+            (t124-t234)*(pp(1)-pp(3))+            (t123-t234)*(pp(1)-pp(4))+ &
                      &                           (t124-t134)*(pp(2)-pp(3))+            (t123-t134)*(pp(2)-pp(4))+            (t123-t124)*(pp(3)-pp(4)))/dd
#endif
                 ! ! Compute new specific internal energy
                 do imat = 1,nmat
                    ee_new(imat) = ee(imat) - ptot * (ff_new(imat) - ff(imat)) / (ff(imat)*gg(imat))
                 end do

                 ! Store new volume fraction
                 do imat = 1,nmat
                    uold(ind_cell(i),imat) = ff_new(imat)
                 end do

!!$                 ! Compute new entropy
!!$                 do imat = 1,nmat
!!$                    s_entry = uold(ind_cell(i),3*nmat+ndim+imat)/uold(ind_cell(i),nmat+imat)
!!$                    g_mat(1) = uold(ind_cell(i),nmat+imat)/uold(ind_cell(i),imat)
!!$                    e_mat(1) = ee_new(imat)*g_mat(1)
!!$                    call eos_s(g_mat,e_mat,s_mat,imat,.false.,1)
!!$                    e_th = s_mat(1)*g_mat(1)**eos_params(imat,1)
!!$                    e_cold = eos_params(imat,3)*(g_mat(1)/eos_params(imat,4))**eos_params(imat,2)
!!$                    if(e_th < 0.001*e_kin .or. e_th < 0.001*e_cold)then
!!$!                       write(*,*)'relaxation: entropy iter=',iter,imat,g_mat(1),e_th,e_cold,s_mat(1),s_entry
!!$                       s_mat(1) = s_entry
!!$                       call eos_s(g_mat,e_mat,s_mat,imat,.true.,1)
!!$                       ee_new(imat) = e_mat(1)/g_mat(1)
!!$                    else
!!$                       uold(ind_cell(i),3*nmat+ndim+imat) = uold(ind_cell(i),nmat+imat)*s_mat(1)
!!$                    endif
!!$                 end do

                 ! Store new partial total energy
                 do imat=1,nmat
                    e_tot = ee_new(imat) + e_kin
                    uold(ind_cell(i),2*nmat+ndim+imat) = uold(ind_cell(i),nmat+imat)*e_tot
                 end do

                 iter=iter+1
                 error=abs(pp(1)-pp(2))/abs(pp(1)+pp(2))

              end do

#if NVAR > NDIM + 3*NMAT
              ! Compute new entropy
              do imat = 1,nmat

                 ! Check that entropy is correct on entry
                 s_entry = unew(ind_cell(i),3*nmat+ndim+imat)/unew(ind_cell(i),nmat+imat)
                 if(s_entry<0)then
                    write(*,*)'end relaxation: negative entropy on entry',imat,s_entry
                 endif

                 s_entry = uold(ind_cell(i),3*nmat+ndim+imat)/uold(ind_cell(i),nmat+imat)
                 g_mat(1) = uold(ind_cell(i),nmat+imat)/uold(ind_cell(i),imat)
                 e_tot = uold(ind_cell(i),2*nmat+ndim+imat)/uold(ind_cell(i),nmat+imat)
                 e_mat(1) = (e_tot-e_kin)*g_mat(1)
                 call eos_s(g_mat,e_mat,s_mat,imat,.false.,1)
                 e_th = s_mat(1)*g_mat(1)**eos_params(imat,1)
                 e_cold = eos_params(imat,3)*(g_mat(1)/eos_params(imat,4))**eos_params(imat,2)
                 if(e_th < 0.001*e_kin .or. e_th < 0.001*e_cold)then
                    ! write(*,*) "Entropy fix: e_th=",e_th,"e_kin=",e_kin,"e_cold=",e_cold
                    s_mat(1) = s_entry
                    call eos_s(g_mat,e_mat,s_mat,imat,.true.,1)
                    ee_new(imat) = e_mat(1)/g_mat(1)
                    e_tot = ee_new(imat) + e_kin
                    uold(ind_cell(i),2*nmat+ndim+imat) = uold(ind_cell(i),nmat+imat)*e_tot
                 else
                    uold(ind_cell(i),3*nmat+ndim+imat) = uold(ind_cell(i),nmat+imat)*s_mat(1)
                 endif

              end do

              ! Check new entropy
              do imat = 1,nmat
                 g_mat(1) = uold(ind_cell(i),nmat+imat)/uold(ind_cell(i),imat)
                 e_tot = uold(ind_cell(i),2*nmat+ndim+imat)/uold(ind_cell(i),nmat+imat)
                 e_mat(1) = (e_tot-e_kin)*g_mat(1)
                 call eos_s(g_mat,e_mat,s_mat,imat,.false.,1)
                 if(s_mat(1)<0)then
                    write(*,*)'end relaxation: negative entropy iter=',iter,imat,g_mat(1),e_mat(1),s_mat(1)
                 endif
              end do
#endif
              if(iter.EQ.iter_max)then
                 write(*,*)'pressure relaxation iter=',iter,ff(1),ff(2),gg(1),gg(2),pp(1),pp(2),error
              endif

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
subroutine phase_transition(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !-------------------------------------------------------------------
  ! Update volume fraction using compressibility source terms
  !-------------------------------------------------------------------
  integer::i,k,ivar,imat,idim,ind,iskip,ncache,igrid,ngrid,iphase,idx,io,trial
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp)::one=1.0_dp,half=0.5_dp, zero=0.0_dp
  integer,parameter::nphase=2
  real(dp),dimension(1:8)::xc
  real(dp)::skip_loc,dx,eps,scale,dx_loc
  integer::ix,iy,iz,nx_loc,iter
  integer::iter_max=200,trial_max=5
  real(dp)::error,error_old,pp_best
  real(dp),dimension(1:nphase)::ff,ee,gg,ff_new,ee_new,p_sat,dedp,dgdp,dfdp,gg_new
  real(dp),dimension(1:nphase)::g_sat
  real(dp),dimension(1:nphase)::gleft,gcen,gright,eleft,ecen,eright
  real(dp),dimension(1:nvector)::gg_mat,ee_mat,cc_mat
  real(dp),dimension(1:nvector)::pp_mat
  real(dp),dimension(1:ndim),save::vv
  real(dp)::pleft,pcen,pright,ftot
  real(dp)::e_kin,e_tot,dtot,ptot,pp_new,f0,e0,p0,g0
  logical,save::read_flag=.false.
  logical::mod_flag=.false.
  integer,parameter::nrows=10000,ncols=5          ! CSV file parameters
  real(dp),dimension(1:nrows, 1:ncols),save::xx   ! Saturation dome values (P, d_v, d_l)
  integer,save::nmax
  real(dp),save::p_crit, d_crit
  real(dp),save::step
  real(dp)::smallgamma,p_0,q,deriv

  ! Ideas for generalising the phase transition routine for an arbitrary combination of
  ! fluids with an arbitrary number of them transitioning:
  ! Here we would need an index list of which phases we are transitioning
  ! The subroutine should be reading in a file with (P_sat12, ... , P_sat1m, P_satl(l+1), ..., P_satln,
  ! d_1, ..., d_m, d_l, ..., d_n, e_1, e_2, ... , e_n) values
  ! Here P_sat12 would be the saturation pressure for the mixture of fluid 1 and fluid 2
  ! d_1 and d_2 the densities at P_sat12 (Analogous for e_1,2)
  ! If phase_params, for example, equals 0 no phase transition considerations
  ! Example phase_params = (/1,2,0/) for NMAT=3, liquid-vapor mixture and another fluid
  ! => Fluid 1,2 are a pair with pre-computed phase transition data available. No phase transition considerations for fluid 3
  ! If P<P_crit for any pair of fluids store their indices in the idx variable and trigger the phase transition routine
  ! After all fluid pairs are relaxed, we perform another pressure relaxation to equalize the pressure before continuing

  ! Read saturation dome values into an array (For now for a single fluid-vapour mixture)
  if (.not. read_flag) then
     xx=0d0
     open (unit=10,file="patch/6eq/psat.csv",action="read",status="old")
     nmax=0
     do
        nmax=nmax+1
        read (10,*,iostat=io) xx(nmax,:)
        if(io.ne.0)exit
     end do
     read_flag = .true.

     ! Critical pressure & density point below which phase transition can happen
     ! File structure (P_sat, d_liquid, d_vapour, e_liquid, e_vapour)
     p_crit = maxval(xx(:,1))
     d_crit = xx(nmax-1,2)
     step   = xx(2,1) - xx(1,1)
     write(*,*)'Saturation dome file read'
     write(*,*)'nmax=',nmax,' P_crit=',p_crit,' d_crit=',d_crit, " step=", step
  end if

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

           ! Only update leaf cells
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
              e_kin=0.0
              do idim = 1,ndim
                 e_kin = e_kin + half*vv(idim)**2
              end do

              ! Calculate the the total pressure
              ptot = 0.0
              do imat=1,nmat
                 ff(imat)   = uold(ind_cell(i),imat)
                 gg_mat(1)  = uold(ind_cell(i),nmat+imat)/ff(imat)
                 e_tot      = uold(ind_cell(i),2*nmat+ndim+imat)/ff(imat)
                 ee_mat(1)  = e_tot - gg_mat(1)*e_kin
                 call eos(gg_mat,ee_mat,pp_mat,cc_mat,imat,.false.,1)
                 ptot = ptot + pp_mat(1)*ff(imat)
              end do

              ! Only check for phase transition below the critical pressure
              if(ptot < p_crit)then

                 mod_flag = .false.

                 ! Compute f0, g0 and e0 for the given mixture
                 idx= 0.0 ! For now
                 f0 = 0.0
                 g0 = 0.0
                 e0 = 0.0
                 do iphase=1,nphase
                    ivar = idx + iphase
                    ff(iphase) = uold(ind_cell(i),ivar)
                    gg(iphase) = uold(ind_cell(i),nmat+ivar)/ff(iphase)
                    e_tot      = uold(ind_cell(i),2*nmat+ndim+ivar)/ff(iphase)
                    ee(iphase) = e_tot - gg(iphase)*e_kin
                    if(ee(iphase)<0)then
                       write(*,*)'phase transition: negative energy',ff(1),ff(2),gg(1),gg(2),ee(iphase),e_tot,e_kin
                    endif
                    f0 = f0 + ff(iphase)
                    g0 = g0 + ff(iphase) * gg(iphase)
                    e0 = e0 + ff(iphase) * ee(iphase)
                 end do
                 g0 = g0/f0
                 e0 = e0/f0

                 ! Set the interpolation kernel
                 k = (ptot/p_crit)*(nmax-1)
                 k = MIN(MAX(k,2),nmax-2)
                 pleft = xx(k-1, 1)
                 pcen  = xx(k  , 1)
                 pright= xx(k+1, 1)

                 do iphase=1,nphase
                    gleft(iphase) = xx(k-1, 1+iphase)
                    gcen(iphase)  = xx(k  , 1+iphase)
                    gright(iphase)= xx(k+1, 1+iphase)
                    eleft(iphase) = xx(k-1, 3+iphase)
                    ecen(iphase)  = xx(k  , 3+iphase)
                    eright(iphase)= xx(k+1, 3+iphase)
                 end do

                 ! Interpolate to the saturation dome values
                 do iphase=1,nphase
                    g_sat(iphase)  = ((ptot - pright)*(ptot - pleft)) /((pcen   - pright)*(pcen   - pleft)) *gcen(iphase)   &
                         &         + ((ptot - pcen)  *(ptot - pleft)) /((pright - pcen)  *(pright - pleft)) *gright(iphase) &
                         &         + ((ptot - pcen)  *(ptot - pright))/((pleft  - pcen)  *(pleft  - pright))*gleft(iphase)
                    ee_new(iphase) = ((ptot - pright)*(ptot - pleft)) /((pcen   - pright)*(pcen   - pleft)) *ecen(iphase)   &
                         &         + ((ptot - pcen)  *(ptot - pleft)) /((pright - pcen)  *(pright - pleft)) *eright(iphase) &
                         &         + ((ptot - pcen)  *(ptot - pright))/((pleft  - pcen)  *(pleft  - pright))*eleft(iphase)
                 end do

                 ! If we are inside the saturation dome, we trigger the mass transfer routine
                 if((gg(1) < g_sat(1) .or. gg(2) > g_sat(2)) .and. g0 < g_sat(1) .and. g0 > g_sat(2))then

                    mod_flag = .true.

                    ! Newton iteration to get the new pressure
                    iter = 0
                    error = 1
                    trial = 0
                    do while(abs(error).GT.1d-12.AND.iter<iter_max)
                       ! Initial guess
                       if(iter==0)then
                          pp_new = ptot !xx(k,1)
                       ! Normal iteration case
                       else
                          deriv  = ff_new(1)*dedp(1)+ff_new(2)*dedp(2) &
                               & + f0*(dgdp(1)*ee_new(2)-dgdp(2)*ee_new(1))/(g_sat(1)-g_sat(2)) &
                               & - (ff_new(1)*ee_new(1)+ff_new(2)*ee_new(2))*(dgdp(1)-dgdp(2))/(g_sat(1)-g_sat(2))
                          pp_new = pp_new - (ff_new(1)*ee_new(1)+ff_new(2)*ee_new(2)-f0*e0)/deriv
                       end if

                       ! Attempt at continuing the iteration when it gets stuck in a loop
                       if(iter == iter_max-1 .and. trial<trial_max)then
                          pp_new = pp_best + pp_best*1e-5
                          iter   = 1
                          trial = trial + 1
                          write(*,*) iter, "Current state: P_n=",pp_new
                          write(*,*) "Retrying..."
                       end if
                       ! Set the interpolation kernel
                       k = (pp_new/p_crit)*(nmax-1)
                       k = MIN(MAX(k,2),nmax-2)
                       pleft = xx(k-1, 1)
                       pcen  = xx(k  , 1)
                       pright= xx(k+1, 1)

                       do iphase=1,nphase
                          gleft(iphase) = xx(k-1, 1+iphase)
                          gcen(iphase)  = xx(k  , 1+iphase)
                          gright(iphase)= xx(k+1, 1+iphase)
                          eleft(iphase) = xx(k-1, 3+iphase)
                          ecen(iphase)  = xx(k  , 3+iphase)
                          eright(iphase)= xx(k+1, 3+iphase)
                       end do

                       ! Interpolate to the saturation dome values
                       do iphase=1,nphase
                          g_sat(iphase)  = ((pp_new - pright)*(pp_new - pleft)) /((pcen   - pright)*(pcen   - pleft)) *gcen(iphase)   &
                               &         + ((pp_new - pcen)  *(pp_new - pleft)) /((pright - pcen)  *(pright - pleft)) *gright(iphase) &
                               &         + ((pp_new - pcen)  *(pp_new - pright))/((pleft  - pcen)  *(pleft  - pright))*gleft(iphase)
                          ee_new(iphase) = ((pp_new - pright)*(pp_new - pleft)) /((pcen   - pright)*(pcen   - pleft)) *ecen(iphase)   &
                               &         + ((pp_new - pcen)  *(pp_new - pleft)) /((pright - pcen)  *(pright - pleft)) *eright(iphase) &
                               &         + ((pp_new - pcen)  *(pp_new - pright))/((pleft  - pcen)  *(pleft  - pright))*eleft(iphase)
                          dgdp  (iphase) = (2*pp_new - (pright+pleft ))/((pcen   - pright)*(pcen   - pleft)) *gcen(iphase)   &
                               &         + (2*pp_new - (pcen  +pleft ))/((pright - pcen)  *(pright - pleft)) *gright(iphase) &
                               &         + (2*pp_new - (pcen  +pright))/((pleft  - pcen)  *(pleft  - pright))*gleft(iphase)
                          dedp  (iphase) = (2*pp_new - (pright+pleft ))/((pcen   - pright)*(pcen   - pleft)) *ecen(iphase)   &
                               &         + (2*pp_new - (pcen  +pleft ))/((pright - pcen)  *(pright - pleft)) *eright(iphase) &
                               &         + (2*pp_new - (pcen  +pright))/((pleft  - pcen)  *(pleft  - pright))*eleft(iphase)
                       end do

                       ! Update volume fractions
                       ff_new(1) = f0 * (g0 - g_sat(2))/(g_sat(1) - g_sat(2))
                       ff_new(2) = f0 * (g_sat(1) - g0)/(g_sat(1) - g_sat(2))

                       ! Update iteration variable
                       iter = iter+1

                       ! Calculate the error
                       error_old = error
                       error  = (ff_new(1)*ee_new(1)+ff_new(2)*ee_new(2)-f0*e0)/(f0*e0)

                       ! Store best guess for possibly restarting the iteration
                       if(abs(error)<abs(error_old) .and. iter>1)then
                          pp_best = pp_new
                       end if
                    end do

                    ! Update the density terms and normalize volume fractions
                    do iphase=1,nphase
                       gg_new(iphase) = g_sat(iphase)
                    end do

                    ! Check for unphysical states landing in the saturation dome when no phase transition is happening (<- Does not do anything so far)
                 else if(gg(1) < g_sat(1) .and. (g0 > g_sat(1) .or. g0 < g_sat(2)))then

                    mod_flag  = .true.

                    ff_new(1) = 1.0 - 1e-2
                    ff_new(2) = 1e-2
                    gg_new(1) = gg(1)
                    gg_new(2) = g_sat(2)
                    ee_new(1) = ee(1)

                 else if(gg(2) > g_sat(2) .and. (g0 > g_sat(1) .or. g0 < g_sat(2)))then

                    mod_flag  = .true.

                    ff_new(1) = 1e-2
                    ff_new(2) = 1.0 - 1e-2
                    gg_new(1) = g_sat(1)
                    gg_new(2) = gg(2)
                    ee_new(2) = ee(2)

                 end if

                 ! If values were modified update existing solutions
                 if(mod_flag)then

                    do iphase = 1,nphase
                       ivar = idx + iphase
                       ! Store new volume fractions & densities
                       uold(ind_cell(i),ivar)      = ff_new(iphase)
                       uold(ind_cell(i),nmat+ivar) = ff_new(iphase) * gg_new(iphase)
                       ! Store new partial total energies
                       e_tot = ee_new(iphase) + gg_new(iphase)*e_kin
                       uold(ind_cell(i),2*nmat+ndim+ivar) = uold(ind_cell(i),ivar)*e_tot
                    end do

                 end if

              end if

           endif

        end do
        ! End loop over octs
     end do
     ! End loop over cells
  end do
  ! End loop over grids
end subroutine phase_transition
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
