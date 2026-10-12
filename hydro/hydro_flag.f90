subroutine hydro_flag(ilevel)
  use amr_commons
  use hydro_commons
#ifdef RT
  use rt_parameters
#endif
  implicit none
  integer::ilevel
  ! -------------------------------------------------------------------
  ! This routine flag for refinement cells that satisfies
  ! some user-defined physical criteria at the level ilevel.
  ! -------------------------------------------------------------------
  integer::i,j,ncache,nok,ix,iy,iz,iskip
  integer::igrid,ind,idim,ngrid,ivar
  integer::nx_loc
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  integer,dimension(1:nvector,0:twondim),save::igridn
  integer,dimension(1:nvector,1:twondim),save::indn

  logical,dimension(1:nvector),save::ok

  real(dp)::dx,dx_loc,scale
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp),dimension(1:nvector,1:nvar_all),save::uug,uum,uud

  if(ilevel==nlevelmax)return
  if(numbtot(1,ilevel)==0)return

  ! Rescaling factors
  dx=0.5d0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  if(    .not. neq_chem  .and.&
       & err_grad_d==-1.0.and.&
       & err_grad_p==-1.0.and.&
       & err_grad_u==-1.0.and.&
#ifdef SOLVERmhd
       & err_grad_A==-1.0.and.&
       & err_grad_B==-1.0.and.&
       & err_grad_C==-1.0.and.&
       & err_grad_B2==-1.0.and.&
#endif
       & jeans_refine(ilevel)==-1.0 )return

#ifdef RT
  if( aexp .lt. rt_refine_aexp) return
  if(    neq_chem        .and.          &
       & err_grad_d==-1.0.and.          &
       & err_grad_p==-1.0.and.          &
       & err_grad_u==-1.0.and.          &
       & jeans_refine(ilevel)==-1.0.and.&
       & rt_err_grad_xHII==-1.0 .and.   &
       & rt_err_grad_xHI==-1.0          &
  & ) &
      return
#endif

  ! Loop over active grids
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector

     ! Gather nvector grids
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do

     ! Gather neighboring offsets
     call getnborgrids(ind_grid,igridn,ngrid)

     ! Loop over cells
     do ind=1,twotondim

        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do

        ! Initialize refinement to false
        do i=1,ngrid
           ok(i)=.false.
        end do

        ! Gather neighboring cells
        call getnborcells(igridn,ind,indn,ngrid)

        ! If a neighbor cell does not exist,
        ! replace it by its father cell
        do j=1,twondim
           do i=1,ngrid
              if(indn(i,j)==0)then
                 indn(i,j)=nbor(ind_grid(i),j)
              end if
           end do
        end do

        ! Loop over dimensions
        do idim=1,ndim
           ! Gather hydro variables
           do ivar=1,nvar_all
              do i=1,ngrid
                 uug(i,ivar)=uold(indn(i,2*idim-1),ivar)
                 uum(i,ivar)=uold(ind_cell(i     ),ivar)
                 uud(i,ivar)=uold(indn(i,2*idim  ),ivar)
              end do
           end do
#ifdef SOLVERmhd
           call hydro_refine(uug,uum,uud,ok,ngrid,ilevel)
#else
           call hydro_refine(uug,uum,uud,ok,ngrid)
#endif
        end do

        if(poisson.and.jeans_refine(ilevel)>0.0)then
           if (collapse_jeans_refine) then
              call collapse_jeans_length_refine(ind_cell,ok,ngrid,ilevel)
           else
              call jeans_length_refine(ind_cell,ok,ngrid,ilevel)
           endif
        endif

        ! Apply geometry-based refinement criteria
        if(r_refine(ilevel)>-1.0)then
           ! Compute cell center in code units
           do idim=1,ndim
              do i=1,ngrid
                 xx(i,idim)=xg(ind_grid(i),idim)+xc(ind,idim)
              end do
           end do
           ! Rescale position from code units to user units
           do idim=1,ndim
              do i=1,ngrid
                 xx(i,idim)=(xx(i,idim)-skip_loc(idim))*scale
              end do
           end do
           call geometry_refine(xx,ok,ngrid,ilevel)
        end if

        ! Count newly flagged cells
        nok=0
        do i=1,ngrid
           if(flag1(ind_cell(i))==0.and.ok(i))then
              nok=nok+1
           end if
        end do

        do i=1,ngrid
           if(ok(i))flag1(ind_cell(i))=1
        end do

        nflag=nflag+nok
     end do
     ! End loop over cells

  end do
  ! End loop over grids

end subroutine hydro_flag
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine jeans_length_refine(ind_cell,ok,ncell,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use constants, only: pi
  implicit none
  integer,intent(in)::ncell,ilevel
  integer,dimension(1:nvector),intent(in)::ind_cell
  logical,dimension(1:nvector),intent(inout)::ok
  !-------------------------------------------------------------------
  ! This routine flags cells for refinement if the Jeans length is
  ! resolved by less than the user-specified number of cells.
  ! Input:
  !   ind_cell: cell indices of the current vector sweep
  ! Output:
  !   ok: updated refinement flag for these cells
  !-------------------------------------------------------------------
  integer::nx_loc,i
  real(dp)::dx,scale,dx_loc,lamb_jeans,n_jeans
  real(dp)::cs2,factG
  real(dp),dimension(1:nvector),save::rho,ekin,erad,emag,etherm

  factG=1
  if(cosmo)factG=3d0/8d0/pi*omega_m*aexp
  n_jeans = jeans_refine(ilevel)

  ! Mesh spacing in this level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  ! Gather density
  do i=1,ncell
     rho(i) = max(uold(ind_cell(i),1),smallr)
  enddo

  ! Compute thermal energy
  call cmp_energy_components(ind_cell,ncell,rho,ekin,erad,emag)
  do i=1,ncell
     etherm(i) = uold(ind_cell(i),neul) - ekin(i) - erad(i) - emag(i)
  enddo

  do i=1,ncell
     ! compute isothermal sound speed (squared)
     cs2 =  etherm(i) / rho(i) * (gamma - 1.0d0)
     ! prevent numerical crash due to negative temperature
     cs2 = max(cs2,smallc**2)
     ! compute the Jeans length (remember G=1)
     lamb_jeans = sqrt( cs2 * pi / (rho(i) * factG) )
     ! the Jeans length must be smaller than n_jeans times the size of the cell
     ok(i) = ok(i) .or. ( n_jeans*dx_loc >= lamb_jeans )
  end do

end subroutine jeans_length_refine
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine collapse_jeans_length_refine(ind_cell,ok,ncell,ilevel)
  use amr_commons
  use hydro_commons
  use constants, only: pi, kB, mH
  implicit none
  integer,intent(in)::ncell,ilevel
  integer,dimension(1:nvector),intent(in)::ind_cell
  logical,dimension(1:nvector),intent(inout)::ok
  !-------------------------------------------------------------------
  ! Jeans-length refinement with a prescribed thermal behaviour: the
  ! sound speed used in the Jeans length follows the floor temperature
  ! collapse_jeans_Tfloor at low density and the actual (adiabatic) gas sound speed
  ! at high density, blending log-linearly in density between collapse_jeans_rho_low
  ! and collapse_jeans_rho_high. See jeans_length_refine for the plain version.
  ! Input:
  !   ind_cell: cell indices of the current vector sweep
  ! Output:
  !   ok: updated refinement flag for these cells
  !-------------------------------------------------------------------
  integer::nx_loc,i
  real(dp)::dx,scale,dx_loc,lamb_jeans,n_jeans,factG
  real(dp)::cs2,cs2_gas,cs2_floor,delta_logrho,interpol,d_iso,d_star
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp),dimension(1:nvector),save::rho,ekin,erad,emag,etherm

  factG=1
  if(cosmo)factG=3d0/8d0/pi*omega_m*aexp
  n_jeans = jeans_refine(ilevel)

  ! Mesh spacing in this level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  ! Convert the density thresholds
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  d_iso     = collapse_jeans_rho_low /scale_d
  d_star    = collapse_jeans_rho_high/scale_d

  ! Compute floor sound speed to code units  (adiabatic cs^2 at collapse_jeans_Tfloor)
  cs2_floor = gamma*kB*collapse_jeans_Tfloor/(mu_gas*mH) / scale_v**2

  ! Gather density
  do i=1,ncell
     rho(i) = max(uold(ind_cell(i),1),smallr)
  enddo

  ! Compute thermal energy
  call cmp_energy_components(ind_cell,ncell,rho,ekin,erad,emag)
  do i=1,ncell
     etherm(i) = uold(ind_cell(i),neul) - ekin(i) - erad(i) - emag(i)
  enddo

  do i=1,ncell
     ! actual (adiabatic) sound speed squared
     cs2_gas = gamma*(gamma-1.0d0) * etherm(i)/rho(i)
     ! prevent numerical crash due to negative temperature
     cs2_gas = max(cs2_gas,smallc**2)

     ! use coldest conditions for setting effective sound speed for the Jeans criterion
     ! This makes refinement more aggressive during the 1st and 2nd protostellar collapse
     if(rho(i) <= d_iso)then
        ! Before 2nd collapse: take coldest temperature (as if it was isothermal at collapse_jeans_Tfloor)
        cs2 = min(cs2_gas, cs2_floor)
     else if(rho(i) >= d_star)then
        ! When protostar has formed: take actual gas temperature
        cs2 = cs2_gas
     else
        ! During the second collapse: transition regime
        ! blending log-linearly in density between collapse_jeans_rho_low and collapse_jeans_rho_high.
        delta_logrho = (log10(rho(i)) - log10(d_iso)) / (log10(d_star) - log10(d_iso))
        ! Interpolate between the two above regimes
        interpol = (1-delta_logrho)*log10(min(cs2_gas, cs2_floor)) + delta_logrho*log10(cs2_gas)
        cs2 = 10.0d0**(interpol)
     endif

     ! compute the Jeans length (remember G=1)
     lamb_jeans = sqrt( cs2 * pi / (rho(i) * factG) )
     ! the Jeans length must be smaller than n_jeans times the size of the cell
     ok(i) = ok(i) .or. ( n_jeans*dx_loc >= lamb_jeans )
  end do

end subroutine collapse_jeans_length_refine
