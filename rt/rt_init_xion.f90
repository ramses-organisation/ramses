!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE rt_init_xion(ilevel)

! Initialize hydrogen ionization state in all cells at given level from
! density and temperature in the cells, assuming chemical equilibrium.
!-------------------------------------------------------------------------
  use amr_commons
  use hydro_commons
  use mpi_mod
  implicit none
  integer:: ilevel
  integer:: ncache,i,igrid,ngrid
  integer,dimension(1:nvector),save:: ind_grid
!-------------------------------------------------------------------------
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel
  ! Do the initialization by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     call rt_init_xion_vsweep(ind_grid, ngrid)
  end do

111 format('   Entering rt_init_xion for level',i2)

END SUBROUTINE rt_init_xion

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE rt_init_xion_vsweep(ind_grid, ngrid)

! Vector sweep initialization of hydrogen ionization state
! ind_grid => Indexes of grids/octs to initialize
! ngrid    => Number of vaid indexes in ind_grid (i.e. number of grids)
!-------------------------------------------------------------------------
  use amr_commons
  use hydro_commons
  use rt_parameters,only:nIons,iIons,isH2,isHe,ixHI,ixHII,ixHeII,ixHeIII
  use cooling_module,only:Y
  implicit none
  integer::ngrid
  integer,dimension(1:nvector)::ind_grid
  integer::i, ind, iskip, idim, nleaf
  real(dp)::scale_nH, scale_T2, scale_l, scale_d, scale_t, scale_v
  integer,dimension(1:nvector),save::ind_cell, ind_leaf
  real(dp)::nH, T2, ekk, err, emag, x, mu, Zsolar
  real(dp),dimension(nIons)::phI_rates       ! Photoionization rates [s-1]
  real(dp),dimension(7)::nSpec               !          Species abundances
#if NENER>0
  integer::irad
#endif
!-------------------------------------------------------------------------
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  pHI_rates(:)=0.0                   ! No UV background for the time being
  ! Loop over cells in each oct
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,ngrid
        ind_cell(i)=iskip+ind_grid(i)
     end do
     ! Gather leaf cells
     nleaf=0
     do i=1,ngrid
        if(son(ind_cell(i))==0)then
           nleaf=nleaf+1
           ind_leaf(nleaf)=ind_cell(i)
        end if
     end do
     if(nleaf .eq. 0) cycle

     do i=1,nleaf
        ! Compute rho
        nH = MAX(uold(ind_leaf(i),1),smallr)   !       Mass density of gas
        Zsolar = z_ave                         ! Metallicity (solar units)
        if(metal) &
             Zsolar=(uold(ind_leaf(i),imetal)) / nH / 0.02
        ! Compute pressure from energy density
        T2 = uold(ind_leaf(i),neul)          ! Energy density (kin+heat)
        ekk = 0.0d0                            !            Kinetic energy
        do idim=2,neul-1
           ekk = ekk+0.5*uold(ind_leaf(i),idim)**2/nH
        end do
        err = 0.0d0
#if NENER>0
        do irad=0,nener-1
           err = err+uold(ind_leaf(i),inener+irad)
        end do
#endif
        emag = 0.0d0
#ifdef SOLVERmhd
        do idim=1,3
           emag=emag+0.125d0*(uold(ind_leaf(i),idim+neul)+uold(ind_leaf(i),idim+nvar))**2
        end do
#endif
        T2 = (gamma-1.0)*(T2-ekk-err-emag)     !     Gamma is ad. exponent
        ! now T2 is pressure (in user units)   !    (relates p and energy)
        ! Compute T2=T/mu in Kelvin from pressure:
        T2 = T2/nH*scale_T2                    !        Ideal gas equation
        ! Compute nH in H/cc (number of H nuclei per cubic centimeter)
        nH = nH*scale_nH

        call cmp_Equilibrium_Abundances(T2,nH,pHI_rates,mu,nSpec,Zsolar)

        ! UPDATE IONIZATION STATES
        if(isH2) then
           if(nrestart==0 .and. cosmo) then !Prevent primordial H2 in cosmo sims
              x = (2.*nSpec(2)+nSpec(3))/(2.*nSpec(2)+nSpec(3)+nspec(4)) ! HI
           else
              x = nSpec(3)/(2.*nSpec(2)+nSpec(3)+nspec(4)) ! HI
           endif
           uold(ind_leaf(i),iIons-1+ixHI) = x*uold(ind_leaf(i),1)
        endif
        x = nSpec(4)/(2.*nSpec(2)+nSpec(3)+nspec(4))        ! HII fraction
        uold(ind_leaf(i),iIons-1+ixHII) = x*uold(ind_leaf(i),1)
        if(Y .gt. 0d0 .and. isHe) then
           x = nSpec(6)/(nSpec(5)+nSpec(6)+nSpec(7))      !  HeII fraction
           uold(ind_leaf(i),iIons-1+ixHeII) = x*uold(ind_leaf(i),1)
           x = nSpec(7)/(nSpec(5)+nSpec(6)+nSpec(7))      ! HeIII fraction
           uold(ind_leaf(i),iIons-1+ixHeIII) = x*uold(ind_leaf(i),1)
        endif
      end do

  end do
  ! End loop over cells
END SUBROUTINE rt_init_xion_vsweep

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE calc_equilibrium_xion(vars, rtvars, xion)

! Calculate and return photoionization equilibrium abundance states for
! a cell
! vars     => Cell variables (rho, v, u, w, etc)
! rtvars   => Cell RT variables (Np1, Fpx1, Fpy1, etc)
! xion     => Equilibrium ionization states of cell
!-------------------------------------------------------------------------
  use amr_commons
  use hydro_commons
  use rt_parameters
  use cooling_module,only:Y
  use rt_cooling_module,only:UVrates,signc
  implicit none
  real(dp),dimension(nvar_all)::vars
  real(dp),dimension(nrtvar)::rtvars
  real(dp),dimension(nIons)::xion
  integer::ip, iI, idim
  real(dp)::scale_nH, scale_T2, scale_l, scale_d, scale_t, scale_v
  real(dp)::scale_Np,scale_Fp,nH,T2,ekk,err,emag,mu,Zsolar,ss_factor
  real(dp),dimension(nIons)::phI_rates       ! Photoionization rates [s-1]
  real(dp),dimension(7)::nSpec               !          Species abundances
#if NENER>0
  integer::irad
#endif
!-------------------------------------------------------------------------
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  call rt_units(scale_Np,scale_Fp)

  ! Calculate photoionization rates:
  phI_rates(:)=0.0
  do ip=1, nGroups
     do iI=1,nIons
        phI_rates(iI) = phI_rates(iI) &
                      + rtvars(iGroups(ip))*scale_Np*signc(ip,iI)
     end do
  end do

  nH = MAX(vars(1),smallr)                  !   Number density of gas [UU]

  Zsolar = z_ave
  if(metal) Zsolar=vars(imetal) / nH / 0.02 !    Metallicity (solar units)

  ! Compute pressure from energy density
  T2 = vars(neul)                         ! Energy dens. (kin+heat) [UU]
  ekk = 0.0d0                             !          Kinetic energy [UU]
  do idim=2,neul-1
     ekk = ekk+0.5*vars(idim)**2/nH
  end do
  err = 0.0d0
#if NENER>0
  do irad=0,nener-1
     err = err+vars(inener+irad)
  end do
#endif
  emag = 0.0d0
#ifdef SOLVERmhd
  do idim=1,3
     emag=emag+0.125d0*(vars(idim+neul)+vars(idim+nvar))**2
  end do
#endif
  T2 = (gamma-1.0)*(T2-ekk-err-emag)        !        Gamma is ad. exponent
                                            !      now T2 is pressure [UU]
  T2 = T2/nH*scale_T2                       !                T/mu [Kelvin]
  nH = nH*scale_nH                          !        Number density [H/cc]

  ! UV background photoionization:
  ss_factor = 1d0
  if(self_shielding) ss_factor = exp(-nH/1d-2)
  if(haardt_madau) phI_rates = phI_rates + UVrates(:,1) * ss_factor

  call cmp_Equilibrium_Abundances(T2, nH, pHI_rates, mu, nSpec, Zsolar)

  if(isH2) xion(ixHI)=nSpec(3)/(2.*nSpec(2)+nSpec(3)+nSpec(4))!    HI frac
  xion(ixHII)=nSpec(4)/(2.*nSpec(2)+nSpec(3)+nSpec(4))        !   HII frac
  if(Y .gt. 0d0 .and. isHe) then
     xion(ixHeII) = nSpec(6)/(nSpec(5)+nSpec(6)+nSpec(7)) !  HeII fraction
     xion(ixHeIII) = nSpec(7)/(nSpec(5)+nSpec(6)+nSpec(7))! HeIII fraction

  endif

END SUBROUTINE calc_equilibrium_xion

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE cmp_Equilibrium_Abundances(T2,nH,phI_rates,mu,nSpec,Zsolar)
!-------------------------------------------------------------------------
  use amr_commons,only:dp
  use rt_cooling_module
  use rt_parameters,only:nIons
  implicit none
  real(dp) ::T2,nH
  real(dp),dimension(nIons)::phI_rates
  real(dp) ::mu,Zsolar
  real(dp),dimension(1:7)::nSpec!-----------------------------------------
  real(dp) ::mu_mid, err_mu, mu_left, mu_right, T, nTot
  integer :: niter
!-------------------------------------------------------------------------
  ! Iteration to find mu                     ! n_E     = n_spec(1) ! e
                                             ! n_H2    = n_spec(2) ! H2
  err_mu=1.                                  ! n_HI    = n_spec(3) ! H
  mu_left=0.5                                ! n_HII   = n_spec(4) ! H+
  mu_right=2.3                               ! n_HEI   = n_spec(5) ! He
  niter=0                                    ! n_HEII  = n_spec(6) ! He+
                                             ! n_HEIII = n_spec(7) ! He++
  do while (err_mu > 1d-4 .and. niter <= 50)
     mu_mid=0.5*(mu_left+mu_right)
     T = T2*mu_mid
     call cmp_chem_eq(T, nH, phI_rates, nSpec, nTot, mu, Zsolar)
     if(mu - mu_mid > 0.)then
        mu_left = mu_mid
     else
        mu_right= mu_mid
     end if
     err_mu = ABS(mu_right-mu_left)/mu_mid
     niter=niter+1
  end do
  mu=mu_mid
  if (niter > 50) then
     write(*,*) 'ERROR in cmp_Equilibrium_Abundances : too many iterations.'
     STOP
  endif

END SUBROUTINE cmp_Equilibrium_Abundances

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE test_cmp_Equilibrium_Abundances()
!-------------------------------------------------------------------------
  use amr_commons,only:dp,z_ave
  use rt_parameters,only:nIons
  implicit none
  integer,parameter::n=500
  real(dp)::T2(n), T20, T21
  real(dp),dimension(nIons)::phI_rates
  real(dp) ::mu
  real(dp),dimension(1:7)::nSpec
  real(dp)::nH, xH2, xHI, xHII, nHe, xHeI, xHeII, xHeIII
  integer::i

! Test cmp_Equilibrium_Abundances routine and write the output to a file
! Run for one gas density with no radiation field and zero metallicity
! (this tends to be the setup most prone to non-convergence fot the
! iterative calculation, due to the weak H2 formation rate at low T)
!-------------------------------------------------------------------------
  T20 = 1d0
  T21 = 6d0
  nH=1e0 ! 1 atom per cm3
  phI_rates(:) = 0d0
  ! Set up evently log-spaced temperatures from 10 to 10^7 K
  do i = 1, n
     T2(i) = T20 + (T21 - T20) * (i - 1) / (n - 1)
  end do
  T2 = 10.d0**T2

  open(unit=10,file='eq-abundances.txt',form='formatted')

  ! Calculate equilbrium ionisation states for all temperatures
  do i = 1, n
    call cmp_Equilibrium_Abundances(T2(i), nH, pHI_rates, mu, nSpec, z_ave)
    xH2    = 2.*nSpec(2)/nH
    xHI    = nSpec(3)/nH
    xHII   = nSpec(4)/nH
    nHe    = nSpec(5)+nSpec(6)+nSpec(7)
    xHeI   = nSpec(5)/nHe
    xHeII  = nSpec(6)/nHe
    xHeIII = nSpec(7)/nHe
    write(10,'(8(ES15.5E3,1x))') T2(i), xH2, xHI, xHII, xHeI, xHeII, xHeIII
  end do
  close(10)

END SUBROUTINE test_cmp_Equilibrium_Abundances
