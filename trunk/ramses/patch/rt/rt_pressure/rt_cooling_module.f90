! RT pressure patch:
! Photon momentum is put into gas momentum. 
! =If rt_isIR=.true., then infrared multiscattering radiation is included.
!  The IR cross sections are set with rt_kIR_sc and rt_kIR_abs, but
!  unless the dust temperature evolution is to be included, only use 
!  rt_kIR_sc.
! =If rt_k_IRabs > 0, then dust cooling/heating by IR is included.
! =If rt_kIR_T is set, the IR cross sections are dependant on T
! =If rt_Trad is set, we use T=Tr, where Tr is the radiation temperature.
! =If rt_isoPress=.true., then a cN is used for the radiation pressure
!  of groups other than IR, to include unresolved flux around sources. 
!  Otherwise F is used.
! =If rt_isOpt=.true., then optical radiation is included, which is
!  reprocessed into IR (using cross section rt_kOpt)
! =UV reprocessing into IR can also be included with rt_kUV.
! =rt_isIRtrap=.true. activates trapping of IR radiation in optically
!  thick regions. This is useful if those regions are badly (or not)
!  resolved.
! ------------------------------------------------------------------------

module rt_cooling_module
  use amr_commons,only:myid  
  use cooling_module,only:X, Y
  use rt_parameters
  use coolrates_module
  implicit none

  private   ! default

  public rt_set_model, rt_solve_cooling, update_UVrates, cmp_chem_eq     &
         , isHe, X, Y, rhoc, kB, mH, T2_min_fix, twopi                   &
         , signc, sigec, PHrate, UVrates, rt_isIR, kappaAbs, kappaSc     &
         , is_kIR_T, iIR, rt_isIRtrap, iIRtrapVar, rt_pressBoost         &
         , rt_isoPress, rt_T_rad, rt_vc, a_r

  ! NOTE: T2=T/mu
  ! Np = photon density, Fp = photon flux, 

  real(dp),parameter::rhoc      = 1.88000d-29    !  Crit. density [g cm-3]
  real(dp),parameter::mH        = 1.66000d-24    !         H atom mass [g]
  real(dp),parameter::kB        = 1.38062d-16    ! Boltzm.const. [erg K-1]
  real(dp),parameter::a_r       = 7.5657d-15   ! Rad.const. [erg cm-3 K-4]
  real(dp),parameter::mu_mol    = 1.2195D0
  real(dp),parameter::T2_min_fix=1.d-2           !     Min temperature [K]
  real(dp),parameter::twopi     = 6.2831853d0    !            Two times pi

  real(dp)::T_min, T_frac, x_min, x_frac, Np_min, Np_frac, Fp_min, Fp_frac

  integer,parameter::iIR=1                           !      IR group index
  integer::iIRtrapVar=1                              ! IRtrap pscalar index
  ! Namelist parameters:
  logical::isHe=.true.
  logical::rt_isoPress=.false.         ! Use cE, not F, for rad. pressure
  real(dp)::rt_pressBoost=1d0          ! Boost on RT pressure            
  logical::rt_isIR=.false.             ! Using IR scattering on dust?    
  logical::rt_isIRtrap=.false.         ! IR trapping in passive scalar?  
  logical::is_kIR_T=.false.            ! k_IR propto T^2?               
  logical::rt_T_rad=.false.            ! Use T_gas = T_rad
  logical::rt_vc=.false.               ! (semi-) relativistic RT
  real(dp),dimension(nGroups)::kappaAbs! Dust absorption opacity    
  real(dp),dimension(nGroups)::kappaSc ! Dust scattering opacity    
  
  ! Cooling constants, updated on SED and c-change [cm3 s-1],[erg cm3 s-1]
  real(dp),dimension(nGroups,nIons)::signc,sigec,PHrate

  real(dp),dimension(nIons, 2)::UVrates     !UV backgr. heating/ion. rates

CONTAINS 

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE rt_set_model(Nmodel, J0in_in, J0min_in, alpha_in, normfacJ0_in, &
     zreioniz_in, correct_cooling, realistic_ne, h, omegab, omega0,      &
     omegaL, astart_sim, T2_sim)
! Initialize cooling. All these parameters are unused at the moment and
! are only there for the original cooling-module.
! Nmodel(integer)     =>     Model for UV background and metals
! J0in_in  (dble)     => Default UV intensity
! J0min_in (dble)     => Minimum UV intensity
! alpha_in (dble)     => Slope of the UV spectrum
! zreioniz_in (dble)  => Reionization redshift
! normfacJ0_in (dble) => Normalization factor fot a Harrdt&Madau UV model
! correct_cooling (integer) => Cooling correction
! realistic_ne (integer) => Use realistic electron density at high z?
! h (dble)            => H0/100
! omegab (dble)       => Omega baryons
! omega0 (dble)       => Omega materal total
! omegaL (dble)       => Omega Lambda
! astart_sim (dble)   => Redshift at which we start the simulation
! T2_sim (dble)      <=  Starting temperature in simulation?
!-------------------------------------------------------------------------
  use UV_module
  real(kind=8) :: J0in_in, zreioniz_in, J0min_in, alpha_in, normfacJ0_in
  real(kind=8) :: astart_sim, T2_sim, h, omegab, omega0, omegaL
  integer  :: Nmodel, correct_cooling, realistic_ne, ig
  real(kind=8) :: astart=0.0001, aend, dasura, T2end, mu, ne
!-------------------------------------------------------------------------
  if(myid==1) write(*,*) &
       '==================RT momentum pressure is turned ON=============='
  if(myid==1 .and. rt_isIRtrap) write(*,*) &
       '=========IR trapping is turned ON=============='
  ! do initialization
  isHe=.true. ; if(Y .le. 0.) isHe=.false.
  T_MIN           = 0.1                  !                      Minimum T2
  T_FRAC          = 0.1            

  x_MIN           = 1.d-6                !    Minimum ionization fractions
  x_FRAC          = 0.1    

  Np_MIN = 1.d-13                        !            Photon density floor
  Np_FRAC = 0.2    

  Fp_MIN  = 1D-13*rt_c_cgs               !           Minimum photon fluxes
  Fp_FRAC = 0.5 !0.2                          !           Fp update restriction    

  if(myid==1) write(*,*) 'IR group index has been set to ',iIR        

  ! Might also put in here filling in of tables of cooling rates, to 
  ! ease the computational load.

  ! Calculate initial temperature
  if (astart_sim < astart) then
     write(*,*) 'ERROR in set_model : astart_sim is too small.'
     write(*,*) 'astart     =',astart
     write(*,*) 'astart_sim =',astart_sim
     STOP
  endif
  aend=astart_sim
  dasura=0.02d0

  call update_rt_c
  call init_UV_background

  if(nrestart==0 .and. cosmo)                                            &
       call rt_evol_single_cell(astart,aend,dasura,h,omegab,omega0       &
                               ,omegaL,-1.0d0,T2end,mu,ne,.false.)
  T2_sim=T2end

END SUBROUTINE rt_set_model

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE update_UVrates
! Set the UV ionization and heating rates according to the given a_exp.
!------------------------------------------------------------------------
  use UV_module
  use amr_parameters,only:aexp
  integer::i
!------------------------------------------------------------------------
  UVrates=0.
  if(.not. rt_UV_hom) RETURN
  
  call inp_UV_rates_table(1./aexp - 1., UVrates)
  if(myid==1) write(*,*) 'UV heating: ',UVrates(:,2)

END SUBROUTINE update_UVrates

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE rt_solve_cooling(T2, xion, Np, Fp, p_gas, dNpdt, dFpdt        &
                           ,nH, c_switch, Zsolar, dt, a_exp, nCell)
! Semi-implicitly solve for new temperature, ionization states, 
! photon density/flux, and gas velocity in a number of cells.
! Parameters: 
! T2     <=> T/mu [K] 
! xion   <=> NION ionization fractions 
! Np     <=> NGROUPS photon number densities [cm-3]
! Fp     <=> NGROUPS * ndim photon number fluxes [cm-2 s-1]
! p_gas  <=> ndim gas momentum densities [cm s-1 g cm-3]
! dNpdt   =>  Op split increment in photon densities during dt
! dFpdt   =>  Op split increment in photon flux magnitudes during dt
! nH      =>  Hydrogen number densities [cm-3]  
! c_switch=>  Cooling switch (1 for cool/heat, 0 for no cool/heat)
! Zsolar  =>  Cell metallicities [solar fraction]
! dt      =>  Timestep size             [s]
! a_exp   =>  Cosmic expansion
! nCell   =>  Number of cells (length of all the above vectors)
!
! We use a slightly modified method of Anninos et al. (1997).
!-------------------------------------------------------------------------
  use amr_commons
  implicit none  
  real(dp),dimension(1:nvector):: T2
  real(dp),dimension(1:nvector, nIons):: xion
  real(dp),dimension(1:nvector, nGroups):: Np, dNpdt
  real(dp),dimension(1:nvector, nGroups, ndim):: Fp, dFpdt
  real(dp),dimension(1:nvector, ndim):: p_gas
  real(dp),dimension(1:nvector):: nH, Zsolar
  logical,dimension(1:nvector):: c_switch
  real(dp)::dt, a_exp
  integer::ncell !--------------------------------------------------------
  real(dp),dimension(1:nvector):: tLeft, ddt
  logical:: dt_ok
  real(dp)::dt_rec
  real(dp):: dT2
  real(dp),dimension(nIons):: dXion
  real(dp),dimension(nGroups):: dNp
  real(dp),dimension(nGroups, ndim):: dFp
  real(dp),dimension(ndim):: dp_gas
  integer::i, ia, ig,  nAct, nAct_next, loopcnt, code
  integer,dimension(1:nvector):: indAct              ! Active cell indexes
!-------------------------------------------------------------------------
  tleft(1:ncell) = dt                !       Time left in dt for each cell
  ddt(1:ncell) = dt                  ! First guess at sub-timestep lengths

  do i=1,ncell
     indact(i) = i                   !      Set up indexes of active cells
     ! Ensure all state vars are legal:
     T2(i) = MAX(T2(i), T2_min_fix)
     xion(i,1:nIons) = MIN(MAX(xion(i,1:nIons), x_MIN),1.d0)
     if(xion(i,2)+xion(i,3) .gt. 1.d0) then
        if(xion(i,2) .gt. xion(i,3)) then
           xion(i,2)=1.d0-xion(i,3)
        else
           xion(i,3)=1.d0-xion(i,2)
        endif
     endif
     if(rt) then
        do ig=1,ngroups
           Np(i,ig) = MAX(smallNp, Np(i,ig))
           call reduce_flux(Fp(i,ig,:),Np(i,ig)*rt_c_cgs)
        end do
     endif
  end do

  ! Loop until all cells have tleft=0
  ! **********************************************
  nAct=nCell                                      ! Currently active cells
  loopcnt=0 ; n_cool_cells=n_cool_cells+nCell     !             Statistics
  do while (nAct .gt. 0)      ! Iterate while there are still active cells
     loopcnt=loopcnt+1   ;   tot_cool_loopcnt=tot_cool_loopcnt+nAct 
     nAct_next=0                     ! Active cells for the next iteration
     do ia=1,nAct                             ! Loop over the active cells
        i = indAct(ia)                        !                 Cell index
        call cool_step(T2(i), xion(i,:), Np(i,:), Fp(i,:,:), p_gas(i,:)  &
                      ,dT2,   dXion(:),  dNp(:),  dFp(:,:),  dp_gas(:)   &
                      ,dNpdt(i,:), dFpdt(i,:,:), nH(i), c_switch(i)      &
                      ,Zsolar(i),  ddt(i), a_exp, dt_ok, dt_rec, code    )
        if(loopcnt .gt. 10000) then
           call display_coolinfo(.true., loopcnt, i, dt-tleft(i), dt     &
                            ,ddt(i), nH(i), T2(i),  xion(i,:),  Np(i,:)  &
                            ,Fp(i,:,:),  p_gas(i,:)                      &
                            ,dT2, dXion, dNp, dFp, dp_gas, code)
        endif
        !if(myid==1 .and. i==1) then 
        !   print*,dt_ok,code,T2(i),Np(i,1),Fp(i,1,:)
        !endif
        if(.not. dt_ok) then
           ddt(i)=ddt(i)/2.                    ! Try again with smaller dt 
           nAct_next=nAct_next+1 ; indAct(nAct_next) = i
           loopCodes(code) = loopCodes(code)+1
           cycle 
        endif
        ! Update the cell state (advance the time by ddt):
        T2(i)     = T2(i) + dT2
        xion(i,:) = xion(i,:) + dXion(:)
        if(nGroups .gt. 0) then 
           Np(i,:)   = Np(i,:) + dNp(:)
           Fp(i,:,:) = Fp(i,:,:) + dFp(:,:)
        endif
        p_gas(i,:)   = p_gas(i,:) + dp_gas(:)

        tleft(i)=tleft(i)-ddt(i)
        if(tleft(i) .gt. 0.) then           ! Not finished with this cell
           nAct_next=nAct_next+1 ; indAct(nAct_next) = i
        else if(tleft(i) .lt. 0.) then        ! Overshot by abs(tleft(i))
           print*,'In rt_solve_cooling: tleft < 0  !!'
           stop
        endif
        ddt(i)=min(dt_rec,tleft(i))    ! Use recommended dt from cool_step
     end do ! end loop over active cells
     nAct=nAct_next
  end do ! end iterative loop
  ! loop statistics
  max_cool_loopcnt=max(max_cool_loopcnt,loopcnt)
END SUBROUTINE rt_solve_cooling


!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE cool_step(T2, xion, Np, Fp, p_gas, dT2, dXion, dNp, dFp       &
                    ,dp_gas, dNpdt, dFpdt, nH, c_switch, Zsolar, dt      &
                    ,a_exp, dt_ok, dt_rec, code)
! Compute change in cell state in timestep dt, or return a recommendation 
! for new timestep-length if dt proves too large.
! T2      => T/mu [K]                                 --- dT2 is new value
! xion    => NION ionization fractions                --- dXion is new
! Np      => NGROUPS photon number densities [cm-3]   --- dNp is new value
! Fp      => NGROUPS * ndim photon fluxes [cm-2 s-1]  --- dFp is new value 
! p_gas   => ndim gas momenta [cm s-1 g cm-3]         --- dp_gas is new
! dNpdt   =>  Op split increment in photon densities during dt
! dFpdt   =>  Op split increment in photon flux magnitudes during dt
! nH      =>  Hydrogen number densities [cm-3]  
! c_switch=>  Cooling switch (1 for cool/heat, 0 for no cool/heat)
! Zsolar  =>  Cell metallicities [solar fraction]
! dt      =>  Timestep size [s]
! a_exp   =>  Cosmic expansion
! dt_ok   <=  .f. if timestep constraints were broken, .t. otherwise
! dt_rec  <=  Recommended timesteps for next iteration
! code    <= Error code in cool step, if dt_ok=.f.
!
! The original values, T2, xion etc, must stay unchanged, while dT2, xion
! etc contain the new values (the difference at the end of the routine).
!-------------------------------------------------------------------------
  use UV_module, ONLY: iUVvars_cool, UV_Nphot_cgs
  use amr_commons
  use const
  implicit none  
  real(dp):: T2, dT2
  real(dp),dimension(nIons):: xion, dXion
  real(dp),dimension(nGroups):: Np, dNp, dNpdt
  real(dp),dimension(nGroups, nDim):: Fp, dFp, dFpdt
  real(dp),dimension(nDim):: p_gas, dp_gas, dmom
  real(dp):: nH, Zsolar, dt, a_exp, dt_rec
  real(dp),dimension(nDim):: u_gas ! Gas velocity
  logical::dt_ok, c_switch!-----------------------------------------------
  real(dp),dimension(nIons),save:: alpha, beta, nN, nI
  real(dp):: dUU, fracMax
  real(dp):: xHeI, mu, TK, nHe, ne, neInit, Hrate, dAlpha, dBeta, s, jac,q
  real(dp):: Crate, dCdT2, X_nHkb, rate, dRate, cr, de, photoRate
  real(dp)::metal_tot,metal_prime
  integer::i,j,code
  real(dp),dimension(nGroups),save:: recRad, phAbs, phSc, dustAbs, dustSc        
  real(dp),dimension(nGroups),save:: kAbs_loc,kSc_loc
  real(dp):: rho
  real(dp):: TR, C_v, E_rad, dE_T, fluxMag
!-------------------------------------------------------------------------
  dt_ok=.false.
  nHe=0.25*nH*Y/X  !         Helium number density
  ! Insert UV background emitted and propagated from void regions:
  if(rt_isDiffuseUVsrc .and. nH .le. rt_UVsrc_nHmax)                     &
                Np(iUVvars_cool-4) = max(Np(iUVvars_cool-4), UV_Nphot_cgs)
  ! U contains the original values, dU the updated ones:
  dT2=T2 ; dXion=xion ; dNp=Np ; dFp=Fp ; dp_gas=p_gas
  ! xHI = MAX(1.-xion(1),0.) ; xHII = xion(1)
  ! xHeII=xion(2) ; xHeIII=xion(3)
  xHeI=MAX(1.-xion(2)-xion(3),0.d0)
  ! nN='neutral' species (pre-ionized), nI=their ionized counterparts
  nN(1)  = nH * (1.d0-xion(1))                                   !     nHI
  nN(2)  = nHe*xHeI                                              !    nHeI
  nN(3)  = nHe*xion(2)                                           !   nHeII
  nI(1)  = nH *xion(1)                                           !    nHII
  nI(2)  = nN(3)                                                 !   nHeII
  nI(3)  = nHe*xion(3)                                           !  nHeIII
  mu= 1./(X*(1.+xion(1)) + 0.25*Y*(1.+xion(2)+2.*xion(3)))   
  if(is_kIR_T) mu=2.33
  TK = dT2 * mu                                        !  Temperature
  if(rt_isTconst) TK=rt_Tconst                         !  Force constant T
  ne= nH*xion(1)+nHE*(xion(2)+2.*xion(3))              !  Electron density
  neInit=ne
  fracMax=0d0    ! Max. fractional update, to check if dt can be increased

  rho = nH / X * mH
  kAbs_loc = kappaAbs
  kSc_loc  = kappaSc
  if(is_kIR_T) then ! k_IR depends on T
     ! Special stuff for Krumholz/Davis experiment
     if(rt_T_rad) then  ! Use radiation temperature for kappa
        E_rad = group_egy(iIR) * ev_to_erg * dNp(iIR)
        TR = max(T2_min_fix,(E_rad*rt_c_cgs/c_cgs/a_r)**0.25)
        dT2 = TR/mu ;   TK = TR
     endif
     kAbs_loc(iIR) = kappaAbs(iIR) * (TK/10d0)**2
     kSc_loc(iIR)  = kappaSc(iIR)  * (TK/10d0)**2
  endif
  ! Set dust absorption and scattering rates [s-1]:
  dustAbs(:)  = kAbs_loc(:) *rho*Zsolar*rt_c_cgs
  dustSc(iIR) = kSc_loc(iIR)*rho*Zsolar*rt_c_cgs

  !(i) UPDATE PHOTON DENSITY AND FLUX ************************************
  if(rt .and. rt_advect) then 
     recRad(1:nGroups)=0. ; phAbs(1:nGroups)=0.              
     ! Scattering rate; reduce the photon flux, but not photon density:
     phSc(1:nGroups)=0.

     ! EMISSION FROM GAS
     if(.not. rt_OTSA .and. rt_advect) then ! ------------- Rec. radiation
        alpha(1) = comp_AlphaA_HII(TK) - comp_AlphaB_HII(TK) 
        ! alpha(2) A-B becomes negative around 1K, hence the max
        alpha(2) = MAX(0.d0,comp_AlphaA_HeII(TK)-comp_AlphaB_HeII(TK))
        alpha(3) = comp_AlphaA_HeIII(TK) - comp_AlphaB_HeIII(TK)
        do i=1,nIons
           if(spec2group(i) .gt. 0)   &     ! Contribution of ion -> group
                recRad(spec2group(i)) = &
                recRad(spec2group(i)) + alpha(i) * nI(i) * ne
        enddo
     endif

     ! ABSORPTION/SCATTERING OF PHOTONS BY GAS
     do i=1,nGroups      ! --------------------------Ionization absorbtion
        phAbs(i) = SUM(nN(:)*signc(i,:)) ! s-1
     end do
     ! IR, optical and UV depletion by dust absorption: ------------------
     if(rt_isIR) & !IR scattering/abs on dust (abs after T update)        
        phSc(iIR)  = phSc(iIR) + dustSc(iIR)                        
     do i=iIR+1,nGroups ! Deplete these photons, since they go into IR 
        phAbs(i) = phAbs(i) + dustAbs(i)                  
     end do

     do i=1,nGroups         ! ------------------- Do the update of N and F
        dNp(i)= MAX(smallNp,                                              &
                   (dt*(recRad(i)+dNpdt(i))+dNp(i)) / (1.d0+dt*phAbs(i)))
        do j=1,nDim
           dFp(i,j) = (dt*dFpdt(i,j)+dFp(i,j))/(1.d0+dt*(phAbs(i)+phSc(i)))
        end do
        call reduce_flux(dFp(i,:),dNp(i)*rt_c_cgs)
        ! ----------------------------------------------------------------
        ! Momentum transfer from photons to gas:                          
        dmom(:)=0d0                                                       
        if(rt_isoPress .and. .not. (rt_isIR .and. i==iIR)) then 
           ! rt_isoPress: assume f=1, where f is reduced flux.
           fluxMag=sqrt(sum(dFp(i,:)))
           if(fluxMag .gt. 0d0)                                          &
                dmom(:) = dmom(:) + dNp(i) * dFp(i,:)/fluxMag * dt       &
                    * (phAbs(i)+phSc(i))  * group_egy(i) * ev_to_erg/c_cgs    
        else ! Use the actual photon flux vector:
           dmom(:) = dmom(:) + dFp(i,:)/rt_c_cgs * dt                    &
                * (phAbs(i)+phSc(i))  * group_egy(i) * ev_to_erg/c_cgs    
        endif                                                             
        ! ----------------------------------------------------------------
     end do
     dp_gas = p_gas + dmom * rt_pressBoost           ! update gas momentum
     dUU=MAXVAL(ABS((dNp-Np)/(Np+Np_MIN))/Np_FRAC)
     if(dUU .gt. 1.) then                                 
        code=1 ;   RETURN                                     ! dt too big
     endif
     fracMax=MAX(fracMax,dUU) ! To check if timestep size can be increased
     dUU=MAXVAL( abs(dFp-Fp) / (abs(Fp)+Fp_MIN)) / Fp_FRAC
     if(dUU .gt. 1.) then                                                 
        code=2 ;   RETURN                                     ! dt too big
     endif
     fracMax=MAX(fracMax,DUU)

     ! Add absorbed UV/optical energy to IR:------------------------------  
     if(rt_isIR) then   
        do i=iIR+1,nGroups
           dNp(iIR) = dNp(iIR) + dustAbs(i) * dt                         &
                               * dNp(i) * group_egy(i) / group_egy(iIR)  
        end do
     endif                                                                
     ! -------------------------------------------------------------------
  endif !if(rt)
  !(ii) UPDATE TEMPERATURE ***********************************************
  if(c_switch .and. cooling .and. .not. rt_T_rad) then
     Hrate=0.                               !  Heating rate [erg cm-3 s-1]
     if(rt .and. rt_advect) then
        do i=1,nGroups                                     !  Photoheating
           Hrate = Hrate + dNp(i) * SUM(nN(:) * PHrate(i,:))
        end do                                                            
     endif                                                                
     if(rt_UV_hom .and. nH .lt. rt_UV_nHSS)                              &
          Hrate = Hrate + SUM(nN(:) * UVrates(:,2))
     Crate = compCoolrate(TK, ne, nN(1), nI(1), nN(2), nN(3), nI(3)      &
                         ,a_exp, dCdT2, RT_OTSA)                !  Cooling
     dCdT2 = dCdT2 * mu                             !  dC/dT2 = mu * dC/dT
     metal_tot=0.d0 ; metal_prime=0.d0                     ! Metal cooling
     if(Zsolar .gt. 0d0) &
          call rt_cmp_metals(T2,nH,mu,metal_tot,metal_prime,a_exp)      
     X_nHkb= X/(1.5 * nH * kB)                     ! Multiplication factor   
     rate  = X_nHkb*(Hrate - Crate - Zsolar*metal_tot)
     dRate = -X_nHkb*(dCdT2 + Zsolar*metal_prime)              ! dRate/dT2
     dUU   = ABS(MAX(T2_min_fix, T2+rate*dt)-T2)     ! 1st order dt constr
     dT2   = MAX(T2_min_fix, T2+rate*dt/(1.-dRate*dt))      ! New T2 value 
     dUU   = MAX(dUU, ABS(dT2-T2)) / (T2+T_MIN) / T_FRAC
     if(dUU .gt. 1.) then                                       ! 10% rule
        code=3 ; RETURN
     endif
     fracMax=MAX(fracMax,dUU)
     TK=dT2*mu
  endif

  if(rt_isIR .and. (kAbs_loc(iIR) .gt. 0d0) .and. .not. rt_T_rad) then
     ! Delta (Cv T) = ( c_red/lambda E - c/lambda a T^4) 
     !              / ( 1/Delta t + 4 c/lambda/C_v a T^3 + c_red/lambda)
     C_v = rho*kb/mh/mu/(gamma-1d0)                                  
     E_rad = group_egy(iIR) * ev_to_erg * dNp(iIR)
     dE_T = (rt_c_cgs * E_rad - c_cgs*a_r*TK**4) &
          / (1d0/kAbs_loc(iIR)/rho/dt + 4d0*c_cgs/C_v*a_r*TK**3+rt_c_cgs)
     dT2 = dT2 + 1d0/mu * 1d0/C_v * dE_T
     dNp(iIR) = dNp(iIR) - dE_T / group_egy(iIR) / ev_to_erg

     dT2 = max(T2_min_fix,dT2)                                   
     dNp(iIR) = max(dNp(iIR), smallNp)                 
     dUU=ABS(dNp(iIR)-Np(iIR))/(Np(iIR)+Np_MIN)/Np_FRAC
     if(dUU .gt. 1.) then                    ! 10% rule for photon density
       code=4 ;   RETURN                          
     endif                                                               
     fracMax=MAX(fracMax,dUU)                                           
                                                                       
     dUU   = ABS(dT2-T2) / (T2+T_MIN) / T_FRAC                         
     if(dUU .gt. 1.) then                                ! 10% rule for T2
        code=5 ; RETURN                                                  
     endif                                                                
     fracMax=MAX(fracMax,dUU)                                             
     TK=dT2*mu                                                          
                                                                        
     call reduce_flux(dFp(iIR,:),dNp(iIR)*rt_c_cgs)           
  endif       

  !(iii) UPDATE xHII******************************************************
  ! First recompute interaction rates since T is updated
  if(rt_OTSA .or. .not. rt_advect) then           !    Recombination rates
     alpha(1) = comp_AlphaB_HII(TK)
     dalpha   = comp_dAlphaB_dT_HII(TK)
  else                               
     alpha(1) = comp_AlphaA_HII(TK)
     dalpha   = comp_dAlphaA_dT_HII(TK)
  endif
  beta(1) = comp_Beta_HI(TK)                      !  Coll. ionization rate
  dBeta   = comp_dBeta_dT_HI(TK)
  cr = beta(1) * ne                               !               Creation
  if(rt) cr = cr + SUM(signc(:,1)*dNp)            !                  [s-1]
  if(rt_UV_hom .and. nH .lt. rt_UV_nHSS) cr = cr + UVrates(1,1) 
  de = alpha(1) * ne                              !            Destruction
  
  ! Not Anninos, but more stable (this IS neccessary, as the one-cell    !
  ! tests oscillate wildly in the Anninos method):                       !
  S  = cr*(1.-dXion(1))-de*dXion(1)
  dUU= ABS(MIN(MAX(dXion(1)+dt*S, x_MIN), 1.)-dXion(1))
  jac=(1.-dXion(1))*(beta(1)*nH-ne*TK*mu*X*dBeta) &          !jac=dS/dxHII
       - cr - de - dXion(1) * (alpha(1)*nH-ne*TK*mu*X*dAlpha)! more Stable
  dXion(1) = xion(1) + dt*(cr*(1.-xion(1))-de*xion(1))/(1.-dt*jac)
  dXion(1) = MIN(MAX(dXion(1), x_MIN),1.d0)
  dUU   = MAX(dUU, ABS(dXion(1)-xion(1))) / (xion(1)+x_MIN) / x_FRAC
  if(dUU .gt. 1.) then
     code=6 ; RETURN
  endif
  fracMax=MAX(fracMax,dUU)
  !End a more stable and accurate integration-----------------------------
  if(isHe) then
     ne= nH*dXion(1)+nHE*(dXion(2)+2.*dXion(3))  ! Upd. bc of changed xhii 
     mu= 1./(X*(1.+dXion(1)) + 0.25*Y*(1.+dXion(2)+2.*dXion(3)))  
     if(.not. rt_isTconst) TK=dT2*mu   !  Update TK because of changed  mu

     !(iv) UPDATE xHeI ***************************************************
     if(rt_OTSA .or. .not. rt_advect) then
        alpha(2) = comp_AlphaB_HeII(TK)
        alpha(3) = comp_AlphaB_HeIII(TK)
     else                               
        alpha(2) = comp_AlphaA_HeII(TK)
        alpha(3) = comp_AlphaA_HeIII(TK)
     endif
     beta(2) = comp_Beta_HeI(TK)
     beta(3) = comp_Beta_HeII(TK)
     ! Creation = recombination of HeII and electrons
     cr = alpha(2) * ne * dXion(2)
     ! Destruction = collisional ionization+photoionization of HeI
     de = beta(2) * ne
     if(rt) de = de + SUM(signc(:,2)*dNp)
     if(rt_UV_hom .and. nH .lt. rt_UV_nHSS) de = de + UVrates(2,1)
     xHeI = (cr*dt+xHeI)/(1.+de*dt)                          !  The update
     xHeI = MIN(MAX(xHeI, 0.),1.)

     !(v) UPDATE xHeII ***************************************************
     ! Creation = coll.- and photo-ionization of HI + rec. of HeIII
     cr = de * xHeI + alpha(3) * ne * dXion(3)
     ! Destruction = rec. of HeII + coll.- and photo-ionization of HeII
     photoRate=0.
     if(rt) photoRate = SUM(signc(:,3)*dNp)
     if(rt_UV_hom .and. nH.lt.rt_UV_nHSS) photoRate=photoRate+UVrates(3,1)
     de = (alpha(2) + beta(3)) * ne + photoRate
     dXion(2) = (cr*dt+dXion(2))/(1.+de*dt)                  !  The update
     dXion(2) = MIN(MAX(dXion(2), x_MIN),1.)

     !(vii) UPDATE xHeIII ************************************************
     ! Creation = coll.- and photo-ionization of HeII
     cr = (beta(3) * ne + photoRate) * dXion(2)            !  xHeII is new
     ! Destruction = rec. of HeIII and e
     de = alpha(3) * ne
     dXion(3) = (cr*dt+dXion(3))/(1.+de*dt)                  !  The update
     dXion(3) = MIN(MAX(dXion(3), x_MIN),1.)

     !(viii) ATOMIC CONSERVATION OF He ***********************************
     if(xHeI .ge. dXion(3)) then   !   Either HeI or HeII is most abundant 
        if(xHeI .le. dXion(2)) dXion(2) = 1.-xHeI-dXion(3) !  HeII most ab
     else                          ! Either HeII or HeIII is most abundant 
        if(dXion(2) .le. dXion(3)) then
           dXion(3) = 1. - xHeI-dXion(2)                          !  HeIII
        else
           dXion(2) = 1. - xHeI-dXion(3)                          !   HeII
        endif
     endif
  endif

  ne = nH*dXion(1)+nHe*(dXion(2)+2.*dXion(3))
  dUU=ABS((ne-neInit)) / (neInit+x_MIN) / x_FRAC
  if(dUU .gt. 1.) then
     code=7 ; RETURN
  endif
  fracMax=MAX(fracMax,dUU)

  if(rt_isTconst) dT2=rt_Tconst/mu

  dT2 = dT2-T2 ; dXion = dXion-xion ; dNp = dNp-Np ; dFp = dFp-Fp
  dp_gas = dp_gas-p_gas ! Now the dUs are really changes, not new values
  !(ix) Check if we are safe to use a bigger timestep in next iteration:
  if(fracMax .lt. 0.5) then
     dt_rec=dt*2.
  else
     dt_rec=dt
  endif
  dt_ok=.true.
  code=0

END SUBROUTINE cool_step

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE display_coolinfo(stopRun, loopcnt, i, dtDone, dt, ddt, nH    &
                            ,T2,  xion,  Np,  Fp,  p_gas                &
                            ,dT2, dXion, dNp, dFp, dp_gas, code)
! Print cooling information to standard output, and maybe stop execution.
!------------------------------------------------------------------------
  use amr_commons
  use rt_parameters
  real(dp),dimension(nIons):: xion, dXion
  real(dp),dimension(nGroups):: Np, dNp
  real(dp),dimension(nGroups, nDim):: Fp, dFp
  real(dp),dimension(nDim):: p_gas, dp_gas
  real(dp)::T2, dT2, dtDone, dt, ddt, nH
  logical::stopRun
  integer::loopcnt,i, code
!------------------------------------------------------------------------
  if(stopRun) write(*, 111) loopcnt
  if(.true.) then
     write(*,900) loopcnt, myid, code, i, dtDone, dt, ddt, rt_c_cgs, nH
     write(*,901) T2,      xion,      Np,      Fp,      p_gas
     write(*,902) dT2,     dXion,     dNp,     dFp,     dp_gas
     write(*,903) dT2/ddt, dXion/ddt, dNp/ddt, dFp/ddt, dp_gas/ddt
     write(*,904) abs(dT2)/(T2+T_MIN), abs(dxion)/(xion+x_MIN),           &
                  abs(dNp)/(Np+Np_MIN), abs(dFp)/(Fp+Fp_MIN)
  endif
  print*,loopcodes
  print*,group_egy(:)
  if(stopRun) then
     print *,'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
     STOP
  endif

111 format(' Stopping because of large number of timestesps in', &
           ' rt_solve_cooling (', I6, ')')
900 format (I3, '  myid=', I2, ' code=', I2, ' i=', I5, ' t=', 1pe12.3,  xs&
            '/', 1pe12.3, ' ddt=', 1pe12.3, ' c=', 1pe12.3, &
            ' nH=', 1pe12.3)
901 format ('  U      =', 20(1pe12.3))
902 format ('  dU     =', 20(1pe12.3))
903 format ('  dU/dt  =', 20(1pe12.3))
904 format ('  dU/U % =', 20(1pe12.3))
END SUBROUTINE display_coolinfo

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE cmp_chem_eq(TK, nH, t_rad_spec, nSpec, nTot, mu)

! Compute chemical equilibrium abundances of e, HI, HII, HeI, HeII, HeIII
! r_rad_spec => photoionization rates [s-1] for HI, HeI, HeII
!------------------------------------------------------------------------
  implicit none
  real(dp),intent(in)::TK, nH
  real(dp),intent(out)::nTot, mu
  real(dp),dimension(1:3),intent(in)::t_rad_spec
  real(dp),dimension(1:6),intent(out)::nSpec!------------------------
  real(dp)::xx, yy
  real(dp)::n_HI, n_HII, n_HEI, n_HEII, n_HEIII, n_E
  real(dp)::t_rad_HI,  t_rad_HEI,  t_rad_HEII
  real(dp)::t_rec_HI,  t_rec_HEI,  t_rec_HEII
  real(dp)::t_ion_HI,  t_ion_HEI,  t_ion_HEII
  real(dp)::t_ion2_HI, t_ion2_HEI, t_ion2_HEII
  real(dp)::x1, err_nE
  integer,parameter::HI=1, HeI=2, HeII=3
!------------------------------------------------------------------------
  xx=(1.-Y)
  yy=Y/(1.-Y)/4.
  
  t_rad_HI   = t_rad_spec(HI)                !      Photoionization [s-1]
  t_rad_HEI  = t_rad_spec(HeI)
  t_rad_HEII = t_rad_spec(HeII)

  if(rt_OTSA) then                           !    Recombination [cm3 s-1]
     t_rec_HI   = comp_AlphaB_HII(TK)        
     t_rec_HEI  = comp_AlphaB_HeII(TK)
     t_rec_HEII = comp_AlphaB_HeIII(TK)
  else 
     t_rec_HI   = comp_AlphaA_HII(TK)        
     t_rec_HEI  = comp_AlphaA_HeII(TK)
     t_rec_HEII = comp_AlphaA_HeIII(TK)
  endif

  t_ion_HI   = comp_Beta_HI(TK)               ! Coll. ionization [cm3 s-1]
  t_ion_HEI  = comp_Beta_HeI(TK)
  t_ion_HEII = comp_Beta_HeII(TK)
  
  n_E = nH        
  err_nE = 1.
  
  do while(err_nE > 1.d-8)
     t_ion2_HI   = t_ion_HI   + t_rad_HI  /MAX(n_E,1e-15*nH)  ! [cm3 s-1]
     t_ion2_HEI  = t_ion_HEI  + t_rad_HEI /MAX(n_E,1e-15*nH)
     t_ion2_HEII = t_ion_HEII + t_rad_HEII/MAX(n_E,1e-15*nH)
     
     n_HI  = t_rec_HI/(t_ion2_HI+t_rec_HI)*nH
     n_HII = t_ion2_HI/(t_ion2_HI+t_rec_HI)*nH
     
     x1=(                                                                &
          t_rec_HEII*t_rec_HEI                                           &
          +t_ion2_HEI*t_rec_HEII                                         &
          +t_ion2_HEII*t_ion2_HEI)                               ! cm6 s-2
     
     n_HEIII = yy*t_ion2_HEII*t_ion2_HEI/x1*nH
     n_HEII  = yy*t_ion2_HEI *t_rec_HEII/x1*nH
     n_HEI   = yy*t_rec_HEII *t_rec_HEI /x1*nH
     
     err_nE = ABS((n_E - (n_HII + n_HEII + 2.*n_HEIII))/nH)
     n_E = 0.5*n_E+0.5*(n_HII + n_HEII + 2.*n_HEIII)
     
  end do
    
  nTOT     = n_E+n_HI+n_HII+n_HEI+n_HEII+n_HEIII
  mu       = nH/xx/nTOT
  nSpec(1) = n_E
  nSpec(2) = n_HI
  nSpec(3) = n_HII
  nSpec(4) = n_HEI
  nSpec(5) = n_HEII
  nSpec(6) = n_HEIII
  
END SUBROUTINE cmp_chem_eq

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE rt_evol_single_cell(astart,aend,dasura,h,omegab,omega0,omegaL   &
                           ,J0min_in,T2end,mu,ne,if_write_result)
!------------------------------------------------------------------------
! Used for initialization of thermal state in cosmological simulations.
!
! astart : valeur du facteur d'expansion au debut du calcul
! aend   : valeur du facteur d'expansion a la fin du calcul
! dasura : la valeur de da/a entre 2 pas de temps
! h      : la valeur de H0/100 
! omegab : la valeur de Omega baryons
! omega0 : la valeur de Omega matiere (total)
! omegaL : la valeur de Omega Lambda
! J0min_in : la valeur du J0min a injecter :
!          Si high_z_realistic_ne alors c'est J0min a a=astart qui
!          est considere
!          Sinon, c'est le J0min habituel.
!          Si J0min_in <=0, les parametres par defaut ou predefinis
!          auparavant sont pris pour le J0min.
! T2end  : Le T/mu en output
! mu     : le poids moleculaire en output
! ne     : le ne en output
! if_write_result : .true. pour ecrire l'evolution de la temperature
!          et de n_e sur l'ecran.
!-------------------------------------------------------------------------
  use amr_commons,only:myid
  use UV_module
  implicit none
  real(kind=8)::astart,aend,T2end,h,omegab,omega0,omegaL,J0min_in,ne,dasura
  logical :: if_write_result
  real(dp)::aexp,daexp,dt_cool,coeff,T2_com, nH_com  
  real(dp),dimension(nIons)::pHI_rates=0., h_rad_spec=0.
  real(kind=8) ::mu
  real(dp) ::cool_tot,heat_tot, mu_dp,diff
  integer::niter
  real(dp) :: n_spec(1:6)
  real(dp),dimension(1:nvector):: T2
  real(dp),dimension(1:nvector, nIons):: xion
  real(dp),dimension(1:nvector, nGroups):: Np, dNpdt
  real(dp),dimension(1:nvector, nGroups, ndim):: Fp, dFpdt
  real(dp),dimension(1:nvector, ndim):: p_gas
  real(dp),dimension(1:nvector)::nH=0., Zsolar=0.
  logical,dimension(1:nvector)::c_switch=.true.
!-------------------------------------------------------------------------
  aexp = astart
  T2_com = 2.726d0 / aexp * aexp**2 / mu_mol
  nH_com = omegab*rhoc*h**2*X/mH

  mu_dp=mu
  call cmp_Equilibrium_Abundances(                                       &
                 T2_com/aexp**2, nH_com/aexp**3, pHI_rates, mu_dp, n_Spec)
  ! Initialize cell state
  T2(1)=T2_com                                          !      Temperature
  xion(1,1)=n_Spec(3)/(nH_com/aexp**3)                  !   HII   fraction
  xion(1,2)=n_Spec(5)/(nH_com/aexp**3)                  !   HeII  fraction
  xion(1,3)=n_Spec(6)/(nH_com/aexp**3)                  !   HeIII fraction
  p_gas(1,:)=0.
  Np(1,:)=0. ; Fp(1,:,:)=0.                  ! Photon densities and fluxes
  dNpdt(1,:)=0. ; dFpdt(1,:,:)=0.                              
  do while (aexp < aend)
     if(rt_UV_hom) call inp_UV_rates_table(1./aexp - 1., UVrates)

     daexp = dasura*aexp
     dt_cool = daexp                                                     &
             / (aexp*100.*h*3.2408608e-20)                               &
             / HsurH0(1.0/dble(aexp)-1.,omega0,omegaL,1.-omega0-omegaL)
     
     nH(1) = nH_com/aexp**3
     T2(1) = T2(1)/aexp**2
     call rt_solve_cooling(T2,xion,Np,Fp,p_gas,dNpdt,dFpdt,nH,c_switch   &
                           ,Zsolar,dt_cool,aexp,1)
     T2(1)=T2(1)*aexp**2
     aexp = aexp + daexp
     if (if_write_result) write(*,'(4(1pe10.3))')                        &
                              aexp,nH(1),T2_com*mu/aexp**2,n_spec(1)/nH(1)
  end do
  T2end=T2(1)/(aexp-daexp)**2
  ne=(n_spec(3)+(n_spec(5)+2.*n_spec(6))*0.25*Y/X)
end subroutine rt_evol_single_cell

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
FUNCTION HsurH0(z,omega0,omegaL,OmegaR)
!------------------------------------------------------------------------
  implicit none
  real(kind=8) :: HsurH0,z,omega0,omegaL,omegaR
!------------------------------------------------------------------------
  HsurH0=sqrt(Omega0*(1.d0+z)**3+OmegaR*(1.d0+z)**2+OmegaL)
END FUNCTION HsurH0

!=======================================================================
subroutine rt_cmp_metals(T2,nH,mu,metal_tot,metal_prime,aexp)
! Taken from the equilibrium cooling_module of RAMSES
! Compute cooling enhancement due to metals
! T2           => Temperature in Kelvin, divided by mu
! nH           => Hydrogen number density (H/cc)
! mu           => Average mass per particle in terms of mH
! metal_tot   <=  Metal cooling contribution to de/dt [erg s-1 cm-3]
! metal_prime <=  d(metal_tot)/dT2 [erg s-1 cm-3 K-1]
!=======================================================================
  implicit none
  real(dp) ::T2,nH,mu,metal_tot,metal_prime,aexp
  ! Cloudy at solar metalicity
  real(dp),dimension(1:91),parameter :: temperature_cc07 = (/ &
       & 3.9684,4.0187,4.0690,4.1194,4.1697,4.2200,4.2703, &
       & 4.3206,4.3709,4.4212,4.4716,4.5219,4.5722,4.6225, &
       & 4.6728,4.7231,4.7734,4.8238,4.8741,4.9244,4.9747, &
       & 5.0250,5.0753,5.1256,5.1760,5.2263,5.2766,5.3269, &
       & 5.3772,5.4275,5.4778,5.5282,5.5785,5.6288,5.6791, &
       & 5.7294,5.7797,5.8300,5.8804,5.9307,5.9810,6.0313, &
       & 6.0816,6.1319,6.1822,6.2326,6.2829,6.3332,6.3835, &
       & 6.4338,6.4841,6.5345,6.5848,6.6351,6.6854,6.7357, &
       & 6.7860,6.8363,6.8867,6.9370,6.9873,7.0376,7.0879, &
       & 7.1382,7.1885,7.2388,7.2892,7.3395,7.3898,7.4401, &
       & 7.4904,7.5407,7.5911,7.6414,7.6917,7.7420,7.7923, &
       & 7.8426,7.8929,7.9433,7.9936,8.0439,8.0942,8.1445, &
       & 8.1948,8.2451,8.2955,8.3458,8.3961,8.4464,8.4967 /)
  real(dp),dimension(1:91),parameter :: excess_cooling_cc07 = (/         &
       & -24.9949,-24.7270,-24.0473,-23.0713,-22.2907,-21.8917,-21.8058, &
       & -21.8501,-21.9142,-21.9553,-21.9644,-21.9491,-21.9134,-21.8559, &
       & -21.7797,-21.6863,-21.5791,-21.4648,-21.3640,-21.2995,-21.2691, &
       & -21.2658,-21.2838,-21.2985,-21.2941,-21.2845,-21.2809,-21.2748, &
       & -21.2727,-21.3198,-21.4505,-21.5921,-21.6724,-21.6963,-21.6925, &
       & -21.6892,-21.7142,-21.7595,-21.7779,-21.7674,-21.7541,-21.7532, &
       & -21.7679,-21.7866,-21.8052,-21.8291,-21.8716,-21.9316,-22.0055, &
       & -22.0800,-22.1600,-22.2375,-22.3126,-22.3701,-22.4125,-22.4353, &
       & -22.4462,-22.4450,-22.4406,-22.4337,-22.4310,-22.4300,-22.4356, &
       & -22.4455,-22.4631,-22.4856,-22.5147,-22.5444,-22.5718,-22.5904, &
       & -22.6004,-22.5979,-22.5885,-22.5728,-22.5554,-22.5350,-22.5159, &
       & -22.4955,-22.4781,-22.4600,-22.4452,-22.4262,-22.4089,-22.3900, &
       & -22.3722,-22.3529,-22.3339,-22.3137,-22.2936,-22.2729,-22.2521 /)
  real(dp),dimension(1:91),parameter :: excess_prime_cc07 = (/           & 
       &   2.0037,  4.7267, 12.2283, 13.5820,  9.8755,  4.8379,  1.8046, &
       &   1.4574,  1.8086,  2.0685,  2.2012,  2.2250,  2.2060,  2.1605, &
       &   2.1121,  2.0335,  1.9254,  1.7861,  1.5357,  1.1784,  0.7628, &
       &   0.1500, -0.1401,  0.1272,  0.3884,  0.2761,  0.1707,  0.2279, &
       &  -0.2417, -1.7802, -3.0381, -2.3511, -0.9864, -0.0989,  0.1854, &
       &  -0.1282, -0.8028, -0.7363, -0.0093,  0.3132,  0.1894, -0.1526, &
       &  -0.3663, -0.3873, -0.3993, -0.6790, -1.0615, -1.4633, -1.5687, &
       &  -1.7183, -1.7313, -1.8324, -1.5909, -1.3199, -0.8634, -0.5542, &
       &  -0.1961, -0.0552,  0.0646, -0.0109, -0.0662, -0.2539, -0.3869, &
       &  -0.6379, -0.8404, -1.1662, -1.3930, -1.6136, -1.5706, -1.4266, &
       &  -1.0460, -0.7244, -0.3006, -0.1300,  0.1491,  0.0972,  0.2463, &
       &   0.0252,  0.1079, -0.1893, -0.1033, -0.3547, -0.2393, -0.4280, &
       &  -0.2735, -0.3670, -0.2033, -0.2261, -0.0821, -0.0754,  0.0634 /)
  real(dp),dimension(1:50),parameter::z_courty=(/                         &
       & 0.00000,0.04912,0.10060,0.15470,0.21140,0.27090,0.33330,0.39880, &
       & 0.46750,0.53960,0.61520,0.69450,0.77780,0.86510,0.95670,1.05300, &
       & 1.15400,1.25900,1.37000,1.48700,1.60900,1.73700,1.87100,2.01300, &
       & 2.16000,2.31600,2.47900,2.64900,2.82900,3.01700,3.21400,3.42100, &
       & 3.63800,3.86600,4.10500,4.35600,4.61900,4.89500,5.18400,5.48800, &
       & 5.80700,6.14100,6.49200,6.85900,7.24600,7.65000,8.07500,8.52100, &
       & 8.98900,9.50000 /)
  real(dp),dimension(1:50),parameter::phi_courty=(/                             &
       & 0.0499886,0.0582622,0.0678333,0.0788739,0.0915889,0.1061913,0.1229119, &
       & 0.1419961,0.1637082,0.1883230,0.2161014,0.2473183,0.2822266,0.3210551, &
       & 0.3639784,0.4111301,0.4623273,0.5172858,0.5752659,0.6351540,0.6950232, &
       & 0.7529284,0.8063160,0.8520859,0.8920522,0.9305764,0.9682031,1.0058810, &
       & 1.0444020,1.0848160,1.1282190,1.1745120,1.2226670,1.2723200,1.3231350, &
       & 1.3743020,1.4247480,1.4730590,1.5174060,1.5552610,1.5833640,1.5976390, &
       & 1.5925270,1.5613110,1.4949610,1.3813710,1.2041510,0.9403100,0.5555344, & 
       & 0.0000000 /)
  real(dp)::TT,lTT,deltaT,lcool,lcool1,lcool2,lcool1_prime,lcool2_prime
  real(dp)::ZZ,deltaZ
  real(dp)::c1=0.4,c2=10.0,TT0=1d5,TTC=1d6,alpha1=0.15
  real(dp)::ux,g_courty,f_courty=1d0,g_courty_prime,f_courty_prime
  integer::iT,iZ
!-------------------------------------------------------------------------
  ZZ=1d0/aexp-1d0
  TT=T2*mu
  lTT=log10(TT)
  ! This is a simple model to take into account the ionization background
  ! on metal cooling (calibrated using CLOUDY). 
  iZ=1+int(ZZ/z_courty(50)*49.)
  iZ=min(iZ,49)
  iZ=max(iZ,1)
  deltaZ=z_courty(iZ+1)-z_courty(iZ)
  ZZ=min(ZZ,z_courty(50))
  ux=1d-4*(phi_courty(iZ+1)*(ZZ-z_courty(iZ))/deltaZ & 
       & + phi_courty(iZ)*(z_courty(iZ+1)-ZZ)/deltaZ )/nH
  g_courty=c1*(TT/TT0)**alpha1+c2*exp(-TTC/TT)
  g_courty_prime=(c1*alpha1*(TT/TT0)**alpha1+c2*exp(-TTC/TT)*TTC/TT)/TT
  f_courty=1d0/(1d0+ux/g_courty)
  f_courty_prime=ux/g_courty/(1d0+ux/g_courty)**2*g_courty_prime/g_courty

  if(lTT.ge.temperature_cc07(91))then
     metal_tot=0d0 !1d-100
     metal_prime=0d0
  else if(lTT.ge.1.0)then
     lcool1=-100d0
     lcool1_prime=0d0
      if(lTT.ge.temperature_cc07(1))then
        iT=1+int((lTT-temperature_cc07(1)) /                             &
             (temperature_cc07(91)-temperature_cc07(1))*90.0)
        iT=min(iT,90)
        iT=max(iT,1)
        deltaT = temperature_cc07(iT+1) - temperature_cc07(iT)
        lcool1 = &
             excess_cooling_cc07(iT+1)*(lTT-temperature_cc07(iT))/deltaT &
           + excess_cooling_cc07(iT)*(temperature_cc07(iT+1)-lTT)/deltaT 
        lcool1_prime  =                                                  &
             excess_prime_cc07(iT+1)*(lTT-temperature_cc07(iT))/deltaT   &
           + excess_prime_cc07(iT)*(temperature_cc07(iT+1)-lTT)/deltaT 
     endif
     ! Fine structure cooling from infrared lines
     lcool2=-31.522879+2.0*lTT-20.0/TT-TT*4.342944d-5
     lcool2_prime=2d0+(20d0/TT-TT*4.342944d-5)*log(10d0)
     ! Total metal cooling and temperature derivative
     metal_tot=10d0**lcool1+10d0**lcool2
     metal_prime=(10d0**lcool1*lcool1_prime+10d0**lcool2*lcool2_prime)/metal_tot
     metal_prime=metal_prime*f_courty+metal_tot*f_courty_prime
     metal_tot=metal_tot*f_courty
  else
     metal_tot=0d0 !1d-100
     metal_prime=0d0
  endif

  metal_tot=metal_tot*nH**2
  metal_prime=           &   ! Convert from DlogLambda/DlogT to DLambda/DT
       metal_prime * metal_tot/TT * mu

end subroutine rt_cmp_metals

END MODULE rt_cooling_module

!************************************************************************
SUBROUTINE updateRTGroups_CoolConstants()
! Update photon group cooling and heating constants, to reflect an update
! in rt_c_cgs and in the cross-sections and energies in the groups.
!------------------------------------------------------------------------
  use rt_cooling_module
  use rt_parameters
  implicit none
  integer::iP, iI
!------------------------------------------------------------------------
  signc=group_csn*rt_c_cgs                                    ! [cm3 s-1]
  sigec=group_cse*rt_c_cgs                                    ! [cm3 s-1]
  do iP=1,nGroups
     do iI=1,nIons               ! Photoheating rates for photons on ions
        PHrate(iP,iI) =  ev_to_erg * &        ! See eq (19) in Aubert(08)
             (sigec(iP,iI) * group_egy(iP) - signc(iP,iI)*ionEvs(iI))
        PHrate(iP,iI) = max(PHrate(iP,iI),0d0) !      No negative heating
     end do
  end do
END SUBROUTINE updateRTGroups_CoolConstants

!************************************************************************
SUBROUTINE reduce_flux(Fp, cNp)
! Make sure the reduced photon flux is less than one
!------------------------------------------------------------------------
  use rt_parameters
  implicit none
  real(dp),dimension(ndim):: Fp
  real(dp):: cNp, fred
!------------------------------------------------------------------------
  fred = sqrt(sum(Fp**2))/cNp
  if(fred .gt. 1.d0) Fp = Fp/fred
END SUBROUTINE reduce_flux
