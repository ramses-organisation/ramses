! Non-equlibrium (in H2, H, He, e, C, N, O, Ne, Mg, Si, S, & Fe)
! cooling module for radiation-hydrodynamics.
! For details, see Rosdahl et al. (2013), Rosdahl & Teyssier (2015),
! Katz et al. (2017), Nickerson, Teyssier, & Rosdahl (2018), 
! Katz (2022), and Katz et al. (2022)
! Joki Rosdahl, Sarah Nickerson, Andreas Bleuler, Romain Teyssier
! and Harley Katz
! NOTE: T2=T/mu, Np = photon density, Fp = photon flux,
module rtz_cooling_module
  use amr_parameters, only: ndim, dp, nvector
  use rt_parameters
  use constants
  use rtz_module
  implicit none

  private   ! default
  public rtz_solve_cooling, rtz_set_model, PHrate, T2_min_fix

  real(dp),parameter::T2_min_fix=1d-2 ! Min temperature [K]
  real(dp),parameter::T_min=0.1, T_frac=0.1
  real(dp),parameter::x_min=1d-20, x_fm=1d-6, x_frac=0.1
  real(dp),parameter::Np_min=1d-13, Np_frac=0.2
  real(dp),parameter::Fp_frac=0.5

CONTAINS

SUBROUTINE rtz_set_model(h, omegab, omega0, omegaL, astart_sim, T2_sim)
  ! Initialize cooling. All these parameters are unused at the moment and
  ! are only there for the original cooling-module.
  ! h (dble)            => H0/100
  ! omegab (dble)       => Omega Baryons
  ! omega0 (dble)       => Omega Matter total
  ! omegaL (dble)       => Omega Lambda
  ! astart_sim (dble)   => Redshift at which we start the simulation
  ! T2_sim (dble)      <=  Starting temperature in simulation?
  !-------------------------------------------------------------------------
  integer::ilevel
  real(dp) :: h, omegab, omega0, omegaL, astart_sim, T2_sim
  real(dp) :: z_decoupling, z_start_sim
  real(dp) :: astart=0.0001, aend, dasura, mu=1., ne
  !-------------------------------------------------------------------------

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

  if(nrestart==0 .and. cosmo) then
      ! Use approximate temperature evolution from 
      ! https://arxiv.org/pdf/astro-ph/0608032
      ! This ignores compton cooling but should be ok
      ! Just don't start the simulation at too high of redshift
      z_decoupling = 150.d0 * ((omegab * h * h /0.023d0)**(2.d0/5.d0))
      z_decoupling = z_decoupling - 1.d0

      z_start_sim = (1.d0 / astart_sim) - 1.d0

      ! If redshift is less than the decoupling redshift, scale by (1+z)^2      
      if (z_start_sim .lt. z_decoupling) then
         T2_sim = 2.725 * (1.d0 + z_decoupling)
         T2_sim = T2_sim * ((1.d0 + z_start_sim)/(1.d0 + z_decoupling))**2.d0
      ! If redshift is greater than the decoupling redshift, use CMB temp.
      else
         T2_sim = 2.725 * (1.d0 + z_start_sim)
      end if
  end if                                   

END SUBROUTINE rtz_set_model

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE rtz_solve_cooling(T2, aexp, xion, nElement, nCO, &
#ifdef RT
     & Np, Fp, p_gas, dNpdt, dFpdt, ilevel, &
#endif
     & dt, nCell, dx_SS_H2)
  ! Semi-implicitly solve for new temperature, ionization states,
  ! photon density/flux, and gas velocity in a number of cells.
  ! Parameters:
  ! T2     <=> T/mu [K]
  ! xion   <=> 27x27 array with all of the mass fractions of each ionization state
  ! Np     <=> nGroups photon number densities [cm-3]
  ! Fp     <=> nGroups * ndim photon number fluxes [cm-2 s-1]
  ! p_gas  <=> ndim gas momentum densities [cm s-1 g cm-3]
  ! dNpdt   =>  Op split increment in photon densities during dt
  ! dFpdt   =>  Op split increment in photon flux magnitudes during dt
  ! nElement=>  Number density of each element [cm^-3] (length 27)
  ! c_switch=>  Cooling switch (1 for cool/heat, 0 for no cool/heat) (OFF)
  ! dt      =>  Timestep size             [s]
  ! nCell   =>  Number of cells (length of all the above vectors)
  ! ilevel  =>  Refinement levels
  ! dx_SS_H2=>  Cell size [cm] used for H2 self shielding 
  !
  ! We use a slightly modified method of Anninos et al. (1997).
  !-------------------------------------------------------------------------
  implicit none
  real(dp):: aexp
  real(dp),dimension(1:nvector):: T2
  real(dp),dimension(1:n_elements, 1:n_elements, 1:nvector):: xion
  real(dp),dimension(1:n_elements, 1:nvector):: nElement 
  real(dp),dimension(1:nvector):: nH
  real(dp),dimension(1:nvector):: nCO
#ifdef RT
  real(dp),dimension(1:ndim, 1:nvector):: p_gas
  real(dp),dimension(1:nGroups, 1:nvector):: Np, dNpdt
  real(dp),dimension(1:ndim, 1:nGroups, 1:nvector):: Fp, dFpdt
  integer::ilevel
  real(dp):: dx_SS_H2 
#endif
!  logical,dimension(1:nvector):: c_switch
  real(dp)::dt
  integer::ncell
  !--------------------------------------------------------
  real(dp),dimension(1:nvector):: tLeft, ddt
  logical:: dt_ok
  real(dp):: dt_rec
  real(dp):: dT2
  real(dp),dimension(1:n_elements,1:n_elements):: dXion
  real(dp),dimension(1:n_elements):: dnElement 
  real(dp)::dCO
  integer::i, ia, nAct, nAct_next, loopcnt, code
  integer,dimension(1:nvector):: indAct              ! Active cell indexes
  real(dp):: one_over_x_FRAC, one_over_T_FRAC
#ifdef RT
  integer::ig
  real(dp):: one_over_rt_c_cgs, one_over_egy_IR_erg
  real(dp):: one_over_Np_FRAC, one_over_Fp_FRAC
  real(dp),dimension(1:ndim):: dp_gas
  real(dp),dimension(nGroups):: dNp
  real(dp),dimension(1:ndim, 1:nGroups):: dFp
  real(dp),dimension(1:nGroups):: group_egy_ratio, group_egy_erg
#endif
  integer*8,dimension(20)::loopCodes=0
  integer::iElement, ion_fracs
  real(dp)::current_mass_frac
  integer:: i_interp, convergence_counter
  integer :: base_unit = 100
  integer :: element_unit, j, iIon
  character(len=50) :: element_filename
  real(dp), allocatable :: ytmp(:)
  real(dp), dimension(1:50)::saved_cooling_rates
  character(len=20), dimension(1:50)::saved_cooling_rates_names
  real(dp)::TK_to_save(1:nvector), mu_to_save(1:nvector)

#ifdef RT
  call rtz_updateRTGroups_CoolConstants(ilevel)
#endif
  ! Store some temporary variables reduce computations
  one_over_T_FRAC = 1d0 / T_FRAC
  one_over_x_FRAC = 1d0 / x_FRAC
#ifdef RT
  one_over_Np_FRAC = 1d0 / Np_FRAC
  one_over_Fp_FRAC = 1d0 / Fp_FRAC
  one_over_rt_c_cgs = 1d0 / rt_c_cgs(ilevel)
  group_egy_erg(1:nGroups) = group_egy(1:nGroups) * eV2erg
  if(rt_isIR) then
     group_egy_ratio(1:nGroups) = group_egy(1:nGroups) / group_egy(iIR)
     one_over_egy_IR_erg = 1d0 / group_egy_erg(iIR)
  endif
#endif
  !-----------------------------------------------------------------------

  ! Check if we are running in equilibrium mode
  if (rtz_equilibrium_test.gt.0) then

      ! Open files for all elements
      do iElement= 1, n_elements
         if (elements(iElement)%atomic_number .gt. 0) then
            element_unit = base_unit + iElement
            write(element_filename, '("element_", I0, "_ions.dat")') iElement
            open(unit=element_unit, file=element_filename, status='unknown')
         end if
      end do

      ! Open file to save cooling and heating rates
      open(unit=base_unit+100, file='coolrates.dat', status='unknown')

      if (rtz_equilibrium_test.eq.1) then 
         !!! USE FOR EQM TESTS WITH COOLING AT CONSTANT RHO
         T2 = 1.d5 ! --> initialize at high temperature
         ! Set the ionization states to neutral
         xion = 0.d0
         do iElement=1,n_elements
            if (elements(iElement)%atomic_number .gt. 0) then
               xion(iElement,1,:) = 1.d0
            end if
         end do
      end if

      do i_interp = 1,300
         ! Initialize the convergence counter
         convergence_counter = 0

         if (rtz_equilibrium_test.eq.2) then
            !!! USE FOR EQM TESTS AT CONSTANT T
            ! Set the temperature
            rt_isTconst = .true.
            rt_Tconst = 10.d0**(((8.d0 - 2.d0) * (real(i_interp,dp) - 1.d0)/(300.d0-1.d0)) + 2.d0)
         end if

         if (rtz_equilibrium_test.eq.1) then 
            !!! USE FOR EQM TESTS WITH COOLING AT CONSTANT RHO
            ! Interpolate over density
            nElement(1:n_elements,1:ncell)  = 0.d0  ! Initialize to zero
            nElement(1,1:ncell)  = 10.d0**(((5.d0 - (-2.d0)) * (real(i_interp,dp) - 1.d0)/(300.d0-1.d0)) + (-2.d0))   
            nElement(2,1:ncell)  = nElement(1,1:ncell) * 8.51d-02 ! Helium
            nElement(6,1:ncell)  = nElement(1,1:ncell) * 2.69d-04 * z_ave ! Carbon
            nElement(7,1:ncell)  = nElement(1,1:ncell) * 6.76d-05 * z_ave ! Nitrogen
            nElement(8,1:ncell)  = nElement(1,1:ncell) * 4.90d-04 * z_ave ! Oxygen
            nElement(10,1:ncell) = nElement(1,1:ncell) * 8.51d-05 * z_ave ! Neon
            nElement(12,1:ncell) = nElement(1,1:ncell) * 3.98d-05 * z_ave ! Magnesium
            nElement(14,1:ncell) = nElement(1,1:ncell) * 3.24d-05 * z_ave ! Silicon
            nElement(16,1:ncell) = nElement(1,1:ncell) * 1.32d-05 * z_ave ! Sulfur
            nElement(26,1:ncell) = nElement(1,1:ncell) * 3.16d-05 * z_ave ! Iron
            nCO(1:ncell) = 0.0 ! CO
         end if

         tleft(1:ncell) = 1.d40             ! Set to an arbitrarily large number
         ddt(1:ncell) = 10000.d0 * 365.25d0 * 60.d0 * 60.d0 ! First guess at sub-timestep lengths

         do i=1,ncell
            indact(i) = i                   !      Set up indexes of active cells

            ! set nH
            nH(i) = nElement(1,i)

            ! Ensure all state vars are legal:
            T2(i) = MAX(T2(i), T2_min_fix)

            ! Make sure the ionization fractions don't go below the min or max
            xion(1:n_elements,1:n_elements,i) = MIN(MAX(xion(1:n_elements,1:n_elements,i), x_MIN),1d0)

            ! Loop over each element and ensure that the ionization fractions
            ! sum to 1
            do iElement=1,n_elements
               ! Check if we are actually using that element
               if (elements(iElement)%atomic_number .gt. 0) then
                  ! REDUCE SO THAT IONIZATION FRACTIONS SUM TO 1
                  ion_fracs = elements(iElement)%n_ions + elements(iElement)%n_mol
                  current_mass_frac = sum(xion(iElement,1:ion_fracs,i))
                  do iIon=1,ion_fracs
                     xion(iElement,iIon,i) = xion(iElement,iIon,i) + ((1.d0 - current_mass_frac) * (xion(iElement,iIon,i) / current_mass_frac))
                  end do
               end if
            end do
         end do

         ! Loop until all cells have tleft=0
         ! **********************************************
         nAct=nCell                                      ! Currently active cells
         loopcnt=0 !; n_cool_cells=n_cool_cells+nCell     !             Statistics
         do while (nAct .gt. 0)      ! Iterate while there are still active cells
            loopcnt=loopcnt+1 !  ;   tot_cool_loopcnt=tot_cool_loopcnt+nAct
            nAct_next=0                     ! Active cells for the next iteration
            do ia=1,nAct                             ! Loop over the active cells
               i = indAct(ia)                        !                 Cell index
               call rtz_cool_step(i)
               if(.not. dt_ok) then
                  convergence_counter = 0 ! Reset the convergence counter
                  ddt(i)=ddt(i)/2.                    ! Try again with smaller dt
                  nAct_next=nAct_next+1 ; indAct(nAct_next) = i
                  loopCodes(code) = loopCodes(code)+1
                  cycle
               endif
               convergence_counter = convergence_counter + 1
               ! Update the cell state (advance the time by ddt):
               T2(i) = T2(i) + dT2
               ! Check for convergence
               xion(:,:,i) = xion(:,:,i) + dXion(:,:)
               nElement(:,i) = nElement(:,i) + dnElement(:)
#ifdef CO
               nCO(i) = nCO(i) + dCO
#endif
               if (convergence_counter .gt. rtz_eqm_min_its) then
                  tleft(i) = 0.0 ! Finish the cell if we have reached convergence
               else
                  tleft(i)=tleft(i)-ddt(i)
                  ! Take at least X iterations
                  if (convergence_counter .lt. rtz_eqm_min_its) then
                     tleft(i)=max(tleft(i),ddt(i))
                  endif
               endif
               if(tleft(i) .gt. 0.) then           ! Not finished with this cell
                  nAct_next=nAct_next+1 ; indAct(nAct_next) = i
               else if(tleft(i) .lt. 0.) then        ! Overshot by abs(tleft(i))
                  print*,'In rtz_solve_cooling: tleft < 0  !!'
                  stop
               endif
               ddt(i)=min(dt_rec,tleft(i))    ! Use recommended dt from rtz_cool_step
            end do ! end loop over active cells
            nAct=nAct_next
         end do ! end iterative loop

         ! Write for debugging
         if (rtz_equilibrium_test.eq.2) then
            if (isH2_rtz) then 
               write(*,*) rt_Tconst, loopcnt, xion(1,1,1), xion(1,2,1), xion(1,3,1)
            else
               write(*,*) rt_Tconst, loopcnt, xion(1,1,1), xion(1,2,1)
            end if
         end if

         if (rtz_equilibrium_test.eq.1) then
            if (isH2_rtz) then 
               if (isCO_rtz) then
                  write(*,*) nH(i), TK_to_save(i), T2(i), mu_to_save(i), loopcnt, nCO(i)/nH(i), (nCO(i)+nElement(6,i))/nH(i), (nCO(i)+nElement(8,i))/nH(i)
               else
                  write(*,*) nH(i), TK_to_save(i), T2(i), mu_to_save(i), loopcnt, xion(1,1,1), xion(1,2,1), xion(1,3,1)
               end if
            else
               write(*,*) nH(i), TK_to_save(i), T2(i), mu_to_save(i), loopcnt, xion(1,1,1), xion(1,2,1)
            end if

            ! Write the cooling and heating rates to file
            if (i_interp.eq.1) write(base_unit+100,'(*(A20, ", "))') 'rho', 'T', 'Tmu', 'mu', saved_cooling_rates_names
            write(base_unit+100,'(*(ES15.6E3, ", "))') nH(i), TK_to_save(i), T2(i), mu_to_save(i), saved_cooling_rates
         end if

         ! Write data to file
         do iElement = 1, n_elements
            if (elements(iElement)%atomic_number .gt. 0) then
               element_unit = base_unit + iElement
               write(element_unit,'(ES15.6, I14, *(ES15.6))') rt_Tconst, loopcnt, &
                  (xion(iElement,j,1), j=1,elements(i)%n_ions)
            end if
         end do

      end do

      ! Close all element files
      do i = iElement, n_elements
         if (elements(iElement)%atomic_number .gt. 0) then
            close(base_unit + iElement)
         end if
      end do

      close(unit=base_unit+100)

      write(*,*) '!************************************************!'
      stop "Program terminated due to equilibrium test"

  ! Otherwise perform the normal loop
  else
      tleft(1:ncell) = dt                !       Time left in dt for each cell
      ddt(1:ncell) = dt                  ! First guess at sub-timestep lengths

      do i=1,ncell
         indact(i) = i                   !      Set up indexes of active cells

         ! set nH
         nH(i) = nElement(1,i)

         ! Ensure all state vars are legal:
         T2(i) = MAX(T2(i), T2_min_fix)

         ! Make sure the ionization fractions don't go below the min or max
         xion(1:n_elements,1:n_elements,i) = MIN(MAX(xion(1:n_elements,1:n_elements,i), x_MIN),1d0)

         ! Loop over each element and ensure that the ionization fractions
         ! sum to 1
         do iElement=1,n_elements
            ! Check if we are actually using that element
            if (elements(iElement)%atomic_number .gt. 0) then
               ! REDUCE SO THAT IONIZATION FRACTIONS SUM TO 1
               ion_fracs = elements(iElement)%n_ions + elements(iElement)%n_mol
               current_mass_frac = sum(xion(iElement,1:ion_fracs,i))
               do iIon=1,ion_fracs
                  xion(iElement,iIon,i) = xion(iElement,iIon,i) + ((1.d0 - current_mass_frac) * (xion(iElement,iIon,i) / current_mass_frac))
               end do
            end if
         end do
#ifdef RT
         do ig=1,nGroups
            Np(ig,i) = MAX(smallNp, Np(ig,i))
            call reduce_flux(Fp(:,ig,i),Np(ig,i)*rt_c_cgs(ilevel))
         end do
#endif
      end do

      ! Loop until all cells have tleft=0
      ! **********************************************
      nAct=nCell                                      ! Currently active cells
      loopcnt=0 !; n_cool_cells=n_cool_cells+nCell     !             Statistics
      do while (nAct .gt. 0)      ! Iterate while there are still active cells
         loopcnt=loopcnt+1 !  ;   tot_cool_loopcnt=tot_cool_loopcnt+nAct
         nAct_next=0                     ! Active cells for the next iteration
         do ia=1,nAct                             ! Loop over the active cells
            i = indAct(ia)                        !                 Cell index
            call rtz_cool_step(i)
            if(.not. dt_ok) then
               ddt(i)=ddt(i)/2.                    ! Try again with smaller dt
               ! ddt(i) = dt_rec              ! Potentially optimized approach
               nAct_next=nAct_next+1 ; indAct(nAct_next) = i
               loopCodes(code) = loopCodes(code)+1
               cycle
            endif
            ! Update the cell state (advance the time by ddt):
            T2(i) = T2(i) + dT2
            xion(:,:,i) = xion(:,:,i) + dXion(:,:)
            nElement(:,i) = nElement(:,i) + dnElement(:)
#ifdef CO
            nCO(i) = nCO(i) + dCO
#endif
#ifdef RT
            Np(:,i) = Np(:,i) + dNp(:)
            Fp(:,:,i) = Fp(:,:,i) + dFp(:,:)
            p_gas(:,i) = p_gas(:,i) + dp_gas(:)
#endif
            tleft(i)=tleft(i)-ddt(i)
            if(tleft(i) .gt. 0.) then           ! Not finished with this cell
               nAct_next=nAct_next+1 ; indAct(nAct_next) = i
            else if(tleft(i) .lt. 0.) then        ! Overshot by abs(tleft(i))
               print*,'In rtz_solve_cooling: tleft < 0  !!'
               stop
            endif
            ddt(i)=min(dt_rec,tleft(i))    ! Use recommended dt from rtz_cool_step
         end do ! end loop over active cells
         nAct=nAct_next
      end do ! end iterative loop

   ! loop statistics
   max_cool_loopcnt=max(max_cool_loopcnt,loopcnt)

  end if

contains

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  SUBROUTINE rtz_cool_step(icell)
    ! Compute change in cell state in timestep ddt(icell), or set in dt_rec
    ! a recommendation for new timestep if ddt(icell) proves too large.
    ! T2      => T/mu [K]                               -- dT2 is new value
    ! xion    => nion ionization fractions              --     dXion is new
    ! Np      => nGroups photon number densities [cm-3]  -- dNp is new value
    ! Fp      => nGroups * ndim photon fluxes [cm-2 s-1] -- dFp is new value
    ! p_gas   => ndim gas momenta [cm s-1 g cm-3]       --    dp_gas is new
    ! dNpdt   =>  Op split increment in photon densities during dt
    ! dFpdt   =>  Op split increment in photon flux magnitudes during dt
    ! nH      =>  Hydrogen number densities [cm-3]
    ! c_switch=>  Cooling switch (1 for cool/heat, 0 for no cool/heat) (OFF)
    ! dt      =>  Timestep size [s]
    ! dt_ok   <=  .f. if timestep constraints were broken, .t. otherwise
    ! dt_rec  <=  Recommended timesteps for next iteration
    ! code    <= Error code in cool step, if dt_ok=.f.
    !
    ! The original values, T2, xion etc, must stay unchanged, while dT2,
    ! dxion etc contain the new values (the difference at the end of the
    ! routine).
    !-----------------------------------------------------------------------
    use amr_commons
    use const
    use collisional_ionization_module
    use recombination_module
    use charge_exchange_module
    use dust_recombination_module
    use photoionization_UVB_module
    use cosmic_ray_ionization_module
    use molecules_module
    use rtz_coolrates_module, only: all_cooling
    implicit none
    integer, intent(in):: icell
    !-----------------------------------------------------------------------
    real(dp):: dUU, fracMax
    real(dp):: mu, TK, ne, neInit
   !  real(dp):: xHI,dxHI, xH2=0d0,dXH2=0d0, xHeI,dxHeI
    real(dp):: Crate, Crate_prime, Crate_prime_a, Crate_prime_b, dCdT2, X_nHkb, rate, dRate, cr, de=0d0
    real(dp):: ss_factor, f_dust
#ifdef RT
    integer::igroup,idim
    real(dp),dimension(ndim):: dmom
    real(dp),dimension(nGroups):: recRad, phAbs, phSc, dustAbs
    real(dp),dimension(nGroups):: dustSc, kAbs_loc, kSc_loc
    real(dp):: TR, one_over_C_v, E_rad, dE_T
   !  real(dp):: G0, eff_peh
    real(dp):: fluxMag, mom_fact
#endif
    real(dp):: rho
    !-----------------------------------------------------------------------
    ! Variables specific to RTZ
    real(dp):: xe
    real(dp):: dust_effective_number_density, dust_to_gas_mass_ratio_over_mw
    real(dp):: HI_number_density, HII_number_density
    real(dp):: paired_ion_number_density
    real(dp):: total_cosmic_ray_ionization_rate, H2_cosmic_ray_ionization_rate
    real(dp):: phi_s, cosmic_ray_scale_factor, primary_cosmic_ray_ionization_rate
    real(dp):: UV_background_G0
    integer:: atomic_number, n_ions, i_other_Element, i_other_Ion, i_current_Element
    integer:: i_current_Ion
    real(dp):: Zsolar, total_G0
    real(dp):: alpha_H2_loc, beta_H2_loc, cr_H2, de_H2, xH2_loc, f_shd
    real(dp):: nElement_dep(n_elements)
#ifdef CO
    real(dp):: cr_CO, de_CO, delta_CO, max_delta_CO, min_delta_CO
    real(dp):: n_C_og, n_O_og, n_CII, n_H2, n_OI, nCII_new, nCO_new, nOI_new
    real(dp):: tot_C, tot_O, x_OI
#endif
    real(dp),dimension(1:27,10):: saved_rates
    !-----------------------------------------------------------------------

    ! RTZ variable initialization

    ! Include dust if we are tracking oxygen
    if (elements(8)%atomic_number .lt. 1) then
       dust_to_gas_mass_ratio_over_mw = 0.d0
       Zsolar = 1.d-40
       nElement_dep(1:n_elements) = nElement(1:n_elements,icell)
    else
       Zsolar = 12.d0 + log10((nElement(8, icell)+1.d-20)/(nElement(1, icell)+1.d-20))
       dust_to_gas_mass_ratio_over_mw = dust_to_gas_scale_RR14(Zsolar)
       Zsolar = Zsolar - 8.69d0
       do iElement=1,n_elements
          nElement_dep(iElement) = nElement(iElement,icell) * (1.d0 - ((1.d0 - elements(iElement)%depletion) * dust_to_gas_mass_ratio_over_mw))
       end do
    end if

    primary_cosmic_ray_ionization_rate = rtz_primary_cosmic_ray_ionization_rate
    UV_background_G0 = rtz_UV_background_G0

    total_G0 = UV_background_G0

    ! END RTZ variable initialization

    dt_ok=.false.
    ! U contains the original values, dU the updated ones
    dT2 = T2(icell) ; dXion(:,:) = xion(:,:,icell)
    dnElement(:) = nElement(:,icell)
#ifdef CO
    dCO = nCO(icell)
#endif

    f_shd = 1.d0
    if (isH2_rtz) then
       f_shd = comp_SH2(nElement_dep(1)*dXion(1,3), dx_SS_H2) * comp_Sd(nElement_dep(1)*dXion(1,1), nElement_dep(1)*dXion(1,3), dx_SS_H2, dust_to_gas_mass_ratio_over_mw)
    end if

#ifdef RT
    dNp(:) = Np(:,icell) ; dFp(:,:) = Fp(:,:,icell)
    dp_gas(:) = p_gas(:,icell)
#endif

    ne = getNe(dXion, nElement_dep(:))
    neInit = ne
    mu = getMu_RTZ(ne, nElement_dep, dXion)
    TK = dT2 * mu                                        !      Temperature
    if(rt_isTconst) TK=rt_Tconst                         ! Force constant T
    fracMax = 0d0 ! Max fractional update, to check if dt can be increased
    ss_factor = 1d0                  ! UV background self_shielding factor
    if (rtz_equilibrium_test.lt.0) then 
       if(self_shielding) ss_factor = exp(-nH(icell)/1d-2)
    end if
    rho = get_rho_rtz(nElement(:,icell))

    ! RTZ -- initialize cosmic ray variables
    ! Measure phi_s for secondary CR ionization
    xe = ne / nElement_dep(1)
    phi_s = secondary_cr_rates(xe)
    ! Now calculate the total CR ionization rate
    total_cosmic_ray_ionization_rate = primary_cosmic_ray_ionization_rate * (1.d0 + phi_s)
    H2_cosmic_ray_ionization_rate = 2.d0 * primary_cosmic_ray_ionization_rate * (1.d0 + phi_s)
    cosmic_ray_scale_factor = H2_cosmic_ray_ionization_rate / 1.d-16
    ! END RTZ

#ifdef RT
    ! Set dust opacities--------------------------------------------------
    kAbs_loc = kappaAbs
    kSc_loc = kappaSc
    if(is_kIR_T) then                            ! k_IR depends on T_rad
          ! For the radiation temperature,  weigh the energy in each group
            ! by its opacity over IR opacity (derived from IR temperature)
       E_rad = group_egy_erg(iIR) * dNp(iIR)
       TR = max(0d0,(E_rad*rt_c_fraction(ilevel)/a_r)**0.25)      ! IR temperature
       kAbs_loc(iIR) = kappaAbs(iIR) * (TR/10d0)**2
       do iGroup=1,nGroups
          if(iGroup .ne. iIR)                                            &
               E_rad = E_rad + kAbs_loc(iGroup) / kAbs_loc(iIR)          &
               * group_egy_erg(iGroup) * dNp(iGroup)
       end do
       TR = max(0d0,(E_rad*rt_c_fraction(ilevel)/a_r)**0.25)  ! Rad. temperature
       if(rt_T_rad) then      ! Use radiation temperature for everything
          dT2 = TR/mu ;   TK = TR
       endif
                 ! Set the IR opacities according to the rad. temperature:
       kAbs_loc(iIR) = kappaAbs(iIR) * (TR/10d0)**2 * exp(-TR/1d3)
       kSc_loc(iIR)  = kappaSc(iIR)  * (TR/10d0)**2 * exp(-TR/1d3)
    endif ! if(is_kIR_T)
                         ! Set dust absorption and scattering rates [s-1]:
    ! TODO(code): update this to harley's dust model
    f_dust = 0.0
    dustAbs(:) =kAbs_loc(:) *rho*Zsolar*f_dust*rt_c_cgs(ilevel)
    dustSc(iIR)=kSc_loc(iIR)*rho*Zsolar*f_dust*rt_c_cgs(ilevel)

    ! UPDATE PHOTON DENSITY AND FLUX *************************************
    if(rt_advect) then
       recRad(1:nGroups)=0. ; phAbs(1:nGroups)=0.
       ! Scattering rate; reduce the photon flux, but not photon density:
       phSc(1:nGroups)=0.

       ! HKnote: OTSA is required with RTZ (for now)

       ! ABSORPTION/SCATTERING OF PHOTONS BY GAS
       do igroup=1,nGroups       ! ----------------Ionization absorbtion
          do i_current_Element=1,n_elements ! loop over elements
             if (elements(i_current_Element)%atomic_number.gt.0) then
                do i_current_Ion=1,elements(i_current_Element)%n_ions-1 ! loop over ions
                   phAbs(igroup) = phAbs(igroup) + nElement_dep(i_current_Element) * dXion(i_current_Element, i_current_Ion) * signc(igroup,i_current_Element,i_current_Ion)  ! s-1
                end do  ! end loop over ions
             end if
          end do ! end loop over elements

          ! Deal with molecules separately
          if (elements(1)%atomic_number.gt.0 .and. isH2_rtz) then
             if (isLW(igroup).eq.1.d0) then 
                phAbs(igroup) = 0.5d0 * nElement_dep(1) * dXion(1, 3) * signc(igroup,1,3) * f_shd  ! s-1
             else
                phAbs(igroup) = 0.5d0 * nElement_dep(1) * dXion(1, 3) * signc(igroup,1,3)
             end if
          end if
       end do
       ! IR, optical and UV depletion by dust absorption: ----------------
       ! IR scattering/abs on dust (abs after T update)
       if(rt_isIR) phSc(iIR)  = phSc(iIR) + dustSc(iIR)
       do igroup=1,nGroups      ! Deplete photons, since they go into IR
          if( .not. (rt_isIR .and. igroup.eq.iIR) ) & ! IR done elsewhere
               phAbs(igroup) = phAbs(igroup) + dustAbs(igroup)
       end do

       dmom(1:ndim)=0d0
       do igroup=1,nGroups     ! ----------------- Do the update of N and F
          dNp(igroup)= MAX(smallNp,                                      &
                        (ddt(icell)*(recRad(igroup)+dNpdt(igroup,icell)) &
                                    +dNp(igroup))                        &
                        / (1d0+ddt(icell)*phAbs(igroup)))
          dUU = ABS(dNp(igroup)-Np(igroup,icell))                        &
                /(Np(igroup,icell)+Np_MIN) * one_over_Np_FRAC
          fracMax=MAX(fracMax,dUU)      ! To check if ddt can be increased
          if(dUU .gt. 1d0) then
             dt_rec = 0.9d0 * ddt(icell) / sqrt(2.d0+fracMax)
             code=1 ;   RETURN                        ! ddt(icell) too big
          endif
          
          ! Update total G0 --> TODO(code) this isn't quite right we sould be checking bounds
          if (group_egy(igroup).gt.5.6d0 .and. group_egy(igroup).lt.13.6d0) then 
             total_G0 = total_G0 + (dNp(igroup) * rt_c_cgs(ilevel) * group_egy(igroup) * eV2erg / (1.6d-3))
          end if

          do idim=1,ndim
             dFp(idim,igroup) = &
                  (ddt(icell)*dFpdt(idim,igroup,icell)+dFp(idim,igroup)) &
                  /(1d0+ddt(icell)*(phAbs(igroup)+phSc(igroup)))
          end do
          call reduce_flux(dFp(:,igroup),dNp(igroup)*rt_c_cgs(ilevel))

          do idim=1,ndim
             dUU = ABS(dFp(idim,igroup)-Fp(idim,igroup,icell))           &
                  / (ABS(Fp(idim,igroup,icell))+Np_MIN*rt_c_cgs(ilevel))         &
                  * one_over_Fp_FRAC
             fracMax=MAX(fracMax,dUU)   ! To check if ddt can be increased
             if(dUU .gt. 1d0) then
                dt_rec = 0.9d0 * ddt(icell) / sqrt(2.d0+fracMax)
                code=2 ;   RETURN                     ! ddt(icell) too big
             endif
          end do

       end do

       do igroup=1,nGroups ! -------Momentum transfer from photons to gas:
          mom_fact = ddt(icell) * (phAbs(igroup) + phSc(igroup)) &
               * group_egy_erg(igroup) * one_over_clight

          if(rt_isoPress .and. .not. (rt_isIR .and. igroup==iIR)) then
             ! rt_isoPress: assume f=1, where f is reduced flux.
             fluxMag=sqrt(sum((dFp(:,igroup))**2))
             if(fluxMag .gt. 0d0) then
                mom_fact = mom_fact * dNp(igroup) / fluxMag
             else
                mom_fact = 0d0
             endif
          else
             mom_fact = mom_fact * one_over_rt_c_cgs
          end if

          do idim = 1, ndim
             dmom(idim) = dmom(idim) + dFp(idim,igroup) * mom_fact
          end do
       end do
       dp_gas = dp_gas + dmom * rt_pressBoost      ! update gas momentum

       ! Add absorbed UV/optical energy to IR:----------------------------
       if(rt_isIR) then
          do igroup=iIR+1,nGroups
             dNp(iIR) = dNp(iIR) + dustAbs(igroup) * ddt(icell)          &
                  * dNp(igroup) * group_egy_ratio(igroup)
          end do
       endif
       ! -----------------------------------------------------------------
    endif !if(rt)
#endif

    ! UPDATE TEMPERATURE *************************************************
    !if(c_switch(icell) .and. .not. rt_isTconst .and. .not. rt_T_rad) then
    if(.not. rt_isTconst .and. .not. rt_T_rad) then
       !HKnote: we call prime first so what we can store the correct cooling rates
       saved_cooling_rates = 0.d0
       call all_cooling(TK + (1.d-5*TK), ne, aexp, nElement_dep(1:n_elements), dXion, nCO(icell), total_G0, dust_to_gas_mass_ratio_over_mw, xe, &
                        primary_cosmic_ray_ionization_rate, H2_cosmic_ray_ionization_rate, & 
                        ss_factor, dNp, ilevel, Crate_prime_a, saved_cooling_rates, saved_cooling_rates_names)
       call all_cooling(TK - (1.d-5*TK), ne, aexp, nElement_dep(1:n_elements), dXion, nCO(icell), total_G0, dust_to_gas_mass_ratio_over_mw, xe, &
                        primary_cosmic_ray_ionization_rate, H2_cosmic_ray_ionization_rate, & 
                        ss_factor, dNp, ilevel, Crate_prime_b, saved_cooling_rates, saved_cooling_rates_names)
       saved_cooling_rates = 0.d0
       call all_cooling(TK, ne, aexp, nElement_dep(1:n_elements), dXion, nCO(icell), total_G0, dust_to_gas_mass_ratio_over_mw, xe, &
                        primary_cosmic_ray_ionization_rate, H2_cosmic_ray_ionization_rate, & 
                        ss_factor, dNp, ilevel, Crate, saved_cooling_rates, saved_cooling_rates_names)
       Crate_prime = (Crate_prime_a - Crate_prime_b) / (2.d-5*TK) ! Central difference should be more stable
       dCdT2 = Crate_prime * mu                            ! dC/dT2 = mu * dC/dT

       X_nHkb = 1.d0/(1.5d0 * (rho/mH) * kB)
       rate  = X_nHkb*Crate
       dRate = -X_nHkb*dCdT2                         ! dRate/dT2
                                                     ! 1st order dt constr
       dUU   = ABS(MAX(T2_min_fix, T2(icell)+rate*ddt(icell))-T2(icell))
                                                            ! New T2 value
       dT2   = MAX(T2_min_fix &
                  ,T2(icell)+rate*ddt(icell)/(1.d0-dRate*ddt(icell)))
       dUU   = MAX(dUU, ABS(dT2-T2(icell))) / (T2(icell)+T_MIN) &
                        *one_over_T_FRAC
     
       fracMax=MAX(fracMax,dUU)
       if(dUU .gt. 1.) then                                     ! 10% rule
         !  write(*,*) "Broken Temperature", T2(icell), dT2, ddt(icell)/1.d12, Crate
          dt_rec = 0.9d0 * ddt(icell) / sqrt(2.d0+fracMax)
          code=3 ; RETURN
       endif
       TK=dT2*mu
       if(rt_isTconst) TK=rt_Tconst                         ! Force constant T
    endif

#ifdef RT
    if(rt_isIR) then
       if(kAbs_loc(iIR) .gt. 0d0 .and. .not. rt_T_rad) then
          ! Evolve IR-Dust equilibrium temperature------------------------
          ! Delta (Cv T)= ( c_red/lambda E - c/lambda a T^4)
          !           / ( 1/Delta t + 4 c/lambda/C_v a T^3 + c_red/lambda)
          one_over_C_v = mH*mu*(gamma-1d0) / (rho*kB)
          E_rad = group_egy_erg(iIR) * dNp(iIR)
          dE_T = (rt_c_cgs(ilevel) * E_rad - c_cgs*a_r*TK**4)                    &
               /(1d0/(kAbs_loc(iIR) * Zsolar * rho * ddt(icell))  &
               +4d0*c_cgs * one_over_C_v *a_r*TK**3+rt_c_cgs(ilevel))
          dT2 = dT2 + 1d0/mu * one_over_C_v * dE_T
          dNp(iIR) = dNp(iIR) - dE_T * one_over_egy_IR_erg

          dT2 = max(T2_min_fix,dT2)
          dNp(iIR) = max(dNp(iIR), smallNp)
          ! 10% rule for photon density:
          dUU = ABS(dNp(iIR)-Np(iIR,icell)) / (Np(iIR,icell)+Np_MIN)     &
                                            * one_over_Np_FRAC
          fracMax=MAX(fracMax,dUU)
          if(dUU .gt. 1.) then
             dt_rec = 0.9d0 * ddt(icell) / sqrt(2.d0+fracMax)
             code=4 ;   RETURN
          endif

          dUU   = ABS(dT2-T2(icell)) / (T2(icell)+T_MIN) * one_over_T_FRAC
          fracMax=MAX(fracMax,dUU)
          if(dUU .gt. 1.) then                           ! 10% rule for T2
             dt_rec = 0.9d0 * ddt(icell) / sqrt(2.d0+fracMax)
             code=5 ; RETURN
          endif
          TK=dT2*mu
          if(rt_isTconst) TK=rt_Tconst                         ! Force constant T
          call reduce_flux(dFp(:,iIR),dNp(iIR)*rt_c_cgs(ilevel))
       endif
    endif
#endif

    !/////////////////////////////////////////
    !//           UPDATE MOLECULES          //
    !/////////////////////////////////////////
    dUU = 0.d0
    alpha_H2_loc = 0.d0
    beta_H2_loc = 0.d0
    cr_H2 = 0.d0
    de_H2 = 0.d0
    if (isH2_rtz) then
       xH2_loc = dXion(1,3) / 2.d0 ! Note that we actually store 2 * xH2

       !! Creation !!

       ! Contains formation on dust and via the primordial channel
       alpha_H2_loc = alpha_H2(TK, dust_to_gas_mass_ratio_over_mw, xe, H2_cosmic_ray_ionization_rate, &
                               total_G0, dXion(1,1), dXion(1,2), nElement_dep(1)) 
       cr_H2 = cr_H2 + alpha_H2_loc

       !! Destruction !!

       ! Collisional dissociation
       if (rtz_include_collisional_ionization) then 
         if (elements(2)%atomic_number.gt.0) then
            beta_H2_loc = beta_H2_krome(TK, dXion(1,1)*nElement_dep(1), ne, xH2_loc*nElement_dep(1), dXion(2,1)*nElement_dep(2)) 
         else
            beta_H2_loc = beta_H2_krome(TK, dXion(1,1)*nElement_dep(1), ne, xH2_loc*nElement_dep(1), 0.d0)
         end if
         de_H2 = de_H2 + beta_H2_loc
       end if

       ! Photodissociation
       if (rtz_include_photoionization) then 
         de_H2 = de_H2 + (UV_background_G0 * 5.68d-11)
       end if 

       ! Cosmic ray destruction
       if (rtz_include_cosmic_ray_ionization) then 
         de_H2 = de_H2 + H2_cosmic_ray_ionization_rate
       end if

#ifdef RT
       ! Photodissociation from the local radiation field
       if (rtz_include_photoionization.and.rt_advect) then
          do igroup=1,nGroups
             if (isLW(igroup).eq.1.d0) then
                de_H2 = de_H2 + (dXion(1,3) * SUM(signc(igroup,1,3) * dNp * f_shd))
             else
                de_H2 = de_H2 + (dXion(1,3) * SUM(signc(igroup,1,3) * dNp * f_shd))
             end if  
          end do
       end if
#endif

       ! Update xH2 and store in dXion
       xH2_loc = (cr_H2*ddt(icell) + xH2_loc)/(1.d0+de_H2*ddt(icell))
       dXion(1,3) = 2.d0 * min(max(xH2_loc,x_MIN),0.5d0)

       ! Check for convergence
       dUU = MAX(dUU,ABS((dXion(1,3)-xion(1,3,icell))/(xion(1,3,icell)+x_FM)))
       dUU = dUU * one_over_x_FRAC
       fracMax=MAX(fracMax,dUU)
       if(dUU .gt. 1.d0) then
         !  write(*,*) "Broken H2", TK, dXion(1,3), xion(1,3,icell), ABS((dXion(1,3)-xion(1,3,icell))/(xion(1,3,icell)+x_FM))
          dt_rec = 0.9d0 * ddt(icell) / sqrt(2.d0+fracMax)
          code=6 !TODO(code) update this code for each ion
          RETURN
       end if

      ! Update mu and T
      mu = getMu_RTZ(ne, nElement_dep, dXion)
      TK = dT2 * mu  
      if(rt_isTconst) TK=rt_Tconst                         ! Force constant T 
    end if

#ifdef CO
    cr_CO = 0.d0
    de_CO = 0.d0
    dUU = 0.d0
    if (isCO_rtz) then
       n_C_og = n_CII + nCO(icell)
       n_O_og = n_OI + nCO(icell)
       n_CII = nElement_dep(6) * dXion(6,2)
       n_OI  = nElement_dep(8) * dXion(8,1)
       n_H2  = 0.5d0 * nElement_dep(1) * dXion(1,3)
       x_OI  = n_OI / nElement_dep(1)

       !! Creation !!
       cr_CO = alpha_CO(total_G0, H2_cosmic_ray_ionization_rate, n_CII, n_H2, x_OI)

       !! Destruction !!
       de_CO = beta_CO(total_G0, H2_cosmic_ray_ionization_rate)

       ! Compute the initial guess of new nCO
       nCO_new = (nCO(icell) + cr_CO*ddt(icell)) / (1.d0 + de_CO*ddt(icell))

       ! Tentative update
       delta_CO = nCO_new - nCO(icell)

       ! Enforce positivity floors: don't destroy more CO than exists,
       ! and don't create more than allowed by available CII and OI
       max_delta_CO = min(n_CII - 1.d-30, n_OI - 1.d-30)
       min_delta_CO = -nCO(icell) + 1.d-30       

       ! Clip delta_CO to enforce physical bounds
       delta_CO = max(min(delta_CO, max_delta_CO), min_delta_CO)

       ! Apply updates
       nCO_new  = nCO(icell) + delta_CO
       nCII_new = n_CII - delta_CO
       nOI_new  = n_OI  - delta_CO

       ! Now update the ion fractions for C and O
       tot_C = sum(nElement_dep(6) * dXion(6,1:elements(6)%n_ions)) - n_CII
       tot_C = tot_C + nCII_new
       do iIon = 1, elements(6)%n_ions
          if (iIon.ne.2) then 
             dXion(6,iIon) = nElement_dep(6) * dXion(6,iIon) / tot_C
          else
             dXion(6,iIon) = nCII_new / tot_C
          end if
       end do

       tot_O = sum(nElement_dep(8) * dXion(8,1:elements(8)%n_ions)) - n_OI
       tot_O = tot_O + nOI_new
       do iIon = 1, elements(8)%n_ions
          if (iIon.ne.1) then 
             dXion(8,iIon) = nElement_dep(8) * dXion(8,iIon) / tot_O
          else
             dXion(8,iIon) = nOI_new / tot_O
          end if
       end do

       ! Check for convergence
       dUU = MAX(dUU,ABS((nCO_new-nCO(icell))/(nCO(icell)+x_FM)))
       dUU = dUU * one_over_x_FRAC
       fracMax=MAX(fracMax,dUU)
       if(dUU .gt. 1.d0) then
         !  write(*,*) "Broken CO", TK, nCO_new, nCO, ABS((nCO_new-nCO(icell))/(nCO(icell)+x_FM))
          dt_rec = 0.9d0 * ddt(icell) / sqrt(2.d0+fracMax)
          code=6 !TODO(code) update this code for each ion
          RETURN
       end if

       ! Update species number densities --> need to fix this in case other species fail
       dCO = nCO_new
       dnElement(6) = dnElement(6) - delta_CO
       nElement_dep(6) = nElement_dep(6) - delta_CO
       dnElement(8) = dnElement(8) - delta_CO
       nElement_dep(8) = nElement_dep(8) - delta_CO
    end if
#endif

    !/////////////////////////////////////////
    !//       UPDATE IONIZATION STATES      //
    !/////////////////////////////////////////

    ! Get the effective dust number density
    dust_effective_number_density = nElement_dep(1) * dust_to_gas_mass_ratio_over_mw

    ! Loop over all elements
    do iElement = 1,n_elements
       if (elements(iElement)%atomic_number > 0) then

          ! Get the atomic number
          atomic_number = elements(iElement)%atomic_number

          ! Get the number of ions
          n_ions = elements(iElement)%n_ions

          ! Loop over the ions and get the relevant rates
          saved_rates = 0.d0

          ! Initialize collisional ionization before the loop
          ! recombination rates are 0 for ground state
          saved_rates(1,2) = collisional_ionization(TK, 1, iElement)

          ! Loop over the number of ions
          do iIon = 1,n_ions
             ! Make sure we have the rates of all next ions
             ! previous ions are already computed
             if (iIon.lt.n_ions) then
                saved_rates(iIon+1,1) = recombination(TK, iIon+1, iElement)
                saved_rates(iIon+1,2) = collisional_ionization(TK, iIon+1, iElement)
                saved_rates(iIon+1,3) = dust_recombination(iIon+1, iElement, TK, UV_background_G0, ne)
             end if

             ! Initialize the change to zero
             dUU = 0.d0

             !/////////////////////////
             !//       Creation      //
             !/////////////////////////
             cr = 0.d0

             ! Account for molecular hydrogen
             if (iElement.eq.1 .and. iIon.eq.1 .and. isH2_rtz) then
                ! Note: no factor of 2 needed since is 2*xH2
                cr = cr + (de_H2 * dXion(1,3))
             end if

             ! Recombination of the more excited ionization state
             if (iIon.lt.n_ions) then  
               !  cr = cr + (recombination(TK, iIon+1, iElement) * ne * dXion(iElement,iIon+1))
                cr = cr + (saved_rates(iIon+1,1) * ne * dXion(iElement,iIon+1))
             end if

             ! Collisional ionization of the less excited state
             if (rtz_include_collisional_ionization) then
               if (iIon.gt.1) then 
                  ! cr = cr + (collisional_ionization(TK, iIon-1, iElement) * ne * dXion(iElement,iIon-1))
                  cr = cr + (saved_rates(iIon-1,2) * ne * dXion(iElement,iIon-1))
               end if
             end if

             ! UVB Photoionization of the less excited state
             if (rtz_include_photoionization) then 
               if (iIon.gt.1) then 
                  cr = cr + (HM12_UVB_z(iElement,iIon-1,1) * ss_factor * dXion(iElement,iIon-1))
               end if
             end if

             ! Photoionization by sub-ionizing ISRF --> only impacts lowest ionization states
             if (rtz_include_HM12_UVB) then 
               if (iIon.eq.2) then 
                  cr = cr + (UV_background_G0 * elements(iElement)%G0_photo_rate * dXion(iElement,iIon-1))
               end if
             end if

             ! Cosmic ray ionization of the less excited state
             if (rtz_include_cosmic_ray_ionization) then 
               if (iIon > 1) then 
                  cr = cr + (cosmic_ray_ionization_rates(iElement,iIon-1) * total_cosmic_ray_ionization_rate * dXion(iElement,iIon-1))
               end if
             end if

             ! Cosmic ray ionization of the less excited state from induced UV
             if (rtz_include_cosmic_ray_ionization) then 
               if (iIon.eq.2) then 
                  cr = cr + (cosmic_ray_ionization_rates_induced_UV(iElement) * cosmic_ray_scale_factor * dXion(iElement,iIon-1))
               end if
             end if

             ! Recombination on dust from the more excited state
             if (rtz_include_dust_recombination) then 
               if (iIon.lt.n_ions) then 
                  ! cr = cr + (dust_recombination(iIon+1, iElement, TK, UV_background_G0, ne) * dust_effective_number_density * dXion(iElement,iIon+1))
                  cr = cr + (saved_rates(iIon+1,3) * dust_effective_number_density * dXion(iElement,iIon+1))
               end if
             end if

#ifdef RT
             ! Photoionization of less excited state from the local radiation field
             if (rtz_include_photoionization.and.rt_advect) then
                if (iIon > 1) then 
                   cr = cr + (dXion(iElement,iIon-1) * SUM(signc(:,iElement,iIon-1)*dNp))
                end if
             end if
#endif

             !/////////////////////////
             !//     Destruction     //
             !/////////////////////////
             de = 0.d0

             ! Account for molecular hydrogen
             if (iElement.eq.1 .and. iIon.eq.1 .and. isH2_rtz) then
               de = de + alpha_H2_loc
             end if

             ! Collisional ionization 
             if (rtz_include_collisional_ionization) then
               if (iIon .lt. n_ions) then 
                  ! de = de + (collisional_ionization(TK, iIon, iElement) * ne)
                  de = de + (saved_rates(iIon,2) * ne)
               end if
             end if

             ! UVB Photoionization
             if (rtz_include_HM12_UVB) then 
               if (iIon .lt. n_ions) then 
                  de = de + (HM12_UVB_z(iElement,iIon,1) * ss_factor)
               end if
             end if

             ! Photoionization by sub-ionizing ISRF --> only impacts lowest ionization states
             if (rtz_include_photoionization) then 
               if (iIon .eq. 1) then
                  de = de + (UV_background_G0 * elements(iElement)%G0_photo_rate)
               end if
             end if

             ! Recombination
             if (iIon .gt. 1) then 
               !  de = de + (recombination(TK, iIon, iElement) * ne)
                de = de + (saved_rates(iIon,1) * ne)
             end if 

             ! Cosmic ray ionization
             if (rtz_include_cosmic_ray_ionization) then 
               if (iIon .lt. n_ions) then
                  de = de + (cosmic_ray_ionization_rates(iElement,iIon) * total_cosmic_ray_ionization_rate)
               end if
             end if

             ! Cosmic ray ionization from induced UV
             if (rtz_include_cosmic_ray_ionization) then 
               if (iIon .eq. 1) then 
                  de = de + (cosmic_ray_ionization_rates_induced_UV(iElement) * cosmic_ray_scale_factor)
               end if
             end if
            
             ! Recombination on dust
             if (rtz_include_dust_recombination) then 
               if (iIon .gt. 1) then 
                  ! de = de + (dust_recombination(iIon, iElement, TK, UV_background_G0, ne) * dust_effective_number_density)
                  de = de + (saved_rates(iIon,3) * dust_effective_number_density)
               end if
             end if

#ifdef RT
             ! Photoionization  from the local radiation field
             if (rtz_include_photoionization.and.rt_advect) then
                if (iIon .lt. n_ions) then 
                   de = de + (dXion(iElement,iIon) * SUM(signc(:,iElement,iIon)*dNp))
                end if
             end if
#endif
             !/////////////////////////
             !//   Charge Transfer   //
             !/////////////////////////
             ! Note, this was split off due to cross species
             ! coupling. Saves us an extra double loop
             if (rtz_include_charge_exchange) then
               if (iElement.eq.1) then !If element is hydrogen
                  ! Loop over all other elements
                  do i_other_Element = 2,n_elements
                     if (elements(i_other_Element)%atomic_number > 0) then
                        ! Loop over all other ionization states
                        do i_other_Ion = 1,elements(i_other_Element)%n_ions
                           ! Get the number density of the other ion
                           paired_ion_number_density = nElement_dep(i_other_Element) * dXion(i_other_Element,i_other_Ion)

                           if (iIon.eq.1) then !H
                              cr = cr + (charge_transfer_ionization(i_other_Ion,i_other_Element,TK) * dXion(iElement,iIon+1) * paired_ion_number_density) ! Example:  O + H+ => O+ + H
                              de = de + (charge_transfer_recombination(i_other_Ion,i_other_Element,TK) * paired_ion_number_density) ! Example:  O+ + H => O + H+
                           else !H+
                              de = de + (charge_transfer_ionization(i_other_Ion,i_other_Element,TK) * paired_ion_number_density) ! Example:  O + H+ => O+ + H
                              cr = cr + (charge_transfer_recombination(i_other_Ion,i_other_Element,TK) * dXion(iElement,iIon-1) * paired_ion_number_density) ! Example:  O+ + H => O + H+
                           end if

                        end do !end loop over ionization states
                     end if
                  end do ! end loop over other elements
               else !All other elements
                  HI_number_density = dXion(1,1) * nElement_dep(1)
                  HII_number_density = dXion(1,2) * nElement_dep(1)

                  if (iIon.gt.1) then
                     ! Ionization from less excited state
                     cr = cr + (charge_transfer_ionization(iIon-1,iElement,TK) * dXion(iElement,iIon-1) * HII_number_density) ! Example:  O + H+ => O+ + H   

                     ! Charge exchange recombination 
                     de = de + (charge_transfer_recombination(iIon,iElement,TK) * HI_number_density) ! Example:  O+ + H => O + H+
                  end if

                  if (iIon.lt.n_ions) then
                     ! Charge exchange ionization
                     de = de + (charge_transfer_ionization(iIon,iElement,TK) * HII_number_density) ! Example:  O + H+ => O+ + H   

                     ! Charge exchange recombination from the more excited state
                     cr = cr + (charge_transfer_recombination(iIon+1,iElement,TK) * dXion(iElement,iIon+1) * HI_number_density) ! Example:  O+ + H => O + H+
                  end if
               end if
             end if

             !/////////////////////////
             !//       Update        //
             !/////////////////////////
             dXion(iElement,iIon) = (cr*ddt(icell) + dXion(iElement,iIon))/(1.d0 + de*ddt(icell))
             dXion(iElement,iIon) = min(max(dXion(iElement,iIon),x_MIN),1.d0)

             ! Get the new electron fraction
             ne = getNe(dXion, nElement_dep(:))
             xe = ne / nElement_dep(1)
             phi_s = secondary_cr_rates(xe)
             total_cosmic_ray_ionization_rate = primary_cosmic_ray_ionization_rate * (1.d0 + phi_s)

             ! Update mu and T --> only H and He (others don't matter)
             if (iElement.lt.3) then 
                mu = getMu_RTZ(ne, nElement_dep, dXion)
                TK = dT2 * mu  
                if(rt_isTconst) TK=rt_Tconst                         ! Force constant T 
             end if

             ! Check for convergence -- Fractional change in ion
             dUU = MAX(dUU,ABS((dXion(iElement,iIon)-xion(iElement,iIon,icell))/(xion(iElement,iIon,icell)+x_FM)))
             dUU = dUU * one_over_x_FRAC
             fracMax=MAX(fracMax,dUU)
             if(dUU .gt. 1.) then
               !  write(*,*) "Broken element/ion", iElement, iIon, dXion(iElement,iIon), xion(iElement,iIon,icell)
                dt_rec = 0.9d0 * ddt(icell) / sqrt(2.d0+fracMax)
                code=6 !TODO(code) update this code for each ion
                RETURN
             end if

             ! Check for convergence -- Fractional change in electrons
             dUU=ABS((ne-neInit)) / (neInit+x_FM) * one_over_x_FRAC
             fracMax=MAX(fracMax,dUU)
             if(dUU .gt. 1.) then
               !  write(*,*) "Broken electron", TK, ABS((ne-neInit)) / (neInit+x_FM)
                dt_rec = 0.9d0 * ddt(icell) / sqrt(2.d0+fracMax)
                code=8
                RETURN
             endif

          end do ! END ION LOOP

          ! REDUCE SO THAT IONIZATION FRACTIONS SUM TO 1
          ion_fracs = elements(iElement)%n_ions + elements(iElement)%n_mol
          current_mass_frac = sum(dXion(iElement,1:ion_fracs))
          do iIon=1,ion_fracs
             dXion(iElement,iIon) = dXion(iElement,iIon) + ((1.d0 - current_mass_frac) * (dXion(iElement,iIon) / current_mass_frac))
          end do
       end if
    end do ! END ELEMENT LOOP

    ! UPDATE CONSTANT T WITH NEW MU **************************************
    if(rt_isTconst)then
       mu = getMu_RTZ(ne, nElement_dep, dXion)
       dT2 = rt_Tconst/mu
    endif

    ! CLEAN UP AND RETURN ************************************************
    dT2 = dT2-T2(icell) ; dXion(:,:) = dXion(:,:)-xion(:,:,icell)
    dnElement(:) = dnElement(:) - nElement(:,icell)
#ifdef CO
    dCO = dCO - nCO(icell)
#endif
#ifdef RT
    dNp(:) = dNp(:)-Np(:,icell) ; dFp(:,:) = dFp(:,:)-Fp(:,:,icell)
    dp_gas(:)= dp_gas(:)-p_gas(:,icell)
#endif
    ! Now the dUs are really changes, not new values
    ! Update the timestep for the next iteration:
   !  dt_rec = 0.5d0 * ddt(icell) / ((0.01d0 + fracMax)**0.5d0)
    dt_rec = 0.9d0 * ddt(icell) / ((0.07d0 + fracMax)**0.3d0)
    dt_rec = min(dt_rec,1E12 * min(TK/100.0,1.0) * min((1.0/nElement_dep(1)),1.0) * min(sqrt(1.0/UV_background_G0),1.0))
    dt_rec = min(dt_rec,rtz_max_cool_timestep)
    dt_ok = .true.
    code=0

    ! Save T and mu
    TK_to_save(icell)=TK
    mu_to_save(icell)=mu

  END SUBROUTINE rtz_cool_step

END SUBROUTINE rtz_solve_cooling

!************************************************************************
FUNCTION getNe(xion,nion) result(ne)
  !Returns the electron number density by looping over all elements 
  !and summing their contributions
  implicit none

  real(dp), intent(in)::xion(1:n_elements,1:n_elements)
  real(dp), intent(in)::nion(1:n_elements)
  real(dp)::ne
  integer::iIons, iElement, n_ions

  ne = 0.d0

  !Loop over all elements
  do iElement=1,n_elements
     if (elements(iElement)%atomic_number .gt. 0) then
        !Get the number of ions
        n_ions = elements(iElement)%n_ions

        !Loop over all ions 
        !Start loop at 2, no electrons in the ground state
        do iIons=2,n_ions
           ne = ne + (nion(iElement) * xion(iElement,iIons) * real(iIons - 1, dp)) 
        end do
     end if
  end do

END FUNCTION getNe

FUNCTION dust_to_gas_scale_RR14(log10_O_over_H) result(ratio)
   ! This returns the dust-to-metal ratio relative to the local value
   !------------------------------------
   ! Zsolar: metallicity in solar units
   ! D2Z_solar ~ D2Z / 0.3
   !------------------------------------
   ! based on Remy-Ruyer, Madden, Galliano et al. (2014)
   ! https://arxiv.org/pdf/1312.3442.pdfR
   ! Broken power law with X_{CO,Z} case (Table 1)
   implicit none

   real(dp),intent(in)::log10_O_over_H
   real(dp)::ratio
   real(dp)::a,alphaH,b,alphaL,xt,x,Xsun
   real(dp)::G2D,G2D_sol,y

   ! y = log (G/D)
   ! x = 12 + log10(O/H)
   ! Xsun = 8.69
   ! 
   ! y = a + alphaH * (Xsun - x) for x>xt
   ! y = b + alphaL * (Xsun - x) for x<=xt
   a      = 2.21d0
   alphaH = 1.00d0 ! MW case
   b      = 0.96d0 ! 0.68
   alphaL = 3.10d0 ! 3.08 
   xt     = 8.10d0 ! 7.96
   Xsun   = 8.69d0
   x = Xsun - log10_O_over_H
   x = max(log10_O_over_H,5.d0) ! Mild extrapolation

   if (log10_O_over_H>xt)then
      y = a + alphaH * (Xsun - log10_O_over_H)
   else
      y = b + alphaL * (Xsun - log10_O_over_H)
   endif

   G2D = 10.d0**y
   G2D_sol = 10.d0**a

   ratio = max(min(G2D_sol / G2D, 1.d0), 0.d0)

END FUNCTION dust_to_gas_scale_RR14

FUNCTION getMu_RTZ(ne, element_number_densities, element_ion_fractions) result(mu)
   implicit none
   real(dp), intent(in):: ne
   real(dp), intent(in):: element_number_densities(27)
   real(dp), intent(in):: element_ion_fractions(27,27)
   real(dp):: mu
   real(dp):: m_bar, n_hat

   integer:: i, j

   m_bar = 0.d0
   n_hat = 0.d0

   do i=1,n_elements
      if (elements(i)%atomic_number.gt.0) then
         do j=1,elements(i)%n_ions
            m_bar = m_bar + (element_number_densities(i) * element_ion_fractions(i,j) * elements(i)%atomic_mass)
            n_hat = n_hat + element_number_densities(i) * element_ion_fractions(i,j)
         end do
      end if
   end do

   ! Include electrons
   n_hat = n_hat + ne

   ! Include contribution from H2
   if (isH2_rtz) then
      m_bar = m_bar + (element_number_densities(1) * element_ion_fractions(1,3) * elements(1)%atomic_mass)
      n_hat = n_hat + (0.5d0 * element_number_densities(1) * element_ion_fractions(1,3))
   end if

   mu = m_bar / n_hat

END FUNCTION getMu_RTZ

FUNCTION get_rho_rtz(element_number_densities) result(rho)
   use constants, only: mH
   implicit none
   real(dp), intent(in):: element_number_densities(27)
   real(dp):: rho

   integer:: i

   rho = 0.d0

   do i=1,n_elements
      if (elements(i)%atomic_number.lt.1) then
         cycle
      end if
      rho = rho + (element_number_densities(i) * elements(i)%atomic_mass) ! gives amu/cm^3
   end do

   rho = mH * rho ! gives g/cm^3
END FUNCTION get_rho_rtz

FUNCTION get_n_rtz(element_number_densities, ne) result(rho_n)
   implicit none
   real(dp), intent(in):: element_number_densities(27)
   real(dp), intent(in):: ne
   real(dp):: rho_n

   integer:: i

   rho_n = 0.d0

   do i=1,n_elements
      if (elements(i)%atomic_number.lt.1) then
         cycle
      end if
      rho_n = rho_n + element_number_densities(i)  ! gives 1/cm^3
   end do

   rho_n = rho_n + ne
END FUNCTION get_n_rtz

#ifdef RT
SUBROUTINE reduce_flux(Fp, cNp)
! Make sure the reduced photon flux is less than one
!------------------------------------------------------------------------
  use rt_parameters
  implicit none
  real(dp),dimension(ndim):: Fp
  real(dp):: cNp, fred
!------------------------------------------------------------------------
  fred = sqrt(sum(Fp**2))/cNp
  if(fred .gt. 1d0) Fp = Fp/fred
END SUBROUTINE reduce_flux
#endif

END MODULE rtz_cooling_module

!************************************************************************
SUBROUTINE rtz_updateRTGroups_CoolConstants(ilevel)
  ! Update photon group cooling and heating constants, to reflect an update
  ! in rt_c_cgs and in the cross-sections and energies in the groups.
  !------------------------------------------------------------------------
  use rt_parameters
  implicit none
#ifdef RT
  !------------------------------------------------------------------------
  integer, intent(in)::ilevel
  integer::iP, iE, iI
  !------------------------------------------------------------------------
  signc(:,:,:) = group_csn*rt_c_cgs(ilevel)        ! [cm3 s-1]
  sigec(:,:,:) = group_cse*rt_c_cgs(ilevel)        ! [cm3 s-1]

  !Photoheating rates for photons on ions
  do iP = 1,nGroups
     do iE = 1,n_elements
        if (elements(iE)%atomic_number.gt.0) then 
           do iI = 1,elements(iE)%n_ions-1
              PHrate(iP,iE,iI) =  eV2erg * &    ! See eq (19) in Aubert(08)
                 (sigec(iP,iE,iI) * group_egy(iP)  &
                 -signc(iP,iE,iI)*ionEvs(iE,iI))
              PHrate(iP,iE,iI) = max(PHrate(iP,iE,iI),0d0)!Heating>0
           end do
        end if
      end do ! End element loop
   end do ! End group loop
#endif
END SUBROUTINE rtz_updateRTGroups_CoolConstants