!! Cooling from frig version (Audit & Hennebelle 2005)
!! TC: some of these processes depend on the composition of the gas which
!!     can be obtained from RT
!! Authors: Patrick Hennebelle, Benjamin Godard
!=======================================================================
subroutine solve_cooling_frig(nH,T2,zsolar,boost,dt,deltaT2,ncell)
!=======================================================================
  implicit none
  ! BRIDGE FUNCTION WITH SAME INTERFACE AS solve_cooling
  ! Input/output variables to this function
  ! nH - hydrogen number density in PHYSICAL units
  ! T2 - temperature / mu in PHYSICAL units
  ! zsolar - Metallicity in solar units (Zphys / 0.02)
  ! boost - raditation boost - exp(-100*nH) if self_shielding=True
  ! dt - cooling timestep in seconds
  ! deltaT2 - temperature change in K/mu (??)
  ! ncell - number of elements in the vector
  integer::ncell
  real(kind=8)::dt
  real(kind=8),dimension(1:ncell)::nH,T2,deltaT2,zsolar,boost
  ! Input/output variables to analytic function calc_temp 
  real(kind=8)::NN,TT, dt_tot_unicode
  ! Temporary variables
  integer::i
  real(kind=8)::TT_ini, mu
  ! Units
  real(kind=8) :: scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  ! HARD-CODED mu TO MAKE TEMPERATURE AGREE WITH HENNEBELLE CODE
  ! TC: this is ugly, we can do better
  ! if RT: get mu from rt_cooling
  mu = 1.4
  scale_T2 = scale_T2 * mu
  do i=1,ncell
     NN = nH(i) ! NOTE!! THE CODE BELOW ASSUMES scale_nH=1 !!
                ! SO WE LEAVE THIS AS IT IS TO KEEP UNITS CONSISTENCY
     TT = T2(i) / scale_T2
     TT_ini = TT
     dt_tot_unicode = dt / scale_t
     call calc_temp(NN,TT,mu,dt_tot_unicode)
     deltaT2(i) = (TT - TT_ini) * scale_T2
  end do
end subroutine solve_cooling_frig
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_temp(NN,TT,mu,dt_tot_unicode)
    use amr_parameters
    use hydro_commons
    use constants, only:kB
    implicit none

    integer :: n,i,j,k,idim, iter,ii
    real(dp) :: dt, dt_tot, temps, dt_max
    real(dp) :: rho,temp,dt_tot_unicode
    real(dp) :: alpha_ct,mu
    real(dp) :: NN,TT, TTold, ref,ref2,dRefdT, eps, vardt,varrel, dTemp
    real(dp) :: rhoutot2
    real(dp) :: scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2

    ! cgs units are used here
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
    scale_T2 = scale_T2 * mu

    if( TT .le. 0.) then
        TT = 50. / scale_T2
        return
    endif

    vardt = 10.**(1./10.); varrel = 0.2

    dt_tot = dt_tot_unicode * scale_t
    TT     = TT * scale_T2

    if (NN .le. smallr) then
        if( NN .le. 0)  write(*,*) 'prob dens',NN
        NN = smallr  
    endif

    alpha_ct = NN*kB/(gamma-1.)

    ! eps - a small offset of T to find gradient in T
    eps = 1d-5

    iter  = 0 ; temps = 0.
    do while ( temps < dt_tot)
        if (TT .lt.0) then
            write(*,*) 'prob Temp',TT, NN
            NN = max(NN,smallr)
            TT = min(4000./NN,8000.)  !2.*4000. / NN
        endif

        TTold = TT

        ! Calculate cooling rate
        !NN is assumed to be in cc and TT in Kelvin
        if (TT < 10035.d0) then
            call cooling_low(TT,NN,ref)
            call cooling_low(TT*(1d0+eps),NN,ref2)
        else
            call cooling_high(TT,NN,ref)
            call cooling_high(TT*(1d0+eps),NN,ref2)
        end if
        
        ! dT = T*(1+eps)-T = eps*T
        dRefdT = (ref2-ref)/(TT*eps)

        ! TODO STG - COPY THIS FUNCTION UP TO HERE, USE ref, drefdT TO 
        !            REPLACE rt_cmp_metals SOMEHOW

        if (iter == 0) then
            if (dRefDT .ne. 0.) then
                dt = abs(1.0E-1 * alpha_ct/dRefDT)
            else
                dt = 1.0E-1 * dt_tot
            endif
            dt_max = dt_tot - temps
            if (dt > 0.7*dt_max) dt = dt_max*(1.+1.0E-12)
        endif

        dTemp = ref/(alpha_ct/dt - dRefdT)

        eps = abs(dTemp/TT)
        if (eps > 0.2) dTemp = 0.2*TTold*dTemp/abs(dTemp)

        TT = TTold + dTemp
        if (TT < 0.) then
            write(*,*) 'Temperature negative !!!'
            write(*,*) 'TTold,TT   = ',TTold,TT
            write(*,*) 'rho   = ',rho
            TT = 100.  !*kelvin
        endif
        iter = iter + 1
        temps = temps + dt
        dt = vardt*varrel*dt/Max(vardt*eps, varrel)
        dt_max = dt_tot - temps
        if (dt > 0.7*dt_max) dt = dt_max*(1.+1.0E-12)
    enddo

    !!now convert temperature in code units
    TT = TT / scale_T2

end subroutine calc_temp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cooling_high(T,n,ref)
    use amr_parameters
    implicit none

    real(dp) :: T,n,ref
    real(dp) :: cold,hot,logT

    ! cooling rate based on Dopita and Sutherland

    logT=log10(T)

    if (logT .LT. 4.0) then
        cold=0.1343*logT**3-1.3906*logT**2+5.1554*logT-31.967
    else if (logT .LT. 4.25) then
        cold=12.64*logT-75.56
    else if (logT .LT. 4.35) then
        cold=-0.3*logT-20.565
    else if (logT .LT. 4.9) then
        cold=1.745*logT-29.463
    else if (logT .LT. 5.4) then
        cold=-20.9125
    else if (logT .LT. 5.9) then
        cold=-1.795*logT-11.219
    else if (logT .LT. 6.2) then
        cold=-21.8095
    else if (logT .LT. 6.7) then
        cold=-1.261*logT-13.991
    else
        cold=-22.44
    endif

    cold=-1.0*10.0**(cold)

    ref=(n**2)*(cold)

end subroutine cooling_high
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cooling_low(T,n,ref)
  use amr_parameters
  implicit none

  real(dp) :: T,n,ref
  real(dp) :: cold,hot,cold_cII,cold_o,cold_h,cold_cII_m,cold_o_m,cold_rec
  real(dp) :: param,G0,epsilon,bet,x,x_ana,ne   ! x is the ionisation rate

  ! cooling and heating function computed from the cooling of 
  ! chemical elements

  ! Carbon abondance 3.5 10-4, depletion 0.4

  !!! We compute the ionisation
  !!! We assume that if x is over 1.d-4 then it is dominated by oxygen
  !!! and that the minimal value is given by the carbon and is 
  !!! 3.5 1.d-4 * depletion * density
  
  !!! For the electrons due to hydrogen we take the formula
  !!! from  Wolfire et al. 2003 appendice C2.
  !!! The ionization rate is set to 1.d-16 G0'=GO/1.7
  !!! Z'd = 1 et phipah=0.5

  ne = 2.4d-3*((T/100d0)**0.25d0)/0.5d0 ! formula C15 of Wolfire et al. 2003

  ! Analytic ionisation in absence of photoionisation
  x_ana = ne / N   ! ionisation
  x_ana = min(x_ana,0.1d0)
  x_ana = max(x_ana,3.5d-4*0.4d0)
  x = x_ana ! (Different variables in case we need photoionised values later)

 ! NOTE - HERE WE USE THE NON-PHOTOIONISED RATES AS THIS MIGHT 
 !        BE TOO HIGH AT x=1
  cold_cII =  92. * 1.38E-16 * 2. * (2.8E-7* ((T/100.)**(-0.5))*x_ana + 8.E-10*((T/100.)**(0.07))) &
       * 3.5E-4 * 0.4 * exp(-92./ T)

  ! oxygen-prompted cooling
  ! abundance 8.6 10-4 depletion 0.8
  cold_o = 1.E-26 * sqrt(T) * (24. * exp(-228./ T) + 7. * exp(-326./ T) )
  ! take oxygen abundance into account
  cold_o = cold_o * 4.5E-4

  ! Hydrogen-prompted cooling
  ! formula from Spitzer 1978
  ! NOTE - RT function handles hydrogen cooling out of equilibrium
  cold_h = 0d0
  if (.not. rt) then
     cold_h = 7.3E-19 * x * exp(-118400./ T )
  endif

  ! cooling from metastables metal lines
  ! formulas from Hollenbach and McKee 1989 (ApJ 342, 306)

  ! Ionised carbon, 1 transition 2P 4P
  ! weight is 1
  ! 2P->4P :
  ! The excitation coefficients depends on the temperature above 10^4 K
  ! les expressions des coefficients d'excitation ont une dependance
  ! abondance 3.5 d-4 depletion 0.4

  cold_cII_m = 6.2d4 * 1.38d-16 * 1.d0 * &    !transition 2P->4P
             ( 2.3d-8* (T/10000.)**(-0.5) * x + 1.d-12 ) *exp(-6.2d4 / T) &
             * 3.5d-4 * 0.4

  if ( T .le. 1.d4 ) then
     cold_o_m = 2.3d4 * 1.38d-16 / 3.d0 * &
             ( 5.1d-9 * (T/10000.)**(0.57) * x + 1.d-12) *exp(-2.3d4/T)
  
     cold_o_m = cold_o_m + &
                4.9d4 * 1.38d-16 / 3.d0  * &
             ( 2.5d-9 * (T/10000.)**(0.57) * x + 1.d-12) *exp(-4.9d4/T)
  
     cold_o_m = cold_o_m + &
                2.6d4 * 1.38d-16 * 1.d0  * &
             ( 5.2d-9 * (T/10000.)**(0.57) * x + 1.d-12) *exp(-2.6d4/T)
  else
     cold_o_m = 2.3d4 * 1.38d-16 / 3.d0 * &
             ( 5.1d-9 * (T/10000.)**(0.17) * x + 1.d-12) *exp(-2.3d4/T)
  
     cold_o_m = cold_o_m + &
                4.9d4 * 1.38d-16 / 3.d0  * &
             ( 2.5d-9 * (T/10000.)**(0.13) * x + 1.d-12) *exp(-4.9d4/T)

     cold_o_m = cold_o_m + &
                2.6d4 * 1.38d-16 * 1.d0  * &
             ( 5.2d-9 * (T/10000.)**(0.15) * x + 1.d-12) *exp(-2.6d4/T)
  endif

  ! oxigen abundance 
  cold_o_m = cold_o_m * 4.5d-4

  ! sum of the cooling terms
  cold = cold_cII  + cold_h  + cold_o  + cold_o_m +  cold_cII_m

  !!!! Computation of the heating term
  !!! Heating on grains is taken into account
  !!! formula 1 et 2  of Wolfire et al. 1995

  !!!! G0 is the UV flux compared to the one given by Habing et Draine
  G0 = 1./1.7

  param = G0 * sqrt(T)/(n*x)
  epsilon = 4.9E-2 / (1. + (param/1925.)**0.73)
  epsilon  = epsilon + 3.7E-2 * (T/1.E4)**0.7 / (1. + (param/5.E3) )

  hot = 1.E-24 * epsilon

  ! for a UV flux of G0/1.7
  hot = hot * G0

  ! recombination cooling on positively charged grains
  bet = 0.74/(T**0.068)
  cold_rec = 4.65E-30*(T**0.94)*(param**bet)*x

  ref = hot*n - (n**2)*(cold + cold_rec)

end subroutine cooling_low