!! Cooling from frig version (Audit & Hennebelle 2005)
!! solve_cooling_ism is used if there is no RT
!! Authors: Patrick Hennebelle, Benjamin Godard
!=======================================================================
subroutine solve_cooling_ism(nH,T2,dt,deltaT2,ncell)
!=======================================================================
  use amr_parameters, only:mu_gas
  implicit none
  ! BRIDGE FUNCTION WITH SAME INTERFACE AS solve_cooling
  ! nH - hydrogen number density in PHYSICAL units
  ! T2 - temperature / mu in PHYSICAL units
  ! dt - cooling timestep in seconds
  ! deltaT2 - temperature change in K/mu (??)
  ! ncell - number of elements in the vector
  integer::ncell
  real(kind=8)::dt
  real(kind=8),dimension(1:ncell)::nH,T2,deltaT2
  real(kind=8)::NN,TT   ! Input/output variables to analytic function calc_temp
  real(kind=8)::TT_ini,mu
  integer::i
  ! mu = 1.4 for Hennebelle code
  mu = mu_gas    ! molecular weight
  do i=1,ncell
     NN = nH(i)      ! H/cc
     TT = T2(i) * mu ! K
     TT_ini = TT
     call calc_temp(NN,TT,dt)
     deltaT2(i) = (TT - TT_ini) / mu
  end do
end subroutine solve_cooling_ism
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_temp(NN,TT,dt_tot)
  ! NN = number densisty collision partners (hydrogen) in H/cm3
  ! TT = temperature in Kelvin (T2*mu)
  ! dt_tot = timestep in s
  use amr_parameters, only:dp
  use hydro_parameters, only:gamma,smallr
  use constants, only:kB
  implicit none

  real(dp) :: NN,TT,mu,dt_tot
  integer  :: iter
  real(dp) :: dt, time, dt_max, vardt
  real(dp) :: TTold, dTemp, eps
  real(dp) :: ref, ref2, varrel
  real(dp) :: dRefdT, alpha_ct

  if( TT .le. 0) then
     TT = 50d0
     return
  endif

  if (NN .le. smallr) then
     if( NN .le. 0)  write(*,*) 'WARNING: problem density in calc_temp',NN
     NN = smallr
  endif

  vardt = 10d0**(0.1d0)
  varrel = 0.2d0

  alpha_ct = NN*kB/(gamma-1d0)

  ! eps - a small offset of T to find gradient in T
  eps = 1d-5

  iter  = 0
  time = 0
  do while (time < dt_tot)
     if (TT .lt. 0) then
        write(*,*) 'WARNING: problem temperature in calc_temp',TT,NN
        NN = max(NN,smallr)
        TT = min(4000d0/NN,8000d0)
     endif

     TTold = TT

     ! Calculate cooling rate
     if (TT < 10035d0) then
        call cooling_low(TT,NN,ref)
        call cooling_low(TT*(1d0+eps),NN,ref2)
     else
        call cooling_high(TT,NN,ref)
        call cooling_high(TT*(1d0+eps),NN,ref2)
     end if

     ! dT = T*(1+eps)-T = eps*T
     dRefdT = (ref2-ref)/(TT*eps)

     if (iter == 0) then
        if (dRefDT .ne. 0) then
           dt = abs(1d-1 * alpha_ct/dRefDT)
        else
           dt = 1d-1 * dt_tot
        endif
        dt_max = dt_tot - time
        if (dt > 0.7d0*dt_max) dt = dt_max*(1.+1d-12)
     endif

     dTemp = ref/(alpha_ct/dt - dRefdT)
     eps = abs(dTemp/TT)
     if (eps > 0.2d0) dTemp = 0.2d0*TTold*dTemp/abs(dTemp)

     TT = TTold + dTemp
     if (TT < 0) then
        write(*,*) 'WARNING: Temperature negative in calc_temp! TTold,TT = ',TTold,TT
        TT = 100
     endif
     iter = iter + 1
     time = time + dt
     dt = vardt*varrel*dt/max(vardt*eps, varrel)
     dt_max = dt_tot - time
     if (dt > 0.7d0*dt_max) dt = dt_max*(1d0+1d-12)
  enddo

end subroutine calc_temp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cooling_high(T,n,ref)
  ! T = physical temperature in K
  ! n = number of Hydrogen atoms for collision / cm3
  ! ref = cooling rate
  use amr_parameters, only:dp
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

  cold=-1*10d0**(cold)

  hot = 0
  ref= hot*n + (n**2)*(cold)

end subroutine cooling_high
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cooling_low(T,n,ref)
  ! T = physical temperature in K
  ! n = number of Hydrogen atoms for collision / cm3
  ! ref = cooling rate
  ! TC: precision errors are significant here...
  use amr_parameters, only:dp,rt
  implicit none

  real(dp) :: T,n,ref
  real(dp) :: cold,hot,cold_cII,cold_o,cold_h,cold_cII_m,cold_o_m,cold_rec
  real(dp) :: param,G0,epsilon,bet,x,ne   ! x is the ionisation rate

  ! Cooling and heating function computed from the cooling of
  ! chemical elements

  ! Carbon abundance 3.5d-4, depletion 0.4

  !!! We compute the ionisation
  !!! We assume that if x is over 1.d-4 then it is dominated by oxygen
  !!! and that the minimal value is given by the carbon and is
  !!! 3.5d-4 * depletion * density

  !!! For the electrons due to hydrogen we take the formula
  !!! from  Wolfire et al. 2003 appendix C2.
  !!! The ionization rate is set to 1d-16 G0'=GO/1.7
  !!! Z'd = 1 et phipah=0.5

  ne = 2.4d-3*((T/100d0)**0.25d0)/0.5d0 ! formula C15 of Wolfire et al. 2003

  ! Analytic ionisation in absence of photoionisation
  x = ne / n   ! ionisation
  x = min(x,0.1d0)
  x = max(x,3.5d-4*0.4d0)

 ! NOTE - HERE WE USE THE NON-PHOTOIONISED RATES AS THIS MIGHT
 !        BE TOO HIGH AT x=1
  cold_cII = 92 * 1.38d-16 * 2 * (2.8d-7* ((T/100)**(-0.5d0))*x + 8d-10*((T/100)**(0.07d0))) &
       * 3.5d-4 * 0.4d0 * exp(-92/T)

  ! Oxygen-prompted cooling
  ! abundance 8.6 10-4 depletion 0.8
  cold_o = 1d-26 * sqrt(T) * (24 * exp(-228/T) + 7 * exp(-326/T))
  ! take oxygen abundance into account
  cold_o = cold_o * 4.5d-4

  ! Hydrogen-prompted cooling
  ! Formula from Spitzer 1978
  ! NOTE - RT function handles hydrogen cooling out of equilibrium
  cold_h = 0
  if (.not. rt) then
     cold_h = 7.3d-19 * x * exp(-118400/T)
  endif

  ! Cooling from metastables metal lines
  ! Formulas from Hollenbach and McKee 1989 (ApJ 342, 306)

  ! Ionised carbon, 1 transition 2P->4P, weight is 1
  ! The excitation coefficients depend on the temperature when T > 10^4 K
  ! abundance 3.5d-4 depletion 0.4

  cold_cII_m = 6.2d4 * 1.38d-16 * 1d0 * &    !transition 2P->4P
             ( 2.3d-8* (T/10000)**(-0.5) * x + 1d-12) *exp(-6.2d4/T) &
             * 3.5d-4 * 0.4

  if (T .le. 1d4) then
     cold_o_m = 2.3d4 * 1.38d-16 / 3d0 * &
             ( 5.1d-9 * (T/10000)**(0.57) * x + 1d-12) *exp(-2.3d4/T)

     cold_o_m = cold_o_m + &
                4.9d4 * 1.38d-16 / 3d0  * &
             ( 2.5d-9 * (T/10000)**(0.57) * x + 1d-12) *exp(-4.9d4/T)

     cold_o_m = cold_o_m + &
                2.6d4 * 1.38d-16 * 1d0  * &
             ( 5.2d-9 * (T/10000)**(0.57) * x + 1d-12) *exp(-2.6d4/T)
  else
     cold_o_m = 2.3d4 * 1.38d-16 / 3d0 * &
             ( 5.1d-9 * (T/10000)**(0.17) * x + 1d-12) *exp(-2.3d4/T)

     cold_o_m = cold_o_m + &
                4.9d4 * 1.38d-16 / 3d0  * &
             ( 2.5d-9 * (T/10000)**(0.13) * x + 1d-12) *exp(-4.9d4/T)

     cold_o_m = cold_o_m + &
                2.6d4 * 1.38d-16 * 1d0  * &
             ( 5.2d-9 * (T/10000)**(0.15) * x + 1d-12) *exp(-2.6d4/T)
  endif

  ! Oxygen abundance
  cold_o_m = cold_o_m * 4.5d-4

  ! Sum of the cooling terms
  cold = cold_cII  + cold_h  + cold_o  + cold_o_m +  cold_cII_m

  !!! Computation of the heating term
  !!! Heating on grains is taken into account
  !!! formula 1 et 2 of Wolfire et al. 1995

  !!!! G0 is the UV flux compared to the one given by Habing et Draine
  G0 = 1d0/1.7d0

  param = G0 * sqrt(T)/(n*x)
  epsilon = 4.9d-2 / (1 + (param/1925)**0.73)
  epsilon  = epsilon + 3.7d-2 * (T/1d4)**0.7 / (1 + (param/5d3) )

  hot = 1d-24 * epsilon

  ! for a UV flux of G0/1.7
  hot = hot * G0

  ! recombination cooling on positively charged grains
  bet = 0.74/(T**0.068)
  cold_rec = 4.65d-30*(T**0.94)*(param**bet)*x

  ref = hot*n - (n**2)*(cold + cold_rec)

end subroutine cooling_low
