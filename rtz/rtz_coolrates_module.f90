! rtz_coolrates_module.f90
module rtz_coolrates_module
  use amr_parameters, only: dp
  implicit none

  private  ! everything is private by default
  public :: initialize_high_temperature_metal_cooling, initialize_fine_structure_tables
  public :: all_cooling

  ! Constants
  real(dp):: MIN_METAL_COOL_DENS = 1.d-10 ! minimum number density to calc cooling
  real(dp):: MIN_COOL_ION = 1.d-10 !//minimum number density to calculate cooling for a particular ion
  real(dp):: KB = 1.380649d-16 ! erg/K
  real(dp):: C_CGS = 2.9979246d10 ! cm/s
  real(dp):: H_PLANCK = 6.626196d-27 ! erg * s
  real(dp):: EV_2_ERG = 1.602d-12 

  ! Initialize the relevant arrays
  integer, parameter:: N_HIGH_T_COOLING_TEMP = 121
  real(dp):: high_t_cooling_temp(N_HIGH_T_COOLING_TEMP)
  real(dp):: high_t_cooling_rates(N_HIGH_T_COOLING_TEMP,27,27)
  logical:: high_t_cooling_rates_tflag(27,27)
  real(dp):: fs_cool_tab(27,160,8)  ! Array for fine structure cooling rates  

  real(dp):: G0_heating_rates(26) = (/ &
    0.d0, & ! 1 - Hydrogen
    0.d0, & ! 2 - Helium
    0.d0, & ! 3
    0.d0, & ! 4
    0.d0, & ! 5
    4.561312122d-22, & ! 6 - Carbon
    0.d0, & ! 7 - Nitrogen
    0.d0, & ! 8 - Oxygen
    0.d0, & ! 9
    0.d0, & ! 10 - Neon
    0.d0, & ! 11
    2.123365613d-22, & ! 12 - Magnesium
    0.d0, & ! 13
    1.31761296d-20, & ! 14 - Silicon
    0.d0, & ! 15
    2.03835276d-21, & ! 16 - Sulfur
    0.d0, & ! 17
    0.d0, & ! 18
    0.d0, & ! 19
    0.d0, & ! 20
    0.d0, & ! 21
    0.d0, & ! 22
    0.d0, & ! 23
    0.d0, & ! 24
    0.d0, & ! 25
    1.451738808d-21 & ! 26 - Iron
    /)

CONTAINS

FUNCTION collisional_ionization_cooling_HI(T) result(rate)
    ! From Cen 1992
    implicit none
    real(dp), intent(in):: T
    real(dp):: rate
    real(dp)::term_1, term_2, term_3

    term_1 = 1.27d-21 * sqrt(T)
    term_2 = 1.d0 /(1.d0 + sqrt(T/1.d5))
    term_3 = exp(-157809.1d0/T)
    rate = term_1 * term_2 * term_3

END FUNCTION collisional_ionization_cooling_HI

FUNCTION collisional_ionization_cooling_HeI(T) result(rate)
    ! From Cen 1992
    implicit none
    real(dp), intent(in):: T
    real(dp):: rate
    real(dp)::term_1, term_2, term_3

    term_1 = 9.38d-22 * sqrt(T)
    term_2 = 1.d0 /(1.d0 + sqrt(T/1.d5))
    term_3 = exp(-285335.4d0/T)
    rate = term_1 * term_2 * term_3

END FUNCTION collisional_ionization_cooling_HeI

FUNCTION collisional_ionization_cooling_HeII(T) result(rate)
    ! From Cen 1992
    implicit none
    real(dp), intent(in):: T
    real(dp):: rate
    real(dp)::term_1, term_2, term_3

    term_1 = 4.95d-22 * sqrt(T)
    term_2 = 1.d0 /(1.d0 + sqrt(T/1.d5))
    term_3 = exp(-631515.d0/T)
    rate = term_1 * term_2 * term_3

END FUNCTION collisional_ionization_cooling_HeII

FUNCTION recombination_cooling_case_B_HII(T) result(rate)
    ! From Hui & Gnedin 1997
    implicit none
    real(dp), intent(in):: T
    real(dp):: rate
    real(dp)::term_1, term_2, term_3, lam_HI

    lam_HI = 315614.d0 / T
    term_1 = 3.435d-30 * T
    term_2 = lam_HI**1.97d0
    term_3 = (1.d0 + ((lam_HI/2.25d0)**0.376d0))**(-3.72d0)
    rate = term_1 * term_2 * term_3

END FUNCTION recombination_cooling_case_B_HII

FUNCTION recombination_cooling_case_B_HeII(T) result(rate)
    ! From Hui & Gnedin 1997
    implicit none
    real(dp), intent(in):: T
    real(dp):: rate
    real(dp)::term_1, term_2, lam_HeI

    lam_HeI = 570670.d0 / T
    term_1 = KB * T * 1.26d-14
    term_2 = lam_HeI**0.75d0
    rate = term_1 * term_2

END FUNCTION recombination_cooling_case_B_HeII

FUNCTION recombination_cooling_case_B_HeIII(T) result(rate)
    ! From Hui & Gnedin 1997
    implicit none
    real(dp), intent(in):: T
    real(dp):: rate
    real(dp)::term_1, term_2, term_3, lam_HeII

    lam_HeII = 1263030.d0 / T
    term_1 = 8.d0 * 3.435d-30 * T
    term_2 = lam_HeII**1.97d0
    term_3 = (1.d0 + ((lam_HeII/2.25d0)**0.376d0))**(-3.72d0)
    rate = term_1 * term_2 * term_3

END FUNCTION recombination_cooling_case_B_HeIII

FUNCTION dielectronic_recombination_cooling_HeII(T) result(rate)
    ! From Black 1981
    implicit none
    real(dp), intent(in):: T
    real(dp):: rate
    real(dp)::term_1, term_2, term_3

    term_1 = 1.24d-13 * (T**(-1.5d0))
    term_2 = exp(-470000.d0/T)
    term_3 = 1.d0 + (0.3d0 * exp(-94000.d0/T))
    rate = term_1 * term_2 * term_3

END FUNCTION dielectronic_recombination_cooling_HeII

FUNCTION collisional_excitation_cooling_HI(T) result(rate)
    ! From Cen 1992
    implicit none
    real(dp), intent(in):: T
    real(dp):: rate
    real(dp)::term_1, term_2, term_3

    term_1 = 7.5d-19
    term_2 = 1.d0 /(1.d0 + sqrt(T/1.d5))
    term_3 = exp(-118348.d0/T)
    rate = term_1 * term_2 * term_3

END FUNCTION collisional_excitation_cooling_HI

FUNCTION collisional_excitation_cooling_HI_seon20(T) result(rate)
    ! From Cen 1992
    implicit none
    real(dp), intent(in):: T
    real(dp):: rate

    rate = 4.13d-19 * exp(-117744.0d0/T)

END FUNCTION collisional_excitation_cooling_HI_seon20

FUNCTION collisional_excitation_cooling_HeII(T) result(rate)
    ! From Cen 1992
    implicit none
    real(dp), intent(in):: T
    real(dp):: rate
    real(dp)::term_1, term_2, term_3

    term_1 = 5.54d-17 * (T**(-0.397d0))
    term_2 = 1.d0 /(1.d0 + sqrt(T/1.d5))
    term_3 = exp(-473638.d0/T)
    rate = term_1 * term_2 * term_3

END FUNCTION collisional_excitation_cooling_HeII

FUNCTION bremmstrahlung(T) result(rate)
    ! Osterbrock and Ferland 2006
    implicit none
    real(dp), intent(in):: T
    real(dp):: rate

    rate = 1.42d-27 * sqrt(T)

END FUNCTION bremmstrahlung

FUNCTION compton_cooling(T, a) result(rate)
    ! Haiman+ 1996
    implicit none
    real(dp), intent(in):: T, a
    real(dp):: rate
    real(dp)::term_1, term_2, term_3

    term_1 = 1.017d-37
    term_2 = (2.727d0/a)**4.d0
    term_3 = T - (2.727d0/a)
    rate = term_1 * term_2 * term_3

END FUNCTION compton_cooling

FUNCTION H2_cooling(nH, nH2, T) result(rate)
    ! From Moseley 2021
    implicit none
    real(dp), intent(in):: nH, nH2, T
    real(dp):: rate
    real(dp):: T3
    real(dp):: n1, n2, n3, n4
    real(dp):: x1, x2, x3, x4
    real(dp):: f1, f2, f3, f4

    T3 = 1d-3 * T
    n1 = 50.d0
    n2 = 450.d0
    n3 = 25.d0
    n4 = 900.d0

    x1 = nH + (5.0d0 * nH2)
    x2 = nH + (4.5d0 * nH2)
    x3 = nH + (0.75d0 * nH2)
    x4 = nH + (0.05d0 * nH2)

    f1 = 1.1d-25 * sqrt(T3) * exp(-0.51d0 / T3)
    f1 = f1 * (((0.7d0 * x1) / (1.d0 + (x1/n1))) + ((0.30d0 * x1)/(1.d0 + (x1/(10.d0*n1)))))

    f2 = 2.0d-25 * T3 * exp(-1.d0 / T3)
    f2 = f2 * (((0.35d0 * x2) / (1.d0 + (x2/n2))) + ((0.65d0 * x2)/(1.d0 + (x2/(10.d0*n2)))))

    f3 = 2.4d-24 * (T3**1.50d0) * exp(-2.0d0 / T3)
    f3 = f3 * (x3 / (1.d0 + (x3/n3)))

    f4 = 1.7d-23 * (T3**1.50d0) * exp(-4.0d0 / T3)
    f4 = f4 * (((0.45d0 * x4) / (1.d0 + (x4/n4))) + ((0.55d0 * x4)/(1.d0 + (x4/(10.d0*n4)))))

    rate = nH2 * (f1 + f2 + f3 + f4)

END FUNCTION H2_cooling

FUNCTION H2_cooling_G15(T) result(rate)
    ! Cooling due to molecular hydrogen from https://ui.adsabs.harvard.edu/abs/2015MNRAS.453.2901G/abstract
    ! see appendix 4.2.1 eqn 30
    ! returns erg/s
    implicit none
    real(dp), intent(in)::T
    real(dp)::rate
    real(dp)::logT3

    logT3 = LOG10(T/1.d3)
    rate = 0.d0
    rate = -20.584225d0
    rate = rate + (5.0194035d0 * logT3)
    rate = rate - (1.5738805d0 * logT3**2.d0)
    rate = rate - (4.7155769d0 * logT3**3.d0)
    rate = rate + (2.4714161d0 * logT3**4.d0)
    rate = rate + (5.4710750d0 * logT3**5.d0)
    rate = rate - (3.9467356d0 * logT3**6.d0)
    rate = rate - (2.2148338d0 * logT3**7.d0)
    rate = rate + (1.8161874d0 * logT3**8.d0)
    rate = (10.d0**(rate))

END FUNCTION H2_cooling_G15

FUNCTION H2_cooling_GA08(T, n_e, n_HI, n_HII, n_HeI, n_H2) result(rate)
    ! Cooling due to molecular hydrogen
    implicit none
    real(dp), intent(in)::T, n_e, n_HI, n_HII, n_HeI, n_H2
    real(dp)::t3, lt, lt3, ltt, loc_T
    real(dp)::gphdl, HDLR, HDLV, gaHI, gaH2, gaHp, gaHe, gael, galdl
    real(dp)::rate

    rate = 0d0
    loc_T = min(T,9999.9d0)
    if (loc_T .le. 10.d0) then
        return
    end if
    
    lt = log10(loc_T)
    t3 = loc_T/1000.d0
    lt3 = log10(t3)

    ! From glover 2011
    gphdl = H2_cooling_G15(loc_T)

    ! Initialise
    gaHI = 0.d0
    gaH2 = 0.d0
    gaHp = 0.d0
    gaHe = 0.d0
    gael = 0.d0

    ! Glover & Abel (2008) ; low-density limit

    ! Excitation by HI
    if (loc_T .lt. 100.d0) then
        gaHI = 10.d0**(-16.818342d0        &
        &         + 37.383713d0*lt3      &
        &         + 58.145166d0*lt3**2.d0   &
        &         + 48.656103d0*lt3**3.d0   &
        &         + 20.159831d0*lt3**4.d0   &
        &         + 3.847961d0*lt3**5.d0)
    else if (loc_T .lt. 1000.d0) then
        gaHI = 10.d0**(-24.311209d0        &
        &         + 3.5692468d0*lt3      &
        &         - 11.33286d0*lt3**2.d0    &
        &         - 27.850082d0*lt3**3.d0   &
        &         - 21.328264d0*lt3**4.d0   &
        &         - 4.2519023d0*lt3**5.d0)
    ! harley checked that extrapolation to 1e4K is ok
    else if (loc_T .lt. 10000.d0) then
        gaHI = 10.d0**(-24.311209d0        &
        &         + 4.6450521d0*lt3       &
        &         - 3.7209846d0*lt3**2.d0    &
        &         + 5.9369081d0*lt3**3.d0    &
        &         - 5.5108047d0*lt3**4.d0    &
        &         + 1.5538288d0*lt3**5.d0)
    end if

    ! Excitation by H2
    ! harley checked that extrapolation to 1e4K is ok
    if (loc_T .lt. 10000.d0) then
        gaH2 = 10.d0**(-23.962112d0         &
        &        + 2.0943374d0*lt3       &
        &        - 0.77151436d0*lt3**2.d0   &
        &        + 0.43693353d0*lt3**3.d0   &
        &        - 0.14913216d0*lt3**4.d0   &
        &        - 0.033638326d0*lt3**5.d0)
    end if

    ! Excitation by He
    ! harley checked that extrapolation to 1e4K is ok
    if (loc_T .lt. 10000.d0) then
        gaHe = 10.d0**(-23.689237d0        &
        &        + 2.1892372d0*lt3      &
        &        - 0.81520438d0*lt3**2.d0   &
        &        + 0.29036281d0*lt3**3.d0   &
        &        - 0.16596184d0*lt3**4.d0   &
        &        + 0.19191375d0*lt3**5.d0)
    end if

    ! Excitation by H+
    ! Update from glover 2015
    if (loc_T .lt. 10000.d0) then
        gaHp = 10.d0**(-22.089523d0      &
        &        + 1.5714711d0*lt3       &
        &        + 0.015391166d0*lt3**2.d0  &
        &        - 0.23619985d0*lt3**3.d0   &
        &        - 0.51002221d0*lt3**4.d0   &
        &        + 0.32168730d0*lt3**5.d0)
    end if

    ! Excitation by electrons
    ! Update from glover 2015
    if (loc_T .lt. 500.d0) then
        gael = 10.d0**(-21.928796d0         &
        &           + 16.815730d0*lt3       &
        &           + 96.743155d0*lt3**2.d0    &
        &           + 343.19180d0*lt3**3.d0    &
        &           + 734.71651d0*lt3**4.d0    &
        &           + 983.67576d0*lt3**5.d0    &
        &           + 801.81247d0*lt3**6.d0    &
        &           + 364.14446d0*lt3**7.d0    &
        &           + 70.609154d0*lt3**8.d0)    
    else if (loc_T .lt. 10000.d0) then
        gael = 10.d0**(-22.921189d0         &
        &           + 1.6802758d0*lt3       &
        &           + 0.93310622d0*lt3**2.d0   &
        &           + 4.0406627d0*lt3**3.d0    &
        &           - 4.7274036d0*lt3**4.d0    &
        &           - 8.8077017d0*lt3**5.d0    &
        &           + 8.9167183d0*lt3**6.d0    &
        &           + 6.4380698d0*lt3**7.d0    &
        &           - 6.3701156d0*lt3**8.d0)    
    end if

    galdl = gaHI*n_HI + gaH2*n_H2 + gaHe*n_HeI + gaHp*n_HII + gael*n_e ! erg/s

    rate = n_H2*gphdl/(1.d0 + gphdl/galdl) ! erg/cm^3/s

END FUNCTION H2_cooling_GA08

FUNCTION cooling_H2GP(nH,nH2,Tgas) result(rate)
    !cooling from Galli&Palla98
    !taken fron krome
    real(dp), intent(in)::nH,nH2,Tgas
    real(dp)::tm, logT, T3
    real(dp)::rate
    real(dp)::LDL,HDLR,HDLV,HDL

    rate = 0.d0

    tm = max(Tgas, 13.0d0)    ! no cooling below 13 Kelvin
    tm = min(Tgas, 1.d5)      ! fixes numerics
    logT = log10(tm)
    T3 = tm * 1.d-3

    !low density limit in erg/s
    LDL = 1.d1**(-103.d0+97.59d0*logT-48.05d0*logT**2 + 10.8d0*logT**3-0.9032d0*logT**4)*nH

    !this will avoid a division by zero and useless calculations
    if(LDL.eq.0d0) then
        rate = 0d0
        return
    end if

    !high density limit
    HDLR = ((9.5e-22*t3**3.76)/(1.+0.12*t3**2.1)*exp(-(0.13/t3)**3)+&
        3.e-24*exp(-0.51/t3)) !erg/s
    HDLV = (6.7e-19*exp(-5.86/t3) + 1.6e-18*exp(-11.7/t3)) !erg/s
    HDL  = HDLR + HDLV !erg/s

    !to avoid division by zero
    if (HDL.eq.0d0) then
        rate = 0d0
    else
        rate = nH2/(1d0/HDL+1d0/LDL) !erg/cm3/s
    endif

END FUNCTION cooling_H2GP

SUBROUTINE initialize_high_temperature_metal_cooling()
    ! Cloudy tables of metal line cooling
    ! which are valid at high temperature
    use amr_commons, only: myid
    implicit none
    integer :: unit_num, ios, i, j

    if(myid.eq.1) write(*,*) 'Initializing High Temperature Metal Cooling Rates'

    ! Load temperatures
    open(newunit=unit_num, file='./data/high_T_cooling/temperatures.dat', status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(*,*) 'Error: Could not open high T temperatures file'
        return
    end if

    do i = 1, N_HIGH_T_COOLING_TEMP
        read(unit_num, *, iostat=ios) high_t_cooling_temp(i)
        if (ios /= 0) exit
    end do
    close(unit_num)

    ! Initialize high-temeprature cooling table to 0
    ! We use -50 because it's log
    high_t_cooling_rates = -50.d0

    ! CARBON
    open(newunit=unit_num, file='./data/high_T_cooling/CARBON/all_cool.dat', status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(*,*) 'Error: Could not open high temperature cooling carbon rates'
        return
    end if

    do i = 1, 7
        read(unit_num, *, iostat=ios) high_t_cooling_rates(:,6,i)
        if (ios /= 0) exit
    end do
    close(unit_num)

    ! NITROGEN
    open(newunit=unit_num, file='./data/high_T_cooling/NITROGEN/all_cool.dat', status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(*,*) 'Error: Could not open high temperature cooling nitrogen rates'
        return
    end if

    do i = 1, 8
        read(unit_num, *, iostat=ios) high_t_cooling_rates(:,7,i)
        if (ios /= 0) exit
    end do
    close(unit_num)

    ! OXYGEN
    open(newunit=unit_num, file='./data/high_T_cooling/OXYGEN/all_cool.dat', status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(*,*) 'Error: Could not open high temperature cooling oxygen rates'
        return
    end if

    do i = 1, 9
        read(unit_num, *, iostat=ios) high_t_cooling_rates(:,8,i)
        if (ios /= 0) exit
    end do
    close(unit_num)

    ! NEON
    open(newunit=unit_num, file='./data/high_T_cooling/NEON/all_cool.dat', status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(*,*) 'Error: Could not open high temperature cooling neon rates'
        return
    end if

    do i = 1, 11
        read(unit_num, *, iostat=ios) high_t_cooling_rates(:,10,i)
        if (ios /= 0) exit
    end do
    close(unit_num)

    ! MAGNESIUM
    open(newunit=unit_num, file='./data/high_T_cooling/MAGNESIUM/all_cool.dat', status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(*,*) 'Error: Could not open high temperature cooling magnesium rates'
        return
    end if

    do i = 1, 13
        read(unit_num, *, iostat=ios) high_t_cooling_rates(:,12,i)
        if (ios /= 0) exit
    end do
    close(unit_num)

    ! SILICON
    open(newunit=unit_num, file='./data/high_T_cooling/SILICON/all_cool.dat', status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(*,*) 'Error: Could not open high temperature cooling silicon rates'
        return
    end if

    do i = 1, 15
        read(unit_num, *, iostat=ios) high_t_cooling_rates(:,14,i)
        if (ios /= 0) exit
    end do
    close(unit_num)

    ! SULFUR
    open(newunit=unit_num, file='./data/high_T_cooling/SULPHUR/all_cool.dat', status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(*,*) 'Error: Could not open high temperature cooling sulfur rates'
        return
    end if

    do i = 1, 17
        read(unit_num, *, iostat=ios) high_t_cooling_rates(:,16,i)
        if (ios /= 0) exit
    end do
    close(unit_num)

    ! IRON
    open(newunit=unit_num, file='./data/high_T_cooling/IRON/all_cool.dat', status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(*,*) 'Error: Could not open high temperature cooling iron rates'
        return
    end if

    do i = 1, 27
        read(unit_num, *, iostat=ios) high_t_cooling_rates(:,26,i)
        if (ios /= 0) exit
    end do
    close(unit_num)

    ! Set all of the fine structure lines to true
    ! prob a better way of doing this but fine for now
    high_t_cooling_rates_tflag = .false.

    ! CI
    high_t_cooling_rates_tflag(6,1) = .true.
    ! CII
    high_t_cooling_rates_tflag(6,2) = .true.
    ! NII
    high_t_cooling_rates_tflag(7,2) = .true.
    ! OI
    high_t_cooling_rates_tflag(8,1) = .true.
    ! OIII
    high_t_cooling_rates_tflag(8,3) = .true.
    ! NeII
    high_t_cooling_rates_tflag(10,2) = .true.
    ! SiI
    high_t_cooling_rates_tflag(14,1) = .true.
    ! SiII
    high_t_cooling_rates_tflag(14,2) = .true.
    ! SI
    high_t_cooling_rates_tflag(16,1) = .true.
    ! FeI
    high_t_cooling_rates_tflag(26,1) = .true.
    ! FeII
    high_t_cooling_rates_tflag(26,2) = .true.

END SUBROUTINE initialize_high_temperature_metal_cooling

FUNCTION get_high_t_cooling_rates(T, ne, element_number_densities, element_ion_fractions, temp_smooth) result(rate)
    ! Computes the high-temperature cooling rates
    ! derived from Harley's custom cloudy models
    use rtz_module, only: elements
    implicit none
    real(dp), intent(in):: T, ne, temp_smooth
    real(dp), intent(in):: element_number_densities(27)
    real(dp), intent(in):: element_ion_fractions(27,27)
    real(dp):: rate
    real(dp):: loc_T, log_T, t_min, t_max, dt, t_scale_fac
    real(dp):: frac_high, frac_low, loc_cooling_rate, loc_temp_smooth
    integer:: idx_low, i, j

    loc_T = max(T,1000.d0)
    log_T = log10(loc_T)
    t_min = 3.d0
    t_max = 9.d0
    dt = 0.05d0
    rate = 0.d0
    t_scale_fac = exp(-1.d0 * ((2000.d0/loc_T)**5.d0))

    ! bounds for temperature --> no cooling
    if (log_T < t_min) then
        return
    end if

    !// Still cool above 10^9 K but set a bound
    if (log_T > t_max) then 
        log_T = t_max
    end if

    !// Prepare for 1D interpolation
    idx_low = floor((log_T - t_min)/dt) + 1

    frac_high = (log_T - high_t_cooling_temp(idx_low)) / (high_t_cooling_temp(idx_low+1) - high_t_cooling_temp(idx_low))
    frac_low = 1.0 - frac_high

    ! Loop over all ions -- not including H and He
    do i = 3, 27
       ! Skip unused elements
       if (elements(i)%atomic_number .lt. 1) then
          cycle
       end if

       ! Skip elements with too low number density
       if (element_number_densities(i) .lt. MIN_COOL_ION) then
          cycle
       end if

       do j = 1, elements(i)%n_ions
          ! Interpolate the cooling rate for the ion
          loc_cooling_rate = frac_low * high_t_cooling_rates(idx_low,i,j)
          loc_cooling_rate = loc_cooling_rate + (frac_high * high_t_cooling_rates(idx_low+1,i,j)) 
          loc_cooling_rate = 10.d0**loc_cooling_rate

          ! Smooth with temeprature if necessary
          loc_temp_smooth = 1.d0
          if (high_t_cooling_rates_tflag(i,j)) then
              loc_temp_smooth = temp_smooth
          end if


         ! Multiply cooling rate by ion and electron number densities
         rate = rate + (element_number_densities(i) * element_ion_fractions(i,j) * ne * loc_cooling_rate * loc_temp_smooth)

       end do
    end do

    rate = rate * t_scale_fac

END FUNCTION get_high_t_cooling_rates

SUBROUTINE initialize_fine_structure_tables()
    ! Initialization for fine structure cooling tables
    use amr_commons, only: myid
    implicit none
    integer :: unit_num, ios, i, j
    integer, parameter::N_LINES = 27
    character(len=19), parameter, dimension(N_LINES) :: file_names = (/ &
        'CII_158um_rates.dat', 'CI_609um_rates.dat ', 'CI_230um_rates.dat ', 'CI_370um_rates.dat ', &
        'NII_205um_rates.dat', 'NII_76um_rates.dat ', 'NII_122um_rates.dat', 'OI_63um_rates.dat  ',  &
        'OI_44um_rates.dat  ', 'OI_145um_rates.dat ', 'NeII_13um_rates.dat', 'SiII_35um_rates.dat', &
        'OIII_88um_rates.dat', 'OIII_33um_rates.dat', 'OIII_52um_rates.dat', 'FeI_14um_rates.dat ', &
        'FeI_24um_rates.dat ', 'FeI_35um_rates.dat ', 'FeII_26um_rates.dat', 'FeII_15um_rates.dat', &
        'FeII_35um_rates.dat', 'SiI_130um_rates.dat', 'SiI_45um_rates.dat ', 'SiI_68um_rates.dat ', &
        'SI_17um_rates.dat  ', 'SI_25um_rates.dat  ', 'SI_56um_rates.dat  ' /)

    if(myid.eq.1) write(*,*) 'Initializing Fine Structure Cooling Tables'

    ! Check that all files exist and load them in
    do i=1,N_LINES
        ! Open the fine structure file
        open(newunit=unit_num, file='./data/fine_structure_data/'//trim(file_names(i)), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*,*) 'Error: Could not open '//trim(file_names(i))//' high temperature cooling iron rates'
            return
        end if

        ! Read the file into the fine structure array
        do j = 1, 160
            read(unit_num, *, iostat=ios) fs_cool_tab(i,j,:)
            if (ios /= 0) exit
        end do
        close(unit_num)

    end do
    
END SUBROUTINE initialize_fine_structure_tables

FUNCTION three_level(g_0, g_1, g_2, lam_10, lam_20, lam_21, &
                     A_10, A_20, A_21, z, T, n_ion, ne, nH, &
                     nHp, nHe, nHep, nHepp, nH2, &
                     idx_0, idx_1, idx_2) result(rate)
    ! Harley's three-level ion solver
    implicit none
    real(dp), intent(in):: g_0, g_1, g_2
    real(dp), intent(in):: lam_10, lam_20, lam_21
    real(dp), intent(in):: A_10, A_20, A_21
    real(dp), intent(in):: z, T, n_ion
    real(dp), intent(in):: ne, nH, nHp, nHe, nHep, nHepp, nH2
    integer, intent(in):: idx_0, idx_1, idx_2
    real(dp):: rate

    real(dp):: T_cmb, logT, tmin, tmax, delta_temp
    real(dp):: nu_10, nu_20, nu_21
    real(dp):: E_10, E_20, E_21
    real(dp):: B_01, B_02, B_12, B_10, B_20, B_21
    real(dp):: B_nu_10, B_nu_20, B_nu_21
    real(dp):: frac_high, frac_low
    real(dp):: q10_oH2, q20_oH2, q21_oH2
    real(dp):: q10_pH2, q20_pH2, q21_pH2
    real(dp):: q10_H, q20_H, q21_H
    real(dp):: q10_Hp, q20_Hp, q21_Hp
    real(dp):: q10_e, q20_e, q21_e
    real(dp):: q10_He, q20_He, q21_He
    real(dp):: q10_Hep, q20_Hep, q21_Hep
    real(dp):: q10_Hepp, q20_Hepp, q21_Hepp
    real(dp):: C_10, C_20, C_21, C_01, C_02, C_12
    real(dp):: n2_n1_top, n2_n1_bot, n2_n1
    real(dp):: n1_n0_top, n1_n0_bot, n1_n0
    real(dp):: n_0, n_1, n_2
    real(dp):: cool_0, cool_1, cool_2
    real(dp):: heat_0, heat_1, heat_2

    integer:: itemp_low, itemp_high

    ! Initialize the result
    rate = 1d-50

    T_cmb = 2.725d0 * (1.d0 + z) ! CMB temperature

    ! Get log T and enforce bounds
    logT = log10(T)

    !// Enforce bounds on temperature
    tmin = 1.d0
    tmax = 4.975d0
    delta_temp = 0.025d0
    logT = max(logT,1.d0)
    logT = min(logT,tmax)

    nu_10 = C_CGS / (lam_10 * 1d-4) ! Hz
    nu_20 = C_CGS / (lam_20 * 1d-4) ! Hz
    nu_21 = C_CGS / (lam_21 * 1d-4) ! Hz

    E_10 = (H_PLANCK * C_CGS / (lam_10 * 1d-4)) / KB ! E/K (K)
    E_20 = (H_PLANCK * C_CGS / (lam_20 * 1d-4)) / KB ! E/K (K)
    E_21 = (H_PLANCK * C_CGS / (lam_21 * 1d-4)) / KB ! E/K (K)

    B_01 = A_10 * ((lam_10 * 1.d-4)**3.d0) * (g_1 / g_0) / (2.d0 * H_PLANCK * C_CGS)
    B_02 = A_20 * ((lam_20 * 1.d-4)**3.d0) * (g_2 / g_0) / (2.d0 * H_PLANCK * C_CGS)
    B_12 = A_21 * ((lam_21 * 1.d-4)**3.d0) * (g_2 / g_1) / (2.d0 * H_PLANCK * C_CGS)

    B_10 = (g_0 / g_1) * B_01
    B_20 = (g_0 / g_2) * B_02
    B_21 = (g_1 / g_2) * B_12

    ! CMB black body spectrum
    B_nu_10 = (2.d0 * H_PLANCK * (nu_10**3.d0) / (C_CGS*C_CGS)) / (exp(H_PLANCK * nu_10 / (KB * T_cmb)) - 1.d0)
    B_nu_20 = (2.d0 * H_PLANCK * (nu_20**3.d0) / (C_CGS*C_CGS)) / (exp(H_PLANCK * nu_20 / (KB * T_cmb)) - 1.d0)
    B_nu_21 = (2.d0 * H_PLANCK * (nu_21**3.d0) / (C_CGS*C_CGS)) / (exp(H_PLANCK * nu_21 / (KB * T_cmb)) - 1.d0)

    ! Find temperature index
    itemp_low = floor((logT - tmin)/delta_temp) + 1
    itemp_high = itemp_low + 1
    frac_high = (logT - (tmin + (real(itemp_low-1,dp) * delta_temp))) / delta_temp
    frac_low = 1.d0 - frac_high

    ! All collison strengths are in cm^3 s^-1
    ! For collisions with o-H2
    q10_oH2 = (fs_cool_tab(idx_0,itemp_low,7) * frac_low) + (fs_cool_tab(idx_0,itemp_high,7) * frac_high)
    q20_oH2 = (fs_cool_tab(idx_1,itemp_low,7) * frac_low) + (fs_cool_tab(idx_1,itemp_high,7) * frac_high)
    q21_oH2 = (fs_cool_tab(idx_2,itemp_low,7) * frac_low) + (fs_cool_tab(idx_2,itemp_high,7) * frac_high)

    ! For collisions with p-H2
    q10_pH2 = (fs_cool_tab(idx_0,itemp_low,8) * frac_low) + (fs_cool_tab(idx_0,itemp_high,8) * frac_high)
    q20_pH2 = (fs_cool_tab(idx_1,itemp_low,8) * frac_low) + (fs_cool_tab(idx_1,itemp_high,8) * frac_high)
    q21_pH2 = (fs_cool_tab(idx_2,itemp_low,8) * frac_low) + (fs_cool_tab(idx_2,itemp_high,8) * frac_high)

    ! For collisions with H
    q10_H = (fs_cool_tab(idx_0,itemp_low,5) * frac_low) + (fs_cool_tab(idx_0,itemp_high,5) * frac_high)
    q20_H = (fs_cool_tab(idx_1,itemp_low,5) * frac_low) + (fs_cool_tab(idx_1,itemp_high,5) * frac_high)
    q21_H = (fs_cool_tab(idx_2,itemp_low,5) * frac_low) + (fs_cool_tab(idx_2,itemp_high,5) * frac_high)

    ! For collisions with H+
    q10_Hp = (fs_cool_tab(idx_0,itemp_low,2) * frac_low) + (fs_cool_tab(idx_0,itemp_high,2) * frac_high)
    q20_Hp = (fs_cool_tab(idx_1,itemp_low,2) * frac_low) + (fs_cool_tab(idx_1,itemp_high,2) * frac_high)
    q21_Hp = (fs_cool_tab(idx_2,itemp_low,2) * frac_low) + (fs_cool_tab(idx_2,itemp_high,2) * frac_high)

    ! For collisions with e
    q10_e = (fs_cool_tab(idx_0,itemp_low,1) * frac_low) + (fs_cool_tab(idx_0,itemp_high,1) * frac_high)
    q20_e = (fs_cool_tab(idx_1,itemp_low,1) * frac_low) + (fs_cool_tab(idx_1,itemp_high,1) * frac_high)
    q21_e = (fs_cool_tab(idx_2,itemp_low,1) * frac_low) + (fs_cool_tab(idx_2,itemp_high,1) * frac_high)

    ! For collisions with He
    q10_He = (fs_cool_tab(idx_0,itemp_low,6) * frac_low) + (fs_cool_tab(idx_0,itemp_high,6) * frac_high)
    q20_He = (fs_cool_tab(idx_1,itemp_low,6) * frac_low) + (fs_cool_tab(idx_1,itemp_high,6) * frac_high)
    q21_He = (fs_cool_tab(idx_2,itemp_low,6) * frac_low) + (fs_cool_tab(idx_2,itemp_high,6) * frac_high)

    ! For collisions with He+
    q10_Hep = (fs_cool_tab(idx_0,itemp_low,3) * frac_low) + (fs_cool_tab(idx_0,itemp_high,3) * frac_high)
    q20_Hep = (fs_cool_tab(idx_1,itemp_low,3) * frac_low) + (fs_cool_tab(idx_1,itemp_high,3) * frac_high)
    q21_Hep = (fs_cool_tab(idx_2,itemp_low,3) * frac_low) + (fs_cool_tab(idx_2,itemp_high,3) * frac_high)

    ! For collisions with He++
    q10_Hepp = (fs_cool_tab(idx_0,itemp_low,4) * frac_low) + (fs_cool_tab(idx_0,itemp_high,4) * frac_high)
    q20_Hepp = (fs_cool_tab(idx_1,itemp_low,4) * frac_low) + (fs_cool_tab(idx_1,itemp_high,4) * frac_high)
    q21_Hepp = (fs_cool_tab(idx_2,itemp_low,4) * frac_low) + (fs_cool_tab(idx_2,itemp_high,4) * frac_high)

    ! Net collision strengths note the 75-25 ortho-para H2
    C_10 = (q10_e * ne) + (q10_H * nH) + (q10_Hp * nHp) + (q10_oH2 * 0.75d0 * nH2) + (q10_pH2 * 0.25d0 * nH2) + (q10_He * nHe) + (q10_Hep * nHep) + (q10_Hepp * nHepp)
    C_20 = (q20_e * ne) + (q20_H * nH) + (q20_Hp * nHp) + (q20_oH2 * 0.75d0 * nH2) + (q20_pH2 * 0.25d0 * nH2) + (q20_He * nHe) + (q20_Hep * nHep) + (q20_Hepp * nHepp)
    C_21 = (q21_e * ne) + (q21_H * nH) + (q21_Hp * nHp) + (q21_oH2 * 0.75d0 * nH2) + (q21_pH2 * 0.25d0 * nH2) + (q21_He * nHe) + (q21_Hep * nHep) + (q21_Hepp * nHepp)

    C_01 = C_10 * (g_1/g_0) * exp(-1.d0 * E_10 / T)
    C_02 = C_20 * (g_2/g_0) * exp(-1.d0 * E_20 / T)
    C_12 = C_21 * (g_2/g_1) * exp(-1.d0 * E_21 / T)

    C_01 = C_01 + (B_01*B_nu_10)
    C_02 = C_02 + (B_02*B_nu_20)
    C_12 = C_12 + (B_12*B_nu_21)

    C_10 = C_10 + (B_10*B_nu_10)
    C_20 = C_20 + (B_20*B_nu_20)
    C_21 = C_21 + (B_21*B_nu_21)

    ! Analytic solution to the 3 level system (see paul goldsmith papers)
    ! Done this way to avoid numerical errors
    n2_n1_top = (C_12 * (C_01 + C_02)) + (C_02 * (A_10 + C_10))
    n2_n1_bot = ((A_21 + C_21 + C_20) * (C_01 + C_02)) - (C_20 * C_02)
    n2_n1 = n2_n1_top / n2_n1_bot

    n1_n0_top = ((A_21 + C_21 + C_20) * (C_01 + C_02)) - (C_20*C_02)
    n1_n0_bot = ((A_21 + C_21 + C_20) * (A_10 + C_10)) + (C_20*C_12)
    n1_n0 = n1_n0_top / n1_n0_bot

    n_0 = 1.d0 / (1.d0 + n1_n0 + (n2_n1 * n1_n0))
    n_1 = 1.d0 / (1.d0 + n2_n1 + (1.d0/n1_n0))
    n_2 = 1.d0 / (1.d0 + (1.d0/n2_n1) + ((1.d0/n2_n1)*(1.d0/n1_n0)))

    ! Cooling and heating rates
    cool_0 = (A_10 + (B_10 * B_nu_10)) * E_10 * KB * n_1 * n_ion
    cool_1 = (A_20 + (B_20 * B_nu_20)) * E_20 * KB * n_2 * n_ion
    cool_2 = (A_21 + (B_21 * B_nu_21)) * E_21 * KB * n_2 * n_ion

    heat_0 = B_01 * B_nu_10 * E_10 * KB * n_0 * n_ion
    heat_1 = B_02 * B_nu_20 * E_20 * KB * n_0 * n_ion
    heat_2 = B_12 * B_nu_21 * E_21 * KB * n_1 * n_ion

    ! Total cooling rate
    rate = (cool_0 + cool_1 + cool_2) - (heat_0 + heat_1 + heat_2)

END FUNCTION three_level

FUNCTION two_level(g_0, g_1, lam_10, A_10, z, T, &
                n_ion, ne, nH,  nHp,  nHe,  nHep, & 
                nHepp,  nH2, idx_0) result(rate)
    !Harley's two-level ion solver
    implicit none
    real(dp), intent(in):: g_0, g_1
    real(dp), intent(in):: lam_10
    real(dp), intent(in):: A_10
    real(dp), intent(in):: z, T, n_ion
    real(dp), intent(in):: ne, nH, nHp, nHe, nHep, nHepp, nH2
    integer, intent(in):: idx_0
    real(dp):: rate

    real(dp):: T_cmb, logT, tmin, tmax, delta_temp
    real(dp):: nu_10
    real(dp):: E_10
    real(dp):: B_01, B_10
    real(dp):: B_nu_10
    real(dp):: frac_high, frac_low
    real(dp):: q10_oH2
    real(dp):: q10_pH2
    real(dp):: q10_H
    real(dp):: q10_Hp
    real(dp):: q10_e
    real(dp):: q10_He
    real(dp):: q10_Hep
    real(dp):: q10_Hepp
    real(dp):: C_10, C_01
    real(dp):: nu_over_nl
    real(dp):: n_0, n_1
    real(dp):: cool_0
    real(dp):: heat_0
    real(dp):: t1, t2

    integer:: itemp_low, itemp_high

    ! Initialize the result
    rate = 1d-50

    T_cmb = 2.725d0 * (1.d0 + z) ! CMB temperature

    ! Get log T and enforce bounds
    logT = log10(T)

    !// Enforce bounds on temperature
    tmin = 1.d0
    tmax = 4.975d0
    delta_temp = 0.025d0
    logT = max(logT,1.d0)
    logT = min(logT,tmax)

    nu_10 = C_CGS / (lam_10 * 1d-4) ! Hz

    E_10 = (H_PLANCK * C_CGS / (lam_10 * 1d-4)) / KB ! E/K (K)

    B_01 = A_10 * ((lam_10 * 1.d-4)**3.d0) * (g_1 / g_0) / (2.d0 * H_PLANCK * C_CGS)

    B_10 = (g_0 / g_1) * B_01

    ! CMB black body spectrum
    B_nu_10 = (2.d0 * H_PLANCK * (nu_10**3.d0) / (C_CGS*C_CGS)) / (exp(H_PLANCK * nu_10 / (KB * T_cmb)) - 1.d0)

    ! Find temperature index.
    itemp_low = floor((logT - tmin)/delta_temp) + 1
    itemp_high = itemp_low + 1
    frac_high = (logT - (tmin + (real(itemp_low-1,dp) * delta_temp))) / delta_temp
    frac_low = 1.d0 - frac_high

    ! All collison strengths are in cm^3 s^-1
    ! For collisions with o-H2
    q10_oH2 = (fs_cool_tab(idx_0,itemp_low,7) * frac_low) + (fs_cool_tab(idx_0,itemp_high,7) * frac_high)

    ! For collisions with p-H2
    q10_pH2 = (fs_cool_tab(idx_0,itemp_low,8) * frac_low) + (fs_cool_tab(idx_0,itemp_high,8) * frac_high)

    ! For collisions with H
    q10_H = (fs_cool_tab(idx_0,itemp_low,5) * frac_low) + (fs_cool_tab(idx_0,itemp_high,5) * frac_high)

    ! For collisions with H+
    q10_Hp = (fs_cool_tab(idx_0,itemp_low,2) * frac_low) + (fs_cool_tab(idx_0,itemp_high,2) * frac_high)

    ! For collisions with e
    q10_e = (fs_cool_tab(idx_0,itemp_low,1) * frac_low) + (fs_cool_tab(idx_0,itemp_high,1) * frac_high)

    ! For collisions with He
    q10_He = (fs_cool_tab(idx_0,itemp_low,6) * frac_low) + (fs_cool_tab(idx_0,itemp_high,6) * frac_high)

    ! For collisions with He+
    q10_Hep = (fs_cool_tab(idx_0,itemp_low,3) * frac_low) + (fs_cool_tab(idx_0,itemp_high,3) * frac_high)

    ! For collisions with He++
    q10_Hepp = (fs_cool_tab(idx_0,itemp_low,4) * frac_low) + (fs_cool_tab(idx_0,itemp_high,4) * frac_high)

    ! Net collision strengths note the 75-25 ortho-para H2
    C_10 = (q10_e * ne) + (q10_H * nH) + (q10_Hp * nHp) + (q10_oH2 * 0.75d0 * nH2) + (q10_pH2 * 0.25d0 * nH2) + (q10_He * nHe) + (q10_Hep * nHep) + (q10_Hepp * nHepp)
    
    C_01 = C_10 * (g_1/g_0) * exp(-1.d0 * E_10 / T)

    ! Analytic solution to the 2 level system (modified version of paul goldsmith papers)
    ! Done this way to avoid numerical errors
    t1 = ( (B_01*B_nu_10) + C_01 )
    t2 = ( A_10 + (B_10*B_nu_10) + C_10 )

    nu_over_nl = t1 / t2
    n_0 = 1.d0 / (nu_over_nl + 1.d0)
    n_1 = 1.d0 / ((1.d0/nu_over_nl) + 1.d0)

    ! Cooling and heating rates
    cool_0 = (A_10 + (B_10 * B_nu_10)) * E_10 * KB * n_1 * n_ion

    heat_0 = B_01 * B_nu_10 * E_10 * KB * n_0 * n_ion

    ! Total cooling rate
    rate = cool_0 - heat_0

END FUNCTION two_level

! OI fine structure calculation function 
FUNCTION OI_fine_structure(T, n_ion, nH, nHp, &
                           ne, nH2, nHe, nHep, &
                           nHepp, z) result(rate)
    implicit none
    real(dp), intent(in):: T, z, n_ion, nH, nHp
    real(dp), intent(in):: ne, nH2, nHe, nHep, nHepp
    real(dp):: rate

    real(dp):: g_0, g_1, g_2 ! Degeneracy variables
    real(dp):: A_10, A_20, A_21 ! Transition rates (s^-1)
    real(dp):: lam_10, lam_20, lam_21 ! Wavelengths (microns)

    g_0 = 5.d0
    g_1 = 3.d0
    g_2 = 1.d0
    
    lam_10 = 63.1679d0
    lam_20 = 44.0453d0
    lam_21 = 145.495d0
    
    A_10 = 8.910d-5
    A_20 = 1.340d-10
    A_21 = 1.750d-5
    
    ! Call three_level function
    rate = three_level(g_0, g_1, g_2, &
                        lam_10, lam_20, lam_21, &
                        A_10, A_20, A_21, &
                        z, T, &
                        n_ion, ne, nH, nHp, nHe, nHep, nHepp, nH2, &
                        8, 9, 10)
END FUNCTION OI_fine_structure

! OIII fine structure calculation function 
FUNCTION OIII_fine_structure(T, n_ion, nH, nHp, &
                           ne, nH2, nHe, nHep, &
                           nHepp, z) result(rate)
    implicit none
    real(dp), intent(in):: T, z, n_ion, nH, nHp
    real(dp), intent(in):: ne, nH2, nHe, nHep, nHepp
    real(dp):: rate

    real(dp):: g_0, g_1, g_2 ! Degeneracy variables
    real(dp):: A_10, A_20, A_21 ! Transition rates (s^-1)
    real(dp):: lam_10, lam_20, lam_21 ! Wavelengths (microns)
    
    g_0 = 1.d0
    g_1 = 3.d0
    g_2 = 5.d0
    
    lam_10 = 88.3323d0
    lam_20 = 32.6523d0
    lam_21 = 51.8004d0
    
    A_10 = 2.597d-5
    A_20 = 3.170d-11
    A_21 = 9.760d-5
    
    ! Call three_level function
    rate = three_level(g_0, g_1, g_2, &
                       lam_10, lam_20, lam_21, &
                       A_10, A_20, A_21, &
                       z, T, &
                       n_ion, ne, nH, nHp, nHe, nHep, nHepp, nH2, &
                       13, 14, 15)
END FUNCTION OIII_fine_structure

FUNCTION CI_fine_structure(T, n_ion, nH, nHp, &
                           ne, nH2, nHe, nHep, &
                           nHepp, z) result(rate)
    implicit none
    real(dp), intent(in):: T, z, n_ion, nH, nHp
    real(dp), intent(in):: ne, nH2, nHe, nHep, nHepp
    real(dp):: rate

    real(dp):: g_0, g_1, g_2 ! Degeneracy variables
    real(dp):: A_10, A_20, A_21 ! Transition rates (s^-1)
    real(dp):: lam_10, lam_20, lam_21 ! Wavelengths (microns)

    g_0 = 1.d0
    g_1 = 3.d0
    g_2 = 5.d0
    
    lam_10 = 609.590d0
    lam_20 = 230.352d0
    lam_21 = 370.269d0
    
    A_10 = 7.930d-8
    A_20 = 1.000d-30
    A_21 = 2.650d-7
    
    ! Call three_level function
    rate = three_level(g_0, g_1, g_2, &
                       lam_10, lam_20, lam_21, &
                       A_10, A_20, A_21, &
                       z, T, &
                       n_ion, ne, nH, nHp, nHe, nHep, nHepp, nH2, &
                       2, 3, 4)
END FUNCTION CI_fine_structure

FUNCTION CII_fine_structure(T, n_ion, nH, nHp, &
                           ne, nH2, nHe, nHep, &
                           nHepp, z) result(rate)
    implicit none
    real(dp), intent(in):: T, z, n_ion, nH, nHp
    real(dp), intent(in):: ne, nH2, nHe, nHep, nHepp
    real(dp):: rate

    real(dp):: g_0, g_1 ! Degeneracy variables
    real(dp):: A_10 ! Transition rates (s^-1)
    real(dp):: lam_10 ! Wavelengths (microns)

    g_0 = 2.d0
    g_1 = 4.d0
    
    lam_10 = 157.636d0
    
    A_10 = 2.290d-6
    
    ! Call two_level function
    rate = two_level(g_0, g_1, &
                    lam_10, &
                    A_10, &
                    z, T, &
                    n_ion, ne, nH, nHp, nHe, nHep, nHepp, nH2, &
                    1)
END FUNCTION CII_fine_structure

FUNCTION NII_fine_structure(T, n_ion, nH, nHp, &
                           ne, nH2, nHe, nHep, &
                           nHepp, z) result(rate)
    implicit none
    real(dp), intent(in):: T, z, n_ion, nH, nHp
    real(dp), intent(in):: ne, nH2, nHe, nHep, nHepp
    real(dp):: rate

    real(dp):: g_0, g_1, g_2 ! Degeneracy variables
    real(dp):: A_10, A_20, A_21 ! Transition rates (s^-1)
    real(dp):: lam_10, lam_20, lam_21 ! Wavelengths (microns)

    g_0 = 1.d0
    g_1 = 3.d0
    g_2 = 5.d0
    
    lam_10 = 205.244d0
    lam_20 = 76.4318d0
    lam_21 = 121.767d0
    
    A_10 = 2.080d-06
    A_20 = 1.000d-30
    A_21 = 7.460d-06
    
    ! Call three_level function
    rate = three_level(g_0, g_1, g_2, &
                       lam_10, lam_20, lam_21, &
                       A_10, A_20, A_21, &
                       z, T, &
                       n_ion, ne, nH, nHp, nHe, nHep, nHepp, nH2, &
                       5, 6, 7)
END FUNCTION NII_fine_structure

FUNCTION SiI_fine_structure(T, n_ion, nH, nHp, &
                           ne, nH2, nHe, nHep, &
                           nHepp, z) result(rate)
    implicit none
    real(dp), intent(in):: T, z, n_ion, nH, nHp
    real(dp), intent(in):: ne, nH2, nHe, nHep, nHepp
    real(dp):: rate

    real(dp):: g_0, g_1, g_2 ! Degeneracy variables
    real(dp):: A_10, A_20, A_21 ! Transition rates (s^-1)
    real(dp):: lam_10, lam_20, lam_21 ! Wavelengths (microns)

    g_0 = 1.d0
    g_1 = 3.d0
    g_2 = 5.d0
    
    lam_10 = 129.641d0
    lam_20 = 44.7993d0
    lam_21 = 68.4548d0
    
    A_10 = 8.250d-6
    A_20 = 3.490d-10
    A_21 = 4.210d-5
    
    ! Call three_level function
    rate = three_level(g_0, g_1, g_2, &
                       lam_10, lam_20, lam_21, &
                       A_10, A_20, A_21, &
                       z, T, &
                       n_ion, ne, nH, nHp, nHe, nHep, nHepp, nH2, &
                       22, 23, 24)
END FUNCTION SiI_fine_structure

FUNCTION SiII_fine_structure(T, n_ion, nH, nHp, &
                           ne, nH2, nHe, nHep, &
                           nHepp, z) result(rate)
    implicit none
    real(dp), intent(in):: T, z, n_ion, nH, nHp
    real(dp), intent(in):: ne, nH2, nHe, nHep, nHepp
    real(dp):: rate

    real(dp):: g_0, g_1 ! Degeneracy variables
    real(dp):: A_10 ! Transition rates (s^-1)
    real(dp):: lam_10 ! Wavelengths (microns)
    
    g_0 = 2.d0
    g_1 = 4.d0
    
    lam_10 = 34.8046d0
    
    A_10 = 2.131d-4
    
    ! Call two_level function
    rate = two_level(g_0, g_1, &
                    lam_10, &
                    A_10, &
                    z, T, &
                    n_ion, ne, nH, nHp, nHe, nHep, nHepp, nH2, &
                    12)
END FUNCTION SiII_fine_structure

FUNCTION NeII_fine_structure(T, n_ion, nH, nHp, &
                           ne, nH2, nHe, nHep, &
                           nHepp, z) result(rate)
    implicit none
    real(dp), intent(in):: T, z, n_ion, nH, nHp
    real(dp), intent(in):: ne, nH2, nHe, nHep, nHepp
    real(dp):: rate

    real(dp):: g_0, g_1 ! Degeneracy variables
    real(dp):: A_10 ! Transition rates (s^-1)
    real(dp):: lam_10 ! Wavelengths (microns)
    
    g_0 = 4.d0
    g_1 = 2.d0
    
    lam_10 = 12.8101d0
    
    A_10 = 8.590d-3
    
    ! Call two_level function
    rate = two_level(g_0, g_1, &
                    lam_10, &
                    A_10, &
                    z, T, &
                    n_ion, ne, nH, nHp, nHe, nHep, nHepp, nH2, &
                    11)
END FUNCTION NeII_fine_structure

FUNCTION FeI_fine_structure(T, n_ion, nH, nHp, &
                           ne, nH2, nHe, nHep, &
                           nHepp, z) result(rate)
    implicit none
    real(dp), intent(in):: T, z, n_ion, nH, nHp
    real(dp), intent(in):: ne, nH2, nHe, nHep, nHepp
    real(dp):: rate

    real(dp):: g_0, g_1, g_2 ! Degeneracy variables
    real(dp):: A_10, A_20, A_21 ! Transition rates (s^-1)
    real(dp):: lam_10, lam_20, lam_21 ! Wavelengths (microns)

    g_0 = 9.d0
    g_1 = 7.d0
    g_2 = 5.d0
    
    lam_10 = 24.0358d0
    lam_20 = 14.2005d0
    lam_21 = 34.7038d0
    
    A_10 = 2.510d-3
    A_20 = 1.000d-30
    A_21 = 1.560d-3
    
    ! Call three_level function
    rate = three_level(g_0, g_1, g_2, &
                       lam_10, lam_20, lam_21, &
                       A_10, A_20, A_21, &
                       z, T, &
                       n_ion, ne, nH, nHp, nHe, nHep, nHepp, nH2, &
                       16, 17, 18)
END FUNCTION FeI_fine_structure

FUNCTION FeII_fine_structure(T, n_ion, nH, nHp, &
                           ne, nH2, nHe, nHep, &
                           nHepp, z) result(rate)
    implicit none
    real(dp), intent(in):: T, z, n_ion, nH, nHp
    real(dp), intent(in):: ne, nH2, nHe, nHep, nHepp
    real(dp):: rate

    real(dp):: g_0, g_1, g_2 ! Degeneracy variables
    real(dp):: A_10, A_20, A_21 ! Transition rates (s^-1)
    real(dp):: lam_10, lam_20, lam_21 ! Wavelengths (microns)

    g_0 = 10.d0
    g_1 = 8.d0
    g_2 = 6.d0
    
    lam_10 = 25.9811d0
    lam_20 = 14.9731d0
    lam_21 = 35.3394d0
    
    A_10 = 2.050d-3
    A_20 = 1.000d-30
    A_21 = 1.560d-3
    
    ! Call three_level function 
    rate = three_level(g_0, g_1, g_2, &
                       lam_10, lam_20, lam_21, &
                       A_10, A_20, A_21, &
                       z, T, &
                       n_ion, ne, nH, nHp, nHe, nHep, nHepp, nH2, &
                       19, 20, 21)
END FUNCTION FeII_fine_structure

FUNCTION SI_fine_structure(T, n_ion, nH, nHp, &
                           ne, nH2, nHe, nHep, &
                           nHepp, z) result(rate)
    implicit none
    real(dp), intent(in):: T, z, n_ion, nH, nHp
    real(dp), intent(in):: ne, nH2, nHe, nHep, nHepp
    real(dp):: rate

    real(dp):: g_0, g_1, g_2 ! Degeneracy variables
    real(dp):: A_10, A_20, A_21 ! Transition rates (s^-1)
    real(dp):: lam_10, lam_20, lam_21 ! Wavelengths (microns)

    g_0 = 5.d0
    g_1 = 3.d0
    g_2 = 1.d0
    
    lam_10 = 25.2421d0
    lam_20 = 17.4278d0
    lam_21 = 56.2957d0
    
    A_10 = 1.400d-3
    A_20 = 7.050d-8
    A_21 = 3.020d-4
    
    ! Call three_level function 
    rate = three_level(g_0, g_1, g_2, &
                       lam_10, lam_20, lam_21, &
                       A_10, A_20, A_21, &
                       z, T, &
                       n_ion, ne, nH, nHp, nHe, nHep, nHepp, nH2, &
                       25, 26, 27)
END FUNCTION SI_fine_structure

FUNCTION dust_recombination_cooling(T, G0, ne, f_dg, nH) result(rate)
    ! Dust recombination cooling
    implicit none
    real(dp), intent(in):: T, G0, ne, f_dg, nH
    real(dp):: rate
    real(dp):: beta_drc

    beta_drc = 0.74d0 / (T**0.068d0)
    rate = (1.5d0 * 4.65d-30) * (T**0.94d0) * ((G0 * sqrt(T) / (0.5*ne))**beta_drc) * ne * 0.5d0 * f_dg * nH
    ! rate = 4.65d-30 * (T**0.94d0) * ((G0 * sqrt(T) / (0.5*ne))**beta_drc) * ne * 0.5d0 * f_dg * nH

END FUNCTION dust_recombination_cooling

FUNCTION dust_recombination_cooling_WD01(T, G0, ne, f_dg, nH) result(rate)
    ! Dust recombination cooling from WD01
    implicit none
    real(dp), intent(in):: T, G0, ne, f_dg, nH
    real(dp):: rate
    real(dp):: Gfac, D0, D1, D2, D3, D4

    Gfac = log(1.7d0 * G0 * sqrt(T) / ne + 50.d0)
    D0 = 0.4535d0
    D1 = 2.234d0
    D2 = -6.266d0
    D3 = 1.442d0
    D4 = 0.05089d0
    rate = 1d-28 * ne * nH * f_dg * (T**(D0 + D1/Gfac)) * exp(D2 + D3*Gfac - D4*Gfac*Gfac)

END FUNCTION dust_recombination_cooling_WD01

FUNCTION dust_gas_collisional_cooling(T, G0, xH2, aexp, nH, f_dg) result(rate)
    implicit none
    real(dp), intent(in):: T, G0, xH2, aexp, nH, f_dg
    real(dp):: rate
    real(dp):: dust_hc_const, T_dust

    dust_hc_const = 1.0d-33 + ((3.8d-33 - 1.0d-33) * max(min(xH2,1.d0),0.d0))

    T_dust = 16.4d0 * ((1.7d0 * G0)**(1.d0/6.d0))
    T_dust = max( T_dust, 2.725d0 * ( (1.d0/aexp) - 1.d0 ) ) ! Limit dust temp minimum to CMB temp

    rate = dust_hc_const * sqrt(T) * (T - T_dust) * ( 1.d0 - ( 0.8d0 * exp(-75.d0/T) ) )
    ! rate = 1.5d0 * rate * nH * nH * f_dg
    rate = rate * nH * nH * f_dg

END FUNCTION dust_gas_collisional_cooling

FUNCTION PE_efficiency(G0, T, ne) result(rate)
    implicit none
    real(dp), intent(in):: G0, T, ne
    real(dp):: rate
    real(dp):: phi_pah, fact

    phi_pah = 0.5d0 ! wolfire+(03)
    fact = G0 * sqrt(T) / (ne * phi_pah)
    rate = (4.9d-2/(1.d0 + 4.d-3 * (fact**0.73d0))) + (3.7d-2 * ((T/1.d4)**0.7d0) / (1.d0 + 2.d-4 * fact))

END FUNCTION PE_efficiency

FUNCTION photoelectric_heating(T, G0, ne, f_dg, nH) result(rate)
    implicit none
    real(dp), intent(in):: T, G0, ne, f_dg, nH
    real(dp):: rate
    real(dp):: eps_PE

    eps_PE = PE_efficiency(G0, T, ne)
    rate = 1.5d0 * 1.3d-24 * eps_PE * G0 * f_dg * nH ! [erg/cm3/s]
    ! rate = 1.3d-24 * eps_PE * G0 * f_dg * nH ! [erg/cm3/s]

END FUNCTION photoelectric_heating

FUNCTION photoelectric_heating_WD01(T, G0, ne, f_dg, nH) result(rate)
    ! This is from weingartner and draine 2001
    implicit none
    real(dp), intent(in):: T, G0, ne, f_dg, nH
    real(dp):: rate
    real(dp):: C0, C1, C2, C3, C4, C5, C6
    real(dp):: Gfac

    C0 = 5.22d0
    C1 = 2.25d0
    C2 = 0.04996d0
    C3 = 0.00430d0
    C4 = 0.147d0
    C5 = 0.431d0
    C6 = 0.692d0

    Gfac = 1.7d0 * G0 * sqrt(T) / ne
    rate = 1.7d-26 * G0 * f_dg * nH * (C0 + C1 * (T**C4))
    rate = rate / (1.d0 + C2 * (Gfac**C5) * (1.d0 + C3 * (Gfac**C6)))

END FUNCTION photoelectric_heating_WD01

FUNCTION cosmic_ray_heating(xe, n_HI, n_HeI, n_H2, ne, xi_h_cr, &
                           element_number_densities, element_ion_fractions) result(rate)
    use rtz_module, only: elements
    use cosmic_ray_ionization_module
    implicit none
    real(dp), intent(in):: xe, n_HI, n_HeI, n_H2, ne, xi_h_cr
    real(dp), intent(in):: element_number_densities(27)
    real(dp), intent(in):: element_ion_fractions(27,27)
    real(dp):: rate
    real(dp):: q_cr, phi_s
    real(dp):: total_cosmic_ray_ionization_rate
    real(dp):: HI_cosmic_ray_ionization_rate_primary
    real(dp):: He_cosmic_ray_ionization_rate_primary
    real(dp):: H2_cosmic_ray_ionization_rate_primary
    real(dp):: induced_UV_heat_scale_fac, x_ion

    integer:: i, j

    q_cr = 6.43d0 * ( 1.d0 + 4.06d0 * sqrt( xe / (0.07d0 + xe) ) ) * EV_2_ERG ! Bialy 2019

    phi_s = secondary_cr_rates(xe)
    total_cosmic_ray_ionization_rate = xi_h_cr * (1.0 + phi_s)

    HI_cosmic_ray_ionization_rate_primary = cosmic_ray_ionization_rates(1,1) * xi_h_cr
    He_cosmic_ray_ionization_rate_primary = cosmic_ray_ionization_rates(2,1) * xi_h_cr
    H2_cosmic_ray_ionization_rate_primary = 2.0d0 * xi_h_cr

    rate = 0.d0

    ! HI
    rate = rate + ( n_HI * HI_cosmic_ray_ionization_rate_primary * q_cr ) ! erg/s/cm^3

    ! HeI
    rate = rate + ( n_HeI * He_cosmic_ray_ionization_rate_primary * q_cr ) ! erg/s/cm^3

    ! H2
    rate = rate + ( n_H2 * H2_cosmic_ray_ionization_rate_primary * q_cr ) ! erg/s/cm^3

    ! Direct heating of the electrons from cosmic rays
    rate = rate + (ne * xi_h_cr * 287.d0 * EV_2_ERG) ! erg/s/cm^3

    ! Cosmic ray heating on metals
    ! --> Here we skip H and He as it's done separately above
    induced_UV_heat_scale_fac = H2_cosmic_ray_ionization_rate_primary / 1d-16

    do i = 3, 27
       ! Skip unused elements
       if (elements(i)%atomic_number .lt. 1) then
          cycle
       end if

       ! Skip elements with too low number density
       if (element_number_densities(i) .lt. MIN_COOL_ION) then
          cycle
       end if

       ! Loop over all ions
       do j = 1,elements(i)%n_ions - 1
          x_ion = element_ion_fractions(i,j)

          ! This is the primary heating term from ionization
          rate = rate + (element_number_densities(i) * x_ion * q_cr * cosmic_ray_ionization_rates(i,j) * total_cosmic_ray_ionization_rate) ! erg/s/cm^3

          ! This is the secondary heating term from induced UV emission
          ! Only happens for the ground state
            if (j .eq. 1) then
                rate = rate + ( element_number_densities(i) * x_ion & 
                                * cosmic_ray_ionization_rates_induced_UV_heat(i) & 
                                * cosmic_ray_ionization_rates_induced_UV(i) * EV_2_ERG & 
                                * induced_UV_heat_scale_fac ) ! erg/s/cm^3
            end if

       end do ! END LOOP OVER IONS

    end do ! END LOOP OVER ELEMENTS

END FUNCTION cosmic_ray_heating

FUNCTION photoheating_UVB(element_number_densities, element_ion_fractions) result(rate)
    ! Photoheating contribution from a UVB
    use rtz_module, only: elements
    use photoionization_UVB_module
    implicit none
    real(dp), intent(in):: element_number_densities(27)
    real(dp), intent(in):: element_ion_fractions(27,27)
    real(dp):: rate

    integer:: i, j
   
    rate = 0.d0

    ! Loop over all elements
    do i = 1, 27
       ! Skip unused elements
       if (elements(i)%atomic_number .lt. 1) then
          cycle
       end if

       ! Skip elements with too low number density
       if (element_number_densities(i) .lt. MIN_COOL_ION) then
          cycle
       end if

       ! Loop over all ions
       do j = 1,elements(i)%n_ions - 1
          rate = rate + (element_number_densities(i) * element_ion_fractions(i,j) * HM12_UVB_z(i,j,2)) ! erg/s/cm^3
       end do ! END LOOP OVER IONS

    end do ! END LOOP OVER ELEMENTS

END FUNCTION photoheating_UVB

FUNCTION photoheating_UVB_G0(G0, element_number_densities, element_ion_fractions) result(rate)
    ! Photoheating from the G0 background
    ! Note that this only impacts the ground state
    use rtz_module, only: elements
    implicit none
    real(dp), intent(in):: G0
    real(dp), intent(in):: element_number_densities(27)
    real(dp), intent(in):: element_ion_fractions(27,27)
    real(dp):: rate

    integer:: i

    rate = 0.d0

    ! Loop over all elements except H and He
    do i = 3, 27
       ! Skip unused elements
       if (elements(i)%atomic_number .lt. 1) then
          cycle
       end if

       ! Skip elements with too low number density
       if (element_number_densities(i) .lt. MIN_COOL_ION) then
          cycle
       end if

       rate = rate + (element_number_densities(i) * element_ion_fractions(i,1) * G0_heating_rates(i) * G0) ! erg/s/cm^3
    
    end do

END FUNCTION photoheating_UVB_G0

FUNCTION Epump(nH, T, xH2, xHI) result(Ep)
    ! Energy converted to heat from UV pumping
    ! see Appendix A http://articles.adsabs.harvard.edu/cgi-bin/nph-iarticle_query?1990ApJ...365..620B&amp;data_type=PDF_HIGH&amp;whole_paper=YES&amp;type=PRINTER&amp;filetype=.pdf
    implicit none
    real(dp), intent(in):: nH, T, xH2, xHI
    real(dp):: Ep
    real(dp):: Crad, Cdex, Cfrac
    Crad = 2.0d-7     !radiation de-excitation rate Burton 1990
    Cdex = (1.0d-12)*((1.4d0*xH2*exp(-18100.d0/(T + 1200.d0))) + (1.d0*xHI*exp(-1000.d0/T)))*sqrt(T)*nH !collisional de-excitation rate Burton 1990
    Cfrac = Cdex / (Cdex + Crad)
    Ep = 2.d0 * Cfrac * EV_2_ERG  !ergs
END FUNCTION Epump

FUNCTION H2_heating_bialy(G0, nH2, nH, T, xH2, xHI, xHII, xe, f_dg, xi_h2_cr) result(rate)
    ! Heating from H2 formation and destruction following Bialy 2018
    use molecules_module, only: alpha_H2_dust, alpha_H2_prim
    implicit none
    real(dp), intent(in):: G0, nH2, nH, T, xH2, xHI, xHII, xe, f_dg, xi_h2_cr
    real(dp):: rate
    real(dp):: D0, I_UV, E_pump, ncrit_factor, n_crit
    real(dp):: Hrate_H2_pump
    real(dp):: E_pd, Hrate_H2_pd
    real(dp):: Hrate_H2_form
    real(dp):: E_form1_dust, E_form2_dust, H2_formation_rate_dust
    real(dp):: E_form1_prim, E_form2_prim, H2_formation_rate_prim

    D0 = 5.68d-11
    I_UV = 1.7d0 * G0

    E_pump = 1.12d0 * EV_2_ERG
    n_crit = 1.1d5 / sqrt(T/1000.d0)
    ncrit_factor = 1.d0 / (1.d0 + (n_crit/nH))

    ! heating from H2 pumping
    Hrate_H2_pump = 9.d0 * D0 * I_UV * E_pump * nH2 * ncrit_factor ! erg/s/cm^3

    ! heating from H2 photodissociation
    E_pd = 0.4d0 * EV_2_ERG
    Hrate_H2_pd = D0 * I_UV * E_pd * nH2 ! erg/s/cm^3

    ! heating from H2 formation
    Hrate_H2_form = 0.d0

    ! heating from H2 formation --> dust channel
    E_form1_dust = 0.2d0 * EV_2_ERG
    E_form2_dust = 4.48d0 * EV_2_ERG * ncrit_factor
    H2_formation_rate_dust = alpha_H2_dust(T, f_dg)

    Hrate_H2_form = Hrate_H2_form + (H2_formation_rate_dust * (E_form1_dust + E_form2_dust) * nH * nH * xHI) ! erg/s/cm^3

    ! heating from H2 formation --> primordial channel
    E_form1_prim = 0.6d0 * EV_2_ERG
    E_form2_prim = 3.13d0 * EV_2_ERG * ncrit_factor
    H2_formation_rate_prim = alpha_H2_prim(T, xe, xi_h2_cr, G0, xHI, xHII)

    Hrate_H2_form = Hrate_H2_form + (H2_formation_rate_prim * (E_form1_prim + E_form2_prim) * nH * nH * xHI) ! erg/s/cm^3
    
    rate = Hrate_H2_pd + Hrate_H2_form + Hrate_H2_pump

END FUNCTION H2_heating_bialy

FUNCTION H2_heating(G0, nH2, nH, T, xH2, xHI, xHII, xe, f_dg, xi_h2_cr) result(rate)
    ! Heating from the formation and destruction of the H2 molecule
    ! This is the same as Katz 2017
    use molecules_module, only: alpha_H2
    implicit none
    real(dp), intent(in):: G0, nH2, nH, T, xH2, xHI, xHII, xe, f_dg, xi_h2_cr
    real(dp):: rate
    real(dp):: fpump, Ebpump, EUV, kUV, HrateLW, H2_formation_rate

    rate = 0.d0

    ! Heating from destruction and pumping
    fpump = 6.94d0 ! Pumping fraction in the ISM: Draine and Bertoldi 1996
    Ebpump = max(Epump(nH, T, xH2, xHI), 0.d0) ! Pumping energy in ergs
    EUV = 0.4d0 * EV_2_ERG ! Energy from photodissociation in ergs (Black and Dalgarno 1977)

    kUV = G0 * 5.68d-11 
    HrateLW = ((kUV*fpump*Ebpump) + (kUV*EUV))*nH2 ! [erg/s/cm3] = [#/s]*[erg/#]*[cm-3]

    rate = rate + max(HrateLW, 0.d0)

    ! Heating from H2 formation
    H2_formation_rate = alpha_H2(T, f_dg, xe, xi_h2_cr, G0, xHI, xHII, nH)
    rate = rate + (2.4d-12 * H2_formation_rate * xHI * nH * nH) ! [cm3 s-1]

END FUNCTION H2_heating

FUNCTION CT_heat_cool(T, element_number_densities, element_ion_fractions) result(rate)
    ! Heating and cooling from charge exchange reactions
    use charge_exchange_module
    implicit none
    real(dp), intent(in):: T
    real(dp), intent(in):: element_number_densities(27)
    real(dp), intent(in):: element_ion_fractions(27,27)
    real(dp):: rate
    real(dp):: loc_nHI, loc_nHII
    
    rate = 0.d0

    loc_nHI = element_number_densities(1) * element_ion_fractions(1,1)
    loc_nHII = element_number_densities(1) * element_ion_fractions(1,2)

    ! Helium
    rate = rate + (charge_transfer_recombination(2, 2, T) * loc_nHI * element_number_densities(2) * element_ion_fractions(2,2) * 10.99d0 * EV_2_ERG) ! H + He+ -> He + H+ 

    ! Carbon
    rate = rate + (charge_transfer_ionization(1, 6, T) * loc_nHII * element_number_densities(6) * element_ion_fractions(6,1) * 2.34d0 * EV_2_ERG) ! H+ + C -> C+ + H
    rate = rate + (charge_transfer_recombination(2, 6, T) * loc_nHI * element_number_densities(6) * element_ion_fractions(6,2) * (-2.34d0) * EV_2_ERG) ! H + C+ -> C + H+
    rate = rate + (charge_transfer_recombination(3, 6, T) * loc_nHI * element_number_densities(6) * element_ion_fractions(6,3) * 4.01d0 * EV_2_ERG) ! H + C++ -> C+ + H+
    rate = rate + (charge_transfer_recombination(4, 6, T) * loc_nHI * element_number_densities(6) * element_ion_fractions(6,4) * 5.73d0 * EV_2_ERG) ! H + C+++ -> C++ + H+
    rate = rate + (charge_transfer_recombination(5, 6, T) * loc_nHI * element_number_densities(6) * element_ion_fractions(6,5) * 11.3d0 * EV_2_ERG) ! H + C++++ -> C+++ + H+

    ! Nitrogen
    rate = rate + (charge_transfer_ionization(1, 7, T) * loc_nHII * element_number_densities(7) * element_ion_fractions(7,1) * (-0.94d0) * EV_2_ERG) ! H+ + N -> N+ + H
    rate = rate + (charge_transfer_recombination(2, 7, T) * loc_nHI * element_number_densities(7) * element_ion_fractions(7,2) * 0.94d0 * EV_2_ERG) ! H + N+ -> N + H+
    rate = rate + (charge_transfer_recombination(3, 7, T) * loc_nHI * element_number_densities(7) * element_ion_fractions(7,3) * 4.56d0 * EV_2_ERG) ! H + N++ -> N+ + H+
    rate = rate + (charge_transfer_recombination(4, 7, T) * loc_nHI * element_number_densities(7) * element_ion_fractions(7,4) * 6.40d0 * EV_2_ERG) ! H + N+++ -> N++ + H+
    rate = rate + (charge_transfer_recombination(5, 7, T) * loc_nHI * element_number_densities(7) * element_ion_fractions(7,5) * 11.0d0 * EV_2_ERG) ! H + N++++ -> N+++ + H+

    ! Oxygen
    rate = rate + (charge_transfer_ionization(1, 8, T) * loc_nHII * element_number_densities(8) * element_ion_fractions(8,1) * (-0.02d0) * EV_2_ERG) ! H+ + O -> O+ + H
    rate = rate + (charge_transfer_recombination(2, 8, T) * loc_nHI * element_number_densities(8) * element_ion_fractions(8,2) * 0.02d0 * EV_2_ERG) ! H + O+ -> O + H+
    rate = rate + (charge_transfer_recombination(3, 8, T) * loc_nHI * element_number_densities(8) * element_ion_fractions(8,3) * 6.65d0 * EV_2_ERG) ! H + O++ -> O+ + H+
    rate = rate + (charge_transfer_recombination(4, 8, T) * loc_nHI * element_number_densities(8) * element_ion_fractions(8,4) * 5.00d0 * EV_2_ERG) ! H + O+++ -> O++ + H+
    rate = rate + (charge_transfer_recombination(5, 8, T) * loc_nHI * element_number_densities(8) * element_ion_fractions(8,5) * 8.47d0 * EV_2_ERG) ! H + O++++ -> O+++ + H+

    ! Neon
    rate = rate + (charge_transfer_recombination(4, 10, T) * loc_nHI * element_number_densities(10) * element_ion_fractions(10,4) * 5.82d0 * EV_2_ERG) ! H + Ne+++ -> Ne++ + H+
    rate = rate + (charge_transfer_recombination(5, 10, T) * loc_nHI * element_number_densities(10) * element_ion_fractions(10,5) * 8.60d0 * EV_2_ERG) ! H + Ne++++ -> Ne+++ + H+

    ! Magnesium
    rate = rate + (charge_transfer_ionization(1, 12, T) * loc_nHII * element_number_densities(12) * element_ion_fractions(12,1) * 1.52d0 * EV_2_ERG) ! H+ + Mg -> Mg+ + H
    rate = rate + (charge_transfer_ionization(2, 12, T) * loc_nHII * element_number_densities(12) * element_ion_fractions(12,2) * (-1.44d0) * EV_2_ERG) ! H+ + Mg+ -> Mg++ + H
    rate = rate + (charge_transfer_recombination(3, 12, T) * loc_nHI * element_number_densities(12) * element_ion_fractions(12,3) * 1.44d0 * EV_2_ERG) ! H + Mg++ -> Mg+ + H+
    rate = rate + (charge_transfer_recombination(4, 12, T) * loc_nHI * element_number_densities(12) * element_ion_fractions(12,4) * 5.73d0 * EV_2_ERG) ! H + Mg+++ -> Mg++ + H+
    rate = rate + (charge_transfer_recombination(5, 12, T) * loc_nHI * element_number_densities(12) * element_ion_fractions(12,5) * 8.60d0 * EV_2_ERG) ! H + Mg++++ -> Mg+++ + H+

    ! Silicon
    rate = rate + (charge_transfer_ionization(1, 14, T) * loc_nHII * element_number_densities(14) * element_ion_fractions(14,1) * 0.12d0 * EV_2_ERG) ! H+ + Si -> Si+ + H
    rate = rate + (charge_transfer_ionization(2, 14, T) * loc_nHII * element_number_densities(14) * element_ion_fractions(14,2) * (-2.72d0) * EV_2_ERG) ! H+ + Si+ -> Si++ + H+
    rate = rate + (charge_transfer_recombination(3, 14, T) * loc_nHI * element_number_densities(14) * element_ion_fractions(14,3) * 2.72d0 * EV_2_ERG) ! H + Si++ -> Si+ + H+
    rate = rate + (charge_transfer_recombination(4, 14, T) * loc_nHI * element_number_densities(14) * element_ion_fractions(14,4) * 4.23d0 * EV_2_ERG) ! H + Si+++ -> Si++ + H+
    rate = rate + (charge_transfer_recombination(5, 14, T) * loc_nHI * element_number_densities(14) * element_ion_fractions(14,5) * 7.49d0 * EV_2_ERG) ! H + Si++++ -> Si+++ + H+

    ! Sulfur
    rate = rate + (charge_transfer_recombination(2, 16, T) * loc_nHI * element_number_densities(16) * element_ion_fractions(16,2) * (-3.24d0) * EV_2_ERG) ! H + S+ -> S + H+
    rate = rate + (charge_transfer_recombination(4, 16, T) * loc_nHI * element_number_densities(16) * element_ion_fractions(16,4) * 5.73d0 * EV_2_ERG) ! H + S+++ -> S++ + H+
    rate = rate + (charge_transfer_recombination(5, 16, T) * loc_nHI * element_number_densities(16) * element_ion_fractions(16,5) * 8.60d0 * EV_2_ERG) ! H + S++++ -> S+++ + H+

    ! Iron
    rate = rate + (charge_transfer_ionization(2, 26, T) * loc_nHII * element_number_densities(26) * element_ion_fractions(26,2) * (-2.56d0) * EV_2_ERG) ! H+ + Fe+ -> Fe++ + H+
    rate = rate + (charge_transfer_recombination(3, 26, T) * loc_nHI * element_number_densities(26) * element_ion_fractions(26,3) * 2.56d0 * EV_2_ERG) ! H + Fe++ -> Fe+ + H+
    rate = rate + (charge_transfer_recombination(4, 26, T) * loc_nHI * element_number_densities(26) * element_ion_fractions(26,4) * 6.30d0 * EV_2_ERG) ! H + Fe+++ -> Fe++ + H+
    rate = rate + (charge_transfer_recombination(5, 26, T) * loc_nHI * element_number_densities(26) * element_ion_fractions(26,5) * 10.0d0 * EV_2_ERG) ! H + Fe++++ -> Fe+++ + H+

END FUNCTION CT_heat_cool

FUNCTION local_photoheating(dNp, element_number_densities, element_ion_fractions, ilevel) result(rate)
    ! Heating due to photoionization from the local radiation field
    use rtz_module, only: elements, n_elements
    use rt_parameters, only: nGroups, PHrate
    implicit none
    real(dp), intent(in):: element_number_densities(27)
    real(dp), intent(in):: element_ion_fractions(27,27)
    real(dp), dimension(nGroups), intent(in):: dNp
    integer, intent(in):: ilevel
    real(dp):: rate
    integer:: ii,jj,kk

    rate = 0.d0
    do ii=1,n_elements
       if (elements(ii)%atomic_number.gt.0) then
          do jj=1,elements(ii)%n_ions - 1
             do kk=1,nGroups
                rate = rate + dNp(kk) * element_number_densities(ii) &
                       * element_ion_fractions(ii,jj) * PHrate(kk,ii,jj)
             end do
          end do
       end if
    end do
END FUNCTION local_photoheating

FUNCTION CO_cooling_koyama_00(n, nH2, nHI, nCO, T) result(rate)
    implicit none
    real(dp), intent(in):: n, nH2, nCO, T, nHI
    real(dp):: rate
    real(dp):: rot, vib_H, vib_H2, k_B, T_pivot

    ! CO cooling

    T_pivot = 3080.d0
    rot = n * nCO * 4.d0*((kB*T)**2.d0)*9.7d-8 / (n * 2.76d0 * kB * (1.d0 + (3.3d6 * (T/1.d3)**(3.d0/4.d0)/n) + 1.5d0*((3.3d6 * (T/1.d3)**(3.d0/4.d0)/n)**0.5d0)))

    vib_H2 = nH2 * nCO * T_pivot * kB * 4.3d-14 * T * EXP(-(3.14d5/T)**0.333d0) * EXP(-T_pivot/T)
    vib_H = nHI * nCO * T_pivot * kB * 3.0d-12 * (T**0.5d0) * EXP(-(2000.d0/T)**3.43d0) * EXP(-T_pivot/T)

    rate = (rot + vib_H2 + vib_H)
END FUNCTION CO_cooling_koyama_00

SUBROUTINE all_cooling(T, ne, aexp, element_number_densities, element_ion_fractions, &
                     nCO, G0, f_dg, xe, xi_h_cr, xi_h2_cr, ss_factor, dNp, ilevel, rate, &
                     saved_cooling_rates, saved_cooling_rates_names)
    ! Main cooling driver
    ! 
    ! T --> Temperature [K]
    ! ne --> Electron number density [cm^-3]
    ! aexp --> Scale factor
    ! element_number_densities --> vector with element number densities [cm^-3]
    ! element_ion_fractions --> 2D array with ion fractions
    ! G0 --> habing band radiation field (MW units)
    ! f_dg --> dust to gas mass ratio normalized by MW value
    ! xe --> electron fraction i.e. ne / (rho/mH)
    ! xi_h_cr --> cosmic ray ionization rate
    use rtz_module, only: elements
    use rt_parameters, only: nGroups, isH2_rtz, rt_advect, rtz_include_charge_exchange, &
                             rtz_include_cosmic_ray_ionization, rtz_include_HM12_UVB
    implicit none
    real(dp), intent(in):: T, ne, aexp, G0, f_dg, xe, xi_h_cr, xi_h2_cr, ss_factor, nCO
    real(dp), intent(in):: element_number_densities(27)
    real(dp), intent(in):: element_ion_fractions(27,27)
    real(dp), dimension(nGroups), intent(in):: dNp
    integer, intent(in):: ilevel
    real(dp), intent(inout):: rate
    real(dp), dimension(50), intent(inout)::saved_cooling_rates
    character(len=20), dimension(50), intent(inout)::saved_cooling_rates_names
    real(dp):: nH_I, nH_II, nH2, nHe_I, nHe_II, nHe_III, nH, xHI, xHII, xH2
    real(dp):: metal_cool_smooth_f1, metal_cool_smooth_f2
    real(dp):: cooling_HI, cooling_HII, cooling_HeI, cooling_HeII
    real(dp):: cooling_HeIII, cooling_bremmstrahlung, cooling_compton
    real(dp):: cooling_H2, total_primordial_cooling
    real(dp):: high_T_metal_cooling
    real(dp):: n_ion_fs, z, total_fine_structure
    real(dp):: cooling_fine_structure_CI, cooling_fine_structure_CII
    real(dp):: cooling_fine_structure_NII
    real(dp):: cooling_fine_structure_OI, cooling_fine_structure_OIII
    real(dp):: cooling_fine_structure_NeII
    real(dp):: cooling_fine_structure_SiI, cooling_fine_structure_SiII
    real(dp):: cooling_fine_structure_SI
    real(dp):: cooling_fine_structure_FeI, cooling_fine_structure_FeII
    real(dp):: dust_cooling, dust_rec_cooling, dust_coll_cooling
    real(dp):: CO_cooling
    real(dp):: photoelectric_heat
    real(dp):: cosmic_ray_heat
    real(dp):: uvb_photoheat
    real(dp):: uvb_photoheat_G0
    real(dp):: h2_heat
    real(dp):: charge_transfer_heat_cool
    real(dp):: photoheating
    real(dp):: total_cooling, total_heating
    integer:: save_cooling_counter

    save_cooling_counter = 1
    saved_cooling_rates_names = ''

    nH_I = element_number_densities(1) * element_ion_fractions(1,1)
    nH_II = element_number_densities(1) * element_ion_fractions(1,2)
    if (isH2_rtz) then 
       nH2 = element_number_densities(1) * element_ion_fractions(1,3) * 0.5d0
    else
       nH2 = 1.d-40
    end if
    nHe_I = element_number_densities(2) * element_ion_fractions(2,1)
    nHe_II = element_number_densities(2) * element_ion_fractions(2,2)
    nHe_III = element_number_densities(2) * element_ion_fractions(2,3)

    nH = nH_I + nH_II + (2.0 * nH2)
    xHI = element_ion_fractions(1,1)
    xHII = element_ion_fractions(1,1)
    if (isH2_rtz) then 
       xH2 = element_ion_fractions(1,2) * 0.5d0 
    else
       xH2 = 1.d-40
    end if

    ! metal cool smoothing parameters
    metal_cool_smooth_f1 = 0.5d0 * (tanh( (5.d-3) * ( T - 1.d4 ) ) + 1.d0 )
    metal_cool_smooth_f2 = 0.5d0 * (tanh( (5.d-3) * ( (-1.d0 * T) + 1.d4 ) ) + 1.d0 )

    ! Cooling from primordial species
    cooling_HI = (collisional_ionization_cooling_HI(T) + collisional_excitation_cooling_HI_seon20(T)) * ne * nH_I
    cooling_HII = recombination_cooling_case_B_HII(T) * ne * nH_II
    cooling_HeI = collisional_ionization_cooling_HeI(T) * ne * nHe_I
    cooling_HeII = (collisional_ionization_cooling_HeII(T) + collisional_excitation_cooling_HeII(T) + recombination_cooling_case_B_HeII(T) + dielectronic_recombination_cooling_HeII(T)) * ne * nHe_II
    cooling_HeIII = recombination_cooling_case_B_HeIII(T) * ne * nHe_III
    cooling_bremmstrahlung = bremmstrahlung(T) * ne * (nH_II + nHe_II + (4.d0 * nHe_III))
    cooling_compton = compton_cooling(T,aexp) * ne
    if (isH2_rtz) then 
       cooling_H2 = H2_cooling(nH_I, nH2, T) 
    !    cooling_H2 = H2_cooling_GA08(T, ne, nH_I, nH_II, nHe_I, nH2)
    !    cooling_H2 = cooling_H2GP(nH_I,nH2,T)
    else
       cooling_H2 = 0.d0
    end if

    total_primordial_cooling = cooling_HI + cooling_HII
    total_primordial_cooling = total_primordial_cooling + cooling_HeI + cooling_HeII + cooling_HeIII
    total_primordial_cooling = total_primordial_cooling + cooling_bremmstrahlung + cooling_compton + cooling_H2

    ! Save cooling rates
    saved_cooling_rates(save_cooling_counter) = cooling_HI; saved_cooling_rates_names(save_cooling_counter) = 'cool_HI'; save_cooling_counter = save_cooling_counter + 1
    saved_cooling_rates(save_cooling_counter) = cooling_HII; saved_cooling_rates_names(save_cooling_counter) = 'cool_HII'; save_cooling_counter = save_cooling_counter + 1
    saved_cooling_rates(save_cooling_counter) = cooling_HeI; saved_cooling_rates_names(save_cooling_counter) = 'cool_HeI'; save_cooling_counter = save_cooling_counter + 1
    saved_cooling_rates(save_cooling_counter) = cooling_HeII; saved_cooling_rates_names(save_cooling_counter) = 'cool_HeII'; save_cooling_counter = save_cooling_counter + 1
    saved_cooling_rates(save_cooling_counter) = cooling_HeIII; saved_cooling_rates_names(save_cooling_counter) = 'cool_HeIII'; save_cooling_counter = save_cooling_counter + 1; 
    saved_cooling_rates(save_cooling_counter) = cooling_bremmstrahlung; saved_cooling_rates_names(save_cooling_counter) = 'cool_bremm'; save_cooling_counter = save_cooling_counter + 1
    saved_cooling_rates(save_cooling_counter) = cooling_compton; saved_cooling_rates_names(save_cooling_counter) = 'cool_compton'; save_cooling_counter = save_cooling_counter + 1
    saved_cooling_rates(save_cooling_counter) = cooling_H2; saved_cooling_rates_names(save_cooling_counter) = 'cool_H2'; save_cooling_counter = save_cooling_counter + 1

    ! Metal line cooling -- high temperature
    high_T_metal_cooling = 0.d0
    high_T_metal_cooling = get_high_t_cooling_rates(T, ne, element_number_densities, &
                                                    element_ion_fractions, metal_cool_smooth_f1)

    ! Save cooling rates
    saved_cooling_rates(save_cooling_counter) = high_T_metal_cooling; saved_cooling_rates_names(save_cooling_counter) = 'high_T_metal_cool'; save_cooling_counter = save_cooling_counter + 1

    ! Metal line cooling -- low temperature (fine structure)
    z = (1.d0 / aexp) - 1.d0
    total_fine_structure = 0.d0

    if (T .lt. 1.1d4) then 
        if (elements(6)%atomic_number .gt. 0) then 
            ! CI cooling
            n_ion_fs = element_number_densities(6) * element_ion_fractions(6,1)
            cooling_fine_structure_CI = 0.d0
            if (n_ion_fs .gt. MIN_COOL_ION) then
                cooling_fine_structure_CI = CI_fine_structure(T, n_ion_fs, nH_I, nH_II, &
                                                            ne, nH2, nHe_I, nHe_II, &
                                                            nHe_III, z)
                total_fine_structure = total_fine_structure + cooling_fine_structure_CI
            end if
            ! Save cooling rates
            saved_cooling_rates(save_cooling_counter) = cooling_fine_structure_CI * metal_cool_smooth_f2; saved_cooling_rates_names(save_cooling_counter) = 'cool_CI'; save_cooling_counter = save_cooling_counter + 1
            
            ! CII cooling
            n_ion_fs = element_number_densities(6) * element_ion_fractions(6,2)
            cooling_fine_structure_CII = 0.d0
            if (n_ion_fs .gt. MIN_COOL_ION) then 
                cooling_fine_structure_CII = CII_fine_structure(T, n_ion_fs, nH_I, nH_II, &
                                                                ne, nH2, nHe_I, nHe_II, &
                                                                nHe_III, z)
                total_fine_structure = total_fine_structure + cooling_fine_structure_CII
            end if
            ! Save cooling rates
            saved_cooling_rates(save_cooling_counter) = cooling_fine_structure_CII * metal_cool_smooth_f2; saved_cooling_rates_names(save_cooling_counter) = 'cool_CII'; save_cooling_counter = save_cooling_counter + 1
        end if
        
        if (elements(7)%atomic_number .gt. 0) then 
            ! NII cooling
            n_ion_fs = element_number_densities(7) * element_ion_fractions(7,2)
            cooling_fine_structure_NII = 0.d0
            if (n_ion_fs .gt. MIN_COOL_ION) then
                cooling_fine_structure_NII = NII_fine_structure(T, n_ion_fs, nH_I, nH_II, &
                                                                ne, nH2, nHe_I, nHe_II, &
                                                                nHe_III, z)
                total_fine_structure = total_fine_structure + cooling_fine_structure_NII
            end if
            ! Save cooling rates
            saved_cooling_rates(save_cooling_counter) = cooling_fine_structure_NII * metal_cool_smooth_f2; saved_cooling_rates_names(save_cooling_counter) = 'cool_NII'; save_cooling_counter = save_cooling_counter + 1
        end if
        
        if (elements(8)%atomic_number .gt. 0) then 
            ! OI cooling
            n_ion_fs = element_number_densities(8) * element_ion_fractions(8,1)
            cooling_fine_structure_OI = 0.d0
            if (n_ion_fs .gt. MIN_COOL_ION) then
                cooling_fine_structure_OI = OI_fine_structure(T, n_ion_fs, nH_I, nH_II, &
                                                            ne, nH2, nHe_I, nHe_II, &
                                                            nHe_III, z)
                total_fine_structure = total_fine_structure + cooling_fine_structure_OI
            end if
            ! Save cooling rates
            saved_cooling_rates(save_cooling_counter) = cooling_fine_structure_OI * metal_cool_smooth_f2; saved_cooling_rates_names(save_cooling_counter) = 'cool_OI'; save_cooling_counter = save_cooling_counter + 1

            ! OIII cooling
            n_ion_fs = element_number_densities(8) * element_ion_fractions(8,3)
            cooling_fine_structure_OIII = 0.d0
            if (n_ion_fs .gt. MIN_COOL_ION) then
                cooling_fine_structure_OIII = OIII_fine_structure(T, n_ion_fs, nH_I, nH_II, &
                                                                ne, nH2, nHe_I, nHe_II, &
                                                                nHe_III, z)
                total_fine_structure = total_fine_structure + cooling_fine_structure_OIII
            end if
            ! Save cooling rates
            saved_cooling_rates(save_cooling_counter) = cooling_fine_structure_OIII * metal_cool_smooth_f2; saved_cooling_rates_names(save_cooling_counter) = 'cool_OIII'; save_cooling_counter = save_cooling_counter + 1
        end if

        if (elements(10)%atomic_number .gt. 0) then 
            ! NeII cooling
            n_ion_fs = element_number_densities(10) * element_ion_fractions(10,2)
            cooling_fine_structure_NeII = 0.d0
            if (n_ion_fs .gt. MIN_COOL_ION) then
                cooling_fine_structure_NeII = NeII_fine_structure(T, n_ion_fs, nH_I, nH_II, &
                                                                ne, nH2, nHe_I, nHe_II, &
                                                                nHe_III, z)
                total_fine_structure = total_fine_structure + cooling_fine_structure_NeII
            end if
            ! Save cooling rates
            saved_cooling_rates(save_cooling_counter) = cooling_fine_structure_NeII * metal_cool_smooth_f2; saved_cooling_rates_names(save_cooling_counter) = 'cool_NeII'; save_cooling_counter = save_cooling_counter + 1
        end if
        
        if (elements(14)%atomic_number .gt. 0) then 
            ! SiI cooling
            n_ion_fs = element_number_densities(14) * element_ion_fractions(14,1)
            cooling_fine_structure_SiI = 0.d0
            if (n_ion_fs .gt. MIN_COOL_ION) then
                cooling_fine_structure_SiI = SiI_fine_structure(T, n_ion_fs, nH_I, nH_II, &
                                                                ne, nH2, nHe_I, nHe_II, &
                                                                nHe_III, z)
                total_fine_structure = total_fine_structure + cooling_fine_structure_SiI
            end if
            ! Save cooling rates
            saved_cooling_rates(save_cooling_counter) = cooling_fine_structure_SiI * metal_cool_smooth_f2; saved_cooling_rates_names(save_cooling_counter) = 'cool_SiI'; save_cooling_counter = save_cooling_counter + 1

            ! SiII cooling
            n_ion_fs = element_number_densities(14) * element_ion_fractions(14,2)
            cooling_fine_structure_SiII = 0.d0
            if (n_ion_fs .gt. MIN_COOL_ION) then
                cooling_fine_structure_SiII = SiII_fine_structure(T, n_ion_fs, nH_I, nH_II, &
                                                                ne, nH2, nHe_I, nHe_II, &
                                                                nHe_III, z)
                total_fine_structure = total_fine_structure + cooling_fine_structure_SiII
            end if
            ! Save cooling rates
            saved_cooling_rates(save_cooling_counter) = cooling_fine_structure_SiII * metal_cool_smooth_f2; saved_cooling_rates_names(save_cooling_counter) = 'cool_SiII'; save_cooling_counter = save_cooling_counter + 1
        end if
        
        if (elements(16)%atomic_number .gt. 0) then 
            ! SI cooling
            n_ion_fs = element_number_densities(16) * element_ion_fractions(16,1)
            cooling_fine_structure_SI = 0.d0
            if (n_ion_fs .gt. MIN_COOL_ION) then
                cooling_fine_structure_SI = SI_fine_structure(T, n_ion_fs, nH_I, nH_II, &
                                                              ne, nH2, nHe_I, nHe_II, &
                                                              nHe_III, z)
                total_fine_structure = total_fine_structure + cooling_fine_structure_SI
            end if
            ! Save cooling rates
            saved_cooling_rates(save_cooling_counter) = cooling_fine_structure_SI * metal_cool_smooth_f2; saved_cooling_rates_names(save_cooling_counter) = 'cool_SI'; save_cooling_counter = save_cooling_counter + 1
        end if
        
        if (elements(26)%atomic_number .gt. 0) then 
            ! FeI cooling
            n_ion_fs = element_number_densities(26) * element_ion_fractions(26,1)
            cooling_fine_structure_FeI = 0.d0
            if (n_ion_fs .gt. MIN_COOL_ION) then 
                cooling_fine_structure_FeI = FeI_fine_structure(T, n_ion_fs, nH_I, nH_II, &
                                                                ne, nH2, nHe_I, nHe_II, &
                                                                nHe_III, z)
                total_fine_structure = total_fine_structure + cooling_fine_structure_FeI
            end if 
            ! Save cooling rates
            saved_cooling_rates(save_cooling_counter) = cooling_fine_structure_FeI * metal_cool_smooth_f2; saved_cooling_rates_names(save_cooling_counter) = 'cool_FeI'; save_cooling_counter = save_cooling_counter + 1
            
            ! FeII cooling
            n_ion_fs = element_number_densities(26) * element_ion_fractions(26,2)
            cooling_fine_structure_FeII = 0.d0
            if (n_ion_fs .gt. MIN_COOL_ION) then
                cooling_fine_structure_FeII = FeII_fine_structure(T, n_ion_fs, nH_I, nH_II, &
                                                                ne, nH2, nHe_I, nHe_II, &
                                                                nHe_III, z)
                total_fine_structure = total_fine_structure + cooling_fine_structure_FeII
            end if  
            ! Save cooling rates
            saved_cooling_rates(save_cooling_counter) = cooling_fine_structure_FeII * metal_cool_smooth_f2; saved_cooling_rates_names(save_cooling_counter) = 'cool_FeII'; save_cooling_counter = save_cooling_counter + 1
        end if   
    end if

    ! Smooth the fine structure cooling with temeprature if needed
    total_fine_structure = total_fine_structure * metal_cool_smooth_f2

    ! Dust cooling
    dust_rec_cooling = dust_recombination_cooling(T, G0, ne, f_dg, element_number_densities(1))
    ! dust_rec_cooling = dust_recombination_cooling_WD01(T, G0, ne, f_dg, element_number_densities(1))
    dust_coll_cooling =  dust_gas_collisional_cooling(T, G0, xH2*2.d0, aexp, element_number_densities(1), f_dg) 
    dust_cooling = dust_rec_cooling + dust_coll_cooling
    
    ! Save cooling rates
    saved_cooling_rates(save_cooling_counter) = dust_rec_cooling; saved_cooling_rates_names(save_cooling_counter) = 'cool_dust_rec'; save_cooling_counter = save_cooling_counter + 1
    saved_cooling_rates(save_cooling_counter) = dust_coll_cooling; saved_cooling_rates_names(save_cooling_counter) = 'cool_dust_col'; save_cooling_counter = save_cooling_counter + 1

#ifdef CO
    CO_cooling = CO_cooling_koyama_00(nH, nH2, nH_I, nCO, T)
    saved_cooling_rates(save_cooling_counter) = CO_cooling; saved_cooling_rates_names(save_cooling_counter) = 'cool_CO'; save_cooling_counter = save_cooling_counter + 1
#endif

    ! Photoelectric heating --> note factor of 1.7 is because IUV
    photoelectric_heat = photoelectric_heating(T, G0, ne, f_dg, element_number_densities(1))
    ! photoelectric_heat = photoelectric_heating_WD01(T, G0, ne, f_dg, element_number_densities(1))

    ! Save cooling rates
    saved_cooling_rates(save_cooling_counter) = photoelectric_heat; saved_cooling_rates_names(save_cooling_counter) = 'heat_PE'; save_cooling_counter = save_cooling_counter + 1

    ! Cosmic ray heating
    cosmic_ray_heat = cosmic_ray_heating(xe, nH_I, nHe_I, nH2, &
                                         ne, xi_h_cr, &
                                         element_number_densities, &
                                         element_ion_fractions)

    ! Save cooling rates
    saved_cooling_rates(save_cooling_counter) = cosmic_ray_heat; saved_cooling_rates_names(save_cooling_counter) = 'heat_CR'; save_cooling_counter = save_cooling_counter + 1

    ! Photoheating from the UV background
    uvb_photoheat = photoheating_UVB(element_number_densities, element_ion_fractions)
    uvb_photoheat = uvb_photoheat * ss_factor ! Account for self-shielding

    ! Save cooling rates
    saved_cooling_rates(save_cooling_counter) = uvb_photoheat; saved_cooling_rates_names(save_cooling_counter) = 'heat_UVB'; save_cooling_counter = save_cooling_counter + 1

    ! Photoheating from the G0 FUV background
    uvb_photoheat_G0 = photoheating_UVB_G0(G0, element_number_densities,element_ion_fractions)

    ! Save cooling rates
    saved_cooling_rates(save_cooling_counter) = uvb_photoheat_G0; saved_cooling_rates_names(save_cooling_counter) = 'heat_G0'; save_cooling_counter = save_cooling_counter + 1

    ! Heating from H2 formation and destruction
    ! h2_heat = H2_heating(G0, nH2, nH, T, xH2, xHI, xHII, xe, f_dg, xi_h2_cr)
    h2_heat = H2_heating_bialy(G0, nH2, nH, T, xH2, xHI, xHII, xe, f_dg, xi_h2_cr)

    ! Save cooling rates
    saved_cooling_rates(save_cooling_counter) = h2_heat; saved_cooling_rates_names(save_cooling_counter) = 'heat_H2'; save_cooling_counter = save_cooling_counter + 1

    ! Heating and cooling from charge transfer reactions
    charge_transfer_heat_cool = CT_heat_cool(T, element_number_densities, element_ion_fractions)

    ! Save cooling rates
    saved_cooling_rates(save_cooling_counter) = charge_transfer_heat_cool; saved_cooling_rates_names(save_cooling_counter) = 'heat_CT'; save_cooling_counter = save_cooling_counter + 1

    ! Photoheating from the local radiation field
    if (rt_advect) then
       photoheating = local_photoheating(dNp, element_number_densities, element_ion_fractions, ilevel)

       ! Save cooling rates
       saved_cooling_rates(save_cooling_counter) = photoheating; saved_cooling_rates_names(save_cooling_counter) = 'heat_PH'; save_cooling_counter = save_cooling_counter + 1
    end if

    !////////////////////////////////////////////////////
    !//           Calculate Heatint & Cooling          //
    !////////////////////////////////////////////////////
    
    !!!!!! Sum all of the cooling rates
    total_cooling = total_primordial_cooling + dust_cooling + high_T_metal_cooling + total_fine_structure 
#ifdef CO
    total_cooling = total_cooling + CO_cooling
#endif

    !!!!!! Sum all of the heating rates
    total_heating = photoelectric_heat + uvb_photoheat_G0

    if (rtz_include_cosmic_ray_ionization) then
       total_heating = total_heating + cosmic_ray_heat
    end if

    if (rtz_include_charge_exchange) then
       total_heating = total_heating + charge_transfer_heat_cool
    end if

    if (rtz_include_HM12_UVB) then 
       total_heating = total_heating + uvb_photoheat
    end if

    if (isH2_rtz) then
       total_heating = total_heating + h2_heat
    end if

    if (rt_advect) then
       total_heating = total_heating + photoheating
    end if

    rate = total_heating - total_cooling

END SUBROUTINE all_cooling

end module rtz_coolrates_module