MODULE coolrates_module
  ! Module for returning cooling and thermochemistry interaction rates.
  ! The temperature dependence is tabulated in rate (nonlog) versus the
  ! log of T, because it is expensive to calculate on the fly.
  ! The rates are interpolated using cubic splines, and extrapolated in
  ! log-log space if temperature is above table boundaries.
  ! Joki Rosdahl and Andreas Bleuler, September 2015.
  ! Nickerson, Teyssier, and Rosdahl, March 2018.

  use amr_parameters,only:dp
  implicit none

  private   ! default
  public init_coolrates_tables, update_coolrates_tables                  &
       , inp_coolrates_table, compCoolrate                               &
       , tbl_alphaZ_H2, tbl_alphaGP_H2, tbl_alphaA_HII, tbl_alphaA_HeII  &
       , tbl_alphaA_HeIII, tbl_alphaB_HII, tbl_alphaB_HeII               &
       , tbl_alphaB_HeIII, tbl_beta_HI, tbl_beta_HeI, tbl_beta_HeII      &
       , tbl_cr_ci_HI, tbl_cr_ci_HeI, tbl_cr_ci_HeII, tbl_cr_ce_HI       &
       , tbl_cr_ce_HeI, tbl_cr_ce_HeII, tbl_cr_r_HII, tbl_cr_r_HeII      &
       , tbl_cr_r_HeIII, tbl_cr_bre, tbl_cr_com, tbl_cr_die              &
       , tbl_beta_H2HI, tbl_beta_H2H2, tbl_cr_H2HI, tbl_cr_H2H2          &
       , tbl_beta_H3B

  ! Default cooling rates table parameters
  integer,parameter     :: nbinT  = 1001
  real(dp),parameter    :: Tmin   = 1d-2
  real(dp),parameter    :: Tmax   = 1d+9
  real(dp)              :: dlogTinv ! Inverse of the bin space (in K)
  real(dp)              :: hTable, h2Table, h3Table   ! Interpol constants
  real(dp)              :: one_over_lnTen, one_over_hTable, one_over_h2Table
  real(dp)              :: three_over_h2Table, two_over_h3Table

  real(dp),dimension(nbinT) :: T_lookup = 0d0 ! Lookup temperature in log K

  type coolrates_table
     ! Cooling and interaction rates (log):
     real(dp),dimension(nbinT)::rates  = 0d0
     ! Temperature derivatives of those rates (drate/dlog(T)):
     real(dp),dimension(nbinT)::primes = 0d0
  end type coolrates_table

  type(coolrates_table),save::tbl_alphaZ_H2  ! H2 dust grain formation
  type(coolrates_table),save::tbl_alphaGP_H2 ! H2 gas-phase formation
  type(coolrates_table),save::tbl_alphaA_HII ! Case A rec. coefficients
  type(coolrates_table),save::tbl_alphaA_HeII
  type(coolrates_table),save::tbl_alphaA_HeIII
  type(coolrates_table),save::tbl_alphaB_HII ! Case B rec. coefficients
  type(coolrates_table),save::tbl_alphaB_HeII
  type(coolrates_table),save::tbl_alphaB_HeIII
  type(coolrates_table),save::tbl_beta_HI ! Collisional ionisation rates
  type(coolrates_table),save::tbl_beta_HeI
  type(coolrates_table),save::tbl_beta_HeII
  type(coolrates_table),save::tbl_cr_ci_HI ! Coll. ionisation cooling
  type(coolrates_table),save::tbl_cr_ci_HeI
  type(coolrates_table),save::tbl_cr_ci_HeII
  type(coolrates_table),save::tbl_cr_ce_HI ! Coll. excitation cooling
  type(coolrates_table),save::tbl_cr_ce_HeI
  type(coolrates_table),save::tbl_cr_ce_HeII
  type(coolrates_table),save::tbl_cr_r_HII ! Recombination cooling
  type(coolrates_table),save::tbl_cr_r_HeII
  type(coolrates_table),save::tbl_cr_r_HeIII
  type(coolrates_table),save::tbl_cr_bre ! Bremsstrahlung cooling rates
  type(coolrates_table),save::tbl_cr_com ! Compton cooling rates
  type(coolrates_table),save::tbl_cr_die ! Dielectronic cooling rates
  type(coolrates_table),save::tbl_beta_H2HI ! Coll. dissociation rates
  type(coolrates_table),save::tbl_beta_H2H2 ! Coll. dissociation rates
  type(coolrates_table),save::tbl_beta_H3B ! Coll. dissociation rates
  type(coolrates_table),save::tbl_cr_H2HI ! Collisional diss. cooling
  type(coolrates_table),save::tbl_cr_H2H2 ! Collisional diss. cooling

CONTAINS

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE init_coolrates_tables(aexp)

! Initialise the cooling rates tables.
!-------------------------------------------------------------------------
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer :: ierr
#endif
  real(dp) :: aexp
  integer :: myid, ncpu, iT
!-------------------------------------------------------------------------
#ifndef WITHOUTMPI
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,ncpu,ierr)
#endif
#ifdef WITHOUTMPI
  myid=0
  ncpu=1
#endif

  ! Initialise the table lookup temperatures -----------------------------
  do iT=1, nbinT
     T_lookup(iT) = log10(Tmin) + (dble(iT)-1d0) / (dble(nbinT)-1d0)     &
                                * (log10(Tmax)-log10(Tmin))
  end do
  dlogTinv = dble(nbinT-1)/(T_lookup(nbinT)-T_lookup(1)) ! (space)^-1
  hTable = 1d0/dlogTinv                             !
  h2Table = hTable*hTable                           ! Constants for table
  h3Table = h2Table*hTable                          ! interpolation

  one_over_lnTen = 1d0/log(10d0)
  one_over_hTable = 1d0/hTable
  one_over_h2Table = 1d0/h2Table
  three_over_h2Table = 3d0/h2Table
  two_over_h3Table = 2d0/h3Table

  do iT = myid+1, nbinT, ncpu ! Loop over TK and assign rates
     call comp_table_rates(iT,aexp)
  end do ! end TK loop

  ! Distribute the complete table between cpus ---------------------------
#ifndef WITHOUTMPI
  call mpi_distribute_coolrates_table(tbl_alphaZ_H2)
  call mpi_distribute_coolrates_table(tbl_alphaGP_H2)
  call mpi_distribute_coolrates_table(tbl_alphaA_HII)
  call mpi_distribute_coolrates_table(tbl_alphaA_HeII)
  call mpi_distribute_coolrates_table(tbl_alphaA_HeIII)

  call mpi_distribute_coolrates_table(tbl_alphaB_HII)
  call mpi_distribute_coolrates_table(tbl_alphaB_HeII)
  call mpi_distribute_coolrates_table(tbl_alphaB_HeIII)

  call mpi_distribute_coolrates_table(tbl_beta_HI)
  call mpi_distribute_coolrates_table(tbl_beta_HeI)
  call mpi_distribute_coolrates_table(tbl_beta_HeII)

  call mpi_distribute_coolrates_table(tbl_cr_ci_HI)
  call mpi_distribute_coolrates_table(tbl_cr_ci_HeI)
  call mpi_distribute_coolrates_table(tbl_cr_ci_HeII)

  call mpi_distribute_coolrates_table(tbl_cr_ce_HI)
  call mpi_distribute_coolrates_table(tbl_cr_ce_HeI)
  call mpi_distribute_coolrates_table(tbl_cr_ce_HeII)

  call mpi_distribute_coolrates_table(tbl_cr_r_HII)
  call mpi_distribute_coolrates_table(tbl_cr_r_HeII)
  call mpi_distribute_coolrates_table(tbl_cr_r_HeIII)

  call mpi_distribute_coolrates_table(tbl_cr_bre)
  call mpi_distribute_coolrates_table(tbl_cr_com)
  call mpi_distribute_coolrates_table(tbl_cr_die)

  call mpi_distribute_coolrates_table(tbl_beta_H2HI)
  call mpi_distribute_coolrates_table(tbl_beta_H2H2)
  call mpi_distribute_coolrates_table(tbl_beta_H3B)
  call mpi_distribute_coolrates_table(tbl_cr_H2HI)
  call mpi_distribute_coolrates_table(tbl_cr_H2H2)

#endif

  if(myid==0) print*,'Coolrates tables initialised '

END SUBROUTINE init_coolrates_tables

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE update_coolrates_tables(aexp)
! Update cooling rates lookup tables which depend on aexp
!-------------------------------------------------------------------------
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer :: ierr
#endif
  real(dp) :: aexp
  integer:: myid, ncpu, iT
!-------------------------------------------------------------------------
#ifndef WITHOUTMPI
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,ncpu,ierr)
#endif
#ifdef WITHOUTMPI
  myid=0
  ncpu=1
#endif
  tbl_cr_com%rates  = 0d0 ; tbl_cr_com%primes = 0d0
  do iT = myid+1, nbinT, ncpu ! Loop over TK and assign rates
     call update_table_rates(iT, aexp)
  end do ! end TK loop

  ! Distribute the complete table between cpus ---------------------------
#ifndef WITHOUTMPI
  call mpi_distribute_coolrates_table(tbl_cr_com)
#endif

  if(myid==0) print*,'Coolrates table updated'
END SUBROUTINE update_coolrates_tables

#ifndef WITHOUTMPI
!PRIVATEXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE mpi_distribute_coolrates_table(table)
! Distribute table between all cpus, assuming table contains only partial
! entries on each cpu, but the whole table is acquired by summing those
! partial tables
!-------------------------------------------------------------------------
  use mpi_mod
  implicit none
  type(coolrates_table)::table
  real(dp),dimension(:),allocatable :: table_mpi_sum
  integer::ierr
!-------------------------------------------------------------------------
  allocate(table_mpi_sum(nbinT))
  call MPI_ALLREDUCE(table%rates,table_mpi_sum                           &
                  ,nbinT,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  table%rates = table_mpi_sum
  call MPI_ALLREDUCE(table%primes,table_mpi_sum                          &
                  ,nbinT,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  table%primes = table_mpi_sum
  deallocate(table_mpi_sum)

END SUBROUTINE mpi_distribute_coolrates_table
#endif

!PRIVATEXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE comp_table_rates(iT, aexp)
! Fill in index iTK in all rates tables.
!-------------------------------------------------------------------------
  use rt_parameters,only:rt_OTSA
  use constants,only:kB
  implicit none
  integer::iT
  real(dp)::aexp, T, T2, Ta, T5, lambda, f, hf, laHII, laHeII, laHeIII
  real(dp)::lowrleft,lowrright,lowr_hi,lowr_h2
  real(dp)::lowvleft_hi,lowvright_hi,lowv_hi,lowv_h2
  real(dp)::TT
!-------------------------------------------------------------------------
  ! Rates are stored in non-log, while temperature derivatives (primes)
  ! are stored in dRate/dlogT (= dRate/dT * T * ln(10)).

  ! The log-log primes are just the normal primes times T/rate,
  ! i.e. dlogL/dlogT = T/L dL/dT

  T = 10d0**T_lookup(iT)

  ! Creation rate of H2 on dust [cm3 s-1] (functio of T from--------------
  ! Hollenback & McKee 1979, average rate between Habart 2004 dense PDRs
  ! and Jura 1974 diffuse ISM)
  ! When using, must multiply by Z*f_dust
  T2 = T/1d2
  f = (1.0+0.4*sqrt(T2)+0.2*T2+0.08*T2**2)
  tbl_alphaZ_H2%rates(iT)   = 9.0d-17 * sqrt(T2) / f
  tbl_alphaZ_H2%primes(iT)  = (-0.12*T2**2 - 0.1*T2 + 0.5) / f           &
                            * log(10d0) * tbl_alphaZ_H2%rates(iT)

  ! Gas phase formation rate rate for H2 on H- [cm3 s-1] assuming---------
  ! equilibrium abundances for H-, as explained in the Appendix of
  ! McKee and Krumholz (2012)
  tbl_alphaGP_H2%rates(iT)  = 8.0d-19*(T/1000.0)**0.88
  tbl_alphaGP_H2%primes(iT) = 0.88                                       &
                            * log(10d0) * tbl_alphaGP_H2%rates(iT)
  ! Case A rec. coefficient [cm3 s-1] for HII (Hui&Gnedin'97)-------------
  lambda = 315614./T                                ! 2d0 * 157807d0 / T
  f = 1d0+(lambda/0.522)**0.47
  tbl_alphaA_HII%rates(iT)  =  1.269d-13 * lambda**1.503 / f**1.923
  tbl_alphaA_HII%primes(iT) = ( 0.90381*(f-1.)/f - 1.503 )               &
                            * log(10d0) * tbl_alphaA_HII%rates(iT)

  ! Case A rec. coefficient [cm3 s-1] for HeII (Hui&Gnedin'97)------------
  lambda = 570670./T
  tbl_alphaA_HeII%rates(iT)  = 3d-14 * lambda**0.654
  tbl_alphaA_HeII%primes(iT) = -0.654                                    &
                             * log(10d0) * tbl_alphaA_HeII%rates(iT)

  ! Case A rec. coefficient [cm3 s-1] for HeIII (Hui&Gnedin'97)-----------
  lambda =  1263030./T
  f= 1d0+(lambda/0.522)**0.47
  tbl_alphaA_HeIII%rates(iT)  =  2.538d-13 * lambda**1.503 / f**1.923
  tbl_alphaA_HeIII%primes(iT) =  ( 0.90381*(f-1.)/f - 1.503 )            &
                              * log(10d0) * tbl_alphaA_HeIII%rates(iT)

  ! Case B rec. coefficient [cm3 s-1] for HII (Hui&Gnedin'97)-------------
  lambda = 315614./T
  f= 1d0+(lambda/2.74)**0.407
  tbl_alphaB_HII%rates(iT)  = 2.753d-14 * lambda**1.5 / f**2.242
  tbl_alphaB_HII%primes(iT) = ( 0.912494*(f-1.)/f - 1.5 )                &
                             * log(10d0) * tbl_alphaB_HII%rates(iT)

  ! Case B rec. coefficient [cm3 s-1] for HeII (Hui&Gnedin'97)------------
  lambda = 570670./T
  tbl_alphaB_HeII%rates(iT)  = 1.26d-14 * lambda**0.75
  tbl_alphaB_HeII%primes(iT) =  -0.75                                    &
                             * log(10d0) * tbl_alphaB_HeII%rates(iT)

  ! Case B rec. coefficient [cm3 s-1] for HeIII (Hui&Gnedin'97)-----------
  lambda = 1263030./T
  f= 1d0+(lambda/2.74)**0.407
  tbl_alphaB_HeIII%rates(iT)  = 5.506d-14 * lambda**1.5 / f**2.242
  tbl_alphaB_HeIII%primes(iT) = ( 0.912494*(f-1.)/f - 1.5 )              &
                              * log(10d0) * tbl_alphaB_HeIII%rates(iT)

  ! Collisional ionization rate [cm3 s-1] of HI (Maselli&'03)-------------
  T5 = T/1d5
  f = 1d0+sqrt(T5) ; hf=0.5d0/f
  tbl_beta_HI%rates(iT)  = 5.85d-11 * sqrt(T) / f * exp(-157809.1d0/T)
  tbl_beta_HI%primes(iT) = (hf+157809.1d0/T)                             &
                         * log(10d0) * tbl_beta_HI%rates(iT)

  ! Collisional ionization rate [cm3 s-1] of HeI (Maselli&'03)------------
  tbl_beta_HeI%rates(iT)  = 2.38d-11 * sqrt(T) / f * exp(-285335.4d0/T)
  tbl_beta_HeI%primes(iT) = (hf+285335.4d0/T)                            &
                          * log(10d0) * tbl_beta_HeI%rates(iT)

  ! Collisional ionization rate [cm3 s-1] of HeII (Maselli&'03)-----------
  tbl_beta_HeII%rates(iT)  = 5.68d-12 * sqrt(T) / f * exp(-631515d0/T)
  tbl_beta_HeII%primes(iT) = (hf+631515d0/T)                            &
                           * log(10d0) * tbl_beta_HeII%rates(iT)

  ! BEGIN COOLING RATES-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  T5 = T/1d5
  f = 1d0+sqrt(T5) ; hf=0.5d0/f

  ! Coll. Ionization Cooling from Cen 1992 (via Maselli et al 2003)
  tbl_cr_ci_HI%rates(iT)  = 1.27d-21 * sqrt(T) / f * exp(-157809.1/T)
  tbl_cr_ci_HI%primes(iT) = (hf+157809.1/T)                              &
                          * log(10d0) * tbl_cr_ci_HI%rates(iT)

  tbl_cr_ci_HeI%rates(iT)  = 9.38d-22 * sqrt(T) / f * exp(-285335.4/T)
  tbl_cr_ci_HeI%primes(iT) = (hf+285335.4/T)                             &
                           * log(10d0) * tbl_cr_ci_HeI%rates(iT)

  tbl_cr_ci_HeII%rates(iT)  = 4.95d-22 * sqrt(T) / f * exp(-631515. /T)
  tbl_cr_ci_HeII%primes(iT) = (hf+631515.0/T)                            &
                            * log(10d0) * tbl_cr_ci_HeII%rates(iT)

  ! Collisional excitation cooling from Cen'92
  tbl_cr_ce_HI%rates(iT)  = 7.5d-19 / f * exp(-118348./T)
  tbl_cr_ce_HI%primes(iT) = (118348./T - 0.5d0 * sqrt(T5) / f )          &
                          * log(10d0) * tbl_cr_ce_HI%rates(iT)

  tbl_cr_ce_HeI%rates(iT)  = 9.10d-27 * T**(-0.1687) / f * exp(-13179./T)
  tbl_cr_ce_HeI%primes(iT) = (13179./T - 0.1687 - 0.5d0 * sqrt(T5) / f)  &
                           * log(10d0) * tbl_cr_ce_HeI%rates(iT)

  tbl_cr_ce_HeII%rates(iT) = 5.54d-17 * T**(-0.397)  / f * exp(-473638./T)
  tbl_cr_ce_HeII%primes(iT) = (473638./T - 0.397 - 0.5d0 * sqrt(T5) / f) &
                            * log(10d0) * tbl_cr_ce_HeII%rates(iT)

  ! Recombination Cooling (Hui&Gnedin'97)
  laHII    = 315614./T
  laHeII   = 570670./T
  laHeIII  = 1263030./T
  if(.not. rt_otsa) then ! Case A
     f = 1d0+(laHII/0.541)**0.502
     tbl_cr_r_HII%rates(iT)    = 1.778d-29 * laHII**1.965 / f**2.697 * T
     tbl_cr_r_HII%primes(iT)   = (-0.965 + 1.35389*(f-1.)/f)             &
                               * log(10d0) * tbl_cr_r_HII%rates(iT)

     tbl_cr_r_HeII%rates(iT)   = 3d-14 * laHeII**0.654 * kB * T
     tbl_cr_r_HeII%primes(iT)  = 0.346                                   &
                               * log(10d0) * tbl_cr_r_HeII%rates(iT)

     f = 1d0+(laHeIII/0.541)**0.502
     tbl_cr_r_HeIII%rates(iT) = 14.224d-29 * laHeIII**1.965 / f**2.697 * T
     tbl_cr_r_HeIII%primes(iT)= (-0.965 + 1.35389*(f-1.)/f)              &
                              * log(10d0) * tbl_cr_r_HeIII%rates(iT)
  else ! Case B
     f = 1d0+(laHII/2.25)**0.376
     tbl_cr_r_HII%rates(iT)    = 3.435d-30 * laHII**1.97 / f**3.72 * T
     tbl_cr_r_HII%primes(iT)   = (-0.97 + 1.39827*(f-1.)/f)              &
                               * log(10d0) * tbl_cr_r_HII%rates(iT)

     tbl_cr_r_HeII%rates(iT)   = 1.26d-14 * laHeII**0.75 * kB * T
     tbl_cr_r_HeII%primes(iT)  = 0.25                                    &
                               * log(10d0) * tbl_cr_r_HeII%rates(iT)

     f = 1d0+(laHeIII/2.25)**0.376
     tbl_cr_r_HeIII%rates(iT)  = 27.48d-30 * laHeIII**1.97 / f**3.72 * T
     tbl_cr_r_HeIII%primes(iT) = (-0.97 + 1.39827*(f-1.)/f)              &
                               * log(10d0) * tbl_cr_r_HeIII%rates(iT)
  endif

  ! Bremsstrahlung from Osterbrock & Ferland 2006
  tbl_cr_bre%rates(iT)  = 1.42d-27 * 1.5 * sqrt(T)
  tbl_cr_bre%primes(iT) = 0.5                                            &
                        * log(10d0) * tbl_cr_bre%rates(iT)

  ! Compton Cooling from Haimann et al. 96, via Maselli et al.
  ! Need to make sure this is done whenever the redshift changes!
  Ta     = 2.727/aexp
  tbl_cr_com%rates(iT)   = 1.017d-37 * Ta**4 * (T-Ta)
  tbl_cr_com%primes(iT)  = T / (T-Ta)                                    &
                         * log(10d0) * tbl_cr_com%rates(iT)

  ! Dielectronic recombination cooling, from Black 1981
  f = 1.24d-13*T**(-1.5d0)*exp(-470000d0/T)
  tbl_cr_die%rates(iT) = f*(1d0+0.3d0*exp(-94000d0/T))
  tbl_cr_die%primes(iT)=0d0
  if(tbl_cr_die%rates(iT) .gt. 0d0) then ! Can simplify w algebra
     tbl_cr_die%primes(iT) = (tbl_cr_die%rates(iT)*(564000.-1.5*T)       &
                              - f*94000.) /T**2 * T * log(10d0)
  endif

  !BEGIN H2 STUFF --------------------------------------------------------

  ! Collisional dissociation ground state (Dove&Mandy 1986) [cm3 s-1]
  tbl_Beta_H2HI%rates(iT) =                                              &
       7.073d-19*(T**2.012)*exp(-5.179d4/T)/(1d0+2.130d-5*T)**3.512
  ! Col dissociation ground state (Martin&Keogh&Mandy 1998)  [cm3 s-1]
  tbl_Beta_H2H2%rates(iT) = &
       5.996d-30*(T**4.1881)*exp(-5.466d4/T)/(1d0+6.761d-6*T)**5.6881
  ! Three body H2 formation 3H -> H2 + H (Forrey 2013) [cm6 s-1]
  ! For H2+2H -> 2H2 use Beta_H3B/8 (Palla 1983)
  tbl_Beta_H3B%rates(iT) = &
       6d-32*T**(-0.25)+2d-31*T**(-0.5)
  ! Prefer to use numerical prime, rather than analytic:
  TT=T*1.0001
  tbl_Beta_H2HI%primes(iT) =                                             &
       7.073d-19*(TT**2.012)*exp(-5.179d4/TT)/(1d0+2.130d-5*TT)**3.512
  tbl_Beta_H2HI%primes(iT) = &
       (tbl_Beta_H2HI%primes(iT) - tbl_Beta_H2HI%rates(iT)) &
       / (0.0001*T) * T * log(10d0) ! last two terms for log-log
  tbl_Beta_H2H2%primes(iT) = &
      5.996d-30*(TT**4.1881)*exp(-5.466d4/TT)/(1d0+6.761d-6*TT)**5.6881
  tbl_Beta_H2H2%primes(iT) = &
       (tbl_Beta_H2H2%primes(iT) - tbl_Beta_H2H2%rates(iT))              &
       / (0.0001*T) * T * log(10d0) ! last two terms for log-log
  tbl_Beta_H3B%primes(iT) = &
       6d-32*TT**(-0.25)+2d-31*TT**(-0.5)
  tbl_Beta_H3B%primes(iT) = &
       (tbl_Beta_H3B%primes(iT) - tbl_Beta_H3B%rates(iT))              &
       / (0.0001*T) * T * log(10d0) ! last two terms for log-log

  ! Collisional H2 cooling************************************************
  ! Hallenbach McKee (1979) + Halle Combes (2012) in cgs
  lowrleft=1.25*exp(-1.70d2/T)*2.35d-14
  lowrright=1.75*exp(-5.05d2/T)*6.97d-14
  lowr_hi=lowrleft*gamma_hi(T,2d0)+lowrright*gamma_hi(T,3d0)
  lowr_h2=lowrleft*gamma_h2(T,2d0)+lowrright*gamma_h2(T,3d0)
  lowvleft_hi=exp(-5860.0/T)*8.09d-13*1.0d-12*sqrt(T)*exp(-1.0d3/T)
  lowvright_hi=exp(-11720.0/T)*1.6d-12*1.6d-12*sqrt(T)*exp(-1.0*(4.0d2/T)**2)
  lowv_hi=lowvleft_hi+lowvright_hi
  lowv_h2=exp(-5860.0/T)*8.09d-13*1.4d-12*sqrt(T)*exp(-1.2d4/(T+1.2d3))

  tbl_cr_H2HI%rates(iT)=lowr_hi+lowv_hi
  tbl_cr_H2H2%rates(iT)=lowr_h2+lowv_h2

  ! Collisional H2 cooling derivatives (calculated numerically)
  TT=T*1.0001
  lowrleft=1.25*exp(-1.70d2/TT)*2.35d-14
  lowrright=1.75*exp(-5.05d2/TT)*6.97d-14
  lowr_hi=lowrleft*gamma_hi(TT,2d0)+lowrright*gamma_hi(TT,3d0)
  lowr_h2=lowrleft*gamma_h2(TT,2d0)+lowrright*gamma_h2(TT,3d0)
  lowvleft_hi=exp(-5860.0/TT)*8.09d-13*1.0d-12*sqrt(TT)*exp(-1.0d3/TT)
  lowvright_hi=exp(-11720.0/TT)*1.6d-12*1.6d-12*sqrt(TT)*exp(-1.0*(4.0d2/TT)**2)
  lowv_hi=lowvleft_hi+lowvright_hi
  lowv_h2=exp(-5860.0/TT)*8.09d-13*1.4d-12*sqrt(TT)*exp(-1.2d4/(TT+1.2d3))

  tbl_cr_H2HI%primes(iT)=lowr_hi+lowv_hi
  tbl_cr_H2HI%primes(iT) = &
       (tbl_cr_H2HI%primes(iT) - tbl_cr_H2HI%rates(iT)) &
       / (0.0001*T) * T * log(10d0) ! last two terms for log-log

  tbl_cr_H2H2%primes(iT)=lowr_h2+lowv_h2
  tbl_cr_H2H2%primes(iT) = &
       (tbl_cr_H2H2%primes(iT) - tbl_cr_H2H2%rates(iT)) &
       / (0.0001*T) * T * log(10d0) ! last two terms for log-log

  !END H2 STUFF ----------------------------------------------------------

END SUBROUTINE comp_table_rates

!PRIVATEXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE update_table_rates(iT,aexp)
! Update index iTK in all compton cooling rates tables due to change in
! aexp
!-------------------------------------------------------------------------
  implicit none
  integer::iT
  real(dp)::aexp, T, Ta
!-------------------------------------------------------------------------
  ! Rates are stored in log, while temperature derivatives (primes) are
  ! stored in non-log
  T = 10d0**T_lookup(iT)
  ! Compton Cooling from Haimann et al. 96, via Maselli et al.
  Ta     = 2.727/aexp
  tbl_cr_com%rates(iT)   = 1.017d-37 * Ta**4 * (T-Ta)
  tbl_cr_com%primes(iT)  = T / (T-Ta)                                    &
                         * log(10d0) * tbl_cr_com%rates(iT)

END SUBROUTINE update_table_rates

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
FUNCTION inp_coolrates_table(rates_table, T, allow_negative, retPrime)
! Returns TABULATED rate value from given table
! rates        => Rates (and primes) table to interpolate
! T            => Temperature [K]
! retPrime     <= Temperature derivative of rate at T (optional).
!-------------------------------------------------------------------------
  implicit none
  type(coolrates_table)::rates_table
  !real(dp),dimension(nbinTK)::rates, ratesPrime
  real(dp),intent(in)::T
  logical,intent(in)::allow_negative
  real(dp),optional::retPrime
  real(dp)::inp_coolrates_table
  integer,save:: iT = 1
  real(dp),save:: facT, yy, yy2, yy3, fa, fb, fprimea, fprimeb, Tlast=-1
  real(dp),save:: alpha, beta, gamma
  logical,save::extrap
!-------------------------------------------------------------------------
  if (.not. (T .eq. Tlast)) then    ! Reuse index if same T from last call
     ! Log of T, snapped to table at the lower boundary, but allowed
     ! to go above upper boundary, in which case we use extrapolation:
     facT = MAX( log10(T), T_lookup(1) )
     extrap=.false.
     if(facT .gt. T_lookup(nbinT)) extrap=.true. ! Above upper limit
     ! Lower closest index in table:
     iT = MIN(MAX(int((facT-T_lookup(1))*dlogTinv)+1, 1), nbinT-1)
     yy=facT-T_lookup(iT)  ! Dist., in log(T), from T to lower table index
     yy2=yy*yy             ! That distance squared
     yy3=yy2*yy            ! ...and cubed
     Tlast = T
  endif

  if(extrap) then ! TK above upper table limit, so extrapolate in log-log:
     alpha = (log10(rates_table%rates(nbinT))             &
              -log10(rates_table%rates(nbinT-1)) )        &
           / (T_lookup(nbinT)-T_lookup(nbinT-1))
     inp_coolrates_table = 10d0**(log10(rates_table%rates(nbinT))        &
                                  + alpha * (facT - T_lookup(nbinT)))
     if( present(retPrime) )                                             &
          retPrime = alpha * inp_coolrates_table / T
     return
  endif

  fa = rates_table%rates(iT)           !
  fb = rates_table%rates(iT+1)         !  Values at neighbouring table
  fprimea=rates_table%primes(iT)       !  indexes
  fprimeb=rates_table%primes(iT+1)     !

  ! Spline interpolation:
  alpha = fprimea
  beta =(fb-fa) * three_over_h2Table - (2d0*fprimea+fprimeb) * one_over_hTable
  gamma = (fprimea+fprimeb) * one_over_h2Table - (fb-fa) * two_over_h3Table
  inp_coolrates_table = fa+alpha*yy+beta*yy2+gamma*yy3
  if(.not. allow_negative) then
     ! Only allow positive rates (spline can go negative):
     inp_coolrates_table = max(0d0,inp_coolrates_table)
  endif
  if( present(retPrime) )                                                &
       retPrime = (alpha+2d0*beta*yy+3d0*gamma*yy2) / T * one_over_lnTen
END FUNCTION inp_coolrates_table

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
FUNCTION compCoolrate(T, ne, nN, nI, dcooldT)

! Compute cooling rate in a cell, using interpolation from cooling rate
! tables
! T        => Cell emperature [K]
! nN       => Neutral and molecular abundances for each ionisation species
! nI       => Ionized and dissociated abundances for each species
! dcooldT  <= Temperature derivative of the rate
! dcooldx  <= Ionized fraction derivative of the rate
! returns:    Cooling rate [erg s-1 cm-3]
!-------------------------------------------------------------------------
  use rt_parameters,only: nIons, isH2, isHe, ixHI, ixHII, ixHeII, ixHeIII
  implicit none
  real(dp)::T, ne
  real(dp),dimension(nIons)::nN,nI
  real(dp)::compCoolrate, dcooldT !---------------------------------------
  real(dp),save::ci_HI,ci_HeI,ci_HeII
  real(dp),save::ci_HI_prime,ci_HeI_prime,ci_HeII_prime
  real(dp),save::cr_H2HI,cr_H2H2,cr_H2HI_prime,cr_H2H2_prime
  real(dp),save::ce_HI,ce_HeI,ce_HeII
  real(dp),save::ce_HI_prime,ce_HeI_prime,ce_HeII_prime
  real(dp),save::r_HII,r_HeII,r_HeIII
  real(dp),save::r_HII_prime,r_HeII_prime,r_HeIII_prime
  real(dp),save::bre, brefac, bre_prime, com, com_prime, die, die_prime
!-------------------------------------------------------------------------
  ! Coll. Ionization Cooling
  ci_HI   = inp_coolrates_table(tbl_cr_ci_HI, T, .false., ci_HI_prime)   &
          * ne * nN(ixHII) ! ne * nHI
  ci_HI_prime = ci_HI_prime * ne * nN(ixHII)
  ! Collisional excitation cooling
  ce_HI   = inp_coolrates_table(tbl_cr_ce_HI, T, .false., ce_HI_prime)   &
          * ne * nN(ixHII) ! ne * nHI
  ce_HI_prime    = ce_HI_prime * ne * nN(ixHII)
  ! Recombination Cooling
  r_HII = inp_coolrates_table(tbl_cr_r_HII, T, .false., r_HII_prime)     &
        * ne * nI(ixHII) ! ne * nHII
  r_HII_prime   = r_HII_prime * ne * nI(ixHII)

  if(isHe) then
     ! Coll. Ionization Cooling
     ci_HeI      = inp_coolrates_table(tbl_cr_ci_HeI, T, .false.         &
                                                       , ci_HeI_prime)   &
                 * ne * nN(ixHeII)  ! ne * nHeI
     ci_HeI_prime= ci_HeI_prime  * ne * nN(ixHeII)
     ci_HeII     = inp_coolrates_table(tbl_cr_ci_HeII, T, .false.        &
                                                        , ci_HeII_prime) &
                 * ne * nN(ixHeIII) ! ne * nHeII
     ci_HeII_prime= ci_HeII_prime * ne * nN(ixHeIII)
     ! Collisional excitation cooling
     ce_HeI      = inp_coolrates_table(tbl_cr_ce_HeI, T, .false.         &
                                                       , ce_HeI_prime)   &
                 * ne * nN(ixHeII) ! ne * nHeI
     ce_HeI_prime= ce_HeI_prime *ne*nN(ixHeII)
     ce_HeII     = inp_coolrates_table(tbl_cr_ce_HeII, T, .false.        &
                                                        , ce_HeII_prime) &
                 * ne * nN(ixHeIII) ! ne * nHeII
     ce_HeII_prime= ce_HeII_prime * ne * nN(ixHeIII)
     ! Recombination Cooling
     r_HeII      = inp_coolrates_table(tbl_cr_r_HeII, T, .false.         &
                                                       , r_HeII_prime)   &
                 * ne * nI(ixHeII) ! ne * nHeII
     r_HeII_prime= r_HeII_prime * ne * nI(ixHeII)
     r_HeIII     = inp_coolrates_table(tbl_cr_r_HeIII, T, .false.        &
                                                        , r_HeIII_prime) &
                 * ne * nI(ixHeIII) ! ne * nHeIII
     r_HeIII_prime= r_HeIII_prime * ne * nI(ixHeIII)
  endif

  ! Bremsstrahlung
  bre    = inp_coolrates_table(tbl_cr_bre, T, .false., bre_prime)
  brefac = ne * nI(ixHII)
  if(isHe) then
     brefac = brefac + ne * ( nI(ixHeII) + 4. * nI(ixHeIII) )
  endif
  bre        = bre * brefac
  bre_prime  = bre_prime * brefac

  ! Compton Cooling
  com       = inp_coolrates_table(tbl_cr_com, T, .true., com_prime) * ne
  com_prime = com_prime * ne

  ! Dielectronic recombination cooling
  if(isHe) then
     die       = inp_coolrates_table(tbl_cr_die, T, .false., die_prime) &
               * ne * nN(ixHeIII) ! ne * nHeII
     die_prime = die_prime * ne * nN(ixHeIII)
  endif

  if(isH2) then
     ! Collisional dissociation (H2) cooling
     cr_H2HI       = inp_coolrates_table(tbl_cr_H2HI, T, .false.         &
                                                       , cr_H2HI_prime)  &
                   * nN(ixHI) * nN(ixHII) ! nH2 * nHI
     cr_H2HI_prime = cr_H2HI_prime * nN(ixHI)*nN(ixHII)

     cr_H2H2       = inp_coolrates_table(tbl_cr_H2H2, T, .false.         &
                                                       , cr_H2H2_prime)  &
                   * nN(ixHI)**2 ! nH2 * nH2
     cr_H2H2_prime = cr_H2H2_prime * nN(ixHI)**2
  endif

  ! Overall Cooling
  compCoolrate          = ci_HI    + r_HII    + ce_HI    + com + bre
  if(isHe) compCoolrate = compCoolrate                                   &
                        + ci_HeI   + r_HeII   + ce_HeI                   &
                        + ci_HeII  + r_HeIII  + ce_HeII  + die
  if(isH2) compCoolrate = compCoolrate                                   &
                        + cr_H2HI  + cr_H2H2

  dCooldT= ci_HI_prime + r_HII_prime + ce_HI_prime + com_prime + bre_prime
  if(isHe) dCooldT = dCooldT                                             &
                   + ci_HeI_prime  + r_HeII_prime  + ce_HeI_prime        &
                   + ci_HeII_prime + r_HeIII_prime + ce_HeII_prime       &
                   + die_prime
  if(isH2) dCooldT = dCooldT                                             &
                   + cr_H2HI_prime + cr_H2H2_prime

END FUNCTION compCoolrate

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION gamma_hi(T,J)
  !-------------------------------------------------------------------------
  implicit none
  real(dp),intent(in)::T,J
  real(dp)::gamma_hi, T3, jconsts
  !-------------------------------------------------------------------------
  T3 = T/1d3
  jconsts=0.33+0.9*exp(-1.0d0*((J-3.5d0)/0.9d0)**2)
  gamma_hi = &
       jconsts*(1.0d-11*sqrt(T3)/(1.0d0+60.0d0*T3**(-4)) + 1.0d-12*T3)

END FUNCTION gamma_hi
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION gamma_h2(T,J)
  !-------------------------------------------------------------------------
  implicit none
  real(dp),intent(in)::T,J
  real(dp)::gamma_h2, T3, jconsts
  !----------------------------------------------------------
  T3 = T/1d3
  jconsts=(0.276*J**2)*exp(-1.0d0*(J/3.18d0)**1.7)
  gamma_h2 = jconsts*(3.3d-12 +6.6d-12*T3)

END FUNCTION gamma_h2

END MODULE coolrates_module
