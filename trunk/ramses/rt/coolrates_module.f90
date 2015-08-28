MODULE coolrates_module
  ! Module for returning cooling and thermochemistry interaction rates.
  ! The temperature dependence is tabulated, because it is expensive to
  ! calculate on the fly.
  
  ! Extrapolate if temperature is above table boundaries.

  use amr_parameters,only:dp
  implicit none

  ! Everything is public for the time being, but good idea to make public
  ! only that which is needed from the outside
  !private   ! default
  !public init_coolrates_table

  ! Default cooling rates table parameters
  integer,parameter     :: nbinT  = 1001
  real(dp),parameter    :: Tmin   = 1d-2
  real(dp),parameter    :: Tmax   = 1d+9
  real(dp)              :: dlogTinv ! Inverse of the bin space (in K)
  real(dp)              :: hTable, h2Table, h3Table   ! Interpol constants
  
  real(dp),dimension(nbinT) :: T_lookup = 0d0 ! Lookup temperature in log K

  type coolrates_table
     ! Cooling and interaction rates (log):
     real(dp),dimension(nbinT)::rates  = 0d0
     ! Temperature derivatives of those rates (dlog/dlog):
     real(dp),dimension(nbinT)::primes = 0d0 
  end type coolrates_table

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
  
CONTAINS

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE init_coolrates_tables(aexp)

! Initialise the cooling rates tables.
!-------------------------------------------------------------------------
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  real(dp) :: aexp, T
  integer :: myid, ncpu, ierr, iT
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
  
  do iT = myid+1, nbinT, ncpu ! Loop over TK and assign values
     call comp_table_rates(iT,aexp)
  end do ! end TK loop
  
  ! Distribute the complete table between cpus ---------------------------
#ifndef WITHOUTMPI
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
#endif

  if(myid==0) print*,'Coolrates tables initialised '
  do iT = 1,nbinT
     !if(myid==0) write(*,901),10d0**T_lookup(iT)     &
         ! , 10d0**tbl_alphaA_HII%rates(iT)         &
         ! , 10d0**tbl_alphaA_HeII%rates(iT)        &
         ! , 10d0**tbl_alphaA_HeIII%rates(iT)       &
         ! , 10d0**tbl_alphaB_HII%rates(iT)         &
         ! , 10d0**tbl_alphaB_HeII%rates(iT)        &
         ! , 10d0**tbl_alphaB_HeIII%rates(iT)       &
         ! , 10d0**tbl_beta_HI%rates(iT)            &
         ! , 10d0**tbl_beta_HeI%rates(iT)           &
         ! , 10d0**tbl_beta_HeII%rates(iT)    &      

         ! , tbl_alphaA_HII%primes(iT)         &
         ! , tbl_alphaA_HeII%primes(iT)        &
         ! , tbl_alphaA_HeIII%primes(iT)       &
         ! , tbl_alphaB_HII%primes(iT)         &
         ! , tbl_alphaB_HeII%primes(iT)        &
         ! , tbl_alphaB_HeIII%primes(iT)       &
         ! , tbl_beta_HI%primes(iT)            &
         ! , tbl_beta_HeI%primes(iT)           &
         ! , tbl_beta_HeII%primes(iT)          

         ! , 10d0**tbl_cr_ci_HI%rates(iT)         &
         ! , 10d0**tbl_cr_ci_HeI%rates(iT)        &
         ! , 10d0**tbl_cr_ci_HeII%rates(iT)       &
         ! , 10d0**tbl_cr_ce_HI%rates(iT)         &
         ! , 10d0**tbl_cr_ce_HeI%rates(iT)        &
         ! , 10d0**tbl_cr_ce_HeII%rates(iT)       &
         ! , 10d0**tbl_cr_r_HII%rates(iT)            &
         ! , 10d0**tbl_cr_r_HeII%rates(iT)           &
         ! , 10d0**tbl_cr_r_HeIII%rates(iT)          

         ! , tbl_cr_ci_HI%primes(iT)         &
         ! , tbl_cr_ci_HeI%primes(iT)        &
         ! , tbl_cr_ci_HeII%primes(iT)       &
         ! , tbl_cr_ce_HI%primes(iT)         &
         ! , tbl_cr_ce_HeI%primes(iT)        &
         ! , tbl_cr_ce_HeII%primes(iT)       &
         ! , tbl_cr_r_HII%primes(iT)            &
         ! , tbl_cr_r_HeII%primes(iT)           &
         ! , tbl_cr_r_HeIII%primes(iT)          

         ! , 10d0**tbl_cr_bre%rates(iT)         &
         ! , 10d0**tbl_cr_com%rates(iT)        &
         ! , 10d0**tbl_cr_die%rates(iT)       &
         ! , tbl_cr_bre%primes(iT)         &
         ! , tbl_cr_com%primes(iT)        &
         ! , tbl_cr_die%primes(iT)       

  end do
901 format (20(1pe12.3))

END SUBROUTINE init_coolrates_tables

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE update_coolrates_tables(aexp)
! Update cooling rates lookup tables which depend on aexp
!-------------------------------------------------------------------------
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  real(dp) :: aexp
  integer:: myid, ncpu, ierr, iT
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
  do iT = myid+1, nbinT, ncpu ! Loop over TK and assign values
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
  implicit none
  include 'mpif.h'
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
  implicit none
  integer::iT
  real(dp)::aexp, T, Ta, T5, lambda, f, hf, laHII, laHeII, laHeIII
  real(dp),parameter::kb=1.3806d-16        ! Boltzmann constant [ergs K-1]
!-------------------------------------------------------------------------
  ! Rates are stored in log, while temperature derivatives (primes) are
  ! stored in non-log.

  ! The log-log primes are just the normal primes times T/rate,
  ! i.e. dlogL/dlogT = T/L dL/dT
  
  T = 10d0**T_lookup(iT)

  ! Case A rec. coefficient [cm3 s-1] for HII (Hui&Gnedin'97)-------------
  lambda = 315614./T                                ! 2.d0 * 157807.d0 / T
  f = 1.d0+(lambda/0.522)**0.47
  tbl_alphaA_HII%rates(iT)  =  1.269d-13 * lambda**1.503 / f**1.923
  tbl_alphaA_HII%primes(iT) = ( 0.90381*(f-1.)/f - 1.503 )

  ! Case A rec. coefficient [cm3 s-1] for HeII (Hui&Gnedin'97)------------
  lambda = 570670./T
  tbl_alphaA_HeII%rates(iT)  =  3.d-14 * lambda**0.654
  tbl_alphaA_HeII%primes(iT) =  -0.654

  ! Case A rec. coefficient [cm3 s-1] for HeIII (Hui&Gnedin'97)-----------
  lambda =  1263030./T
  f= 1.d0+(lambda/0.522)**0.47
  tbl_alphaA_HeIII%rates(iT)  =  2.538d-13 * lambda**1.503 / f**1.923 
  tbl_alphaA_HeIII%primes(iT) =  ( 0.90381*(f-1.)/f - 1.503 )

  ! Case B rec. coefficient [cm3 s-1] for HII (Hui&Gnedin'97)-------------
  lambda = 315614./T
  f= 1.d0+(lambda/2.74)**0.407
  tbl_alphaB_HII%rates(iT)  = 2.753d-14 * lambda**1.5 / f**2.242
  tbl_alphaB_HII%primes(iT) =  ( 0.912494*(f-1.)/f - 1.5 )

  ! Case B rec. coefficient [cm3 s-1] for HeII (Hui&Gnedin'97)------------
  lambda = 570670./T
  tbl_alphaB_HeII%rates(iT)  = 1.26d-14 * lambda**0.75
  tbl_alphaB_HeII%primes(iT) =  -0.75

  ! Case B rec. coefficient [cm3 s-1] for HeIII (Hui&Gnedin'97)-----------
  lambda = 1263030./T
  f= 1.d0+(lambda/2.74)**0.407
  tbl_alphaB_HeIII%rates(iT) = 5.506d-14 * lambda**1.5 / f**2.242
  tbl_alphaB_HeIII%primes(iT) = ( 0.912494*(f-1.)/f - 1.5 )

  ! Collisional ionization rate [cm3 s-1] of HI (Maselli&'03)-------------
  T5 = T/1d5
  f = 1d0+sqrt(T5) ; hf=0.5d0/f
  tbl_beta_HI%rates(iT)  = 5.85d-11 * sqrt(T) / f * exp(-157809.1d0/T)
  tbl_beta_HI%primes(iT) = (hf+157809.1d0/T)

  ! Collisional ionization rate [cm3 s-1] of HeI (Maselli&'03)------------
  tbl_beta_HeI%rates(iT)  = 2.38d-11 * sqrt(T) / f * exp(-285335.4d0/T)
  tbl_beta_HeI%primes(iT) = (hf+285335.4d0/T)

  ! Collisional ionization rate [cm3 s-1] of HeII (Maselli&'03)-----------
  tbl_beta_HeII%rates(iT)  = 5.68d-12 * sqrt(T) / f * exp(-631515.d0/T)
  tbl_beta_HeII%primes(iT) = (hf+631515.d0/T)

  ! BEGIN COOLING RATES-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  T5 = T/1d5
  f = 1d0+sqrt(T5) ; hf=0.5d0/f

  ! Coll. Ionization Cooling from Cen 1992 (via Maselli et al 2003)
  tbl_cr_ci_HI%rates(iT)  = 1.27d-21 * sqrt(T) / f * exp(-157809.1/T)
  tbl_cr_ci_HI%primes(iT) = hf+157809.1/T

  tbl_cr_ci_HeI%rates(iT)  = 9.38d-22 * sqrt(T) / f * exp(-285335.4/T)
  tbl_cr_ci_HeI%primes(iT) = hf+285335.4/T

  tbl_cr_ci_HeII%rates(iT)  = 4.95d-22 * sqrt(T) / f * exp(-631515. /T)
  tbl_cr_ci_HeII%primes(iT) = hf+631515.0/T

  ! Collisional excitation cooling from Cen'92
  tbl_cr_ce_HI%rates(iT)  = 7.5d-19 / f * exp(-118348./T)
  tbl_cr_ce_HI%primes(iT) = 118348./T - 0.5d0 * sqrt(T5) / f

  tbl_cr_ce_HeI%rates(iT)  = 9.10d-27 * T**(-0.1687) / f * exp(-13179./T) 
  tbl_cr_ce_HeI%primes(iT) = 13179./T - 0.1687 - 0.5d0 * sqrt(T5) / f

  tbl_cr_ce_HeII%rates(iT) = 5.54d-17 * T**(-0.397)  / f * exp(-473638./T)
  tbl_cr_ce_HeII%primes(iT) = 473638./T - 0.397 - 0.5d0 * sqrt(T5) / f
  
  ! Recombination Cooling (Hui&Gnedin'97)
  laHII    = 315614./T                                             
  laHeII   = 570670./T                                            
  laHeIII  = 1263030./T
  if(.not. rt_otsa) then ! Case A
     f = 1.d0+(laHII/0.541)**0.502                         
     tbl_cr_r_HII%rates(iT)    = 1.778d-29 * laHII**1.965 / f**2.697 * T
     tbl_cr_r_HII%primes(iT)   = - 0.965 + 1.35389*(f-1.)/f
 
     tbl_cr_r_HeII%rates(iT)   = 3.d-14 * laHeII**0.654 * kb * T
     tbl_cr_r_HeII%primes      = 0.346

     f = 1.d0+(laHeIII/0.541)**0.502                      
     tbl_cr_r_HeIII%rates(iT) = 14.224d-29 * laHeIII**1.965 / f**2.697 * T
     tbl_cr_r_HeIII%primes(iT)= - 0.965 + 1.35389*(f-1.)/f
  else ! Case B
     f = 1.d0+(laHII/2.25)**0.376                        
     tbl_cr_r_HII%rates(iT)    = 3.435d-30 * laHII**1.97 / f**3.72 * T
     tbl_cr_r_HII%primes(iT)   = - 0.97 + 1.39827*(f-1.)/f

     tbl_cr_r_HeII%rates(iT)   = 1.26d-14 * laHeII**0.75 * kb * T
     tbl_cr_r_HeII%primes(iT)  = 0.25

     f = 1.d0+(laHeIII/2.25)**0.376                       
     tbl_cr_r_HeIII%rates(iT)  = 27.48d-30 * laHeIII**1.97 / f**3.72 * T
     tbl_cr_r_HeIII%primes(iT) = - 0.97 + 1.39827*(f-1.)/f
  endif

  ! Bremsstrahlung from Osterbrock & Ferland 2006
  tbl_cr_bre%rates(iT)  = 1.42d-27 * 1.5 * sqrt(T)
  tbl_cr_bre%primes(iT) = 0.5

  ! Compton Cooling from Haimann et al. 96, via Maselli et al.
  ! Need to make sure this is done whenever the redshift changes!!!!!!!!!!
  Ta     = 2.727/aexp                          
  tbl_cr_com%rates(iT)   = 1.017d-37 * Ta**4 * (T-Ta)
  tbl_cr_com%primes(iT)  = T / (T-Ta)   

  ! Dielectronic recombination cooling, from Black 1981
  f = 1.24d-13*T**(-1.5d0)*exp(-470000.d0/T)
  tbl_cr_die%rates(iT) = f*(1.D0+0.3d0*exp(-94000.d0/T))
  tbl_cr_die%primes(iT) = (tbl_cr_die%rates(iT)*(564000.-1.5*T) &
       - f*94000.) /T**2 * T/ tbl_cr_die%rates(iT)  ! Can fix/simplify w algebra

  call log_table(tbl_alphaA_HII, iT)
  call log_table(tbl_alphaA_HeII, iT) ; call log_table(tbl_alphaA_HeIII, iT)
  call log_table(tbl_alphaB_HII, iT)
  call log_table(tbl_alphaB_HeII, iT) ; call log_table(tbl_alphaB_HeIII, iT)
  call log_table(tbl_beta_HI, iT)
  call log_table(tbl_beta_HeI, iT)    ; call log_table(tbl_beta_HeII, iT)
  call log_table(tbl_cr_ci_HI, iT)
  call log_table(tbl_cr_ci_HeI, iT)   ; call log_table(tbl_cr_ci_HeII, iT)
  call log_table(tbl_cr_ce_HI, iT)
  call log_table(tbl_cr_ce_HeI, iT)   ; call log_table(tbl_cr_ce_HeII, iT)
  call log_table(tbl_cr_r_HII, iT)
  call log_table(tbl_cr_r_HeII, iT)   ; call log_table(tbl_cr_r_HeIII, iT)
  call log_table(tbl_cr_bre, iT)
  call log_table(tbl_cr_com, iT)      ; call log_table(tbl_cr_die, iT)

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
  tbl_cr_com%primes(iT)  = T / (T-Ta)   
  call log_table(tbl_cr_com, iT)

END SUBROUTINE update_table_rates

!PRIVATEXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE log_table(table, iT)
  ! Convert rates table entries to logarithm and make corresponding
  ! arrangements for the log-log derivative  
  !---------------------------------------------------------------------
  implicit none
  type(coolrates_table)::table
  integer::iT
  !---------------------------------------------------------------------
  table%rates(iT) = log10(MAX(table%rates(iT),1d-99))
  if( table%rates(iT) .le. -99 ) table%primes(iT)=0d0  
END SUBROUTINE log_table

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
FUNCTION inp_coolrates_table(rates_table, T, retPrime)
! Returns TABULATED rate value from given table
! rates        => Rates (and primes) table to interpolate
! T            => Temperature [K]
! retPrime     <= Temperature derivative of rate at T (optional).  
!-------------------------------------------------------------------------
  implicit none
  type(coolrates_table)::rates_table
  !real(dp),dimension(nbinTK)::rates, ratesPrime
  real(dp),intent(in)::T
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

  if(extrap) then ! TK is above upper table limit, so we extrapolate:
     alpha = (rates_table%rates(nbinT)-rates_table%rates(nbinT-1))        &
           / (T_lookup(nbinT)-T_lookup(nbinT-1))
     inp_coolrates_table = &
          10d0**(rates_table%rates(nbinT)                                 &
          + alpha * (facT - T_lookup(nbinT)))
     if( present(retPrime) )                                              &
          retPrime = alpha * inp_coolrates_table / T
     return
  endif 
  
  fa = rates_table%rates(iT)   !
  fb = rates_table%rates(iT+1) !  Values at neighbouring table
  fprimea=rates_table%primes(iT)       !  indexes
  fprimeb=rates_table%primes(iT+1)     !

  ! Spline interpolation:
  alpha = fprimea
  beta = 3d0*(fb-fa)/h2Table-(2d0*fprimea+fprimeb)/hTable
  gamma = (fprimea+fprimeb)/h2Table-2d0*(fb-fa)/h3Table
  inp_coolrates_table = 10d0**(fa+alpha*yy+beta*yy2+gamma*yy3)
  if( present(retPrime) )                                                &
       retPrime = inp_coolrates_table / T                                &
                                       * (alpha+2d0*beta*yy+3d0*gamma*yy2)
END FUNCTION inp_coolrates_table

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
FUNCTION compCoolrate(T, ne, nHI, nHII, nHeI, nHeII, nHeIII, aexp,       &
                       dcooldT, rt_OTSA)

! Compute cooling rate in a cell, using interpolation from cooling rate
! tables
! T        => Cell emperature [K]
! xion     => Cell ionization fraction
! nH       => Hydrogen number density [cm-3]  
! aexp     => Cosmic expansion
! dcooldT <=  Temperature derivative of the rate
! dcooldx <=  Ionized fraction derivative of the rate
! returns:  Resulting cooling rate [erg s-1 cm-3]  
!-------------------------------------------------------------------------
  implicit none  
  real(dp)::T, ne, nHI, nHII, nHeI, nHeII, nHeIII, aexp
  real(dp)::compCoolrate, dcooldT
  logical::RT_OTSA!-------------------------------------------------------
  real(dp),save::ci_HI,ci_HeI,ci_HeII
  real(dp),save::ci_HI_prime,ci_HeI_prime,ci_HeII_prime
  real(dp),save::ce_HI,ce_HeI,ce_HeII
  real(dp),save::ce_HI_prime,ce_HeI_prime,ce_HeII_prime
  real(dp),save::r_HII,r_HeII,r_HeIII
  real(dp),save::r_HII_prime,r_HeII_prime,r_HeIII_prime
  real(dp),save::bre, bre_prime, com, com_prime, die, die_prime
!-------------------------------------------------------------------------
  ! Coll. Ionization Cooling
  ci_HI   = inp_coolrates_table(tbl_cr_ci_HI, T, ci_HI_prime)            &
          * ne * nHI
  ci_HeI  = inp_coolrates_table(tbl_cr_ci_HeI, T, ci_HeI_prime)          &
          * ne * nHeI
  ci_HeII = inp_coolrates_table(tbl_cr_ci_HeII, T, ci_HeII_prime)        &
          * ne * nHeII
  ci_HI_prime    = ci_HI_prime   * ne * nHI
  ci_HeI_prime   = ci_HeI_prime  * ne * nHeI
  ci_HeII_prime  = ci_HeII_prime * ne * nHeII

  ! Collisional excitation cooling
  ce_HI   = inp_coolrates_table(tbl_cr_ce_HI, T, ce_HI_prime)            &
          * ne * nHI
  ce_HeI  = inp_coolrates_table(tbl_cr_ce_HeI, T, ce_HeI_prime)          &
          * ne *   nHeI
  ce_HeII = inp_coolrates_table(tbl_cr_ce_HeII, T, ce_HeII_prime)        &
          * ne *  nHeII
  ce_HI_prime    = ce_HI_prime   * ne * nHI
  ce_HeI_prime   = ce_HeI_prime  * ne * nHeI
  ce_HeII_prime  = ce_HeII_prime * ne * nHeII

  ! Recombination Cooling
  r_HII   = inp_coolrates_table(tbl_cr_r_HII, T, r_HII_prime)            &
          * ne * nHII
  r_HeII  = inp_coolrates_table(tbl_cr_r_HeII, T, r_HeII_prime)          &
          * ne * nHeII
  r_HeIII = inp_coolrates_table(tbl_cr_r_HeIII, T, r_HeIII_prime)        &
          * ne * nHeIII
  r_HII_prime   = r_HII_prime   * ne * nHII
  r_HeII_prime  = r_HeII_prime  * ne * nHeII
  r_HeIII_prime = r_HeIII_prime * ne * nHeIII

  ! Bremsstrahlung
  bre  = inp_coolrates_table(tbl_cr_bre, T, bre_prime)                   &
       * ne * (nHII + nHeII + 4. * nHeIII)
  bre_prime = bre_prime  * ne * (nHII + nHeII + 4. * nHeIII)

  ! Compton Cooling
  com       = inp_coolrates_table(tbl_cr_com, T, com_prime) * ne    
  com_prime = com_prime * ne

  ! Dielectronic recombination cooling
  die = inp_coolrates_table(tbl_cr_die, T, die_prime) * ne * nHeII
  die_prime = die_prime * ne * nHeII
 
  ! Overall Cooling
  compCoolrate  = ci_HI    + r_HII    + ce_HI    + com                   &
                + ci_HeI   + r_HeII   + ce_HeI   + die                   &
                + ci_HeII  + r_HeIII  + ce_HeII  + bre

  dCooldT       = ci_HI_prime   + r_HII_prime   + ce_HI_prime            &
                + ci_HeI_prime  + r_HeII_prime  + ce_HeI_prime           &
                + ci_HeII_prime + r_HeIII_prime + ce_HeII_prime          &
                + bre_prime     + com_prime     + die_prime
  

END FUNCTION compCoolrate


END MODULE coolrates_module
