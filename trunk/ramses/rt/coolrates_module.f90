MODULE coolrates_module
  ! Add an init routine that tabulates cooling rates and reaction rates
  ! Add a routine to interpolate the table and return those rates
  ! Table.T
  ! Use actual calculations if temperature is outside table boundaries.

  use amr_parameters,only:dp
  implicit none

  !private   ! default

  !public init_coolrates_table, get_AlphaA_HII

  ! Default cooling rates table parameters
  integer,parameter     :: nbinTK = 101
  real(dp),parameter    :: TKmin  = 1d-2
  real(dp),parameter    :: TKmax  = 1d+9
  real(dp)              :: dlogTinv ! Inverse of the bin space (in K)
  real(dp)              :: hTable, h2Table, h3Table   ! Interpol constants
  
  type coolrates_table
     real(kind=8),dimension(:)    ,pointer::TK
     real(kind=8),dimension(:)    ,pointer::alphaA_HII
     real(kind=8),dimension(:)    ,pointer::alphaA_HII_prime
  end type coolrates_table

  type(coolrates_table)::table, table_mpi
  
CONTAINS

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE init_coolrates_table

! Initialise the cooling rates table.
!-------------------------------------------------------------------------
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  real(kind=8):: TK
  integer:: myid, ncpu, ierr
  integer:: iTK
  logical,save:: first=.true.
!-------------------------------------------------------------------------
#ifndef WITHOUTMPI
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,ncpu,ierr)
#endif
#ifdef WITHOUTMPI
  myid=0
  ncpu=1
#endif

  ! Allocate the tables --------------------------------------------------
  if(.not.first)then
     deallocate(table%TK)
     deallocate(table%alphaA_HII)
  else
     first=.false.
  endif

  allocate(table%TK(1:nbinTK))
  allocate(table%alphaA_HII(1:nbinTK))
  allocate(table%alphaA_HII_prime(1:nbinTK))

  ! Initialise the table temperature indexes -----------------------------
  do iTK=1, nbinTK
     table%TK(iTK) = log10(TKmin)    &
                    + (dble(iTK)-1d0) / (dble(nbinTK)-1d0) &
                    * (log10(TKmax)-log10(TKmin))
  end do
  dlogTinv = dble(nbinTK-1)/(table%TK(nbinTK)-table%TK(1)) ! (space)^-1
  hTable = 1d0/dlogTinv                     !
  h2Table = hTable*hTable                   ! Constants for table
  h3Table = h2Table*hTable                  ! interpolation
  
  ! Initialise the table values ------------------------------------------
  table%alphaA_HII = 0d0
  table%alphaA_HII_prime = 0d0
  do iTK = myid+1, nbinTK, ncpu ! Loop over TK and assign values
  !do iTK=1, nbinTK
     TK = 10d0**table%TK(iTK)
     table%alphaA_HII(iTK) = log10(MAX(comp_AlphaA_HII(TK),1d-100))
     table%alphaA_HII_prime(iTK) = log10(MAX(comp_dAlphaA_dT_HII(TK),1d-100))
  end do ! end TK loop
  
  ! Distribute the complete table between cpus ---------------------------
#ifndef WITHOUTMPI
  allocate(table_mpi%alphaA_HII(1:nbinTK))
  allocate(table_mpi%alphaA_HII_prime(1:nbinTK))

  call MPI_ALLREDUCE(table%alphaA_HII,table_mpi%alphaA_HII               &
       ,nbinTK,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  table%alphaA_HII = table_mpi%alphaA_HII
  call MPI_ALLREDUCE(table%alphaA_HII_prime,table_mpi%alphaA_HII_prime   &
       ,nbinTK,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  table%alphaA_HII_prime = table_mpi%alphaA_HII_prime

  deallocate(table_mpi%alphaA_HII)
  deallocate(table_mpi%alphaA_HII_prime)
#endif

  if(myid==0) print*,'Coolrates table initialised '
  do iTK = 1,nbinTK
     if(myid==0) print*,10d0**table%TK(iTK)     &
          , 10d0**table%alphaA_HII(iTK)         &
          , 10d0**table%alphaA_HII_prime(iTK)
     !if(myid==1) print*,table%alphaA_HII
  end do

  !stop
  
END SUBROUTINE init_coolrates_table

!PRIVATE XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION comp_AlphaA_HII(T)

! Returns case A rec. coefficient [cm3 s-1] for HII (Hui&Gnedin'97)
! T           => Temperature [K]
!-------------------------------------------------------------------------
  implicit none  
  real(dp),intent(in)::T
  real(dp)::comp_AlphaA_HII, lambda
!-------------------------------------------------------------------------
  lambda=  315614./T                                ! 2.d0 * 157807.d0 / T
  comp_AlphaA_HII =  1.269d-13 * lambda**1.503 &
                     / ( ( 1.d0+(lambda/0.522)**0.47 )**1.923 )
END FUNCTION comp_AlphaA_HII

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
FUNCTION get_AlphaA_HII(T, useLastIndex, prime)

! Returns TABULATED case A rec. coefficient [cm3 s-1] for HII
! T            => Temperature [K]
! useLastIndex => Use interpolation indexes from last coolrates call
! prime        <= Temperature derivative of returned rate
!-------------------------------------------------------------------------
  implicit none
  real(dp),intent(in)::T
  logical, optional, intent (in) :: useLastIndex
  real(dp),optional::prime
  real(dp)::get_AlphaA_HII
!-------------------------------------------------------------------------
  get_AlphaA_HII = inp_coolrates_table(table%alphaA_HII                  &
                          ,table%alphaA_HII_prime, T, useLastIndex, prime)
END FUNCTION get_AlphaA_HII

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
FUNCTION inp_coolrates_table(rates, ratesPrime, T, useLastIndex, retPrime)
! Returns TABULATED rate value from given table
! rates        => Rates table to interpolate
! ratesPrime   => Prime rates table to interpolate
! T            => Temperature [K]
! useLastIndex => Re-use table index and weights from last call to
!                 this function.
! retPrime     <= Temperature derivative of rate at T (optional).  
!-------------------------------------------------------------------------
  implicit none
  real(dp),dimension(nbinTK)::rates, ratesPrime
  real(dp),intent(in)::T
  logical, optional, intent (in) :: useLastIndex
  real(dp),optional::retPrime

  real(dp)::inp_coolrates_table

  integer,save:: iT = 1
  real(dp),save:: facT, yy, yy2, yy3, fa, fb, fprimea, fprimeb
  real(dp),save:: alpha, beta, gamma
!-------------------------------------------------------------------------
  if (.not. useLastIndex) then
     ! Log of T, snapped to table:
     facT = MIN( MAX( log10(T), table%TK(1) ), table%TK(nbinTK) )
     ! Lower closest index in table:
     iT = MIN(MAX(int((facT-table%TK(1))*dlogTinv)+1, 1), nbinTK-1) 
     yy=facT-table%TK(iT)  ! Dist., in log(T), from T to lower table index
     yy2=yy*yy             ! That distance squared
     yy3=yy2*yy            ! ...and cubed
  endif
  
  fa = rates(iT)                !
  fb = rates(iT+1)              !  Values at neighbouring table
  fprimea=ratesPrime(iT)       !  indexes
  fprimeb=ratesPrime(iT+1)     !

  ! Spline interpolation:
  alpha = fprimea
  beta = 3d0*(fb-fa)/h2Table-(2d0*fprimea+fprimeb)/hTable
  gamma = (fprimea+fprimeb)/h2Table-2d0*(fb-fa)/h3Table
  inp_coolrates_table = 10d0**(fa+alpha*yy+beta*yy2+gamma*yy3)
  if( present(retPrime) ) &
       retPrime = inp_coolrates_table / T &
                * (alpha+2d0*beta*yy+3d0*gamma*yy2)

END FUNCTION inp_coolrates_table

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION comp_dAlphaA_dT_HII(T)

! Returns Temperature derivative of the case A rec. coeff. of HII 
! T           => Temperature [K]
!-------------------------------------------------------------------------
  implicit none  
  real(dp),intent(in)::T
  real(dp)::comp_dAlphaA_dT_HII, lambda, f
!-------------------------------------------------------------------------
  lambda=  315614./T
  f= 1.d0+(lambda/0.522)**0.47
  comp_dAlphaA_dT_HII = 1.269d-13 * lambda**1.503 / f**1.923             &
                        / T * ( 0.90381*(f-1.)/f - 1.503 )
END FUNCTION comp_dAlphaA_dT_HII

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION comp_AlphaA_HeII(T)

! Returns case A rec. coefficient [cm3 s-1] for HeII (Hui&Gnedin'97)
! T           => Temperature [K]
!-------------------------------------------------------------------------
  implicit none  
  real(dp),intent(in)::T
  real(dp)::comp_AlphaA_HeII, lambda
!-------------------------------------------------------------------------
  lambda=  570670./T
  comp_AlphaA_HeII =  3.d-14 * lambda**0.654
END FUNCTION comp_AlphaA_HeII

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION comp_AlphaA_HeIII(T)

! Returns case A rec. coefficient [cm3 s-1] for HeIII (Hui&Gnedin'97)
! T           => Temperature [K]
!-------------------------------------------------------------------------
  implicit none  
  real(dp),intent(in)::T
  real(dp)::comp_AlphaA_HeIII, lambda
!-------------------------------------------------------------------------
  lambda=  1263030./T
  comp_AlphaA_HeIII =  2.538d-13 * lambda**1.503                         &
                     / ( ( 1.d0+(lambda/0.522)**0.47 )**1.923 )
END FUNCTION comp_AlphaA_HeIII

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION comp_AlphaB_HII(T)

! Returns case B rec. coefficient [cm3 s-1] for HII (Hui&Gnedin'97)
! T           => Temperature [K]
!-------------------------------------------------------------------------
  implicit none  
  real(dp),intent(in)::T
  real(dp)::comp_AlphaB_HII, lambda
!-------------------------------------------------------------------------
  lambda = 315614./T
  comp_AlphaB_HII =  &
       2.753d-14 * lambda**1.5 / ( (1.d0+(lambda/2.74)**0.407)**2.242 )
END FUNCTION comp_AlphaB_HII

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION comp_AlphaB_HeII(T)

! Returns case B rec. coefficient [cm3 s-1] for HeII (Hui&Gnedin'97)
! T           => Temperature [K]
!-------------------------------------------------------------------------
  implicit none  
  real(dp),intent(in)::T
  real(dp)::comp_AlphaB_HeII, lambda
!-------------------------------------------------------------------------
  lambda = 570670./T
  comp_AlphaB_HeII = 1.26d-14 * lambda**0.75
END FUNCTION comp_AlphaB_HeII

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION comp_AlphaB_HeIII(T)

! Returns case B rec. coefficient [cm3 s-1] for HeIII (Hui&Gnedin'97)
! T           => Temperature [K]
!-------------------------------------------------------------------------
  implicit none  
  real(dp),intent(in)::T
  real(dp)::comp_AlphaB_HeIII, lambda
!-------------------------------------------------------------------------
  lambda=  1263030./T
  comp_AlphaB_HeIII =  5.506d-14 * lambda**1.5                           &
                     / ( ( 1.d0+(lambda/2.74)**0.407 )**2.242 )
END FUNCTION comp_AlphaB_HeIII

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION comp_dAlphaB_dT_HII(T)

! Returns temperature derivative of the case B recombination rate
! Gnedin 1997
! T           => Temperature [K]
!-------------------------------------------------------------------------
  implicit none  
  real(dp),intent(in)::T
  real(dp)::comp_dAlphaB_dT_HII, lambda, f
!-------------------------------------------------------------------------
  lambda = 315614./T
  f= 1.d0+(lambda/2.74)**0.407
  comp_dAlphaB_dT_HII = 2.753d-14 * lambda**1.5 / f**2.242               &
                  / T * ( 0.912494*(f-1.)/f - 1.5 )
END FUNCTION comp_dAlphaB_dT_HII

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION comp_Beta_HI(T)

! Returns collisional ionization rate [cm3 s-1] of HI (Maselli&'03)
! T           => Temperature [K]
!-------------------------------------------------------------------------
  implicit none  
  real(dp),intent(in)::T
  real(dp)::comp_Beta_HI, T5
!-------------------------------------------------------------------------
  T5 = T/1.d5
  comp_Beta_HI = &
                5.85d-11 * sqrt(T) / (1.d0+sqrt(T5)) * exp(-157809.1d0/T)
END FUNCTION comp_Beta_HI

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION comp_dBeta_dT_HI(T)

! Returns T derivative of collisional ionization rate [cm3 s-1] of HI 
! (Maselli&'03)
! T           => Temperature [K]
!-------------------------------------------------------------------------
  implicit none  
  real(dp),intent(in)::T
  real(dp)::f, T5, comp_dBeta_dT_HI
!-------------------------------------------------------------------------
  T5 = T/1.d5
  f = 1.+sqrt(T5)
  comp_dBeta_dT_HI = 5.85d-11 * sqrt(T) / f * exp(-157809.1d0/T)         &
                * (0.5/T + 157809.1d0/T**2 - .5/(sqrt(1.d5*T)*f))
END FUNCTION comp_dBeta_dT_HI

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION comp_Beta_HeI(T)

! Returns collisional ionization rate [cm3 s-1] of HeI (Maselli&'03)
! T           => Temperature [K]
!-------------------------------------------------------------------------
  implicit none  
  real(dp),intent(in)::T
  real(dp)::comp_Beta_HeI, T5
!-------------------------------------------------------------------------
  T5 = T/1.d5
  comp_Beta_HeI = &
                2.38d-11 * sqrt(T) / (1.d0+sqrt(T5)) * exp(-285335.4d0/T)
END FUNCTION comp_Beta_HeI

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION comp_Beta_HeII(T)

! Returns collisional ionization rate [cm3 s-1] of HeII (Maselli&'03)
! T           => Temperature [K]
!-------------------------------------------------------------------------
  implicit none  
  real(dp),intent(in)::T
  real(dp)::comp_Beta_HeII, T5
!-------------------------------------------------------------------------
  T5 = T/1.d5
  comp_Beta_HeII = &
                5.68d-12 * sqrt(T) / (1.d0+sqrt(T5)) * exp(-631515.d0/T)
END FUNCTION comp_Beta_HeII

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION comp_collExrate_HI(TK)

! Gives Collisional Excitation rate coefficient for the 1->2
! transition in hydrogen atom energies, in [cm3 s-1], according to
! Callaway, Unnikrishnan & Oza 1987.
! TK      => Temperatue in Kelvin
!-------------------------------------------------------------------------
  real(dp),intent(in)::TK
  real(dp)::comp_collExrate_HI
  real(dp),parameter ::kB    = 1.38062d-16       ! Boltzm.const. [erg K-1]
!-------------------------------------------------------------------------
  comp_collExrate_HI =                                                   &
       2.41d-6/sqrt(TK) * (TK/1.d4)**0.22 * exp(-1.63d-11/(kb*TK))
END FUNCTION comp_collExrate_HI

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
FUNCTION compCoolrate(T, ne, nHI, nHII, nHeI, nHeII, nHeIII, aexp,     &
                       dcooldT, rt_OTSA)

! Compute cooling rate in a cell
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
  real(dp)::compCoolrate, dcooldT!----------------------------------------
  real(dp)::T2, sqrtT, sqrt5T, fT5, sqrtT_fT5, f, Ta
  real(dp)::laHII, laHeII, laHeIII 
  real(dp)::ci_HI, ci_HeI, ci_HeII, dci_HI, dci_HeI, dci_HeII
  real(dp)::r_HII, r_HeII, r_HeIII, dr_HII, dr_HeII, dr_HeIII
  real(dp)::ce_HI, ce_HeI, ce_HeII, dce_HI, dce_HeI, dce_HeII
  real(dp)::bre, dbre, com, dcom, die, ddie
  real(dp),parameter::kb=1.3806d-16        ! Boltzmann constant [ergs K-1]
  logical::RT_OTSA
!-------------------------------------------------------------------------
  T2       = T**2
  sqrtT    = sqrt(T)
  sqrt5T   = sqrt(1.d5*T)
  fT5      = 1.+sqrt(T/1.d5)

  ! Coll. Ionization Cooling from Cen 1992 (via Maselli et al 2003)
  sqrtT_fT5 = sqrtT/fT5
  ci_HI   = 1.27d-21 * sqrtT_fT5 * exp(-157809.1/T)           * ne *   nHI
  ci_HeI  = 9.38d-22 * sqrtT_fT5 * exp(-285335.4/T)           * ne *  nHeI
  ci_HeII = 4.95d-22 * sqrtT_fT5 * exp(-631515. /T)           * ne * nHeII
  f    = 0.5/T - 0.5/(sqrt5T+T)
  dci_HI  = ci_HI *   ( f + 157809.1/T2 )                             !dHI
  dci_HeI = ci_HeI *  ( f + 285335.4/T2 )                            !dHeI
  dci_HeII= ci_HeII * ( f + 631515.0/T2 )                           !dHeII

  ! Recombination Cooling (Hui&Gnedin'97)
  laHII  = 315614./T                                             
  laHeII  = 570670./T                                            
  laHeIII  = 1263030./T
  if(.not. rt_otsa) then ! Case A
     f       = 1.d0+(laHII/0.541)**0.502                         
     r_HII   = 1.778d-29 * laHII**1.965 / f**2.697 * T      * ne *   nHII
     dr_HII  = r_HII/T * ( - 0.965 + 1.35389*(f-1.)/f )             !dHII
 
     r_HeII  = 3.d-14 * laHeII**0.654 * kb * T              * ne *  nHeII
     dr_HeII = r_HeII / T * 0.346                                  !dHeII

     f       = 1.d0+(laHeIII/0.541)**0.502                      
     r_HeIII = 14.224d-29 * laHeIII**1.965 / f**2.697 * T   * ne * nHeIII
     dr_HeIII= r_HeIII/T * ( - 0.965 + 1.35389*(f-1.)/f )         !dHeIII     
  else ! Case B
     f      = 1.d0+(laHII/2.25)**0.376                        
     r_HII  = 3.435d-30 * laHII**1.97 / f**3.72 * T         * ne *   nHII
     dr_HII = r_HII/T * ( - 0.97 + 1.39827*(f-1.)/f )               !dHII

     r_HeII  = 1.26d-14 * laHeII**0.75 * kb * T             * ne *  nHeII
     dr_HeII = r_HeII / T * 0.25                                   !dHeII

     f       = 1.d0+(laHeIII/2.25)**0.376                       
     r_HeIII = 27.48d-30 * laHeIII**1.97 / f**3.72 * T      * ne * nHeIII
     dr_HeIII =r_HeIII/T * ( - 0.97 + 1.39827*(f-1.)/f )          !dHeIII
  endif

  ! Collisional excitation cooling from Cen'92
  ce_HI  = 7.5d-19 * exp(-118348./T) / fT5                  * ne *    nHI
  ce_HeI = 9.10d-27 * exp(-13179./T) / fT5 * T**(-0.1687)   * ne *   nHeI
  ce_HeII = 5.54d-17 * exp(-473638./T) / fT5 * T**(-0.397)  * ne *  nHeII

  f = 0.5/(sqrt5T+T)
  dce_HI  = ce_HI * ( 118348./T2 - f)                                !dHI
  dce_HeI  = ce_HeI * ( 13179./T2 - f - 0.1687/T)                   !dHeI
  dce_HeII  = ce_HeII * ( 473638./T2 - f - 0.397/T)                !dHeII

  ! Bremsstrahlung from Osterbrock & Ferland 2006
  bre  = 1.42d-27 * 1.5 * sqrtT       * ne * (nHII + nHeII + 4. * nHeIII)
  dbre = bre * 0.5 / T

  ! Compton Cooling from Haimann et al. 96, via Maselli et al.
  Ta     = 2.727/aexp                          
  com     = 1.017d-37 * Ta**4 * (T-Ta)                               * ne    
  dcom  = com / (T-Ta)   

  ! Dielectronic recombination cooling, from Black 1981
  f = 1.24D-13*T**(-1.5d0)*exp(-470000.d0/T)                 * ne * nHeII
  die = f*(1.D0+0.3d0*exp(-94000.d0/T))
  ddie = die/T2*(564000.-1.5*T) - f*94000./T2
 
  ! Overall Cooling
  compCoolrate  = ci_HI    + r_HII    + ce_HI    + com    &
                + ci_HeI   + r_HeII   + ce_HeI   + die    &
                + ci_HeII  + r_HeIII  + ce_HeII  + bre

  dCooldT       = dci_HI   + dr_HII   + dce_HI   + dcom   &
                + dci_HeI  + dr_HeII  + dce_HeI  + ddie   &
                + dci_HeII + dr_HeIII + dce_HeII + dbre

END FUNCTION compCoolrate


END MODULE coolrates_module
