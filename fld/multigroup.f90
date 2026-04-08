!> \file
!! Contains module coeff_xi, functions artheta4(), xi(),
!! deriv_artheta4(), deriv_xi() and cal_Teg(), Pr_nu(),
!! Qr_nu , av(), BPlanck() and Div_BPlanck() and
!! subroutines create_groups(), tabulate_art4(), cal_hr()
!! and read_omegas()
!<

!###########################################################
!###########################################################
!###########################################################
!###########################################################

module coeff_xi
  use amr_parameters, only : dp
  real(dp), parameter :: C0 = -6.4939394022668291494e+00_dp
  real(dp), parameter :: C1 = -7.2459729120796762965e+00_dp
  real(dp), parameter :: C2 = -4.0425479966955722625e+00_dp
  real(dp), parameter :: C3 = -1.1702323487246738515e+00_dp
  real(dp), parameter :: C4 = -0.1738347185981371462e+00_dp
  real(dp), parameter :: C5 = -0.7276729479633296956e-02_dp
  real(dp), parameter :: C6 =  6.4939394022668291494e+00_dp
  real(dp), parameter :: C7 =  1.1158054400000000000e+00_dp
  real(dp), parameter :: limhigh = 1.0e+02_dp, limlow = 2.0e-05_dp
end module coeff_xi

!###########################################################
!###########################################################
!###########################################################
!###########################################################

!  Function RADIATION_SOURCE
!
!> Computes radiation source
!<
function radiation_source(T,igrp)
  use const
  use radiation_parameters, only : ngrp,grey_rad_transfer,Tray_min,stellar_photon
  use constants, only:aR
  implicit none

  real(dp), intent(in) :: T
  integer , intent(in) :: igrp
  real(dp)             :: radiation_source,artheta4

  if(grey_rad_transfer)then
     radiation_source = aR*T**4
  else
     if(T.gt.zero) then
        radiation_source = artheta4(T,igrp)
     else
        radiation_source = zero
     end if
  end if

  if(stellar_photon .and. igrp==1)radiation_source=aR*Tray_min**4 
  
  return

end function radiation_source

!###########################################################
!###########################################################
!###########################################################
!###########################################################

!  Function DERIV_RADIATION_SOURCE
!
!> Computes radiation source derivative
!<
function deriv_radiation_source(T,igrp)
  use const
  use radiation_parameters, only : ngrp,grey_rad_transfer,Tray_min,stellar_photon
  use constants, only:aR
  implicit none

  real(dp), intent(in) :: T
  integer , intent(in) :: igrp
  real(dp)              :: deriv_radiation_source,deriv_artheta4

  if(grey_rad_transfer)then
     deriv_radiation_source = four*aR*T**3
  else
     if(T.gt.zero) then
        deriv_radiation_source = deriv_artheta4(T,igrp)
     else
        deriv_radiation_source = zero
     end if
  end if

  if(stellar_photon .and. igrp==1)deriv_radiation_source=four*aR*Tray_min**3

  return

end function deriv_radiation_source

!###########################################################
!###########################################################
!###########################################################
!###########################################################

!  Function ARTHETA4
!
!> Computes the energy of a Planck black body distribution
!! inside a given group.
!<
function artheta4(Tray,igrp)
  use const
  use radiation_parameters, only : nu_min_hz,nu_max_hz,eray_min
  use constants           , only : pi,hplanck,kb,clight

  implicit none

  real(dp), intent(in) :: Tray
  integer , intent(in) :: igrp
  real(dp)             :: constant,xmin,xmax,xsimin,xsimax,xi,artheta4,BPlanck

  constant = (eight*pi*kb**4)/(clight*hplanck)**3

  xmin = hplanck*nu_min_hz(igrp)/(kb*Tray)
  xmax = hplanck*nu_max_hz(igrp)/(kb*Tray)

  xsimin = xi(xmin) ; xsimax = xi(xmax)

  if(xsimin==xsimax)then
     artheta4 = max(BPlanck(half*(nu_min_hz(igrp)+nu_max_hz(igrp)),Tray)*(nu_max_hz(igrp)-nu_min_hz(igrp)),eray_min)
  else
     artheta4 = max(constant*(Tray**4)*(xsimax-xsimin),eray_min)
  endif

  return

end function artheta4
!###########################################################
!###########################################################
!###########################################################
!###########################################################

!  Function XI
!
!> Used by artheta4() to compute the energy of a Planckian
!! distribution between two frequencies.
!<
function xi(nu)
  use const
  use coeff_xi
  implicit none

  real(dp),intent(in) :: nu
  real(dp)            :: xi

  if(nu >= limhigh)then
     xi = c6
  elseif(nu <= limlow)then
     xi = zero
  else
     xi = exp(-c7*nu) * ( c0 + c1*nu + c2*(nu**2) + c3*(nu**3) + c4*(nu**4) + c5*(nu**5) ) + c6
  endif

  if(xi < zero)then
     write(*,*)'negative xi!',xi,nu
     read(*,*)
  endif

end function xi

!###########################################################
!###########################################################
!###########################################################
!###########################################################

!  Function DERIV_ARTHETA4
!
!> Derivative of the artheta4() function which is used
!! to compute the radiative energy source term.
!<
function deriv_artheta4(Tray,igrp)

  use const
  use radiation_parameters, only : nu_min_hz,nu_max_hz,eray_min,deray_min
  use constants           , only : pi,kb,clight,hplanck

  implicit none

  real(dp), intent(in) :: Tray
  integer , intent(in) :: igrp
  real(dp)             :: xi,deriv_xi,deriv_artheta4,nu,dnu,Div_BPlanck
  real(dp)             :: constant,xmin,xmax,xsimin,xsimax,v1,v2,v3

  constant = (eight*pi*kb**4)/(clight*hplanck)**3

  xmin = hplanck*nu_min_hz(igrp)/(kb*Tray)
  xmax = hplanck*nu_max_hz(igrp)/(kb*Tray)

  xsimin = xi(xmin) ; xsimax = xi(xmax)

  if(xsimin==xsimax)then
     deriv_artheta4 = max(Div_BPlanck(half*(nu_min_hz(igrp)+nu_max_hz(igrp)),Tray)*(nu_max_hz(igrp)-nu_min_hz(igrp)),deray_min)
  else
     v1 = four*constant*(Tray**3)*(xsimax-xsimin)
     v2 = deriv_xi(xmin,Tray)
     v3 = deriv_xi(xmax,Tray)
     deriv_artheta4 = v1 + constant*(Tray**4)*(v3-v2)
     deriv_artheta4 = max(deriv_artheta4,deray_min)
  endif

  return

end function deriv_artheta4

!###########################################################
!###########################################################
!###########################################################
!###########################################################

!  Function DERIV_XI
!
!> Derivative of the xi() function which is used in the
!! artheta4() computations.
!<
function deriv_xi(x,T)

  use amr_parameters, only : dp
  use coeff_xi
  use const

  implicit none

  real(dp),intent(in) :: x,T
  real(dp)            :: deriv_xi,p0,p1,p2

  if(x >= limhigh)then
     deriv_xi = zero
  elseif(x <= limlow)then
     deriv_xi = zero
  else
     p0 = exp(-c7*x)
     p1 = c7*(x/T) * p0 * (c0 + c1*x + c2*(x**2) + c3*(x**3) + c4*(x**4) + c5*(x**5))
     p2 = -p0*( c1*x/T + 2.0d0*c2*(x**2)/T + 3.0d0*c3*(x**3)/T + 4.0d0*c4*(x**4)/T + 5.0d0*c5*(x**5)/T )
     deriv_xi = p1 + p2
  endif

end function deriv_xi

!###########################################################
!###########################################################
!###########################################################
!###########################################################

!  Function CAL_TEG
!
!> Computes the temperature of a black body which would
!! provide the same energy as Eg inside a given group igrp.
!<
function cal_Teg(Eg,igrp)

  use amr_parameters      , only : dp
  use radiation_parameters, only : inverse_art4_E,Ninv_art4,dEr_inv_art4

  implicit none

  real(dp), intent(in) :: Eg
  integer , intent(in) :: igrp
  integer              :: iEg
  real(dp)             :: cal_Teg,m,x1,x2,y1,y2,lEg

  lEg = log10(Eg)

!!$  print*, 'begin cal_Teg *****'
!!$  write(*,'(4x,15(ES17.8))') Eg,lEg,inverse_art4_E(1,igrp,1),inverse_art4_E(1,igrp,Ninv_art4)

  if(lEg < inverse_art4_E(1,igrp,1))then
     cal_Teg = inverse_art4_E(2,igrp,1)
  elseif(lEg >= inverse_art4_E(1,igrp,Ninv_art4))then
     cal_Teg = inverse_art4_E(2,igrp,Ninv_art4)
  else
     iEg = int((lEg-inverse_art4_E(1,igrp,1))/dEr_inv_art4(igrp)) + 1
     x1 = inverse_art4_E(1,igrp,iEg  )
     x2 = inverse_art4_E(1,igrp,iEg+1)
     y1 = inverse_art4_E(2,igrp,iEg  )
     y2 = inverse_art4_E(2,igrp,iEg+1)
     ! compute gradient
     m = (y2-y1)/(x2-x1)
     cal_Teg = (m*(lEg-x1))+y1
!!$     write(*,'(8x,i4,15(ES17.8))') iEg,x1,lEg,x2,y1,cal_Teg,y2
  endif

  cal_Teg = 10.0_dp**(cal_Teg)

!!$  write(*,'(8x,ES17.8)') cal_Teg

end function cal_Teg

!###########################################################
!###########################################################
!###########################################################
!###########################################################

!  Function CAL_TEG_SLOW
!
!> Computes the temperature of a black body which would
!! provide the same energy as Eg inside a given group igrp.
!! Slow version of cal_Teg but needed when increment in
!! radiative energy is not constant.
!<
function cal_Teg_slow(Eg,igrp)

  use amr_parameters      , only : dp
  use radiation_parameters, only : inverse_art4_T,Ninv_art4,ngrp

  implicit none

  real(dp), intent(in) :: Eg
  integer , intent(in) :: igrp
  integer              :: i
  real(dp)             :: cal_Teg_slow,m,x1,x2,y1,y2
  logical              :: q

  q            = .true.
  cal_Teg_slow =  inverse_art4_T(ngrp+1,1)
  m            =  0.0_dp

  i = 1

!!$  write(*,*) 'Teg_slow: Eg',Eg,cal_Teg_slow
!!$  write(*,*) 'inv',inverse_art4_T(igrp,i)

  ! search through grid and locate temp
  do while(q)
     if((Eg < inverse_art4_T(igrp,i)).or.(i == Ninv_art4))then
        q = .false.
     elseif(Eg == inverse_art4_T(igrp,i)) then
        cal_Teg_slow = inverse_art4_T(ngrp+1,i)
        q = .false.
     else
        if(Eg < inverse_art4_T(igrp,i+1)) then
!!$           write(*,*) 'Teg_slow: i,cal_Teg_slow',i,igrp
!!$           write(*,*) 'cal_Teg_slow',cal_Teg_slow
!!$           write(*,*) 'Eg',Eg
!!$           write(*,*) 'inverse',inverse_art4_T(igrp,i)
           ! first order linear interpolation
           x1 = inverse_art4_T(igrp  ,i  )
           x2 = inverse_art4_T(igrp  ,i+1)
           y1 = inverse_art4_T(ngrp+1,i  )
           y2 = inverse_art4_T(ngrp+1,i+1)
           ! compute gradient
           m = (y2-y1)/(x2-x1)
           cal_Teg_slow = (m*(Eg-x1))+y1
           q = .false.
           !write(*,*) 'Teg_slow: i,cal_Teg_slow',i,cal_Teg_slow,Eg,inverse_art4_T(igrp,i)
        else
           i = i + 1
        endif
     endif
     !write(*,*) 'Teg_slow: i,cal_Teg_slow',i,cal_Teg_slow,Eg,inverse_art4_T(igrp,i)
  enddo

  !write(*,*) 'Teg_slow: i,cal_Teg_slow',i,cal_Teg_slow,Eg,inverse_art4_T(igrp,i)

end function cal_Teg_slow

!###########################################################
!###########################################################
!###########################################################
!###########################################################

!  Subroutine CREATE_GROUPS
!
!> Creates the frequency groups used in a multigroup
!! radiative hydrodynamics simulation.
!!
!! - If the number of groups is 1, then the group lower and
!!   upper limits are set to numin and numax which are read
!!   in from the namelist.
!!
!! - If ngrp > 1, the routine either reads a list of
!!   group boundaries from a file named 'groups.dat' or
!!   splits the groups automatically.
!!
!! - The groups are either read in Hz or eV units by
!!   specifying the keyword 'freqs_in_Hz' to .true. or
!!   .false.
!!
!! - If the groups are split automatically, this can be done
!!   either logarithmically or lineraly by setting the
!!   keyword 'split_groups_log' to .true. or .false. The
!!   splitting is performed between the numin and numax
!!   frequencies read in from the namelist.
!!
!! - An end group holding frequencies from the last group
!!   boundary to +infinity can be created by setting the
!!   keyword 'extra_end_group' to .true. in the namelist.
!!   This will NOT create an additional group, it will
!!   split the range (numax-numin) into (ngrp-1)
!!   groups and make the last group from numax to
!!   +infinity.
!<
subroutine create_groups

  use amr_parameters, only : dp
  use amr_commons,only: myid
  use hydro_parameters,only: ngrp
  use radiation_parameters

  implicit none

  real(dp)    :: fstep,hz2ev,ev2hz
  integer :: igrp,nstep

  if(ngrp == 1)then
     if(freqs_in_Hz)then
        nu_min_hz(1) = numin
        nu_max_hz(1) = numax
     else
        nu_min_ev(1) = numin
        nu_max_ev(1) = numax
     endif
  else

     ! read groups from file?
     if(read_groups)then

        open(19,file='groups.dat',status='old')
        do igrp = 1,ngrp
           if(freqs_in_Hz)then
              read(19,*)nu_min_hz(igrp),nu_max_hz(igrp)
           else
              read(19,*)nu_min_ev(igrp),nu_max_ev(igrp)
           endif
        enddo

     else ! by default split groups evenly

        if(extra_end_group)then
           nstep = ngrp-1
        else
           nstep = ngrp
        endif

        if(freqs_in_Hz)then

           ! Frequencies in Hz
           if(split_groups_log)then
              fstep = (log10(numax)-log10(numin))/float(nstep)
              nu_min_hz(1) = numin
              nu_max_hz(1) = 10.0d0**(log10(numin) + fstep)
              do igrp = 2,nstep
                 nu_min_hz(igrp) = nu_max_hz(igrp-1)
                 nu_max_hz(igrp) = 10.0d0**(log10(nu_min_hz(igrp)) + fstep)
              enddo
           else
              fstep = (numax-numin)/float(nstep)
              nu_min_hz(1) = numin
              nu_max_hz(1) = fstep + numin
              do igrp = 2,nstep
                 nu_min_hz(igrp) = nu_max_hz(igrp-1)
                 nu_max_hz(igrp) = nu_min_hz(igrp) + fstep
              enddo
           endif

           if(extra_end_group)then
              nu_min_hz(ngrp) = nu_max_hz(ngrp-1)
              nu_max_hz(ngrp) = frequency_upperlimit
           endif

        else

           ! Frequencies in eV
           if(split_groups_log)then
              fstep = (log10(numax)-log10(numin))/float(nstep)
              nu_min_ev(1) = numin
              nu_max_ev(1) = 10.0d0**(log10(numin) + fstep)
              do igrp = 2,nstep
                 nu_min_ev(igrp) = nu_max_ev(igrp-1)
                 nu_max_ev(igrp) = 10.0d0**(log10(nu_min_ev(igrp)) + fstep)
              enddo
           else
              fstep = (numax-numin)/float(nstep)
              nu_min_ev(1) = numin
              nu_max_ev(1) = fstep + numin
              do igrp = 2,nstep
                 nu_min_ev(igrp) = nu_max_ev(igrp-1)
                 nu_max_ev(igrp) = nu_min_ev(igrp) + fstep
              enddo
           endif

           if(extra_end_group)then
              nu_min_ev(ngrp) = nu_max_ev(ngrp-1)
              nu_max_ev(ngrp) = frequency_upperlimit
           endif

        endif
     endif
  endif

  do igrp = 1,ngrp
     if(freqs_in_Hz)then
        nu_min_ev(igrp) = hz2ev(nu_min_hz(igrp))  ! convert from Hz to eV
        nu_max_ev(igrp) = hz2ev(nu_max_hz(igrp))  ! convert from Hz to eV
     else
        nu_min_hz(igrp) = ev2hz(nu_min_ev(igrp))  ! convert from eV to Hz
        nu_max_hz(igrp) = ev2hz(nu_max_ev(igrp))  ! convert from eV to Hz
     endif
  enddo

  if(myid==1) then
     write(*,*) ' '
     write(*,*) 'Number of groups: ',ngrp
     write(*,*) ' '
     write(*,*) 'Group frequencies in Hz:'
     do igrp = 1,ngrp
        write(*,998) igrp,nu_min_hz(igrp),nu_max_hz(igrp)
     enddo
     write(*,*) ' '
     write(*,*) 'Group frequencies in eV:'
     do igrp = 1,ngrp
        write(*,998) igrp,nu_min_ev(igrp),nu_max_ev(igrp)
     enddo
     write(*,*) ' '
  endif

  998 format('igrp, numin, numax = ',i4,2(2x,es12.4))

  return

end subroutine create_groups

!###########################################################
!###########################################################
!###########################################################
!###########################################################

!  Function HZ2EV
!
!> Converts frequency from Hz to eV.
!<
function hz2ev(nu)

  use amr_parameters, only : dp
  use constants     , only : hplanck,eV2erg

  implicit none
  
  real(dp), intent(in) :: nu
  real(dp)             :: hz2eV
  
  hz2eV = nu*hplanck/eV2erg

end function hz2ev

!###########################################################
!###########################################################
!###########################################################
!###########################################################

!  Function EV2HZ
!
!> Converts frequency from eV to Hz.
!<
function ev2hz(nu)

  use amr_parameters, only : dp
  use constants     , only : hplanck,eV2erg

  implicit none
  
  real(dp), intent(in) :: nu
  real(dp)             :: ev2hz
  
  ev2hz = nu*eV2erg/hplanck

end function ev2hz

!###########################################################
!###########################################################
!###########################################################
!###########################################################

!  Subroutine TABULATE_ART4
!
!> Tabulates the artheta4 function to find group interface
!! values for Eray in the comoving frame matter/radiation
!! coupling terms.
!<
subroutine tabulate_art4

  use amr_parameters      , only : dp
  use radiation_parameters, only : ngrp,Ninv_art4,inverse_art4_T,inverse_art4_E,dEr_inv_art4

  implicit none

  integer  :: i,igrp,j
  real(dp) :: T,T1,T2,artheta4,cal_Teg_slow,dTinv_art4

  allocate(inverse_art4_T(ngrp+1,Ninv_art4),inverse_art4_E(2,ngrp,Ninv_art4),dEr_inv_art4(ngrp))

  T1 = log10(1.0_dp) ; T2 = log10(1.0e+07_dp)

  dTinv_art4 = (T2 - T1)/real(Ninv_art4-1,dp)

  ! First pass: log-regular in temperature
  do i = 1,Ninv_art4
     T = real(i-1,dp)*dTinv_art4 + T1
     do igrp = 1,ngrp
        inverse_art4_T(igrp,i) = log10(artheta4(10.0_dp**(T),igrp))
     enddo
     inverse_art4_T(ngrp+1,i) = T
  enddo

  ! Second pass: re-sample curves with regular dEr
  do igrp = 1,ngrp
     dEr_inv_art4(igrp) = (inverse_art4_T(igrp,Ninv_art4) - inverse_art4_T(igrp,1))/real(Ninv_art4-1,dp)
     do i = 1,Ninv_art4
        inverse_art4_E(1,igrp,i) = real(i-1,dp)*dEr_inv_art4(igrp)+inverse_art4_T(igrp,1)
     enddo
  enddo
  ! Warning: do NOT merge this loop with the previous one (does not work with ifort -O3)
  do igrp=1,ngrp
     do i = 1,Ninv_art4
        inverse_art4_E(2,igrp,i) = cal_Teg_slow(inverse_art4_E(1,igrp,i),igrp)
     enddo
  enddo

  !deallocate(inverse_art4_T)

  return

end subroutine tabulate_art4

!###########################################################
!###########################################################
!###########################################################
!###########################################################

!  Function BPLANCK
!
!> Computes the Planck Black Body distribution function.
!<
function BPlanck(nu,T)

  use amr_parameters, only : dp
  use constants     , only : pi,kb,clight,hplanck
  use coeff_xi      , only : limhigh
  use const

  implicit none

  real(dp), intent(in) :: nu,T
  real(dp)             :: BPlanck

  if((hplanck*nu/(kb*T)) > limhigh)then
     BPlanck = (eight*pi*hplanck*nu**3)/clight**3 * exp(-hplanck*nu/(kb*T))
  else
     BPlanck = (eight*pi*hplanck*nu**3)/clight**3 / ( exp(hplanck*nu/(kb*T)) - one )
  endif

end function BPlanck

!###########################################################
!###########################################################
!###########################################################
!###########################################################

!  Function DIV_BPLANCK
!
!> Computes the derivative of the Planck Black Body
!! distribution function.
!<
function Div_BPlanck(nu,T)

  use amr_parameters, only : dp
  use cooling_module, only : kb,hplanck
  use coeff_xi  , only : limhigh
  use const

  implicit none

  real(dp), intent(in) :: nu,T
  real(dp)             :: Div_BPlanck,x,BPlanck,y,ee

  x = hplanck*nu/(kb*T)
  if(x > limhigh)then
     Div_BPlanck = BPlanck(nu,T) * (x/T)
  else
     ee = exp(x)
     y = ee / (ee - one)
     Div_BPlanck = BPlanck(nu,T) * (x/T) * y
  endif

end function Div_BPlanck

!###########################################################
!###########################################################
!###########################################################
!###########################################################

!  Function PR_NU
!
!> Computes the group interface values used to calculate the
!! Doppler terms Pr_nu in the coupling of the radiation to
!! the fluid in the comoving frame.
!<
function Pr_nu(Pr,Er,igrp,nu,scale_d,scale_v)

  use amr_parameters, only : dp
  use radiation_parameters,only:nu_min_hz,nu_max_hz

  implicit none

  integer , intent(in) :: igrp
  real(dp), intent(in) :: Pr,Er,nu,scale_d,scale_v
  real(dp)             :: Te,Dedd,cal_Teg,BPlanck,Pr_nu

  Dedd  = Pr/Er
  Te    = cal_Teg(Er*scale_d*scale_v**2,igrp)
  Pr_nu = Dedd*BPlanck(nu,Te)/(scale_d*scale_v**2)

end function Pr_nu

!###########################################################
!###########################################################
!###########################################################
!###########################################################

!  Function QR_NU
!
!> Computes the group interface values used to calculate the
!! Doppler terms Qr_nu in the coupling of the radiation to
!! the fluid in the comoving frame.
!<
function Qr_nu(Qr,Er,igrp,nu,scale_d,scale_v)

  use amr_parameters, only : dp

  implicit none

  integer , intent(in) :: igrp
  real(dp), intent(in) :: Qr,Er,nu,scale_d,scale_v
  real(dp)             :: Qr_nu,Hedd,Te,cal_Teg,BPlanck

  Hedd  = Qr/Er
  Te    = cal_Teg(Er*scale_d*scale_v**2,igrp)
  Qr_nu = Hedd*BPlanck(nu,Te)/(scale_d*scale_v**2)

end function Qr_nu

!###########################################################
!###########################################################
!###########################################################
!###########################################################

!  Function AV
!
!> Returns minmod average.
!<
function av(a,b)

  use amr_parameters, only : dp
  use const

  implicit none

  real(dp) :: av,a,b

  ! minmod averaging
  if(abs(a) > abs(b))then
     av = b
  else 
     av = a
  endif
  if(a*b < zero)then
     av = zero
  endif

end function av

!###########################################################
!###########################################################
!###########################################################
!###########################################################

!  Subroutine CAL_HR
!
!> Computes the H 3D tensor used to calculate the Doppler
!! terms Q in the coupling of the radiation to the fluid in
!! the comoving frame.
!<
subroutine cal_Hr(E,F,Hr)

  use amr_parameters      , only : ndim,dp
  use constants           , only : clight
  use radiation_parameters, only : irad_trans_model,irad_trans_model_p1,irad_trans_model_m1
  use const

  implicit none

  real(dp)                  :: E
  real(dp),dimension(1:3  ) :: F
  real(dp),dimension(3,3,3) :: Hr

  integer :: i
  real(dp)    :: normef,fx,fy,fz,w1,w2,omega

  Hr = zero

  select case(irad_trans_model)

  case(irad_trans_model_p1) ! 'P1'

     return

  case(irad_trans_model_m1) ! 'M1'

     fx = zero ; fy = zero ; fz = zero

                   fx = F(1)/(clight*E) 
     if(ndim.gt.1) fy = F(2)/(clight*E)
     if(ndim.gt.2) fz = F(3)/(clight*E)

     normef = 0.
     do i = 1,ndim
        normef = normef + (F(i)/clight/E)**2
     enddo
     normef = sqrt(normef)

     if(normef > one) normef = one

     w1 = omega(1,normef)
     w2 = omega(2,normef)

     Hr(1,1,1) = w1*fx*three + w2*fx*fx*fx ; Hr(2,1,1) = w1*fy       + w2*fy*fx*fx ; Hr(3,1,1) = w1*fz       + w2*fz*fx*fx 
     Hr(1,2,1) = w1*fy       + w2*fx*fy*fx ; Hr(2,2,1) = w1*fx       + w2*fy*fy*fx ; Hr(3,2,1) =               w2*fz*fy*fx
     Hr(1,3,1) = w1*fz       + w2*fx*fz*fx ; Hr(2,3,1) =               w2*fy*fz*fx ; Hr(3,3,1) = w1*fx       + w2*fz*fz*fx

     Hr(1,1,2) = w1*fy       + w2*fx*fx*fy ; Hr(2,1,2) = w1*fx       + w2*fy*fx*fy ; Hr(3,1,2) =               w2*fz*fx*fy
     Hr(1,2,2) = w1*fx       + w2*fx*fy*fy ; Hr(2,2,2) = w1*fy*three + w2*fy*fy*fy ; Hr(3,2,2) = w1*fz       + w2*fz*fy*fy
     Hr(1,3,2) =               w2*fx*fz*fy ; Hr(2,3,2) = w1*fz       + w2*fy*fz*fy ; Hr(3,3,2) = w1*fy       + w2*fz*fz*fy

     Hr(1,1,3) = w1*fz       + w2*fx*fx*fz ; Hr(2,1,3) =               w2*fy*fx*fz ; Hr(3,1,3) = w1*fx       + w2*fz*fx*fz
     Hr(1,2,3) =               w2*fx*fy*fz ; Hr(2,2,3) = w1*fz       + w2*fy*fy*fz ; Hr(3,2,3) = w1*fy       + w2*fz*fy*fz
     Hr(1,3,3) = w1*fx       + w2*fx*fz*fz ; Hr(2,3,3) = w1*fy       + w2*fy*fz*fz ; Hr(3,3,3) = w1*fz*three + w2*fz*fz*fz

  end select

  return

end subroutine cal_Hr

!###########################################################
!###########################################################
!###########################################################
!###########################################################

!  Function OMEGA
!
!> Performs a Hermitian cubic spline interpolation between
!! tabulated points on the omega curves to find a specific
!! value of omega_1 or omega_2.
!<
function omega(p,fr)

  use amr_parameters, only : dp
  use radiation_parameters,only:Nomega,f_array,w1_array,dw1_array,w2_array,dw2_array
  use const

  implicit none

  real(dp)   , intent(in) :: fr
  integer, intent(in) :: p
  real(dp)                :: x0,x1,y0,y1,m0,m1,omega,fnew
  real(dp)                :: h00,h10,h01,h11,h,t
  integer             :: i0,i1

  if((p.ne.1).and.(p.ne.2)) stop 'invalid p'

  fnew = fr

  if(fnew > one) fnew = one

  if(fnew==one)then
     omega = real(p,dp) - one
  else
     ! locate between which two array elements fr falls
     i0 = int(fnew*real(Nomega,dp))
     ! define points at x+1
     i1 = i0 + 1
     ! get x and y values for array points
     x0 = f_array(i0)
     x1 = f_array(i1)
     if(p==1)then
        y0 =  w1_array(i0)
        y1 =  w1_array(i1)
        m0 = Dw1_array(i0)
        m1 = Dw1_array(i1)
     else
        y0 =  w2_array(i0)
        y1 =  w2_array(i1)
        m0 = Dw2_array(i0)
        m1 = Dw2_array(i1)
     endif

     h = x1 - x0
     t = (fnew-x0)/h

     h00 = (one+two*t)*(one-t)*(one-t)
     h10 = t*(one-t)*(one-t)
     h01 = t*t*(three-two*t)
     h11 = t*t*(t-one)

     omega = h00*y0 + h10*h*m0 + h01*y1 + h11*h*m1

  endif

end function omega

!###########################################################
!###########################################################
!###########################################################
!###########################################################

!  Subroutine READ_OMEGAS
!
!> Creates a tabulated array of omega values in order to
!! compute the H tensor.
!<
subroutine read_omegas

  use amr_parameters, only : dp
  use radiation_parameters, only : f_array,w1_array,dw1_array,w2_array,dw2_array

  implicit none

     f_array(  0) = 0.000000000000e+00 ; w1_array(  0) = 0.200000000000E+00 ; dw1_array(  0) =  0.000000000000E+00 ; w2_array(  0) = 0.482142857143E+00 ; dw2_array(  0) = 0.000000000000E+00
     f_array(  1) = 0.100000000000e-01 ; w1_array(  1) = 0.199990356721E+00 ; dw1_array(  1) = -0.192874186607E-02 ; w2_array(  1) = 0.482163952013E+00 ; dw2_array(  1) = 0.421920237875E-02
     f_array(  2) = 0.200000000000e-01 ; w1_array(  2) = 0.199961421821E+00 ; dw1_array(  2) = -0.385849623214E-02 ; w2_array(  2) = 0.482227250059E+00 ; dw2_array(  2) = 0.844109215605E-02
     f_array(  3) = 0.300000000000e-01 ; w1_array(  3) = 0.199913180114E+00 ; dw1_array(  3) = -0.579027559821E-02 ; w2_array(  3) = 0.482332791593E+00 ; dw2_array(  3) = 0.126683567304E-01
     f_array(  4) = 0.400000000000e-01 ; w1_array(  4) = 0.199845606286E+00 ; dw1_array(  4) = -0.772509246429E-02 ; w2_array(  4) = 0.482480643799E+00 ; dw2_array(  4) = 0.169036835004E-01
     f_array(  5) = 0.500000000000e-01 ; w1_array(  5) = 0.199758664900E+00 ; dw1_array(  5) = -0.966415666507E-02 ; w2_array(  5) = 0.482670900736E+00 ; dw2_array(  5) = 0.211503373143E-01
     f_array(  6) = 0.600000000000e-01 ; w1_array(  6) = 0.199652309346E+00 ; dw1_array(  6) = -0.116079935741E-01 ; w2_array(  6) = 0.482903686402E+00 ; dw2_array(  6) = 0.254095801429E-01
     f_array(  7) = 0.700000000000e-01 ; w1_array(  7) = 0.199526484436E+00 ; dw1_array(  7) = -0.135581199196E-01 ; w2_array(  7) = 0.483179147144E+00 ; dw2_array(  7) = 0.296855749031E-01
     f_array(  8) = 0.800000000000e-01 ; w1_array(  8) = 0.199381123250E+00 ; dw1_array(  8) = -0.155154276734E-01 ; w2_array(  8) = 0.483497460882E+00 ; dw2_array(  8) = 0.339806563647E-01
     f_array(  9) = 0.900000000000e-01 ; w1_array(  9) = 0.199216148691E+00 ; dw1_array(  9) = -0.174809756129E-01 ; w2_array(  9) = 0.483858832569E+00 ; dw2_array(  9) = 0.382976474651E-01
     f_array( 10) = 0.100000000000e+00 ; w1_array( 10) = 0.199031473009E+00 ; dw1_array( 10) = -0.194558358817E-01 ; w2_array( 10) = 0.484263495574E+00 ; dw2_array( 10) = 0.426394103351E-01
     f_array( 11) = 0.110000000000e+00 ; w1_array( 11) = 0.198826997657E+00 ; dw1_array( 11) = -0.214410958009E-01 ; w2_array( 11) = 0.484711712108E+00 ; dw2_array( 11) = 0.470088516363E-01
     f_array( 12) = 0.120000000000e+00 ; w1_array( 12) = 0.198602613132E+00 ; dw1_array( 12) = -0.234378597313E-01 ; w2_array( 12) = 0.485203773696E+00 ; dw2_array( 12) = 0.514089280577E-01
     f_array( 13) = 0.130000000000e+00 ; w1_array( 13) = 0.198358198794E+00 ; dw1_array( 13) = -0.254472509947E-01 ; w2_array( 13) = 0.485740001699E+00 ; dw2_array( 13) = 0.558426519930E-01
     f_array( 14) = 0.140000000000e+00 ; w1_array( 14) = 0.198093622668E+00 ; dw1_array( 14) = -0.274704138590E-01 ; w2_array( 14) = 0.486320747909E+00 ; dw2_array( 14) = 0.603130974217E-01
     f_array( 15) = 0.150000000000e+00 ; w1_array( 15) = 0.197808741222E+00 ; dw1_array( 15) = -0.295085155969E-01 ; w2_array( 15) = 0.486946395188E+00 ; dw2_array( 15) = 0.648234060165E-01
     f_array( 16) = 0.160000000000e+00 ; w1_array( 16) = 0.197503399126E+00 ; dw1_array( 16) = -0.315627486249E-01 ; w2_array( 16) = 0.487617358180E+00 ; dw2_array( 16) = 0.693767935024E-01
     f_array( 17) = 0.170000000000e+00 ; w1_array( 17) = 0.197177428994E+00 ; dw1_array( 17) = -0.336343327304E-01 ; w2_array( 17) = 0.488334084086E+00 ; dw2_array( 17) = 0.739765562942E-01
     f_array( 18) = 0.180000000000e+00 ; w1_array( 18) = 0.196830651093E+00 ; dw1_array( 18) = -0.357245173976E-01 ; w2_array( 18) = 0.489097053501E+00 ; dw2_array( 18) = 0.786260784397E-01
     f_array( 19) = 0.190000000000e+00 ; w1_array( 19) = 0.196462873039E+00 ; dw1_array( 19) = -0.378345842398E-01 ; w2_array( 19) = 0.489906781330E+00 ; dw2_array( 19) = 0.833288388997E-01
     f_array( 20) = 0.200000000000e+00 ; w1_array( 20) = 0.196073889458E+00 ; dw1_array( 20) = -0.399658495499E-01 ; w2_array( 20) = 0.490763817772E+00 ; dw2_array( 20) = 0.880884191967E-01
     f_array( 21) = 0.210000000000e+00 ; w1_array( 21) = 0.195663481630E+00 ; dw1_array( 21) = -0.421196669783E-01 ; w2_array( 21) = 0.491668749390E+00 ; dw2_array( 21) = 0.929085114678E-01
     f_array( 22) = 0.220000000000e+00 ; w1_array( 22) = 0.195231417102E+00 ; dw1_array( 22) = -0.442974303525E-01 ; w2_array( 22) = 0.492622200254E+00 ; dw2_array( 22) = 0.977929269588E-01
     f_array( 23) = 0.230000000000e+00 ; w1_array( 23) = 0.194777449265E+00 ; dw1_array( 23) = -0.465005766481E-01 ; w2_array( 23) = 0.493624833180E+00 ; dw2_array( 23) = 0.102745605000E+00
     f_array( 24) = 0.240000000000e+00 ; w1_array( 24) = 0.194301316916E+00 ; dw1_array( 24) = -0.487305891273E-01 ; w2_array( 24) = 0.494677351057E+00 ; dw2_array( 24) = 0.107770622510E+00
     f_array( 25) = 0.250000000000e+00 ; w1_array( 25) = 0.193802743772E+00 ; dw1_array( 25) = -0.509890006592E-01 ; w2_array( 25) = 0.495780498273E+00 ; dw2_array( 25) = 0.112872204071E+00
     f_array( 26) = 0.260000000000e+00 ; w1_array( 26) = 0.193281437958E+00 ; dw1_array( 26) = -0.532773972369E-01 ; w2_array( 26) = 0.496935062244E+00 ; dw2_array( 26) = 0.118054732630E+00
     f_array( 27) = 0.270000000000e+00 ; w1_array( 27) = 0.192737091462E+00 ; dw1_array( 27) = -0.555974217108E-01 ; w2_array( 27) = 0.498141875054E+00 ; dw2_array( 27) = 0.123322760888E+00
     f_array( 28) = 0.280000000000e+00 ; w1_array( 28) = 0.192169379537E+00 ; dw1_array( 28) = -0.579507777554E-01 ; w2_array( 28) = 0.499401815213E+00 ; dw2_array( 28) = 0.128681023427E+00
     f_array( 29) = 0.290000000000e+00 ; w1_array( 29) = 0.191577960084E+00 ; dw1_array( 29) = -0.603392340922E-01 ; w2_array( 29) = 0.500715809539E+00 ; dw2_array( 29) = 0.134134449649E+00
     f_array( 30) = 0.300000000000e+00 ; w1_array( 30) = 0.190962472967E+00 ; dw1_array( 30) = -0.627646289896E-01 ; w2_array( 30) = 0.502084835176E+00 ; dw2_array( 30) = 0.139688177610E+00
     f_array( 31) = 0.310000000000e+00 ; w1_array( 31) = 0.190322539304E+00 ; dw1_array( 31) = -0.652288750665E-01 ; w2_array( 31) = 0.503509921753E+00 ; dw2_array( 31) = 0.145347568807E+00
     f_array( 32) = 0.320000000000e+00 ; w1_array( 32) = 0.189657760692E+00 ; dw1_array( 32) = -0.677339644253E-01 ; w2_array( 32) = 0.504992153695E+00 ; dw2_array( 32) = 0.151118224032E+00
     f_array( 33) = 0.330000000000e+00 ; w1_array( 33) = 0.188967718389E+00 ; dw1_array( 33) = -0.702819741452E-01 ; w2_array( 33) = 0.506532672702E+00 ; dw2_array( 33) = 0.157006000376E+00
     f_array( 34) = 0.340000000000e+00 ; w1_array( 34) = 0.188251972429E+00 ; dw1_array( 34) = -0.728750721684E-01 ; w2_array( 34) = 0.508132680403E+00 ; dw2_array( 34) = 0.163017029488E+00
     f_array( 35) = 0.350000000000e+00 ; w1_array( 35) = 0.187510060691E+00 ; dw1_array( 35) = -0.755155236152E-01 ; w2_array( 35) = 0.509793441195E+00 ; dw2_array( 35) = 0.169157737223E+00
     f_array( 36) = 0.360000000000e+00 ; w1_array( 36) = 0.186741497885E+00 ; dw1_array( 36) = -0.782056975680E-01 ; w2_array( 36) = 0.511516285294E+00 ; dw2_array( 36) = 0.175434864794E+00
     f_array( 37) = 0.370000000000e+00 ; w1_array( 37) = 0.185945774484E+00 ; dw1_array( 37) = -0.809480743679E-01 ; w2_array( 37) = 0.513302611997E+00 ; dw2_array( 37) = 0.181855491583E+00
     f_array( 38) = 0.380000000000e+00 ; w1_array( 38) = 0.185122355565E+00 ; dw1_array( 38) = -0.837452534730E-01 ; w2_array( 38) = 0.515153893188E+00 ; dw2_array( 38) = 0.188427059764E+00
     f_array( 39) = 0.390000000000e+00 ; w1_array( 39) = 0.184270679583E+00 ; dw1_array( 39) = -0.865999619315E-01 ; w2_array( 39) = 0.517071677096E+00 ; dw2_array( 39) = 0.195157400923E+00
     f_array( 40) = 0.400000000000e+00 ; w1_array( 40) = 0.183390157045E+00 ; dw1_array( 40) = -0.895150635295E-01 ; w2_array( 40) = 0.519057592331E+00 ; dw2_array( 40) = 0.202054764868E+00
     f_array( 41) = 0.410000000000e+00 ; w1_array( 41) = 0.182480169098E+00 ; dw1_array( 41) = -0.924935686797E-01 ; w2_array( 41) = 0.521113352221E+00 ; dw2_array( 41) = 0.209127850847E+00
     f_array( 42) = 0.420000000000e+00 ; w1_array( 42) = 0.181540066006E+00 ; dw1_array( 42) = -0.955386451247E-01 ; w2_array( 42) = 0.523240759474E+00 ; dw2_array( 42) = 0.216385841427E+00
     f_array( 43) = 0.430000000000e+00 ; w1_array( 43) = 0.180569165520E+00 ; dw1_array( 43) = -0.986536295356E-01 ; w2_array( 43) = 0.525441711190E+00 ; dw2_array( 43) = 0.223838439298E+00
     f_array( 44) = 0.440000000000e+00 ; w1_array( 44) = 0.179566751131E+00 ; dw1_array( 44) = -0.101842040099E+00 ; w2_array( 44) = 0.527718204262E+00 ; dw2_array( 44) = 0.231495907309E+00
     f_array( 45) = 0.450000000000e+00 ; w1_array( 45) = 0.178532070182E+00 ; dw1_array( 45) = -0.105107590192E+00 ; w2_array( 45) = 0.530072341193E+00 ; dw2_array( 45) = 0.239369112083E+00
     f_array( 46) = 0.460000000000e+00 ; w1_array( 46) = 0.177464331845E+00 ; dw1_array( 46) = -0.108454203261E+00 ; w2_array( 46) = 0.532506336365E+00 ; dw2_array( 46) = 0.247469571592E+00
     f_array( 47) = 0.470000000000e+00 ; w1_array( 47) = 0.176362704942E+00 ; dw1_array( 47) = -0.111886029035E+00 ; w2_array( 47) = 0.535022522809E+00 ; dw2_array( 47) = 0.255809507113E+00
     f_array( 48) = 0.480000000000e+00 ; w1_array( 48) = 0.175226315593E+00 ; dw1_array( 48) = -0.115407461206E+00 ; w2_array( 48) = 0.537623359523E+00 ; dw2_array( 48) = 0.264401900068E+00
     f_array( 49) = 0.490000000000e+00 ; w1_array( 49) = 0.174054244681E+00 ; dw1_array( 49) = -0.119023156752E+00 ; w2_array( 49) = 0.540311439372E+00 ; dw2_array( 49) = 0.273260554271E+00
     f_array( 50) = 0.500000000000e+00 ; w1_array( 50) = 0.172845525117E+00 ; dw1_array( 50) = -0.122738057074E+00 ; w2_array( 50) = 0.543089497652E+00 ; dw2_array( 50) = 0.282400164218E+00
     f_array( 51) = 0.510000000000e+00 ; w1_array( 51) = 0.171599138880E+00 ; dw1_array( 51) = -0.126557411158E+00 ; w2_array( 51) = 0.545960421369E+00 ; dw2_array( 51) = 0.291836390107E+00
     f_array( 52) = 0.520000000000e+00 ; w1_array( 52) = 0.170314013817E+00 ; dw1_array( 52) = -0.130486800991E+00 ; w2_array( 52) = 0.548927259295E+00 ; dw2_array( 52) = 0.301585940376E+00
     f_array( 53) = 0.530000000000e+00 ; w1_array( 53) = 0.168989020175E+00 ; dw1_array( 53) = -0.134532169495E+00 ; w2_array( 53) = 0.551993232903E+00 ; dw2_array( 53) = 0.311666662670E+00
     f_array( 54) = 0.540000000000e+00 ; w1_array( 54) = 0.167622966842E+00 ; dw1_array( 54) = -0.138699851272E+00 ; w2_array( 54) = 0.555161748251E+00 ; dw2_array( 54) = 0.322097644240E+00
     f_array( 55) = 0.550000000000e+00 ; w1_array( 55) = 0.166214597260E+00 ; dw1_array( 55) = -0.142996606518E+00 ; w2_array( 55) = 0.558436408921E+00 ; dw2_array( 55) = 0.332899322958E+00
     f_array( 56) = 0.560000000000e+00 ; w1_array( 56) = 0.164762584990E+00 ; dw1_array( 56) = -0.147429658466E+00 ; w2_array( 56) = 0.561821030135E+00 ; dw2_array( 56) = 0.344093610271E+00
     f_array( 57) = 0.570000000000e+00 ; w1_array( 57) = 0.163265528873E+00 ; dw1_array( 57) = -0.152006734828E+00 ; w2_array( 57) = 0.565319654156E+00 ; dw2_array( 57) = 0.355704027626E+00
     f_array( 58) = 0.580000000000e+00 ; w1_array( 58) = 0.161721947763E+00 ; dw1_array( 58) = -0.156736113724E+00 ; w2_array( 58) = 0.568936567132E+00 ; dw2_array( 58) = 0.367755858121E+00
     f_array( 59) = 0.590000000000e+00 ; w1_array( 59) = 0.160130274772E+00 ; dw1_array( 59) = -0.161626674693E+00 ; w2_array( 59) = 0.572676317545E+00 ; dw2_array( 59) = 0.380276315401E+00
     f_array( 60) = 0.600000000000e+00 ; w1_array( 60) = 0.158488850976E+00 ; dw1_array( 60) = -0.166687955452E+00 ; w2_array( 60) = 0.576543736441E+00 ; dw2_array( 60) = 0.393294732145E+00
     f_array( 61) = 0.610000000000e+00 ; w1_array( 61) = 0.156795918523E+00 ; dw1_array( 61) = -0.171930215181E+00 ; w2_array( 61) = 0.580543959649E+00 ; dw2_array( 61) = 0.406842770837E+00
     f_array( 62) = 0.620000000000e+00 ; w1_array( 62) = 0.155049613072E+00 ; dw1_array( 62) = -0.177364505238E+00 ; w2_array( 62) = 0.584682452242E+00 ; dw2_array( 62) = 0.420954659961E+00
     f_array( 63) = 0.630000000000e+00 ; w1_array( 63) = 0.153247955482E+00 ; dw1_array( 63) = -0.183002748329E+00 ; w2_array( 63) = 0.588965035496E+00 ; dw2_array( 63) = 0.435667459291E+00
     f_array( 64) = 0.640000000000e+00 ; w1_array( 64) = 0.151388842667E+00 ; dw1_array( 64) = -0.188857827376E+00 ; w2_array( 64) = 0.593397916664E+00 ; dw2_array( 64) = 0.451021358533E+00
     f_array( 65) = 0.650000000000e+00 ; w1_array( 65) = 0.149470037497E+00 ; dw1_array( 65) = -0.194943685479E+00 ; w2_array( 65) = 0.597987721926E+00 ; dw2_array( 65) = 0.467060014337E+00
     f_array( 66) = 0.660000000000e+00 ; w1_array( 66) = 0.147489157651E+00 ; dw1_array( 66) = -0.201275438648E+00 ; w2_array( 66) = 0.602741532927E+00 ; dw2_array( 66) = 0.483830931566E+00
     f_array( 67) = 0.670000000000e+00 ; w1_array( 67) = 0.145443663261E+00 ; dw1_array( 67) = -0.207869503270E+00 ; w2_array( 67) = 0.607666927370E+00 ; dw2_array( 67) = 0.501385895766E+00
     f_array( 68) = 0.680000000000e+00 ; w1_array( 68) = 0.143330843199E+00 ; dw1_array( 68) = -0.214743740607E+00 ; w2_array( 68) = 0.612772024233E+00 ; dw2_array( 68) = 0.519781465078E+00
     f_array( 69) = 0.690000000000e+00 ; w1_array( 69) = 0.141147799830E+00 ; dw1_array( 69) = -0.221917621087E+00 ; w2_array( 69) = 0.618065534232E+00 ; dw2_array( 69) = 0.539079531379E+00
     f_array( 70) = 0.700000000000e+00 ; w1_array( 70) = 0.138891432003E+00 ; dw1_array( 70) = -0.229412411632E+00 ; w2_array( 70) = 0.623556816301E+00 ; dw2_array( 70) = 0.559347962331E+00
     f_array( 71) = 0.710000000000e+00 ; w1_array( 71) = 0.136558416035E+00 ; dw1_array( 71) = -0.237251389940E+00 ; w2_array( 71) = 0.629255940939E+00 ; dw2_array( 71) = 0.580661338366E+00
     f_array( 72) = 0.720000000000e+00 ; w1_array( 72) = 0.134145184409E+00 ; dw1_array( 72) = -0.245460090386E+00 ; w2_array( 72) = 0.635173761485E+00 ; dw2_array( 72) = 0.603101801480E+00
     f_array( 73) = 0.730000000000e+00 ; w1_array( 73) = 0.131647901825E+00 ; dw1_array( 73) = -0.254066587191E+00 ; w2_array( 73) = 0.641321994515E+00 ; dw2_array( 73) = 0.626760036270E+00
     f_array( 74) = 0.740000000000e+00 ; w1_array( 74) = 0.129062438200E+00 ; dw1_array( 74) = -0.263101821699E+00 ; w2_array( 74) = 0.647713310813E+00 ; dw2_array( 74) = 0.651736408089E+00
     f_array( 75) = 0.750000000000e+00 ; w1_array( 75) = 0.126384338154E+00 ; dw1_array( 75) = -0.272599982107E+00 ; w2_array( 75) = 0.654361438630E+00 ; dw2_array( 75) = 0.678142288700E+00
     f_array( 76) = 0.760000000000e+00 ; w1_array( 76) = 0.123608786384E+00 ; dw1_array( 76) = -0.282598945829E+00 ; w2_array( 76) = 0.661281281298E+00 ; dw2_array( 76) = 0.706101606882E+00
     f_array( 77) = 0.770000000000e+00 ; w1_array( 77) = 0.120730568255E+00 ; dw1_array( 77) = -0.293140797121E+00 ; w2_array( 77) = 0.668489051650E+00 ; dw2_array( 77) = 0.735752670260E+00
     f_array( 78) = 0.780000000000e+00 ; w1_array( 78) = 0.117744024773E+00 ; dw1_array( 78) = -0.304272435550E+00 ; w2_array( 78) = 0.676002426259E+00 ; dw2_array( 78) = 0.767250316117E+00
     f_array( 79) = 0.790000000000e+00 ; w1_array( 79) = 0.114643000941E+00 ; dw1_array( 79) = -0.316046294831E+00 ; w2_array( 79) = 0.683840723117E+00 ; dw2_array( 79) = 0.800768463596E+00
     f_array( 80) = 0.800000000000e+00 ; w1_array( 80) = 0.111420786277E+00 ; dw1_array( 80) = -0.328521196582E+00 ; w2_array( 80) = 0.692025107214E+00 ; dw2_array( 80) = 0.836503158934E+00
     f_array( 81) = 0.810000000000e+00 ; w1_array( 81) = 0.108070045991E+00 ; dw1_array( 81) = -0.341763370139E+00 ; w2_array( 81) = 0.700578829512E+00 ; dw2_array( 81) = 0.874676230561E+00
     f_array( 82) = 0.820000000000e+00 ; w1_array( 82) = 0.104582740977E+00 ; dw1_array( 82) = -0.355847678308E+00 ; w2_array( 82) = 0.709527506120E+00 ; dw2_array( 82) = 0.915539704460E+00
     f_array( 83) = 0.830000000000e+00 ; w1_array( 83) = 0.100950034301E+00 ; dw1_array( 83) = -0.370859100587E+00 ; w2_array( 83) = 0.718899446213E+00 ; dw2_array( 83) = 0.959381175256E+00
     f_array( 84) = 0.840000000000e+00 ; w1_array( 84) = 0.971621813085E-01 ; dw1_array( 84) = -0.386894541136E+00 ; w2_array( 84) = 0.728726039481E+00 ; dw2_array( 84) = 0.100653038994E+01
     f_array( 85) = 0.850000000000e+00 ; w1_array( 85) = 0.932083996661E-01 ; dw1_array( 85) = -0.404065050437E+00 ; w2_array( 85) = 0.739042216844E+00 ; dw2_array( 85) = 0.105736738594E+01
     f_array( 86) = 0.860000000000e+00 ; w1_array( 86) = 0.890767146385E-01 ; dw1_array( 86) = -0.422498579680E+00 ; w2_array( 86) = 0.749887002174E+00 ; dw2_array( 86) = 0.111233264427E+01
     f_array( 87) = 0.870000000000e+00 ; w1_array( 87) = 0.847537735118E-01 ; dw1_array( 87) = -0.442343429609E+00 ; w2_array( 87) = 0.761304178149E+00 ; dw2_array( 87) = 0.117193988857E+01
     f_array( 88) = 0.880000000000e+00 ; w1_array( 88) = 0.802246211610E-01 ; dw1_array( 88) = -0.463772617233E+00 ; w2_array( 88) = 0.773343096843E+00 ; dw2_array( 88) = 0.123679240851E+01
     f_array( 89) = 0.890000000000e+00 ; w1_array( 89) = 0.754724261056E-01 ; dw1_array( 89) = -0.486989474897E+00 ; w2_array( 89) = 0.786059676142E+00 ; dw2_array( 89) = 0.130760415588E+01
     f_array( 90) = 0.900000000000e+00 ; w1_array( 90) = 0.704781426104E-01 ; dw1_array( 90) = -0.512234934468E+00 ; w2_array( 90) = 0.799517638208E+00 ; dw2_array( 90) = 0.138522742777E+01
     f_array( 91) = 0.910000000000e+00 ; w1_array( 91) = 0.652200888646E-01 ; dw1_array( 91) = -0.539797165468E+00 ; w2_array( 91) = 0.813790068442E+00 ; dw2_array( 91) = 0.147068984712E+01
     f_array( 92) = 0.920000000000e+00 ; w1_array( 92) = 0.596734129822E-01 ; dw1_array( 92) = -0.570024586660E+00 ; w2_array( 92) = 0.828961407241E+00 ; dw2_array( 92) = 0.156524482125E+01
     f_array( 93) = 0.930000000000e+00 ; w1_array( 93) = 0.538094056736E-01 ; dw1_array( 93) = -0.603343865088E+00 ; w2_array( 93) = 0.845130040041E+00 ; dw2_array( 93) = 0.167044218495E+01
     f_array( 94) = 0.940000000000e+00 ; w1_array( 94) = 0.475945975265E-01 ; dw1_array( 94) = -0.640285581678E+00 ; w2_array( 94) = 0.862411738759E+00 ; dw2_array( 94) = 0.178823032434E+01
     f_array( 95) = 0.950000000000e+00 ; w1_array( 95) = 0.409895430072E-01 ; dw1_array( 95) = -0.681522287652E+00 ; w2_array( 95) = 0.880944360071E+00 ; dw2_array( 95) = 0.192111002367E+01
     f_array( 96) = 0.960000000000e+00 ; w1_array( 96) = 0.339471274188E-01 ; dw1_array( 96) = -0.727927975827E+00 ; w2_array( 96) = 0.900894491000E+00 ; dw2_array( 96) = 0.207237939241E+01
     f_array( 97) = 0.970000000000e+00 ; w1_array( 97) = 0.264100996931E-01 ; dw1_array( 97) = -0.780678218773E+00 ; w2_array( 97) = 0.922467320156E+00 ; dw2_array( 97) = 0.224655548053E+01
     f_array( 98) = 0.980000000000e+00 ; w1_array( 98) = 0.183072236399E-01 ; dw1_array( 98) = -0.841439339546E+00 ; w2_array( 98) = 0.945922409205E+00 ; dw2_array( 98) = 0.245019224521E+01
     f_array( 99) = 0.990000000000e+00 ; w1_array( 99) = 0.954653858989E-02 ; dw1_array( 99) = -0.912806423112E+00 ; w2_array( 99) = 0.971602190644E+00 ; dw2_array( 99) = 0.269383781672E+01
     f_array(100) = 1.000000000000e+00 ; w1_array(100) = 0.000000000000E+00 ; dw1_array(100) = -0.100000000000E+01 ; w2_array(100) = 0.100000000000E+01 ; dw2_array(100) = 0.300000000000E+01

  return

end subroutine read_omegas
