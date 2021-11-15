!> \file
!! Contains subroutines compute_valp(), bubble_sort()
!! and function interpol_valp()
!<

!###########################################################
!###########################################################
!###########################################################
!###########################################################

!  Subroutine COMPUTE_VALP
!
!> Computes the tabulated eigenvalues for the M1 radiative
!! transfer solver.
!<
subroutine compute_valp
  
#if USE_M_1==1

  use constants, only : c_cgs
  use radiation_parameters, only : valp,n_points
  use const

  implicit none

  integer :: i,j,k
  real(dp) :: normef,theta,eps,pi

  real(dp), dimension(4,4) :: mat
  real(dp), dimension(4) :: WR,WI
  integer,parameter :: LDVL=1,LDVR=1
  real(dp), dimension(LDVL,4) :: VL
  real(dp), dimension(LDVR,4) :: VR
  integer :: info,lda=4
  integer,parameter :: lwork=20
  real(dp),dimension(lwork) :: work

  real(dp)                    :: E
  real(dp),dimension(1:3    ) :: F
  real(dp),dimension(3,3    ) :: Dedd,Dedd_dE
  real(dp),dimension(3,3,1:3) :: Dedd_dF

  pi=acos(-one)

  E=one
  F=zero

  n_points=50

  allocate(valp(0:n_points,0:n_points,0:n_points,4))

  do i=0,n_points
     normef=float(i)/n_points

     do j=0,n_points
        theta=j*pi/n_points
        
        !Warning : F is normalized by c for the cal_Dedd call
        F(1)=E*cos(theta)*normef
        F(2)=E*sin(theta)*normef

        call cal_Dedd(E,F,Dedd,Dedd_dE,Dedd_dF)

        do k=0,n_points
           eps=float(k)/n_points

           mat=zero

           ! Warning the matrix is normalized to get eigenvalues in units of clight
           mat(1,1) = zero
           mat(1,2) = one
           mat(2,1) = Dedd(1,1)+Dedd_dE(1,1  )
           mat(2,2) =           Dedd_dF(1,1,1) *eps

           mat(1,3) = zero
           mat(2,3) =           Dedd_dF(1,1,2) *eps

           mat(3,1) = Dedd(1,2)+Dedd_dE(1,2  )
           mat(3,2) =           Dedd_dF(1,2,1) *eps
           mat(3,3) =           Dedd_dF(1,2,2) *eps

           mat(1,4) = zero
           mat(2,4) =           Dedd_dF(1,1,3) *eps
           mat(3,4) =           Dedd_dF(1,2,3) *eps

           mat(4,1) = Dedd(1,3)+Dedd_dE(1,3  )
           mat(4,2) =           Dedd_dF(1,3,1) *eps
           mat(4,3) =           Dedd_dF(1,3,2) *eps
           mat(4,4) =           Dedd_dF(1,3,3) *eps

           !Warning : dsyev works only for symmetric matrix -> use dgeev instead
           call dgeev('N','N',4,mat,lda,WR,WI,VL,LDVL,VR,LDVR,WORK,LWORK,info)
           
           ! sort the real part of the eigenvalues in ascending order
           !
           ! warning : the eigenvalues may be not real but complex...
           !           in particular it happens when eps->0 and f->1
           call bubble_sort(4,WR)
           valp(i,j,k,1:4)=WR(1:4)
        enddo
     enddo
  enddo

  return

#endif

end subroutine compute_valp

!###########################################################
!###########################################################
!###########################################################
!###########################################################

!  Subroutine BUBBLE_SORT
!
!> Sorts array ARR of length n in ascending order.
!<
subroutine bubble_sort(n,arr)

  use amr_parameters, only : dp

  implicit none

  integer                :: n,i,j
  real(dp)               :: a
  real(dp), dimension(n) :: arr

  do j = 1,n
     do i = 2,n
        if(arr(i) < arr(i-1))then
           a        = arr(i  )
           arr(i  ) = arr(i-1)
           arr(i-1) = a
        endif
     enddo
  enddo

  return

end subroutine bubble_sort

!###########################################################
!###########################################################
!###########################################################
!###########################################################

!  Function INTERPOL_VALP
!
!> Performs interpolation between tabulated point on the
!! eigenvalues curves to find a particular eigenvalue
!! during the simulation.
!<
function interpol_valp(f,t,e,ivalp)

  use radiation_parameters, only : n_points,valp
  use const

  implicit none
  
  integer :: ivalp
  integer :: i_f,i_t,i_e
  integer :: i_f_p,i_t_p,i_e_p
  
  real(dp) :: f,t,e,pi
  real(dp) :: interpol_valp
  real(dp) :: fmin,tmin,emin,df,dt,de,lf,lt,le,dd1,dd2,dd3,de1,de2,de3

  pi=acos(-one)

  fmin = zero
  tmin = zero
  emin = zero

  df = one/n_points
  dt = pi /n_points
  de = one/n_points

  lf = (f -fmin )/df   ! No +0.5 since we calculate them at the cell interfaces and not the centres
  lt = (t -tmin )/dt
  le = (e -emin )/de

  i_f = int(lf) ; i_t = int(lt) ; i_e = int(le)

  dd1 = lf - float(i_f) 
  dd2 = lt - float(i_t)
  dd3 = le - float(i_e)

  de1 = one - dd1
  de2 = one - dd2
  de3 = one - dd3

  ! Boundary conditions: at the boundaries, dd = 0 so value of b is unimportant as long as b < 1
  i_f_p = min(n_points,i_f+1)
  i_t_p = min(n_points,i_t+1)
  i_e_p = min(n_points,i_e+1)

  interpol_valp = zero
  interpol_valp = interpol_valp + de1*de2*de3*valp(i_f  ,i_t  ,i_e  ,ivalp)
  interpol_valp = interpol_valp + dd1*de2*de3*valp(i_f_p,i_t  ,i_e  ,ivalp)
  interpol_valp = interpol_valp + de1*dd2*de3*valp(i_f  ,i_t_p,i_e  ,ivalp)
  interpol_valp = interpol_valp + dd1*dd2*de3*valp(i_f_p,i_t_p,i_e  ,ivalp)
  interpol_valp = interpol_valp + de1*de2*dd3*valp(i_f  ,i_t  ,i_e_p,ivalp)
  interpol_valp = interpol_valp + dd1*de2*dd3*valp(i_f_p,i_t  ,i_e_p,ivalp)
  interpol_valp = interpol_valp + de1*dd2*dd3*valp(i_f  ,i_t_p,i_e_p,ivalp)
  interpol_valp = interpol_valp + dd1*dd2*dd3*valp(i_f_p,i_t_p,i_e_p,ivalp)

end function interpol_valp

!###########################################################
!###########################################################
!###########################################################
!###########################################################

!  Subroutine CAL_DEDD
!
!> Computes the Eddington tensor Dedd and its derivatives
!! with respect to Er and Fr, Dedd_dE and Dedd_dF,
!! respectively.
!!
!! WARNING: The input variable F is normalized by clight,
!! therefore, reduced_flux = F / E
!<
subroutine cal_Dedd(E,F,Dedd,Dedd_dE,Dedd_dF)

  use amr_parameters      , only : ndim
  use constants           , only : c_cgs
  use radiation_parameters, only : irad_trans_model,irad_trans_model_p1,irad_trans_model_m1
  use const

  implicit none

  real(dp)                       :: E
  real(dp),dimension(1:3       ) :: F
  real(dp),dimension(1:3       ) :: unit_vector,reduced_flux
  real(dp),dimension(3,3       ) :: Dedd,Dedd_dE
  real(dp),dimension(3,3,3     ) :: Dedd_dF

  integer :: i
  real(dp)    :: normeF,chi,chiprime
  real(dp)    :: g, h, gprime, hprime
 

  Dedd    = zero
  Dedd_dE = zero
  Dedd_dF = zero

  select case(irad_trans_model)

  case(irad_trans_model_p1) ! 'P1'

     do i = 1,3
        Dedd(i,i) = one/three
     enddo
     return

  case(irad_trans_model_m1) ! 'M1'

     reduced_flux=zero
     reduced_flux(1:ndim)=F(1:ndim)/E

     normef = zero
     do i=1,ndim
        normef = normef + reduced_flux(i)**2
     enddo
     normef = sqrt(normef)

     if(normef.gt.one) normef = one

     chi      = (three+four*normef**2)/(five+two*sqrt(four-three*normef**2))
     chiprime = two*normef/sqrt(four-three*normef**2)

     !other form from Bruno Dubroca
     !chi      = (one+normef**2*(one+normef**2))/three
     !chiprime = two*normef*(one+two*normef**2))/three

     g      = (one-chi)/two
     h      = (three*chi-one)/two
     gprime = -chiprime/two
     hprime = three*chiprime/two

     do i=1,3
        Dedd(i,i)=g
     enddo

     ! recompute the norm for unit_vector
     normef = zero
     do i=1,ndim
        normef = normef + reduced_flux(i)**2
     enddo
     normef = sqrt(normef)

     if(normef.gt.1.e-5) then
        unit_vector = reduced_flux/normef

        ! Eddington tensor
        Dedd(1,1) = g + h*unit_vector(1)*unit_vector(1) ; Dedd(1,2) =     h*unit_vector(1)*unit_vector(2) ; Dedd(1,3) =     h*unit_vector(1)*unit_vector(3)
        Dedd(2,1) =     h*unit_vector(2)*unit_vector(1) ; Dedd(2,2) = g + h*unit_vector(2)*unit_vector(2) ; Dedd(2,3) =     h*unit_vector(2)*unit_vector(3)
        Dedd(3,1) =     h*unit_vector(3)*unit_vector(1) ; Dedd(3,2) =     h*unit_vector(3)*unit_vector(2) ; Dedd(3,3) = g + h*unit_vector(3)*unit_vector(3)

        ! Eddington tensor derivative with respect to Ei * Ei
        Dedd_dE(1,1) = - normef * (gprime + hprime*unit_vector(1)*unit_vector(1))
        Dedd_dE(1,2) = - normef * (         hprime*unit_vector(1)*unit_vector(2))
        Dedd_dE(1,3) = - normef * (         hprime*unit_vector(1)*unit_vector(3))
        Dedd_dE(2,1) = - normef * (         hprime*unit_vector(2)*unit_vector(1))
        Dedd_dE(2,2) = - normef * (gprime + hprime*unit_vector(2)*unit_vector(2))
        Dedd_dE(2,3) = - normef * (         hprime*unit_vector(2)*unit_vector(3))
        Dedd_dE(3,1) = - normef * (         hprime*unit_vector(3)*unit_vector(1))
        Dedd_dE(3,2) = - normef * (         hprime*unit_vector(3)*unit_vector(2))
        Dedd_dE(3,3) = - normef * (gprime + hprime*unit_vector(3)*unit_vector(3))

        ! Eddington tensor derivative with respect to Fx * Ei
        Dedd_dF(1,1,1) = ( gprime*unit_vector(1) + hprime*unit_vector(1)*unit_vector(1)*unit_vector(1) + h*(-two*unit_vector(1)*unit_vector(1)*unit_vector(1) /normef+two*unit_vector(1)/normef) )
        Dedd_dF(1,2,1) = (                         hprime*unit_vector(1)*unit_vector(1)*unit_vector(2) + h*(-two*unit_vector(1)*unit_vector(1)*unit_vector(2) /normef+    unit_vector(2)/normef) )
        Dedd_dF(1,3,1) = (                         hprime*unit_vector(1)*unit_vector(1)*unit_vector(3) + h*(-two*unit_vector(1)*unit_vector(1)*unit_vector(3) /normef+    unit_vector(3)/normef) )

        Dedd_dF(2,1,1) = (                         hprime*unit_vector(1)*unit_vector(2)*unit_vector(1) + h*(-two*unit_vector(1)*unit_vector(2)*unit_vector(1) /normef+    unit_vector(2)/normef) )
        Dedd_dF(2,2,1) = ( gprime*unit_vector(1) + hprime*unit_vector(1)*unit_vector(2)*unit_vector(2) + h*(-two*unit_vector(1)*unit_vector(2)*unit_vector(2) /normef                          ) )
        Dedd_dF(2,3,1) = (                         hprime*unit_vector(1)*unit_vector(2)*unit_vector(3) + h*(-two*unit_vector(1)*unit_vector(2)*unit_vector(3) /normef                          ) )

        Dedd_dF(3,1,1) = (                         hprime*unit_vector(1)*unit_vector(3)*unit_vector(1) + h*(-two*unit_vector(1)*unit_vector(3)*unit_vector(1) /normef+    unit_vector(3)/normef) )
        Dedd_dF(3,2,1) = (                         hprime*unit_vector(1)*unit_vector(3)*unit_vector(2) + h*(-two*unit_vector(1)*unit_vector(3)*unit_vector(2) /normef                          ) )
        Dedd_dF(3,3,1) = ( gprime*unit_vector(1) + hprime*unit_vector(1)*unit_vector(3)*unit_vector(3) + h*(-two*unit_vector(1)*unit_vector(3)*unit_vector(3) /normef                          ) )

        if(ndim.gt.1) then
           ! Eddington tensor derivative with respect to Fy * Ei
           Dedd_dF(1,1,2) = ( gprime*unit_vector(2) + hprime*unit_vector(2)*unit_vector(1)*unit_vector(1) + h*(-two*unit_vector(2)*unit_vector(1)*unit_vector(1) /normef                          ) )
           Dedd_dF(1,2,2) = (                         hprime*unit_vector(2)*unit_vector(1)*unit_vector(2) + h*(-two*unit_vector(2)*unit_vector(1)*unit_vector(2) /normef+    unit_vector(1)/normef) )
           Dedd_dF(1,3,2) = (                         hprime*unit_vector(2)*unit_vector(1)*unit_vector(3) + h*(-two*unit_vector(2)*unit_vector(1)*unit_vector(3) /normef                          ) ) 

           Dedd_dF(2,1,2) = (                         hprime*unit_vector(2)*unit_vector(2)*unit_vector(1) + h*(-two*unit_vector(2)*unit_vector(2)*unit_vector(1) /normef+    unit_vector(1)/normef) )
           Dedd_dF(2,2,2) = ( gprime*unit_vector(2) + hprime*unit_vector(2)*unit_vector(2)*unit_vector(2) + h*(-two*unit_vector(2)*unit_vector(2)*unit_vector(2) /normef+two*unit_vector(2)/normef) )
           Dedd_dF(2,3,2) = (                         hprime*unit_vector(2)*unit_vector(2)*unit_vector(3) + h*(-two*unit_vector(2)*unit_vector(2)*unit_vector(3) /normef+    unit_vector(3)/normef) )

           Dedd_dF(3,1,2) = (                         hprime*unit_vector(2)*unit_vector(3)*unit_vector(1) + h*(-two*unit_vector(2)*unit_vector(3)*unit_vector(1) /normef                          ) )
           Dedd_dF(3,2,2) = (                         hprime*unit_vector(2)*unit_vector(3)*unit_vector(2) + h*(-two*unit_vector(2)*unit_vector(3)*unit_vector(2) /normef+    unit_vector(3)/normef) )
           Dedd_dF(3,3,2) = ( gprime*unit_vector(2) + hprime*unit_vector(2)*unit_vector(3)*unit_vector(3) + h*(-two*unit_vector(2)*unit_vector(3)*unit_vector(3) /normef                          ) )

           if(ndim.gt.2) then
              ! Eddington tensor derivative with respect to Fz * Ei
              Dedd_dF(1,1,3) = ( gprime*unit_vector(3) + hprime*unit_vector(3)*unit_vector(1)*unit_vector(1) + h*(-two*unit_vector(3)*unit_vector(1)*unit_vector(1) /normef                          ) )
              Dedd_dF(1,2,3) = (                         hprime*unit_vector(3)*unit_vector(1)*unit_vector(2) + h*(-two*unit_vector(3)*unit_vector(1)*unit_vector(2) /normef                          ) )
              Dedd_dF(1,3,3) = (                         hprime*unit_vector(3)*unit_vector(1)*unit_vector(3) + h*(-two*unit_vector(3)*unit_vector(1)*unit_vector(3) /normef+    unit_vector(1)/normef) ) 

              Dedd_dF(2,1,3) = (                         hprime*unit_vector(3)*unit_vector(2)*unit_vector(1) + h*(-two*unit_vector(3)*unit_vector(2)*unit_vector(1) /normef                          ) )
              Dedd_dF(2,2,3) = ( gprime*unit_vector(3) + hprime*unit_vector(3)*unit_vector(2)*unit_vector(2) + h*(-two*unit_vector(3)*unit_vector(2)*unit_vector(2) /normef                          ) )
              Dedd_dF(2,3,3) = (                         hprime*unit_vector(3)*unit_vector(2)*unit_vector(3) + h*(-two*unit_vector(3)*unit_vector(2)*unit_vector(3) /normef+    unit_vector(2)/normef) )

              Dedd_dF(3,1,3) = (                         hprime*unit_vector(3)*unit_vector(3)*unit_vector(1) + h*(-two*unit_vector(3)*unit_vector(3)*unit_vector(1) /normef+    unit_vector(1)/normef) )
              Dedd_dF(3,2,3) = (                         hprime*unit_vector(3)*unit_vector(3)*unit_vector(2) + h*(-two*unit_vector(3)*unit_vector(3)*unit_vector(2) /normef+    unit_vector(2)/normef) )
              Dedd_dF(3,3,3) = ( gprime*unit_vector(3) + hprime*unit_vector(3)*unit_vector(3)*unit_vector(3) + h*(-two*unit_vector(3)*unit_vector(3)*unit_vector(3) /normef+two*unit_vector(3)/normef) )

           endif
        endif
     endif

  end select

  return

end subroutine cal_Dedd
