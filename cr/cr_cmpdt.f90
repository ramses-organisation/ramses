!#########################################################
!#########################################################
subroutine cr_cmpdt_a2(a2,rho,crgather,ncell)
  ! Add the CR pressure contribution to the squared gas sound speed,
  ! a2 += sum_g E_cr,g * gamma_g*(gamma_g-1)/rho, in cells dense enough
  ! for the CRs to stay coupled to the gas.
  use amr_parameters, only: dp,nvector
  use hydro_parameters, only: smallr
  use cr_parameters, only: ncr_groups,gamma_cr,cr_smallr_decouple
  implicit none
  integer,intent(in)::ncell
  real(dp),dimension(1:nvector),intent(inout)::a2
  real(dp),dimension(1:nvector),intent(in)::rho
  real(dp),dimension(1:nvector,1:ncr_groups),intent(in)::crgather
  real(dp),dimension(1:nvector),save::cr_cs
  integer::k,icr

  do k = 1, ncell
     cr_cs(k)=0.0d0
  end do
  do icr = 1,ncr_groups
     do k = 1, ncell
        cr_cs(k)=cr_cs(k) + crgather(k,icr) * gamma_cr(icr)*(gamma_cr(icr)-1.0d0)
     end do
  end do
  do k = 1, ncell
     cr_cs(k)=cr_cs(k)/rho(k)
     ! Only consider CR sound speed where rho is not tiny
     if(rho(k) .gt. smallr*cr_smallr_decouple) &
          a2(k) = a2(k) + cr_cs(k)
  end do

end subroutine cr_cmpdt_a2
!#########################################################
!#########################################################
subroutine cr_cmpdt_vamax(uu,ncell)
  ! Track the maximum Alfven speed for the adaptive reduced light speed.
  use amr_parameters, only: dp,nvector
  use hydro_parameters, only: nvar
  use cr_parameters, only: cr_va_max
  implicit none
  integer,intent(in)::ncell
  real(dp),dimension(1:nvector,1:nvar+3),intent(in)::uu
  real(dp)::BNva
  integer::k,idim

  do k = 1, ncell
     BNva=0d0
     do idim=1,3
        BNva = BNva + (0.5d0*(uu(k,5+idim)+uu(k,nvar+idim)))**2
     end do
     cr_va_max=max(cr_va_max, sqrt(BNva/uu(k,1)))
  end do

end subroutine cr_cmpdt_vamax
