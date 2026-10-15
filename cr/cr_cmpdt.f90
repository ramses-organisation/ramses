!#########################################################
!#########################################################
subroutine cr_cmpdt(uu,gg,ecr,dx,dt,ncell)
  ! Gas Courant time step with the CR pressure in the sound speed:
  ! mhd/courant_fine.f90 cmpdt plus a2 += sum_g E_cr,g*gamma_g*(gamma_g-1)/rho
  ! in cells dense enough for the CRs to stay coupled to the gas.
  use amr_parameters
  use hydro_parameters
  use const
  use cr_parameters, only: ncr_groups,gamma_cr,cr_smallr_decouple
  implicit none
  integer::ncell
  real(dp)::dx,dt
  real(dp),dimension(1:nvector,1:nvar+3)::uu
  real(dp),dimension(1:nvector,1:ndim)::gg
  real(dp),dimension(1:nvector,1:ncr_groups)::ecr
  real(dp),dimension(1:nvector),save::a2,B2,rho,ctot,cr_cs

  real(dp)::dtcell,smallp,cf,cc,bc,bn
  integer::k,idim,icr
#if NENER>0
  integer::irad
#endif

  smallp = smallr*smallc**2/gamma

  ! Convert to primitive variables
  do k = 1,ncell
     uu(k,1)=max(uu(k,1),smallr)
     rho(k)=uu(k,1)
  end do
  do idim = 1,3
     do k = 1, ncell
        uu(k,idim+1) = uu(k,idim+1)/rho(k)
     end do
  end do

  do k = 1,ncell
     B2(k)=zero
  end do
  do idim = 1,3
     do k = 1, ncell
        Bc = half*(uu(k,neul+idim)+uu(k,nvar+idim))
        B2(k)=B2(k)+Bc**2
        uu(k,neul) = uu(k,neul)-half*uu(k,1)*uu(k,idim+1)**2-half*Bc**2
     end do
  end do
#if NENER>0
  do irad = 1,nener
     do k = 1, ncell
        uu(k,neul) = uu(k,neul)-uu(k,nhydro+irad)
     end do
  end do
#endif

  ! Compute thermal sound speed
  do k = 1, ncell
     uu(k,neul) = max((gamma-one)*uu(k,neul),smallp)
     a2(k)=gamma*uu(k,neul)/uu(k,1)
  end do
#if NENER>0
  do irad = 1,nener
     do k = 1, ncell
        a2(k) = a2(k) + gamma_rad(irad)*(gamma_rad(irad)-1)*uu(k,nhydro+irad)/uu(k,1)
     end do
  end do
#endif

  ! Add the CR sound speed where rho is not tiny
  do k = 1, ncell
     cr_cs(k)=0.0d0
  end do
  do icr = 1,ncr_groups
     do k = 1, ncell
        cr_cs(k)=cr_cs(k) + ecr(k,icr) * gamma_cr(icr)*(gamma_cr(icr)-1.0d0)
     end do
  end do
  do k = 1, ncell
     cr_cs(k)=cr_cs(k)/rho(k)
     if(rho(k) .gt. smallr*cr_smallr_decouple) &
          a2(k) = a2(k) + cr_cs(k)
  end do

  ! Compute maximum wave speed (fast magnetosonic)
  do k = 1, ncell
     ctot(k)=zero
  end do
  if(ischeme.eq.1)then
     do idim = 1,ndim   ! WARNING: ndim instead of 3
        do k = 1, ncell
           ctot(k)=ctot(k)+abs(uu(k,idim+1))
        end do
     end do
  else
     do idim = 1,ndim   ! WARNING: ndim instead of 3
        do k = 1, ncell
           cc=half*(B2(k)/rho(k)+a2(k))
           BN=half*(uu(k,5+idim)+uu(k,nvar+idim))
           cf=sqrt(cc+sqrt(cc**2-a2(k)*BN**2/rho(k)))
           ctot(k)=ctot(k)+abs(uu(k,idim+1))+cf
        end do
     end do
  endif

  ! Compute gravity strength ratio
  do k = 1, ncell
     rho(k)=zero
  end do
  do idim = 1,ndim
     do k = 1, ncell
        rho(k)=rho(k)+abs(gg(k,idim))
     end do
  end do
  do k = 1, ncell
     rho(k)=rho(k)*dx/ctot(k)**2
     rho(k)=MAX(rho(k),0.0001_dp)
  end do

  ! Compute maximum time step for each authorized cell
  dt = courant_factor*dx/smallc
  do k = 1,ncell
     dtcell=dx/ctot(k)*(sqrt(one+two*courant_factor*rho(k))-one)/rho(k)
     dt = min(dt,dtcell)
  end do

end subroutine cr_cmpdt
!#########################################################
!#########################################################
subroutine cr_cmpdt_vamax(uu,ncell)
  ! Track the maximum Alfven speed for the adaptive reduced light speed.
  use amr_parameters, only: dp,nvector
  use hydro_parameters, only: nvar,smallr
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
     cr_va_max=max(cr_va_max, sqrt(BNva/max(uu(k,1),smallr)))
  end do

end subroutine cr_cmpdt_vamax
