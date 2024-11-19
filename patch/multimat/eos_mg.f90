!###########################################################
!###########################################################
!###########################################################
!###########################################################

subroutine eos(f,g,q,p,c,kappa_mat,kappa_hat,ncell)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none
  integer::ncell
  real(dp),dimension(1:nvector,1:nmat)::f,g,kappa_mat
  real(dp),dimension(1:nvector,1:npri)::q
  real(dp),dimension(1:nvector)::p,c,kappa_hat
  ! Compute total pressure and sound speed from total internal energy
  ! On entry:
  ! f is the volume fraction of each fluid
  ! g is the true density of each fluid
  ! q are the total primitive variables (d, u, eint)
  ! On exit:
  ! p is the total pressure
  ! c is the total sound speed
  integer::k,imat
  real(dp)::p0,a0
  real(dp)::smallgamma,biggamma,e_0, rho_0,e_c,P_c, delpc
  real(dp),dimension(1:nvector),save::ec_hat,pc_hat, alpha_hat
  real(dp),dimension(1:nvector),save::gamma_hat,delpc_hat

  kappa_hat(1:ncell)=zero
  alpha_hat(1:ncell)=zero
  pc_hat (1:ncell)=zero
  ec_hat(1:ncell)=zero
  delpc_hat(1:ncell)=zero
  gamma_hat(1:ncell)=zero
  ! Mie-Gruneisen implementation trial
  do imat = 1,nmat
    ! Get Mie-Grueneisen EOS parameters
    smallgamma=eos_params(imat,1);biggamma=eos_params(imat,2);e_0=eos_params(imat,3);rho_0=eos_params(imat,4)
    ! Test case for checking whether the model works
    ! smallgamma=eos_params(imat,1); biggamma=0 ; p0=eos_params(imat,2)
    ! e_0 =  p0
    ! P - P_c = (gamma - one) * (e - e_c) ; e = P/(gamma-1) + (e_c-P_c/(gamma-1))
    a0 = one / (smallgamma-one)
    do k = 1,ncell
      ! Update Mie-Gruneisen terms for each material
      ! print *, g(k,imat)
      e_c = e_0 * (g(k,imat)/rho_0)**biggamma
      P_c = e_0 * (biggamma-one) * (g(k,imat)/rho_0)**biggamma
      delpc = (e_0/rho_0) * biggamma * (biggamma-one) * (g(k,imat)/rho_0)**(biggamma-one)
      ! Update total values
      alpha_hat(k) = alpha_hat(k) + f(k,imat) * a0
      ec_hat(k) = ec_hat(k) + f(k,imat) * e_c
      pc_hat(k) = pc_hat(k) + (f(k,imat)*a0*P_c)
      delpc_hat(k) = delpc_hat(k) + (f(k,imat)*a0*delpc)
    end do
  end do
  pc_hat(1:ncell) = pc_hat(1:ncell)/alpha_hat(1:ncell)
  delpc_hat(1:ncell) = delpc_hat(1:ncell)/alpha_hat(1:ncell)
  gamma_hat(1:ncell)=one/alpha_hat(1:ncell)+one
  do k = 1,ncell
    ! Calculate the pressure for given internal energy
    p(k) = (q(k,npri) - ec_hat(k)) / alpha_hat(k) + pc_hat(k)
  end do

  do imat = 1,nmat
    smallgamma=eos_params(imat,1);biggamma=eos_params(imat,2);e_0=eos_params(imat,3);rho_0=eos_params(imat,4)
    ! smallgamma=eos_params(imat,1); biggamma=0 ; p0=eos_params(imat,2)
    ! e_0 =  p0
    ! P - P_c = (gamma - one) * (e - e_c) ; e = P/(gamma-1) + (e_c-P_c/(gamma-1))
    a0 = one / (smallgamma-one)
    do k = 1,ncell
      ! Calculate the bulk moduli
      ! c_mat**2 = P_c' + smallgamma/rho * (P-P_c)
      kappa_mat(k,imat) = g(k,imat) * MAX(smallc**2, delpc+(smallgamma/g(k,imat))*(p(k)-P_c))
      kappa_hat(k) = kappa_hat(k) + f(k,imat) / kappa_mat(k,imat)
    end do
  end do
  ! We have calculated 1/kappa_hat so we invert it
  kappa_hat(1:ncell) = one / kappa_hat(1:ncell)
  ! Calculate the speed of sound (old method)
  do k = 1,ncell
    c(k) = delpc_hat(k)+(gamma_hat(k)/q(k,1))*(p(k)-pc_hat(k))
    c(k) = sqrt(max(c(k),smallc**2))
  end do
  ! Calculate the total speed of sound(new method)
  ! do k = 1,ncell
  !    c(k) = kappa_hat(k) / q(k,1)
  !    c(k) = sqrt(max(c(k),smallc**2))
  ! end do

end subroutine eos

!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine eosinv(f,g,q,e,c,kappa_mat,kappa_hat,ncell)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none
  integer::ncell
  real(dp),dimension(1:nvector,1:nmat)::f,g,c_mat,kappa_mat
  real(dp),dimension(1:nvector,1:npri)::q
  real(dp),dimension(1:nvector)::e,c,kappa_hat
  ! Compute total internal energy and sound speed from total pressure
  ! On entry:
  ! f is the volume fraction of each fluid
  ! g is the true density of each fluid
  ! q are the total primitive variables (d, u, P)
  ! On exit:
  ! eint is the total internal energy
  ! c is the sound speed of each fluid
  ! c_hat is the total sound speed
  integer::k,imat
  real(dp)::p0,a0
  real(dp)::smallgamma,biggamma,e_0,rho_0,e_c,P_c, delpc
  real(dp),dimension(1:nvector),save::ec_hat,pc_hat, alpha_hat
  real(dp),dimension(1:nvector),save::gamma_hat,delpc_hat

  kappa_hat(1:ncell)=zero
  alpha_hat(1:ncell)=zero
  pc_hat (1:ncell)=zero
  ec_hat(1:ncell)=zero
  delpc_hat(1:ncell)=zero
  gamma_hat(1:ncell)=zero
  ! Mie-Gruneisen implementation trial
  do imat = 1,nmat
    ! Get Mie-Grueneisen EOS parameters
    smallgamma=eos_params(imat,1);biggamma=eos_params(imat,2);e_0=eos_params(imat,3);rho_0=eos_params(imat,4)
    ! P - P_c = (gamma - one) * (e - e_c) ; e = P/(gamma-1) + (e_c-P_c/(gamma-1))
    ! Test case for checking whether the model works
    ! smallgamma=eos_params(imat,1); biggamma=0 ; p0=eos_params(imat,2)
    ! e_0 =  p0
    ! P - P_c = (gamma - one) * (e - e_c) ; e = P/(gamma-1) + (e_c-P_c/(gamma-1))
    a0 = one / (smallgamma-one)
    do k = 1,ncell
      ! Update Mie-Gruneisen terms for each material
      e_c = e_0 * (g(k,imat)/rho_0)**biggamma
      P_c = e_0 * (biggamma-one) * (g(k,imat)/rho_0)**biggamma
      delpc = (e_0/rho_0) * biggamma * (biggamma-one) * (g(k,imat)/rho_0)**(biggamma-one)
      ! Update total values
      alpha_hat(k) = alpha_hat(k) + f(k,imat) * a0
      ec_hat(k) = ec_hat(k) + f(k,imat) * e_c
      pc_hat(k) = pc_hat(k) + (f(k,imat)*a0*P_c)
      delpc_hat(k) = delpc_hat(k) + (f(k,imat)*a0*delpc)
      end do
  end do
  pc_hat(1:ncell) = pc_hat(1:ncell)/alpha_hat(1:ncell)
  delpc_hat(1:ncell) = delpc_hat(1:ncell)/alpha_hat(1:ncell)
  gamma_hat(1:ncell)=one/alpha_hat(1:ncell)+one
  ! Calculate the internal energy for given pressure
  do k=1,ncell
    e(k) = alpha_hat(k) * (q(k,npri)-pc_hat(k)) + ec_hat(k)
  end do

  do imat = 1,nmat
    smallgamma=eos_params(imat,1);biggamma=eos_params(imat,2);e_0=eos_params(imat,3);rho_0=eos_params(imat,4)
    ! smallgamma=eos_params(imat,1); biggamma=0 ; p0=eos_params(imat,2)
    ! e_0 =  p0
    ! P - P_c = (gamma - one) * (e - e_c) ; e = P/(gamma-1) + (e_c-P_c/(gamma-1))
    a0 = one / (smallgamma-one)
    do k = 1,ncell
      ! Calculate the bulk moduli
      ! c_mat**2 = P_c' + smallgamma/rho * (P-P_c)
      kappa_mat(k,imat) = g(k,imat) * MAX(smallc**2, delpc+(smallgamma/g(k,imat))*(q(k,npri)-P_c))
      kappa_hat(k) = kappa_hat(k) + f(k,imat) / kappa_mat(k,imat)
    end do
  end do
  ! We have calculated 1/kappa_hat so we invert it
  kappa_hat(1:ncell) = one / kappa_hat(1:ncell)
  ! Calculate the speed of sound (old method)
  do k = 1,ncell
    c(k) = delpc_hat(k)+(gamma_hat(k)/q(k,1))*(q(k,npri)-pc_hat(k))
    c(k) = sqrt(max(c(k),smallc**2))
  end do
  ! Calculate the total speed of sound (new method)
  ! do k = 1,ncell
  !    c(k) = kappa_hat(k) / q(k,1)
  !    c(k) = sqrt(max(c(k),smallc**2))
  ! end do
end subroutine eosinv
