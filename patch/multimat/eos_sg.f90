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
  real(dp)::g0,p0,a0,b0
  real(dp),dimension(1:nvector),save::alpha_tot,beta_tot
  real(dp),dimension(1:nvector),save::gamma_tot,pinf_tot

  ! Initialize variables as zero and update them iteratively
  alpha_tot(1:ncell)=zero
  beta_tot (1:ncell)=zero
  do imat = 1,nmat
    ! Get stiffened gas EOS parameters
    g0=eos_params(imat,1); p0=eos_params(imat,2)
    a0=one/(g0-one); b0=p0*g0/(g0-one)
    ! p=(gamma-1)*e-gamma*pinf; e=beta+alpha*p
    do k = 1,ncell
        ! Update alpha_tot and beta_tot
        alpha_tot(k) = alpha_tot(k) + f(k,imat) * a0
        beta_tot (k) = beta_tot (k) + f(k,imat) * b0
    end do
  end do
  gamma_tot(1:ncell)=one/alpha_tot(1:ncell)+one
  pinf_tot (1:ncell)=beta_tot(1:ncell)/alpha_tot(1:ncell)/gamma_tot(1:ncell)
  ! Calculate the pressure for given internal energy
  do k = 1,ncell
     p(k) = (q(k,npri) - beta_tot(k)) / alpha_tot(k)
  end do
  ! Calculate the bulk modulus
  do imat = 1,nmat
    ! Get stiffened gas EOS parameters
    g0=eos_params(imat,1); p0=eos_params(imat,2)
    a0=one/(g0-one); b0=p0*g0/(g0-one)
    ! p=(gamma-1)*e-gamma*pinf; e=beta+alpha*p
    do k = 1,ncell
      kappa_mat(k,imat) = MAX(g(k,imat) * smallc**2, g0  * (p(k) + p0))
      kappa_hat(k) = kappa_hat(k) + f(k,imat) / kappa_mat(k,imat)
    end do
  end do
  kappa_hat(1:ncell)=one/kappa_hat(1:ncell)
  ! Calculate the speed of sound (old method)
  do k = 1,ncell
    c(k) = gamma_tot(k)*(p(k)+pinf_tot(k))/q(k,1)
    c(k) = sqrt(max(c(k),smallc**2))
  end do
  ! Calculate the total speed of sound (new method)
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
  real(dp)::g0,p0,a0,b0
  real(dp),dimension(1:nvector),save::alpha_tot,beta_tot
  real(dp),dimension(1:nvector),save::gamma_tot,pinf_tot

  ! Initialize variables as zero and update them iteratively
  alpha_tot(1:ncell)=zero
  beta_tot (1:ncell)=zero
  do imat = 1,nmat
    ! Get stiffened gas EOS parameters
    g0=eos_params(imat,1); p0=eos_params(imat,2)
    a0=one/(g0-one); b0=p0*g0/(g0-one)
    ! p=(gamma-1)*e-gamma*pinf; e=beta+alpha*p
    do k = 1,ncell
      ! Update alpha_tot and beta_tot
        alpha_tot(k) = alpha_tot(k) + f(k,imat) * a0
        beta_tot (k) = beta_tot (k) + f(k,imat) * b0
    end do
  end do
  gamma_tot(1:ncell)=one/alpha_tot(1:ncell)+one
  pinf_tot (1:ncell)=beta_tot(1:ncell)/alpha_tot(1:ncell)/gamma_tot(1:ncell)
  ! Calculate the interal energy for given pressure
  do k = 1,ncell
    e(k) = alpha_tot(k) * q(k,npri) + beta_tot(k)
  end do
  ! Calculate the bulk modulus
  do imat = 1,nmat
    ! Get stiffened gas EOS parameters
    g0=eos_params(imat,1); p0=eos_params(imat,2)
    a0=one/(g0-one); b0=p0*g0/(g0-one)
    ! p=(gamma-1)*e-gamma*pinf; e=beta+alpha*p
    do k = 1,ncell
      kappa_mat(k,imat) = MAX(g(k,imat) * smallc**2, g0  * (q(k,npri) + p0))
      kappa_hat(k) = kappa_hat(k) + f(k,imat) / kappa_mat(k,imat)
    end do
  end do
  kappa_hat(1:ncell)=one/kappa_hat(1:ncell)
  ! Calculate the speed of sound (old method)
  do k = 1,ncell
    c(k) = gamma_tot(k)*(q(k,npri)+pinf_tot(k))/q(k,1)
    c(k) = sqrt(max(c(k),smallc**2))
  end do
  ! Calculate the total speed of sound (new method)
  !   do k = 1,ncell
  !      c(k) = kappa_hat(k) / q(k,1)
  !      c(k) = sqrt(max(c(k),smallc**2))
  !   end do
end subroutine eosinv
