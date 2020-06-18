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
  real(dp)::smallgamma,biggamma,e_0,e_c,P_c, delpc
  real(dp)::rho_0,C_v,T_0,E_1,E_2,A_1,A_2
  real(dp),dimension(1:nvector),save::ec_hat,pc_hat, alpha_hat
  real(dp),dimension(1:nvector),save::gamma_hat,delpc_hat

  if(eos_name == 'stiffened gas')then       
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

   else if(eos_name == 'mie-grueneisen')then
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
            e_c = e_0 * (max(g(k,imat),smallr)/rho_0)**biggamma
            P_c = e_0 * (biggamma-one) * (max(g(k,imat),smallr)/rho_0)**biggamma
            delpc = (e_0/rho_0) * biggamma * (biggamma-one) * (max(g(k,imat),smallr)/rho_0)**(biggamma-one)
            ! Update total values
            alpha_hat(k) = alpha_hat(k) + f(k,imat) * a0
            ec_hat(k) = ec_hat(k) + f(k,imat) * e_c
            pc_hat(k) = pc_hat(k) + f(k,imat) * P_c * a0 
            delpc_hat(k) = delpc_hat(k) + f(k,imat) * delpc * a0 
         end do
      end do
      pc_hat(1:ncell) = pc_hat(1:ncell)/alpha_hat(1:ncell)
      delpc_hat(1:ncell) = delpc_hat(1:ncell)/alpha_hat(1:ncell)
      gamma_hat(1:ncell)=one/alpha_hat(1:ncell)+one
      do k = 1,ncell
         ! Calculate the pressure for given internal energy
         p(k) = max(q(k,npri) - ec_hat(k),0.0) / alpha_hat(k) + pc_hat(k)
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
            P_c = e_0 * (biggamma-one) * (max(g(k,imat),smallr)/rho_0)**biggamma
            delpc = (e_0/rho_0) * biggamma * (biggamma-one) * (max(g(k,imat),smallr)/rho_0)**(biggamma-one)
            kappa_mat(k,imat) = max(g(k,imat),smallr) * max(smallc**2,delpc+(smallgamma/max(g(k,imat),smallr))*max(p(k)-P_c,0.0))
            kappa_hat(k) = kappa_hat(k) + f(k,imat) / kappa_mat(k,imat)
         end do
      end do
      ! We have calculated 1/kappa_hat so we invert it
      kappa_hat(1:ncell) = one / kappa_hat(1:ncell)
      ! Calculate the speed of sound (old method)
      do k = 1,ncell
      c(k) = delpc_hat(k)+(gamma_hat(k)/max(q(k,1),smallr)*max(p(k)-pc_hat(k),0.0)
      c(k) = sqrt(max(c(k),smallc**2))
      end do
      ! Calculate the total speed of sound(new method)
      ! do k = 1,ncell 
      !    c(k) = kappa_hat(k) / q(k,1)
      !    c(k) = sqrt(max(c(k),smallc**2))
      ! end do

   else if(eos_name == 'cochran-chan')then
      kappa_hat(1:ncell)=zero
      alpha_hat(1:ncell)=zero
      pc_hat (1:ncell)=zero
      ec_hat(1:ncell)=zero
      delpc_hat(1:ncell)=zero
      gamma_hat(1:ncell)=zero
      ! Cochran-Chan EOS written in terms of the Mie-Grueneisen EOS 
      do imat = 1,nmat
         ! Get Mie-Grueneisen EOS parameters
         rho_0=eos_params(imat,1);C_v=eos_params(imat,2);T_0=eos_params(imat,3)
         E_1=eos_params(imat,4);E_2=eos_params(imat,5);A_1=eos_params(imat,6);A_2=eos_params(imat,7)
         smallgamma=eos_params(imat,8)
         ! P - P_c = (gamma - one) * (e - e_c) ; e = P/(gamma-1) + (e_c-P_c/(gamma-1))
         a0 = one / (smallgamma-one) 
         do k = 1,ncell
            ! Update Mie-Gruneisen terms for each material
            e_c = -A_1 / (rho_0*(1. - E_1))  * (((rho_0/g(k,imat))**(1.-E_1))-1.) + A_2 / (rho_0*(1. - E_2)) * (((rho_0/g(k,imat))**(1.-E_2))-1.) - C_v * T_0
            P_c = A_1 * (rho_0/g(k,imat))**(-E_1) - A_2 * (rho_0/g(k,imat))**(-E_2)
            delpc = A_1 * E_1 * (g(k,imat)**(E_1-1)/rho_0**E_1) - A_2 * E_2 * (g(k,imat)**(E_2-1)/rho_0**E_2)
            ! Update total values
            alpha_hat(k) = alpha_hat(k) + f(k,imat) * a0
            ec_hat(k) = ec_hat(k) + f(k,imat) * e_c
            pc_hat(k) = pc_hat(k) + (f(k,imat)*a0*P_c) 
            delpc_hat(k) = delpc_hat(k) + (f(k,imat)*a0*delpc) 
         end do
      end do
      pc_hat(1:ncell) = pc_hat(1:ncell)/alpha_hat(1:ncell)
      delpc_hat(1:ncell) = delpc_hat(1:ncell)/alpha_hat(1:ncell)
      gamma_hat(1:ncell) = one/alpha_hat(1:ncell)+one
      do k = 1,ncell
         ! Calculate the pressure for given internal energy
         p(k) = max(q(k,npri) - ec_hat(k),0.0) / alpha_hat(k) + pc_hat(k)
      end do 

      do imat = 1,nmat
         rho_0=eos_params(imat,1);C_v=eos_params(imat,2);T_0=eos_params(imat,3)
         E_1=eos_params(imat,4);E_2=eos_params(imat,5);A_1=eos_params(imat,6);A_2=eos_params(imat,7)
         smallgamma=eos_params(imat,8)
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
   end if
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
  real(dp)::smallgamma,biggamma,e_0,e_c,P_c, delpc
  real(dp)::rho_0,C_v,T_0,E_1,E_2,A_1,A_2
  real(dp),dimension(1:nvector),save::ec_hat,pc_hat, alpha_hat
  real(dp),dimension(1:nvector),save::gamma_hat,delpc_hat

   if(eos_name=='stiffened gas')then 
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

   else if(eos_name=='mie-grueneisen')then
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
            e_c = e_0 * (max(g(k,imat),smallr)/rho_0)**biggamma
            P_c = e_0 * (biggamma-one) * (max(g(k,imat),smallr)/rho_0)**biggamma
            delpc = (e_0/rho_0) * biggamma * (biggamma-one) * (max(g(k,imat),smallr)/rho_0)**(biggamma-one)
            ! Update total values
            alpha_hat(k) = alpha_hat(k) + f(k,imat) * a0
            ec_hat(k) = ec_hat(k) + f(k,imat) * e_c
            pc_hat(k) = pc_hat(k) + f(k,imat) * P_c * a0 
            delpc_hat(k) = delpc_hat(k) + f(k,imat) * delpc * a0 
         end do
      end do
      pc_hat(1:ncell) = pc_hat(1:ncell)/alpha_hat(1:ncell)
      delpc_hat(1:ncell) = delpc_hat(1:ncell)/alpha_hat(1:ncell)
      gamma_hat(1:ncell) = one/alpha_hat(1:ncell)+one
      ! Calculate the internal energy for given pressure
      do k=1,ncell
         e(k) = alpha_hat(k) * max(q(k,npri)-pc_hat(k),0.0) + ec_hat(k)
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
            P_c = e_0 * (biggamma-one) * (max(g(k,imat),smallr)/rho_0)**biggamma
            delpc = (e_0/rho_0) * biggamma * (biggamma-one) * (max(g(k,imat),smallr)/rho_0)**(biggamma-one)
            kappa_mat(k,imat) = g(k,imat) * MAX(smallc**2, delpc+(smallgamma/max(g(k,imat),smallr))*max(q(k,npri)-P_c,0.0))
            kappa_hat(k) = kappa_hat(k) + f(k,imat) / kappa_mat(k,imat)
         end do
      end do
      ! We have calculated 1/kappa_hat so we invert it
      kappa_hat(1:ncell) = one / kappa_hat(1:ncell)
      ! Calculate the speed of sound (old method)
      do k = 1,ncell
      c(k) = delpc_hat(k)+(gamma_hat(k)/q(k,1))*max(q(k,npri)-pc_hat(k),0.0)
      c(k) = sqrt(max(c(k),smallc**2))
      end do
      ! Calculate the total speed of sound (new method)
      ! do k = 1,ncell 
      !    c(k) = kappa_hat(k) / q(k,1)
      !    c(k) = sqrt(max(c(k),smallc**2))
      ! end do

   else if(eos_name == 'cochran-chan')then
      kappa_hat(1:ncell)=zero
      alpha_hat(1:ncell)=zero
      pc_hat (1:ncell)=zero
      ec_hat(1:ncell)=zero
      delpc_hat(1:ncell)=zero
      gamma_hat(1:ncell)=zero
      ! Cochran-Chan EOS written in terms of the Mie-Grueneisen EOS 
      do imat = 1,nmat
         ! Get Mie-Grueneisen EOS parameters
         rho_0=eos_params(imat,1);C_v=eos_params(imat,2);T_0=eos_params(imat,3)
         E_1=eos_params(imat,4);E_2=eos_params(imat,5);A_1=eos_params(imat,6);A_2=eos_params(imat,7)
         smallgamma=eos_params(imat,8)
         ! P - P_c = (gamma - one) * (e - e_c) ; e = P/(gamma-1) + (e_c-P_c/(gamma-1))
         a0 = one / (smallgamma-one) 
         do k = 1,ncell
            ! Update Mie-Gruneisen terms for each material
            e_c = -A_1 / (rho_0*(1. - E_1))  * (((rho_0/g(k,imat))**(1.-E_1))-1.) + A_2 / (rho_0*(1. - E_2)) * (((rho_0/g(k,imat))**(1.-E_2))-1.) - C_v * T_0
            P_c = A_1 * (rho_0/g(k,imat))**(-E_1) - A_2 * (rho_0/g(k,imat))**(-E_2)
            delpc = A_1 * E_1 * (g(k,imat)**(E_1-1)/rho_0**E_1) - A_2 * E_2 * (g(k,imat)**(E_2-1)/rho_0**E_2)
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
         rho_0=eos_params(imat,1);C_v=eos_params(imat,2);T_0=eos_params(imat,3)
         E_1=eos_params(imat,4);E_2=eos_params(imat,5);A_1=eos_params(imat,6);A_2=eos_params(imat,7)
         smallgamma=eos_params(imat,8)
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
   end if
end subroutine eosinv

!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmpdt(uu,grav,rr,dx,dt,ncell)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none
  integer::ncell
  real(dp)::dx,dt
  real(dp),dimension(1:nvector,1:nvar)::uu
  real(dp),dimension(1:nvector,1:ndim)::grav
  real(dp),dimension(1:nvector)::rr
  
  real(dp),dimension(1:nvector,1:npri),save::qq
  real(dp),dimension(1:nvector,1:nmat),save::ff,gg,kappa_matt
  real(dp),dimension(1:nvector),save::ekin,dtot,cc,st,pp,kappa_hatt
  real(dp)::dtcell,eps
  integer::k,idim,imat

  ! Convert to primitive variable
     
  ! Volume fraction and fluid density
  do imat = 1,nmat
     do k = 1,ncell
        ff(k,imat) = uu(k,imat+npri)
        gg(k,imat) = uu(k,imat+npri+nmat)
     end do
  end do
     
  ! Compute density
  do k = 1,ncell
     qq(k,1) = uu(k,1)
  end do
  
  ! Compute velocity and specific kinetic energy
  ekin(1:ncell)=0.0
  do idim = 1,ndim
     do k = 1,ncell
        qq(k,idim+1) = uu(k,idim+1)/uu(k,1)
        ekin(k) = ekin(k) + half*qq(k,idim+1)**2
     end do
  end do
  
  ! Compute total internal energy
  do k = 1,ncell
     qq(k,npri) = uu(k,npri) - uu(k,1)*ekin(k)
  end do
  
  ! Call eos routine
  call eos(ff,gg,qq,pp,cc,kappa_matt,kappa_hatt,ncell)
  
  ! Compute wave speed
  if(geom==3)then
     do k = 1,ncell
        eps = dx/two/rr(k)
        cc(k) = (abs(qq(k,2))+cc(k))*(one+eps)**2/(one+third*eps**2)
     end do
  else if(geom==2)then
     do k = 1,ncell
        eps = dx/two/rr(k)
        cc(k) = (abs(qq(k,2))+cc(k))*(one+eps)
     end do
  else
     do k = 1,ncell
        cc(k) = abs(qq(k,2))+cc(k)
     end do
  endif
  do idim = 2,ndim
     do k = 1,ncell 
        cc(k) = cc(k) + abs(qq(k,idim+1))+cc(k)
     end do
  end do

  ! Compute gravity strength ratio
  do k = 1,ncell
     st(k) = zero
  end do
  do idim = 1,ndim
     do k = 1,ncell 
        st(k) = st(k) + abs(grav(k,idim))
     end do
  end do
  do k = 1,ncell
     st(k) = st(k)*dx/cc(k)**2
     st(k) = MAX(st(k),0.0001_dp)
  end do

  ! Compute maximum time step for each authorized cell
  dt = courant_factor*dx/smallc
  do k = 1,ncell
     dtcell = dx/cc(k)*(sqrt(one+two*courant_factor*st(k))-one)/st(k)
     dt = min(dt,dtcell)
  end do

end subroutine cmpdt
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine hydro_refine(ug,um,ud,ok,current_dim,ncell)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none
  ! dummy arguments
  integer::ncell,current_dim
  real(dp),dimension(1:nvector,1:nvar)::ug,um,ud
  logical ,dimension(1:nvector)       ::ok
  
  integer::k,idim,imat
  real(dp),dimension(1:nvector,1:npri),save::qg,qm,qd
  real(dp),dimension(1:nvector,1:nmat),save::fg,fm,fd
  real(dp),dimension(1:nvector,1:nmat),save::gg,gm,gd
  real(dp),dimension(1:nvector,1:nmat),save::kappa_matg,kappa_matm,kappa_matd
  real(dp),dimension(1:nvector),save::eking,ekinm,ekind
  real(dp),dimension(1:nvector),save::pg,pm,pd
  real(dp),dimension(1:nvector),save::cg,cm,cd
  real(dp),dimension(1:nvector),save::kappa_hatg,kappa_hatm,kappa_hatd
  logical ,dimension(1:nvector),save::wg,wd,bg,bm,bd
  real(dp)::ffg,ffm,ffd,ddg,ddm,ddd
  real(dp)::ppg,ppm,ppd,vvg,vvm,vvd
  real(dp)::ccg,ccm,ccd,error
  
  ! Convert to primitive variables
     
  ! Volume fraction and fluid density
  do imat = 1,nmat
     do k = 1,ncell
        fg(k,imat) = ug(k,imat+npri)
        fm(k,imat) = um(k,imat+npri)
        fd(k,imat) = ud(k,imat+npri)
        gg(k,imat) = ug(k,imat+npri+nmat)
        gm(k,imat) = um(k,imat+npri+nmat)
        gd(k,imat) = ud(k,imat+npri+nmat)
     end do
  end do
     
  ! Detect embedded body
  if(static)then
     do k=1,ncell
        bg(k)=fg(k,1)>0.01
        bm(k)=fm(k,1)>0.01
        bd(k)=fd(k,1)>0.01
     end do
     do k=1,ncell
        if(bm(k))then
           fg(k,1:nmat)=fm(k,1:nmat)
           gg(k,1:nmat)=gm(k,1:nmat)
           ug(k,1:nvar)=um(k,1:nvar)
           fd(k,1:nmat)=fm(k,1:nmat)
           gd(k,1:nmat)=gm(k,1:nmat)
           ud(k,1:nvar)=um(k,1:nvar)
        else
           if(bg(k))then
              fg(k,1:nmat)=fm(k,1:nmat)
              gg(k,1:nmat)=gm(k,1:nmat)
              ug(k,1:nvar)=um(k,1:nvar)
              ug(k,current_dim+1)=-um(k,current_dim+1)
           endif
           if(bd(k))then
              fd(k,1:nmat)=fm(k,1:nmat)
              gd(k,1:nmat)=gm(k,1:nmat)
              ud(k,1:nvar)=um(k,1:nvar)
              ud(k,current_dim+1)=-um(k,current_dim+1)
           endif
        endif
     enddo
  endif

  ! Compute total density
  do k = 1,ncell
     qg(k,1) = ug(k,1)
     qm(k,1) = um(k,1)
     qd(k,1) = ud(k,1)
  end do
  
  ! Compute velocity and specific kinetic energy
  eking(1:ncell)=0.0
  ekinm(1:ncell)=0.0
  ekind(1:ncell)=0.0
  do idim = 1,ndim
     do k = 1,ncell
        qg(k,idim+1) = ug(k,idim+1)/ug(k,1)
        qm(k,idim+1) = um(k,idim+1)/um(k,1)
        qd(k,idim+1) = ud(k,idim+1)/ud(k,1)
        eking(k) = eking(k) + half*qg(k,idim+1)**2
        ekinm(k) = ekinm(k) + half*qm(k,idim+1)**2
        ekind(k) = ekind(k) + half*qd(k,idim+1)**2
     end do
  end do
  
  ! Compute total internal energy
  do k = 1,ncell
     qg(k,npri) = ug(k,npri) - qg(k,1)*eking(k)
     qm(k,npri) = um(k,npri) - qm(k,1)*ekinm(k)
     qd(k,npri) = ud(k,npri) - qd(k,1)*ekind(k)
  end do
  
  ! Call eos routine
  call eos(fg,gg,qg,pg,cg,kappa_matg,kappa_hatg,ncell)
  call eos(fm,gm,qm,pm,cm,kappa_matm,kappa_hatm,ncell)
  call eos(fd,gd,qd,pd,cd,kappa_matd,kappa_hatd,ncell)
  
  ! Compute errors
  if(err_grad_d >= 0.)then
     do k=1,ncell
        ddg=abs(qg(k,1)); ddm=abs(qm(k,1)); ddd=abs(qd(k,1))
        error=2.0d0*MAX( &
             & ABS((ddd-ddm)/(ddd+ddm+floor_d)) , &
             & ABS((ddm-ddg)/(ddm+ddg+floor_d)) )
        ok(k) = ok(k) .or. error > err_grad_d
     end do
  end if

  if(err_grad_f >= 0.)then
     do imat=1,nmat
        do k=1,ncell
           ffg=fg(k,imat); ffm=fm(k,imat); ffd=fd(k,imat)
           error=2.0d0*MAX( &
                & ABS((ffd-ffm)/(ffd+ffm+floor_f)) , &
                & ABS((ffm-ffg)/(ffm+ffg+floor_f)) )
           ok(k) = ok(k) .or. error > err_grad_f
        end do
     end do
  end if

  if(err_grad_p > -1.0)then
     do k=1,ncell
        ppg=pg(k); ppm=pm(k); ppd=pd(k)
        error=2.0d0*MAX( &
             & ABS((ppd-ppm)/(ppd+ppm+floor_p)), &
             & ABS((ppm-ppg)/(ppm+ppg+floor_p)) )
        ok(k) = ok(k) .or. error > err_grad_p
     end do
  end if

  if(err_grad_u >= 0.)then
     do idim = 1,ndim
        do k=1,ncell
           vvg=qg(k,idim+1); vvm=qm(k,idim+1); vvd=qd(k,idim+1)
           ccg=cg(k)       ; ccm=cm(k)       ; ccd=cd(k)
           error=2.0d0*MAX( &
                & ABS((vvd-vvm)/(ccd+ccm+ABS(vvd)+ABS(vvm)+floor_u)) , &
                & ABS((vvm-vvg)/(ccm+ccg+ABS(vvm)+ABS(vvg)+floor_u)) )
           ok(k) = ok(k) .or. error > err_grad_u
        end do
     end do
  end if

!!$  if(static)then
!!$     do k=1,ncell
!!$        if(wg(k).or.wd(k))then
!!$           ddg=abs(qg(k,1)); ddm=abs(qm(k,1)); ddd=abs(qd(k,1))
!!$           error=2.0d0*MAX( &
!!$                & ABS((ddd-ddm)/(ddd+ddm+floor_d)) , &
!!$                & ABS((ddm-ddg)/(ddm+ddg+floor_d)) )
!!$           write(*,*)wg(k),wd(k)
!!$           write(*,*)bg(k),bd(k)           
!!$           write(*,*)ddg,ddm,ddd
!!$           write(*,*)error,ok(k)
!!$        endif
!!$     end do
!!$  endif

end subroutine hydro_refine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine riemann_acoustic(fl,fr,gl,gr,ql,qr,cl,cr,fgdnv,ggdnv,qgdnv,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  ! dummy arguments
  integer::ngrid
  real(dp),dimension(1:nvector,1:npri)::ql,qr,qgdnv
  real(dp),dimension(1:nvector,1:nmat)::fl,fr,fgdnv
  real(dp),dimension(1:nvector,1:nmat)::gl,gr,ggdnv
  real(dp),dimension(1:nvector,1:nmat)::kappa_matl,kappa_matr
  real(dp),dimension(1:nvector)::cl,cr
  real(dp),dimension(1:nvector)::kappa_hatl,kappa_hatr

  ! local variables
  integer::i,n,ir,ie,imat
  real(dp)::smallp,delp_p

  ! local arrays
  real(dp),dimension(1:nvector,1:npri),save::qo,qstar
  real(dp),dimension(1:nvector,1:nmat),save::fo,fstar
  real(dp),dimension(1:nvector,1:nmat),save::go,gstar
  real(dp),dimension(1:nvector,1:nmat),save::kappa_matstar
  real(dp),dimension(1:nvector),save::cstar,estar,pstar,ustar
  real(dp),dimension(1:nvector),save::kappa_hatstar
  real(dp),dimension(1:nvector),save::sgnm ,spin ,spout,ushock
  real(dp),dimension(1:nvector),save::frac,co,el,er
  logical ,dimension(1:nvector),save::wall
  real(dp)::wl,wr,ul,ur,pl,pr,dl,dr,dstar

  ! Sound speed
  call eosinv(fl,gl,ql,el,cl,kappa_matl,kappa_hatl,ngrid)
  call eosinv(fr,gr,qr,er,cr,kappa_matr,kappa_hatr,ngrid)

  ! Acoustic star state
  do i=1,ngrid
     dl = ql(i,1); dr = qr(i,1)
     ul = ql(i,2); ur = qr(i,2)
     pl = ql(i,npri); pr = qr(i,npri)
     wl = cl(i)*dl
     wr = cr(i)*dr
     pstar(i) = ( (wr*pl+wl*pr)+wl*wr*(ul-ur) ) / (wl+wr)
     ustar(i) = ( (wr*ur+wl*ul)+      (pl-pr) ) / (wl+wr)
  end do

  ! Left going or right going contact wave
  do i=1,ngrid   
     sgnm(i) = sign(one,ustar(i))
  end do

  ! Left or right unperturbed state
  do i = 1,ngrid
     if(sgnm(i)==one)then
        fo(i,1:nmat) = fl(i,1:nmat)
        go(i,1:nmat) = gl(i,1:nmat)
        qo(i,1:npri) = ql(i,1:npri)
        co(i) = cl(i)
     else
        fo(i,1:nmat) = fr(i,1:nmat)
        go(i,1:nmat) = gr(i,1:nmat)
        qo(i,1:npri) = qr(i,1:npri)
        co(i) = cr(i)
     end if
  end do

  ! Star region density, internal energy and sound speed
  do i = 1,ngrid
     dstar = qo(i,1) + (pstar(i)-qo(i,npri))/co(i)**2
     dstar = max(dstar,smallr)
     fstar(i,1:nmat) = fo(i,1:nmat)
     gstar(i,1:nmat) = go(i,1:nmat)
     qstar(i,1) = dstar
     qstar(i,2) = ustar(i)
     qstar(i,npri) = pstar(i)
#if NDIM>1
     qstar(i,3) = qo(i,3)
#endif
#if NDIM>2
     qstar(i,4) = qo(i,4)
#endif
  end do
  call eosinv(fstar,gstar,qstar,estar,cstar,kappa_matstar,kappa_hatstar,ngrid)

  ! Head and tail speed of rarefaction
  do i=1,ngrid
     spout(i) = co   (i)-sgnm(i)*qo   (i,2)
     spin (i) = cstar(i)-sgnm(i)*qstar(i,2)
  end do

  ! Shock speed
  do i=1,ngrid
     ushock(i) = half*(spin(i)+spout(i))
     ushock(i) = max(ushock(i),-sgnm(i)*qstar(i,2))
  end do
  do i=1,ngrid
     if(pstar(i)>=qo(i,npri))then
        spout(i)=ushock(i)
        spin (i)=spout (i)
     end if
  end do

  ! Sample the solution at x/t=0
  do i=1,ngrid
     if(spout(i)<zero)then      ! Initial state
        qgdnv(i,1:npri) = qo(i,1:npri)
        fgdnv(i,1:nmat) = fo(i,1:nmat)
        ggdnv(i,1:nmat) = go(i,1:nmat)
     else if(spin(i)>=zero)then  ! Star region
        qgdnv(i,1:npri) = qstar(i,1:npri)
        fgdnv(i,1:nmat) = fstar(i,1:nmat)
        ggdnv(i,1:nmat) = gstar(i,1:nmat)
     else                        ! Rarefaction
        frac(i) = spout(i)/(spout(i)-spin(i))
        qgdnv(i,1:npri) = frac(i)*qstar(i,1:npri) + (one - frac(i))*qo(i,1:npri)
        fgdnv(i,1:nmat) = frac(i)*fstar(i,1:nmat) + (one - frac(i))*fo(i,1:nmat)
        ggdnv(i,1:nmat) = frac(i)*gstar(i,1:nmat) + (one - frac(i))*go(i,1:nmat)
     end if
  end do

end subroutine riemann_acoustic
!###########################################################
!###########################################################
!###########################################################
!###########################################################

