module resistivity_table
   use amr_parameters, ONLY:dp
   use nimhd_parameters
   implicit none

   ! number of steps in the table for
   integer :: bchimie       ! B field (currently only used for reading table)
   integer :: nchimie       ! density 
   integer :: tchimie       ! temperature

   ! step sizes in table variables for
   real(dp) :: dbchimie     ! B field
   real(dp) :: dnchimie     ! density 
   real(dp) :: dtchimie     ! temperature

   ! minimum values for extrapolation variables in table for
   real(dp) :: bminchimie   ! B field
   real(dp) :: nminchimie   ! density  
   real(dp) :: tminchimie   ! temperature

   real(dp),allocatable,dimension(:,:,:,:)::resistivite_chimie ! resistivites chimie
   real(dp),allocatable,dimension(:,:,:)::resistivite_chimie_x ! to read in resistivity table

contains

subroutine read_resistivities
   use nimhd_parameters
   use constants, ONLY:mH
   !----------------------------------------------------------------
   ! Read non-ideal MHD resistivities from file and construct table
   !----------------------------------------------------------------
   integer::i,j,k
   integer   :: Nvarchimie        ! number of chemical species = ions + 3*grains (neutral,+,-) 
   real(dp)::bmaxchimie,dummy
   real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
   call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

   ! read header: get the number of bins for rho and temperature, as well as number of chemical species included
   open(42,file='resnh.dat', status='old')
   read(42,*) nchimie, tchimie, nvarchimie
   read(42,*)
   read(42,*)
   allocate(resistivite_chimie_x(-1:nvarchimie,nchimie,tchimie))

   ! read table values
   do i=1,tchimie
      do j=1,nchimie
         read(42,*)resistivite_chimie_x(0:nvarchimie,j,i),dummy,dummy,dummy,dummy,resistivite_chimie_x(-1,j,i)
      end do
      read(42,*)
   end do
   close(42)

   rho_threshold=max(rho_threshold,resistivite_chimie_x(0,1,1)*(mu_gas*mH)/scale_d) ! input in part/cc, output in code units
   nminchimie=(resistivite_chimie_x(0,1,1))
   dnchimie=(log10(resistivite_chimie_x(0,nchimie,1))-log10(resistivite_chimie_x(0,1,1)))/&
            &(nchimie-1)
   tminchimie=(resistivite_chimie_x(-1,1,1))
   dtchimie=(log10(resistivite_chimie_x(-1,1,tchimie))-log10(resistivite_chimie_x(-1,1,1)))/&
            &(tchimie-1)

   ! values for Btable
   bminchimie=1d-10
   bmaxchimie=1d5               ! ok for first core in nimhd. maybe not enough for second core.
   bchimie=100
   dbchimie=(log10(bmaxchimie)-log10(bminchimie))/real((bchimie-1),dp)
   allocate(resistivite_chimie(0:3,1:nchimie,1:tchimie,1:bchimie)) !memory leak

   ! construct table
   call construct_resistivity_table(nvarchimie)

   deallocate(resistivite_chimie_x)


end subroutine read_resistivities
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine construct_resistivity_table(nvarchimie)

   use amr_commons, only : myid
   use constants, only:pi,c_cgs, kB
   implicit none

   ! ---- chemistry ---------
   integer   :: nion=7            ! number of ions
   integer   :: nbins_grains      ! number of grain sizes
   real(dp)  :: nbins_real
   integer   :: Nvarchimie        ! number of chemical species = ions + 3*grains (neutral,+,-) 
   ! species charge and mass
   real(dp), allocatable, dimension(:) :: q
   real(dp), allocatable, dimension(:) :: m
   ! grain radius and mass for each bin
   real(dp), allocatable, dimension(:) :: r_g
   real(dp), allocatable, dimension(:) :: m_g
   ! grain density in g/cc
   real(dp), parameter :: rho_s=2.3_dp
   ! grain size distribution parameters, cf Kunz & Mouschovias 2009
   real(dp), parameter :: a_0=0.0375d-4      ! cm
   real(dp), parameter :: a_min=0.0181d-4    ! cm
   real(dp), parameter :: a_max=0.9049d-4    ! cm
   real(dp), parameter :: zeta=a_min/a_max   ! a_min/a_max
   real(dp), parameter :: lambda_pow=-3.5d0  ! Coeff power law
   ! other constants
   real(dp), parameter :: me=9.1094d-28 ! masse de l'electron, en g
   real(dp), parameter :: mp=1.6726d-24 ! masse du proton, en g
   real(dp), parameter :: e=4.803204d-10 ! charge de l'electron, en cgs
   !--------------------------

   integer  :: iB,iH,iT,i
   real(dp) :: B,nH,T,sigH,sigO,sigP

   ! resistivites (cf Kunz & Mouschovias 2009)
   real(dp), allocatable, dimension(:)  :: sigma         ! sigma_s
   real(dp), allocatable, dimension(:)  :: tau_sn
   real(dp), allocatable, dimension(:)  :: omega

   ! table contains values for neutral, positive and negatively charged grains for each size bin
   ! so we divide by 3 to get the actual amount of size bins
   nbins_real=real(nvarchimie-nion,dp)/3.0_dp
   nbins_grains=floor(nbins_real)
   if (nbins_real.ne.real(nbins_grains,dp)) then
      print*, 'issue in number of species'
      stop
   endif

   allocate(r_g(nbins_grains)) ! grain sizes
   allocate(m_g(nbins_grains)) ! grain masses
   allocate(m(nvarchimie))     ! particle (ions + grains) masses
   allocate(q(nvarchimie))     ! particle (ions + grains) charges
   allocate(sigma(nvarchimie))
   allocate(tau_sn(nvarchimie))
   allocate(omega(nvarchimie))

   ! determine grain sizes
   if(nbins_grains==1) then
     ! if we have only 1 bin, we take this average value
     r_g(1)=a_0
   else
     do  i=1,nbins_grains    ! cf Kunz & Mouschovias 2009
         r_g(i)=a_min*zeta**(-(i-1.d0)/nbins_grains)*(5.d0*(1.d0-zeta**(0.5/nbins_grains))/(1.d0-zeta**(2.5/nbins_grains)))**0.5
     end do
   end if

   ! set particle charges
   q(:)=1d0*e    ! cations
   q(1)=-1d0*e   ! electron
   do  i=nion+1,Nvarchimie
      if (mod(i-nion,3)==0) q(i)=0d0     ! neutral grains
      if (mod(i-nion,3)==1) q(i)=1d0*e   ! positively charged grains
      if (mod(i-nion,3)==2) q(i)=-1d0*e  ! negatively charged grains
   end do

   ! set particle and grain masses
   m(:) = 0d0
   m(1) = me               ! e-
   m(2) = 23.5d0*mp        ! ions metalliques
   m(3) = 29d0*mp          ! ions moleculaires
   m(4) = 3*mp             ! H3+
   m(5) = mp               ! H+
   m(6) = 12d0*mp          ! C+
   m(7) = 4d0*mp           ! He+
   do i=1,nbins_grains     ! masse des grains
      m_g(i)=4d0/3d0*pi*r_g(i)**3*rho_s
      m(nion+3*(i-1)+1:nion+3*i)=m_g(i) !mass for neutral, postive and negative grain is the same
   end do


   ! Compute table values
   if(myid==1) write(*,*) 'Computing 3D resistivities table'
   tau_sn      = 0.0_dp
   omega       = 0.0_dp
   sigma       = 0.0_dp
      
   do  iB=1,bchimie
   do  iT=1,tchimie
   do  iH=1,nchimie

      nh=resistivite_chimie_x(0,iH,iT)  ! density (.cc) of current point
      T =resistivite_chimie_x(-1,iH,iT) ! temperature
      B =10.0d0**(log10(bminchimie)+dble(iB-1)*dbchimie)

      do i=1,nion
         if (i==1) then
            ! electron
            tau_sn(i) = 1d0/1.16d0*(m(i)+2d0*mp)/(2d0*mp)*1d0/(nH/2d0*1.3d-9)
        else
            ! ions
            tau_sn(i) = 1d0/1.14d0*(m(i)+2d0*mp)/(2d0*mp)*1d0/(nH/2d0*1.69d-9)
         end if
         omega(i) = q(i)*B/(m(i)*c_cgs)
         sigma(i) = resistivite_chimie_x(i,iH,iT)*nH*(q(i))**2*tau_sn(i)/m(i)
      end do

      do  i=1,nbins_grains
         ! positively charged grains
         tau_sn(nion+1+3*(i-1)) = 1d0/1.28d0*(m_g(i)+2d0*mp)/(2d0*mp)*1d0/(nH/2d0*(pi*r_g(i)**2*(8d0*kB*T/(pi*2d0*mp))**0.5))
         omega( nion+1+3*(i-1)) = q(nion+1+3*(i-1))*B/(m_g(i)*c_cgs)
         sigma( nion+1+3*(i-1)) = resistivite_chimie_x(nion+1+3*(i-1),iH,iT)*nH*(q(nion+1+3*(i-1)))**2*tau_sn(nion+1+3*(i-1))/m_g(i)
         ! negatively charged grains
         tau_sn(nion+2+3*(i-1)) = tau_sn(nion+1+3*(i-1))
         omega( nion+2+3*(i-1)) = q(nion+2+3*(i-1))*B/(m_g(i)*c_cgs)
         sigma( nion+2+3*(i-1)) = resistivite_chimie_x(nion+2+3*(i-1),iH,iT)*nH*(q(nion+2+3*(i-1)))**2*tau_sn(nion+2+3*(i-1))/m_g(i)
         ! neutral grains
         tau_sn(nion+3+3*(i-1)) = tau_sn(nion+1+3*(i-1))
         omega( nion+3+3*(i-1)) = q(nion+3+3*(i-1))*B/(m_g(i)*c_cgs)
         sigma( nion+3+3*(i-1)) = resistivite_chimie_x(nion+3+3*(i-1),iH,iT)*nH*(q(nion+3+3*(i-1)))**2*tau_sn(nion+3+3*(i-1))/m_g(i)

      end do

      ! sum contributions from individual ions and grains
      sigP=0d0
      sigO=0d0
      sigH=0d0
      do i=1,nvarchimie
         sigP=sigP+sigma(i)
         sigO=sigO+sigma(i)/(1d0+(omega(i)*tau_sn(i))**2)
         sigH=sigH-sigma(i)*omega(i)*tau_sn(i)/(1d0+(omega(i)*tau_sn(i))**2)
      end do

      resistivite_chimie(1,iH,iT,iB)=log10(sigP)
      resistivite_chimie(2,iH,iT,iB)=log10(sigO)
      resistivite_chimie(3,iH,iT,iB)=log10(abs(sigH))
      resistivite_chimie(0,iH,iT,iB)=sign(1.0d0,sigH)

   end do
   end do
   end do
   
   deallocate(r_g,m_g,q,m)
   deallocate(omega,sigma,tau_sn)

   if(myid==1) write(*,*) '3D resistivities table complete'

end subroutine construct_resistivity_table
!###########################################################
!###########################################################
!###########################################################
subroutine interpolate_table(rho_cell,temp_cell,mag_cell,sigO,sigH,sigP)

   use amr_parameters,    only : dp
   use constants,         only : pi
   implicit none
   ! input: density, temperature, magnetic field strength and ionisation rate in the cell
   real(dp), intent(in)::rho_cell,temp_cell,mag_cell
   ! output: interpolated resistivity
   real(dp), intent(out)::sigO,sigH,sigP

   real(dp)::i_n,i_t,i_b
   integer::j,k,l
   real(dp)::j_dp,k_dp,l_dp,r_dp
   real(dp)::BBcgs,sigav
   real(dp), dimension(3)::x

   real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
   call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

   ! determine indices in table (as reals) 
   i_n=(1d0+(log10(rho_cell)-log10(nminchimie))/dnchimie)
   i_t=(1d0+(log10(temp_cell)-log10(tminchimie))/dtchimie)
   BBcgs=sqrt(mag_cell*(4d0*pi*scale_d*(scale_v)**2)) ! change to Gauss
   i_b=(1d0+(log10(BBcgs)-log10(bminchimie))/dbchimie)
   ! don't go outside of table
   i_n = max(1d0, min(i_n,real(nchimie)))
   i_t = max(1d0, min(i_t,real(tchimie)))
   i_b = max(1d0, min(i_b,real(bchimie)))

   ! convert to integers, lower bound
   j = min(floor(i_n),nchimie-1)
   k = min(floor(i_t),tchimie-1)
   l = min(floor(i_b),bchimie-1)
   j_dp = real(j,dp)
   k_dp = real(k,dp)
   l_dp = real(l,dp)

   x(1:3)=(1d0-(i_n-j_dp))*(1d0-(i_t-k_dp))*(1d0-(i_b-l_dp)) * (resistivite_chimie(1:3,j,  k  ,l))+&
             &((i_n-j_dp))*(1d0-(i_t-k_dp))*(1d0-(i_b-l_dp)) * (resistivite_chimie(1:3,j+1,k  ,l))+&
         &(1d0-(i_n-j_dp))*(    (i_t-k_dp))*(1d0-(i_b-l_dp)) * (resistivite_chimie(1:3,j,  k+1,l))+&
             &((i_n-j_dp))*(    (i_t-k_dp))*(1d0-(i_b-l_dp)) * (resistivite_chimie(1:3,j+1,k+1,l))+&
         &(1d0-(i_n-j_dp))*(1d0-(i_t-k_dp))*(    (i_b-l_dp)) * (resistivite_chimie(1:3,j,  k,  l+1))+&
             &((i_n-j_dp))*(1d0-(i_t-k_dp))*(    (i_b-l_dp)) * (resistivite_chimie(1:3,j+1,k,  l+1))+&
         &(1d0-(i_n-j_dp))*(    (i_t-k_dp))*(    (i_b-l_dp)) * (resistivite_chimie(1:3,j,  k+1,l+1))+&
             &((i_n-j_dp))*(    (i_t-k_dp))*(    (i_b-l_dp)) * (resistivite_chimie(1:3,j+1,k+1,l+1))
 
   sigP= 10d0**x(1)
   sigO= 10d0**x(2)

   ! modification since x(3) can be negative we simply use the sign of the leftmost
   ! point. If there is a sign inversion, we set it to zero.
   ! If you are using Hall resisitvities, this could be improved by using a linear
   ! interpolation instead of log.
   sigH=(10d0**x(3))*resistivite_chimie(0,j,k,l)
   sigav = sum(resistivite_chimie(0,j:j+1,k:k+1,l:l+1)) / 8d0
   if(sigav .ne. resistivite_chimie(0,j,k,l))then
      sigH = 0
   endif

end subroutine interpolate_table

end module resistivity_table