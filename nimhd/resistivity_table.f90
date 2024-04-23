module resistivity_table
   use amr_parameters, ONLY:dp
   use nimhd_commons
   use nimhd_parameters
   implicit none

   ! values at which the resistivities where evaluated
   !real(dp),allocatable,dimension(:)::table_values_b  ! magnetic field
   !real(dp),allocatable,dimension(:)::table_values_n  ! density
   !real(dp),allocatable,dimension(:)::table_values_t  ! temperature
   !real(dp),allocatable,dimension(:)::table_values_xi ! ionisation rate

   ! number of steps in the table for
   integer :: bchimie       ! B field (currently only used for reading table)
   integer :: nchimie       ! density 
   integer :: tchimie       ! temperature
   integer :: xichimie      ! ionisation rate

   ! step sizes in table variables for
   real(dp) :: dbchimie     ! B field
   real(dp) :: dnchimie     ! density 
   real(dp) :: dtchimie     ! temperature
   real(dp) :: dxichimie    ! ionisation rate

   ! minimum values for extrapolation variables in table for
   real(dp) :: bminchimie   ! B field
   real(dp) :: nminchimie   ! density  
   real(dp) :: tminchimie   ! temperature
   real(dp) :: ximinchimie  ! ionisation rate

   integer   :: nion=9            ! number of ions
   integer   :: nbins_grains      ! number of grain sizes
   !integer,parameter   :: Nvarchimie=7+(3*nbins_grains)   ! nombre d'especes ioniques
   integer   :: Nvarchimie   ! nombre d'especes ioniques
   !integer,parameter   :: Ntot=Nvar+6  ! nombre total d'especes considerees


   ! grain variables for each bin
   real(dp), allocatable, dimension(:) :: r_g    ! radius
   real(dp), allocatable, dimension(:) :: m_g    ! mass


   real(dp),allocatable,dimension(:,:,:,:)::resistivite_chimie_x ! to read in resistivity table


   !!!! Pour les collisions avec les grains
   real(dp), allocatable, dimension(:) :: q           ! charge 
   real(dp), allocatable, dimension(:) :: m           ! masse

   real(dp), parameter :: me=9.1094d-28 ! masse de l'electron, en g
   real(dp), parameter :: mp=1.6726d-24 ! masse du proton, en g
   real(dp), parameter :: e=4.803204d-10 ! charge de l'electron, en cgs


contains

subroutine read_resistivities
   use nimhd_parameters
   use nimhd_commons
   use constants, ONLY:mH

   integer::i,j,k
   real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
   call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
   !------------------------------------------
   ! Read resistivity tables for non-ideal MHD
   !------------------------------------------
   ! this should be refactored
   ! add parameter to choose how to compute the resistivity
   ! 1) fixed resistivity
   ! 2) analytical model resitivity(rho,T), Shu? Default option? Can be easily patched by the user
   ! 3) a table (like from Marchand 2016)
   ! make separate routine that computes the resistivity based on the chosen option
   ! move use_x3d to patch. We don't want to include a large table in ramses

   if(resistivity_table_ndim==4)then
      open(42,file=res_table_name,status='old')
      read(42,*) nchimie, tchimie, xichimie, nvarchimie
      read(42,*)
      read(42,*)
      allocate(resistivite_chimie_x(-2:nvarchimie+4,nchimie,tchimie,xichimie))

      do k=1,xichimie
         do i=1,tchimie
            do j=1,nchimie
               ! TC: resistivite_chimie_x(-2:0,....) contain the values of n, T and Xi?
               ! the rest contain the contributions to the resists from each ion and grain size?
               read(42,*)resistivite_chimie_x(-2:nvarchimie+4,j,i,k)
            end do
            read(42,*)
         end do
      end do
      close(42)

      rho_threshold=max(rho_threshold,resistivite_chimie_x(-2,1,1,1)*(mu_gas*mH)/scale_d) ! input in part/cc, output in code units
      nminchimie=(resistivite_chimie_x(-2,1,1,1))
      dnchimie=(log10(resistivite_chimie_x(-2,nchimie,1,1))-log10(resistivite_chimie_x(-2,1,1,1)))/&
               &(nchimie-1)
      tminchimie=(resistivite_chimie_x(-1,1,1,1))
      dtchimie=(log10(resistivite_chimie_x(-1,1,tchimie,1))-log10(resistivite_chimie_x(-1,1,1,1)))/&
               &(tchimie-1)
      ximinchimie=(resistivite_chimie_x(0,1,1,1))
      dxichimie=(log10(resistivite_chimie_x(0,1,1,xichimie))-log10(resistivite_chimie_x(0,1,1,1)))/&
               &(xichimie-1)
      call rq_3d
      call nimhd_4dtable
   else
      print*, 'must choose an input for abundances or resistivities'
      stop
   endif

end subroutine read_resistivities

!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine nimhd_3dtable
   use nimhd_commons
   use amr_commons, only : myid
   use constants, only:pi,c_cgs, kB
   implicit none

   integer  :: iB,iH,iT,i
   real(dp) :: B,bmaxchimie,nH,T,sigH,sigO,sigP
   real(dp), dimension(nvarchimie) :: x

   real(dp), allocatable, dimension(:)  :: sigma         ! sigma_s
   real(dp), allocatable, dimension(:)  :: tau_sn
   real(dp), allocatable, dimension(:)  :: omega

   if(myid==1) write(*,*) 'Computing 3D resistivities table'
   
   ! values for Btable
   bminchimie=1d-10
   bmaxchimie=1d5               ! ok for first core in nimhd. maybe not enough for second core.
   bchimie=100
   dbchimie=(log10(bmaxchimie)-log10(bminchimie))/real((bchimie-1),dp)

   allocate(resistivite_chimie(0:3,1:nchimie,1:tchimie,1,1:bchimie))
   allocate(sigma(nvarchimie))
   allocate(tau_sn(nvarchimie))
   allocate(omega(nvarchimie))

   tau_sn      = 0.0_dp
   omega       = 0.0_dp
   sigma       = 0.0_dp
      
   do  iB=1,bchimie
   do  iT=1,tchimie
   do  iH=1,nchimie

      nh=resistivite_chimie_x(0,iH,iT,1)  ! density (.cc) of current point
      B =10.0d0**(log10(bminchimie)+dble(iB-1)*dbchimie)
      T =resistivite_chimie_x(-1,iH,iT,1)
      
      !inp=nh/2d0/H2_fraction     ! inp is neutrals.cc, to fit densionbis
      !print *, q(1),m(1)
      !stop
      do  i=1,nvarchimie-3*nbins_grains
         if  (i==1) then  ! electron
            x(i) = resistivite_chimie_x(i,iH,iT,1)
            tau_sn(i) = 1d0/1.16d0*(m(i)+2d0*mp)/(2d0*mp)*1d0/(nH/2d0*1.3d-9)
            omega(i) = q(i)*B/(m(i)*c_cgs)
            sigma(i) = x(i)*nH*(q(i))**2*tau_sn(i)/m(i)
         else if (i>=2 .and. i<=7) then ! ions
            x(i) = resistivite_chimie_x(i,iH,iT,1)
            tau_sn(i) = 1d0/1.14d0*(m(i)+2d0*mp)/(2d0*mp)*1d0/(nH/2d0*1.69d-9)
            omega(i) = q(i)*B/(m(i)*c_cgs)
            sigma(i) = x(i)*nH*(q(i))**2*tau_sn(i)/m(i)
         end if
      end do
      do  i=1,nbins_grains   ! grains
         x(8+3*(i-1)) = resistivite_chimie_x(8+3*(i-1),iH,iT,1)
         tau_sn(8+3*(i-1))=1d0/1.28d0*(m_g(i)+2d0*mp)/(2d0*mp)*1d0/(nH/2d0*(pi*r_g(i)**2*(8d0*kB*T/(pi*2d0*mp))**0.5))    ! g+
         omega(8+3*(i-1)) = q(8+3*(i-1))*B/(m_g(i)*c_cgs)
         sigma(8+3*(i-1)) = x(8+3*(i-1))*nH*(q(8+3*(i-1)))**2*tau_sn(8+3*(i-1))/m_g(i)

         x(9+3*(i-1)) = resistivite_chimie_x(9+3*(i-1),iH,iT,1)
         tau_sn(9+3*(i-1))=tau_sn(8+3*(i-1))          ! g-
         omega(9+3*(i-1)) = q(9+3*(i-1))*B/(m_g(i)*c_cgs)
         sigma(9+3*(i-1)) = x(9+3*(i-1))*nH*(q(9+3*(i-1)))**2*tau_sn(9+3*(i-1))/m_g(i)

         x(10+3*(i-1)) = resistivite_chimie_x(10+3*(i-1),iH,iT,1)
         tau_sn(10+3*(i-1))=tau_sn(8+3*(i-1))          ! g0
         omega(10+3*(i-1)) = q(10+3*(i-1))*B/(m_g(i)*c_cgs)
         sigma(10+3*(i-1)) = x(10+3*(i-1))*nH*(q(10+3*(i-1)))**2*tau_sn(10+3*(i-1))/m_g(i)

      end do

      sigP=0d0
      sigO=0d0
      sigH=0d0

      do i=1,nvarchimie
         sigP=sigP+sigma(i)
         sigO=sigO+sigma(i)/(1d0+(omega(i)*tau_sn(i))**2)
         sigH=sigH-sigma(i)*omega(i)*tau_sn(i)/(1d0+(omega(i)*tau_sn(i))**2)
      end do
      resistivite_chimie(1,iH,iT,1,iB)=log10(sigP)
      resistivite_chimie(2,iH,iT,1,iB)=log10(sigO)
      resistivite_chimie(3,iH,iT,1,iB)=log10(abs(sigH))
      resistivite_chimie(0,iH,iT,1,iB)=sign(1.0d0,sigH)
   end do
   end do
   end do
   
   deallocate(r_g,m_g,q,m,sigma,tau_sn)
   deallocate(omega)
   deallocate(resistivite_chimie_x)

   if(myid==1) write(*,*) '3D resistivities table complete'

end subroutine nimhd_3dtable
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine rq_3d    ! just use it once if you change the grains distribution to update the variables q(:) and r_g(:)
   use nimhd_commons
   use constants, only:pi
   implicit none
   integer :: i
   real(dp):: Lp1,Lp3,Lp4,nbins_real

   ! grain density in g/cc
   real(dp), parameter :: rho_s=2.3_dp

   ! grain size distribution parameters, cf Kunz & Mouschovias 2009
   real(dp), parameter :: a_0=0.0375d-4      ! cm
   real(dp), parameter :: a_min=0.0181d-4    ! cm
   real(dp), parameter :: a_max=0.9049d-4    ! cm
   real(dp), parameter :: zeta=a_min/a_max   ! a_min/a_max
   real(dp), parameter :: lambda_pow=-3.5d0  ! Coeff power law

   allocate(r_g(nbins_grains)) ! grain sizes
   allocate(m_g(nbins_grains)) ! grain masses
   allocate(m(nvarchimie))     ! particle (ions + grains) masses
   allocate(q(nvarchimie))     ! particle (ions + grains) charges

   ! determine grain sizes

   ! table contains values for neutral, positive and negatively charged grains for each size bin
   ! so we divide by 3 to get the actual amount of size bins
   nbins_real=real(nvarchimie-nion,dp)/3.0_dp
   nbins_grains=floor(nbins_real)
   if (nbins_real.ne.real(nbins_grains,dp)) then
      print*, 'issue in number of species'
      stop
   endif
   
   Lp1=dble(lambda_pow+1)
   Lp3=dble(lambda_pow+3)
   Lp4=dble(lambda_pow+4)

   if(nbins_grains==1) then
     ! if we have only 1 bin, we take this average value
     r_g(1)=a_0
   else
     do  i=1,nbins_grains    ! cf Kunz & Mouschovias 2009
       r_g(i)=a_min*zeta**(-dble(i)/dble(nbins_grains)) * &
            & dsqrt( Lp1/Lp3* (1d0-zeta**(Lp3/dble(nbins_grains)))/(1d0-zeta**(Lp1/dble(nbins_grains))))
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
   m(8) = 39.098*mp        ! K+
   m(9) = 22.99d0*mp       ! Na+

   do i=1,nbins_grains     ! masse des grains
      m_g(i)=4d0/3d0*pi*r_g(i)**3*rho_s
      m(nion+3*(i-1)+1:nion+3*i)=m_g(i) !mass for neutral, postive and negative grain is the same
   end do

end subroutine rq_3d
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine nimhd_4dtable
   use nimhd_commons
   use amr_commons, only : myid
   use constants, only:pi,c_cgs,kB
   implicit none

   integer  :: iB,iH,iT,iX,i
   real(dp) :: B,bmaxchimie,nH,T,xi,sigH,sigO,sigP
   real(dp) :: sigv, muuu


   ! resistivites (cf Kunz & Mouschovias 2009)
   real(dp), allocatable, dimension(:)  :: sigma         ! sigma_s
   real(dp), allocatable, dimension(:)  :: tau_sn
   real(dp), allocatable, dimension(:)  :: omega

   if(myid==1) write(*,*) 'Computing 3D resistivities table'
   
   ! values for Btable
   bminchimie=1d-10
   bmaxchimie=1d10               ! ok for first core in nimhd. maybe not enough for second core.
   bchimie=150
   dbchimie=(log10(bmaxchimie)-log10(bminchimie))/real((bchimie-1),dp)

   allocate(resistivite_chimie(0:3,1:nchimie,1:tchimie,1:xichimie,1:bchimie))
   allocate(sigma(nvarchimie))
   allocate(tau_sn(nvarchimie))
   allocate(omega(nvarchimie))

   tau_sn      = 0.0_dp
   omega       = 0.0_dp
   sigma       = 0.0_dp
      
   do  iB=1,bchimie
   do  iX=1,xichimie
   do  iT=1,tchimie
   do  iH=1,nchimie

      nh=resistivite_chimie_x(-2,iH,iT,iX)  ! density (.cc) of current point
      B =10.0d0**(log10(bminchimie)+dble(iB-1)*dbchimie)
      T =resistivite_chimie_x(-1,iH,iT,iX)
      xi =resistivite_chimie_x(0,iH,iT,iX)
      
!      write(*,*) ih,it,ib,nh,b,t
!      read(*,*)
      !inp=nh/2d0/H2_fraction     ! inp is neutrals.cc, to fit densionbis
      !print *, q(1),m(1)
      !stop

      ! TC: I don't get why you need to read a table if you are still going to throw a bunch of additional formulas at it

      do i=1,nion
        if  (i==1) then 
          ! electron
          sigv=3.16d-11 * (dsqrt(8d0*kB*1d-7*T/(pi*me*1d-3))*1d-3)**1.3d0
          tau_sn(i) = 1d0/1.16d0*(m(i)+2d0*mp)/(2d0*mp)*1d0/(nH/2d0*sigv)
        else
          ! ions
          muuu=m(i)*2d0*mp/(m(i)+2d0*mp)
          if (i==2 .or. i==3) then
            sigv=2.4d-9 *(dsqrt(8d0*kB*1d-7*T/(pi*muuu*1d-3))*1d-3)**0.6d0
          else if (i==4) then
            sigv=2d-9 * (dsqrt(8d0*kB*1d-7*T/(pi*muuu*1d-3))*1d-3)**0.15d0
          else if (i==5) then
            sigv=3.89d-9 * (dsqrt(8d0*kB*1d-7*T/(pi*muuu*1d-3))*1d-3)**(-0.02d0)
          else
            sigv=1.69d-9
          end if
          tau_sn(i) = 1d0/1.14d0*(m(i)+2d0*mp)/(2d0*mp)*1d0/(nH/2d0*sigv)
        end if
        omega(i) = q(i)*B/(m(i)*c_cgs)
        sigma(i) = resistivite_chimie_x(i,iH,iT,iX)*nH*(q(i))**2*tau_sn(i)/m(i)
      end do
      
      do  i=1,nbins_grains
        ! grains TC: neutral? Shouldn't this be charged? In any case, something is missing
        tau_sn(nion+1+3*(i-1))= 1d0/1.28d0*(m_g(i)+2d0*mp)/(2d0*mp)*1d0/(nH/2d0*(pi*r_g(i)**2*(8d0*kB*T/(pi*2d0*mp))**0.5))
        omega(nion+1+3*(i-1)) = q(nion+1+3*(i-1))*B/(m_g(i)*c_cgs)
        sigma(nion+1+3*(i-1)) = resistivite_chimie_x(nion+1+3*(i-1),iH,iT,iX)*nH*(q(nion+1+3*(i-1)))**2*tau_sn(nion+1+3*(i-1))/m_g(i)
      
        tau_sn(nion+2+3*(i-1))= tau_sn(nion+1+3*(i-1))
        omega(nion+2+3*(i-1)) = q(nion+2+3*(i-1))*B/(m_g(i)*c_cgs)
        sigma(nion+2+3*(i-1)) = resistivite_chimie_x(nion+2+3*(i-1),iH,iT,iX)*nH*(q(nion+2+3*(i-1)))**2*tau_sn(nion+2+3*(i-1))/m_g(i)
      
      end do

      ! TC: sum contributions from individual ions and grains
      sigP=0d0
      sigO=0d0
      sigH=0d0
      do i=1,nvarchimie
         sigP=sigP+sigma(i)
         sigO=sigO+sigma(i)/(1d0+(omega(i)*tau_sn(i))**2)
         sigH=sigH-sigma(i)*omega(i)*tau_sn(i)/(1d0+(omega(i)*tau_sn(i))**2)
      end do

      if(sigH==0d0) sigH=1d-30

      ! TC: what is 0, 1, 2 and 3?
      resistivite_chimie(0,iH,iT,iX,iB)=sign(1.0d0,sigH)
      resistivite_chimie(1,iH,iT,iX,iB)=log10(sigP)
      resistivite_chimie(2,iH,iT,iX,iB)=log10(sigO)
      resistivite_chimie(3,iH,iT,iX,iB)=log10(abs(sigH))

   end do
   end do
   end do
   end do

   deallocate(r_g,m_g,q,m)
   deallocate(omega,sigma,tau_sn)
   deallocate(resistivite_chimie_x)
   
   if(myid==1) write(*,*) '3D resistivities table complete'

end subroutine nimhd_4dtable

!###########################################################
!###########################################################
!###########################################################
!subroutine sig_x2d(ll,ii,j,k,lb,ib,sigO,sigH,sigP)
subroutine interpolate_table(rho_cell,temp_cell,mag_cell,ionisrate_cell,sigO,sigH,sigP)

   use amr_parameters,    only : dp
   use constants,         only : pi
   implicit none
   ! input: density, temperature, magnetic field strength and ionisation rate in the cell
   real(dp), intent(in)::rho_cell,temp_cell,mag_cell,ionisrate_cell
   ! output: interpolated resistivity
   real(dp), intent(out)::sigO,sigH,sigP

   real(dp)::i_n,i_t,i_b,i_r
   integer::j,k,l,r
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

   if(resistivity_table_ndim==3)then
      r = 1
      x(1:3)=(1d0-(i_n-j_dp))*(1d0-(i_t-k_dp))*(1d0-(i_b-l_dp)) * (resistivite_chimie(1:3,j,  k  ,r,l))+&
                &((i_n-j_dp))*(1d0-(i_t-k_dp))*(1d0-(i_b-l_dp)) * (resistivite_chimie(1:3,j+1,k  ,r,l))+&
            &(1d0-(i_n-j_dp))*(    (i_t-k_dp))*(1d0-(i_b-l_dp)) * (resistivite_chimie(1:3,j,  k+1,r,l))+&
                &((i_n-j_dp))*(    (i_t-k_dp))*(1d0-(i_b-l_dp)) * (resistivite_chimie(1:3,j+1,k+1,r,l))+&
            &(1d0-(i_n-j_dp))*(1d0-(i_t-k_dp))*(    (i_b-l_dp)) * (resistivite_chimie(1:3,j,  k,  r,l+1))+&
                &((i_n-j_dp))*(1d0-(i_t-k_dp))*(    (i_b-l_dp)) * (resistivite_chimie(1:3,j+1,k,  r,l+1))+&
            &(1d0-(i_n-j_dp))*(    (i_t-k_dp))*(    (i_b-l_dp)) * (resistivite_chimie(1:3,j,  k+1,r,l+1))+&
                &((i_n-j_dp))*(    (i_t-k_dp))*(    (i_b-l_dp)) * (resistivite_chimie(1:3,j+1,k+1,r,l+1))
 
      sigav = sum(resistivite_chimie(0,j:j+1,k:k+1,r,l:l+1)) / 8d0

   else if(resistivity_table_ndim==4)then
      i_r=(1d0+(log10(ionisrate_cell)-log10(ximinchimie))/dxichimie)
      i_r = max(1d0, min(i_r,real(xichimie)))
      r = min(floor(i_r),xichimie-1)
      r_dp = real(r,dp)
   
      x(1:3)=(1d0-(i_n-j_dp))*(1d0-(i_t-k_dp))*(1d0-(i_r-r_dp))*(1d0-(i_b-l_dp))*(resistivite_chimie(1:3,j,  k,  r,  l))+&
                &((i_n-j_dp))*(1d0-(i_t-k_dp))*(1d0-(i_r-r_dp))*(1d0-(i_b-l_dp))*(resistivite_chimie(1:3,j+1,k,  r,  l))+&
            &(1d0-(i_n-j_dp))*(    (i_t-k_dp))*(1d0-(i_r-r_dp))*(1d0-(i_b-l_dp))*(resistivite_chimie(1:3,j,  k+1,r,  l))+&
                &((i_n-j_dp))*(    (i_t-k_dp))*(1d0-(i_r-r_dp))*(1d0-(i_b-l_dp))*(resistivite_chimie(1:3,j+1,k+1,r,  l))+&
            &(1d0-(i_n-j_dp))*(1d0-(i_t-k_dp))*(    (i_r-r_dp))*(1d0-(i_b-l_dp))*(resistivite_chimie(1:3,j,  k,  r+1,l))+&
                &((i_n-j_dp))*(1d0-(i_t-k_dp))*(    (i_r-r_dp))*(1d0-(i_b-l_dp))*(resistivite_chimie(1:3,j+1,k,  r+1,l))+&
            &(1d0-(i_n-j_dp))*(    (i_t-k_dp))*(    (i_r-r_dp))*(1d0-(i_b-l_dp))*(resistivite_chimie(1:3,j,  k+1,r+1,l))+&
                &((i_n-j_dp))*(    (i_t-k_dp))*(    (i_r-r_dp))*(1d0-(i_b-l_dp))*(resistivite_chimie(1:3,j+1,k+1,r+1,l))+&
            &(1d0-(i_n-j_dp))*(1d0-(i_t-k_dp))*(1d0-(i_r-r_dp))*(    (i_b-l_dp))*(resistivite_chimie(1:3,j,  k,  r,  l+1))+&
                &((i_n-j_dp))*(1d0-(i_t-k_dp))*(1d0-(i_r-r_dp))*(    (i_b-l_dp))*(resistivite_chimie(1:3,j+1,k,  r,  l+1))+&
            &(1d0-(i_n-j_dp))*(    (i_t-k_dp))*(1d0-(i_r-r_dp))*(    (i_b-l_dp))*(resistivite_chimie(1:3,j,  k+1,r,  l+1))+&
                &((i_n-j_dp))*(    (i_t-k_dp))*(1d0-(i_r-r_dp))*(    (i_b-l_dp))*(resistivite_chimie(1:3,j+1,k+1,r,  l+1))+&
            &(1d0-(i_n-j_dp))*(1d0-(i_t-k_dp))*(    (i_r-r_dp))*(    (i_b-l_dp))*(resistivite_chimie(1:3,j,  k,  r+1,l+1))+&
                &((i_n-j_dp))*(1d0-(i_t-k_dp))*(    (i_r-r_dp))*(    (i_b-l_dp))*(resistivite_chimie(1:3,j+1,k,  r+1,l+1))+&
            &(1d0-(i_n-j_dp))*(    (i_t-k_dp))*(    (i_r-r_dp))*(    (i_b-l_dp))*(resistivite_chimie(1:3,j,  k+1,r+1,l+1))+&
                &((i_n-j_dp))*(    (i_t-k_dp))*(    (i_r-r_dp))*(    (i_b-l_dp))*(resistivite_chimie(1:3,j+1,k+1,r+1,l+1))

      sigav = sum(resistivite_chimie(0,j:j+1,k:k+1,r:r+1,l:l+1)) / 16d0

   endif

   sigP= 10d0**x(1)
   sigO= 10d0**x(2)

   ! modification since x(3) can be negative we simply use the sign of the leftmost
   ! point. If there is a sign inversion, we set it to zero.
   ! If you are using Hall resisitvities, this could be improved by using a linear
   ! interpolation instead of log.
   sigH=(10d0**x(3))*resistivite_chimie(0,j,k,r,l)
   if(sigav .ne. resistivite_chimie(0,j,k,r,l))then
      sigH = 0
   endif

end subroutine interpolate_table

end module resistivity_table