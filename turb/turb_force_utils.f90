#if USE_TURB==1
! TURB_FORCE_UTILS.F90
! A.McLeod - 24/11/2012
! Utility routines for turbulent forcing
! ============================================================================

#if NDIM==3
#define HERMITIAN_FIELD 1
#endif

#ifdef NPRE
#if !(NPRE==4)
#define DOUBLE_PRECISION 1
#endif
#endif

#define PI 3.141592653589793_dp

subroutine find_conj_pair(i,j,k,ii,jj,kk)
   use turb_commons
   implicit none
   integer, intent(in)     :: i,j,k         ! Grid coordinates
   integer, intent(out)    :: ii,jj,kk      ! Hermitian pair

   ii = TURB_GS - i
   if (ii >= TURB_GS) ii = ii - TURB_GS
   jj = TURB_GS - j
   if (jj >= TURB_GS) jj = jj - TURB_GS
   kk = TURB_GS - k
   if (kk >= TURB_GS) kk = kk - TURB_GS

end subroutine find_conj_pair

subroutine find_unitk(i,j,k,limit,unitk)
   use turb_commons
   implicit none
   integer, intent(in)         :: i,j,k,limit
   integer                     :: ii,jj,kk
   real (kind=dp), intent(out) :: unitk(1:NDIM)

   ii = i
   if (i > limit) ii = TURB_GS - i
#if NDIM>1
   jj = j
   if (j > limit) jj = TURB_GS - j
#endif
#if NDIM>2
   kk = k
   if (k > limit) kk = TURB_GS - k
#endif

#if NDIM==1
   unitk = (/real(ii,dp)/)
#elif NDIM==2
   unitk = (/real(ii,dp),real(jj,dp)/)
#else
   unitk = (/real(ii,dp),real(jj,dp),real(kk,dp)/)
#endif
#if NDIM>1
   unitk = unitk / sqrt(sum(unitk**2))
#endif

end subroutine find_unitk

subroutine calc_power_spectrum(k, power_spectrum)
   use turb_commons
   implicit none
   integer, intent(in)        :: k(1:3)         ! Wavevector
   real(kind=dp), intent(out) :: power_spectrum ! Power value
   real(kind=dp)              :: k_mag          ! Wavevector magnitude

   ! Remark that the components of k are between -TURB_GS and TURB_GS
   ! with k=1 corresponding to the box size
   select case(forcing_power_spectrum)
      case('power_law')
         ! alpha^-2 power spectrum
         if (all(k==0)) then
            power_spectrum = 0
            return
         end if

         k_mag = sqrt(real(sum(k**2),dp))
         if (k_mag > (TURB_GS/2)) then
            power_spectrum = 0
            return
         end if
         power_spectrum = k_mag**(-2)

      case('parabolic')
         ! 'parabola' large-scale modes power spectrum
         power_spectrum = 0._dp
         k_mag = sqrt(real(sum(k**2),dp))
         if ((k_mag > 1.0_dp) .AND. (k_mag < 3.0_dp)) then
             power_spectrum = 1.0 - (k_mag-2.0)**2
         end if

      case('konstandin')
         ! forcing between k=1 (max) and k=1 (zero) as in Konstandin 2015
         power_spectrum = 0._dp
         k_mag = sqrt(real(sum(k**2),dp))
         if ((k_mag >= 0.999999999999999_dp) .AND. (k_mag < 2.0_dp)) then
             power_spectrum = 2.0 - (k_mag)
         end if

      case('test')
         power_spectrum = 0._dp
         if (k(1)==1 .AND. k(2)==0 .AND. k(3)==0) then
            power_spectrum = 1.0
         end if

      case default
         write (6,*) "Unknown forcing_power_spectrum!"
         write (6,*) "Use 'power_law', 'parabolic', 'konstandin' or 'test'"
         stop
      end select

end subroutine calc_power_spectrum

subroutine gaussian_cmplx(G)
   !use constants, only:pi
   use turb_commons
   implicit none
   ! Generates 3 complex random variates, where each variate has a
   ! complex value drawn from a Gaussian distribution EITHER
   ! by assigning a Gaussian magnitude and uniformly distributed argument OR
   ! by assigning Gaussian real and imaginary components
   complex(kind=cdp), intent(out) :: G(1:ndim) ! Random gaussian values
   integer              :: d           ! Dimension counter
   real(kind=dp)        :: Rnd(1:3)    ! Random numbers
   real(kind=dp)        :: w           ! Aux. variable for random Gaussian
   real(kind=dp)        :: mag(1:ndim), arg(1:ndim) ! Magnitude and argument
   ! Create Gaussian distributed random numbers
   ! Marsaglia version of the Box-Muller algorithm
   ! This creates a Gaussian with a spread of 1
   ! and centred on 0 (hopefully)
   ! You get two independent variates each time

   ! Method 1: Gaussian magnitude and uniformly-random argument
   do d=1,ndim
      do
         call kiss64_double(2, kiss64_state, Rnd(1:2))
         Rnd(1) = 2.0 * Rnd(1) - 1.0
         Rnd(2) = 2.0 * Rnd(2) - 1.0
         w = Rnd(1)**2 + Rnd(2)**2
         if (w<1.0) exit
      end do
      w = sqrt( (-2.0 * log( w ) ) / w )

      mag(d) = Rnd(1) * w
      ! Throwing away second random variate (something of a waste)
   end do

   call kiss64_double(ndim, kiss64_state, Rnd(1:ndim))
   arg = (Rnd(1:ndim) * 2.0 * PI) - 1.0

   G = cmplx(mag * cos(arg), mag * sin(arg), kind=cdp)

   ! Method 2: Gaussian real and imaginary components
!    do d=1,ndim
!       do
!          call kiss64_double(ndim, kiss64_state, Rnd(1:ndim))
!          Rnd(1) = 2.0 * Rnd(1) - 1.0
!          Rnd(2) = 2.0 * Rnd(2) - 1.0
!          w = Rnd(1)**2 + Rnd(2)**2
!          if (w<1.0) exit
!       end do
!       w = sqrt( (-2.0 * log( w ) ) / w )
!
!       G(d) = cmplx(Rnd(1) * w, Rnd(2) * w)
!    end do

end subroutine gaussian_cmplx

subroutine add_turbulence(turb_field, dt)
   use turb_commons
   implicit none
   ! Take a complex Fourier field of turbulence, and add a random component
   ! from a Wiener process (in this case Gaussian perturbations) multiplied
   ! by an initial distribution of power and then projected by Helmholtz
   ! decomposition.
   ! Optionally, this field can be Hermitian i.e. purely real after the
   ! transform, but I'm not sure if this produces a purely even field...
   complex(kind=cdp), intent(inout) :: turb_field(1:NDIM, 0:TGRID_X,&
                                                 &0:TGRID_Y, 0:TGRID_Z)
                                           ! Complex field to add to
   real(kind=dp), intent(in)       :: dt   ! Width of Gaussian which drives
                                           ! Wiener process

   integer              :: i, j, k         ! Loop variables
   integer              :: halfway         ! Half of grid size

#if defined(HERMITIAN_FIELD)
   integer              :: ii,jj,kk        ! Hermitian pair
   logical              :: hermitian_pair  ! Do we need to take conjugate?
   logical              :: own_conjg       ! Are we are own conjugate pair?
#endif

   complex(kind=cdp)    :: complex_vec(1:NDIM) ! Complex Fourier vector
#if NDIM>1
   real(kind=dp)        :: unitk(1:NDIM)       ! Unit vector parallel to k
   complex(kind=cdp)    :: unitk_cmplx(1:NDIM) ! Unit vector parallel to k
   complex(kind=cdp)    :: comp_cmplx(1:NDIM)  ! Compressive component
   complex(kind=cdp)    :: sol_cmplx(1:NDIM)   ! Solenoidal components
#endif

   halfway = TURB_GS / 2 ! integer division

   do k=0,TGRID_Z
      do j=0,TGRID_Y
         do i=0,TGRID_X
#if defined(HERMITIAN_FIELD)
            hermitian_pair = .FALSE.
#endif
            ! Check there is any power in this mode, else cycle
            if (power_spec(i,j,k) == 0.0_dp) cycle

#if defined(HERMITIAN_FIELD)
            ! Test if we are own conjugate or need to use another conjugate
            own_conjg = .FALSE.
            call find_conj_pair(i,j,k,ii,jj,kk)
            if (i==ii .AND. j==jj .AND. k==kk) own_conjg = .TRUE.
            if (k > halfway) then
               hermitian_pair = .TRUE.
            else if (k == 0 .OR. 2*k == TURB_GS) then
               if (j > halfway) then
                  hermitian_pair = .TRUE.
               else if (j == 0 .OR. 2*j == TURB_GS) then
                  if (i > halfway) hermitian_pair = .TRUE.
               end if
            end if
            if (hermitian_pair) then
               turb_field(1:3,i,j,k) = &
                  & conjg(turb_field(1:3,ii,jj,kk))
               cycle
            end if
#endif

            ! Ornstein-Uhlenbeck process
            ! dF(k,t) = F_0(k) P(k) dW - F(k,t) dt / T
            ! where F(k,t) is the vector fourier amplitude,
            ! F_0 is the power spectrum distribution of amplitudes,
            ! P(k) is the projection (Helmholtz decomposition),
            ! dt is the timestep of integration,
            ! and dW is the Wiener process, described by a random variate
            ! from a normal distribution with a mean of zero and a
            ! variance of dt i.e. a standard deviation of sqrt(dt)

            ! In this we only calculate the first term - the random
            ! growth term.

            ! Random power/phase (complex Gaussian) multiplied by power
            ! spectrum and multiplied by width dt
            call gaussian_cmplx(complex_vec)
            complex_vec = complex_vec * sqrt(dt) * power_spec(i,j,k)

#if defined(HERMITIAN_FIELD)
            if (own_conjg) then
               ! To be own conjugate, phase must be zero
               ! (in order to be real)
               complex_vec = abs(complex_vec)
            end if
#endif

#if NDIM>1
            ! Helmholtz decomposition!
            call find_unitk(i, j, k, halfway, unitk)
            unitk_cmplx = cmplx(unitk, kind=cdp)
            comp_cmplx = unitk_cmplx *  dot_product(unitk_cmplx,complex_vec)
            sol_cmplx = complex_vec - comp_cmplx
            complex_vec = comp_cmplx*comp_frac + sol_cmplx*sol_frac
            ! note that there are two degrees of freedom for
            ! solenoidal/transverse modes, and only one for
            ! longitudinal/compressive modes,
            ! and so for comp_frac == sol_frac == 0.5, you get
            ! F(solenoidal) == 2 * F(compressive)
            ! as per Federrath et al 2010
            turb_field(1:NDIM,i,j,k) = turb_field(1:NDIM,i,j,k) + &
               & complex_vec
#else
            turb_field(1:NDIM,i,j,k) = turb_field(1:NDIM,i,j,k) + &
               & complex_vec
#endif

         end do
      end do
   end do

end subroutine add_turbulence

subroutine decay_turbulence(old_turb_field, new_turb_field)
   use turb_commons
   implicit none
   ! Take a old complex Fourier field of turbulence, and find
   ! decay. Remove this decay from new turbulent field.
   complex(kind=cdp), intent(in)    :: old_turb_field(1:NDIM, 0:TGRID_X,&
                                                     &0:TGRID_Y, 0:TGRID_Z)
                                           ! Old field for reference
   complex(kind=cdp), intent(inout) :: new_turb_field(1:NDIM, 0:TGRID_X,&
                                                     &0:TGRID_Y, 0:TGRID_Z)
                                           ! Complex field to subtract from

  ! Ornstein-Uhlenbeck process
  ! dF(k,t) = F_0(k) P(k) dW - F(k,t) dt / T
  ! where F(k,t) is the vector fourier amplitude,
  ! F_0 is the power spectrum distribution of amplitudes,
  ! P(k) is the projection (Helmholtz decomposition),
  ! dt is the timestep of integration,
  ! and dW is the Wiener process, described by a random variate
  ! from a normal distribution with a mean of zero and a
  ! variance of dt i.e. a standard deviation of sqrt(dt)

  ! In this we only calculate the second term - the exponential
  ! decay term.

   new_turb_field = new_turb_field - (old_turb_field * turb_decay_frac)

end subroutine decay_turbulence

subroutine FFT_1D(complex_field, real_field)
   use turb_commons
   implicit none
#include "fftw3.f"
   ! Transform complex field into purely real field for 1D vector field

   integer (kind=ILP)   :: plan            ! FFTW plan pointer
   complex(kind=cdp), intent(in)  :: complex_field(0:TGRID_X)
                                           ! Complex field to transform
   real(kind=dp), intent(out)     :: real_field(0:TGRID_X)
                                           ! Result of transforms

   complex(kind=cdp), allocatable :: fftfield(:) ! Memory for FFT

   ! Allocate storage for performing FFTs
   allocate(fftfield(0:TGRID_X))

#if defined(DOUBLE_PRECISION)
   call dfftw_plan_dft_1d(plan, TURB_GS, fftfield, &
      & fftfield, FFTW_BACKWARD, FFTW_ESTIMATE)
#else
   call sfftw_plan_dft_1d(plan, TURB_GS, fftfield, &
      & fftfield, FFTW_BACKWARD, FFTW_ESTIMATE)
#endif

   fftfield = complex_field(:)

#if defined(DOUBLE_PRECISION)
   call dfftw_execute_dft(plan, fftfield, fftfield)
#else
   call sfftw_execute_dft(plan, fftfield, fftfield)
#endif
   real_field(:) = real(fftfield, dp) / (turb_gs_real)

#if defined(DOUBLE_PRECISION)
   call dfftw_destroy_plan(plan)
#else
   call sfftw_destroy_plan(plan)
#endif

   deallocate(fftfield)

end subroutine FFT_1D

subroutine FFT_2D(complex_field, real_field)
   use turb_commons
   implicit none
#include "fftw3.f"
   ! Transform complex field into purely real field for 2D vector field

   integer              :: d               ! Dimension counter
   integer (kind=ILP)   :: plan            ! FFTW plan pointer
   complex(kind=cdp), intent(in)  :: complex_field(1:2, 0:TGRID_X,&
                                                  &0:TGRID_Y)
                                           ! Complex field to transform
   real(kind=dp), intent(out)     :: real_field(1:2, 0:TGRID_X,&
                                               &0:TGRID_Y)
                                           ! Result of transforms

   complex(kind=cdp), allocatable :: fftfield(:,:) ! Memory for FFT

   ! Allocate storage for performing FFTs
   allocate(fftfield(0:TGRID_X,0:TGRID_Y))

#if defined(DOUBLE_PRECISION)
   call dfftw_plan_dft_2d(plan, TURB_GS, TURB_GS, fftfield, &
      & fftfield, FFTW_BACKWARD, FFTW_ESTIMATE)
#else
   call sfftw_plan_dft_2d(plan, TURB_GS, TURB_GS, fftfield, &
      & fftfield, FFTW_BACKWARD, FFTW_ESTIMATE)
#endif

   do d=1,2
      fftfield = complex_field(d,:,:)
#if defined(DOUBLE_PRECISION)
      call dfftw_execute_dft(plan, fftfield, fftfield)
#else
      call sfftw_execute_dft(plan, fftfield, fftfield)
#endif
      real_field(d,:,:) = real(fftfield, dp) / (turb_gs_real**2)
   end do

#if defined(DOUBLE_PRECISION)
   call dfftw_destroy_plan(plan)
#else
   call sfftw_destroy_plan(plan)
#endif

   deallocate(fftfield)

end subroutine FFT_2D

subroutine FFT_3D(complex_field, real_field)
   use turb_commons
   implicit none
#include "fftw3.f"
   ! Transform complex field into purely real field for 3D vector field

   integer              :: d               ! Dimension counter
   integer (kind=ILP)   :: plan            ! FFTW plan pointer
   complex(kind=cdp), intent(in)  :: complex_field(1:3, 0:TGRID_X,&
                                                  &0:TGRID_Y, 0:TGRID_Z)
                                           ! Complex field to transform
   real(kind=dp), intent(out)     :: real_field(1:3, 0:TGRID_X,&
                                               &0:TGRID_Y, 0:TGRID_Z)
                                           ! Result of transforms

   complex(kind=cdp), allocatable :: fftfield(:,:,:) ! Memory for FFT

   ! Allocate storage for performing FFTs
   allocate(fftfield(0:TGRID_X,0:TGRID_Y,0:TGRID_Z))

#if defined(DOUBLE_PRECISION)
   call dfftw_plan_dft_3d(plan, TURB_GS, TURB_GS, TURB_GS, fftfield, &
      & fftfield, FFTW_BACKWARD, FFTW_ESTIMATE)
#else
   call sfftw_plan_dft_3d(plan, TURB_GS, TURB_GS, TURB_GS, fftfield, &
      & fftfield, FFTW_BACKWARD, FFTW_ESTIMATE)
#endif

   do d=1,3
      fftfield = complex_field(d,:,:,:)
#if defined(DOUBLE_PRECISION)
      call dfftw_execute_dft(plan, fftfield, fftfield)
#else
      call sfftw_execute_dft(plan, fftfield, fftfield)
#endif
      real_field(d,:,:,:) = real(fftfield, dp) / (turb_gs_real**3)
   end do

#if defined(DOUBLE_PRECISION)
   call dfftw_destroy_plan(plan)
#else
   call sfftw_destroy_plan(plan)
#endif

   deallocate(fftfield)

end subroutine FFT_3D

subroutine proj_rms_norm(sol_frac_in, P)
   use turb_commons
   implicit none
   ! Calculate the (empirically measured and fitted) normalization for
   ! projection of random vectors with solenoidal fraction 'sol_frac'
   real(kind=dp), intent(in)      :: sol_frac_in ! Solenoidal fraction
   real(kind=dp), intent(out)     :: P           ! Normalization constant

#if NDIM==3
   P = (0.797d0 * sol_frac_in**2) - (0.529d0 * sol_frac_in) + 0.568d0

   ! for reference, to maintain magnitude of vectors
   !P = (0.563 * sol_frac_in**2) - (0.258 * sol_frac_in) + 0.487
#else
   ! Not tested for NDIM /= 3
   P = 1.0d0
#endif

end subroutine proj_rms_norm

subroutine power_rms_norm(power_in, P)
   use turb_commons
   implicit none
   ! Calculate the normalizations for FFT of initial power spectrum
   real(kind=dp), intent(in)  :: power_in(0:TGRID_X,&
                                         &0:TGRID_Y,0:TGRID_Z)
   real(kind=dp), intent(out) :: P  ! Normalization constant
   complex(kind=cdp)           :: complex_field(1:NDIM, 0:TGRID_X,&
                                               &0:TGRID_Y, 0:TGRID_Z)
                                                 ! Complex field to transform
   real(kind=dp)              :: real_field(1:NDIM, 0:TGRID_X,&
                                           &0:TGRID_Y, 0:TGRID_Z)
                                                 ! Result of transforms
   integer                    :: d               ! Dimension counter

   do d=1,NDIM
      complex_field(d,:,:,:) = cmplx(power_in, kind=cdp)
   end do

#if NDIM==1
   call FFT_1D(complex_field, real_field)
#elif NDIM==2
   call FFT_2D(complex_field, real_field)
#else
   call FFT_3D(complex_field, real_field)
#endif

   P = sqrt(sum(real_field**2)/size(real_field))

end subroutine power_rms_norm

subroutine turb_force_calc(ncache, x_cell, rho, aturb)
   use turb_commons
   implicit none
   integer, intent(in)        :: ncache
   real(kind=dp), intent(in)  :: x_cell(1:ndim,1:nvector) ! Positions
   real(kind=dp), intent(in)  :: rho(1:nvector)           ! Densities
   real(kind=dp), intent(out) :: aturb(1:ndim, 1:nvector) ! Turbulent forcing

   integer                    :: nok                 ! no. of OK cells
   integer                    :: ok_cell(1:nvector)  ! 'ok' cells
   integer                    :: i                   ! cell counter
   real(kind=dp)              :: r(1:ndim,1:nvector) ! Position in turb grid
   real(kind=dp)              :: dr1(1:ndim,1:nvector)
                                              ! Position in cell
   real(kind=dp)              :: dr2(1:ndim,1:nvector)
                                              ! 1 - position in cell
   integer                    :: bmin(1:ndim,1:nvector)
                                              ! 'top-left' corner of box
   integer                    :: bmax(1:ndim,1:nvector)
                                              ! 'bottom-right' corner of box
   real(kind=dp)              :: cube_vals(1:ndim,1:nvector,1:twotondim)
                                              ! Values of turbulent forcing
                                              ! from top left to bottom right
   real(kind=dp)              :: interp(1:ndim,1:nvector,1:twotondim)
                                              ! Interpolation weights
#ifdef ANALYTIC_FORCING_TEST
   real(kind=dp)              :: c            ! y-x position
#endif

   aturb = 0.0

   nok = 0
   do i=1,ncache
   ! Position of particle in 'grid' space
      if (rho(i) < turb_min_rho) then
         continue ! Less than minimum density for adding turbulence
      else
         nok = nok + 1
         ok_cell(nok) = i
         r(:,nok) = x_cell(:,i)
      end if
   end do

   ! Find ids of top left and bottom right corner of encompassing cube
   ! These can be the same if a particle is right on the boundary, but won't
   ! affect the result.
   bmin(:,1:nok) = floor(r(:,1:nok))
   bmax(:,1:nok) = ceiling(r(:,1:nok))

   ! Right-hand edge is equal to left-hand edge (periodic),
   ! which gives TURB_GS + 1 interpolation points
   where (bmin(:,1:nok)==TURB_GS) bmin(:,1:nok)=0
                                    ! this only happens for r==TURB_GS
   where (bmax(:,1:nok)==TURB_GS) bmax(:,1:nok)=0

   dr1(:,1:nok) = r(:,1:nok) - real(floor(r(:,1:nok)),dp)
   dr2(:,1:nok) = 1.0_dp - dr1(:,1:nok)

   ! Find cube values
#if NDIM==1
   do i=1,nok
      cube_vals(:,i,1) = afield_now(:, bmin(1,i), 0, 0)
      cube_vals(:,i,2) = afield_now(:, bmax(1,i), 0, 0)
   end do

   do i=1,nok
      interp(:,i,1) = dr2(1,i)
      interp(:,i,2) = dr1(1,i)
   end do
#elif NDIM==2
   do i=1,nok
      cube_vals(:,i,1) = afield_now(:, bmin(1,i), bmin(2,i), 0)
      cube_vals(:,i,2) = afield_now(:, bmax(1,i), bmin(2,i), 0)
      cube_vals(:,i,3) = afield_now(:, bmin(1,i), bmax(2,i), 0)
      cube_vals(:,i,4) = afield_now(:, bmax(1,i), bmax(2,i), 0)
   end do

   do i=1,nok
      interp(:,i,1) = dr2(1,i) * dr2(2,i)
      interp(:,i,2) = dr1(1,i) * dr2(2,i)
      interp(:,i,3) = dr2(1,i) * dr1(2,i)
      interp(:,i,4) = dr1(1,i) * dr1(2,i)
   end do
#else
   do i=1,nok
      cube_vals(:,i,1) = afield_now(:, bmin(1,i), bmin(2,i), bmin(3,i))
      cube_vals(:,i,2) = afield_now(:, bmax(1,i), bmin(2,i), bmin(3,i))
      cube_vals(:,i,3) = afield_now(:, bmin(1,i), bmax(2,i), bmin(3,i))
      cube_vals(:,i,4) = afield_now(:, bmax(1,i), bmax(2,i), bmin(3,i))
      cube_vals(:,i,5) = afield_now(:, bmin(1,i), bmin(2,i), bmax(3,i))
      cube_vals(:,i,6) = afield_now(:, bmax(1,i), bmin(2,i), bmax(3,i))
      cube_vals(:,i,7) = afield_now(:, bmin(1,i), bmax(2,i), bmax(3,i))
      cube_vals(:,i,8) = afield_now(:, bmax(1,i), bmax(2,i), bmax(3,i))
   end do

   ! Find interpolation values
   do i=1,nok
      interp(:,i,1) = dr2(1,i) * dr2(2,i) * dr2(3,i)
      interp(:,i,2) = dr1(1,i) * dr2(2,i) * dr2(3,i)
      interp(:,i,3) = dr2(1,i) * dr1(2,i) * dr2(3,i)
      interp(:,i,4) = dr1(1,i) * dr1(2,i) * dr2(3,i)
      interp(:,i,5) = dr2(1,i) * dr2(2,i) * dr1(3,i)
      interp(:,i,6) = dr1(1,i) * dr2(2,i) * dr1(3,i)
      interp(:,i,7) = dr2(1,i) * dr1(2,i) * dr1(3,i)
      interp(:,i,8) = dr1(1,i) * dr1(2,i) * dr1(3,i)
   end do
#endif

   ! Find interpolated value
   dr1(1:ndim,1:nok) = sum(interp(:,1:nok,:) * cube_vals(:,1:nok,:), dim=3)
#ifdef ANALYTIC_FORCING_TEST
   ! Test with angled periodic forcing; y=x line are diverging lines,
   ! y=x+0.5 and y=x-0.5 are converging lines.
   do i=1,nok
      c = (r(2,i) - r(1,i)) / turb_gs_real
      if (c==-0.5 .OR. c==0.0 .OR. c==0.5) then
         ! on a line; no forcing
         aturb(1:2, ok_cell(i)) = 0.0
      else if (c > 0.5) then
         ! top/left corner
         aturb(1:2, ok_cell(i)) = 0.01 * (/1.0, -1.0/)
      else if (c > 0.0) then
         ! upper left
         aturb(1:2, ok_cell(i)) = 0.01 * (/-1.0, 1.0/)
      else if (c > -0.5) then
         ! lower left
         aturb(1:2, ok_cell(i)) = 0.01 * (/1.0, -1.0/)
      else
         ! bottom/right corner
         aturb(1:2, ok_cell(i)) = 0.01 * (/-1.0, 1.0/)
      end if
      !aturb(1,ok_cell(i)) = sign(1.0_dp, 0.5_dp - (r(1,i)/turb_gs_real))
   end do
#else
   do i=1,nok
      aturb(:,ok_cell(i)) = dr1(:,i)
   end do
#endif

end subroutine turb_force_calc

subroutine current_turb_rms(rms_val)
   use turb_commons
   implicit none
   real (kind=dp), intent(out) :: rms_val
   integer                     :: i, j, k
   real (kind=dp)              :: asqd

   rms_val = 0.0
   do k=0,TGRID_Z
      do j=0,TGRID_Y
         do i=0,TGRID_X
            asqd = sum(afield_now(1:ndim, i, j, k)**2)
            rms_val = rms_val + asqd
         end do
      end do
   end do

   rms_val = sqrt(rms_val / (turb_gs_real**NDIM))

end subroutine current_turb_rms

! PRNG of Marsaglia
subroutine spin_up(s)
   use turb_parameters
   integer(ILP), intent(inout)   :: s(4)
   integer                       :: i
   integer(ILP)                  :: r
   do i=1,10000
      call kiss64_core(s, r)
   end do
end subroutine spin_up

subroutine kiss64_core(s, r)
   use turb_parameters
   integer(ILP), intent(inout) :: s(4)
   integer(ILP), intent(out)   :: r
   integer(ILP) :: t
   t = ishft(s(1), 58) + s(4)
   if (ishft(s(1), -63) .eq. ishft(t, -63)) then
      s(4) = ishft(s(1), -6) + ishft(s(1), -63)
   else
      s(4) = ishft(s(1), -6) + 1 - ishft(s(1) + t, -63)
   endif
   s(1) = t + s(1)
   s(2) = ieor(s(2), ishft(s(2),13))
   s(2) = ieor(s(2), ishft(s(2),-17))
   s(2) = ieor(s(2), ishft(s(2),43))
   s(3)= 6906969069_ILP * s(3) + 1234567_ILP
   r = s(1) + s(2) + s(3)
end subroutine kiss64_core

subroutine kiss64_double(N, s, out_array)
   use turb_parameters
   integer, intent(in)           :: N
   integer(ILP), intent(inout)   :: s(4)
   real(dp), intent(out)         :: out_array(1:N)
   integer                       :: i
   integer(ILP)                  :: int_array(1:N)
   integer(ILP), parameter       :: randmax = 9223372036854775807_ILP
#ifndef DOUBLE_PRECISION
   double precision              :: dbl_array(1:N)
#endif

   do i=1,N
      call kiss64_core(s, int_array(i))
   end do

   int_array = iand(int_array, z'FFFFFFFFFFFFF')
   int_array = ieor(int_array, z'3FF0000000000000')
#ifdef DOUBLE_PRECISION
   out_array = transfer(int_array, out_array) - 1.0_dp
#else
   dbl_array = transfer(int_array, out_array)
   out_array = real(dbl_array, dp)
#endif

end subroutine kiss64_double
#endif
