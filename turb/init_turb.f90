subroutine init_turb
  use turb_commons
  implicit none
  !--------------------------------------------------
  ! Local variables
  !--------------------------------------------------

   integer       :: i, j, k                ! Loop variables
   integer       :: k_vec(1:3)             ! Wavevector
   integer       :: all_stat(1:4)          ! Allocation statuses

   integer              :: n_seed=4        ! Length of random seed, 4 for KISS64
   integer              :: clock           ! Integer clock time

   real(kind=dp)        :: power_norm      ! Normalization from power spectrum
   real(kind=dp)        :: proj_norm       ! Normalization from projection
   real(kind=dp)        :: OU_norm         ! Normalization for OU process

   real(kind=dp) :: turb_rms_now           ! Realised rms of the static field

   integer, parameter :: instant_turb_mult=5
                                         ! Number of autocorrelation times
                                         ! to evolve with instant turbulence
   real(kind=dp), parameter :: instant_turb_percent=10.0
                                         ! Display percentage at this interval
   real(kind=dp)      :: cur_percent     ! Percentage currently tracking

   ! Tasks to always be done (including MPI non-root tasks)
   ! ---------------------------------------------------------------------------

   all_stat = 0

   ! Allocate turbulent force storage
   allocate(fturb(1:ncoarse+twotondim*ngridmax,1:ndim), stat=all_stat(1))

   ! Allocate grids
   allocate(afield_last(1:NDIM,0:TGRID_X,0:TGRID_Y,0:TGRID_Z),&
            &stat=all_stat(2))
   allocate(afield_next(1:NDIM,0:TGRID_X,0:TGRID_Y,0:TGRID_Z),&
            &stat=all_stat(3))
   allocate(afield_now(1:NDIM,0:TGRID_X,0:TGRID_Y,0:TGRID_Z),&
            &stat=all_stat(4))

   if (any(all_stat /= 0)) stop 'Out of memory in init_turb!'

   ! Set grid spacing
   turb_space = (/1.d0, 1.d0, 1.d0/) / turb_gs_real

   ! Set turbulence update time from autocorrelation time and number of substeps
   turb_dt = turb_T / real(turb_Ndt,dp)

   ! Tasks that do not need to be done by MPI non-root tasks
   ! ---------------------------------------------------------------------------
   ! Only the root task holds the turbulent spectra (turb_last, turb_next), the
   ! power spectrum and the PRNG state; the other tasks receive the resulting
   ! acceleration fields at the rendezvous below. Note myid==1 also in serial
   ! builds, so this branch needs no preprocessor guard.

   if (myid == 1) then

      if (nrestart == 0) then
         ! Set up random seed (modified from gfortran docs)
         if (turb_seed == -1) then
             call system_clock(count=clock)
             kiss64_state = clock + 37 * (/(i-1,i=1,n_seed)/)
         else
             kiss64_state = turb_seed
         end if
         call spin_up(kiss64_state)
      end if

      ! Allocate grids
      allocate(turb_last(1:NDIM,0:TGRID_X,0:TGRID_Y,0:TGRID_Z), stat=all_stat(1))
      allocate(turb_next(1:NDIM,0:TGRID_X,0:TGRID_Y,0:TGRID_Z), stat=all_stat(2))
      allocate(power_spec(0:TGRID_X,0:TGRID_Y,0:TGRID_Z), stat=all_stat(3))

      if (any(all_stat /= 0)) stop 'Out of memory in init_turb!'

      ! Set decay fraction per timestep dt
      turb_decay_frac = turb_dt / turb_T ! == 1 / turbNdt

      ! Set solenoidal fraction from compressive fraction
      sol_frac = 1.0_dp - comp_frac

      ! Set turbulent field time
      turb_next_time = t ! will be updated later

      ! Set initial power distribution (should be parameterised in some fashion)
      do k=0,TGRID_Z
         if (k > TURB_GS / 2) then
            k_vec(3) = k - TURB_GS
         else
            k_vec(3) = k
         end if
         do j=0,TGRID_Y
            if (j > TURB_GS / 2) then
               k_vec(2) = j - TURB_GS
            else
               k_vec(2) = j
            end if
            do i=0,TGRID_X
               if (i > TURB_GS / 2) then
                  k_vec(1) = i - TURB_GS
               else
                  k_vec(1) = i
               end if
               call calc_power_spectrum(k_vec, power_spec(i,j,k))
            end do
         end do
      end do

      ! Calculate turbulent normalization
      ! Power normalization comes from FFT of initial power spectrum
      call power_rms_norm(power_spec, power_norm)
      ! Projection normalization was empirically estimated and fitted
      call proj_rms_norm(sol_frac, proj_norm)
      ! OU norm comes from standard deviation of OU process
      OU_norm = sqrt(turb_T / 2.0_dp)
      ! Combination of all factored (reciprocal for easy multiplication)
      turb_norm = 1.0_dp / (power_norm * proj_norm * OU_norm)

      if (nrestart > 0) then
         ! Restart - load turbulent fields from files and perform FFTs
         call read_turb_fields
#if NDIM==1
         call FFT_1D(turb_last(1,:,0,0), afield_last(1,:,0,0))
         call FFT_1D(turb_next(1,:,0,0), afield_next(1,:,0,0))
#elif NDIM==2
         call FFT_2D(turb_last(:,:,:,0), afield_last(:,:,:,0))
         call FFT_2D(turb_next(:,:,:,0), afield_next(:,:,:,0))
#else
         call FFT_3D(turb_last, afield_last)
         call FFT_3D(turb_next, afield_next)
#endif
         afield_last = afield_last * turb_norm * turb_rms
         afield_next = afield_next * turb_norm * turb_rms
      else
         ! Not a restart - set up initial field and perform FFT
         turb_next = cmplx(0, 0, kind=cdp)
         call add_turbulence(turb_next, turb_dt)

         ! Fourier transform
#if NDIM==1
         call FFT_1D(turb_next(1,:,0,0), afield_next(1,:,0,0))
#elif NDIM==2
         call FFT_2D(turb_next(:,:,:,0), afield_next(:,:,:,0))
#else
         call FFT_3D(turb_next, afield_next)
#endif
         afield_next = afield_next * turb_norm * turb_rms

         ! Call turb_next_field to create second field
         call turb_next_field
         if (instant_turb) then
            write (6,'(A,$)') " Evolving initial turbulent field..."
            call flush(6)
            cur_percent = instant_turb_percent
            do i=1,instant_turb_mult*turb_Ndt
               call turb_next_field
               turb_next_time = turb_last_time
               turb_last_time = turb_last_time - turb_dt
               write (6,'(A,$)') '.'
               if (100.0*real(i)/real(instant_turb_mult*turb_Ndt)&
                   & >= cur_percent) then
                  write (6, '(A,I0,A,$)') '(',int(cur_percent+0.5),'%)'
                  cur_percent = cur_percent + instant_turb_percent
                  if (cur_percent >= 99.99999999999) cur_percent = 200.0
               end if
               call flush(6)
            end do
            write (6,'(A,$)') " done."
            write (6,*)
            call flush(6)
            ! The spin-up loop above holds the clock at its fixed point
            ! (turb_last_time, turb_next_time) = (t, t + turb_dt) while it
            ! regenerates the fields, so it already leaves the correct bracket
            ! for the current time. Do NOT advance the clock again here: doing
            ! so used to leave turb_last_time one turb_dt ahead of t, so the
            ! initial afield_now came out as the extrapolation
            ! 2*afield_last - afield_next rather than a field on the segment.
         end if
      end if

   end if

   ! Tasks to always be done (including MPI non-root tasks)
   ! ---------------------------------------------------------------------------

#ifndef WITHOUTMPI
   ! Rendezvous: for the root task this call sends the fields built above, for
   ! every other task it receives them into the arrays allocated at the top.
   call mpi_share_turb_fields(.TRUE.)
#endif

   ! Set up afield_now
   select case (turb_type)
   case (1)
      ! Evolving forcing: interpolate between the two bracketing fields. The time
      ! fraction is 0 on a fresh start (with or without instant_turb) and lands
      ! part way through the interval on a restart.
      call turb_interpolate_now
   case (2)
      ! Fixed forcing: a single static field, held for the whole run by
      ! turb_check_time. It must NOT be interpolated - blending two independent
      ! fields gives an rms below turb_rms, which is what broke type 2 on
      ! restart. afield_last pairs with turb_last, which is what gets written to
      ! and reloaded from the restart files, so this is stable across a restart.
      afield_now = afield_last

      ! Renormalise onto the requested amplitude. afield_last is a single draw
      ! of the Ornstein-Uhlenbeck process, whose rms fluctuates by a few percent
      ! about sqrt(ndim)*turb_rms. For turb_type=1 that scatter is intended and
      ! averages out over the run, but here the field is frozen, so whatever
      ! deviation this particular draw happens to have would bias the driving
      ! for the whole simulation.
      ! current_turb_rms works on afield_now, and afield_last is rebuilt
      ! identically from turb_last after a restart, so the same factor is
      ! recovered and the field stays stable across a restart. It is also
      ! computed identically on every task, so no communication is needed.
      call current_turb_rms(turb_rms_now)
      if (turb_rms_now > 0.0_dp) then
         afield_now = afield_now * &
              & (sqrt(real(ndim,dp))*turb_rms / turb_rms_now)
      end if
   case (3)
      ! Decaying: no forcing is ever applied. This field is consumed once, by
      ! init_flow_fine, as the initial velocity field, then zeroed by
      ! turb_check_time. Any single field will do, but it must be a single
      ! field and it must be the same one on every task.
      afield_now = afield_last

      ! Renormalise for the same reason as turb_type=2: this single draw of the
      ! process would otherwise be a few percent off the requested amplitude.
      ! Here that amplitude is a velocity rather than a forcing, so this fixes
      ! the initial velocity dispersion, and hence the initial Mach number, at
      ! exactly sqrt(ndim)*turb_rms instead of leaving it to the draw.
      call current_turb_rms(turb_rms_now)
      if (turb_rms_now > 0.0_dp) then
         afield_now = afield_now * &
              & (sqrt(real(ndim,dp))*turb_rms / turb_rms_now)
      end if
   end select

end subroutine init_turb
