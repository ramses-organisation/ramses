subroutine init_turb
  use turb_commons
  implicit none
  !-------------------------------------------------------------------
  ! Set up the turbulence module. Two independent features live in it:
  !
  !   driven_turb  -> driving, a forcing field applied at every timestep
  !   initial_turb -> a turbulent velocity field applied once, at t=0
  !
  ! They share only the grid dimensions and the random number generator;
  ! everything else belongs to one or the other and is set up by
  ! init_turb_driving / init_turb_initial.
  !-------------------------------------------------------------------

   integer       :: i                      ! Loop variable
   integer       :: n_seed=4               ! Length of random seed, 4 for KISS64
   integer       :: clock                  ! Integer clock time

   ! Only the root task holds the spectra and the PRNG state; the other tasks
   ! receive the resulting fields at the rendezvous below. Note myid==1 also in
   ! serial builds, so these branches need no preprocessor guard.
   if (myid == 1 .AND. nrestart == 0) then
      ! Set up random seed (modified from gfortran docs). On a restart the state
      ! is restored from file by read_turb_fields instead.
      if (turb_seed == -1) then
          call system_clock(count=clock)
          kiss64_state = clock + 37 * (/(i-1,i=1,n_seed)/)
      else
          kiss64_state = turb_seed
      end if
      call spin_up(kiss64_state)
   end if

   ! The driving is set up first, so that its random sequence - and hence a
   ! driven run - is unaffected by whether initial_turb is switched on.
   if (driven_turb) call init_turb_driving

   ! The initial field is only ever used on a fresh start, by init_flow_fine.
   ! On a restart it must not be drawn at all: that would consume random numbers
   ! after read_turb_fields has restored the PRNG state, and the driving would
   ! then diverge from the run being continued.
   if (initial_turb .AND. nrestart == 0) call init_turb_initial

end subroutine init_turb
!#####################################################################
!#####################################################################
!#####################################################################
subroutine init_turb_driving
  use turb_commons
  implicit none
  !-------------------------------------------------------------------
  ! Set up the turbulent driving: the storage for the force, the two
  ! bracketing forcing fields, their power spectrum and normalisation.
  !-------------------------------------------------------------------

   integer       :: i                      ! Loop variable
   integer       :: all_stat(1:4)          ! Allocation statuses

   real(kind=dp)        :: power_norm      ! Normalization from power spectrum
   real(kind=dp)        :: proj_norm       ! Normalization from projection
   real(kind=dp)        :: OU_norm         ! Normalization for OU process

   integer, parameter :: instant_turb_mult=5
                                         ! Number of autocorrelation times
                                         ! to evolve with instant turbulence
   real(kind=dp), parameter :: instant_turb_percent=10.0
                                         ! Display percentage at this interval
   real(kind=dp)      :: cur_percent     ! Percentage currently tracking

   ! Tasks to always be done (including MPI non-root tasks)
   ! ---------------------------------------------------------------------------

   all_stat = 0

   ! Allocate turbulent force storage and the fields it is built from
   allocate(fturb(1:ncoarse+twotondim*ngridmax,1:ndim), stat=all_stat(1))
   allocate(afield_last(1:NDIM,0:TGRID_X,0:TGRID_Y,0:TGRID_Z), stat=all_stat(2))
   allocate(afield_next(1:NDIM,0:TGRID_X,0:TGRID_Y,0:TGRID_Z), stat=all_stat(3))
   allocate(afield_now(1:NDIM,0:TGRID_X,0:TGRID_Y,0:TGRID_Z), stat=all_stat(4))

   if (any(all_stat /= 0)) stop 'Out of memory in init_turb_driving!'

   ! Set turbulence update time from autocorrelation time and number of substeps
   turb_dt = turb_T / real(turb_Ndt,dp)

   ! Tasks that do not need to be done by MPI non-root tasks
   ! ---------------------------------------------------------------------------
   ! Only the root task holds the turbulent spectra (turb_last, turb_next), the
   ! power spectrum and the PRNG state; the other tasks receive the resulting
   ! acceleration fields at the rendezvous below. Note myid==1 also in serial
   ! builds, so this branch needs no preprocessor guard.

   if (myid == 1) then

      ! Allocate the spectra, held by the root task only
      allocate(turb_last(1:NDIM,0:TGRID_X,0:TGRID_Y,0:TGRID_Z), stat=all_stat(1))
      allocate(turb_next(1:NDIM,0:TGRID_X,0:TGRID_Y,0:TGRID_Z), stat=all_stat(2))
      allocate(power_spec(0:TGRID_X,0:TGRID_Y,0:TGRID_Z), stat=all_stat(3))

      if (any(all_stat(1:3) /= 0)) stop 'Out of memory in init_turb_driving!'

      ! Set decay fraction per timestep dt
      turb_decay_frac = turb_dt / turb_T ! == 1 / turbNdt

      ! Set solenoidal fraction from compressive fraction
      sol_frac = 1.0_dp - comp_frac

      ! Set turbulent field time
      turb_next_time = t ! will be updated later

      ! Set the power distribution of the driving
      call build_power_spectrum(forcing_power_spectrum, power_spec)

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
         call add_turbulence(turb_next, power_spec, comp_frac, sol_frac, turb_dt)

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
         end if
      end if
   end if

   ! Tasks to always be done (including MPI non-root tasks)
   ! ---------------------------------------------------------------------------

#ifndef WITHOUTMPI
   ! Rendezvous: for the root task this sends the fields built above, for
   ! every other task it receives them into the arrays allocated at the top.
   call mpi_share_turb_fields(.TRUE.)
#endif

   ! Set up afield_now, i.e. the field the driving actually applies.
   ! To be done after the rendezvous.
   if (turb_evolving) then
      ! Interpolate between the two bracketing fields. The time fraction is 0
      ! on a fresh start (with or without instant_turb) and lands part way
      ! through the interval on a restart.
      call turb_interpolate_now
   else
      ! A single static field, held for the whole run by turb_check_time. It
      ! must NOT be interpolated - blending two independent fields gives an
      ! rms below turb_rms.
      afield_now = afield_last
   end if

   ! Renormalise onto the requested amplitude. afield_last is a single draw
   ! of the Ornstein-Uhlenbeck process, whose rms fluctuates by a few percent
   ! about sqrt(ndim)*turb_rms. Forced on for a static field.
   if (turb_exact_rms) call turb_normalise_rms

end subroutine init_turb_driving
!#####################################################################
!#####################################################################
!#####################################################################
subroutine init_turb_initial
  use turb_commons
  implicit none
  !-------------------------------------------------------------------
  ! Draw the turbulent field used for the initial velocity.
  !
  ! This is a single draw, not a process: there is no time evolution,
  ! no autocorrelation time, and no analytic normalisation. The field
  ! only has to have the right shape - add_turb_init_velocity measures
  ! its rms and scales it onto initial_turb_vrms - so none of the
  ! driving parameters enter here.
  !-------------------------------------------------------------------

   integer :: all_stat                    ! Allocation status

   real(kind=dp), allocatable     :: power_spec_init(:,:,:)
                                         ! Power distribution of the draw, which
                                         ! may differ from the driving spectrum
   complex(kind=cdp), allocatable :: turb_ic(:,:,:,:)
                                         ! Scratch spectrum for the draw

   ! Allocated on every task: the root task fills it, the others receive it
   allocate(afield_init(1:NDIM,0:TGRID_X,0:TGRID_Y,0:TGRID_Z), stat=all_stat)
   if (all_stat /= 0) stop 'Out of memory in init_turb_initial!'

   if (myid == 1) then

      allocate(power_spec_init(0:TGRID_X,0:TGRID_Y,0:TGRID_Z), stat=all_stat)
      if (all_stat /= 0) stop 'Out of memory in init_turb_initial!'
      call build_power_spectrum(initial_turb_spectrum, power_spec_init)

      ! Draw into a scratch spectrum of its own, so that with driving switched on
      ! the Ornstein-Uhlenbeck state in turb_next is left untouched and the two
      ! fields stay uncorrelated.
      allocate(turb_ic(1:NDIM,0:TGRID_X,0:TGRID_Y,0:TGRID_Z), stat=all_stat)
      if (all_stat /= 0) stop 'Out of memory in init_turb_initial!'

      ! add_turbulence takes the variance of the Wiener increment as its last
      ! argument. Starting from an empty field, that increment is the whole field
      ! and the variance is a pure scale factor, so a unit variance is the natural
      ! choice: the amplitude is set by initial_turb_vrms, not here.
      turb_ic = cmplx(0, 0, kind=cdp)
      call add_turbulence(turb_ic, power_spec_init, initial_turb_comp_frac,&
                           & 1.0_dp - initial_turb_comp_frac, 1.0_dp)
#if NDIM==1
      call FFT_1D(turb_ic(1,:,0,0), afield_init(1,:,0,0))
#elif NDIM==2
      call FFT_2D(turb_ic(:,:,:,0), afield_init(:,:,:,0))
#else
      call FFT_3D(turb_ic, afield_init)
#endif

      deallocate(turb_ic)
      deallocate(power_spec_init)
   end if

#ifndef WITHOUTMPI
   ! Rendezvous: for the root task this sends the field drawn above, for every
   ! other task it receives it into afield_init.
   call mpi_share_turb_init_field
#endif

end subroutine init_turb_initial
!#####################################################################
!#####################################################################
!#####################################################################
subroutine add_turb_init_velocity(ilevel)
  use amr_commons
  use hydro_commons
  use turb_commons
  implicit none
  integer::ilevel
  !-------------------------------------------------------------------
  ! Add the turbulent initial velocity field to the gas and update the
  ! total energy accordingly. Called once, from init_flow_fine, on a
  ! fresh start. This is not a force: the field is scaled by
  ! turb_init_vscale so that the resulting velocity dispersion is
  ! initial_turb_vrms, and added directly to the momentum.
  !-------------------------------------------------------------------
  integer::ncache,ngrid,i,igrid,iskip,ind,idim
  integer::ix,iy,iz,nx_loc
#ifdef SOLVERmhd
  integer::nndim=3
#else
  integer::nndim=ndim
#endif
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  real(kind=dp) :: x_cell(1:ndim,1:nvector)     ! Cell positions
  real(kind=dp) :: rho(1:nvector)               ! Cell densities
  real(kind=dp) :: vturb(1:ndim,1:nvector)      ! Turbulent velocity field
  real(dp),dimension(1:3,1:twotondim)::xc
  real(dp),dimension(1:3)::skip_loc
  real(dp)::dx,dx_loc,scale,d,turb_init_vscale,rms_now

!$omp threadprivate(ind_grid,ind_cell)

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Scale the field onto the requested velocity dispersion. Measuring the rms
  ! here keeps this correct whatever amplitude the field was generated at.
  call current_turb_rms(afield_init, rms_now)
  if (rms_now > 0.0_dp) then
     turb_init_vscale = initial_turb_vrms / rms_now
  else
     turb_init_vscale = 0.0_dp
  end if

  ! Mesh size at level ilevel in coarse cell units
  dx=0.5D0**ilevel

  ! Rescaling factors
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=turb_gs_real/dble(nx_loc)
  dx_loc=dx*scale

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(1,ind)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(2,ind)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(3,ind)=(dble(iz)-0.5D0)*dx
  end do

  ! Loop over active grids by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do

     ! Loop over cells
     do ind=1,twotondim

        ! Gather cell indices
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        ! Gather cell centre positions
        do i=1,ngrid
           do idim=1,ndim
              x_cell(idim,i)=xg(ind_grid(i),idim)+xc(idim,ind)
           end do
        end do
        ! Rescale position from code units to 0->turb_gs_real units
        do i=1,ngrid
           do idim=1,ndim
              x_cell(idim,i)=(x_cell(idim,i)-skip_loc(idim))*scale
           end do
        end do

        ! Gather cell densities
        do i=1,ngrid
           rho(i) = uold(ind_cell(i), 1)
        end do

        ! Interpolate the turbulent field onto the cells
        call turb_force_calc(afield_init, ngrid, x_cell, rho, vturb)

        do i=1,ngrid
           d = max(uold(ind_cell(i),1),smallr)

           ! Remove the kinetic energy, leaving internal (+ magnetic) energy
           do idim=1,nndim
              uold(ind_cell(i),neul) = uold(ind_cell(i),neul) &
                   & - 0.5d0*uold(ind_cell(i),idim+1)**2/d
           end do

           ! Add the turbulent velocity to the momentum
           do idim=1,ndim
              uold(ind_cell(i),idim+1) = uold(ind_cell(i),idim+1) &
                   & + d*vturb(idim,i)*turb_init_vscale
           end do

           ! Restore the total energy with the new kinetic energy
           do idim=1,nndim
              uold(ind_cell(i),neul) = uold(ind_cell(i),neul) &
                   & + 0.5d0*uold(ind_cell(i),idim+1)**2/d
           end do
        end do

     end do
     ! End loop over cells

  end do
  ! End loop over grids

111 format('   Entering add_turb_init_velocity for level',i2)

end subroutine add_turb_init_velocity
