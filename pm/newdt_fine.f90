subroutine newdt_fine(ilevel)
  use pm_commons
  use amr_commons
  use hydro_commons
  use poisson_commons, ONLY: self_gravity, gravity_rho_ana_type
#ifdef RT
  use rt_parameters, ONLY: rt_advect, rt_nsubcycle
#endif
#if USE_TURB==1
  use turb_commons
#endif
  use constants, ONLY: pi
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::info
#endif
  integer::ilevel
  !-----------------------------------------------------------
  ! This routine compute the time step using 4 constraints:
  ! 1- a Courant-type condition using particle velocity
  ! 2- the gravity free-fall time
  ! 3- 10% maximum variation for aexp
  ! 4- maximum step time for ATON
  ! This routine also computes the particle kinetic energy.
  !-----------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart
  integer::npart1,ip,isink,idim
  integer,dimension(1:nvector),save::ind_part
  real(kind=8)::dt_loc,dt_all,ekin_loc,ekin_all
  real(dp)::tff,fourpi,threepi2
  real(dp)::dx,dx_loc,nx_loc,scale
  real(dp)::vsink2,vsink_max
#ifdef ATON
  real(dp)::aton_time_step,dt_aton
#endif
#ifdef RT
  real(dp)::dt_rt
#endif
  integer :: iout_frw
  real(dp) :: t_target, a_target

  logical :: ok

!$omp threadprivate(ind_part)

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  threepi2=3.0d0*pi**2

  ! Save old time step
  dtold(ilevel)=dtnew(ilevel)

  ! Maximum time step
  dtnew(ilevel)=boxlen/smallc
  if(poisson.and.(self_gravity.or.gravity_rho_ana_type>0))then
     fourpi=4.0d0*pi
     if(cosmo)fourpi=1.5d0*omega_m*aexp
     if (sink)then
        tff=sqrt(threepi2/8/fourpi/(rho_max(ilevel)+rho_sink_tff(ilevel)))
     else
        tff=sqrt(threepi2/8/fourpi/(rho_max(ilevel)+smallr))
     end if
     dtnew(ilevel)=MIN(dtnew(ilevel),courant_factor*tff)
  end if
  if(cosmo)then
     dtnew(ilevel)=MIN(dtnew(ilevel),0.1d0/hexp)
  end if

  ! Check sink velocity
  if(sink)then
     dx=0.5d0**ilevel
     nx_loc=dble(icoarse_max-icoarse_min+1)
     scale=boxlen/nx_loc
     dx_loc=dx*scale
     vsink_max=0d0
     do isink=1,nsink
        if(.not. new_born(isink))then
           vsink2=0d0
           do idim=1,ndim
              vsink2=vsink2+vsink(isink,idim)**2
           end do
          if(sink_descent)then
             vsink_max=MAX(vsink_max,sqrt(vsink2)+graddescent_over_dt(isink))
          else
             vsink_max=MAX(vsink_max,sqrt(vsink2))
          endif
        endif
     end do
     if(vsink_max.GT.0d0)then
        dtnew(ilevel)=MIN(dtnew(ilevel),courant_factor*dx_loc/vsink_max)
     endif
  endif

#ifdef ATON
  ! Maximum time step for ATON
  if(aton)then
     dt_aton = aton_time_step()
     if(dt_aton>0d0)then
        dtnew(ilevel)=MIN(dtnew(ilevel),dt_aton)
     end if
  end if
#endif

#ifdef RT
  ! Maximum time step for radiative transfer
  if(rt_advect)then
     call get_rt_courant_dt(dt_rt, ilevel)
     dtnew(ilevel) = 0.99999d0 * MIN(dtnew(ilevel), dt_rt * rt_nsubcycle)
     if(static) RETURN
  endif
#endif

#if USE_TURB==1
  ! Maximum time step from turbulent forcing
  if (driven_turb) then
     dtnew(ilevel) = min(dtnew(ilevel), turb_dt)
  end if
#endif

  if(pic) then

     dt_all=dtnew(ilevel); dt_loc=dt_all
     ekin_all=0; ekin_loc=0

     ! Compute maximum time step on active region
     if(numbl(myid,ilevel)>0)then
!$omp parallel private(jgrid,igrid,npart1,ipart,jpart,ok,ip) reduction(MIN:dt_loc) reduction(+:ekin_loc)
        ! Loop over grids
        ip=0
!$omp do schedule(dynamic,10)
        do jgrid=1,active(ilevel)%ngrid
           igrid=active(ilevel)%igrid(jgrid)
           npart1=numbp(igrid)   ! Number of particles in the grid
           if(npart1>0)then
              ! Loop over particles
              ipart=headp(igrid)
              do jpart=1,npart1
                 ! MC_tracer ================
                 ! skip tracer particles here
                 if (tracer) then
                    ok = is_not_tracer(typep(ipart))
                 else
                    ok = .true.
                 end if
                 ! End MC Tracer ============

                 if (ok) then
                    ip=ip+1
                    ind_part(ip)=ipart
                    if(ip==nvector)then
                       call newdt2(ind_part,dt_loc,ekin_loc,ip,ilevel)
                       ip=0
                    end if
                 end if
                 ipart=nextp(ipart)    ! Go to next particle
              end do
              ! End loop over particles
           end if
        end do
!$omp end do nowait
        ! End loop over grids
        if(ip>0)call newdt2(ind_part,dt_loc,ekin_loc,ip,ilevel)
!$omp end parallel
     end if

     ! Minimize time step over all cpus
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(dt_loc,dt_all,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
          & MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(ekin_loc,ekin_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
          & MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
     dt_all=dt_loc
     ekin_all=ekin_loc
#endif
     ekin_tot=ekin_tot+ekin_all
     dtnew(ilevel)=MIN(dtnew(ilevel),dt_all)

  end if

  if(hydro)call courant_fine(ilevel)

  !------------------------------------------------------------
  ! Adjust timestep to ensure outputs occur exactly at target a
  ! Target time is obtained via interpolation of FRW table
  !------------------------------------------------------------
  if (exact_output_time) then
     if (ilevel == levelmin) then
        if (iout <= noutput) then
           if (cosmo) then
              ! --- cosmological case (use aout → interpolate t_target) ---
              a_target = aout(iout)

              ! Find bracketing indices in aexp_frw
              iout_frw = 1
              do while (iout_frw <= n_frw)
                 if (aexp_frw(iout_frw-1) >= a_target .and. a_target >= aexp_frw(iout_frw)) exit
                 iout_frw = iout_frw + 1
              end do

              if (iout_frw <= n_frw) then
                 ! Linear interpolation of target time
                 t_target = tau_frw(iout_frw  ) * (a_target - aexp_frw(iout_frw-1)) / &
                            (aexp_frw(iout_frw  ) - aexp_frw(iout_frw-1)) + &
                            tau_frw(iout_frw-1) * (a_target - aexp_frw(iout_frw  )) / &
                            (aexp_frw(iout_frw-1) - aexp_frw(iout_frw  ))

                 ! Adjust timestep to hit the target output time exactly
                 if (t < t_target - eps_t .and. t + dtnew(ilevel) >= t_target - eps_t) then
                    dtnew(ilevel) = t_target - t
                 end if

              end if
           else
              ! --- non-cosmological case (use tout directly) ---
              t_target = tout(iout)

              if (t < t_target - eps_t .and. t + dtnew(ilevel) >= t_target - eps_t) then
                 dtnew(ilevel) = t_target - t
              end if

           end if
        end if
     end if
  end if

111 format('   Entering newdt_fine for level ',I2)

end subroutine newdt_fine
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine newdt2(ind_part,dt_loc,ekin_loc,nn,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
  real(kind=8)::dt_loc,ekin_loc
  integer::nn,ilevel
  integer,dimension(1:nvector)::ind_part

  integer::i,idim,nx_loc
  real(dp)::dx,dx_loc,scale,dtpart
  real(dp),dimension(1:nvector),save::v2,mmm
  real(dp),dimension(1:nvector,1:ndim)::vvv

!$omp threadprivate(v2,mmm)

  ! Compute time step
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  v2(1:nn)=0.0D0
  do idim=1,ndim
     do i=1,nn
        vvv(i, idim) = vp(ind_part(i), idim)
        v2(i)=max(v2(i),vvv(i, idim)**2)
        ! v2(i)=v2(i)+vp(ind_part(i),idim)**2
     end do
  end do
  do i=1,nn
     if(v2(i)>0.0D0)then
        dtpart=courant_factor*dx_loc/sqrt(v2(i))
        dt_loc=MIN(dt_loc,dtpart)
     end if
  end do

  ! Fetch mass
  do i = 1, nn
     mmm(i) = mp(ind_part(i))
  end do

  ! Compute kinetic energy
  do idim=1,ndim
     do i=1,nn
        ekin_loc=ekin_loc+0.5D0*mmm(i)*vvv(i, idim)**2
     end do
  end do

end subroutine newdt2
