subroutine make_sn_stellar
  use pm_commons
  use amr_commons
  use hydro_commons
  use sink_feedback_parameters
  use constants, only:pi,pc2cm,mH,M_sun
  use mpi_mod
  implicit none

  integer:: ivar, ilevel, ind, ix, iy, iz, ngrid, iskip, idim
  integer:: i, nx_loc, igrid, ncache
  integer, dimension(1:nvector), save:: ind_grid, ind_cell
  real(dp):: dx, scale, dx_loc, vol_loc
  real(dp), dimension(1:3):: skip_loc
  real(dp), dimension(1:twotondim, 1:3):: xc
  logical, dimension(1:nvector), save:: ok
  real(dp), dimension(1:nvector, 1:ndim), save:: xx
  real(dp):: sn_r, sn_m, sn_p, sn_e, sn_d, sn_ed
  real(dp):: rr,pgas,dgas,ekin
  integer:: info
  real(dp),dimension(1:nvector,1:ndim)::x
  real(dp),dimension(1:3):: xshift, x_sn
  logical, save:: first = .true.
  real(dp), save:: xseed
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  logical, dimension(1:nstellarmax):: mark_del
  integer:: istellar,isink
  real(dp)::T_sn,sn_ed_lim,pnorm_sn,vol_sn
  real(dp)::mass_sn, mass_sn_all,dens_moy,r_cooling, Tsat_local
  real(dp)::pgas_check,pgas_check_all,egas_check,egas_check_all
  integer, parameter:: navg = 3
  real(dp), dimension(1:3):: avg_center
  real(dp):: avg_radius
  real(dp), dimension(1:navg):: avg_rpow
  real(dp), dimension(1:navg, 1:nvar+3):: avg_upow
  real(dp), dimension(1:navg):: avg
  real(dp):: norm, distance_sn, ekin_before, ekin_after
  logical::r_cooling_resolved=.false.

  if(.not. hydro)return
  if(ndim .ne. 3)return

  if(verbose)write(*,*)'Entering make_sn_stellar'

  if (first) then
     xseed = 0.5
     call random_number(xseed)
     first = .false.
  endif

  ! Mesh spacing in that level
  nx_loc = icoarse_max - icoarse_min + 1
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale = boxlen / dble(nx_loc)

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! radius of sphere where we dump the SN energy/momentum/mass (numerical, we choose 3 coarse cells)
  sn_r = 3.0d0*(0.5d0**levelmin)*scale
  ! namelist option to set a minimum radius for the SN remnant
  if(sn_r_sat .ne. 0) sn_r = max(sn_r, sn_r_sat * pc2cm / scale_l)

  ! we assume the SN energy and momentum is always the same
  sn_p = sn_p_ref
  sn_e = sn_e_ref
  Tsat_local = Tsat

  ! Loop over stellar objects to determine whether one is turning supernovae
  ! After the explosion, the object is removed from the list
  mark_del = .false.
  do istellar = 1, nstellar
    if(t - tstellar(istellar) < ltstellar(istellar)) cycle
    mark_del(istellar) = .true.

    ! find corresponding sink
    isink = 1
    do while ((isink.le.nsink) .and. (id_stellar(istellar) .ne. idsink(isink)))
      isink = isink + 1
    end do
    if (isink.gt.nsink) then
      write(*,*)"BUG: COULD NOT FIND SINK"
      call clean_stop
    endif

    ! the mass of the massive star
    sn_m = mstellar(istellar)

    !remove the mass that is dumped in the grid from the sink
    msink(isink) = msink(isink) - sn_m

    ! maximum distance the massive star could have traveled from the sink
    ! = the velocity dispersion (namelist) times the life time of the object
    distance_sn = ltstellar(istellar)*Vdisp

    ! determine random location within that distance
    ! find a random point within a sphere of radius < 1
    norm=2.
    do while (norm .gt. 1)
       norm=0.
       do idim = 1, ndim
          call random_number(xseed)
          xshift(idim) = (xseed-0.5)*2.
          norm = norm + xshift(idim)**2
       end do
    end do
    do idim = 1, ndim
        ! shift should not be more than half the box size
        xshift(idim) = xshift(idim) * min(distance_sn,0.5d0*boxlen)
    end do

    ! place the supernovae around sink particles
    x_sn(:) = xsink(isink, :) + xshift(:)
    ! if SN went outside the box, just invert the shift (this avoids having to check boundary conditions)
    do idim = 1,ndim
       if( (x_sn(idim) .lt. 0) .or. (x_sn(idim) .gt. boxlen)) x_sn(idim) = xsink(isink, idim) - xshift(idim)
    end do

    ! Now that we have the location of the SN, we check the properties of the surroundings
    avg_center = x_sn
    avg_radius = sn_r
    avg_rpow = 0.0d0
    avg_upow = 0.0d0
    ! avg_rpow(1) = 0 ; avg_upow(1, :) = 0 -> integrand(1) = 1
    ! avg_rpow(2) = 0 ; avg_upow(2, 1) = 1 -> integrand(2) = density
    ! avg_rpow(3) = 1 ; avg_upow(3, :) = 0 -> integrand(3) = radius
    avg_upow(2, 1) = 1.0d0
    avg_rpow(3) = 1.0d0
    call sphere_average(navg, avg_center, avg_radius, avg_rpow, avg_upow, avg)
    vol_sn = avg(1)
    mass_sn = avg(2) + sn_m ! region average + ejecta
    pnorm_sn = avg(3)

    ! density of the gas ejected by the SN, assumed to be distributed uniformly over the SN sphere
    sn_d = sn_m / vol_sn
    ! energy density of SN, assume uniform over SN sphere
    sn_ed = sn_e / vol_sn

    ! average density of gas in SN radius, after explosion (includes ejecta)
    dens_moy = mass_sn / vol_sn
    dens_moy = dens_moy*scale_d/mH !H/cc

    ! estimate cooling radius (Martizzi et al 2015)
    r_cooling = 6.3d0 * pc2cm/scale_l * (dens_moy/100d0)**(-0.42)
    ! if cooling radius resolved -> thermal feedback
    ! else -> momentum feedback
    r_cooling_resolved = (r_cooling > sn_r)
    ! determine SN momentum, see Iffrig and Hennebelle 2015
    ! NOT USED FOR NOW
    !sn_p = sn_p_ref * (dens_moy / 10.)**(-0.117647)

    pgas_check=0
    egas_check=0

    !now loop over cells again and dump energies, mass and momentum
    !loop over levels
    do ilevel = levelmin, nlevelmax
      ! Computing local volume (important for averaging hydro quantities)
      dx = 0.5d0**ilevel
      dx_loc = dx * scale
      vol_loc = dx_loc**ndim

      ! Cell center position relative to grid center position
      do ind=1,twotondim
        iz = (ind - 1) / 4
        iy = (ind - 1 - 4 * iz) / 2
        ix = (ind - 1 - 2 * iy - 4 * iz)
        xc(ind,1) = (dble(ix) - 0.5d0) * dx
        xc(ind,2) = (dble(iy) - 0.5d0) * dx
        xc(ind,3) = (dble(iz) - 0.5d0) * dx
      end do

      ! Loop over grids
      ncache=active(ilevel)%ngrid
      do igrid = 1, ncache, nvector
        ngrid = min(nvector, ncache - igrid + 1)
        do i = 1, ngrid
          ind_grid(i) = active(ilevel)%igrid(igrid + i - 1)
        end do

        ! Loop over cells
        do ind = 1, twotondim
          ! Gather cell indices
          iskip = ncoarse + (ind - 1) * ngridmax
          do i = 1, ngrid
            ind_cell(i) = iskip + ind_grid(i)
          end do

          ! Gather cell center positions
          do i = 1, ngrid
            xx(i, :) = xg(ind_grid(i), :) + xc(ind, :)
          end do
          ! Rescale position from coarse grid units to code units
          do idim=1,ndim
             do i=1,ngrid
                xx(i,idim)=(xx(i,idim)-skip_loc(idim))*scale
             end do
          end do

          ! Flag leaf cells
          do i = 1, ngrid
            ok(i) = (son(ind_cell(i)) == 0)
          end do

          do i = 1, ngrid
            if(ok(i)) then
                rr = 0.
                do idim=1,ndim
                   rr = rr + ((xx(i,idim) - x_sn(idim)) / sn_r)**2
                enddo
                rr = sqrt(rr)

                ! if cell inside SN sphere, dump fraction of energy and momentum
                if(rr < 1.) then
                  ! add ejecta density
                  uold(ind_cell(i), 1) = uold(ind_cell(i), 1) + sn_d
                  dgas = uold(ind_cell(i), 1)

                  ! momemtum feedback

                  ! compute velocity of the gas within this cell assuming energy equipartition
                  ! limit the velocity using Vsat
                  pgas = min(sn_p / pnorm_sn * rr, Vsat * dgas)
                  pgas_check = pgas_check + pgas * vol_loc

                  ! kinetic energy before SN
                  ekin_before=0.
                  do idim=1,ndim
                    ekin_before = ekin_before + ( (uold(ind_cell(i),idim+1))**2 ) / dgas / 2.
                  enddo

                  ! add momemtum
                  do idim=1,ndim
                     uold(ind_cell(i),idim+1) = uold(ind_cell(i),idim+1) + pgas * (xx(i,idim) - x_sn(idim)) / (rr * sn_r)
                  enddo

                  ! kinetic energy after
                  ekin_after=0.
                  do idim=1,ndim
                    ekin_after = ekin_after + ( (uold(ind_cell(i),idim+1))**2 ) / dgas / 2.
                  enddo

                  ! add extra kinetic energy
                  uold(ind_cell(i), 2+ndim) = uold(ind_cell(i), 2+ndim) + (ekin_after - ekin_before)
                  egas_check = egas_check + (ekin_after - ekin_before) * vol_loc

                  ! Thermal feedback

                  !before adding thermal energy make sure the temperature is not too high (too small timesteps otherwise)
                  T_sn = (sn_ed / dgas * (gamma-1.) ) * scale_T2
                  T_sn = min( T_sn , Tsat_local) / scale_T2
                  sn_ed_lim = T_sn * dgas / (gamma-1.)
                  egas_check = egas_check + sn_ed_lim * vol_loc
                  uold(ind_cell(i), 2+ndim) = uold(ind_cell(i), 2+ndim) + sn_ed_lim

                end if

            end if
          end do
          ! End loop over sublist of cells
        end do
        ! End loop over cells
      end do
      ! End loop over grids
    end do
    ! End loop over levels

#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(pgas_check,pgas_check_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    call MPI_ALLREDUCE(egas_check,egas_check_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#else
    pgas_check_all = pgas_check
    egas_check_all = egas_check
#endif

    if(myid == 1) write(*, *) "SN EVENT at", x_sn(1),x_sn(2),x_sn(3), "from sink", idsink(isink),"t=",t
    if(myid == 1) write(*, *) "Resolved?",r_cooling_resolved,"r_cool=",r_cooling,"r_sn=",sn_r,"density=",dens_moy,"starmass=",sn_m*(scale_d*scale_l**3)/M_sun
    if(myid == 1) write(*, *) "momentum (injected, expected)=", pgas_check_all * scale_d * scale_l**3 * scale_v, sn_p * scale_d * scale_l**3 * scale_v
    if(myid == 1) write(*, *) "energy (injected, expected)=", egas_check_all*(scale_d * scale_v**2 * scale_l**3), sn_e*(scale_d * scale_v**2 * scale_l**3)

  end do ! end of the loop over stellar objects

  call delete_stellar(mark_del)

  ! Update hydro quantities for split cells
  do ilevel = nlevelmax, levelmin, -1
    call upload_fine(ilevel)
    do ivar = 1, nvar
      call make_virtual_fine_dp(uold(1, ivar), ilevel)
    enddo
  enddo

end subroutine make_sn_stellar
!################################################################
!################################################################
!################################################################
!################################################################
subroutine sphere_average(navg, center, radius, rpow, upow, avg)
    use amr_parameters, only: boxlen, dp, hydro, icoarse_max, icoarse_min &
        & , jcoarse_min, kcoarse_min, levelmin, ndim, ngridmax, nlevelmax &
        & , nvector, twotondim, verbose
    use amr_commons, only: active, ncoarse, son, xg, myid
    use hydro_parameters, only: nvar
    use hydro_commons, only: uold
    use mpi_mod
    implicit none

    ! Integrate quantities over spheres
    ! The integrand is (r / radius)**rpow(iavg) * product(u(ivar)**upow(iavg, ivar), ivar=1:nvar)

    integer, intent(in):: navg                               ! Number of quantities
    real(dp), dimension(1:ndim), intent(in):: center         ! Sphere center
    real(dp), intent(in):: radius                            ! Sphere radius
    real(dp), dimension(1:navg), intent(in):: rpow           ! Power of radius in the integral
#ifdef SOLVERmhd
    real(dp), dimension(1:navg, 1:nvar+3), intent(in):: upow ! Power of hydro variables in the integral
#else
    real(dp), dimension(1:navg, 1:nvar), intent(in):: upow   ! Power of hydro variables in the integral
#endif
    real(dp), dimension(1:navg), intent(out):: avg           ! Averages

    integer:: i, ivar, ilevel, igrid, ind, ix, iy, iz, iskip, idim
    integer:: nx_loc, ncache, ngrid
    integer, dimension(1:nvector):: ind_grid, ind_cell
    logical, dimension(1:nvector):: ok

    real(dp):: scale, dx, dx_loc, vol_loc, rr
    real(dp), dimension(1:3):: skip_loc
    real(dp), dimension(1:twotondim, 1:3):: xc
    real(dp), dimension(1:nvector, 1:ndim):: xx

    integer:: info
    real(dp), dimension(1:navg):: avg_loc
    real(dp), dimension(1:navg):: integrand
    real(dp), dimension(1:navg):: utemp

    if(.not. hydro)return
    if(ndim .ne. 3)return

    if(verbose .and. myid == 1) write(*, *) 'Entering sphere_average'

    ! Mesh spacing in that level
    nx_loc = icoarse_max - icoarse_min + 1
    skip_loc=(/0.0d0,0.0d0,0.0d0/)
    if(ndim>0)skip_loc(1)=dble(icoarse_min)
    if(ndim>1)skip_loc(2)=dble(jcoarse_min)
    if(ndim>2)skip_loc(3)=dble(kcoarse_min)
    scale = boxlen / dble(nx_loc)

    avg_loc = 0.0d0

    do ilevel = levelmin, nlevelmax
        ! Computing local volume (important for averaging hydro quantities)
        dx = 0.5d0**ilevel
        dx_loc = dx * scale
        vol_loc = dx_loc**ndim

        ! Cell center position relative to grid center position
        do ind = 1, twotondim
            iz = (ind - 1) / 4
            iy = (ind - 1 - 4 * iz) / 2
            ix = (ind - 1 - 2 * iy - 4 * iz)
            if(ndim>0) xc(ind, 1) = (dble(ix) - 0.5d0) * dx
            if(ndim>1) xc(ind, 2) = (dble(iy) - 0.5d0) * dx
            if(ndim>2) xc(ind, 3) = (dble(iz) - 0.5d0) * dx
        end do

        ! Loop over grids
        ncache = active(ilevel)%ngrid
        do igrid = 1, ncache, nvector
            ngrid = min(nvector, ncache - igrid + 1)
            do i = 1, ngrid
                ind_grid(i) = active(ilevel)%igrid(igrid + i - 1)
            end do

            ! Loop over cells
            do ind = 1, twotondim
                ! Gather cell indices
                iskip = ncoarse + (ind - 1) * ngridmax
                do i = 1, ngrid
                    ind_cell(i) = iskip + ind_grid(i)
                end do

                ! Gather cell center positions
                do i = 1, ngrid
                    xx(i, :) = xg(ind_grid(i), :) + xc(ind, :)
                end do

                ! Rescale position from coarse grid units to code units
                do idim=1,ndim
                   do i=1,ngrid
                      xx(i,idim)=(xx(i,idim)-skip_loc(idim))*scale
                   end do
                end do

                ! Flag leaf cells
                do i = 1, ngrid
                    ok(i) = (son(ind_cell(i)) == 0)
                end do

                do i = 1, ngrid
                    if(ok(i)) then
                        rr = sqrt(sum(((xx(i, :) - center) / radius)**2))
                        if(rr < 1.) then
                            integrand = rr**rpow
                            where(abs(rpow) < 1.0d-10) ! Avoid NaNs of the form 0**0
                                integrand = 1.0d0
                            end where
#ifdef SOLVERmhd
                            do ivar = 1, nvar + 3
#else
                            do ivar = 1, nvar
#endif
                                utemp(:) = uold(ind_cell(i), ivar)
                                where(abs(upow(:, ivar)) < 1.0d-10) ! Avoid NaNs of the form 0**0
                                    utemp = 1.0d0
                                end where
                                integrand = integrand * utemp**upow(:, ivar)
                            end do
                            avg_loc = avg_loc + vol_loc * integrand
                        endif
                        ! End test on radius
                    endif
                    ! End test on leaf cells
                end do
                ! End loop over sublist of cells
            end do
            ! End loop over cells
        end do
        ! End loop over grids
    end do
    ! End loop over levels

#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(avg_loc, avg, navg, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
#else
    avg = avg_loc
#endif
end subroutine sphere_average
