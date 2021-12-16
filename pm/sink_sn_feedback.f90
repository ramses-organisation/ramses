subroutine make_sn_stellar
  use pm_commons
  use amr_commons
  use hydro_commons
  use sink_feedback_parameters
  use constants, only:pi,pc2cm
  use mpi_mod
  implicit none

  integer:: ivar
  integer:: ilevel, ind, ix, iy, iz, ngrid, iskip, idim
  integer:: i, nx_loc, igrid, ncache
  integer, dimension(1:nvector), save:: ind_grid, ind_cell
  real(dp):: dx, scale, dx_loc, vol_loc
  real(dp), dimension(1:3):: skip_loc
  real(dp), dimension(1:twotondim, 1:3):: xc
  logical, dimension(1:nvector), save:: ok
  real(dp), dimension(1:nvector, 1:ndim), save:: xx
  real(dp):: sn_r, sn_m, sn_p, sn_e, sn_vol, sn_d, sn_ed
  real(dp):: rr,pgas,dgas,ekin,mass_sn_tot,mass_sn_tot_all
  integer:: info
  real(dp),dimension(1:nvector,1:ndim)::x
  real(dp),dimension(1:3):: xshift, x_sn
  real(dp), save:: xseed
  real(dp) ::dens_max_loc,dens_max_loc_all
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  logical, dimension(1:nstellarmax):: mark_del
  integer:: istellar,isink
  real(dp)::T_sn,sn_ed_lim,pnorm_sn,vol_sn,vol_rap
  real(dp)::pgas_check,pgas_check_all
  integer, parameter:: navg = 3
  integer, parameter:: nsph = 1
  real(dp), dimension(1:nsph, 1:3):: avg_center
  real(dp), dimension(1:nsph):: avg_radius
  real(dp), dimension(1:navg):: avg_rpow
  real(dp), dimension(1:navg, 1:nvar+3):: avg_upow
  real(dp), dimension(1:navg, 1:nsph):: avg
  real(dp):: norm, rad_sn

  if(.not. hydro)return
  if(ndim .ne. 3)return

  if(verbose)write(*,*)'Entering make_sn_stellar'

  ! TC: random number for direction should be looked at more carefully
  ! RNG should be seeded first

  ! Mesh spacing in that level
  nx_loc = icoarse_max - icoarse_min + 1
  skip_loc=(/dble(icoarse_min),dble(jcoarse_min),dble(kcoarse_min)/)
  scale = boxlen / dble(nx_loc)

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  sn_r = 3.0d0*(0.5d0**levelmin)*scale
  if(sn_r_sat .ne. 0) sn_r = max(sn_r, sn_r_sat * pc2cm / scale_l) !impose a minimum size of 12 pc for the radius
  sn_m = sn_mass_ref !note this is replaced later
  sn_p = sn_p_ref
  sn_e = sn_e_ref
  sn_vol = 4. / 3. * pi * sn_r**3

  !we loop over stellar objects to determine whether one is turning supernovae
  !after it happens, the object is removed from the list
  mark_del = .false.
  do istellar = 1, nstellar
    if(t - tstellar(istellar) < ltstellar(istellar)) cycle
    mark_del(istellar) = .true.

    ! find correct index in sink array which will be equal or lower than id_sink due to sink merging
    isink = id_stellar(istellar)
    do while (id_stellar(istellar) .ne. idsink(isink))
      isink = isink - 1
    end do

    ! the mass of the massive stars 
    sn_m = mstellar(istellar) 

    !remove the mass that is dumped in the grid
    msink(isink) = msink(isink) - sn_m

    !the velocity dispersion times the life time of the object
    rad_sn = ltstellar(istellar)*Vdisp

    !find a random point within a sphere of radius < 1
    !TC: can this be done more efficiently?
    norm=2.
    do while (norm .gt. 1)
       norm=0.
       do idim = 1, ndim
          call random_number(xseed)
          xshift(idim) = (xseed-0.5)*2.
          norm = norm + xshift(idim)**2
       end do
    end do
    do idim = 1, ndim - 1 
        xshift(idim) = xshift(idim) * rad_sn
    end do
    !special treatment for the z-coordinates to maintain it at low altitude
    xshift(3) = xshift(3) * min(rad_sn,100.)
    !TC: this 100 doesn't have units, while rad_sn is in c.u.

    ! place the supernovae around sink particles
    x_sn(:) = xsink(isink, :) + xshift(:)

    !apply periodic boundary conditions (only along x and y)
    !PH note that this should also be modified for the shearing box 24/01/2017
    ! TODO: check which BC!
    ! if BC -> aply, else: move SN back in the box
    if( x_sn(1) .lt. 0) x_sn(1) = boxlen - x_sn(1) 
    if( x_sn(2) .lt. 0) x_sn(2) = boxlen - x_sn(2) 
    if( x_sn(1) .gt. boxlen) x_sn(1) = - boxlen + x_sn(1) 
    if( x_sn(2) .gt. boxlen) x_sn(2) = - boxlen + x_sn(2) 

    avg_center(1, :) = x_sn(:)
    avg_radius = sn_r
    avg_rpow = 0.0d0
    avg_upow = 0.0d0
    ! avg_rpow(1) = 0 ; avg_upow(1, :) = 0 -> integrand(1) = 1
    ! avg_rpow(2) = 0 ; avg_upow(2, 1) = 1 -> integrand(2) = density
    avg_upow(2, 1) = 1.0d0
    ! avg_rpow(3) = 1 ; avg_upow(3, :) = 0 -> integrand(3) = radius
    avg_rpow(3) = 1.0d0
    call sphere_average(navg, nsph, avg_center, avg_radius, avg_rpow, avg_upow, avg)
    vol_sn = avg(1, 1)
    pnorm_sn = avg(3, 1)

    !compute energy and mass density
    sn_d = sn_m / vol_sn
    sn_ed = sn_e / vol_sn

    dens_max_loc = 0.
    mass_sn_tot = 0.
    dens_max_loc_all = 0.
    mass_sn_tot_all = 0.
    pgas_check=0.

    !now loop over cells again and damp energies, mass and momentum
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
!                rr = sqrt(sum(((xx(i,:) - x_sn(:)) / sn_r)**2))
                rr = 0.
                do idim=1,ndim
                   rr = rr + ((xx(i,idim) - x_sn(idim)) / sn_r)**2
                enddo
                rr = sqrt(rr)

                if(rr < 1.) then
                  uold(ind_cell(i), 1) = uold(ind_cell(i), 1) + sn_d
                  dgas = uold(ind_cell(i), 1)

                  if(dgas .gt. dens_max_loc) dens_max_loc = dgas

                  mass_sn_tot = mass_sn_tot + dgas*vol_loc


                  !compute velocity of the gas within this cell assuming 
                  !energy equipartition
                  pgas = min(sn_p / pnorm_sn * rr /  dgas , Vsat) * dgas
                  pgas_check = pgas_check + pgas * vol_loc

                  ekin=0.
                  do idim=1,ndim
                     ekin = ekin + ( (uold(ind_cell(i),idim+1))**2 ) / dgas / 2.
                  enddo

                  uold(ind_cell(i), 2+ndim) = uold(ind_cell(i), 2+ndim) - ekin

                  do idim=1,ndim
                     uold(ind_cell(i),idim+1) = uold(ind_cell(i),idim+1) + pgas * (xx(i,idim) - x_sn(idim)) / (rr * sn_r)
                  enddo

                  ekin=0.
                  do idim=1,ndim
                     ekin = ekin + ( (uold(ind_cell(i),idim+1))**2 ) / dgas / 2.
                  enddo

                  !before adding thermal energy make sure the temperature is not too high (too small timesteps otherwise)
                  T_sn = (sn_ed / dgas * (gamma-1.) ) * scale_T2
                  T_sn = min( T_sn , Tsat) / scale_T2
                  sn_ed_lim = T_sn * dgas / (gamma-1.)

                  uold(ind_cell(i), 2+ndim) = uold(ind_cell(i), 2+ndim) + ekin + sn_ed_lim

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
    call MPI_ALLREDUCE(mass_sn_tot,mass_sn_tot_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    call MPI_ALLREDUCE(dens_max_loc,dens_max_loc_all,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
    call MPI_ALLREDUCE(pgas_check,pgas_check_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#else
    mass_sn_tot_all = mass_sn_tot
    dens_max_loc_all = dens_max_loc
    pgas_check_all = pgas_check
#endif

    if(myid == 1) write(*, *) "SN event: momentum (injected, expected)=", pgas_check_all, sn_p
    if(myid == 1) write(*, *) "Physical units:", pgas_check_all * scale_d * scale_l**3 * scale_v, sn_p * scale_d * scale_l**3 * scale_v

    !calculate grid effect
    vol_rap = vol_sn / sn_vol

    !TC: should be outputted to the log
    if(myid .eq. 1) then 
       open(103,file='supernovae2.txt',form='formatted',status='unknown',access='append')
         write(103,112) t,x_sn(1),x_sn(2),x_sn(3),dens_max_loc_all,mass_sn_tot_all,vol_rap,pgas_check_all,sn_p
       close(103)
    endif

112 format(9e12.4)

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
subroutine sphere_average(navg, nsph, center, radius, rpow, upow, avg)
    use amr_parameters, only: boxlen, dp, hydro, icoarse_max, icoarse_min &
        & , jcoarse_min, kcoarse_min, levelmin, ndim, ngridmax, nlevelmax &
        & , nvector, twotondim, verbose
    use amr_commons, only: active, ncoarse, son, xg, myid
    use hydro_parameters, only: nvar
    use hydro_commons, only: uold
    use mpi_mod
    implicit none

    ! Integrate quantities over spheres
    ! The integrand is (r / radius(isph))**rpow(iavg) * product(u(ivar)**upow(iavg, ivar), ivar=1:nvar)

    integer, intent(in):: navg                               ! Number of quantities
    integer, intent(in):: nsph                               ! Number of spheres
    real(dp), dimension(1:nsph, 1:ndim), intent(in):: center ! Sphere centers
    real(dp), dimension(1:nsph), intent(in):: radius         ! Sphere radii
    real(dp), dimension(1:navg), intent(in):: rpow           ! Power of radius in the integral
#ifdef SOLVERmhd
    real(dp), dimension(1:navg, 1:nvar+3), intent(in):: upow ! Power of hydro variables in the integral
#else
    real(dp), dimension(1:navg, 1:nvar), intent(in):: upow   ! Power of hydro variables in the integral
#endif
    real(dp), dimension(1:navg, 1:nsph), intent(out):: avg   ! Averages

    integer:: i, ivar, ilevel, igrid, ind, ix, iy, iz, iskip, isph, idim
    integer:: nx_loc, ncache, ngrid
    integer, dimension(1:nvector):: ind_grid, ind_cell
    logical, dimension(1:nvector):: ok

    real(dp):: scale, dx, dx_loc, vol_loc, rr
    real(dp), dimension(1:3):: skip_loc
    real(dp), dimension(1:twotondim, 1:3):: xc
    real(dp), dimension(1:nvector, 1:ndim):: xx

    integer:: info
    real(dp), dimension(1:navg, 1:nsph):: avg_loc
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
                        do isph = 1, nsph
                            rr = sqrt(sum(((xx(i, :) - center(isph, :)) / radius(isph))**2))

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
                                avg_loc(:, isph) = avg_loc(:, isph) + vol_loc * integrand
                            endif
                            ! End test on radius
                        end do
                        ! End loop over spheres
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
    call MPI_ALLREDUCE(avg_loc, avg, navg * nsph, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
#else
    avg = avg_loc
#endif
end subroutine sphere_average
