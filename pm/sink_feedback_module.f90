!module sink_feedback_parameters
!end module sink_feedback_parameters

module sink_feedback_module
  use amr_parameters,only:dp,ndim
  implicit none

  public


  logical::sn_feedback_sink = .false. !SN feedback emanates from the sink
  logical::sn_feedback_cr=.false.     !Add CR component to the SN feedback

  !mass, energy and momentum of supernova for forcing by sinks
  ! sn_e in erg, sn_p in g cm/s , sn_mass in g
  real(dp):: sn_e=1.d51 , sn_p=4.d43 , sn_mass=2.e33
  real(dp):: sn_e_ref , sn_p_ref , sn_mass_ref

  !limit speed and temperatue in supernova remnants, minimum radius for the SN remnant
  real(dp):: Tsat=1.d99 , Vsat=1.d99, sn_r_sat=0.

  !dispersion velocity of the stellar objects in km/s
  real(dp):: Vdisp=1. 

  !fraction of the SN energy put in the CR (WARNING, this is currently not set to sn_e at maxiume, so if fcr=1, then the energy inptu can ne 2*sn_e)
  real(dp)::fcr=0.1

  ! Stellar object related arrays, those parameters are read in  read_stellar_params 
  logical:: sn_direct = .false.
  character(LEN=100)::stellar_strategy='local' ! local: create stellar particles from each sink
                                               ! global: create when the total mass in sinks exceeds stellar_msink_th

  integer:: nstellarmax=2000 ! maximum number of stellar objects
  integer:: nstellar = 0 ! current number of stellar objects
  real(dp):: imf_index, imf_low, imf_high ! power-law IMF model: PDF index, lower and higher mass bounds (Msun)
  real(dp):: lt_t0, lt_m0, lt_a, lt_b ! Stellar lifetime model: t(M) = lt_t0 * exp(lt_a * (log(lt_m0 / M))**lt_b)

  ! Stellar ionizing flux model: S(M) = stf_K * (M / stf_m0)**stf_a / (1 + (M / stf_m0)**stf_b)**stf_c
  ! This is a fit from Vacca et al. 1996
  ! Corresponding routine : vaccafit
  real(dp)::stf_K=9.634642584812752d48 ! s**(-1) then normalised in code units in read_stellar
  real(dp)::stf_m0=2.728098824280431d1 ! Msun then normalised in code units in read_stellar
  real(dp)::stf_a=6.840015602892084d0
  real(dp)::stf_b=4.353614230584390d0
  real(dp)::stf_c=1.142166657042991d0 

  !     hii_w: density profile exponent (n = n_0 * (r / r_0)**(-hii_w))
  !     hii_alpha: recombination constant in code units
  !     hii_c: HII region sound speed
  !     hii_t: fiducial HII region lifetime, it is normalised in code units in read_stellar 
  !     hii_T2: HII region temperature
  real(dp):: hii_w, hii_alpha, hii_c, hii_t, hii_T2
  real(dp):: stellar_msink_th ! sink mass threshold for stellar object creation (Msun)
  real(dp), allocatable, dimension(:, :):: xstellar ! stellar object position
  real(dp), allocatable, dimension(:):: mstellar, tstellar, ltstellar ! stellar object mass, birth time, life time
  integer, allocatable, dimension(:):: id_stellar !the id  of the sink to which it belongs
  ! Allow users to pre-set stellar mass selection for physics comparison runs, etc
  ! Every time mstellar is added to, instead of a random value, use mstellarini
  integer,parameter::nstellarini=5000
  real(dp),dimension(nstellarini)::mstellarini ! List of stellar masses to use

  ! Use the supernova module?
  logical::FB_on = .false.

  !series of supernovae specified by "hand"
  ! Number of supernovae (max limit and number active in namelist)
  integer,parameter::NSNMAX=1000
  integer::FB_nsource=0

  ! Type of source ('supernova', 'wind')
  ! NOTE: 'supernova' is a single dump, 'wind' assumes these values are per year
  character(LEN=10),dimension(1:NSNMAX)::FB_sourcetype='supernova'

  ! Feedback start and end times (NOTE: supernova ignores FB_end)
  real(dp),dimension(1:NSNMAX)::FB_start = 1d10
  real(dp),dimension(1:NSNMAX)::FB_end = 1d10
  ! Book-keeping for whether the SN has happened
  logical,dimension(1:NSNMAX)::FB_done = .false.
  
  ! Source position in units from 0 to 1
  real(dp),dimension(1:NSNMAX)::FB_pos_x = 0.5d0
  real(dp),dimension(1:NSNMAX)::FB_pos_y = 0.5d0
  real(dp),dimension(1:NSNMAX)::FB_pos_z = 0.5d0

  ! Ejecta mass in solar masses (/year for winds)
  real(dp),dimension(1:NSNMAX)::FB_mejecta = 1.d0

  ! Energy of source in ergs (/year for winds)
  real(dp),dimension(1:NSNMAX)::FB_energy = 1.d51

  ! Use a thermal dump? (Otherwise add kinetic energy)
  ! Note that if FB_radius is 0 it's thermal anyway
  logical,dimension(1:NSNMAX)::FB_thermal = .false.

  ! Radius to deposit energy inside in number of cells (at highest level)
  real(dp),dimension(1:NSNMAX)::FB_radius = 12d0

CONTAINS


subroutine vaccafit(M,S)
  ! This routine is called in sink_RT_feedback
  ! perform a fit of the Vacca et al. 96 ionising flux
  ! M - stellar mass / solar masses
  ! S - photon emission rate in / s

  real(dp),intent(in)::M
  real(dp),intent(out)::S
  
  S = stf_K * (M / stf_m0)**stf_a / (1. + (M / stf_m0)**stf_b)**stf_c

end subroutine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine make_sn_stellar
  use pm_commons
  use amr_commons
  use hydro_commons
  use mpi_mod
  implicit none

  integer:: ivar
  integer:: ilevel, ind, ix, iy, iz, ngrid, iskip, idim
  integer:: i, nx_loc, igrid, ncache
  integer, dimension(1:nvector), save:: ind_grid, ind_cell
  real(dp):: dx
  real(dp):: scale, dx_min, dx_loc, vol_loc
  real(dp), dimension(1:3):: xbound, skip_loc
  real(dp), dimension(1:twotondim, 1:3):: xc
  logical, dimension(1:nvector), save:: ok

  real(dp), dimension(1:nvector, 1:ndim), save:: xx
  real(dp):: sn_r, sn_m, sn_p, sn_e, sn_vol, sn_d, sn_ed, sn_rp
  real(dp):: rr, pi,dens_max,pgas,dgas,ekin,mass_sn_tot,dens_max_all,mass_sn_tot_all
  integer:: n_sn,n_sn_all,info
  integer ,dimension(1:nvector)::cc
  real(dp),dimension(1:nvector,1:ndim)::x
  real(dp),dimension(1:3):: xshift, x_sn

  logical, save:: first = .true.
  real(dp), save:: xseed

  real(dp) ::dens_max_loc,dens_max_loc_all

  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m, pc

  logical, dimension(1:nstellarmax):: mark_del
  integer:: istellar

  real(dp)::T_sn,sn_ed_lim,pnorm_sn,pnorm_sn_all,vol_sn,vol_sn_all,vol_rap, mass_sn, mass_sn_all,dens_moy

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

  pi = acos(-1.0)

  if (first) then
     xseed = 0.5
     call random_number(xseed)
     first = .false.
  endif

  ! Mesh spacing in that level
  xbound(1:3) = (/ dble(nx), dble(ny), dble(nz) /)
  nx_loc = icoarse_max - icoarse_min + 1
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale = boxlen / dble(nx_loc)
  dx_min = scale * 0.5d0**nlevelmax

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_m=scale_d*(scale_l**3)
  pc = 3.08d18 / scale_l

  sn_r = 3.0d0*(0.5d0**levelmin)*scale
  if(sn_r_sat .ne. 0) sn_r = max(sn_r, sn_r_sat * pc) !impose a minimum size of 12 pc for the radius
!  sn_r = 2.*(0.5**levelmin)*scale
  sn_m = sn_mass_ref !note this is replaced later
  sn_p = sn_p_ref
  sn_e = sn_e_ref
  sn_rp = 0.

    sn_vol = 4. / 3. * pi * sn_r**3

  !we loop over stellar objets to determine whether one is turning supernovae
  !after it happens, the object is removed from the list
  mark_del = .false.
  do istellar = 1, nstellar
    if(t - tstellar(istellar) < ltstellar(istellar)) cycle
    mark_del(istellar) = .true.

    ! the mass of the massive stars 
    sn_m = mstellar(istellar) 

    !remove the mass that is dumped in the grid
    msink(id_stellar(istellar)) = msink(id_stellar(istellar)) - sn_m

    !the velocity dispersion times the life time of the object
    rad_sn = ltstellar(istellar)*Vdisp

    !find a random point within a sphere of radius < 1
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

!    x_sn(:) = xstellar(istellar, :) + xshift(:)
    !place the supernovae around sink particles
    x_sn(:) = xsink(id_stellar(istellar), :) + xshift(:)

    !apply periodic boundary conditions (only along x and y)
    !PH note that this should also be modified for the shearing box 24/01/2017
    if( x_sn(1) .lt. 0) x_sn(1) = boxlen - x_sn(1) 
    if( x_sn(2) .lt. 0) x_sn(2) = boxlen - x_sn(2) 
    if( x_sn(1) .gt. boxlen) x_sn(1) = - boxlen + x_sn(1) 
    if( x_sn(2) .gt. boxlen) x_sn(2) = - boxlen + x_sn(2) 



    if(.true.) then
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
      mass_sn = avg(2, 1) + vol_sn * sn_d ! region average + ejecta
      pnorm_sn = avg(3, 1)
    else
    !do a first path to compute the volume of the cells that are enclosed in the supernovae radius
    !this is to correct for the grid effects
    vol_sn = 0. ; vol_sn_all=0.
    mass_sn = 0. ; mass_sn_all=0.
    pnorm_sn = 0. ; pnorm_sn_all=0.

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
        if(ndim>0) xc(ind,1) = (dble(ix) - 0.5d0) * dx
        if(ndim>1) xc(ind,2) = (dble(iy) - 0.5d0) * dx
        if(ndim>2) xc(ind,3) = (dble(iz) - 0.5d0) * dx
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
                   rr = rr + ( (xx(i,idim) - x_sn(idim)) / sn_r)**2
                enddo

                if(rr < 1.) then
                  vol_sn = vol_sn + vol_loc
                  mass_sn = mass_sn + vol_loc * (uold(ind_cell(i), 1) + sn_d)
                  pnorm_sn = pnorm_sn + vol_loc * sqrt(rr)
                endif
            endif
          end do
          !  End loop over sublist of cells
        end do
        ! End loop over cells
      end do
      ! End loop over grids
    end do
    ! End loop over levels

#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(vol_sn,vol_sn_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    call MPI_ALLREDUCE(mass_sn,mass_sn_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    call MPI_ALLREDUCE(pnorm_sn,pnorm_sn_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)

    vol_sn  = vol_sn_all
    mass_sn = mass_sn_all
    pnorm_sn = pnorm_sn_all
#endif
    end if

    !compute energy and mass density
    sn_d = sn_m / vol_sn
    sn_ed = sn_e / vol_sn

    dens_moy = mass_sn / vol_sn

    dens_max_loc = 0.
    mass_sn_tot = 0.
    dens_max_loc_all = 0.
    mass_sn_tot_all = 0.
    n_sn=0.

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
                  !pgas = sqrt(eff_sn*sn_ed / max(dens_moy,dgas) ) * dgas 
                  !for cells where dgas < dens_moy, take the same velocity as if
                  !the density were equal to dens_moy
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
#if NCR>0
                  if(sn_feedback_cr)then
                     do ivar=1,ncr
                        uold(ind_cell(i),firstindex_ent+ivar)=uold(ind_cell(i),firstindex_ent+ivar) &
                             & + sn_ed*fcr
                        uold(ind_cell(i), 2+ndim) = uold(ind_cell(i), 2+ndim) + sn_ed*fcr
                     end do
                  end if
#endif
                  n_sn = n_sn + 1

                  !write(*,*) 'put SN, myid , n_sn ',myid, n_sn, 'x,y,z: ',x_sn(1),x_sn(2),x_sn(3), 'density ',dgas

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
    call MPI_ALLREDUCE(n_sn,n_sn_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
    call MPI_ALLREDUCE(mass_sn_tot,mass_sn_tot_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    call MPI_ALLREDUCE(dens_max_loc,dens_max_loc_all,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
    call MPI_ALLREDUCE(pgas_check,pgas_check_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#else
    n_sn_all = n_sn
    mass_sn_tot_all = mass_sn_tot
    dens_max_loc_all = dens_max_loc
    pgas_check_all = pgas_check
#endif

    if(myid == 1) write(*, *) "SN momentum (injected, expected):", pgas_check_all, sn_p
    if(myid == 1) write(*, *) "Physical units:", pgas_check_all * scale_d * scale_l**3 * scale_v, sn_p * scale_d * scale_l**3 * scale_v

  !write(*,*) '4 n_sn ', n_sn_all

  !  if(myid .eq. cc(1)) write(*,*) '4 n_sn ', n_sn_all

    !calculate grid effect
    vol_rap = vol_sn / sn_vol

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
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

! THIS SECTION DEALS WITH INDIVIDUAL FIXED SOURCES IN THE NAMELIST

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
subroutine feedback_fixed(ilevel)
  use amr_commons
  use hydro_parameters
  implicit none
  integer::ilevel,isn
  !---------------------
  ! Dump a supernova/winds into cells
  ! NOTE: variables internally still labelled "sn", but valid for winds too
  !---------------------
  ! This check should be done already in amr_step, but best to be sure I suppose
  if(FB_on .and. ilevel == levelmin) then
     do isn=1,FB_nsource
        ! Deal with supernovae
        if(FB_sourcetype(isn) .eq. 'supernova') then
          if(t >= FB_start(isn) .and. .not. FB_done(isn)) then
             ! NOTE: FB_done NOT SAVED BY OUTPUT, SO RELAUNCHING AFTER FB_start WILL LAUNCH SUPERNOVA
             if(myid == 1) write (*,*) 'Supernova ',isn,' @ t = ', t, &
                  & 'FB_start =', FB_start(isn)
             call make_fb_fixed(ilevel,isn)
             FB_done(isn) = .true.
          endif
        endif
        ! Deal with winds
        if(FB_sourcetype(isn) .eq. 'wind') then
          if(t >= FB_start(isn) .and. t <= FB_end(isn)) then
            ! Inject winds
            call make_fb_fixed(ilevel,isn)
            ! Log file book-keeping, wind starts
            if(.not. FB_done(isn)) then
              if(myid == 1) write (*,*) 'Wind started ',isn,' @ t = ', t, &
                    & 'FB_start =', FB_start(isn)
              FB_done(isn) = .true.
            endif
          endif
          ! Log file book-keeping, wind ends
          if(t > FB_end(isn) .and. FB_done(isn)) then
            if(myid == 1) write (*,*) 'Wind stopped ',isn,' @ t = ', t, &
                  & 'FB_end =', FB_end(isn)
            FB_done(isn) = .false.
          endif
        endif
     end do
  endif
end subroutine feedback_fixed
!################################################################
!################################################################
!################################################################
!################################################################
subroutine make_fb_fixed(currlevel,isn)
  ! Adapted from O. Iffrig's make_sn_blast
  use amr_commons
  use hydro_commons
  implicit none

  integer, intent(in) :: currlevel,isn

  integer:: ivar
  integer:: ilevel, ind, ix, iy, iz, ngrid, iskip, idim
  integer:: i, nx_loc, igrid, ncache
  integer, dimension(1:nvector), save:: ind_grid, ind_cell
  real(dp):: dx, dt
  real(dp):: scale, dx_min, dx_loc, vol_loc

  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::scale_msun, scale_ecgs

  real(dp), dimension(1:3):: xbound, skip_loc
  real(dp), dimension(1:twotondim, 1:3):: xc
  logical, dimension(1:nvector), save:: ok

  real(dp),dimension(1:ndim):: sn_cent
  real(dp), dimension(1:nvector, 1:ndim), save:: xx
  real(dp):: sn_r, sn_m, sn_e, sn_vol, sn_d, sn_ed, dx_sel, sn_p, sn_v
  real(dp):: rr, pi
  real(dp), dimension(1:ndim)::rvec
  logical:: sel = .false.
  real(dp),parameter::m_sun=1.9891d33  ! Solar mass [g]
  real(dp),parameter::year2=3.154d7  ! 1 year [s]

  if(.not. hydro)return
  if(ndim .ne. 3)return

  if(verbose)write(*,*)'Entering make_fb_fixed'

  pi = acos(-1.0)

  ! Mesh spacing in that level
  xbound(1:3) = (/ dble(nx), dble(ny), dble(nz) /)
  nx_loc = icoarse_max - icoarse_min + 1
  skip_loc = (/ 0.0d0, 0.0d0, 0.0d0 /)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale = boxlen / dble(nx_loc)
  dx_min = scale * 0.5d0**nlevelmax

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_msun = scale_d * scale_l**ndim / m_sun
  scale_ecgs = scale_d * scale_v**2 * scale_l**ndim

  ! Hard-code sn properties to centre of box
  sn_r = FB_radius(isn)*(0.5**nlevelmax)*scale

  !Set up injection energy, mass, velocity and position
  sn_m = FB_mejecta(isn) / scale_msun ! Put in 10 solar masses
  sn_e = FB_energy(isn) / scale_ecgs
  sn_v = sqrt(2.0*(sn_e*scale_ecgs)/(sn_m*scale_msun*m_sun))
  sn_v = sn_v / scale_v
  sn_cent(1)= FB_pos_x(isn)*boxlen
#if NDIM>1
  sn_cent(2)= FB_pos_y(isn)*boxlen
#endif
#if NDIM>2
  sn_cent(3)= FB_pos_z(isn)*boxlen
#endif

  ! HACK - force Courant condition on winds before first hydro step
  ! TODO: think about how this is synchronised more carefully
  dt = min(dtnew(currlevel),dt)
  dt = dtnew(currlevel)

  ! If this is a wind, scale the luminosity by the timestep
  if (FB_sourcetype(isn) .eq. 'wind') then
    sn_m = sn_m * dt*scale_t/year2 ! dt in years
    sn_e = sn_e * dt*scale_t/year2
  endif

  ! HACK !!! - KINETIC BLAST ONLY WORKS FOR sn_r > 0.0 !!!
  if(sn_r /= 0.0) then
     sn_vol = 4. / 3. * pi * sn_r**3
     sn_d = sn_m / sn_vol
     sn_ed = sn_e / sn_vol
     sn_p = sn_d*sn_v ! uniform momentum of blast ejecta
  end if
     
  if(myid .eq. 1 .and. FB_sourcetype(isn) .eq. 'supernova') then
     write(*,*) 'Supernova blast! Wow!'
#if NDIM==3
     write(*,*) 'x_sn, y_sn, z_sn, ',sn_cent(1),sn_cent(2),sn_cent(3)
#elif NDIM==2
     write(*,*) 'x_sn, y_sn, ',sn_cent(1),sn_cent(2)
#elif NDIM==1
     write(*,*) 'x_sn, ',sn_cent(1)
#endif
  endif

  ! Loop over levels
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
      if(ndim>0) xc(ind,1) = (dble(ix) - 0.5d0) * dx
      if(ndim>1) xc(ind,2) = (dble(iy) - 0.5d0) * dx
      if(ndim>2) xc(ind,3) = (dble(iz) - 0.5d0) * dx
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
            if(sn_r == 0.0) then
              sn_d = sn_m / vol_loc ! XXX
              sn_ed = sn_e / vol_loc ! XXX
              rr = 1.0
              do idim = 1, ndim
                !rr = rr * max(1.0 - abs(xx(i, idim) - sn_center(sn_i, idim)) / dx_sel, 0.0)
                rr = rr * max(1.0 - abs(xx(i, idim) - sn_cent(idim)) / dx_loc, 0.0)
              end do
              !if(rr > 0.0) then
                !if(.not. sel) then
                  !! We found a leaf cell near the supernova center
                  !sel = .true.
                  !sn_d = sn_m / sn_vol
                  !sn_ed = sn_e / sn_vol
                !end if
                uold(ind_cell(i), 1) = uold(ind_cell(i), 1) + sn_d * rr
                uold(ind_cell(i), 2+ndim) = uold(ind_cell(i), 2+ndim) + sn_ed * rr
              !end if
            else
               ! Get direction to point the explosion in
               do idim = 1, ndim
                  rvec(idim) = (xx(i, idim) - sn_cent(idim)) / sn_r
               enddo
               rr = sqrt(sum(rvec**2))
               rvec = rvec/rr
   
               if(rr < 1.) then
                  uold(ind_cell(i), 1) = uold(ind_cell(i), 1) + sn_d
                  ! If not entirely thermal injection, add some velocity
                  if (.not.FB_thermal(isn)) then
                     do idim=1,ndim 
                        uold(ind_cell(i),1+idim) = uold(ind_cell(i),1+idim) + &
                          & sn_p * rvec(idim)
                     enddo
                  end if
                  uold(ind_cell(i), 2+ndim) = uold(ind_cell(i), 2+ndim) + sn_ed
               endif
            endif
          endif
        end do
      end do
      ! End loop over cells
    end do
    ! End loop over grids
  end do
  ! End loop over levels

  ! Update hydro quantities for split cells
  do ilevel = nlevelmax, levelmin, -1
    call upload_fine(ilevel)
    do ivar = 1, nvar
      call make_virtual_fine_dp(uold(1, ivar), ilevel)
    enddo
  enddo
end subroutine make_fb_fixed

END MODULE sink_feedback_module
