
module deltaE_module

  use amr_commons

  implicit none

  ! DeltaE params
  logical::deltaE_enable=.true. ! Whether to enable deltaE output
  logical::deltaE_correct_pressure_fix=.true. ! Whether to be "pressure_fix" aware when computing energy
  logical::deltaE_debug=.false. ! print intermediate energies
  logical::deltaE_force_all_levels=.false. ! Force computation on all levels for all processes
  integer::deltaE_level_turb=0 ! Level at wich the turbulence energy is computed (0 = compute on leaf cells.)

  ! Arrays
  integer, parameter :: nb_energy_kind = 13
  integer, parameter :: iekin = 1, iepot = 2, ieint = 3, iemag = 4, iekin_gas = 5, &
                        & iekin_part = 6, iepot_gas = 7, iepot_part = 8, &
                        & iekin_gas_turb = 9, iekin_part_star = 10, &
                        & iekin_part_dm = 11, iepot_part_star = 12, iepot_part_dm = 13
  character(len=14), dimension(1:nb_energy_kind) :: energy_names = [character(len=14) :: &
                        & "ekin", "epot", "eint", &
                        & "emag", "ekin_gas", "ekin_part",&
                        & "epot_gas", "epot_part", &
                        & "ekin_gas_turb", "ekin_part_star", "ekin_part_dm", &
                        & "epot_part_star", "epot_part_dm"]

  type process
    real(dp), dimension(1:nb_energy_kind) :: v
    character(len=20) :: name
  end type 

  type processes
    type(process) :: flux_gas, flux_part, cooling, gravity_gas, gravity_part, &
                      &  star_formation, feedback, turb_driving, magnetic_diffusion, corrections
  contains
    procedure :: initialize_processes
    procedure :: print_processes
  end type

  type(processes) :: deltaE = processes( &
    process(0.0d0, "flux_gas"), &
    process(0.0d0, "flux_part"), &
    process(0.0d0, "cooling"), &
    process(0.0d0, "gravity_gas"), &
    process(0.0d0, "gravity_part"), &
    process(0.0d0, "star_formation"), &
    process(0.0d0, "feedback"), &
    process(0.0d0, "turb_driving"), &
    process(0.0d0, "magnetic_diffusion"), &
    process(0.0d0, "corrections"))

contains

  subroutine read_deltaE_params(namelist_unit,nml_ok)
   
    integer,intent(in)::namelist_unit
    logical,intent(inout)::nml_ok
    integer::nml_err
   
    namelist/deltaE_params/deltaE_correct_pressure_fix,deltaE_debug,deltaE_enable,deltaE_force_all_levels,deltaE_level_turb
   ! Go to the beginning of the file
    rewind(namelist_unit)

    ! Read namelist
    read(namelist_unit,NML=deltaE_params,IOSTAT=nml_err)

    if(nml_err>0)then
      if(myid==1)write(*,*)'Error reading namelist &DELTAE_PARAMS. Check formatting.'
      nml_ok=.false.
    end if

    call deltaE%initialize_processes

  end subroutine read_deltaE_params

  subroutine initialize_processes(this)

    implicit none
    class(processes), intent(inout) :: this
   
    this%flux_gas%v = 0.0d0
    this%flux_part%v = 0.0d0
    this%cooling%v = 0.0d0
    this%gravity_gas%v = 0.0d0
    this%gravity_part%v = 0.0d0
    this%star_formation%v = 0.0d0
    this%feedback%v = 0.0d0
    this%turb_driving%v = 0.0d0
    this%magnetic_diffusion%v = 0.0d0
    this%corrections%v = 0d0

  end subroutine

  subroutine print_energies(energies)
    implicit none
    real(dp), dimension(1:nb_energy_kind), intent(in) :: energies
    integer :: i

    write (*, *) "Energy output"
    do i = 1, nb_energy_kind
      write (*, 998) " ", energy_names(i), energies(i)
    end do
    write (*, *) "End energy output"

998 format(A1, A14, 1pe20.7)

  end subroutine

  subroutine print_processes(this)

    implicit none
    class(processes), intent(in) :: this
    integer :: i

    write (*, *) "DeltaE output"
    write (*, 997) " ", this%flux_gas%name, this%flux_part%name, this%cooling%name, this%gravity_gas%name, &
                        this%gravity_part%name, this%star_formation%name, this%feedback%name,  this%turb_driving%name, this%magnetic_diffusion%name, &
                        & this%corrections%name
    do i = 1, nb_energy_kind
      write (*, 998) energy_names(i), this%flux_gas%v(i), this%flux_part%v(i), this%cooling%v(i), this%gravity_gas%v(i), &
                   this%gravity_part%v(i), this%star_formation%v(i), this%feedback%v(i), this%turb_driving%v(i),  this%magnetic_diffusion%v(i), &
                   & this%corrections%v(i)
    end do
    write (*, *) "End deltaE output"

997 format(A20, *(A20))
998 format(A20, *(1pe20.5))

  end subroutine

  subroutine compute_energies(ilevel, use_unew, energies)
    implicit none

    integer, intent(in):: ilevel
    logical, intent(in):: use_unew
    real(dp), dimension(1:nb_energy_kind), intent(inout) :: energies

    energies = 0.0d0
    if (hydro) then
      call compute_energy_gas(ilevel, energies, use_unew)
    end if
    if (pic) then
      !call make_tree_fine(ilevel)
      call compute_energy_part(ilevel, energies)
    end if
    energies(iekin) = energies(iekin_gas) + energies(iekin_part)
    energies(iepot) = energies(iepot_gas) +  energies(iepot_part)

  end subroutine

  subroutine compute_transfer(levelstart, levelend, use_unew, deltaE_process, step)

    implicit none
    integer, intent(in):: levelstart, levelend
    integer :: act_levelstart, act_levelend
    logical, intent(in):: use_unew
    type(process), intent(inout) :: deltaE_process
    integer, intent(in) :: step

    integer :: ilevel

    real(dp), dimension(1:nb_energy_kind), save :: energy_before
    real(dp), dimension(1:nb_energy_kind) ::energy_level, energy_after

    if (deltaE_force_all_levels) then
      act_levelstart = levelmin
      act_levelend = nlevelmax
    else 
      act_levelstart = levelstart
      act_levelend = levelend
    end if

    if (step == 1) then
      energy_before = 0.0d0
    else
      energy_after = 0.0d0
    end if

    do ilevel = act_levelstart, act_levelend
      call compute_energies(ilevel, use_unew, energy_level)
      if (step == 1) then
        energy_before = energy_before + energy_level
      else
        energy_after = energy_after + energy_level
      end if
    end do


    if (step == 2) then

      if (isnan(energy_after(iepot_part)) .or. isnan(energy_before(iepot_part))) then
        if (myid == 1) write(*,*) "DeltaE Warning: nan epot found in ", deltaE_process%name
        energy_after(iepot_part) = 0.0d0
        energy_before(iepot_part) = 0.0d0
        energy_after(iepot) = energy_after(iepot_gas) 
        energy_before(iepot) = energy_before(iepot_gas) 
      end if


      if (isnan(energy_after(iepot_gas)) .or. isnan(energy_before(iepot_gas))) then
        if (myid == 1) write(*,*) "DeltaE Warning: nan epot found in ", deltaE_process%name
        energy_after(iepot_gas) = 0.0d0
        energy_before(iepot_gas) = 0.0d0
        energy_after(iepot) = energy_after(iepot_part) 
        energy_before(iepot) = energy_before(iepot_part) 
      end if

      deltaE_process%v = deltaE_process%v + energy_after - energy_before

      if (deltaE_debug .and. myid == 1) then
        write(*,*) "DeltaE debug before ", deltaE_process%name
        call print_energies(energy_before)
        write(*,*) "DeltaE debug after ", deltaE_process%name
        call print_energies(energy_after)
      end if
    end if

  end subroutine

  subroutine compute_energy_gas(ilevel, energies, use_unew)
    use amr_commons
    use hydro_commons
    use poisson_commons
    use mpi_mod
    implicit none
#ifndef WITHOUTMPI
    integer::info
    real(kind=8), dimension(4)::comm_buffin, comm_buffout
#endif
    integer, intent(in)::ilevel
    real(dp), dimension(1:nb_energy_kind), intent(inout) :: energies
    logical, intent(in)::use_unew

    integer::i, ivar,  ind, ncache, igrid, iskip
    integer::nleaf, ngrid, nx_loc
    integer, dimension(1:nvector), save::ind_grid, ind_cell, ind_leaf

    real(dp):: dx, vol, scale
    real(kind=8)::mass_loc, ekin_loc, eint_loc, emag_loc, epot_loc, ekin_leaf, emag_leaf, epot_leaf
    real(kind=8)::mass_all
    real(dp), dimension(1:nvector, 1:nvar_all), save::uu
    real(dp), dimension(1:nvector, 1:ndim), save::gg

    real(dp):: e_cons, e_prim, e_trunc, div

    mass_loc = 0.0d0
    ekin_loc = 0.0d0; ekin_leaf = 0.0d0
    emag_loc = 0.0d0; emag_leaf = 0.0d0
    eint_loc = 0.0d0
    epot_loc = 0.0d0; epot_leaf = 0.0d0

    energies(iekin_gas) = 0.0d0
    energies(ieint) = 0.0d0
    energies(iemag) = 0.0d0
    energies(iepot) = 0.0d0

    if (numbtot(1, ilevel) == 0) return

    ! Mesh spacing at that level
    nx_loc = icoarse_max - icoarse_min + 1
    scale = boxlen/dble(nx_loc)
    dx = 0.5D0**ilevel*scale
    vol = dx**ndim

    ! Loop over active grids by vector sweeps
    ncache = active(ilevel)%ngrid
    do igrid = 1, ncache, nvector
      ngrid = MIN(nvector, ncache - igrid + 1)
      do i = 1, ngrid
        ind_grid(i) = active(ilevel)%igrid(igrid + i - 1)
      end do

      ! Loop over cells
      do ind = 1, twotondim
        iskip = ncoarse + (ind - 1)*ngridmax
        do i = 1, ngrid
          ind_cell(i) = ind_grid(i) + iskip
        end do

        ! Gather leaf cells
        nleaf = 0
        do i = 1, ngrid
          if (son(ind_cell(i)) == 0) then
            nleaf = nleaf + 1
            ind_leaf(nleaf) = ind_cell(i)
          end if
        end do

        ! Gather hydro variables
        do ivar = 1, nvar_all
          if (use_unew) then
            do i = 1, nleaf
              uu(i, ivar) = unew(ind_leaf(i), ivar)
            end do
          else
            do i = 1, nleaf
              uu(i, ivar) = uold(ind_leaf(i), ivar)
            end do
          end if
        end do

        if (poisson) then 
          do i = 1, nleaf
            epot_leaf = 0.5*vol*uu(i, 1)*phi_old(ind_leaf(i))
            epot_loc = epot_loc + epot_leaf
          end do
        end if

        ! Compute total internal energy step 1
        do i = 1, nleaf
          eint_loc = eint_loc + uu(i, neul)*vol
        end do

        ! Compute total energies
        do ivar = 1, ndim
          do i = 1, nleaf
            ! Compute total kinetic energy
            ekin_leaf = 0.5d0*vol*uu(i, 1 + ivar)**2/uu(i, 1)
            ekin_loc = ekin_loc + ekin_leaf
#ifdef SOLVERmhd
            ! Compute total magnetic energy
            emag_leaf = 0.125d0*(uu(i, neul + ivar) + uu(i, nvar + ivar))**2*vol
            emag_loc = emag_loc + emag_leaf
#endif

            ! Compute total internal energy step 2
            eint_loc = eint_loc - ekin_leaf - emag_leaf
          end do
        end do

        ! Compute total internal energy step 3 
#if NENER>0
        do ivar = 1, nener
          do i = 1, nleaf
            eint_loc = eint_loc - uu(i, nhydro + ivar)*vol
          end do
        end do
#endif

        if(pressure_fix .and. deltaE_correct_pressure_fix)then
          ! Correct internal energy if too small
          do i=1, nleaf
            ekin_leaf = 0.
            do ivar = 1, ndim
              ekin_leaf = ekin_leaf +  0.5d0*uu(i, 1 + ivar)**2/uu(i, 1)
            end do
#if NENER>0
            do ivar=1,nener
              ekin_leaf = ekin_leaf + uu(i, nhydro + ivar)
            end do
#endif
#ifdef SOLVERmhd
            do ivar = 1, ndim
              ekin_leaf = ekin_leaf + 0.125d0*(uu(i, neul + ivar) + uu(i, nvar + ivar))**2
            end do
#endif

            e_cons = uu(i, neul) - ekin_leaf
            e_prim = enew(ind_leaf(i))
            ! Note: here divu=-div.u*dt
            div = abs(divu(ind_leaf(i)))*dx/dtnew(ilevel)
            e_trunc = beta_fix*uu(i, 1)*max(div,3.0d0*hexp*dx)**2
            if(e_cons<e_trunc)then
                eint_loc = eint_loc - e_cons*vol + e_prim*vol
            end if
          end do
        end if
      end do
      ! End loop over cells

    end do
    ! End loop over grids

    ! Compute global quantities
#ifndef WITHOUTMPI
    comm_buffin(1) = ekin_loc
    comm_buffin(2) = eint_loc
    comm_buffin(3) = emag_loc
    comm_buffin(4) = epot_loc
    call MPI_ALLREDUCE(comm_buffin, comm_buffout, 4, MPI_DOUBLE_PRECISION, MPI_SUM,&
            &MPI_COMM_WORLD, info)
    energies(iekin_gas) = comm_buffout(1)
    energies(ieint) = comm_buffout(2)
    energies(iemag) = comm_buffout(3)
    energies(iepot_gas) = comm_buffout(4)
#endif
#ifdef WITHOUTMPI
    energies(iekin_gas) = ekin_loc
    energies(ieint) = eint_loc
    energies(iemag) = emag_loc
    energies(iepot_gas) = epot_loc
#endif

    if (deltaE_level_turb == 0 .or. deltaE_level_turb == ilevel) then
      if (use_unew) then
        energies(iekin_gas_turb) = compute_ekin_turb(ilevel, unew)
      else 
        energies(iekin_gas_turb) = compute_ekin_turb(ilevel, uold)
      end if
    end if

  end subroutine compute_energy_gas

  function compute_ekin_turb(ilevel, uarray) result(ekin_turb)
    use amr_commons
    use hydro_commons
    use mpi_mod
    use amr_constants, only: i1min, i1max, j1min, j1max, k1min, k1max

    implicit none

    ! Parameters and outputs
    integer, intent(in) :: ilevel
    real(dp),dimension(:,:), intent(in):: uarray
    real(dp) :: ekin_turb

    ! AMR Variables
    integer :: i, ind, ncache, igrid, iskip, idim
    integer :: ncell_turb, ngrid, nx_loc
    integer, dimension(nvector) :: ind_grid, ind_cell, ind_cell_turb, index_current_grid
    integer, dimension(1:nvector, 1:threetondim)::nbor_cells
    integer::i1, j1, k1, ind_father

    ! Physical variables
    real(dp) :: scale, dx, vol, ekin_turb_loc
    real(dp), dimension(nvector, ndim):: velocity_mean ! Mean velocity in the grid
    real(dp), dimension(nvector):: mass
    real(dp), dimension(ndim) :: mass_velocity ! Mass* velocity of the current leaf cell
    real(dp) :: vel

    ! MPI variables
#ifndef WITHOUTMPI
    integer::info
#endif

    ! Initialize arrays
    ekin_turb_loc = 0.0d0

    ! Mesh spacing at that level
    nx_loc = icoarse_max - icoarse_min + 1
    scale = boxlen/dble(nx_loc)
    dx = 0.5D0**ilevel*scale
    vol = dx**ndim

    ! Loop over active grids by vector sweeps
    ncache = active(ilevel)%ngrid
    do igrid = 1, ncache, nvector

      velocity_mean = 0.0d0
      mass = 0.0d0

      ngrid = MIN(nvector, ncache - igrid + 1)

      do i = 1, ngrid
        ind_grid(i) = active(ilevel)%igrid(igrid + i - 1)
        ind_cell(i) = father(ind_grid(i)) ! also gather father cells
      end do

      ! Collect neighbor cells
      call get3cubefather(ind_cell, nbor_cells, ngrid, ilevel - 1)

      ! Loop over 3x3x3 neighboring father cells
      do k1 = k1min, k1max
        do j1 = j1min, j1max
          do i1 = i1min, i1max

            ! Get neighbor cell index
            ind_father = 1 + i1 + 3*j1 + 9*k1

            do i = 1, ngrid
              ! Get mass and velocity and aggregate
              ! Always use uold for the father cells (TODO do average ??)
              mass_velocity = 2**(ndim)*vol*uold(nbor_cells(i, ind_father), 2:1 + ndim)
              velocity_mean(i, :) = velocity_mean(i, :) + mass_velocity
              mass(i) = mass(i) + 2**(ndim)*vol*uold(nbor_cells(i, ind_father), 1)
            end do
          end do
        end do
      end do

      ! Normalize
      do i = 1, ngrid
        velocity_mean(i, :) = velocity_mean(i, :)/mass(i)
      end do

      ! Loop over cells
      do ind = 1, twotondim
        iskip = ncoarse + (ind - 1)*ngridmax
        do i = 1, ngrid
          ind_cell(i) = ind_grid(i) + iskip
        end do

        ! Gather cells on which the computation will be done (depending on whether deltaE_level_turb is set)
        ncell_turb = 0
        do i = 1, ngrid
          if (son(ind_cell(i)) == 0 .or. deltaE_level_turb == ilevel) then
            ncell_turb = ncell_turb + 1
            ind_cell_turb(ncell_turb) = ind_cell(i)
            index_current_grid(ncell_turb) = i
          end if
        end do

        ! Compute variance
        do i = 1, ncell_turb
          do idim = 1, ndim
            vel = uarray(ind_cell_turb(i), 1 + idim)/max(uarray(ind_cell_turb(i), 1), smallr)
            ekin_turb_loc = ekin_turb_loc + 0.5*(vol*uarray(ind_cell_turb(i), 1)*(vel - velocity_mean(index_current_grid(i), idim))**2)
          end do
        end do
      end do
    end do

#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(ekin_turb_loc, ekin_turb, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
#endif
#ifdef WITHOUTMPI
    ekin_turb = ekin_turb_loc
#endif

  end function

  subroutine compute_energy_part(ilevel, energies)
    use amr_commons
    use hydro_commons
    use mpi_mod
    use pm_commons

    implicit none

    ! Parameters & outputs
    integer, intent(in):: ilevel
    real(dp), dimension(1:nb_energy_kind), intent(inout) :: energies

    real(dp), dimension(-NFAMILIES:NFAMILIES) :: ekin_families_loc
    real(dp), dimension(-NFAMILIES:NFAMILIES) :: epot_families_loc

    real(kind=8)::ekin_part_loc, epot_part_loc

    ! AMR variables
    integer :: igrid, jgrid

    ! Particle variable
    integer::ipart, jpart, next_part, ip, npart1, local_counter

    integer, dimension(1:nvector), save:: ind_part

    ! MPI variables
#ifndef WITHOUTMPI
    integer::info
#endif

    energies(iepot_part) = 0.d0
    energies(iepot_part_dm) = 0.d0
    energies(iepot_part_star) = 0.d0
    energies(iekin_part) = 0.d0
    energies(iekin_part_dm) = 0.d0
    energies(iekin_part_star) = 0.d0

    ekin_part_loc = 0.d0
    epot_part_loc = 0.d0
    ekin_families_loc = 0.d0
    epot_families_loc = 0.d0

    if (pic) then
      ip = 0
      ! Loop over particles that are not tracers
      do jgrid = 1, active(ilevel)%ngrid
        igrid = active(ilevel)%igrid(jgrid)
        npart1 = numbp(igrid)  ! Number of particles in the grid
        if (npart1 > 0) then
          ipart = headp(igrid)
          ! Loop over particles
          do jpart = 1, npart1
            ! Save next particle  <---- Very important !!!
            next_part = nextp(ipart)
            ! Skip tracers
            if (.not. is_tracer(typep(ipart))) then
              ip = ip + 1
              ind_part(ip) = ipart
              if (ip == nvector) then
                call ekin_part_helper(ind_part, ekin_part_loc, ekin_families_loc, epot_part_loc, epot_families_loc, ip, ilevel)
                ip = 0
              end if
            end if

            ipart = next_part  ! Go to next particle
          end do
          ! End loop over particles
        end if
      end do
      ! End loop over grids
      if (ip > 0) then
        call ekin_part_helper(ind_part, ekin_part_loc, ekin_families_loc, epot_part_loc, epot_families_loc,  ip, ilevel)
      end if

#ifndef WITHOUTMPI
      call MPI_ALLREDUCE(ekin_part_loc, energies(iekin_part), 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
      call MPI_ALLREDUCE(ekin_families_loc(FAM_DM), energies(iekin_part_dm), 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
  call MPI_ALLREDUCE(ekin_families_loc(FAM_STAR), energies(iekin_part_star), 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
      call MPI_ALLREDUCE(epot_part_loc, energies(iepot_part), 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
      call MPI_ALLREDUCE(epot_families_loc(FAM_DM), energies(iepot_part_dm), 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
  call MPI_ALLREDUCE(epot_families_loc(FAM_STAR), energies(iepot_part_star), 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
#endif
#ifdef WITHOUTMPI
      energies(iekin_part) = ekin_part_loc
      energies(iekin_part_dm) = ekin_families_loc(FAM_DM)
      energies(iekin_part_star) = ekin_families_loc(FAM_STAR)
      energies(iepot_part) = epot_part_loc
      energies(iepot_part_dm) = epot_families_loc(FAM_DM)
      energies(iepot_part_star) = epot_families_loc(FAM_STAR)
#endif

    end if

  end subroutine compute_energy_part

  subroutine ekin_part_helper(ind_part, ekin_loc, ekin_families_loc, epot_loc, epot_families_loc, nn, ilevel)
    use amr_commons
    use pm_commons
    use hydro_commons
    implicit none

    integer::nn, ilevel
    integer, dimension(1:nvector)::ind_part

    integer::i, idim, nx_loc
    real(dp)::scale
    real(dp), dimension(1:nvector), save:: mass_part, phi_part
    integer, dimension(1:nvector) :: type_part
    real(dp), dimension(1:nvector, 1:ndim)::vel_part
    real(dp), dimension(-NFAMILIES:NFAMILIES), intent(inout) :: ekin_families_loc
    real(dp), dimension(-NFAMILIES:NFAMILIES), intent(inout) :: epot_families_loc
    real(dp), intent(inout):: ekin_loc, epot_loc

    real(dp)::dx

    ! Compute time step
    dx = 0.5D0**ilevel
    nx_loc = (icoarse_max - icoarse_min + 1)
    scale = boxlen/dble(nx_loc)

    do idim = 1, ndim
      do i = 1, nn
        vel_part(i, idim) = vp(ind_part(i), idim)
      end do
    end do

    ! Fetch mass and type
    do i = 1, nn
      mass_part(i) = mp(ind_part(i))
      type_part(i) = typep(ind_part(i))%family
#ifdef OUTPUT_PARTICLE_POTENTIAL
      phi_part(i) = ptcl_phi(ind_part(i))
      epot_loc = epot_loc + 0.5 * mass_part(i) * phi_part(i)
      epot_families_loc(type_part(i)) = epot_families_loc(type_part(i)) + 0.5 * mass_part(i) *  phi_part(i)
#endif
    end do

    ! Compute kinetic energy
    do idim = 1, ndim
      do i = 1, nn
        ekin_loc = ekin_loc + 0.5D0*mass_part(i)*vel_part(i, idim)**2
        ekin_families_loc(type_part(i)) = ekin_families_loc(type_part(i)) + 0.5D0*mass_part(i)*vel_part(i, idim)**2
      end do
    end do

  end subroutine ekin_part_helper

  !#########################################################################
  !#########################################################################
  !#########################################################################
  !#########################################################################
!   subroutine epot_part_helper(ind_grid, ind_part, ind_grid_part, epot_loc, epot_families_loc, ng, np, ilevel)
!     use amr_commons
!     use pm_commons
!     use poisson_commons
!     use hydro_commons, ONLY: uold, smallr
!     use, intrinsic :: ieee_arithmetic

!     implicit none
!     integer::ng, np, ilevel
!     integer, dimension(1:nvector)::ind_grid
!     integer, dimension(1:nvector)::ind_grid_part, ind_part
!     !------------------------------------------------------------
!     ! This routine computes the potential energy of each particle by
!     ! inverse CIC.
!     ! If particle sits entirely in fine level, then CIC is performed
!     ! at level ilevel. Otherwise, it is performed at level ilevel-1.
!     ! This routine is called by compute_epot_part.
!     !------------------------------------------------------------

!     ! Results
!     real(dp), dimension(-NFAMILIES:NFAMILIES), intent(inout) :: epot_families_loc
!     real(dp), intent(inout):: epot_loc

!     ! Temporary scalar
!     real(dp) :: epot_part
!     integer :: family_part

!     logical::error
!     integer::i, j, ind, idim, nx_loc, isink
!     real(dp)::dx, dx_loc, scale, vol_loc
!     ! Grid-based arrays
!     integer, dimension(1:nvector), save::father_cell
!     real(dp), dimension(1:nvector, 1:ndim), save::x0
!     integer, dimension(1:nvector, 1:threetondim), save::nbors_father_cells
!     ! Particle-based arrays
!     logical, dimension(1:nvector), save::ok
!     real(dp), dimension(1:nvector, 1:ndim), save::x, ff, new_xp, new_vp, dd, dg
!     integer, dimension(1:nvector, 1:ndim), save::ig, id, igg, igd, icg, icd
!     real(dp), dimension(1:nvector, 1:twotondim), save::vol
!     integer, dimension(1:nvector, 1:twotondim), save::igrid, icell, indp, kg
!     real(dp), dimension(1:3)::skip_loc

!     if (isnan(epot_loc)) then 
!       return 
!     end if

!     ! Mesh spacing in that level
!     dx = 0.5D0**ilevel
!     nx_loc = (icoarse_max - icoarse_min + 1)
!     skip_loc = (/0.0d0, 0.0d0, 0.0d0/)
!     if (ndim > 0) skip_loc(1) = dble(icoarse_min)
!     if (ndim > 1) skip_loc(2) = dble(jcoarse_min)
!     if (ndim > 2) skip_loc(3) = dble(kcoarse_min)
!     scale = boxlen/dble(nx_loc)
!     dx_loc = dx*scale
!     vol_loc = dx_loc**3

!     ! Lower left corner of 3x3x3 grid-cube
!     do idim = 1, ndim
!       do i = 1, ng
!         x0(i, idim) = xg(ind_grid(i), idim) - 3.0D0*dx
!       end do
!     end do

!     ! Gather neighboring father cells (should be present anytime !)
!     do i = 1, ng
!       father_cell(i) = father(ind_grid(i))
!     end do
!     call get3cubefather(father_cell, nbors_father_cells, &
!          & ng, ilevel)

!     ! Rescale particle position at level ilevel
!     do idim = 1, ndim
!       do j = 1, np
!         x(j, idim) = xp(ind_part(j), idim)/scale + skip_loc(idim)
!       end do
!     end do
!     do idim = 1, ndim
!       do j = 1, np
!         x(j, idim) = x(j, idim) - x0(ind_grid_part(j), idim)
!       end do
!     end do
!     do idim = 1, ndim
!       do j = 1, np
!         x(j, idim) = x(j, idim)/dx
!       end do
!     end do

!     ! Check for illegal moves
!     error = .false.
!     do idim = 1, ndim
!       do j = 1, np
!         if (x(j, idim) < 0.5D0 .or. x(j, idim) > 5.5D0) error = .true.
!       end do
!     end do
!     if (error) then
!       epot_loc =  ieee_value(epot_loc, ieee_quiet_nan)
!       epot_families_loc(:) = ieee_value(epot_loc, ieee_quiet_nan)
!       return
!     end if

!     ! CIC at level ilevel (dd: right cloud boundary; dg: left cloud boundary)
!     do idim = 1, ndim
!       do j = 1, np
!         dd(j, idim) = x(j, idim) + 0.5D0
!         id(j, idim) = int(dd(j, idim))
!         dd(j, idim) = dd(j, idim) - id(j, idim)
!         dg(j, idim) = 1.0D0 - dd(j, idim)
!         ig(j, idim) = id(j, idim) - 1
!       end do
!     end do

!     ! Compute parent grids
!     do idim = 1, ndim
!       do j = 1, np
!         igg(j, idim) = ig(j, idim)/2
!         igd(j, idim) = id(j, idim)/2
!       end do
!     end do
! #if NDIM==1
!     do j = 1, np
!       kg(j, 1) = 1 + igg(j, 1)
!       kg(j, 2) = 1 + igd(j, 1)
!     end do
! #endif
! #if NDIM==2
!     do j = 1, np
!       kg(j, 1) = 1 + igg(j, 1) + 3*igg(j, 2)
!       kg(j, 2) = 1 + igd(j, 1) + 3*igg(j, 2)
!       kg(j, 3) = 1 + igg(j, 1) + 3*igd(j, 2)
!       kg(j, 4) = 1 + igd(j, 1) + 3*igd(j, 2)
!     end do
! #endif
! #if NDIM==3
!     do j = 1, np
!       kg(j, 1) = 1 + igg(j, 1) + 3*igg(j, 2) + 9*igg(j, 3)
!       kg(j, 2) = 1 + igd(j, 1) + 3*igg(j, 2) + 9*igg(j, 3)
!       kg(j, 3) = 1 + igg(j, 1) + 3*igd(j, 2) + 9*igg(j, 3)
!       kg(j, 4) = 1 + igd(j, 1) + 3*igd(j, 2) + 9*igg(j, 3)
!       kg(j, 5) = 1 + igg(j, 1) + 3*igg(j, 2) + 9*igd(j, 3)
!       kg(j, 6) = 1 + igd(j, 1) + 3*igg(j, 2) + 9*igd(j, 3)
!       kg(j, 7) = 1 + igg(j, 1) + 3*igd(j, 2) + 9*igd(j, 3)
!       kg(j, 8) = 1 + igd(j, 1) + 3*igd(j, 2) + 9*igd(j, 3)
!     end do
! #endif
!     do ind = 1, twotondim
!       do j = 1, np
!         igrid(j, ind) = son(nbors_father_cells(ind_grid_part(j), kg(j, ind)))
!       end do
!     end do

!     ! Check if particles are entirely in level ilevel
!     ok(1:np) = .true.
!     do ind = 1, twotondim
!       do j = 1, np
!         ok(j) = ok(j) .and. igrid(j, ind) > 0
!       end do
!     end do

!     ! If not, rescale position at level ilevel-1
!     do idim = 1, ndim
!       do j = 1, np
!         if (.not. ok(j)) then
!           x(j, idim) = x(j, idim)/2.0D0
!         end if
!       end do
!     end do
!     ! If not, redo CIC at level ilevel-1
!     do idim = 1, ndim
!       do j = 1, np
!         if (.not. ok(j)) then
!           dd(j, idim) = x(j, idim) + 0.5D0
!           id(j, idim) = int(dd(j, idim))
!           dd(j, idim) = dd(j, idim) - id(j, idim)
!           dg(j, idim) = 1.0D0 - dd(j, idim)
!           ig(j, idim) = id(j, idim) - 1
!         end if
!       end do
!     end do

!     ! Compute parent cell position
!     do idim = 1, ndim
!       do j = 1, np
!         if (ok(j)) then
!           icg(j, idim) = ig(j, idim) - 2*igg(j, idim)
!           icd(j, idim) = id(j, idim) - 2*igd(j, idim)
!         else
!           icg(j, idim) = ig(j, idim)
!           icd(j, idim) = id(j, idim)
!         end if
!       end do
!     end do
! #if NDIM==1
!     do j = 1, np
!       icell(j, 1) = 1 + icg(j, 1)
!       icell(j, 2) = 1 + icd(j, 1)
!     end do
! #endif
! #if NDIM==2
!     do j = 1, np
!       if (ok(j)) then
!         icell(j, 1) = 1 + icg(j, 1) + 2*icg(j, 2)
!         icell(j, 2) = 1 + icd(j, 1) + 2*icg(j, 2)
!         icell(j, 3) = 1 + icg(j, 1) + 2*icd(j, 2)
!         icell(j, 4) = 1 + icd(j, 1) + 2*icd(j, 2)
!       else
!         icell(j, 1) = 1 + icg(j, 1) + 3*icg(j, 2)
!         icell(j, 2) = 1 + icd(j, 1) + 3*icg(j, 2)
!         icell(j, 3) = 1 + icg(j, 1) + 3*icd(j, 2)
!         icell(j, 4) = 1 + icd(j, 1) + 3*icd(j, 2)
!       end if
!     end do
! #endif
! #if NDIM==3
!     do j = 1, np
!       if (ok(j)) then
!         icell(j, 1) = 1 + icg(j, 1) + 2*icg(j, 2) + 4*icg(j, 3)
!         icell(j, 2) = 1 + icd(j, 1) + 2*icg(j, 2) + 4*icg(j, 3)
!         icell(j, 3) = 1 + icg(j, 1) + 2*icd(j, 2) + 4*icg(j, 3)
!         icell(j, 4) = 1 + icd(j, 1) + 2*icd(j, 2) + 4*icg(j, 3)
!         icell(j, 5) = 1 + icg(j, 1) + 2*icg(j, 2) + 4*icd(j, 3)
!         icell(j, 6) = 1 + icd(j, 1) + 2*icg(j, 2) + 4*icd(j, 3)
!         icell(j, 7) = 1 + icg(j, 1) + 2*icd(j, 2) + 4*icd(j, 3)
!         icell(j, 8) = 1 + icd(j, 1) + 2*icd(j, 2) + 4*icd(j, 3)
!       else
!         icell(j, 1) = 1 + icg(j, 1) + 3*icg(j, 2) + 9*icg(j, 3)
!         icell(j, 2) = 1 + icd(j, 1) + 3*icg(j, 2) + 9*icg(j, 3)
!         icell(j, 3) = 1 + icg(j, 1) + 3*icd(j, 2) + 9*icg(j, 3)
!         icell(j, 4) = 1 + icd(j, 1) + 3*icd(j, 2) + 9*icg(j, 3)
!         icell(j, 5) = 1 + icg(j, 1) + 3*icg(j, 2) + 9*icd(j, 3)
!         icell(j, 6) = 1 + icd(j, 1) + 3*icg(j, 2) + 9*icd(j, 3)
!         icell(j, 7) = 1 + icg(j, 1) + 3*icd(j, 2) + 9*icd(j, 3)
!         icell(j, 8) = 1 + icd(j, 1) + 3*icd(j, 2) + 9*icd(j, 3)
!       end if
!     end do
! #endif

!     ! Compute parent cell adresses
!     do ind = 1, twotondim
!       do j = 1, np
!         if (ok(j)) then
!           indp(j, ind) = ncoarse + (icell(j, ind) - 1)*ngridmax + igrid(j, ind)
!         else
!           indp(j, ind) = nbors_father_cells(ind_grid_part(j), icell(j, ind))
!         end if
!       end do
!     end do

!     ! Compute cloud volumes
! #if NDIM==1
!     do j = 1, np
!       vol(j, 1) = dg(j, 1)
!       vol(j, 2) = dd(j, 1)
!     end do
! #endif
! #if NDIM==2
!     do j = 1, np
!       vol(j, 1) = dg(j, 1)*dg(j, 2)
!       vol(j, 2) = dd(j, 1)*dg(j, 2)
!       vol(j, 3) = dg(j, 1)*dd(j, 2)
!       vol(j, 4) = dd(j, 1)*dd(j, 2)
!     end do
! #endif
! #if NDIM==3
!     do j = 1, np
!       vol(j, 1) = dg(j, 1)*dg(j, 2)*dg(j, 3)
!       vol(j, 2) = dd(j, 1)*dg(j, 2)*dg(j, 3)
!       vol(j, 3) = dg(j, 1)*dd(j, 2)*dg(j, 3)
!       vol(j, 4) = dd(j, 1)*dd(j, 2)*dg(j, 3)
!       vol(j, 5) = dg(j, 1)*dg(j, 2)*dd(j, 3)
!       vol(j, 6) = dd(j, 1)*dg(j, 2)*dd(j, 3)
!       vol(j, 7) = dg(j, 1)*dd(j, 2)*dd(j, 3)
!       vol(j, 8) = dd(j, 1)*dd(j, 2)*dd(j, 3)
!     end do
! #endif

!     ! Gather potential energ

!     if (poisson) then
!       do ind = 1, twotondim
!         do j = 1, np
!           family_part = typep(ind_part(j))%family
!           epot_part = mp(ind_part(j))*phi(indp(j, ind))*vol(j, ind)
!           epot_loc = epot_loc + epot_part
!           epot_families_loc(family_part) = epot_families_loc(family_part) + epot_part
!         end do
!       end do
!     end if

!   end subroutine epot_part_helper
  !#########################################################################
  !#########################################################################
  !#########################################################################
  !#########################################################################

end module deltaE_module
