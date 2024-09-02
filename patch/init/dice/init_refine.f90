!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_refine
    use amr_commons
    use pm_commons
    use dice_commons
    implicit none
    ! -------------------------------------------
    ! This routine builds the initial AMR grid
    ! -------------------------------------------
    integer :: ilevel

    if (myid == 1) write(*, *) 'Building initial AMR grid'
    init = .true.

    ! Base refinement
    do ilevel = 1, levelmin
        call flag
        call refine
    end do

    ! Further refinements if necessary
    do ilevel = levelmin + 1, nlevelmax
        if (initfile(levelmin) /= ' '.and. initfile(ilevel) == ' ') exit
        if (hydro) call init_flow
#ifdef RT
        if (rt) call rt_init_flow
#endif
        if (ivar_refine == 0) call init_refmap
        call flag
        call refine
        if (nremap > 0) call load_balance
        if (numbtot(1, ilevel) == 0) exit
    end do

    ! Final pass to initialize the flow
    init = .false.
    if (hydro) call init_flow
#ifdef RT
    if (rt) call rt_init_flow
#endif

end subroutine init_refine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_refine_2
    ! --------------------------------------------------------------
    ! This routine builds additional refinements to the
    ! the initial AMR grid for filetype ne 'grafic'
    ! DICE patch: It is ensured that all the particles are
    ! transfered down to level 1 before initialising the grid
    ! --------------------------------------------------------------
    use amr_commons
    use hydro_commons
#ifdef RT
    use rt_hydro_commons
#endif
    use pm_commons
    use poisson_commons
    use dice_commons
    implicit none
    integer :: ilevel, i, ivar
    real(dp) :: eps_star2

    if (filetype /= 'grafic') then

        if (myid == 1 .and. amr_struct) then
            write(*, *) "Initial conditions with AMR data structure"
            write(*, '(A50)') "__________________________________________________"
        end if
        do i = levelmin, nlevelmax + 1
            ! DICE------
            do ilevel = levelmin - 1, 1, - 1
                if (pic) call merge_tree_fine(ilevel)
            end do
            ! ----------
            call refine_coarse
            do ilevel = 1, nlevelmax
                call build_comm(ilevel)
                call make_virtual_fine_int(cpu_map(1), ilevel)
                call refine_fine(ilevel)
                ! DICE------
                if (pic) call make_tree_fine(ilevel)
                ! ----------
                if (hydro) call init_flow_fine(ilevel)
                ! DICE------
                if (pic) then
                    call kill_tree_fine(ilevel)
                    call virtual_tree_fine(ilevel)
                end if
                ! ----------
#ifdef RT
                if (rt) call rt_init_flow_fine(ilevel)
#endif
            end do

            ! DICE------
            do ilevel = nlevelmax - 1, levelmin, - 1
                if (pic) call merge_tree_fine(ilevel)
            end do
            ! ----------
            if (nremap > 0) call load_balance

            do ilevel = levelmin, nlevelmax
                if (pic) call make_tree_fine(ilevel)
                if (poisson) call rho_fine(ilevel, 2)
                if (hydro) call init_flow_fine(ilevel)
                if (pic) then
                    call kill_tree_fine(ilevel)
                    call virtual_tree_fine(ilevel)
                end if
            end do

            do ilevel = nlevelmax, levelmin, - 1
                if (pic) call merge_tree_fine(ilevel)
                if (hydro) then
                    call upload_fine(ilevel)
#ifdef SOLVERmhd
                    do ivar = 1, nvar + 3
#else
                        do ivar = 1, nvar
#endif
                            call make_virtual_fine_dp(uold(1, ivar), ilevel)
#ifdef SOLVERmhd
                        end do
#else
                    end do
#endif
                    if (simple_boundary) call make_boundary_hydro(ilevel)
                end if
#ifdef RT
                if (rt) then
                    call rt_upload_fine(ilevel)
                    do ivar = 1, nrtvar
                        call make_virtual_fine_dp(rtuold(1, ivar), ilevel)
                    end do
                    if (simple_boundary) call rt_make_boundary_hydro(ilevel)
                end if
#endif
            end do

            do ilevel = nlevelmax, 1, - 1
                call flag_fine(ilevel, 2)
            end do
            call flag_coarse

        end do
        ! DICE------
        do ilevel = levelmin - 1, 1, - 1
            if (pic) call merge_tree_fine(ilevel)
        end do
        call kill_gas_part(1)
        do ilevel = 1, nlevelmax
            if (pic) then
                call make_tree_fine(ilevel)
                call kill_tree_fine(ilevel)
                call virtual_tree_fine(ilevel)
            end if
        end do
        do ilevel = nlevelmax, levelmin, - 1
            call merge_tree_fine(ilevel)
        end do
        deallocate(up)
        if (sf_virial) then
            eps_star2 = eps_star
            eps_star = 0d0
            do ilevel = nlevelmax, levelmin, - 1
                call star_formation(ilevel)
            end do
            eps_star = eps_star2
        end if
        dice_init = .false.
        ! ----------
    end if ! if .not. 'grafic'


#ifdef RT
    if (neq_chem .and. rt_is_init_xion) then
        if (myid == 1) write(*, *) 'Initializing ionization states from T profile'
        do ilevel = nlevelmax, 1, - 1
            call rt_init_xion(ilevel)
            call upload_fine(ilevel)
#ifdef SOLVERmhd
            do ivar = 1, nvar + 3
#else
                do ivar = 1, nvar
#endif
                    call make_virtual_fine_dp(uold(1, ivar), ilevel)
#ifdef SOLVERmhd
                end do
#else
            end do
#endif
            if (simple_boundary) call make_boundary_hydro(ilevel)
        end do
    end if
#endif

end subroutine init_refine_2
!################################################################
!################################################################
!################################################################
!################################################################
subroutine kill_gas_part(ilevel)
    use pm_commons
    use amr_commons
    use mpi_mod
    implicit none
    integer :: ilevel
    ! --------------------------------------------------------
    ! This subroutine removes the gas particles
    ! initially present in the gadget1 DICE output
    ! --------------------------------------------------------
    integer :: igrid, jgrid, ipart, jpart, next_part
    integer :: ig, ip, npart1, npart2, icpu, info
    integer, dimension(1:nvector) :: ind_grid, ind_part, ind_grid_part
    logical, dimension(1:nvector) :: ok=.true.
    integer :: npart_all
    integer, dimension(1:ncpu) :: npart_cpu, npart_cpu_all

    npart_cpu = 0
    npart_all = 0

    if (numbtot(1, ilevel) == 0) return
    ! Gather gas particles.
    ! Loop over cpus
    do icpu = 1, ncpu
        igrid = headl(icpu, ilevel)
        ig = 0
        ip = 0
        ! Loop over grids
        do jgrid = 1, numbl(icpu, ilevel)
            npart1 = numbp(igrid)  ! Number of particles in the grid
            npart2 = 0
            ! Count gas particles
            if (npart1 > 0) then
                ipart = headp(igrid)
                ! Loop over particles
                do jpart = 1, npart1
                    ! Save next particle   <--- Very important !!!
                    next_part = nextp(ipart)
                    if (idp(ipart) == 1) then
                        npart2 = npart2 + 1
                    end if
                    ipart = next_part  ! Go to next particle
                end do
                npart_cpu(myid) = npart_cpu(myid) + npart2
            end if
            ! Gather gas particles
            if (npart2 > 0) then
                ig = ig + 1
                ind_grid(ig) = igrid
                ipart = headp(igrid)
                ! Loop over particles
                do jpart = 1, npart1
                    ! Save next particle   <--- Very important !!!
                    next_part = nextp(ipart)
                    ! Select only gas particles
                    if (idp(ipart) == 1) then
                        if (ig == 0) then
                            ig = 1
                            ind_grid(ig) = igrid
                        end if
                        ip = ip + 1
                        ind_part(ip) = ipart
                        ind_grid_part(ip) = ig
                    end if
                    if (ip == nvector) then
                        call remove_list(ind_part, ind_grid_part, ok, ip)
                        call add_free_cond(ind_part, ok, ip)
                        ip = 0
                        ig = 0
                    end if
                    ipart = next_part  ! Go to next particle
                end do
                ! End loop over particles
            end if

            igrid = next(igrid)   ! Go to next grid
        end do

        ! End loop over grids
        if (ip > 0) then
            call remove_list(ind_part, ind_grid_part, ok, ip)
            call add_free_cond(ind_part, ok, ip)
        end if
    end do

#ifndef WITHOUTMPI
    ! Give an array of number of gas on each cpu available to all cpus
    call MPI_ALLREDUCE(npart_cpu, npart_cpu_all, ncpu, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, info)
#endif
    npart_all = sum(npart_cpu_all(1:ncpu))
    if (npart_all > 0) then
        if (myid == 1) then
            write(*, '(A50)') "__________________________________________________"
            write(*, '(A,I15)') ' Gas particles deleted ->', npart_all
            write(*, '(A50)') "__________________________________________________"
        end if
    end if
    npart_cpu(myid) = npart
#ifndef WITHOUTMPI
    ! Give an array of number of gas on each cpu available to all cpus
    call MPI_ALLREDUCE(npart_cpu, npart_cpu_all, ncpu, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, info)
#endif
    npart_all = sum(npart_cpu_all(1:ncpu))
    if (npart_all > 0) then
        if (myid == 1) then
            write(*, '(A50)') "__________________________________________________"
            write(*, '(A,I15)') ' Remaining particles ->', npart_all
            write(*, '(A50)') "__________________________________________________"
        end if
    end if

    do icpu = 1, ncpu
        igrid = headl(icpu, ilevel)
        ig = 0
        ip = 0
        ! Loop over grids
        do jgrid = 1, numbl(icpu, ilevel)
            npart1 = numbp(igrid)  ! Number of particles in the grid
            npart2 = 0
            ! Count gas particles
            if (npart1 > 0) then
                ipart = headp(igrid)
                ! Loop over particles
                do jpart = 1, npart1
                    ! Save next particle   <--- Very important !!!
                    next_part = nextp(ipart)
                    if (idp(ipart) > 0) idp(ipart) = idp(ipart) - 1
                    ipart = next_part  ! Go to next particle
                end do
                npart_cpu(myid) = npart_cpu(myid) + npart2
            end if
        end do
    end do


111 format('   Entering kill_gas_part for level ', I2)
    ! ---------------------------------------------
end subroutine
