subroutine read_hydro_params(nml_ok)
    use amr_commons
    use hydro_commons
    use mpi_mod
    implicit none
    logical :: nml_ok
    ! --------------------------------------------------
    ! Local variables
    ! --------------------------------------------------
    integer :: i, idim, nboundary_true=0
    integer , dimension(1:MAXBOUND) :: bound_type
    real(dp) :: scale, ek_bound, lor, h

    ! --------------------------------------------------
    ! Namelist definitions
    ! --------------------------------------------------
    namelist / init_params / filetype, initfile, multiple, nregion, region_type &
        & , x_center, y_center, z_center, aexp_ini &
        & , length_x, length_y, length_z, exp_region &
        & , d_region, u_region, v_region, w_region, p_region
    namelist / hydro_params / gamma, eos, courant_factor, smallr, smallc &
        & , niter_riemann, slope_type &
        & , pressure_fix, beta_fix, scheme, riemann
    namelist / refine_params / x_refine, y_refine, z_refine, r_refine &
        & , a_refine, b_refine, exp_refine, jeans_refine &
        & , m_refine, mass_sph, err_grad_d, err_grad_p, err_grad_u, err_grad_lor &
        & , floor_d, floor_u, floor_p &
        & , interpol_var, interpol_type
    namelist / boundary_params / nboundary, bound_type &
        & , ibound_min, ibound_max, jbound_min, jbound_max &
        & , kbound_min, kbound_max &
        & , d_bound, u_bound, v_bound, w_bound, p_bound
    namelist / physics_params / cooling, haardt_madau, metal, isothermal &
        & , n_star, T2_star, g_star, del_star, eps_star &
        & , eta_sn, yield, rbubble, f_ek, ndebris, f_w &
        & , J21, a_spec, z_ave, z_reion, n_sink

    ! Read namelist file
    rewind(1)
    read(1, NML = init_params, END = 101)
    goto 102
101 write(*, *) ' You need to set up namelist &INIT_PARAMS in parameter file'
    call clean_stop
102 rewind(1)
    if (nlevelmax > levelmin) read(1, NML = refine_params)
    rewind(1)
    if (hydro) read(1, NML = hydro_params)
    rewind(1)
    read(1, NML = boundary_params, END = 103)
    simple_boundary = .true.
    goto 104
103 simple_boundary = .false.
104 if (nboundary > MAXBOUND) then
        write(*, *) 'Error: nboundary>MAXBOUND'
        call clean_stop
    end if
    rewind(1)
    read(1, NML = physics_params, END = 105)
105 continue

    !!! check for inconsitencies in namelist
    if (nlevelmax > levelmin) then
        if (eos == 'TM') then
            if (interpol_var == 1) then
                write(*, *) ,'TM does only works with interpol_var=0'
                stop
            end if
        else
            if (interpol_var /= 1) then
                write(*, *) ,'ID EOS does only works with interpol_var=1'
                stop
            end if
        end if
    end if

    ! -------------------------------------------------
    ! This section deals with hydro boundary conditions
    ! -------------------------------------------------
    if (simple_boundary .and. nboundary == 0) then
        simple_boundary = .false.
    end if

    if (simple_boundary) then

        ! Compute new coarse grid boundaries
        do i = 1, nboundary
            if (ibound_min(i) * ibound_max(i) == 1 .and. ndim > 0 .and. bound_type(i) > 0) then
                nx = nx + 1
                if (ibound_min(i) == - 1) then
                    icoarse_min = icoarse_min + 1
                    icoarse_max = icoarse_max + 1
                end if
                nboundary_true = nboundary_true + 1
            end if
        end do
        do i = 1, nboundary
            if (jbound_min(i) * jbound_max(i) == 1 .and. ndim > 1 .and. bound_type(i) > 0) then
                ny = ny + 1
                if (jbound_min(i) == - 1) then
                    jcoarse_min = jcoarse_min + 1
                    jcoarse_max = jcoarse_max + 1
                end if
                nboundary_true = nboundary_true + 1
            end if
        end do
        do i = 1, nboundary
            if (kbound_min(i) * kbound_max(i) == 1 .and. ndim > 2 .and. bound_type(i) > 0) then
                nz = nz + 1
                if (kbound_min(i) == - 1) then
                    kcoarse_min = kcoarse_min + 1
                    kcoarse_max = kcoarse_max + 1
                end if
                nboundary_true = nboundary_true + 1
            end if
        end do

        ! Compute boundary geometry
        do i = 1, nboundary
            if (ibound_min(i) * ibound_max(i) == 1 .and. ndim > 0 .and. bound_type(i) > 0) then
                if (ibound_min(i) == - 1) then
                    ibound_min(i) = icoarse_min + ibound_min(i)
                    ibound_max(i) = icoarse_min + ibound_max(i)
                    if (bound_type(i) == 1) boundary_type(i) = 1
                    if (bound_type(i) == 2) boundary_type(i) = 11
                    if (bound_type(i) == 3) boundary_type(i) = 21
                else
                    ibound_min(i) = icoarse_max + ibound_min(i)
                    ibound_max(i) = icoarse_max + ibound_max(i)
                    if (bound_type(i) == 1) boundary_type(i) = 2
                    if (bound_type(i) == 2) boundary_type(i) = 12
                    if (bound_type(i) == 3) boundary_type(i) = 22
                end if
                if (ndim > 1) jbound_min(i) = jcoarse_min + jbound_min(i)
                if (ndim > 1) jbound_max(i) = jcoarse_max + jbound_max(i)
                if (ndim > 2) kbound_min(i) = kcoarse_min + kbound_min(i)
                if (ndim > 2) kbound_max(i) = kcoarse_max + kbound_max(i)
            else if (jbound_min(i) * jbound_max(i) == 1 .and. ndim > 1 .and. bound_type(i) > 0) then
                ibound_min(i) = icoarse_min + ibound_min(i)
                ibound_max(i) = icoarse_max + ibound_max(i)
                if (jbound_min(i) == - 1) then
                    jbound_min(i) = jcoarse_min + jbound_min(i)
                    jbound_max(i) = jcoarse_min + jbound_max(i)
                    if (bound_type(i) == 1) boundary_type(i) = 3
                    if (bound_type(i) == 2) boundary_type(i) = 13
                    if (bound_type(i) == 3) boundary_type(i) = 23
                else
                    jbound_min(i) = jcoarse_max + jbound_min(i)
                    jbound_max(i) = jcoarse_max + jbound_max(i)
                    if (bound_type(i) == 1) boundary_type(i) = 4
                    if (bound_type(i) == 2) boundary_type(i) = 14
                    if (bound_type(i) == 3) boundary_type(i) = 24
                end if
                if (ndim > 2) kbound_min(i) = kcoarse_min + kbound_min(i)
                if (ndim > 2) kbound_max(i) = kcoarse_max + kbound_max(i)
            else if (kbound_min(i) * kbound_max(i) == 1 .and. ndim > 2 .and. bound_type(i) > 0) then
                ibound_min(i) = icoarse_min + ibound_min(i)
                ibound_max(i) = icoarse_max + ibound_max(i)
                jbound_min(i) = jcoarse_min + jbound_min(i)
                jbound_max(i) = jcoarse_max + jbound_max(i)
                if (kbound_min(i) == - 1) then
                    kbound_min(i) = kcoarse_min + kbound_min(i)
                    kbound_max(i) = kcoarse_min + kbound_max(i)
                    if (bound_type(i) == 1) boundary_type(i) = 5
                    if (bound_type(i) == 2) boundary_type(i) = 15
                    if (bound_type(i) == 3) boundary_type(i) = 25
                else
                    kbound_min(i) = kcoarse_max + kbound_min(i)
                    kbound_max(i) = kcoarse_max + kbound_max(i)
                    if (bound_type(i) == 1) boundary_type(i) = 6
                    if (bound_type(i) == 2) boundary_type(i) = 16
                    if (bound_type(i) == 3) boundary_type(i) = 26
                end if
            end if
        end do
        do i = 1, nboundary
            ! Check for errors
            if ( (ibound_min(i) < 0 .or. ibound_max(i) > (nx - 1)) .and. (ndim > 0) .and. bound_type(i) > 0 ) then
                if (myid == 1) write(*, *) 'Error in the namelist'
                if (myid == 1) write(*, *) 'Check boundary conditions along X direction', i
                nml_ok = .false.
            end if
            if ( (jbound_min(i) < 0 .or. jbound_max(i) > (ny - 1)) .and. (ndim > 1) .and. bound_type(i) > 0) then
                if (myid == 1) write(*, *) 'Error in the namelist'
                if (myid == 1) write(*, *) 'Check boundary conditions along Y direction', i
                nml_ok = .false.
            end if
            if ( (kbound_min(i) < 0 .or. kbound_max(i) > (nz - 1)) .and. (ndim > 2) .and. bound_type(i) > 0) then
                if (myid == 1) write(*, *) 'Error in the namelist'
                if (myid == 1) write(*, *) 'Check boundary conditions along Z direction', i
                nml_ok = .false.
            end if
        end do
    end if
    nboundary = nboundary_true
    if (simple_boundary .and. nboundary == 0) then
        simple_boundary = .false.
    end if

    ! --------------------------------------------------
    ! Compute boundary conservative variables
    ! --------------------------------------------------
    do i = 1, nboundary
        lor = (1d0 - u_bound(i) ** 2 - v_bound(i) ** 2 - w_bound(i) ** 2) ** (- 1./ 2.)
        boundary_var(i, 1) = MAX(lor * d_bound(i), smallr)
        h = 1d0 + p_bound(i) / d_bound(i) * gamma / (gamma - 1)
        if (eos == 'TM') h = 5d0 / 2d0 * p_bound(i) / d_bound(i) + sqrt(9d0 / 4d0 * (p_bound(i) / d_bound(i)) ** 2 + 1)
        boundary_var(i, 2) = lor ** 2 * d_bound(i) * u_bound(i) * h
        boundary_var(i, 3) = lor ** 2 * d_bound(i) * v_bound(i) * h
        boundary_var(i, 4) = lor ** 2 * d_bound(i) * w_bound(i) * h
        boundary_var(i, 5) = lor ** h * d_bound(i) * h - p_bound(i)
    end do

    ! -----------------------------------
    ! Rearrange level dependent arrays
    ! -----------------------------------
    do i = nlevelmax, levelmin, - 1
        jeans_refine(i) = jeans_refine(i - levelmin + 1)
    end do
    do i = 1, levelmin - 1
        jeans_refine(i) = - 1.0
    end do

end subroutine read_hydro_params
