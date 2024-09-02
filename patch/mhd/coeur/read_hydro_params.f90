subroutine read_hydro_params(nml_ok)
    use amr_commons
    use hydro_commons
    implicit none
#ifndef WITHOUTMPI
    include 'mpif.h'
#endif
    logical :: nml_ok
    ! --------------------------------------------------
    ! Local variables
    ! --------------------------------------------------
    integer :: i, idim, nboundary_true=0
    integer , dimension(1:MAXBOUND) :: bound_type
    real(dp) :: scale, ek_bound, em_bound

    ! --------------------------------------------------
    ! Namelist definitions
    ! --------------------------------------------------
    namelist / init_params / filetype, initfile, multiple, nregion, region_type &
        & , x_center, y_center, z_center, aexp_ini &
        & , length_x, length_y, length_z, exp_region &
        & , d_region, u_region, v_region, w_region, p_region &
        & , alpha_dense_core, beta_dense_core, crit_dense_core, delta_rho &
#if NENER > 0
        & , prad_region &
#endif
        & , A_region, B_region, C_region
    namelist / hydro_params / gamma, courant_factor, smallr, smallc &
        & , niter_riemann, slope_type, slope_mag_type &
#if NENER > 0
        & , gamma_rad &
#endif
        & , pressure_fix, beta_fix, scheme, riemann, riemann2d
    namelist / refine_params / x_refine, y_refine, z_refine, r_refine &
        & , a_refine, b_refine, exp_refine, jeans_refine, mass_cut_refine &
        & , m_refine, mass_sph, err_grad_d, err_grad_p, err_grad_u &
        & , err_grad_A, err_grad_B, err_grad_C, err_grad_B2 &
        & , floor_d, floor_u, floor_p, ivar_refine, var_cut_refine &
        & , floor_A, floor_B, floor_C, floor_B2 &
        & , interpol_var, interpol_type, interpol_mag_type, sink_refine
    namelist / boundary_params / nboundary, bound_type &
        & , ibound_min, ibound_max, jbound_min, jbound_max &
        & , kbound_min, kbound_max &
        & , d_bound, u_bound, v_bound, w_bound, p_bound &
        & , A_bound, B_bound, C_bound
    namelist / physics_params / cooling, haardt_madau, metal, isothermal &
        & , m_star, n_star, T2_star, g_star, del_star, eps_star, jeans_ncells &
        & , eta_sn, yield, rbubble, f_ek, ndebris, f_w, mass_gmc, kappa_IR &
        & , J21, a_spec, z_ave, z_reion, eta_mag, delayed_cooling &
        & , self_shielding, smbh, agn, B_ave, t_diss

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

    ! ------------------------------------------------
    ! set ischeme
    ! ------------------------------------------------
    SELECT CASE (scheme)
        CASE ('muscl')
        ischeme = 0
        CASE ('induction')
        ischeme = 1

        CASE DEFAULT
        write(*, *) 'unknown scheme'
        call clean_stop
    END SELECT
    ! ------------------------------------------------
    ! set iriemann
    ! ------------------------------------------------
    SELECT CASE (riemann)
        CASE ('llf')
        iriemann = 0
        CASE ('roe')
        iriemann = 1
        CASE ('hll')
        iriemann = 2
        CASE ('hlld')
        iriemann = 3
        CASE ('upwind')
        iriemann = 4
        CASE ('hydro')
        iriemann = 5

        CASE DEFAULT
        write(*, *) 'unknown riemann solver'
        call clean_stop
    END SELECT
    ! ------------------------------------------------
    ! set iriemann
    ! ------------------------------------------------
    SELECT CASE (riemann2d)
        CASE ('llf')
        iriemann2d = 0
        CASE ('roe')
        iriemann2d = 1
        CASE ('upwind')
        iriemann2d = 2
        CASE ('hll')
        iriemann2d = 3
        CASE ('hlla')
        iriemann2d = 4
        CASE ('hlld')
        iriemann2d = 5
        CASE DEFAULT
        write(*, *) 'unknown 2D riemann solver'
        call clean_stop
    END SELECT

    ! --------------------------------------------------
    ! Make sure virtual boundaries are expanded to
    ! account for staggered mesh representation
    ! --------------------------------------------------
    nexpand_bound = 2

    ! --------------------------------------------------
    ! Check for star formation
    ! --------------------------------------------------
    if (eps_star > 0) then
        ! For historical reference:
        ! t_star=0.1635449d0*(n_star/0.1d0)**(-0.5d0)/eps_star
        star = .true.
        pic = .true.
    end if

    ! --------------------------------------------------
    ! Check for metal
    ! --------------------------------------------------
    if (metal .and. nvar < (ndim + 6)) then
        if (myid == 1) write(*, *) 'Error: metals need nvar >= ndim+6'
        if (myid == 1) write(*, *) 'Modify hydro_parameters.f90 and recompile'
        nml_ok = .false.
    end if

    ! --------------------------------------------------
    ! Check for non-thermal energies
    ! --------------------------------------------------
#if NENER > 0
    if (nvar < (8 + nener)) then
        if (myid == 1) write(*, *) 'Error: non-thermal energy need nvar >= ndim+2+nener'
        if (myid == 1) write(*, *) 'Modify NENER and recompile'
        nml_ok = .false.
    end if
#endif

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
        boundary_var(i, 1) = MAX(d_bound(i), smallr)
        boundary_var(i, 2) = d_bound(i) * u_bound(i)
        boundary_var(i, 3) = d_bound(i) * v_bound(i)
        boundary_var(i, 4) = d_bound(i) * w_bound(i)
        boundary_var(i, 6) = A_bound(i)
        boundary_var(i, 7) = B_bound(i)
        boundary_var(i, 8) = C_bound(i)
        boundary_var(i, nvar + 1) = A_bound(i)
        boundary_var(i, nvar + 2) = B_bound(i)
        boundary_var(i, nvar + 3) = C_bound(i)
        ek_bound = 0.5d0 * d_bound(i) * (u_bound(i) ** 2 + v_bound(i) ** 2 + w_bound(i) ** 2)
        em_bound = 0.5d0 * (A_bound(i) ** 2 + B_bound(i) ** 2 + C_bound(i) ** 2)
        boundary_var(i, 5) = ek_bound + em_bound + P_bound(i) / (gamma - 1.0d0)
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

    ! -----------------------------------
    ! Sort out passive variable indices
    ! -----------------------------------
    imetal = 9 + nener
    idelay = imetal
    if (metal) idelay = imetal + 1
    ixion = idelay
    if (delayed_cooling) ixion = idelay + 1
    ichem = ixion
    if (aton) ichem = ixion + 1

    ! -----------------------------------
    ! Set magnetic slope limiters
    ! -----------------------------------
    if (slope_mag_type == - 1) then
        slope_mag_type = slope_type
    end if
    if (interpol_mag_type == - 1) then
        interpol_mag_type = interpol_type
    end if

end subroutine read_hydro_params
