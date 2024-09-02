!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmpdt(uu, gg, dx, dt, ncell)
    use amr_parameters
    use hydro_parameters
    use const
    implicit none
    integer :: ncell
    real(dp) :: dx, dt
    real(dp), dimension(1:nvector, 1:nvar) :: uu
    real(dp), dimension(1:nvector, 1:ndim) :: gg

    real(dp) :: dtcell, smallp
    integer :: k, idim
#if NENER > 0
    integer :: irad
#endif

    smallp = smallc ** 2 / gamma

    ! Convert to primitive variables
    do k = 1, ncell
        uu(k, 1) = max(uu(k, 1), smallr)
    end do
    ! Velocity
    do idim = 1, ndim
        do k = 1, ncell
            uu(k, idim + 1) = uu(k, idim + 1) / uu(k, 1)
        end do
    end do
    ! Internal energy
    do idim = 1, ndim
        do k = 1, ncell
            uu(k, ndim + 2) = uu(k, ndim + 2) - half * uu(k, 1) * uu(k, idim + 1) ** 2
        end do
    end do
#if NENER > 0
    do irad = 1, nener
        do k = 1, ncell
            uu(k, ndim + 2) = uu(k, ndim + 2) - uu(k, ndim + 2 + irad)
        end do
    end do
#endif

    ! Debug
    if (debug) then
        do k = 1, ncell
            if (uu(k, ndim + 2) <= 0 .or. uu(k, 1) <= smallr) then
                write(*, *) 'stop in cmpdt'
                write(*, *) 'dx   =', dx
                write(*, *) 'ncell=', ncell
                write(*, *) 'rho  =', uu(k, 1)
                write(*, *) 'P    =', uu(k, ndim + 2)
                write(*, *) 'vel  =', uu(k, 2:ndim + 1)
                stop
            end if
        end do
    end if

    ! Compute pressure
    do k = 1, ncell
        uu(k, ndim + 2) = max((gamma - one) * uu(k, ndim + 2), uu(k, 1) * smallp)
    end do
#if NENER > 0
    do irad = 1, nener
        do k = 1, ncell
            uu(k, ndim + 2 + irad) = (gamma_rad(irad) - one) * uu(k, ndim + 2 + irad)
        end do
    end do
#endif

    ! Compute sound speed
    do k = 1, ncell
        uu(k, ndim + 2) = gamma * uu(k, ndim + 2)
    end do
#if NENER > 0
    do irad = 1, nener
        do k = 1, ncell
            uu(k, ndim + 2) = uu(k, ndim + 2) + gamma_rad(irad) * uu(k, ndim + 2 + irad)
        end do
    end do
#endif
    do k = 1, ncell
        uu(k, ndim + 2) = sqrt(uu(k, ndim + 2) / uu(k, 1))
    end do

    ! Compute wave speed
    do k = 1, ncell
        uu(k, ndim + 2) = dble(ndim) * uu(k, ndim + 2)
    end do
    do idim = 1, ndim
        do k = 1, ncell
            uu(k, ndim + 2) = uu(k, ndim + 2) + abs(uu(k, idim + 1))
        end do
    end do

    ! Compute gravity strength ratio
    do k = 1, ncell
        uu(k, 1) = zero
    end do
    do idim = 1, ndim
        do k = 1, ncell
            uu(k, 1) = uu(k, 1) + abs(gg(k, idim))
        end do
    end do
    do k = 1, ncell
        uu(k, 1) = uu(k, 1) * dx / uu(k, ndim + 2) ** 2
        uu(k, 1) = MAX(uu(k, 1), 0.0001_dp)
    end do

    ! Compute maximum time step for each authorized cell
    dt = courant_factor * dx / smallc
    do k = 1, ncell
        dtcell = dx / uu(k, ndim + 2) * (sqrt(one + two * courant_factor * uu(k, 1)) - one) / uu(k, 1)
        dt = min(dt, dtcell)
    end do

end subroutine cmpdt
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine hydro_refine(ug, um, ud, ok, nn)
    use amr_parameters
    use hydro_parameters
    use const
#ifdef RT
    use rt_parameters
#endif
    implicit none
    ! dummy arguments
    integer nn
    real(dp) :: ug(1:nvector, 1:nvar)
    real(dp) :: um(1:nvector, 1:nvar)
    real(dp) :: ud(1:nvector, 1:nvar)
    logical :: ok(1:nvector)

    integer :: k, idim
    real(dp), dimension(1:nvector), save :: eking, ekinm, ekind
    real(dp) :: dg, dm, dd, pg, pm, pd, vg, vm, vd, cg, cm, cd, error
#if NENER > 0
    integer :: irad
#endif

    ! Convert to primitive variables
    do k = 1, nn
        ug(k, 1) = max(ug(k, 1), smallr)
        um(k, 1) = max(um(k, 1), smallr)
        ud(k, 1) = max(ud(k, 1), smallr)
    end do
    ! Velocity
    do idim = 1, ndim
        do k = 1, nn
            ug(k, idim + 1) = ug(k, idim + 1) / ug(k, 1)
            um(k, idim + 1) = um(k, idim + 1) / um(k, 1)
            ud(k, idim + 1) = ud(k, idim + 1) / ud(k, 1)
        end do
    end do
    ! Pressure
    do k = 1, nn
        eking(k) = zero
        ekinm(k) = zero
        ekind(k) = zero
    end do
    do idim = 1, ndim
        do k = 1, nn
            eking(k) = eking(k) + half * ug(k, 1) * ug(k, idim + 1) ** 2
            ekinm(k) = ekinm(k) + half * um(k, 1) * um(k, idim + 1) ** 2
            ekind(k) = ekind(k) + half * ud(k, 1) * ud(k, idim + 1) ** 2
        end do
    end do
#if NENER > 0
    do irad = 1, nener
        do k = 1, nn
            eking(k) = eking(k) + ug(k, ndim + 2 + irad)
            ekinm(k) = ekinm(k) + um(k, ndim + 2 + irad)
            ekind(k) = ekind(k) + ud(k, ndim + 2 + irad)
        end do
    end do
#endif
    do k = 1, nn
        ug(k, ndim + 2) = (gamma - one) * (ug(k, ndim + 2) - eking(k))
        um(k, ndim + 2) = (gamma - one) * (um(k, ndim + 2) - ekinm(k))
        ud(k, ndim + 2) = (gamma - one) * (ud(k, ndim + 2) - ekind(k))
    end do
    ! Passive scalars
#if NVAR > NDIM + 2
    do idim = ndim + 3, nvar
        do k = 1, nn
            ug(k, idim) = ug(k, idim) / ug(k, 1)
            um(k, idim) = um(k, idim) / um(k, 1)
            ud(k, idim) = ud(k, idim) / ud(k, 1)
        end do
    end do
#endif

    ! Compute errors
    if (err_grad_d >= 0.) then
        do k = 1, nn
            dg = ug(k, 1); dm = um(k, 1); dd = ud(k, 1)
            error = 2.0d0 * MAX( &
                & ABS((dd - dm) / (dd + dm + floor_d)) , &
                & ABS((dm - dg) / (dm + dg + floor_d)) )
            ok(k) = ok(k) .or. error > err_grad_d
        end do
    end if

    if (err_grad_p >= 0.) then
        do k = 1, nn
            pg = ug(k, ndim + 2); pm = um(k, ndim + 2); pd = ud(k, ndim + 2)
            error = 2.0d0 * MAX( &
                & ABS((pd - pm) / (pd + pm + floor_p)), &
                & ABS((pm - pg) / (pm + pg + floor_p)) )
            ok(k) = ok(k) .or. error > err_grad_p
        end do
    end if

    if (err_grad_u >= 0.) then
        do idim = 1, ndim
            do k = 1, nn
                vg = ug(k, idim + 1); vm = um(k, idim + 1); vd = ud(k, idim + 1)
                cg = sqrt(max(gamma * ug(k, ndim + 2) / ug(k, 1), floor_u ** 2))
                cm = sqrt(max(gamma * um(k, ndim + 2) / um(k, 1), floor_u ** 2))
                cd = sqrt(max(gamma * ud(k, ndim + 2) / ud(k, 1), floor_u ** 2))
                error = 2.0d0 * MAX( &
                    & ABS((vd - vm) / (cd + cm + ABS(vd) + ABS(vm) + floor_u)) , &
                    & ABS((vm - vg) / (cm + cg + ABS(vm) + ABS(vg) + floor_u)) )
                ok(k) = ok(k) .or. error > err_grad_u
            end do
        end do
    end if

#ifdef RT
    ! Ionization state (only Hydrogen)
    if (rt_err_grad_xHII >= 0.) then ! ---------------------------------------
        do k = 1, nn
            dg = min(1d0, max(0d0, ug(k, iIons)))
            dm = min(1d0, max(0d0, um(k, iIons)))
            dd = min(1d0, max(0d0, ud(k, iIons)))
            error = 2.0d0 * MAX( &
                & ABS((dd - dm) / (dd + dm + rt_floor_xHII)) , &
                & ABS((dm - dg) / (dm + dg + rt_floor_xHII)) )
            ok(k) = ok(k) .or. error > rt_err_grad_xHII
        end do
    end if

    ! Neutral state (only Hydrogen)
    if (rt_err_grad_xHI  >= 0.) then ! ---------------------------------------
        do k = 1, nn
            dg = min(1d0, max(0d0, 1d0 - ug(k, iIons)))
            dm = min(1d0, max(0d0, 1d0 - um(k, iIons)))
            dd = min(1d0, max(0d0, 1d0 - ud(k, iIons)))
            error = 2.0d0 * MAX( &
                & ABS((dd - dm) / (dd + dm + rt_floor_xHI)) , &
                & ABS((dm - dg) / (dm + dg + rt_floor_xHI)) )
            ok(k) = ok(k) .or. error > rt_err_grad_xHI
        end do
    end if
#endif

end subroutine hydro_refine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine riemann_approx(qleft, qright, fgdnv, ngrid)
    use amr_parameters
    use hydro_parameters
    use const
    implicit none

    ! dummy arguments
    integer :: ngrid
    real(dp), dimension(1:nvector, 1:nvar) :: qleft, qright
    real(dp), dimension(1:nvector, 1:nvar + 1) :: fgdnv

    ! local arrays
    real(dp), dimension(1:nvector, 1:nvar + 1), save :: qgdnv
    real(dp), dimension(1:nvector), save :: rl   , ul   , pl   , cl
    real(dp), dimension(1:nvector), save :: rr   , ur   , pr   , cr
    real(dp), dimension(1:nvector), save :: ro   , uo   , po   , co
    real(dp), dimension(1:nvector), save :: rstar, ustar, pstar, cstar
    real(dp), dimension(1:nvector), save :: wl   , wr   , wo
    real(dp), dimension(1:nvector), save :: sgnm , spin , spout, ushock
    real(dp), dimension(1:nvector), save :: frac , delp , pold
    integer , dimension(1:nvector), save :: ind  , ind2

    ! local variables
    real(dp) :: smallp, gamma6, ql, qr, usr, usl, wwl, wwr, smallpp, entho, etot
    integer :: i, j, n, iter, n_new

    ! constants
    smallp = smallc ** 2 / gamma
    smallpp = smallr * smallp
    gamma6 = (gamma + one) / (two * gamma)
    entho = one / (gamma - one)

    ! Pressure, density and velocity
    do i = 1, ngrid
        rl(i) = MAX(qleft (i, 1), smallr)
        ul(i) =    qleft (i, 2)
        pl(i) = MAX(qleft (i, 3), rl(i) * smallp)
        rr(i) = MAX(qright(i, 1), smallr)
        ur(i) =    qright(i, 2)
        pr(i) = MAX(qright(i, 3), rr(i) * smallp)
    end do

    ! Lagrangian sound speed
    do i = 1, ngrid
        cl(i) = gamma * pl(i) * rl(i)
        cr(i) = gamma * pr(i) * rr(i)
    end do

    ! First guess
    do i = 1, ngrid
        wl(i) = sqrt(cl(i)); wr(i) = sqrt(cr(i))
        pstar(i) = ((wr(i) * pl(i) + wl(i) * pr(i)) + wl(i) * wr(i) * (ul(i) - ur(i))) / (wl(i) + wr(i))
        pstar(i) = MAX(pstar(i), 0.0_dp)
        pold (i) = pstar(i)
    end do
    n = ngrid
    do i = 1, n
        ind(i) = i
    end do

    ! Newton-Raphson iterations to find pstar at the required accuracy
    ! for a two-shock Riemann problem
    do iter = 1, niter_riemann
        do i = 1, n
            wwl = sqrt(cl(ind(i)) * (one + gamma6 * (pold(i) - pl(ind(i))) / pl(ind(i))))
            wwr = sqrt(cr(ind(i)) * (one + gamma6 * (pold(i) - pr(ind(i))) / pr(ind(i))))
            ql = two * wwl ** 3 / (wwl ** 2 + cl(ind(i)))
            qr = two * wwr ** 3 / (wwr ** 2 + cr(ind(i)))
            usl = ul(ind(i)) - (pold(i) - pl(ind(i))) / wwl
            usr = ur(ind(i)) + (pold(i) - pr(ind(i))) / wwr
            delp(i) = MAX(qr * ql / (qr + ql) * (usl - usr),- pold(i))
        end do
        do i = 1, n
            pold(i) = pold(i) + delp(i)
        end do
        ! Convergence indicator
        do i = 1, n
            uo(i) = ABS(delp(i) / (pold(i) + smallpp))
        end do
        n_new = 0
        do i = 1, n
            if (uo(i) > 1d-06) then
                n_new = n_new + 1
                ind2(n_new) = ind (i)
                po  (n_new) = pold(i)
            end if
        end do
        j = n_new
        do i = 1, n
            if (uo(i) <= 1d-06) then
                n_new = n_new + 1
                ind2(n_new) = ind (i)
                po  (n_new) = pold(i)
            end if
        end do
        ind (1:n) = ind2(1:n)
        pold(1:n) = po  (1:n)
        n = j
    end do

    ! Star region pressure
    ! for a two-shock Riemann problem
    do i = 1, ngrid
        pstar(ind(i)) = pold(i)
    end do
    do i = 1, ngrid
        wl(i) = sqrt(cl(i) * (one + gamma6 * (pstar(i) - pl(i)) / pl(i)))
        wr(i) = sqrt(cr(i) * (one + gamma6 * (pstar(i) - pr(i)) / pr(i)))
    end do

    ! Star region velocity
    ! for a two shock Riemann problem
    do i = 1, ngrid
        ustar(i) = half * (ul(i) + (pl(i) - pstar(i)) / wl(i) + &
            &           ur(i) - (pr(i) - pstar(i)) / wr(i) )
    end do

    ! Left going or right going contact wave
    do i = 1, ngrid
        sgnm(i) = sign(one, ustar(i))
    end do

    ! Left or right unperturbed state
    do i = 1, ngrid
        if (sgnm(i) == one) then
            ro(i) = rl(i)
            uo(i) = ul(i)
            po(i) = pl(i)
            wo(i) = wl(i)
        else
            ro(i) = rr(i)
            uo(i) = ur(i)
            po(i) = pr(i)
            wo(i) = wr(i)
        end if
    end do
    do i = 1, ngrid
        co(i) = max(smallc, sqrt(abs(gamma * po(i) / ro(i))))
    end do

    ! Star region density
    do i = 1, ngrid
        if (pstar(i) >= po(i)) then
            ! Shock
            rstar(i) = ro(i) / (one + ro(i) * (po(i) - pstar(i)) / wo(i) ** 2)
        else
            ! Rarefaction
            rstar(i) = ro(i) * (pstar(i) / po(i)) ** (one / gamma)
        end if
    end do

    do i = 1, ngrid
        ! Prevent vacuum formation in star region
        rstar(i) = max(rstar(i), smallr)
        ! Star region sound speed
        cstar(i) = sqrt(abs(gamma * pstar(i) / rstar(i)))
        cstar(i) = max(cstar(i), smallc)
        ! Compute rarefaction head and tail speed
        spout(i) = co   (i) - sgnm(i) * uo   (i)
        spin (i) = cstar(i) - sgnm(i) * ustar(i)
        ! Compute shock speed
        ushock(i) = wo(i) / ro(i) - sgnm(i) * uo(i)
    end do

    do i = 1, ngrid
        if (pstar(i) >= po(i)) then
            spout(i) = ushock(i)
            spin (i) = spout (i)
        end if
    end do

    ! Sample the solution at x/t=0
    do i = 1, ngrid
        if (spout(i) <= zero) then
            qgdnv(i, 1) = ro(i)
            qgdnv(i, 2) = uo(i)
            qgdnv(i, 3) = po(i)
        else if (spin(i) >= zero) then
            qgdnv(i, 1) = rstar(i)
            qgdnv(i, 2) = ustar(i)
            qgdnv(i, 3) = pstar(i)
        else
            frac(i) = spout(i) / (spout(i) - spin(i))
            qgdnv(i, 2) = frac(i) * ustar(i) + (one - frac(i)) * uo(i)
            qgdnv(i, 3) = frac(i) * pstar(i) + (one - frac(i)) * po(i)
            qgdnv(i, 1) = ro(i) * (qgdnv(i, 3) / po(i)) ** (one / gamma)
        end if
    end do

    ! Passive scalars
#if NVAR > 3
    do n = 4, nvar
        do i = 1, ngrid
            if (sgnm(i) == one) then
                qgdnv(i, n) = qleft(i, n)
            else
                qgdnv(i, n) = qright(i, n)
            end if
        end do
    end do
#endif
    ! Specific internal energy
    do i = 1, ngrid
        qgdnv(i, nvar + 1) = po(i) / ro(i) * entho
    end do

    ! Compute fluxes
    do i = 1, ngrid
        fgdnv(i, 1) = qgdnv(i, 1) * qgdnv(i, 2)  ! Mass density
        fgdnv(i, 2) = qgdnv(i, 3) + qgdnv(i, 1) * qgdnv(i, 2) ** 2  ! Normal momentum
        etot = qgdnv(i, 3) * entho + half * qgdnv(i, 1) * qgdnv(i, 2) ** 2
#if NDIM > 1
        etot = etot + half * qgdnv(i, 1) * qgdnv(i, 4) ** 2
#endif
#if NDIM > 2
        etot = etot + half * qgdnv(i, 1) * qgdnv(i, 5) ** 2
#endif
        fgdnv(i, 3) = qgdnv(i, 2) * (etot + qgdnv(i, 3))     ! Total energy
    end do
    ! Other advected quantities
    do n = 4, nvar + 1
        do i = 1, ngrid
            fgdnv(i, n) = fgdnv(i, 1) * qgdnv(i, n)
        end do
    end do


end subroutine riemann_approx
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine riemann_acoustic(qleft, qright, fgdnv, ngrid)
    use amr_parameters
    use hydro_parameters
    use const
    implicit none

    ! dummy arguments
    integer :: ngrid
    real(dp), dimension(1:nvector, 1:nvar) :: qleft, qright
    real(dp), dimension(1:nvector, 1:nvar + 1) :: fgdnv

    ! local variables
    integer :: i, n
    real(dp) :: smallp, entho, etot

    ! local arrays
    real(dp), dimension(1:nvector, 1:nvar + 1), save :: qgdnv
    real(dp), dimension(1:nvector), save :: rl   , ul   , pl   , cl
    real(dp), dimension(1:nvector), save :: rr   , ur   , pr   , cr
    real(dp), dimension(1:nvector), save :: ro   , uo   , po   , co
    real(dp), dimension(1:nvector), save :: rstar, ustar, pstar, cstar
    real(dp), dimension(1:nvector), save :: wl   , wr   , wo
    real(dp), dimension(1:nvector), save :: sgnm , spin , spout, ushock
    real(dp), dimension(1:nvector), save :: frac

    ! constants
    smallp = smallc ** 2 / gamma
    entho = one / (gamma - one)

    ! Initial states pressure, density and velocity
    do i = 1, ngrid
        rl(i) = max(qleft (i, 1), smallr)
        ul(i) =    qleft (i, 2)
        pl(i) = max(qleft (i, 3), rl(i) * smallp)
        rr(i) = max(qright(i, 1), smallr)
        ur(i) =    qright(i, 2)
        pr(i) = max(qright(i, 3), rr(i) * smallp)
    end do

    ! Acoustic star state
    do i = 1, ngrid
        cl(i) = sqrt(gamma * pl(i) / rl(i))
        cr(i) = sqrt(gamma * pr(i) / rr(i))
        wl(i) = cl(i) * rl(i)
        wr(i) = cr(i) * rr(i)
        pstar(i) = ((wr(i) * pl(i) + wl(i) * pr(i)) + wl(i) * wr(i) * (ul(i) - ur(i))) &
            &   / (wl(i) + wr(i))
        ustar(i) = ((wr(i) * ur(i) + wl(i) * ul(i)) + (pl(i) - pr(i))) &
            &   / (wl(i) + wr(i))
        !!$ pstar(i) = MAX(pstar(i),zero)
    end do

    ! Left going or right going contact wave
    do i = 1, ngrid
        sgnm(i) = sign(one, ustar(i))
    end do

    ! Left or right unperturbed state
    do i = 1, ngrid
        if (sgnm(i) == one) then
            ro(i) = rl(i)
            uo(i) = ul(i)
            po(i) = pl(i)
            wo(i) = wl(i)
            co(i) = cl(i)
        else
            ro(i) = rr(i)
            uo(i) = ur(i)
            po(i) = pr(i)
            wo(i) = wr(i)
            co(i) = cr(i)
        end if
    end do

    ! Star region density and sound speed
    do i = 1, ngrid
        rstar(i) = ro(i) + (pstar(i) - po(i)) / co(i) ** 2
        rstar(i) = max(rstar(i), smallr)
        cstar(i) = sqrt(abs(gamma * pstar(i) / rstar(i)))
        cstar(i) = max(cstar(i), smallc)
    end do

    ! Head and tail speed of rarefaction
    do i = 1, ngrid
        spout(i) = co   (i) - sgnm(i) * uo   (i)
        spin (i) = cstar(i) - sgnm(i) * ustar(i)
    end do

    ! Shock speed
    do i = 1, ngrid
        ushock(i) = half * (spin(i) + spout(i))
        ushock(i) = max(ushock(i),- sgnm(i) * ustar(i))
    end do
    do i = 1, ngrid
        if (pstar(i) >= po(i)) then
            spout(i) = ushock(i)
            spin (i) = spout (i)
        end if
    end do

    ! Sample the solution at x/t=0
    do i = 1, ngrid
        if (spout(i) < zero) then      ! Initial state
            qgdnv(i, 1) = ro(i)
            qgdnv(i, 2) = uo(i)
            qgdnv(i, 3) = po(i)
        else if (spin(i) >= zero) then  ! Star region
            qgdnv(i, 1) = rstar(i)
            qgdnv(i, 2) = ustar(i)
            qgdnv(i, 3) = pstar(i)
        else                        ! Rarefaction
            frac(i) = spout(i) / (spout(i) - spin(i))
            qgdnv(i, 1) = frac(i) * rstar(i) + (one - frac(i)) * ro(i)
            qgdnv(i, 2) = frac(i) * ustar(i) + (one - frac(i)) * uo(i)
            qgdnv(i, 3) = frac(i) * pstar(i) + (one - frac(i)) * po(i)
        end if
    end do

    ! Passive scalars
#if NVAR > 3
    do n = 4, nvar
        do i = 1, ngrid
            if (sgnm(i) == one) then
                qgdnv(i, n) = qleft (i, n)
            else
                qgdnv(i, n) = qright(i, n)
            end if
        end do
    end do
#endif
    ! Specific internal energy
    do i = 1, ngrid
        qgdnv(i, nvar + 1) = po(i) / ro(i) * entho
    end do

    ! Compute fluxes
    do i = 1, ngrid
        fgdnv(i, 1) = qgdnv(i, 1) * qgdnv(i, 2)  ! Mass density
        fgdnv(i, 2) = qgdnv(i, 3) + qgdnv(i, 1) * qgdnv(i, 2) ** 2  ! Normal momentum
        etot = qgdnv(i, 3) * entho + half * qgdnv(i, 1) * qgdnv(i, 2) ** 2
#if NDIM > 1
        etot = etot             + half * qgdnv(i, 1) * qgdnv(i, 4) ** 2
#endif
#if NDIM > 2
        etot = etot             + half * qgdnv(i, 1) * qgdnv(i, 5) ** 2
#endif
        fgdnv(i, 3) = qgdnv(i, 2) * (etot + qgdnv(i, 3))     ! Total energy
    end do
    ! Other advected quantities
    do n = 4, nvar + 1
        do i = 1, ngrid
            fgdnv(i, n) = fgdnv(i, 1) * qgdnv(i, n)
        end do
    end do

end subroutine riemann_acoustic
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine riemann_llf(qleft, qright, fgdnv, ngrid)
    use amr_parameters
    use hydro_parameters
    use const
    implicit none

    ! dummy arguments
    integer :: ngrid
    real(dp), dimension(1:nvector, 1:nvar) :: qleft, qright
    real(dp), dimension(1:nvector, 1:nvar + 1) :: fgdnv

    ! local arrays
    real(dp), dimension(1:nvector, 1:nvar + 1), save :: fleft, fright
    real(dp), dimension(1:nvector, 1:nvar + 1), save :: uleft, uright
    real(dp), dimension(1:nvector), save :: cmax

    ! local variables
    integer :: i, n
    real(dp) :: smallp, entho
    real(dp) :: rl   , ul   , pl   , cl
    real(dp) :: rr   , ur   , pr   , cr

    ! Constants
    smallp = smallc ** 2 / gamma
    entho = one / (gamma - one)

    ! ===========================
    ! Compute maximum wave speed
    ! ===========================
    do i = 1, ngrid
        ! Left states
        rl = max(qleft (i, 1), smallr)
        ul =     qleft (i, 2)
        pl = max(qleft (i, 3), rl * smallp)
        cl = gamma * pl
#if NENER > 0
        do n = 1, nener
            cl = cl + gamma_rad(n) * qleft(i, ndim + 2 + n)
        end do
#endif
        cl = sqrt(cl / rl)
        ! Right states
        rr = max(qright(i, 1), smallr)
        ur =     qright(i, 2)
        pr = max(qright(i, 3), rr * smallp)
        cr = gamma * pr
#if NENER > 0
        do n = 1, nener
            cr = cr + gamma_rad(n) * qright(i, ndim + 2 + n)
        end do
#endif
        cr = sqrt(cr / rr)
        ! Local max. wave speed
        cmax(i) = max(abs(ul) + cl, abs(ur) + cr)
    end do

    ! ===============================
    ! Compute conservative variables
    ! ===============================
    do i = 1, ngrid
        ! Mass density
        uleft (i, 1) = qleft (i, 1)
        uright(i, 1) = qright(i, 1)
        ! Normal momentum
        uleft (i, 2) = qleft (i, 1) * qleft (i, 2)
        uright(i, 2) = qright(i, 1) * qright(i, 2)
        ! Total energy
        uleft (i, 3) = qleft (i, 3) * entho + half * qleft (i, 1) * qleft (i, 2) ** 2
        uright(i, 3) = qright(i, 3) * entho + half * qright(i, 1) * qright(i, 2) ** 2
#if NDIM > 1
        uleft (i, 3) = uleft (i, 3)       + half * qleft (i, 1) * qleft (i, 4) ** 2
        uright(i, 3) = uright(i, 3)       + half * qright(i, 1) * qright(i, 4) ** 2
#endif
#if NDIM > 2
        uleft (i, 3) = uleft (i, 3)       + half * qleft (i, 1) * qleft (i, 5) ** 2
        uright(i, 3) = uright(i, 3)       + half * qright(i, 1) * qright(i, 5) ** 2
#endif
#if NENER > 0
        do n = 1, nener
            uleft (i, 3) = uleft (i, 3) + qleft (i, ndim + 2 + n) / (gamma_rad(n) - one)
            uright(i, 3) = uright(i, 3) + qright(i, ndim + 2 + n) / (gamma_rad(n) - one)
        end do
#endif
    end do
    ! Transverse velocities
#if NDIM > 1
    do n = 4, ndim + 2
        do i = 1, ngrid
            uleft (i, n) = qleft (i, 1) * qleft (i, n)
            uright(i, n) = qright(i, 1) * qright(i, n)
        end do
    end do
#endif
    ! Non-thermal energies
#if NENER > 0
    do n = 1, nener
        do i = 1, ngrid
            uleft (i, ndim + 2 + n) = qleft (i, ndim + 2 + n) / (gamma_rad(n) - one)
            uright(i, ndim + 2 + n) = qright(i, ndim + 2 + n) / (gamma_rad(n) - one)
        end do
    end do
#endif
    ! Other passively advected quantities
#if NVAR > 2 + NDIM + NENER
    do n = 3 + ndim + nener, nvar
        do i = 1, ngrid
            uleft (i, n) = qleft (i, 1) * qleft (i, n)
            uright(i, n) = qright(i, 1) * qright(i, n)
        end do
    end do
#endif
    ! Thermal energy
    do i = 1, ngrid
        uleft (i, nvar + 1) = qleft (i, 3) * entho
        uright(i, nvar + 1) = qright(i, 3) * entho
    end do

    ! ==============================
    ! Compute left and right fluxes
    ! ==============================
    do i = 1, ngrid
        ! Mass density
        fleft (i, 1) = qleft (i, 2) * uleft (i, 1)
        fright(i, 1) = qright(i, 2) * uright(i, 1)
        ! Normal momentum
        fleft (i, 2) = qleft (i, 2) * uleft (i, 2) + qleft (i, 3)
        fright(i, 2) = qright(i, 2) * uright(i, 2) + qright(i, 3)
#if NENER > 0
        do n = 1, nener
            fleft (i, 2) = fleft (i, 2) + qleft (i, ndim + 2 + n)
            fright(i, 2) = fright(i, 2) + qright(i, ndim + 2 + n)
        end do
#endif
        ! Total energy
        fleft (i, 3) = qleft (i, 2) * (uleft (i, 3) + qleft (i, 3))
        fright(i, 3) = qright(i, 2) * (uright(i, 3) + qright(i, 3))
#if NENER > 0
        do n = 1, nener
            fleft (i, 3) = fleft (i, 3) + qleft (i, 2) * qleft (i, ndim + 2 + n)
            fright(i, 3) = fright(i, 3) + qright(i, 2) * qright(i, ndim + 2 + n)
        end do
#endif
    end do
    ! Other passively advected quantities
    do n = 4, nvar + 1
        do i = 1, ngrid
            fleft (i, n) = qleft (i, 2) * uleft (i, n)
            fright(i, n) = qright(i, 2) * uright(i, n)
        end do
    end do

    ! =============================
    ! Compute Lax-Friedrich fluxes
    ! =============================
    do n = 1, nvar + 1
        do i = 1, ngrid
            fgdnv(i, n) = half * (fleft(i, n) + fright(i, n) - cmax(i) * (uright(i, n) - uleft(i, n)))
        end do
    end do

end subroutine riemann_llf
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine riemann_hll(qleft, qright, fgdnv, ngrid)
    USE amr_parameters
    USE const
    USE hydro_parameters
    ! 1D HLL Riemann solver
    IMPLICIT NONE
    integer :: ngrid
    real(dp), dimension(1:nvector, 1:nvar) :: qleft, qright
    real(dp), dimension(1:nvector, 1:nvar + 1) :: fgdnv

    real(dp), dimension(1:nvector, 1:nvar + 1), save :: fleft, fright
    real(dp), dimension(1:nvector, 1:nvar + 1), save :: uleft, uright
    real(dp), dimension(1:nvector), save :: SL, SR
    integer :: i, n
    real(dp) :: smallp, entho
    real(dp) :: rl   , ul   , pl   , cl
    real(dp) :: rr   , ur   , pr   , cr

    ! Constants
    smallp = smallc ** 2 / gamma
    entho = one / (gamma - one)

    ! ===========================
    ! Compute maximum wave speed
    ! ===========================
    do i = 1, ngrid
        ! Left states
        rl = max(qleft (i, 1), smallr)
        ul =     qleft (i, 2)
        pl = max(qleft (i, 3), rl * smallp)
        cl = gamma * pl
#if NENER > 0
        do n = 1, nener
            cl = cl + gamma_rad(n) * qleft(i, ndim + 2 + n)
        end do
#endif
        cl = sqrt(cl / rl)
        ! Right states
        rr = max(qright(i, 1), smallr)
        ur =     qright(i, 2)
        pr = max(qright(i, 3), rr * smallp)
        cr = gamma * pr
#if NENER > 0
        do n = 1, nener
            cr = cr + gamma_rad(n) * qright(i, ndim + 2 + n)
        end do
#endif
        cr = sqrt(cr / rr)
        ! Left and right max. wave speed
        SL(i) = min(min(ul, ur) - max(cl, cr), zero)
        SR(i) = max(max(ul, ur) + max(cl, cr), zero)
    end do

    ! ===============================
    ! Compute conservative variables
    ! ===============================
    do i = 1, ngrid
        ! Mass density
        uleft (i, 1) = qleft (i, 1)
        uright(i, 1) = qright(i, 1)
        ! Normal momentum
        uleft (i, 2) = qleft (i, 1) * qleft (i, 2)
        uright(i, 2) = qright(i, 1) * qright(i, 2)
        ! Total energy
        uleft (i, 3) = qleft (i, 3) * entho + half * qleft (i, 1) * qleft (i, 2) ** 2
        uright(i, 3) = qright(i, 3) * entho + half * qright(i, 1) * qright(i, 2) ** 2
#if NDIM > 1
        uleft (i, 3) = uleft (i, 3)       + half * qleft (i, 1) * qleft (i, 4) ** 2
        uright(i, 3) = uright(i, 3)       + half * qright(i, 1) * qright(i, 4) ** 2
#endif
#if NDIM > 2
        uleft (i, 3) = uleft (i, 3)       + half * qleft (i, 1) * qleft (i, 5) ** 2
        uright(i, 3) = uright(i, 3)       + half * qright(i, 1) * qright(i, 5) ** 2
#endif
#if NENER > 0
        do n = 1, nener
            uleft (i, 3) = uleft (i, 3) + qleft (i, ndim + 2 + n) / (gamma_rad(n) - one)
            uright(i, 3) = uright(i, 3) + qright(i, ndim + 2 + n) / (gamma_rad(n) - one)
        end do
#endif
    end do
    ! Transverse velocities
#if NDIM > 1
    do n = 4, ndim + 2
        do i = 1, ngrid
            uleft (i, n) = qleft (i, 1) * qleft (i, n)
            uright(i, n) = qright(i, 1) * qright(i, n)
        end do
    end do
#endif
    ! Non-thermal energies
#if NENER > 0
    do n = 1, nener
        do i = 1, ngrid
            uleft (i, ndim + 2 + n) = qleft (i, ndim + 2 + n) / (gamma_rad(n) - one)
            uright(i, ndim + 2 + n) = qright(i, ndim + 2 + n) / (gamma_rad(n) - one)
        end do
    end do
#endif
    ! Other passively advected quantities
#if NVAR > 2 + NDIM + NENER
    do n = 3 + ndim + nener, nvar
        do i = 1, ngrid
            uleft (i, n) = qleft (i, 1) * qleft (i, n)
            uright(i, n) = qright(i, 1) * qright(i, n)
        end do
    end do
#endif
    ! Thermal energy
    do i = 1, ngrid
        uleft (i, nvar + 1) = qleft (i, 3) * entho
        uright(i, nvar + 1) = qright(i, 3) * entho
    end do

    ! ==============================
    ! Compute left and right fluxes
    ! ==============================
    do i = 1, ngrid
        ! Mass density
        fleft (i, 1) = uleft (i, 2)
        fright(i, 1) = uright(i, 2)
        ! Normal momentum
        fleft (i, 2) = qleft (i, 3) + uleft (i, 2) * qleft (i, 2)
        fright(i, 2) = qright(i, 3) + uright(i, 2) * qright(i, 2)
#if NENER > 0
        do n = 1, nener
            fleft (i, 2) = fleft (i, 2) + qleft (i, ndim + 2 + n)
            fright(i, 2) = fright(i, 2) + qright(i, ndim + 2 + n)
        end do
#endif
        ! Total energy
        fleft (i, 3) = qleft (i, 2) * (uleft (i, 3) + qleft (i, 3))
        fright(i, 3) = qright(i, 2) * (uright(i, 3) + qright(i, 3))
#if NENER > 0
        do n = 1, nener
            fleft (i, 3) = fleft (i, 3) + qleft (i, 2) * qleft (i, ndim + 2 + n)
            fright(i, 3) = fright(i, 3) + qright(i, 2) * qright(i, ndim + 2 + n)
        end do
#endif
    end do
    ! Other advected quantities
    do n = 4, nvar + 1
        do i = 1, ngrid
            fleft (i, n) = qleft (i, 2) * uleft (i, n)
            fright(i, n) = qright(i, 2) * uright(i, n)
        end do
    end do

    ! ===================
    ! Compute HLL fluxes
    ! ===================
    do n = 1, nvar + 1
        do i = 1, ngrid
            fgdnv(i, n) = (SR(i) * fleft(i, n) - SL(i) * fright(i, n) &
                & + SR(i) * SL(i) * (uright(i, n) - uleft(i, n))) / (SR(i) - SL(i))
        end do
    end do

end subroutine riemann_hll
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine riemann_hllc(qleft, qright, fgdnv, ngrid)
    use amr_parameters
    use hydro_parameters
    use const
    implicit none

    ! HLLC Riemann solver (Toro)
    integer :: ngrid
    real(dp), dimension(1:nvector, 1:nvar) :: qleft, qright
    real(dp), dimension(1:nvector, 1:nvar + 1) :: fgdnv

    real(dp) :: SL, SR
    real(dp) :: entho
    real(dp) :: rl, pl, ul, ecinl, etotl, el, ptotl
    real(dp) :: rr, pr, ur, ecinr, etotr, er, ptotr
    real(dp) :: cfastl, rcl, rstarl, estarl
    real(dp) :: cfastr, rcr, rstarr, estarr
    real(dp) :: etotstarl, etotstarr
    real(dp) :: ustar, ptotstar
    real(dp) :: ro, uo, ptoto, etoto, eo
    real(dp) :: smallp
    integer :: ivar, i
#if NENER > 0
    REAL(dp), dimension(1:nener) :: eradl, eradr, erado
    REAL(dp), dimension(1:nener) :: eradstarl, eradstarr
    integer :: irad
#endif

    ! constants
    smallp = smallc ** 2 / gamma
    entho = one / (gamma - one)

    do i = 1, ngrid

        ! Left variables
        rl = max(qleft (i, 1), smallr)
        Pl = max(qleft (i, 3), rl * smallp)
        ul =    qleft (i, 2)

        el = Pl * entho
        ecinl = half * rl * ul * ul
#if NDIM > 1
        ecinl = ecinl + half * rl * qleft(i, 4) ** 2
#endif
#if NDIM > 2
        ecinl = ecinl + half * rl * qleft(i, 5) ** 2
#endif
        etotl = el + ecinl
#if NENER > 0
        do irad = 1, nener
            eradl(irad) = qleft(i, 2 + ndim + irad) / (gamma_rad(irad) - one)
            etotl = etotl + eradl(irad)
        end do
#endif
        Ptotl = Pl
#if NENER > 0
        do irad = 1, nener
            Ptotl = Ptotl + qleft(i, 2 + ndim + irad)
        end do
#endif

        ! Right variables
        rr = max(qright(i, 1), smallr)
        Pr = max(qright(i, 3), rr * smallp)
        ur =    qright(i, 2)

        er = Pr * entho
        ecinr = half * rr * ur * ur
#if NDIM > 1
        ecinr = ecinr + half * rr * qright(i, 4) ** 2
#endif
#if NDIM > 2
        ecinr = ecinr + half * rr * qright(i, 5) ** 2
#endif
        etotr = er + ecinr
#if NENER > 0
        do irad = 1, nener
            eradr(irad) = qright(i, 2 + ndim + irad) / (gamma_rad(irad) - one)
            etotr = etotr + eradr(irad)
        end do
#endif
        Ptotr = Pr
#if NENER > 0
        do irad = 1, nener
            Ptotr = Ptotr + qright(i, 2 + ndim + irad)
        end do
#endif

        ! Find the largest eigenvalues in the normal direction to the interface
        cfastl = gamma * Pl
#if NENER > 0
        do irad = 1, nener
            cfastl = cfastl + gamma_rad(irad) * qleft(i, ndim + 2 + irad)
        end do
#endif
        cfastl = sqrt(max(cfastl / rl, smallc ** 2))

        cfastr = gamma * Pr
#if NENER > 0
        do irad = 1, nener
            cfastr = cfastr + gamma_rad(irad) * qright(i, ndim + 2 + irad)
        end do
#endif
        cfastr = sqrt(max(cfastr / rr, smallc ** 2))

        ! Compute HLL wave speed
        SL = min(ul, ur) - max(cfastl, cfastr)
        SR = max(ul, ur) + max(cfastl, cfastr)

        ! Compute lagrangian sound speed
        rcl = rl * (ul - SL)
        rcr = rr * (SR - ur)

        ! Compute acoustic star state
        ustar   = (rcr * ur   + rcl * ul   +  (Ptotl - Ptotr)) / (rcr + rcl)
        Ptotstar = (rcr * Ptotl + rcl * Ptotr + rcl * rcr * (ul - ur)) / (rcr + rcl)

        ! Left star region variables
        rstarl = rl * (SL - ul) / (SL - ustar)
        etotstarl = ((SL - ul) * etotl - Ptotl * ul + Ptotstar * ustar) / (SL - ustar)
        estarl = el * (SL - ul) / (SL - ustar)
#if NENER > 0
        do irad = 1, nener
            eradstarl(irad) = eradl(irad) * (SL - ul) / (SL - ustar)
        end do
#endif

        ! Right star region variables
        rstarr = rr * (SR - ur) / (SR - ustar)
        etotstarr = ((SR - ur) * etotr - Ptotr * ur + Ptotstar * ustar) / (SR - ustar)
        estarr = er * (SR - ur) / (SR - ustar)
#if NENER > 0
        do irad = 1, nener
            eradstarr(irad) = eradr(irad) * (SR - ur) / (SR - ustar)
        end do
#endif

        ! Sample the solution at x/t=0
        if (SL > 0d0) then
            ro = rl
            uo = ul
            Ptoto = Ptotl
            etoto = etotl
            eo = el
#if NENER > 0
            do irad = 1, nener
                erado(irad) = eradl(irad)
            end do
#endif
        else if (ustar > 0d0) then
            ro = rstarl
            uo = ustar
            Ptoto = Ptotstar
            etoto = etotstarl
            eo = estarl
#if NENER > 0
            do irad = 1, nener
                erado(irad) = eradstarl(irad)
            end do
#endif
        else if (SR > 0d0) then
            ro = rstarr
            uo = ustar
            Ptoto = Ptotstar
            etoto = etotstarr
            eo = estarr
#if NENER > 0
            do irad = 1, nener
                erado(irad) = eradstarr(irad)
            end do
#endif
        else
            ro = rr
            uo = ur
            Ptoto = Ptotr
            etoto = etotr
            eo = er
#if NENER > 0
            do irad = 1, nener
                erado(irad) = eradr(irad)
            end do
#endif
        end if

        ! =========================
        ! Compute the Godunov flux
        ! =========================
        fgdnv(i, 1) = ro * uo
        fgdnv(i, 2) = ro * uo * uo + Ptoto
        fgdnv(i, 3) = (etoto + Ptoto) * uo
        ! Transverse velocities
#if NDIM > 1
        do ivar = 4, ndim + 2
            if (ustar > 0) then
                fgdnv(i, ivar) = ro * uo * qleft (i, ivar)
            else
                fgdnv(i, ivar) = ro * uo * qright(i, ivar)
            end if
        end do
#endif
        ! Non-thermal energies
#if NENER > 0
        do irad = 1, nener
            fgdnv(i, ndim + 2 + irad) = uo * erado(irad)
        end do
#endif
        ! Other passively advected quantities
#if NVAR > 2 + NDIM + NENER
        do ivar = 3 + ndim + nener, nvar
            if (ustar > 0) then
                fgdnv(i, ivar) = ro * uo * qleft (i, ivar)
            else
                fgdnv(i, ivar) = ro * uo * qright(i, ivar)
            end if
        end do
#endif
        ! Thermal energy
        fgdnv(i, nvar + 1) = uo * eo

    end do

end subroutine riemann_hllc
!###########################################################
!###########################################################
!###########################################################
!###########################################################
