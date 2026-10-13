subroutine add_viscosity_source_terms(ilevel)
    use amr_commons
    use hydro_commons
    use poisson_commons
    use pm_commons

    implicit none
    integer::ilevel
    !--------------------------------------------------------------------------
    ! This routine adds to unew the viscosity terms
    ! Only the momentum and the
    ! total energy are modified in array unew.
    !--------------------------------------------------------------------------
    integer::i, ind, iskip, nx_loc, j, iskip_son
    integer::ncache, igrid, ngrid, idim, jdim
    integer, dimension(1:8)::i3l, j3l, k3l, i3r, j3r, k3r, i3c, j3c, k3c
    real(kind=8)::scale, dx, dx_loc
    real(dp)::scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2

    integer, dimension(1:nvector), save::ind_grid, ind_cell
    integer, dimension(1:nvector, 0:twondim), save::igridn
    integer, dimension(1:nvector, 1:ndim), save::ind_left, ind_right
    real(dp), dimension(1:nvector, -1:1, -1:1, -1:1, 1:ndim), save::vel_neigh
    real(dp) :: v_ileft_jtop, v_iright_jtop, v_ileft_jbottom, v_iright_jbottom
    integer, dimension(1:nvector, 1:threetondim), save::nbors_father_cells
    real(dp), dimension(1:nvector, 1:ndim), save::den_left, den_right
    real(dp), dimension(1:nvector, 1:ndim), save::viscosity_term
    real(dp), dimension(1:nvector) :: energy_visc, e_prev
    real(dp) :: nu_viscosity
    real(dp), dimension(1:ndim) :: vel
    real(dp) :: den
    real(dp) ::vel_left, vel_right
    real(dp) :: den_dvel_left, den_dvel_right ! density times velocity derivative left and right
     real(dp):: den_En_left, den_En_right ! density times velocity times velocity derivative left and right
    real(dp) :: d = 0, u = 0, v = 0, w = 0, du = 0, dv = 0, dw = 0, e_kin = 0, e_nokin = 0
#ifdef SOLVERMHD
    real(dp) :: A, B, C
#endif

    real(dp), dimension(1:twotondim, 1:3)::xc
    real(dp), dimension(1:3)::skip_loc

    real(dp), dimension(1:ndim):: nu_viscosity_left, nu_viscosity_right

    integer, dimension(1:nvector, 0:twondim), save::ibuffer_father
    real(dp), dimension(1:nvector, 0:twondim, 1:nvar), save::u1
    real(dp), dimension(1:nvector, 1:twotondim, 1:nvar), save::u2
    real(dp), dimension(1:nvector, iu1:iu2, ju1:ju2, ku1:ku2, 1:nvar), save::uloc
    real(dp), dimension(1:nvector, iu1:iu2, ju1:ju2, ku1:ku2, 1:ndim), save::x_diag

    integer::i1, j1, k1, i2, j2, k2, i3, j3, k3, nexist, nbuffer, ind_son, ind_father, ivar
    integer, dimension(1:nvector), save::igrid_nbor, ind_buffer, ind_exist, ind_nexist
    integer::i1min, i1max, j1min, j1max, k1min, k1max
    integer::i2min, i2max, j2min, j2max, k2min, k2max
    integer::i3min, i3max, j3min, j3max, k3min, k3max
    integer :: e(3,3)

    if (numbtot(1, ilevel) == 0) return
    if (verbose) write (*, 111) ilevel

    nx_loc = icoarse_max - icoarse_min + 1
    skip_loc = (/0.0d0, 0.0d0, 0.0d0/)
    if (ndim > 0) skip_loc(1) = dble(icoarse_min)
    if (ndim > 1) skip_loc(2) = dble(jcoarse_min)
    if (ndim > 2) skip_loc(3) = dble(kcoarse_min)
    scale = boxlen/dble(nx_loc)
    dx = 0.5D0**ilevel
    dx_loc = dx*scale

    call units(scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2)

  ! These arrays are neccesary to obtain the correct neighbors according to each cell
    i3l(1:8) = (/0, 1, 0, 1, 0, 1, 0, 1/); j3l(1:8) = (/0, 0, 1, 1, 0, 0, 1, 1/); k3l(1:8) = (/0, 0, 0, 0, 1, 1, 1, 1/)
    i3c(1:8) = (/1, 2, 1, 2, 1, 2, 1, 2/); j3c(1:8) = (/1, 1, 2, 2, 1, 1, 2, 2/); k3c(1:8) = (/1, 1, 1, 1, 2, 2, 2, 2/)
    i3r(1:8) = (/2, 3, 2, 3, 2, 3, 2, 3/); j3r(1:8) = (/2, 2, 3, 3, 2, 2, 3, 3/); k3r(1:8) = (/2, 2, 2, 2, 3, 3, 3, 3/)

    ! Unit offset, useful to access diagonal neighbors
    e = reshape([ &
        1,0,0, &
        0,1,0, &
        0,0,1 ], shape(e))

    ! Loop over myid grids by vector sweeps
    ncache = active(ilevel)%ngrid
    do igrid = 1, ncache, nvector

      ! Gather nvector grids
      ngrid = MIN(nvector, ncache - igrid + 1)
      do i = 1, ngrid
        ind_grid(i) = active(ilevel)%igrid(igrid + i - 1)
      end do

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Here begins the diagonal neighbors block of code !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !This was taken from the godfine1 subroutine in this same file

  ! TODO : just call the stencil method

      ! Compute father cell index
      do i = 1, ngrid
        ind_cell(i) = father(ind_grid(i))
      end do

      call get3cubefather(ind_cell, nbors_father_cells, ngrid, ilevel)

      ! Integer constants
      i1min = 0; i1max = 0; i2min = 0; i2max = 0; i3min = 1; i3max = 1
      j1min = 0; j1max = 0; j2min = 0; j2max = 0; j3min = 1; j3max = 1
      k1min = 0; k1max = 0; k2min = 0; k2max = 0; k3min = 1; k3max = 1
      if (ndim > 0) then
        i1max = 2; i2max = 1; i3max = 2
      end if
      if (ndim > 1) then
        j1max = 2; j2max = 1; j3max = 2
      end if
      if (ndim > 2) then
        k1max = 2; k2max = 1; k3max = 2
      end if

      !---------------------------
      ! Gather 6x6x6 cells stencil
      !---------------------------
      ! Loop over 3x3x3 neighboring father cells
      do k1 = k1min, k1max
      do j1 = j1min, j1max
      do i1 = i1min, i1max

        ! Check if neighboring grid exists
        nbuffer = 0
        nexist = 0
        ind_father = 1 + i1 + 3*j1 + 9*k1
        do i = 1, ngrid
          igrid_nbor(i) = son(nbors_father_cells(i, ind_father))
          if (igrid_nbor(i) > 0) then
            nexist = nexist + 1
            ind_exist(nexist) = i
          else
            nbuffer = nbuffer + 1
            ind_nexist(nbuffer) = i
            ind_buffer(nbuffer) = nbors_father_cells(i, ind_father)
          end if
        end do

        ! If not, interpolate hydro variables from parent cells
        if (nbuffer > 0) then
          call getnborfather(ind_buffer, ibuffer_father, nbuffer, ilevel)
          do j = 0, twondim
            do ivar = 1, nvar
              do i = 1, nbuffer
                u1(i, j, ivar) = uold(ibuffer_father(i, j), ivar)
              end do
            end do
          end do
          call interpol_hydro(u1, u2, nbuffer)
        end if

        ! Loop over 2x2x2 cells
        do k2 = k2min, k2max
        do j2 = j2min, j2max
        do i2 = i2min, i2max

          ind_son = 1 + i2 + 2*j2 + 4*k2
          iskip = ncoarse + (ind_son - 1)*ngridmax
          do i = 1, nexist
            ind_cell(i) = iskip + igrid_nbor(ind_exist(i))
          end do

          i3 = 1; j3 = 1; k3 = 1
          if (ndim > 0) then
            i3 = 1 + 2*(i1 - 1) + i2
          end if
          if (ndim > 1) then
            j3 = 1 + 2*(j1 - 1) + j2
          end if
          if (ndim > 2) then
            k3 = 1 + 2*(k1 - 1) + k2
          end if

          ! Gather hydro variables
          do ivar = 1, nvar
            do i = 1, nexist
              uloc(ind_exist(i), i3, j3, k3, ivar) = uold(ind_cell(i), ivar)
            end do
            do i = 1, nbuffer
              uloc(ind_nexist(i), i3, j3, k3, ivar) = u2(i, ind_son, ivar)
            end do
          end do
          ! Gather neighboring positions for the dx and neighboring viscosities for the laplacian
          do idim = 1, ndim
            do i = 1, ngrid
              x_diag(i, i3, j3, k3, idim) = (xg(igrid_nbor(i), idim) + xc(ind_son, idim) - skip_loc(idim))*scale
            end do
          end do

        end do
        end do
        end do
        ! End loop over cells

      end do
      end do
      end do
      ! End loop over neighboring grids
      ! Loop over grid cells
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!End of diagonal indices block !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Loop over cells
      do ind = 1, twotondim

        iskip_son = ncoarse + (ind - 1)*ngridmax

        ! Compute central cell index
        iskip = ncoarse + (ind - 1)*ngridmax
        do i = 1, ngrid
          ind_cell(i) = iskip + ind_grid(i)
        end do

        do i = 1, ngrid


          ! Central value
          vel_neigh(i, 0, 0, 0, 1:ndim) = uloc(i, i3c(ind), j3c(ind), k3c(ind), 2:ndim + 1)/max(uloc(i, i3c(ind), j3c(ind), k3c(ind), 1), smallr)

          ! left and right neighbors
          ! x coordinate
          vel_neigh(i, -1, 0, 0, 1:ndim) = uloc(i, i3l(ind), j3c(ind), k3c(ind), 2:ndim + 1)/max(uloc(i, i3l(ind), j3c(ind), k3c(ind), 1), smallr)
          vel_neigh(i, 1, 0, 0, 1:ndim)  = uloc(i, i3r(ind), j3c(ind), k3c(ind), 2:ndim + 1)/max(uloc(i, i3r(ind), j3c(ind), k3c(ind), 1), smallr)
          den_left(i, 1)  = max(uloc(i, i3l(ind), j3c(ind), k3c(ind), 1), smallr)
          den_right(i, 1) = max(uloc(i, i3r(ind), j3c(ind), k3c(ind), 1), smallr)

#if NDIM > 1
          ! y coordinate
          vel_neigh(i, 0, -1, 0, 1:ndim) = uloc(i, i3c(ind), j3l(ind), k3c(ind), 2:ndim + 1)/max(uloc(i, i3c(ind), j3l(ind), k3c(ind), 1), smallr)
          vel_neigh(i, 0, 1, 0, 1:ndim)  = uloc(i, i3c(ind), j3r(ind), k3c(ind), 2:ndim + 1)/max(uloc(i, i3c(ind), j3r(ind), k3c(ind), 1), smallr)
          den_left(i, 2)  = max(uloc(i, i3c(ind), j3l(ind), k3c(ind), 1), smallr)
          den_right(i, 2) = max(uloc(i, i3c(ind), j3r(ind), k3c(ind), 1), smallr)

#endif
#if NDIM > 2
          ! z coordinate
          vel_neigh(i, 0, 0, -1, 1:ndim) = uloc(i, i3c(ind), j3c(ind), k3l(ind), 2:ndim + 1)/max(uloc(i, i3c(ind), j3c(ind), k3l(ind), 1), smallr)
          vel_neigh(i, 0, 0, 1, 1:ndim)  = uloc(i, i3c(ind), j3c(ind), k3r(ind), 2:ndim + 1)/max(uloc(i, i3c(ind), j3c(ind), k3r(ind), 1), smallr)
          den_left(i, 3)  = max(uloc(i, i3c(ind), j3c(ind), k3l(ind), 1), smallr)
          den_right(i, 3) = max(uloc(i, i3c(ind), j3c(ind), k3r(ind), 1), smallr)
#endif

          ! Diagonal neighbors
          ! same z coordinate
          vel_neigh(i, -1, -1, 0, 1:ndim) = uloc(i, i3l(ind), j3l(ind), k3c(ind), 2:ndim + 1)/max(uloc(i, i3l(ind), j3l(ind), k3c(ind), 1), smallr)
          vel_neigh(i, -1, 1, 0, 1:ndim)  = uloc(i, i3l(ind), j3r(ind), k3c(ind), 2:ndim + 1)/max(uloc(i, i3l(ind), j3r(ind), k3c(ind), 1), smallr)
          vel_neigh(i, 1, -1, 0, 1:ndim)  = uloc(i, i3r(ind), j3l(ind), k3c(ind), 2:ndim + 1)/max(uloc(i, i3r(ind), j3l(ind), k3c(ind), 1), smallr)
          vel_neigh(i, 1, 1,  0, 1:ndim)  = uloc(i, i3r(ind), j3r(ind), k3c(ind), 2:ndim + 1)/max(uloc(i, i3r(ind), j3r(ind), k3c(ind), 1), smallr)

          ! same y coordinate
          vel_neigh(i, -1, 0, -1, 1:ndim) = uloc(i, i3l(ind), j3c(ind), k3l(ind), 2:ndim + 1)/max(uloc(i, i3l(ind), j3c(ind), k3l(ind), 1), smallr)
          vel_neigh(i, -1, 0, 1, 1:ndim)  = uloc(i, i3l(ind), j3c(ind), k3r(ind), 2:ndim + 1)/max(uloc(i, i3l(ind), j3c(ind), k3r(ind), 1), smallr)
          vel_neigh(i, 1, 0, -1, 1:ndim)  = uloc(i, i3r(ind), j3c(ind), k3l(ind), 2:ndim + 1)/max(uloc(i, i3r(ind), j3c(ind), k3l(ind), 1), smallr)
          vel_neigh(i, 1, 0, 1, 1:ndim)   = uloc(i, i3r(ind), j3c(ind), k3r(ind), 2:ndim + 1)/max(uloc(i, i3r(ind), j3c(ind), k3r(ind), 1), smallr)

          ! same x coordinate
          vel_neigh(i, 0, -1, -1, 1:ndim) = uloc(i, i3c(ind), j3l(ind), k3l(ind), 2:ndim + 1)/max(uloc(i, i3c(ind), j3l(ind), k3l(ind), 1), smallr)
          vel_neigh(i, 0, -1, 1, 1:ndim)  = uloc(i, i3c(ind), j3l(ind), k3r(ind), 2:ndim + 1)/max(uloc(i, i3c(ind), j3l(ind), k3r(ind), 1), smallr)
          vel_neigh(i, 0, 1, -1, 1:ndim)  = uloc(i, i3c(ind), j3r(ind), k3l(ind), 2:ndim + 1)/max(uloc(i, i3c(ind), j3r(ind), k3l(ind), 1), smallr)
          vel_neigh(i, 0, 1, 1, 1:ndim)   = uloc(i, i3c(ind), j3r(ind), k3r(ind), 2:ndim + 1)/max(uloc(i, i3c(ind), j3r(ind), k3r(ind), 1), smallr)

        end do


  ! Compute the viscosity term
        viscosity_term(1:ngrid, 1:ndim) = 0.0d0
        energy_visc(1:ngrid) = 0.0d0
        do i = 1, ngrid

          den = max(uold(ind_cell(i), 1), smallr) ! density of the cell
          vel(1:ndim) = uold(ind_cell(i), 2:ndim + 1)/ den! velocity of the cell

          ! Viscosity of the cell
          select case (viscosity_kind)
          case ('constant_uniform')

            nu_viscosity          = nu_viscosity_shear
            nu_viscosity_left(:)  = nu_viscosity_shear
            nu_viscosity_right(:) = nu_viscosity_shear

          end select

          ! Add non crossed terms (F1_j = d_i (rho nu d_i v_j ))
          ! + (1 - 2/dim) d_i (rho nu d_i v_i) from the non-crossed terms F2 and F3 when idim = jdim
          ! Same for energy : E1 = d_i (rho nu v_j d_i v_j )) 
          ! + (1 - 2/dim) d_i (rho nu v_i d_i v_i) from the non-crossed terms E2 and E3 when idim = jdim

          do jdim = 1, ndim ! component of the laplacian and the velocity
            do idim = 1, ndim ! direction for derivatives

              ! e(i, j)  is the krocknecker symbole delta_ij
              vel_left  = vel_neigh(i, - e(1, idim), - e(2, idim), - e(3, idim), jdim)
              vel_right = vel_neigh(i,   e(1, idim),   e(2, idim),   e(3, idim), jdim)

              ! first derivative on the cell faces
              den_dvel_left  =  ((den + den_left(i,idim)) / 2.0 ) * (nu_viscosity + nu_viscosity_left(idim)) / 2.0 * ((vel(jdim) - vel_left ) / (dx_loc))
              den_dvel_right =  ((den + den_right(i,idim)) / 2.0) * (nu_viscosity + nu_viscosity_right(idim)) / 2.0 * ((vel_right - vel(jdim) ) /  (dx_loc))

              ! second derivative at the cell center
              viscosity_term(i, jdim) = viscosity_term(i, jdim) + (1 + e(idim, jdim) * (1.0 - 2.0 / ndim)) * (den_dvel_right - den_dvel_left) / (dx_loc)

              ! First derivative for the energy term
              den_En_left  =  ((vel(jdim) + vel_left) / 2.0 ) * den_dvel_left
              den_En_right =  ((vel(jdim) + vel_right) / 2.0) * den_dvel_right

              ! Second derivative for energy term
              energy_visc(i) = energy_visc(i) + (1 + e(idim, jdim) * (1.0 - 2.0 / ndim)) * (den_En_right - den_En_left)  / dx_loc

              ! Bulk viscosity
              if (idim == jdim) then
                den_dvel_left  =  ((den + den_left(i,idim)) / 2.0 ) * nu_viscosity_bulk * ((vel(jdim) - vel_left ) / (dx_loc))
                den_dvel_right =  ((den + den_right(i,idim)) / 2.0) * nu_viscosity_bulk * ((vel_right - vel(jdim) ) /  (dx_loc))
                viscosity_term(i, jdim) = viscosity_term(i, jdim) + (den_dvel_right - den_dvel_left) / (dx_loc)

                den_En_left  =  ((vel(jdim) + vel_left) / 2.0 ) * den_dvel_left
                den_En_right =  ((vel(jdim) + vel_right) / 2.0) * den_dvel_right
                energy_visc(i) = energy_visc(i) +  (den_En_right - den_En_left)  / dx_loc
              end if
            end do
          end do

          ! Add crossed terms to the laplacian

          den = max(uold(ind_cell(i), 1), smallr) ! density of the cell

          ! Crossed term F2_j = d_i (rho nu d_j v_i )
          ! For the energy: E2 = d_i (rho nu v_j d_j v_i )
          do jdim = 1, ndim
            do idim = 1, ndim
                if (idim == jdim) cycle ! this case corresponds to a non crossed term, which is already included in the previous loop

                v_ileft_jtop     = vel_neigh(i, - e(1, idim) - e(1, jdim), - e(2, idim) - e(2, jdim), - e(3, idim) - e(3, jdim), idim)
                v_ileft_jbottom  = vel_neigh(i, - e(1, idim) + e(1, jdim), - e(2, idim) + e(2, jdim), - e(3, idim) + e(3, jdim), idim)

                v_iright_jtop    = vel_neigh(i,   e(1, idim) - e(1, jdim),   e(2, idim) - e(2, jdim),   e(3, idim) - e(3, jdim), idim)
                v_iright_jbottom = vel_neigh(i,   e(1, idim) + e(1, jdim),   e(2, idim) + e(2, jdim),   e(3, idim) + e(3, jdim), idim)

                ! First derivative of the velocity at the left and right neighbor center in the direction jdim
                den_dvel_left  = den_left(i, idim)*nu_viscosity_left(idim) * (v_ileft_jbottom - v_ileft_jtop)/(2 * dx_loc)
                den_dvel_right = den_right(i, idim)*nu_viscosity_right(idim) * (v_iright_jbottom - v_iright_jtop)/(2 * dx_loc)

                ! Second derivative at the cell center
                viscosity_term(i, jdim) = viscosity_term(i, jdim) + (den_dvel_right - den_dvel_left) / (2* dx_loc)

                ! First derivative for the energy term
                den_En_left  =  ((vel(jdim) + vel_left) / 2.0 ) * den_dvel_left
                den_En_right =  ((vel(jdim) + vel_right) / 2.0) * den_dvel_right

                ! Second derivative for energy term
                energy_visc(i) = energy_visc(i) + (den_En_right - den_En_left) / (2* dx_loc)

            end do
          end do


          ! Crossed term F3_i = - (2/ndim) d_i rho nu d_j v_j
          ! For the energy: E3 = - (2/ndim) d_i rho nu v_i d_j v_j

          do jdim = 1, ndim
            do idim = 1, ndim
                if (idim == jdim) cycle ! this case corresponds to a non crossed term

                v_ileft_jtop     = vel_neigh(i, - e(1, idim) - e(1, jdim), - e(2, idim) - e(2, jdim), - e(3, idim) - e(3, jdim), jdim)
                v_ileft_jbottom  = vel_neigh(i, - e(1, idim) + e(1, jdim), - e(2, idim) + e(2, jdim), - e(3, idim) + e(3, jdim), jdim)

                v_iright_jtop    = vel_neigh(i,   e(1, idim) - e(1, jdim),   e(2, idim) - e(2, jdim),   e(3, idim) - e(3, jdim), jdim)
                v_iright_jbottom = vel_neigh(i,   e(1, idim) + e(1, jdim),   e(2, idim) + e(2, jdim),   e(3, idim) + e(3, jdim), jdim)

                ! First derivative of the velocity at the left and right neighbor center in the direction jdim
                den_dvel_left = den_left(i, idim) * nu_viscosity_left(idim)*(v_ileft_jbottom - v_ileft_jtop)/(2 * dx_loc)
                den_dvel_right = den_right(i, idim) * nu_viscosity_right(idim)*(v_iright_jbottom - v_iright_jtop)/(2 * dx_loc)

                !  - (2/ndim) d_i rho nu d_j v_j
                viscosity_term(i, idim) = viscosity_term(i, idim) - (2.0/ndim)*(den_dvel_right - den_dvel_left) / (2 * dx_loc)

                ! First derivative for the energy term
                den_En_left  =  ((vel(idim) + vel_left) / 2.0 ) * den_dvel_left
                den_En_right =  ((vel(idim) + vel_right) / 2.0) * den_dvel_right

                ! Second derivative for energy term
                energy_visc(i) = energy_visc(i) - (2.0/ndim) * (den_En_right - den_En_left) / (2* dx_loc)

                ! Bulk viscosity contribution
                den_dvel_left = den_left(i, idim) * nu_viscosity_bulk*(v_ileft_jbottom - v_ileft_jtop)/(2 * dx_loc)
                den_dvel_right = den_right(i, idim) * nu_viscosity_bulk*(v_iright_jbottom - v_iright_jtop)/(2 * dx_loc)

                viscosity_term(i, idim) = viscosity_term(i, idim) + (den_dvel_right - den_dvel_left) / (2 * dx_loc)

                den_En_left  =  ((vel(idim) + vel_left) / 2.0 ) * den_dvel_left
                den_En_right =  ((vel(idim) + vel_right) / 2.0) * den_dvel_right
                energy_visc(i) = energy_visc(i) + (den_En_right - den_En_left) / (2* dx_loc)
            end do
          end do

        end do

        do i = 1, ngrid
          ! Add viscosity term at time t
          d = max(unew(ind_cell(i), 1), smallr)
          du = 0.0d0; dv = 0.0d0; dw = 0.0d0
          du = unew(ind_cell(i), 2)
          dv = unew(ind_cell(i), 3)
#if NDIM > 2
          dw = unew(ind_cell(i), 4)
#endif
          du = du + viscosity_term(i, 1)*dtnew(ilevel)
          unew(ind_cell(i), 2) = du
          dv = dv + viscosity_term(i, 2)*dtnew(ilevel)
          unew(ind_cell(i), 3) = dv
#if NDIM > 2
          dw = dw + viscosity_term(i, 3)*dtnew(ilevel)
          unew(ind_cell(i), 4) = dw
#endif

          ! Update the energy
          e_prev = unew(ind_cell(i), ndim + 2)
          unew(ind_cell(i), ndim + 2) = e_prev + energy_visc(i)*dtnew(ilevel) 

        end do

      end do
      ! End loop over cells
    end do
    ! End loop over sweeps

  111 format('    Entering add_viscosity_terms for level ', i2)

  end subroutine add_viscosity_source_terms
