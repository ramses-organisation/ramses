MODULE cr_flux_module

  ! Two-moment cosmic-ray closure and intercell-flux module.
  ! This is a PURE flux/closure module: no source terms, no NENER.

  use amr_commons
  use hydro_commons
  use hydro_parameters         ! iu1..ku2, if1..kf2, nvar, smallr
  use cr_parameters            ! ncr_groups, ncrvar, iCRu, gamma_cr, cr_vmax, ...
  implicit none

  private   ! default

  ! CR-specific stencil-loop bounds used by cmp_cr_wavespeeds.
  integer,parameter::ifcr1=0
  integer,parameter::jfcr1=1-ndim/2
  integer,parameter::kfcr1=1-ndim/3

  public cmp_cr_faces, rotatevec, invrotatevec

CONTAINS

!************************************************************************
subroutine cmp_cr_flux_tensors(uin_gas, uin_cr, iGrp, nGrid, ftens, vmax, bfield)

! Compute central fluxes for a CR group, for each cell in a vector
! of grids.
! The flux tensor is a three by four tensor (2*3 and 1*2 in 1D and 2D,
! respectively) where the first row is CR flux (x,y,z) and
! the other three rows compose the Eddington tensor (see Yiang&Peng 2017)
! input/output:
! uin_gas   => gas (rho,vel,B) variables of all cells in a vector of grids
! uin_cr    => CR (E,F) variables of all cells in a vector of grids
! igrp      => CR group number
! ngrid     => Number of 'valid' grids in uin.
! ftens     <=  Group flux tensors for all the cells.
!------------------------------------------------------------------------
  real(dp),dimension(nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3)::uin_gas
  real(dp),dimension(nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ncrvar)::uin_cr
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nDim+1,1:ndim)::ftens
  integer::iGrp, nGrid!---------------------------------------------------
  real(dp),dimension(1:ndim)::crflux
  real(dp)::Ecr, vmax
  integer::i, j, k, idim, jdim, n, icrE, nedge
  real(dp)::mu1,mu2,chi,b_norm2,crflux_norm
  real(dp),dimension(nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim)::bfield
  real(dp)::vmax2,Ecr2,aniso_term
  !------------------------------------------------------------------------

  icrE = iCRu+(ndim+1)*(iGrp-1) ! starting index of cr variables (base 1)
  ! Loop 6X6X6 cells in grid, from -1 to 4.
  do k = ku1, ku2
  do j = ju1, ju2
  do i = iu1, iu2

     nedge = 0              ! Check if we're at a corner and if so, cycle
     if(i.lt.1 .or. i.gt.2) nedge=nedge+1
     if(ndim.gt.1 .and. (j.lt.1 .or. j.gt.2)) nedge=nedge+1
     if(ndim.gt.2 .and. (k.lt.1 .or. k.gt.2)) nedge=nedge+1
     if(nedge.ge.2) cycle

     do n=1,ngrid
        Ecr =   uin_cr(n, i, j, k, icrE)            ! CR density in cell
        crflux = uin_cr(n,i,j,k,icrE+1 : icrE+ndim) ! CR flux vector
        if(Ecr .lt. 0d0) then
          write(*,*)'negative CR density in cmp_flux_tensors. -EXITING-'
          call clean_stop
        endif
        ftens(n,i,j,k,1,1:ndim)= crflux  !   First row is CR flux
        ! Rest is Eddington tensor
        if(isotropic_pressure) then
           ftens(n,i,j,k,2:ndim+1,1:ndim) = 0d0
           do idim = 1, ndim
              ftens(n,i,j,k,idim+1,idim) = &
                   Ecr*vmax**2/3d0
           enddo
        else
           ! M1 closure
           vmax2=vmax**2
           Ecr2=Ecr**2
           b_norm2    =SUM(bfield(n,i,j,k,1:ndim)**2)
           crflux_norm=SUM(crflux**2)
           mu1=MIN(crflux_norm/(vmax2*Ecr2),1.0d0)
           mu2=(3d0+4d0*mu1)/(5d0+2d0*sqrt(4d0-3d0*mu1))
           chi=0.5d0*(1d0-mu2)
           aniso_term=(1d0-3.0d0*chi)/b_norm2
           do idim = 1, ndim
              do jdim = 1, ndim
                 ftens(n,i,j,k,idim+1,jdim) = bfield(n,i,j,k,idim)*bfield(n,i,j,k,jdim)
              end do
           end do
           ftens(n,i,j,k,2:ndim+1,1:ndim) = ftens(n,i,j,k,2:ndim+1,1:ndim)*aniso_term
           do idim = 1, ndim
              ftens(n,i,j,k,idim+1,idim) = ftens(n,i,j,k,idim+1,idim) + chi
           end do
           ftens(n,i,j,k,2:ndim+1,1:ndim) = ftens(n,i,j,k,2:ndim+1,1:ndim)*Ecr*vmax2
        endif
    enddo
  enddo
  enddo
  enddo

end subroutine cmp_cr_flux_tensors

!************************************************************************
SUBROUTINE cmp_cr_wavespeeds(uin_gas, uin_cr, iGrp, ngrid, lmax, ilevel, dt)

  !  Compute CR wavespeeds for given vector of sub-grids.
  !
  !  inputs/outputs
  !  uin_gas     => input gas (rho,vel,B) cell states
  !  uin_cr      => input CR (E,F) cell states
  !  iGrp        => CR group number
  !  ngrid       => number of sub-grids of 3^ndim cells
  !  lmax       <=  return maximum cell wavespeeds
  !  ilevel     <=  current refinement level
  !  dt         <=  current CR timestep length
  !
  !  other vars
  !  iu1,iu2     |first and last index of input array,
  !  ju1,ju2     |cell centered,
  !  ku1,ku2     |including buffer cells.
  !------------------------------------------------------------------------
    real(dp),dimension(nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3), &
                                                            intent(in)::uin_gas
    real(dp),dimension(nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ncrvar), &
                                                            intent(in)::uin_cr
    real(dp),dimension(nvector,iu1:iu2,ju1:ju2,ku1:ku2,ndim)::   lmax
    integer,intent(in)::iGrp, ngrid, ilevel
    real(dp),intent(in)::dt !----------------------------------------------
    real(dp)::scale_kappa, dx, dx_loc, scale, Ecr, va
    integer::icrE, i, j, k, n, nx_loc, idim, nedge
    real(dp),dimension(1:3)::B_field, gradEcr, Dcr_vec
    real(dp)::norm,bdotgradE,cosp,sinp,cost,sint,bxby,Dcr_dir
  !------------------------------------------------------------------------
    icrE = iCRu+(ndim+1)*(iGrp-1) ! starting index of cr variables (base 1)

    ! Cell width in ilevel
    dx=0.5D0**ilevel
    nx_loc=(icoarse_max-icoarse_min+1)
    scale=boxlen/dble(nx_loc)
    dx_loc=dx*scale
    do k=kfcr1,kf2                                !
    do j=jfcr1,jf2                                !  Loop each cell in grid
    do i=ifcr1,if2                                !

    nedge=0                   ! Check if we're at a corner and if so, cycle
    if(mod(i,3).eq.0) nedge=nedge+1
    if(ndim.gt.1 .and. mod(j,3).eq.0) nedge=nedge+1
    if(ndim.gt.2 .and. mod(k,3).eq.0) nedge=nedge+1
    if(nedge.ge.2) cycle

    do n=1,ngrid                                     ! Loop buffer of grids

       Ecr    = uin_cr(n, i, j, k, icrE)

       ! Magnetic field, needed to rotate Dcr and for bdotgradEcr
       norm=0.
       do idim=1,3
          B_field(idim) = 0.5 * &
            (uin_gas(n, i, j, k, 5+idim) + uin_gas(n, i, j, k, nvar+idim) )
          norm = norm + B_field(idim)**2
       end do
       norm = max(sqrt(norm), 1d-30)

       bxby = sqrt(B_field(1)**2+B_field(2)**2)
       if(norm.gt.1e-10) then
          sint = bxby/norm
          cost = B_field(3)/norm
       else
          sint = 1d0
          cost = 0d0
       endif
       if(bxby.gt.1e-10) then
          sinp = B_field(2)/bxby
          cosp = B_field(1)/bxby
       else
          sinp = 0d0
          cosp = 1d0
       endif
       B_field = B_field/norm

       va=0.
       if(mom_streaming_diffusion) va=norm/sqrt(uin_gas(n,i,j,k,1))
       if(mom_streaming_diffusion .and. v_alfven.gt.0.0) va = v_alfven

       ! Calculate grad Pcr
       gradEcr(1) = (uin_cr(n, i+1, j, k, icrE)-uin_cr(n, i-1, j, k, icrE)) &
                      / (2d0*dx_loc)
#if NDIM>1
       gradEcr(2) = (uin_cr(n, i, j+1, k, icrE)-uin_cr(n, i, j-1, k, icrE)) &
                      / (2d0*dx_loc)
#endif
#if NDIM>2
       gradEcr(3) = (uin_cr(n, i, j, k+1, icrE)-uin_cr(n, i, j, k-1, icrE)) &
                      / (2d0*dx_loc)
#endif

       ! Calculate B dot grad Pcr
       bdotgradE = 0.
       do idim=1,ndim
          bdotgradE = bdotgradE + B_field(idim) * gradEcr(idim)
       enddo

       ! Diffusion, eq 10 in JO17
       Dcr_vec = (/ DCR_code(igrp), DCR_code(igrp)*Dcr_perp_factor(iGrp), DCR_code(igrp)*Dcr_perp_factor(iGrp) /)
       if(mom_streaming_diffusion) &
          Dcr_vec(1) = Dcr_vec(1) + &
              min(DCRmax_code, 3./max(abs(bdotgradE),1d-50) * va * gamma_cr(iGrp) * max(Ecr,smallcr))

       ! Rotate Dcr_vec so it is parallel with B, hence
       ! describing Dcr in the simulation coordinate system
       ! (instead of the B coordinate system)
       call invrotatevec(sint, cost, sinp, cosp, Dcr_vec(1), Dcr_vec(2), Dcr_vec(3))

       ! Calculate wavespeeds
       Dcr_dir = abs(Dcr_vec(1)) ! x component of rotated Dcr
       lmax(n,i,j,k,1) = &
          cmp_cr_lmax(dx_loc, Dcr_dir, cr_vmax(ilevel), dt)

#if NDIM>1
       Dcr_dir = abs(Dcr_vec(2)) ! y component of rotated Dcr
       lmax(n,i,j,k,2) = &
          cmp_cr_lmax(dx_loc, Dcr_dir, cr_vmax(ilevel), dt)
#endif

#if NDIM>2
       Dcr_dir = abs(Dcr_vec(3)) ! z component of rotated Dcr
       lmax(n,i,j,k,3) = &
          cmp_cr_lmax(dx_loc, Dcr_dir, cr_vmax(ilevel), dt)
#endif

    end do
    end do
    end do
    end do

END SUBROUTINE cmp_cr_wavespeeds

!************************************************************************
FUNCTION cmp_cr_lmax(dx_loc, dcoeff, vmax, dt)

! Compute maximum local wavespeed
!------------------------------------------------------------------------
  real(dp)::dx_loc, cmp_cr_lmax, vmax, dt
  real(dp)::tau, dcoeff, r_factor
!------------------------------------------------------------------------
  !tau = cr_f_taucell * dx_loc / dcoeff * vmax * sqrt(1.5)
  tau = cr_f_taucell * 0.5d0 * dx_loc**2 / dcoeff / dt
  if(tau.lt.1e-3) then
    r_factor = sqrt((1.0 - 0.5*tau**2))
  else
    r_factor = sqrt((1.-exp(-min(tau,10.)**2))/min(tau,1e8)**2) ! Capital R on p 6 in YP17
  endif
  if(isotropic_pressure) then
     cmp_cr_lmax = r_factor * vmax / sqrt(3d0)
  else
     cmp_cr_lmax = r_factor * vmax
  endif

END FUNCTION cmp_cr_lmax

!************************************************************************
FUNCTION cmp_cr_face(fdn, fup, udn, uup, lminus, lplus)

! Compute HLLE intercell fluxes for all (four) CR variables.
! fdn    => flux function in the cell downwards from the border
! fup    => flux function in the cell upwards from the border
! udn    => state (cr density and flux) downwards from the border
! uup    => state (cr density and flux) upwards from the border
! lminus => minimum intercell wavespeed
! lplus  => maximum intercell wavespeed
! returns      flux vector for the given state variables, i.e. line nr dim
!              in the 3*4 flux function tensor
!------------------------------------------------------------------------
  real(dp),dimension(nDim+1)::fdn, fup, udn, uup, cmp_cr_face
  real(dp)::lminus, lplus, coeff, llmax
!------------------------------------------------------------------------
  if (cr_HLLE) then
    coeff = 0D0
    if (abs(lplus - lminus) > 1D-20) coeff = 0.5D0 * (lplus + lminus) / (lplus - lminus)
    cmp_cr_face = 0.5D0 * (fdn + fup - lminus * udn - lplus * uup) &
                  + coeff * (fdn - fup - lminus * udn + lplus * uup)
  else ! Lax Friedrich
    llmax = max(abs(lplus), abs(lminus))
    cmp_cr_face = 0.5 * (fdn + fup - llmax * (uup - udn)) !LLF flux
  endif
END FUNCTION cmp_cr_face

!************************************************************************
FUNCTION minmod(left, right)

! Minmod interpolation function for intercell fluxes
  real(dp),dimension(ndim+1)::left,right,minmod
  integer::i
!------------------------------------------------------------------------
  do i=1,ndim+1
    if (abs(left(i) ) .gt. abs(right(i))) minmod(i)=right(i)
    if (abs(right(i)) .ge. abs(left(i) )) minmod(i)=left(i)
    if (left(i)*right(i) .le. 0.0d0     ) minmod(i)=0.0d0
  end do
END FUNCTION minmod

!************************************************************************
FUNCTION vminmod(left, right)

! Minmod interpolation function for velocity
  real(dp)::left,right,vminmod
!------------------------------------------------------------------------
  if (abs(left ) .gt. abs(right)) vminmod=right
  if (abs(right) .ge. abs(left )) vminmod=left
  if (left*right .le. 0.0d0     ) vminmod=0.0d0
END FUNCTION vminmod

!************************************************************************
SUBROUTINE cmp_cr_faces(uin_gas,uin_cr,iFlx,dx,dt,iGrp,ngrid,ilevel)

!  Compute intercell fluxes for one CR group in all dimensions,
!  using the Eddington tensor with the Yiang+Peng'17 closure relation.
!  The intercell fluxes are the right-hand sides of the equation:
!      dq/dt = - nabla dot f   (eq 12 in YP17)
!  where q=[Ecr, Fx/ccr^2, Fy/ccr^2, Fz/ccr^2], ccr the reduced wavespeed
!  and f the Eddington pressure tensor. A flux at index i,j,k represents
!  flux across the lower faces of that cell, i.e. at i-1/2 etc.
!
!  inputs/outputs
!  uin_gas     => input gas (rho,vel,B) states
!  uin_cr      => input CR (E,F) states
!  iFlx       <=  return fluxes in the 3 coord directions.
!  dx          => cell width
!  dt          => time step
!  iGrp        => CR group number
!  ngrid       => number of sub-grids
!  ilevel      => level being updated
!
!  other vars
!  iu1,iu2     |First and last index of input array,
!  ju1,ju2     |cell centered,
!  ku1,ku2     |including buffer cells (6x6x6).
!  if1,if2     |First and last index of output array,
!  jf1,jf2     |edge centered, for active
!  kf1,kf2     |cells only (3x3x3).
!------------------------------------------------------------------------
  use amr_parameters
  use amr_commons
  use const
  implicit none
  real(dp),dimension(nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3)::uin_gas
  real(dp),dimension(nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ncrvar)::uin_cr
  real(dp),dimension(nvector,if1:if2,jf1:jf2,kf1:kf2,1:ncrvar,1:ndim)::iFlx
  real(dp)::dx, dt
  integer,intent(in)::iGrp, nGrid, ilevel
  integer::iP0, iP1
  real(dp),save, &                                     !   Central fluxes
           dimension(nvector,iu1:iu2,ju1:ju2,ku1:ku2, ndim+1, ndim)::cFlx
  real(dp),save, &                                     ! Cell wavespeeds
           dimension(nvector,iu1:iu2,ju1:ju2,ku1:ku2,ndim)::  lmax
  ! Upwards and downwards fluxes and states of the group
  real(dp),dimension(nDim+1),save:: fdn, fup, udn, uup
  real(dp):: lminus, lplus                        ! Intercell wavespeeds
  real(dp)::dtdx, prod(ndim+1)
  integer ::i, j, k, n
  real(dp),dimension(ndim+1),save::slopeLM,slopeRM,slopeM
  real(dp),dimension(ndim+1),save::slopeLL,slopeL
  real(dp):: vup,vdn,meanadv,meandiffv,aup,adn
  real(dp),dimension(nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim)::bfield
  integer::idim
!------------------------------------------------------------------------
  iP0 = 1+(iGrp-1)*(ndim+1)            ! Index of CR group energy density
  iP1 = iP0+nDim
  ! Magnetic field, needed for M1; in the default P1 path (isotropic_pressure)
  ! cmp_cr_flux_tensors never reads bfield, so skip building the full stencil.
  if(.not.isotropic_pressure) then
     do idim=1,ndim
        bfield(:,:,:,:,idim) = 0.5*(uin_gas(:,:,:,:,5+idim)+uin_gas(:,:,:,:,nvar+idim))
     end do
  endif

  ! compute flux tensors for all the cells with correction
  call cmp_cr_flux_tensors(uin_gas, uin_cr, iGrp, ngrid, cFlx, cr_vmax(ilevel),bfield)!  flux tensors

  ! Wavespeeds in each cell
  call cmp_cr_wavespeeds(uin_gas, uin_cr, iGrp, ngrid, lmax, ilevel, dt)

  ! Solve for 1D flux in X direction
  !----------------------------------------------------------------------
  dtdx=dt/dx
  do i=if1,if2                                 !
  do j=jf1,jf2                                 !        each cell in grid
  do k=kf1,kf2                                 !
     if(ndim.gt.1 .and. j.eq.jf2) cycle
     if(ndim.gt.2 .and. k.eq.kf2) cycle
     do n=1,ngrid                              ! <- buffer of grids
        fdn = cFlx(n,  i-1, j, k, :, 1    )    !
        fup = cFlx(n,  i,   j, k, :, 1    )    !  upwards and downwards
        udn = uin_cr( n,  i-1, j, k, iP0:iP1 ) !  conditions
        uup = uin_cr( n,  i,   j, k, iP0:iP1 ) !
        vdn  = uin_gas( n,  i-1, j, k, 2) / uin_gas(n,i-1,j,k,1) ! left velocity
        vup  = uin_gas( n,  i,   j, k, 2) / uin_gas(n,i  ,j,k,1) ! right velocity

        if(cr_interpolation) then ! Second-order interpolation with using Van-Leer slope Limiter
        ! interpolation of U
        slopeLM = (fup-fdn)/dx
        slopeRM = (cFlx(n, i+1, j, k, :, 1) - fup)/dx
        prod = slopeLM*slopeRM
        slopeM=0.
        where(prod.gt.0.) slopeM=2.*prod/(slopeLM+slopeRM)
        slopeLL = (fdn - cFlx(n, i-2, j, k, :, 1))/dx
        prod = slopeLL*slopeLM
        slopeL=0.
        where(prod.gt.0) slopeL=2.*prod/(slopeLL+slopeLM)
        fdn = fdn+slopeL*0.5d0*dx
        fup = fup-slopeM*0.5d0*dx

        ! interpolation of F
        slopeLM = (uup-udn)/dx
        slopeRM = (uin_cr(n, i+1, j, k, iP0:iP1) - uup)/dx
        prod = slopeLM*slopeRM
        slopeM=0.
        where(prod.gt.0) slopeM=2.*prod/(slopeLM+slopeRM)
        slopeLL = (udn - uin_cr(n, i-2, j, k, iP0:iP1 ))/dx
        prod = slopeLL*slopeLM
        slopeL=0.
        where(prod.gt.0.) slopeL=2.*prod/(slopeLL+slopeLM)
        udn = udn+slopeL*0.5d0*dx
        uup = uup-slopeM*0.5d0*dx

        else if(cr_use_minmod) then ! Second-order interpolation with using minmod slope Limiter
            slopeLM = (fup-fdn)/dx
            slopeRM = (cFlx(n, i+1, j, k, :, 1) - fup)/dx
            slopeM  = minmod(slopeLM,slopeRM)
            slopeLL = (fdn - cFlx(n, i-2, j, k, :, 1))/dx
            slopeL  = minmod(slopeLL,slopeLM)
            fdn = fdn+slopeL*0.5d0*dx
            fup = fup-slopeM*0.5d0*dx

            slopeLM = (uup-udn)/dx
            slopeRM = (uin_cr(n, i+1, j, k, iP0:iP1) - uup)/dx
            slopeM  = minmod(slopeLM,slopeRM)
            slopeLL = (udn - uin_cr(n, i-2, j, k, iP0:iP1 ))/dx
            slopeL  = minmod(slopeLL,slopeLM)
            udn = udn+slopeL*0.5d0*dx
            uup = uup-slopeM*0.5d0*dx
        endif
        meanadv = 0.5*(vdn+vup)
        meandiffv = 0.5*( lmax(n,i-1,j,k,1) + lmax(n,i,j,k,1) )
        adn = min(meanadv-meandiffv, vdn-lmax(n,i-1,j,k,1))
        adn = max(adn,-cr_vmax(ilevel)*sqrt(1./3.))
        aup = max(meanadv+meandiffv, vup+lmax(n,i,j,k,1))
        aup = min(aup,cr_vmax(ilevel)*sqrt(1./3.))
        lminus = min(adn,0.)
        lplus = max(aup,0.)

        iFlx(n, i, j, k, iP0:iP1, 1)=&
             cmp_cr_face( fdn, fup, udn, uup, lminus, lplus)*dtdx
     end do
  end do
  end do
  end do

  ! Solve for 1D flux in Y direction
  !----------------------------------------------------------------------
#if NDIM>1
  dtdx=dt/dx
  do i=if1,if2-1
  do j=jf1,jf2
  do k=kf1,kf2
     if(ndim.gt.2 .and. k.eq.kf2) cycle
     do n=1,ngrid
           fdn = cFlx(n,  i, j-1, k, :, 2    )
           fup = cFlx(n,  i, j,   k, :, 2    )
           udn = uin_cr( n,  i, j-1, k, iP0:iP1 )
           uup = uin_cr( n,  i, j,   k, iP0:iP1 )
           vdn  = uin_gas( n,  i, j-1, k,3) / uin_gas(n,i,j-1,k,1) ! left velocity
           vup  = uin_gas( n,  i ,j,   k,3) / uin_gas(n,i,j,  k,1) ! right velocity

           if(cr_interpolation) then ! Second-order interpolation with using Van-Leer slope Limiter
           ! interpolation of U
           slopeLM = (fup-fdn)/dx
           slopeRM = (cFlx(n, i, j+1, k, :, 2) - fup)/dx
           prod = slopeLM*slopeRM
           slopeM=0.
           where(prod.gt.0.) slopeM=2.*prod/(slopeLM+slopeRM)
           slopeLL = (fdn - cFlx(n, i, j-2, k, :, 2))/dx
           prod = slopeLL*slopeLM
           slopeL=0.
           where(prod.gt.0) slopeL=2.*prod/(slopeLL+slopeLM)
           fdn = fdn+slopeL*0.5d0*dx
           fup = fup-slopeM*0.5d0*dx

           ! interpolation of F
           slopeLM = (uup-udn)/dx
           slopeRM = (uin_cr(n, i, j+1, k, iP0:iP1) - uup)/dx
           prod = slopeLM*slopeRM
           slopeM=0.
           where(prod.gt.0) slopeM=2.*prod/(slopeLM+slopeRM)
           slopeLL = (udn - uin_cr(n, i, j-2, k, iP0:iP1 ))/dx
           prod = slopeLL*slopeLM
           slopeL=0.
           where(prod.gt.0.) slopeL=2.*prod/(slopeLL+slopeLM)
           udn = udn+slopeL*0.5d0*dx
           uup = uup-slopeM*0.5d0*dx

           else if(cr_use_minmod) then ! Second-order interpolation with using minmod slope Limiter
              slopeLM = (fup-fdn)/dx
              slopeRM = (cFlx(n, i, j+1, k, :, 2) - fup)/dx
              slopeM  = minmod(slopeLM,slopeRM)
              slopeLL = (fdn - cFlx(n, i, j-2, k, :, 2))/dx
              slopeL  = minmod(slopeLL,slopeLM)
              fdn = fdn+slopeL*0.5d0*dx
              fup = fup-slopeM*0.5d0*dx

              slopeLM = (uup-udn)/dx
              slopeRM = (uin_cr(n, i, j+1, k, iP0:iP1) - uup)/dx
              slopeM  = minmod(slopeLM,slopeRM)
              slopeLL = (udn - uin_cr(n, i, j-2, k, iP0:iP1))/dx
              slopeL  = minmod(slopeLL,slopeLM)
              udn = udn+slopeL*0.5d0*dx
              uup = uup-slopeM*0.5d0*dx
           endif
           meanadv = 0.5*(vdn+vup)
           meandiffv = 0.5*( lmax(n,i,j-1,k,2) + lmax(n,i,j,k,2) )
           adn = min(meanadv-meandiffv, vdn-lmax(n,i,j-1,k,2))
           adn = max(adn,-cr_vmax(ilevel)*sqrt(1./3.))
           aup = max(meanadv+meandiffv, vup+lmax(n,i,j,k,2))
           aup = min(aup,cr_vmax(ilevel)*sqrt(1./3.))
           lminus = min(adn,0.)
           lplus = max(aup,0.)

           iFlx(n, i, j, k, iP0:iP1, 2)=&
                cmp_cr_face( fdn, fup, udn, uup, lminus, lplus)*dtdx
     end do
  end do
  end do
  end do
#endif

  ! Solve for 1D flux in Z direction
  !----------------------------------------------------------------------
#if NDIM>2
  dtdx=dt/dx
  do i=if1,if2-1
  do j=jf1,jf2-1
  do k=kf1,kf2
     do n=1,ngrid
           fdn = cFlx(n,  i, j, k-1, :, 3    )
           fup = cFlx(n,  i, j, k,   :, 3    )
           udn = uin_cr( n,  i, j, k-1, iP0:iP1 )
           uup = uin_cr( n,  i, j, k,   iP0:iP1 )
           vdn  = uin_gas( n,  i, j, k-1, 4) / uin_gas(n,i,  j,k-1,1) ! left velocity
           vup  = uin_gas( n,  i ,j, k,   4) / uin_gas(n,i  ,j,k,  1) ! right velocity

           if(cr_interpolation) then ! Second-order interpolation with using Van-Leer slope Limiter
           ! interpolation of U
           slopeLM = (fup-fdn)/dx
           slopeRM = (cFlx(n, i, j, k+1, :, 3) - fup)/dx
           prod = slopeLM*slopeRM
           slopeM=0.
           where(prod.gt.0.) slopeM=2.*prod/(slopeLM+slopeRM)
           slopeLL = (fdn - cFlx(n, i, j, k-2, :, 3))/dx
           prod = slopeLL*slopeLM
           slopeL=0.
           where(prod.gt.0) slopeL=2.*prod/(slopeLL+slopeLM)
           fdn = fdn+slopeL*0.5d0*dx
           fup = fup-slopeM*0.5d0*dx

           ! interpolation of F
           slopeLM = (uup-udn)/dx
           slopeRM = (uin_cr(n, i, j, k+1, iP0:iP1) - uup)/dx
           prod = slopeLM*slopeRM
           slopeM=0.
           where(prod.gt.0) slopeM=2.*prod/(slopeLM+slopeRM)
           slopeLL = (udn - uin_cr(n, i, j, k-2, iP0:iP1 ))/dx
           prod = slopeLL*slopeLM
           slopeL=0.
           where(prod.gt.0.) slopeL=2.*prod/(slopeLL+slopeLM)
           udn = udn+slopeL*0.5d0*dx
           uup = uup-slopeM*0.5d0*dx

           else if(cr_use_minmod) then ! Second-order interpolation with using minmod slope Limiter
              slopeLM = (fup-fdn)/dx
              slopeRM = (cFlx(n, i, j, k+1, :, 3) - fup)/dx
              slopeM  = minmod(slopeLM,slopeRM)
              slopeLL = (fdn - cFlx(n, i, j, k-2, :, 3))/dx
              slopeL  = minmod(slopeLL,slopeLM)
              fdn = fdn+slopeL*0.5d0*dx
              fup = fup-slopeM*0.5d0*dx

              slopeLM = (uup-udn)/dx
              slopeRM = (uin_cr(n, i, j, k+1, iP0:iP1) - uup)/dx
              slopeM  = minmod(slopeLM,slopeRM)
              slopeLL = (udn - uin_cr(n, i, j, k-2, iP0:iP1))/dx
              slopeL  = minmod(slopeLL,slopeLM)
              udn = udn+slopeL*0.5d0*dx
              uup = uup-slopeM*0.5d0*dx
           endif

           meanadv = 0.5*(vdn+vup)
           meandiffv = 0.5*( lmax(n,i,j,k-1,3) + lmax(n,i,j,k,3) )
           adn = min(meanadv-meandiffv, vdn-lmax(n,i,j,k-1,3))
           adn = max(adn,-cr_vmax(ilevel)*sqrt(1./3.))
           aup = max(meanadv+meandiffv, vup+lmax(n,i,j,k,3))
           aup = min(aup,cr_vmax(ilevel)*sqrt(1./3.))
           lminus = min(adn,0.)
           lplus  = max(aup,0.)

           iFlx(n, i, j, k, iP0:iP1, 3)=&
                cmp_cr_face( fdn, fup, udn, uup, lminus, lplus)*dtdx
      end do
  end do
  end do
  end do
#endif

end subroutine cmp_cr_faces

!************************************************************************
PURE SUBROUTINE rotatevec(sint, cost, sinp, cosp, v1, v2, v3)
  !  Rotate vector v by t=theta and p=phi
  !  i.e. rotate to the local coordinate system from theta, phi.
  !  Hence the x-component of the result is the component of v parallel
  !  with the theta,phi vector.
  !------------------------------------------------------------------------
    implicit none
    real(dp),intent(in):: sint, cost, sinp, cosp
    real(dp),intent(inout)::v1,v2,v3
    real(dp)::newv1, newv3
  !------------------------------------------------------------------------
    ! First apply R1, then apply R2
    newv1 =  cosp * v1 + sinp * v2
    v2 = -sinp * v1 + cosp * v2

    ! Now apply R2
    v1 =  sint * newv1 + cost * v3
    newv3 = -cost * newv1 + sint * v3
    v3 = newv3
END SUBROUTINE rotatevec

!************************************************************************
PURE SUBROUTINE invrotatevec(sint, cost, sinp, cosp, v1, v2, v3)
  !  Inverse-rotate vector v by t=theta and p=phi,
  !  i.e. rotate v onto theta, pi
  !
  !------------------------------------------------------------------------
    implicit none
    real(dp),intent(in):: sint, cost, sinp, cosp
    real(dp),intent(inout)::v1,v2,v3
    real(dp)::newv1, newv2
  !------------------------------------------------------------------------
    ! First apply R2^-1, then apply R1^-1
    newv1 = sint * v1 - cost * v3
    v3 = cost * v1 + sint * v3

    ! now apply R1^-1
    v1 = cosp * newv1 - sinp * v2
    newv2 = sinp * newv1 + cosp * v2
    v2 = newv2
END SUBROUTINE invrotatevec

END MODULE cr_flux_module
