MODULE cr_flux_module

  ! Two-moment cosmic-ray closure and intercell-flux module.
  ! This is a PURE flux/closure module: no source terms, no NENER.

  use amr_commons
  use hydro_commons
  use hydro_parameters         ! iu1..ku2, if1..kf2, nvar, smallr
  use cr_parameters            ! ncr_groups, ncrvar, Ecr_idx, gamma_cr, cr_vmax, ...
  implicit none

  private   ! default

  ! CR-specific stencil-loop bounds used by cmp_cr_wavespeeds.
  integer,parameter::ifcr1=0
  integer,parameter::jfcr1=1-ndim/2
  integer,parameter::kfcr1=1-ndim/3

  public cmp_cr_faces, cr_limit_flux, rotatevec, invrotatevec

CONTAINS

!************************************************************************
subroutine cmp_cr_flux_tensors(uin_cr, iGrp, nGrid, ftens, vmax, bfield)

! Compute central fluxes for a CR group, for each cell in a vector
! of grids.
! The flux tensor is a three by four tensor (2*3 and 1*2 in 1D and 2D,
! respectively) where the first row is CR flux (x,y,z) and
! the other three rows compose the Eddington tensor (Jiang & Oh 2018)
! input/output:
! uin_cr    => CR (E,F) variables of all cells in a vector of grids
! igrp      => CR group number
! ngrid     => Number of 'valid' grids in uin.
! ftens     <=  Group flux tensors for all the cells.
! vmax      => reduced light speed at this level
! bfield    => cell-centred B field (only referenced for the M1 closure)
!------------------------------------------------------------------------
  real(dp),dimension(nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ncrvar),intent(in)::uin_cr
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nDim+1,1:ndim),intent(inout)::ftens
  real(dp),dimension(nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3),intent(in)::bfield
  integer,intent(in)::iGrp, nGrid
  real(dp),intent(in)::vmax !---------------------------------------------
  real(dp),dimension(1:ndim)::crflux
  real(dp)::Ecr, pdiag
  integer::i, j, k, idim, jdim, n, icrE, nedge
  real(dp)::mu1,mu2,chi,b_norm2,crflux_norm
  real(dp)::vmax2,Ecr2,aniso_term
  !------------------------------------------------------------------------

  icrE = Ecr_idx(iGrp) ! starting index of cr variables (base 1)
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
        Ecr =   max(uin_cr(n, i, j, k, icrE), cr_efloor)   ! CR density (floored so 0 background is safe in M1)
        crflux = uin_cr(n,i,j,k,icrE+1 : icrE+ndim) ! CR flux vector
        ftens(n,i,j,k,1,1:ndim)= crflux  !   First row is CR flux
        ! Rest is Eddington tensor
        if(cr_isotropic_pressure) then
           ! P1 closure: isotropic pressure, diagonal tensor E*c^2/3
           ftens(n,i,j,k,2:ndim+1,1:ndim) = 0d0
           pdiag = Ecr*vmax**2/3d0
           do idim = 1, ndim
              ftens(n,i,j,k,idim+1,idim) = pdiag
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
SUBROUTINE cmp_cr_wavespeeds(uin_gas, uin_cr, bfield, iGrp, ngrid, lmax, ilevel, dt)

!  Compute CR wavespeeds for given vector of sub-grids.
!
!  inputs/outputs
!  uin_gas     => input gas (rho,vel) cell states
!  uin_cr      => input CR (E,F) cell states
!  bfield      => cell-centred B field (all 3 components)
!  iGrp        => CR group number
!  ngrid       => number of sub-grids of 3^ndim cells
!  ilevel      => current refinement level
!  dt          => current CR timestep length
!  lmax       <=  return maximum cell wavespeeds
!
!  other vars
!  iu1,iu2     |first and last index of input array,
!  ju1,ju2     |cell centered,
!  ku1,ku2     |including buffer cells.
!------------------------------------------------------------------------
    real(dp),dimension(nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar_all), &
                                                            intent(in)::uin_gas
    real(dp),dimension(nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ncrvar), &
                                                            intent(in)::uin_cr
    real(dp),dimension(nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3), &
                                                            intent(in)::bfield
    real(dp),dimension(nvector,iu1:iu2,ju1:ju2,ku1:ku2,ndim), &
                                                            intent(inout)::lmax
    integer,intent(in)::iGrp, ngrid, ilevel
    real(dp),intent(in)::dt !----------------------------------------------
    real(dp)::dx, dx_loc, scale, Ecr, va
    integer::icrE, i, j, k, n, nx_loc, idim, nedge
    real(dp),dimension(1:3)::B_field, gradEcr, Dcr_vec
    real(dp)::norm,bdotgradE,cosp,sinp,cost,sint,bxby
    real(dp)::Dcr_par,Dcr_perp
  !------------------------------------------------------------------------
    icrE = Ecr_idx(iGrp) ! starting index of cr variables (base 1)

    ! Diffusion coefficients along/across B, eq 10 in JO18
    Dcr_par  = DCR_code(iGrp)
    Dcr_perp = DCR_code(iGrp)*Dcr_perp_factor(iGrp)

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

       ! Magnetic field, needed to rotate Dcr and for bdotgradEcr
       norm=0.
       do idim=1,3
          B_field(idim) = bfield(n,i,j,k,idim)
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

       Dcr_vec = (/ Dcr_par, Dcr_perp, Dcr_perp /)

       if(cr_streaming_diffusion) then
          ! Streaming as an effective diffusivity 3*va*gamma_cr*Ecr/|b.gradEcr|
          va = norm/sqrt(uin_gas(n,i,j,k,1))
          if(cr_v_alfven.gt.0.0) va = cr_v_alfven
          B_field = B_field/norm

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

          Ecr = uin_cr(n, i, j, k, icrE)
          Dcr_vec(1) = Dcr_vec(1) + &
              min(DCRmax_code, 3./max(abs(bdotgradE),1d-50) * va * gamma_cr(iGrp) * max(Ecr,cr_efloor))
       endif

       ! Rotate Dcr_vec so it is parallel with B, hence
       ! describing Dcr in the simulation coordinate system
       ! (instead of the B coordinate system)
       call invrotatevec(sint, cost, sinp, cosp, Dcr_vec(1), Dcr_vec(2), Dcr_vec(3))

       ! Calculate wavespeeds
       do idim=1,ndim
          lmax(n,i,j,k,idim) = &
             cmp_cr_lmax(dx_loc, abs(Dcr_vec(idim)), cr_vmax(ilevel), dt)
       end do

    end do
    end do
    end do
    end do

END SUBROUTINE cmp_cr_wavespeeds

!************************************************************************
PURE FUNCTION cmp_cr_lmax(dx_loc, dcoeff, vmax, dt)

! Compute maximum local wavespeed
!------------------------------------------------------------------------
  real(dp),intent(in)::dx_loc, vmax, dt, dcoeff
  real(dp)::cmp_cr_lmax
  real(dp)::tau, r_factor
!------------------------------------------------------------------------
  tau = cr_f_taucell * 0.5d0 * dx_loc**2 / dcoeff / dt
  if(tau.lt.1e-3) then
    r_factor = sqrt((1.0 - 0.5*tau**2))
  else
    r_factor = sqrt((1.-exp(-min(tau,10.)**2))/min(tau,1e8)**2) ! Capital R on p 6 in YP17
  endif
  if(cr_isotropic_pressure) then
     cmp_cr_lmax = r_factor * vmax / sqrt(3d0)
  else
     cmp_cr_lmax = r_factor * vmax
  endif

END FUNCTION cmp_cr_lmax

!************************************************************************
PURE FUNCTION cmp_cr_face(fdn, fup, udn, uup, lminus, lplus)

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
  real(dp),dimension(nDim+1),intent(in)::fdn, fup, udn, uup
  real(dp),intent(in)::lminus, lplus
  real(dp),dimension(nDim+1)::cmp_cr_face
  real(dp)::coeff, llmax
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
PURE SUBROUTINE cr_recon_face(qmm, qdn, qup, qpp, dx, rdn, rup)

! Slope-limited (harmonic-mean) linear reconstruction at the qdn|qup
! interface, from the four cells qmm,qdn,qup,qpp centred at i-2,i-1,i,i+1
! relative to the face at i-1/2.
! rdn/rup <= the reconstructed states just below/above the face.
!------------------------------------------------------------------------
  real(dp),dimension(nDim+1),intent(in)::qmm, qdn, qup, qpp
  real(dp),intent(in)::dx
  real(dp),dimension(nDim+1),intent(out)::rdn, rup
  real(dp),dimension(nDim+1)::slopeLL,slopeLM,slopeRM,slopeL,slopeM,prod
!------------------------------------------------------------------------
  slopeLM = (qup-qdn)/dx
  slopeRM = (qpp-qup)/dx
  prod = slopeLM*slopeRM
  slopeM = 0.
  where(prod.gt.0.) slopeM = 2.*prod/(slopeLM+slopeRM)
  slopeLL = (qdn-qmm)/dx
  prod = slopeLL*slopeLM
  slopeL = 0.
  where(prod.gt.0.) slopeL = 2.*prod/(slopeLL+slopeLM)
  rdn = qdn + slopeL*0.5d0*dx
  rup = qup - slopeM*0.5d0*dx
END SUBROUTINE cr_recon_face

!************************************************************************
SUBROUTINE cmp_cr_faces(uin_gas,uin_cr,iFlx,dx,dt,iGrp,ngrid,ilevel)

!  Compute intercell fluxes for one CR group in all dimensions,
!  using the Eddington tensor with the Jiang & Oh 2018 closure relation.
!  The intercell fluxes are the right-hand sides of the equation:
!      dq/dt = - nabla dot f   (eq 12 in JO18)
!  where q=[Ecr, Fx/ccr^2, Fy/ccr^2, Fz/ccr^2], ccr the reduced wavespeed
!  and f the Eddington pressure tensor. A flux at index i,j,k represents
!  flux across the lower faces of that cell, i.e. at i-1/2 etc.
!  All directions share one loop: the sweep direction spans faces
!  if1:if2, transverse directions span active cells only.
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
  real(dp),dimension(nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar_all),intent(in)::uin_gas
  real(dp),dimension(nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ncrvar),intent(in)::uin_cr
  real(dp),dimension(nvector,if1:if2,jf1:jf2,kf1:kf2,1:ncrvar,1:ndim),intent(inout)::iFlx
  real(dp),intent(in)::dx, dt
  integer,intent(in)::iGrp, nGrid, ilevel
  real(dp),save, &                                     !   Central fluxes
           dimension(nvector,iu1:iu2,ju1:ju2,ku1:ku2, ndim+1, ndim)::cFlx
  real(dp),save, &                                     ! Cell wavespeeds
           dimension(nvector,iu1:iu2,ju1:ju2,ku1:ku2,ndim)::  lmax
  real(dp),save, &                                     ! Cell-centred B
           dimension(nvector,iu1:iu2,ju1:ju2,ku1:ku2,3)::  bfield
  ! Reconstructed fluxes and states below/above the face
  real(dp),dimension(nDim+1):: fdn, fup, udn, uup
  real(dp):: lminus, lplus                        ! Intercell wavespeeds
  real(dp):: vup,vdn,meanadv,meandiffv,aup,adn
  real(dp):: dtdx,vclamp
  integer :: iP0,iP1,i,j,k,n,idim,i0,j0,k0,ihi,jhi,khi
!------------------------------------------------------------------------
  iP0 = Ecr_idx(iGrp)
  iP1 = iP0+nDim
  ! Cell-centred B field, all 3 components: cmp_cr_wavespeeds needs the full
  ! vector for the Dcr rotation angles even in 1D/2D.
  do idim=1,3
     bfield(:,:,:,:,idim) = 0.5*(uin_gas(:,:,:,:,5+idim)+uin_gas(:,:,:,:,nvar+idim))
  end do

  ! Central flux tensors for all cells
  call cmp_cr_flux_tensors(uin_cr, iGrp, ngrid, cFlx, cr_vmax(ilevel), bfield)

  ! Wavespeeds in each cell
  call cmp_cr_wavespeeds(uin_gas, uin_cr, bfield, iGrp, ngrid, lmax, ilevel, dt)

  ! Solve for the 1D flux at the lower cell faces in each direction
  !----------------------------------------------------------------------
  dtdx=dt/dx
  vclamp=cr_vmax(ilevel)*sqrt(1./3.) ! single-precision literal kept for bit-compatibility
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     ihi=if2; if(idim/=1)ihi=if2-1
     jhi=jf2; if(ndim>1 .and. idim/=2)jhi=jf2-1
     khi=kf2; if(ndim>2 .and. idim/=3)khi=kf2-1
     do k=kf1,khi
     do j=jf1,jhi
     do i=if1,ihi
        do n=1,ngrid                              ! <- buffer of grids
           ! Limited linear reconstruction of central fluxes and CR states
           call cr_recon_face(cFlx(n,i-2*i0,j-2*j0,k-2*k0,:,idim), &
                &             cFlx(n,i-i0  ,j-j0  ,k-k0  ,:,idim), &
                &             cFlx(n,i     ,j     ,k     ,:,idim), &
                &             cFlx(n,i+i0  ,j+j0  ,k+k0  ,:,idim), dx, fdn, fup)
           call cr_recon_face(uin_cr(n,i-2*i0,j-2*j0,k-2*k0,iP0:iP1), &
                &             uin_cr(n,i-i0  ,j-j0  ,k-k0  ,iP0:iP1), &
                &             uin_cr(n,i     ,j     ,k     ,iP0:iP1), &
                &             uin_cr(n,i+i0  ,j+j0  ,k+k0  ,iP0:iP1), dx, udn, uup)

           vdn = uin_gas(n,i-i0,j-j0,k-k0,1+idim) / uin_gas(n,i-i0,j-j0,k-k0,1)
           vup = uin_gas(n,i   ,j   ,k   ,1+idim) / uin_gas(n,i   ,j   ,k   ,1)

           ! HLLE wavespeeds: advection -+ diffusion speed, capped at vmax/sqrt(3)
           meanadv   = 0.5*(vdn+vup)
           meandiffv = 0.5*( lmax(n,i-i0,j-j0,k-k0,idim) + lmax(n,i,j,k,idim) )
           adn = min(meanadv-meandiffv, vdn-lmax(n,i-i0,j-j0,k-k0,idim))
           adn = max(adn,-vclamp)
           aup = max(meanadv+meandiffv, vup+lmax(n,i,j,k,idim))
           aup = min(aup,vclamp)
           lminus = min(adn,0.)
           lplus  = max(aup,0.)

           iFlx(n, i, j, k, iP0:iP1, idim)=&
                cmp_cr_face( fdn, fup, udn, uup, lminus, lplus)*dtdx
        end do
     end do
     end do
     end do
  end do

end subroutine cmp_cr_faces

!************************************************************************
PURE SUBROUTINE cr_limit_flux(Ecr, Fcr, vmax)

! Rescale a superluminal CR flux back to |F| <= vmax*Ecr (M1 closure)
! or vmax*Ecr/sqrt(3) (P1 closure).
!------------------------------------------------------------------------
  real(dp),intent(in)::Ecr, vmax
  real(dp),dimension(ndim),intent(inout)::Fcr
  real(dp)::fred
!------------------------------------------------------------------------
  fred = sqrt(sum(Fcr**2))/(vmax*Ecr)
  if(cr_isotropic_pressure) fred = fred*sqrt(3d0)
  if(fred>1.0) Fcr = Fcr/fred
END SUBROUTINE cr_limit_flux

!************************************************************************
PURE SUBROUTINE rotatevec(sint, cost, sinp, cosp, v1, v2, v3)
  !  Rotate vector v by t=theta and p=phi
  !  i.e. rotate to the local coordinate system from theta, phi.
  !  Hence the x-component of the result is the component of v parallel
  !  with the theta,phi vector.
  !------------------------------------------------------------------------
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
