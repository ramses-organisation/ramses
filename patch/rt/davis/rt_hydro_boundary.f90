!################################################################
!################################################################
!################################################################
!################################################################
subroutine rt_make_boundary_hydro(ilevel)
  use amr_commons
  use rt_hydro_commons
  use hydro_commons
  use poisson_commons                                                       !DAVIS
  use rt_cooling_module,only:a_r,iIR,is_kIR_T,kappaSc, T2_min_fix
  implicit none
  integer::ilevel
  ! -------------------------------------------------------------------
  ! This routine set up boundary conditions for fine levels.
  ! -------------------------------------------------------------------
  integer::ibound,boundary_dir,idim,inbor
  integer::i,ncache,ivar,igrid,ngrid,ind
  integer::iskip,iskip_ref,gdim,nx_loc,ix,iy,iz
  integer,dimension(1:8)::ind_ref,alt
  integer,dimension(1:nvector),save::ind_grid,ind_grid_ref
  integer,dimension(1:nvector),save::ind_cell,ind_cell_ref

  real(dp)::switch,dx,dx_loc,scale
  real(dp),dimension(1:3)::gs,skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp),dimension(1:nvector,1:nrtvar),save::uu

  integer::rtType !RT-- Type of rt variable =0,1,2,3
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2 &
       ,scale_Np,scale_Fp
  real(dp):: F_star, cE_up, Fy_limit, Fx_up, Fy_up, Fmag_up, fred_up, chi_up  &
       ,nbold, cPy_up, eps, Fy_down, cE_down, Fy_half, fred_down       &
       ,cPy_down, cPy_half, factor, Fmag, EIR, kir_sc, TR, tau, gas_rho
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  call rt_units(scale_Np, scale_Fp)

  if(.not. simple_boundary)return
  if(verbose)write(*,111)ilevel

  ! Mesh size at level ilevel
  dx=0.5D0**ilevel

  ! Rescaling factors
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  ! Loop over boundaries
  do ibound=1,nboundary
     ! Compute direction of reference neighbors
     boundary_dir=boundary_type(ibound)-10*(boundary_type(ibound)/10)
     if(boundary_dir==1)inbor=2
     if(boundary_dir==2)inbor=1
     if(boundary_dir==3)inbor=4
     if(boundary_dir==4)inbor=3
     if(boundary_dir==5)inbor=6
     if(boundary_dir==6)inbor=5

     ! Compute index of reference cells
     ! Reflexive boundary
     if(boundary_type(ibound)== 1)ind_ref(1:8)=(/2,1,4,3,6,5,8,7/)
     if(boundary_type(ibound)== 2)ind_ref(1:8)=(/2,1,4,3,6,5,8,7/)
     if(boundary_type(ibound)== 3)ind_ref(1:8)=(/3,4,1,2,7,8,5,6/)
     if(boundary_type(ibound)== 4)ind_ref(1:8)=(/3,4,1,2,7,8,5,6/)
     if(boundary_type(ibound)== 5)ind_ref(1:8)=(/5,6,7,8,1,2,3,4/)
     if(boundary_type(ibound)== 6)ind_ref(1:8)=(/5,6,7,8,1,2,3,4/)
     ! Free boundary
     if(boundary_type(ibound)==11)ind_ref(1:8)=(/1,1,3,3,5,5,7,7/)
     if(boundary_type(ibound)==12)ind_ref(1:8)=(/2,2,4,4,6,6,8,8/)
     if(boundary_type(ibound)==13)ind_ref(1:8)=(/1,2,1,2,5,6,5,6/)
     if(boundary_type(ibound)==14)ind_ref(1:8)=(/3,4,3,4,7,8,7,8/)
     if(boundary_type(ibound)==15)ind_ref(1:8)=(/1,2,3,4,1,2,3,4/)
     if(boundary_type(ibound)==16)ind_ref(1:8)=(/5,6,7,8,5,6,7,8/)
     ! Imposed boundary (used only for flag1)
     if(boundary_type(ibound)==21)ind_ref(1:8)=(/1,1,3,3,5,5,7,7/)
     if(boundary_type(ibound)==22)ind_ref(1:8)=(/2,2,4,4,6,6,8,8/)
     if(boundary_type(ibound)==23)ind_ref(1:8)=(/1,2,1,2,5,6,5,6/)
     if(boundary_type(ibound)==24)ind_ref(1:8)=(/3,4,3,4,7,8,7,8/)
     if(boundary_type(ibound)==25)ind_ref(1:8)=(/1,2,3,4,1,2,3,4/)
     if(boundary_type(ibound)==26)ind_ref(1:8)=(/5,6,7,8,5,6,7,8/)

     ! Velocity sign switch for reflexive boundary conditions
     gs=(/1,1,1/)
     if(boundary_type(ibound)==1.or.boundary_type(ibound)==2)gs(1)=-1
     if(boundary_type(ibound)==3.or.boundary_type(ibound)==4)gs(2)=-1
     if(boundary_type(ibound)==5.or.boundary_type(ibound)==6)gs(3)=-1

     ! Loop over grids by vector sweeps
     ncache=boundary(ibound,ilevel)%ngrid
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=boundary(ibound,ilevel)%igrid(igrid+i-1)
        end do

        ! Gather neighboring reference grid
        do i=1,ngrid
           ind_grid_ref(i)=son(nbor(ind_grid(i),inbor))
        end do

        ! Loop over cells
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do

           ! Gather neighboring reference cell
           iskip_ref=ncoarse+(ind_ref(ind)-1)*ngridmax
           do i=1,ngrid
              ind_cell_ref(i)=iskip_ref+ind_grid_ref(i)
           end do

           ! Wall and free boundary conditions
           if((boundary_type(ibound)/10).ne.2)then

              ! Gather reference hydro variables
              do ivar=1,nrtvar
                 do i=1,ngrid
                    uu(i,ivar)=rtuold(ind_cell_ref(i),ivar)
                 end do
              end do
              ! BEGIN DAVIS-----------------------------------------------
              if(boundary_dir .eq. 3) then ! Bottom boundary
                 do i=1,ngrid

!#define DAVIS_FZERO
#ifdef DAVIS_FZERO
                    ! Simple and flux works, but pressure not correct
                    F_star = 1.03d4/scale_Fp
                    factor=1d0
                    if(rt_flux_scheme .eq. 'hll') factor=sqrt(3d0)
                    !uu(i,1)=uu(i,1) + factor/rt_c * (2d0*F_star - uu(i,3))
                    !uu(i,2)=0d0
                    !uu(i,3)=0d0

                    gas_rho = uold(ind_cell_ref(i),1)
                    kIR_sc  = kappaSc(iIR)
                    if(is_kIR_T) then         !                   k_IR depends on T
                       EIR = group_egy(iGroupIR) * ev_to_erg * uu(i,1) *scale_Np
                       TR = max(T2_min_fix,(EIR*rt_c_cgs/c_cgs/a_r)**0.25)
                       kIR_sc  = kIR_sc  * (min(TR,820d0)/10d0)**2
                       if(gas_rho*scale_d .le. 7.22d-26) kIR_sc=0d0
                    endif
                    kIR_sc = kIR_sc*scale_d*scale_l ! Convert to code units

                    tau = gas_rho * kIR_sc * dx_loc
                    !uu(i,1) = uu(i,1) + factor/rt_c * (F_star*max(2d0,3d0*tau) - uu(i,3))
                    uu(i,1) = uu(i,1) + factor/rt_c * (2d0*F_star - uu(i,3))
                    uu(i,2) = 0d0
                    uu(i,3) = 0d0

#endif

#define DAVIS_FNONZERO
#ifdef DAVIS_FNONZERO
                    !! Fy_down=F_star, and adjust cE_down accordingly
                    F_star = 1.03d4/scale_Fp
                    cE_up = uu(i,1)*rt_c
                    Fx_up = uu(i,2)
                    Fy_up = uu(i,3)
                    Fy_down = F_star
                    cE_down = Fy_down - Fy_up + cE_up
                    uu(i,1) = cE_down/rt_c
                    uu(i,2) = 0d0
                    uu(i,3) = Fy_down

!                    ! Interpolate betweeen free-streaming and diffusion limits:
!                    gas_rho = uold(ind_cell_ref(i),1)
!                    kIR_sc  = rt_kIR_sc
!                    if(is_kIR_T) then         !                   k_IR depends on T
!                       EIR = group_egy(iGroupIR) * ev_to_erg * uu(i,1) *scale_Np
!                       TR = max(T2_min_fix,(EIR*rt_c_cgs/c_cgs/a_r)**0.25)
!                       kIR_sc  = kIR_sc  * (min(TR,820d0)/10d0)**2
!                       if(gas_rho*scale_d .le. 7.22d-26) kIR_sc=0d0
!                    endif
!                    kIR_sc = kIR_sc*scale_d*scale_l ! Convert to code units
!
!                    tau = gas_rho * kIR_sc * dx_loc
!                    cE_down = cE_up - Fy_up + Fy_down * ( 1d0/(tau+1d0) + 3d0*tau )
!                    uu(i,1) = cE_down/rt_c
!                    uu(i,2) = 0d0
!                    uu(i,3) = Fy_down

!                    ! Try to do as close to Fx_down = -Fx_up as possible:
!                    uu(i,2) = -Fx_up
!                    Fmag = sqrt(Fy_down**2+Fx_up**2)
!                    if( cE_down .lt. Fmag ) then
!                       uu(i,2) = -Fx_up/abs(Fx_up) * sqrt(cE_down**2 - Fy_down**2)
!                    endif
#endif

                 end do
              endif
              if(boundary_dir .eq. 4) then ! Top boundary = zero valued
                 do i=1,ngrid
                    uu(i,1) = 1d-50
                    uu(i,2) = 0d0
                    uu(i,3) = 0d0
                 end do
              endif
              ! END DAVIS-------------------------------------------------
              ! Scatter to boundary region
              do ivar=1,nrtvar
                 switch=1
                 ! ------------------------------------------------------------
                 ! need to switch photon flux in RT...
                 ! rtVar= ivar = 1,2,3,4,.. = Np, Fx, Fy, Fz, Np, Fx,..
                 rtType = mod(ivar-1, ndim+1)
                 !if(rtType .ne. 0) switch=gs(rtType) !0=Np, 1=Fx, 2=Fy, 3=Fz
                 if(rtType .ne. 0 .and. boundary_dir .ne. 3) &
                      switch=gs(rtType) !0=Np, 1=Fx, 2=Fy, 3=Fz
                 ! ------------------------------------------------------------
                 do i=1,ngrid
                    rtuold(ind_cell(i),ivar)=uu(i,ivar)*switch
                 end do

              end do

              ! Imposed boundary conditions
           else

              ! Compute cell center in code units
              do idim=1,ndim
                 do i=1,ngrid
                    xx(i,idim)=xg(ind_grid(i),idim)+xc(ind,idim)
                 end do
              end do

              ! Rescale position from code units to user units
              do idim=1,ndim
                 do i=1,ngrid
                    xx(i,idim)=(xx(i,idim)-skip_loc(idim))*scale
                 end do
              end do

              call rt_boundana(xx,uu,dx_loc,ibound,ngrid)

              ! Scatter variables
              do ivar=1,nrtvar
                 do i=1,ngrid
                    rtuold(ind_cell(i),ivar)=uu(i,ivar)
                 end do
              end do

           end if

        end do
        ! End loop over cells

     end do
     ! End loop over grids

  end do
  ! End loop over boundaries

111 format('   Entering rt_make_boundary_hydro for level ',I2)

end subroutine rt_make_boundary_hydro
