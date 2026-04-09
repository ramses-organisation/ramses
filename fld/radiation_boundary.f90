!################################################################
!################################################################
!################################################################ 
!################################################################
subroutine make_boundary_diffusion(ilevel,igroup)
  use amr_commons
  use hydro_commons
  use radiation_parameters
  use units_commons
  implicit none
  ! -------------------------------------------------------------------
  ! This routine set up boundary conditions for fine levels.
  ! -------------------------------------------------------------------
  integer,intent(IN)::ilevel,igroup
  integer::ibound,boundary_dir,idim,inbor
  integer::i,ncache,ivar,igrid,ngrid,ind,ht
  integer::iskip,iskip_ref,nx_loc,ix,iy,iz,igrp
  
#if NDUST>0  
  integer::idust
#endif
  real(dp)::sum_dust
  
  integer,dimension(1:8)::ind_ref
  integer,dimension(1:nvector),save::ind_grid,ind_grid_ref
  integer,dimension(1:nvector),save::ind_cell,ind_cell_ref

  real(dp)::dx,dx_loc,scale
  real(dp)::rosseland_ana
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp),dimension(1:nvector,1:nvar+3),save::uu
  real(dp)::dd,t2,t2r,cal_Teg,usquare,emag,erad_loc,eps,ekin,Cv,rho

#if USE_FLD==1
  real(dp)::scale_nH,scale_T2,scale_t,scale_v,scale_d,scale_l
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
#endif

  if(.not. simple_boundary)return

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
     ! Zero flux
     if(boundary_type(ibound)== 1)ind_ref(1:8)=(/2,1,4,3,6,5,8,7/)
     if(boundary_type(ibound)== 2)ind_ref(1:8)=(/2,1,4,3,6,5,8,7/)
     if(boundary_type(ibound)== 3)ind_ref(1:8)=(/3,4,1,2,7,8,5,6/)
     if(boundary_type(ibound)== 4)ind_ref(1:8)=(/3,4,1,2,7,8,5,6/)
     if(boundary_type(ibound)== 5)ind_ref(1:8)=(/5,6,7,8,1,2,3,4/)
     if(boundary_type(ibound)== 6)ind_ref(1:8)=(/5,6,7,8,1,2,3,4/)
     ! Zero flux
     if(boundary_type(ibound)==11)ind_ref(1:8)=(/1,1,3,3,5,5,7,7/)
     if(boundary_type(ibound)==12)ind_ref(1:8)=(/2,2,4,4,6,6,8,8/)
     if(boundary_type(ibound)==13)ind_ref(1:8)=(/1,2,1,2,5,6,5,6/)
     if(boundary_type(ibound)==14)ind_ref(1:8)=(/3,4,3,4,7,8,7,8/)
     if(boundary_type(ibound)==15)ind_ref(1:8)=(/1,2,3,4,1,2,3,4/)
     if(boundary_type(ibound)==16)ind_ref(1:8)=(/5,6,7,8,5,6,7,8/)
     ! Imposed boundary
     if(boundary_type(ibound)==21)ind_ref(1:8)=(/1,1,3,3,5,5,7,7/)
     if(boundary_type(ibound)==22)ind_ref(1:8)=(/2,2,4,4,6,6,8,8/)
     if(boundary_type(ibound)==23)ind_ref(1:8)=(/1,2,1,2,5,6,5,6/)
     if(boundary_type(ibound)==24)ind_ref(1:8)=(/3,4,3,4,7,8,7,8/)
     if(boundary_type(ibound)==25)ind_ref(1:8)=(/1,2,3,4,1,2,3,4/)
     if(boundary_type(ibound)==26)ind_ref(1:8)=(/5,6,7,8,5,6,7,8/)

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

           ! Zero flux boundary conditions
           if((boundary_type(ibound)/10).ne.2)then

              ! Gather reference variables and  scatter to boundary region
              do i=1,ngrid
                 if(son(ind_cell(i)) == 0)then
                    uold(ind_cell(i),8+igroup)  = uold(ind_cell_ref(i),8+igroup)
                    unew(ind_cell(i),8+igroup)  = unew(ind_cell_ref(i),8+igroup)
                    unew(ind_cell(i),5)     = unew(ind_cell_ref(i),5)
                    unew(ind_cell(i),nvar+3)= unew(ind_cell_ref(i),5)
                    enew(ind_cell(i))       = unew(ind_cell_ref(i),8+igroup)
                    divu(ind_cell(i))       = divu(ind_cell_ref(i))
                    unew(ind_cell(i),2)     = unew(ind_cell_ref(i),2)
                 end if
              end do

              ! Imposed boundary conditions
           else

              ! Compute cell center in code units and rescale position from code units to user units
              do idim=1,ndim
                 do i=1,ngrid
                    if(son(ind_cell(i)) == 0)then
                       xx(i,idim)=(xg(ind_grid(i),idim)+xc(ind,idim)-skip_loc(idim))*scale
                    end if
                 end do
              end do
              
              call boundana(xx,uu,dx_loc,ibound,ngrid)
              
              ! Scatter variables
              do i=1,ngrid 
                 if(son(ind_cell(i)) == 0)then
                    dd=max(uu(i,1),smallr)
                    
                    usquare=0.0_dp
                    do idim=1,ndim
                       usquare=usquare+(uu(i,idim+1)/uu(i,1))**2
                    end do
                    ! Compute total magnetic energy
                    emag = 0.0_dp
                    do ivar=1,3
                       emag = emag + 0.125_dp*(uu(i,5+ivar) + uu(i,nvar+ivar))**2
                    end do
                    ! Compute total non-thermal+radiative energy
                    erad_loc=0.0_dp
                    do igrp=1,nener
                       erad_loc=erad_loc+uu(i,8+igrp)
                    enddo

                    rho   = uu(i,1)
                    ekin  = rho*usquare*0.5_dp
                    eps   = (uu(i,5)-ekin-emag-erad_loc)
                    sum_dust =0.0d0
#if NDUST>0
                    do idust = 1, ndust
                       sum_dust = sum_dust + uu(i,firstindex_ndust+idust)/rho
                    end do
#endif           
                    call temperature_eos((1.0d0-sum_dust)*rho,eps,t2,ht)                    
!                     t2    = Tr_floor ! comment this for radiative shock

                    unew(ind_cell(i),nvar+3) = t2
                    unew(ind_cell(i),5)      = t2
                    
                    uold(ind_cell(i),firstindex_er+igroup)= uu(i,firstindex_er+igroup)*scale_d*scale_v**2/(scale_E0)
                    unew(ind_cell(i),firstindex_er+igroup)= uold(ind_cell(i),firstindex_er+igroup)
                    enew(ind_cell(i)         )= uold(ind_cell(i),firstindex_er+igroup)
                    unew(ind_cell(i),2       )= 0.0_dp
                    
                    ! Compute Rosseland opacity
                    t2r = cal_Teg(unew(ind_cell(i),firstindex_er+igroup)*scale_E0,igroup)
                    divu(ind_cell(i))= rosseland_ana(dd*scale_d,t2,t2r,igroup,in_sink(ind_cell(i)))/scale_kappa
                    if(divu(ind_cell(i))*dx_loc .lt. min_optical_depth) divu(ind_cell(i))=min_optical_depth/dx_loc
                    
                 end if
              end do
           end if
              
        end do
        ! End loop over cells
           
     end do
     ! End loop over grids

  end do
  ! End loop over boundaries


111 format('   Entering make_boundary_diffusion for level ',I2)

end subroutine make_boundary_diffusion

!################################################################
!################################################################
!################################################################ 
!################################################################

subroutine make_boundary_diffusion_tot(ilevel)
  use amr_commons,only:boundary,son,ncoarse,nbor,xg
  use hydro_commons
  use radiation_parameters
  use const
  use units_commons
  implicit none
  ! -------------------------------------------------------------------
  ! This routine set up boundary conditions for fine levels.
  ! -------------------------------------------------------------------
  integer,intent(IN)::ilevel
  integer::ibound,boundary_dir,idim,inbor,igroup,ht
  integer::i,ncache,ivar,igrid,ngrid,ind
  integer::iskip,iskip_ref,gdim,nx_loc,ix,iy,iz,igrp,irad
  integer,dimension(1:8)::ind_ref
  integer,dimension(1:nvector),save::ind_grid,ind_grid_ref
  integer,dimension(1:nvector),save::ind_cell,ind_cell_ref

  real(dp)::dx,dx_loc,scale
  real(dp)::rosseland_ana,planck_ana
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp),dimension(1:nvector,1:nvar+3),save::uu
  real(dp),dimension(1:nvector)::cond,relax
  real(dp)::dd,t2,t2r,cal_Teg,usquare,emag,erad_loc,eps,ekin,Cv,rho

#if NDUST>0
  integer::idust
#endif
  real(dp)::sum_dust
  
#if USE_FLD==1
  real(dp)::scale_nH,scale_T2,scale_t,scale_v,scale_d,scale_l
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
#endif

  If(.not. simple_boundary)return

  ! Mesh size at level ilevel
  dx=half**ilevel

  ! Rescaling factors
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/zero,zero,zero/)
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
     if(ndim>0)xc(ind,1)=(dble(ix)-half)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-half)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-half)*dx
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
     ! Zero flux
     if(boundary_type(ibound)== 1)ind_ref(1:8)=(/2,1,4,3,6,5,8,7/)
     if(boundary_type(ibound)== 2)ind_ref(1:8)=(/2,1,4,3,6,5,8,7/)
     if(boundary_type(ibound)== 3)ind_ref(1:8)=(/3,4,1,2,7,8,5,6/)
     if(boundary_type(ibound)== 4)ind_ref(1:8)=(/3,4,1,2,7,8,5,6/)
     if(boundary_type(ibound)== 5)ind_ref(1:8)=(/5,6,7,8,1,2,3,4/)
     if(boundary_type(ibound)== 6)ind_ref(1:8)=(/5,6,7,8,1,2,3,4/)
     ! Zero flux
     if(boundary_type(ibound)==11)ind_ref(1:8)=(/1,1,3,3,5,5,7,7/)
     if(boundary_type(ibound)==12)ind_ref(1:8)=(/2,2,4,4,6,6,8,8/)
     if(boundary_type(ibound)==13)ind_ref(1:8)=(/1,2,1,2,5,6,5,6/)
     if(boundary_type(ibound)==14)ind_ref(1:8)=(/3,4,3,4,7,8,7,8/)
     if(boundary_type(ibound)==15)ind_ref(1:8)=(/1,2,3,4,1,2,3,4/)
     if(boundary_type(ibound)==16)ind_ref(1:8)=(/5,6,7,8,5,6,7,8/)
     ! Imposed boundary
     if(boundary_type(ibound)==21)ind_ref(1:8)=(/1,1,3,3,5,5,7,7/)
     if(boundary_type(ibound)==22)ind_ref(1:8)=(/2,2,4,4,6,6,8,8/)
     if(boundary_type(ibound)==23)ind_ref(1:8)=(/1,2,1,2,5,6,5,6/)
     if(boundary_type(ibound)==24)ind_ref(1:8)=(/3,4,3,4,7,8,7,8/)
     if(boundary_type(ibound)==25)ind_ref(1:8)=(/1,2,3,4,1,2,3,4/)
     if(boundary_type(ibound)==26)ind_ref(1:8)=(/5,6,7,8,5,6,7,8/)

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

           ! Zero flux boundary conditions
           if((boundary_type(ibound)/10).ne.2)then

              ! Gather reference variables and  scatter to boundary region
              do i=1,ngrid
                 if(son(ind_cell(i)) == 0)then
#if USE_FLD==1
                    do igroup=1,ngrp
                       kappaR_bicg(ind_cell(i),igroup)       = kappaR_bicg(ind_cell_ref(i),igroup)
                    enddo
#endif
                    do irad=1,nvar_trad
                       unew    (ind_cell(i),ind_trad(irad)) = unew    (ind_cell_ref(i),ind_trad(irad))
                       uold    (ind_cell(i),ind_trad(irad)) = uold    (ind_cell_ref(i),ind_trad(irad))
                    enddo
                    do irad = 1,nvar_bicg
                       if(bicg_to_cg) var_bicg(ind_cell(i),irad, 2) = var_bicg(ind_cell_ref(i),irad, 2)
                       var_bicg(ind_cell(i),irad, 5) = var_bicg(ind_cell_ref(i),irad, 5)
                       if(.not.bicg_to_cg) var_bicg(ind_cell(i),irad, 6) = var_bicg(ind_cell_ref(i),irad, 6)
                    enddo

                 end if
              end do

              ! Imposed boundary conditions
           else

              ! Compute cell center in code units and rescale position from code units to user units
              do idim=1,ndim
                 do i=1,ngrid
                    if(son(ind_cell(i)) == 0)then
                       xx(i,idim)=(xg(ind_grid(i),idim)+xc(ind,idim)-skip_loc(idim))*scale
                    end if
                 end do
              end do

              call boundana(xx,uu,dx_loc,ibound,ngrid)

              ! Scatter variables
              do i=1,ngrid 
                 if(son(ind_cell(i)) == 0)then
                    dd=max(uu(i,1),smallr)

                    usquare=zero
                    do idim=1,ndim
                       usquare=usquare+(uu(i,idim+1)/uu(i,1))**2
                    end do
                    ! Compute total magnetic energy
                    emag = zero
                    do ivar=1,3
                       emag = emag + 0.125_dp*(uu(i,5+ivar) + uu(i,nvar+ivar))**2
                    end do
                    ! Compute total non-thermal+radiative energy
                    erad_loc=zero
                    do igrp=1,nener
                       erad_loc=erad_loc+uu(i,8+igrp)
                    enddo

                    rho   = uu(i,1)
                    ekin  = rho*usquare*half
                    eps   = (uu(i,5)-ekin-emag-erad_loc)

                    sum_dust =0.0d0
#if NDUST>0
                    do idust = 1, ndust
                       sum_dust = sum_dust + uu(i,firstindex_ndust+idust)/rho
                    end do
#endif       
                    call temperature_eos((1.0_dp-sum_dust)*rho,eps,t2,ht)
                    
#if NGRP>0
                    uu(i,ind_trad(1)) = t2
!                    uu(i,ind_trad(1)) = Tr_floor ! comment this for radiative shock
#endif

#if USE_FLD==1
                    ! Compute Rosseland opacity
                    do igroup=1,ngrp
                       t2r = cal_Teg(uu(i,firstindex_er+igroup)*scale_d*scale_v**2,igroup)
                       kappaR_bicg(ind_cell(i),igroup)= rosseland_ana(dd*scale_d,uu(i,ind_trad(1)),t2r,igroup,in_sink(ind_cell(i)))/scale_kappa
                       if( kappaR_bicg(ind_cell(i),igroup)*dx_loc .lt. min_optical_depth)  kappaR_bicg(ind_cell(i),igroup)=min_optical_depth/dx_loc
                    enddo
#endif
                    do irad=1,nvar_trad
                       uold(ind_cell(i),ind_trad(irad)) = uu(i,ind_trad(irad)) / norm_trad(irad)
                       unew(ind_cell(i),ind_trad(irad)) = uold(ind_cell(i),ind_trad(irad))
                    enddo

                   do irad = 1,nvar_bicg
                      if(bicg_to_cg) var_bicg(ind_cell(i),irad, 2) = zero
                      var_bicg(ind_cell(i),irad, 5) = zero
                      if(.not.bicg_to_cg) var_bicg(ind_cell(i),irad, 6) = zero
                   enddo

                 end if
              end do

           end if

        end do
        ! End loop over cells

     end do
     ! End loop over grids

  end do
  ! End loop over boundaries


111 format('   Entering make_boundary_diffusion for level ',I2)

end subroutine make_boundary_diffusion_tot
