! Gather routines related to DICE initial conditions

subroutine dice_reset_uold(ilevel)
  use amr_commons
  use hydro_commons
  use dice_commons
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine sets array uold to zero before calling
  ! the hydro scheme. uold is set to zero in virtual boundaries as well.
  ! Used for DICE initial conditions.
  !--------------------------------------------------------------------------
  integer::i,ivar,ind,icpu,iskip

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,112)ilevel

  ! Set uold to uold for myid cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do ivar=1,nvar_all
        do i=1,active(ilevel)%ngrid
           uold(active(ilevel)%igrid(i)+iskip,ivar)=0D0
        end do
     end do
  end do

  ! Set uold to 0 for virtual boundary cells
  do icpu=1,ncpu
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do ivar=1,nvar_all
        do i=1,reception(icpu,ilevel)%ngrid
           uold(reception(icpu,ilevel)%igrid(i)+iskip,ivar)=0D0
        end do
     end do
  end do
  end do

112 format('   Entering dice_reset_uold for level ',i2)

end subroutine dice_reset_uold


subroutine dice_init_uold(ilevel)
  use amr_commons
  use hydro_commons
  use dice_commons
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine sets array unew to its initial value uold before calling
  ! the hydro scheme. unew is set to zero in virtual boundaries.
  !--------------------------------------------------------------------------
  integer::i,ivar,ind,icpu,iskip,idim
  real(dp)::u,e
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,113)ilevel

  ! Set uold to namelist values for myid cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do ivar=nvar_all,1,-1
        do i=1,active(ilevel)%ngrid
           if(uold(active(ilevel)%igrid(i)+iskip,1).lt.IG_rho/scale_nH) then
              uold(active(ilevel)%igrid(i)+iskip,ivar)                      = 0D0
              if(ivar.eq.1) uold(active(ilevel)%igrid(i)+iskip,ivar)        = max(IG_rho/scale_nH,smallr)
              if(ivar.eq.ndim+2)then
                 uold(active(ilevel)%igrid(i)+iskip,ivar) = IG_T2/scale_T2/(gamma-1)*max(IG_rho/scale_nH,smallr)
              endif
              if(metal) then
                if(ivar.eq.imetal) uold(active(ilevel)%igrid(i)+iskip,ivar) = max(IG_rho/scale_nH,smallr)*IG_metal
              endif
           endif
        end do
     end do
  end do
  ! Set cell averaged kinetic energy
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,active(ilevel)%ngrid
        ! Initialisation of the refinement mask
        if((ic_mask_ivar.gt.0).and.(ivar_refine.gt.0).and.(ic_mask_ivar.le.nvar).and.(ivar_refine.le.nvar))then
     ! Switch to K/mu for ic_mask_ivar=ndim+2 case
           if(ic_mask_ivar.eq.ndim+2) then
        u = uold(active(ilevel)%igrid(i)+iskip,ic_mask_ivar)*scale_T2*(gamma-1)
           else
        u = uold(active(ilevel)%igrid(i)+iskip,ic_mask_ivar)
           endif
           if(ic_mask_ivar.gt.1)then
              u = u/uold(active(ilevel)%igrid(i)+iskip,1)
           endif
           if((u.ge.ic_mask_min).and.(u.le.ic_mask_max))then
              uold(active(ilevel)%igrid(i)+iskip,ivar_refine) = 1.0*uold(active(ilevel)%igrid(i)+iskip,1)
           endif
        endif
        e = 0d0
        do idim=1,ndim
           e = e+0.5*uold(active(ilevel)%igrid(i)+iskip,idim+1)**2/uold(active(ilevel)%igrid(i)+iskip,1)
        enddo
        uold(active(ilevel)%igrid(i)+iskip,ndim+2) = uold(active(ilevel)%igrid(i)+iskip,ndim+2)+e
     end do
  end do

#ifdef SOLVERmhd
  ! set constant magnetic field
  CALL dice_mag_constant(ilevel)
  ! toroidal field
  CALL dice_mag_compute(ilevel)
#endif

  ! Set uold to 0 for virtual boundary cells
  do icpu=1,ncpu
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do ivar=1,nvar_all
        do i=1,reception(icpu,ilevel)%ngrid
           uold(reception(icpu,ilevel)%igrid(i)+iskip,ivar)=0.0
        end do
     end do
  end do
  end do

113 format('   Entering dice_init_uold for level ',i2)

end subroutine dice_init_uold

subroutine condinit_dice(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  use dice_commons
  implicit none
  integer::ilevel
  !------------------------------------------------------------------
  ! This routine computes the initial density field at level ilevel using
  ! the CIC scheme from particles that are not entirely in
  ! level ilevel (boundary particles).
  ! Arrays flag1 and flag2 are used as temporary work space.
  !------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,idim,icpu,next_part
  integer::i,ig,ip,npart1,npart2
  real(dp)::dx

  integer,dimension(1:nvector),save::ind_grid,ind_cell
  integer,dimension(1:nvector),save::ind_part,ind_grid_part
  real(dp),dimension(1:nvector,1:ndim),save::x0

!$omp threadprivate(ind_grid,ind_cell,ind_part,ind_grid_part,x0)

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,114)ilevel

  ! Mesh spacing in that level
  dx=0.5D0**ilevel

  ! Loop over cpus
  do icpu=1,ncpu
     ! Loop over grids
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0

        ! Count gas particles
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if(ic_mask_ptype.eq.-1)then
                 if(idp(ipart).eq.1)then
                    npart2=npart2+1
                 endif
              else
                 npart2=npart2+1
              endif
              ipart=next_part  ! Go to next particle
           end do
        endif

        ! Gather gas particles
        if(npart2>0)then
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)

           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if(ic_mask_ptype.eq.-1)then
                 if(idp(ipart).eq.1)then
                    if(ig==0)then
                       ig=1
                       ind_grid(ig)=igrid
                    end if
                    ip=ip+1
                    ind_part(ip)=ipart
                    ind_grid_part(ip)=ig
                 endif
              else
                 if(ig==0)then
                    ig=1
                    ind_grid(ig)=igrid
                 end if
                 ip=ip+1
                 ind_part(ip)=ipart
                 ind_grid_part(ip)=ig
              endif
              if(ip==nvector)then
                 ! Lower left corner of 3x3x3 grid-cube
                 do idim=1,ndim
                    do i=1,ig
                       x0(i,idim)=xg(ind_grid(i),idim)-3.0D0*dx
                    end do
                 end do
                 do i=1,ig
                    ind_cell(i)=father(ind_grid(i))
                 end do
#if NDIM==3
                 if(amr_struct) then
                    call init_gas_ngp(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
                 else
#endif
                    call init_gas_cic(ind_cell,ind_part,ind_grid_part,x0,ig,ip,ilevel)
#if NDIM==3
                 endif
#endif
                 ip=0
                 ig=0
              end if
              ipart=next_part  ! Go to next particle
           end do
           ! End loop over particles

        end if

        igrid=next(igrid)   ! Go to next grid
     end do
     ! End loop over grids

     if(ip>0)then
        ! Lower left corner of 3x3x3 grid-cube
        do idim=1,ndim
           do i=1,ig
              x0(i,idim)=xg(ind_grid(i),idim)-3.0D0*dx
           end do
        end do
        do i=1,ig
           ind_cell(i)=father(ind_grid(i))
        end do
#if NDIM==3
        if(amr_struct) then
           call init_gas_ngp(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
        else
#endif
           call init_gas_cic(ind_cell,ind_part,ind_grid_part,x0,ig,ip,ilevel)
#if NDIM==3
        endif
#endif
     end if

  end do

114 format('   Entering condinit_dice for level ',I2)

end subroutine condinit_dice
!==================================================================================
!==================================================================================
!==================================================================================

subroutine init_gas_cic(ind_cell,ind_part,ind_grid_part,x0,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use dice_commons
  use cooling_module
  implicit none
  integer::ng,np,ilevel
  integer ,dimension(1:nvector)::ind_cell,ind_grid_part,ind_part
  real(dp),dimension(1:nvector,1:ndim)::x0
  !------------------------------------------------------------------
  ! This routine computes the initial density field at level ilevel using
  ! the CIC scheme. Only cells that are in level ilevel
  ! are updated by the input particle list.
  !------------------------------------------------------------------
  logical::error
  integer::j,ind,idim,nx_loc
  real(dp)::dx,dx_loc,scale
  ! Grid-based arrays
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  ! Particle-based arrays
  logical ,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector,1:ndim),save::xx,dd,dg
  integer ,dimension(1:nvector,1:ndim),save::ig,id,igg,igd,icg,icd
  real(dp),dimension(1:nvector,1:twotondim),save::vol
  integer ,dimension(1:nvector,1:twotondim),save::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:nvector),save::ethermal
  real(dp),dimension(1:nvector),save::vol_loc
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v

!$omp threadprivate(nbors_father_cells,ok,xx,dd,dg,ig,id,igg,igd,icg,icd,vol,igrid,icell,indp,kg,ethermal,vol_loc)


  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc(1:nvector)=dx_loc**ndim

  ! Gather neighboring father cells (should be present anytime !)
  call get3cubefather(ind_cell,nbors_father_cells,ng,ilevel)

  ! Rescale particle position at level ilevel
  do idim=1,ndim
     do j=1,np
        xx(j,idim)=xp(ind_part(j),idim)/scale+skip_loc(idim)
        xx(j,idim)=xx(j,idim)-x0(ind_grid_part(j),idim)
        xx(j,idim)=xx(j,idim)/dx
     end do
  end do

  ! Check for illegal moves
  error=.false.
  do idim=1,ndim
     do j=1,np
        if(xx(j,idim)<0.5D0.or.xx(j,idim)>5.5D0)error=.true.
     end do
  end do
  if(error)then
     write(*,*)'problem in cic'
     do idim=1,ndim
        do j=1,np
           if(xx(j,idim)<0.5D0.or.xx(j,idim)>5.5D0)then
              write(*,*)xx(j,1:ndim)
           endif
        end do
     end do
     stop
  end if

  ! CIC at level ilevel (dd: right cloud boundary; dg: left cloud boundary)
  do idim=1,ndim
     do j=1,np
        dd(j,idim)=xx(j,idim)+0.5D0
        id(j,idim)=dd(j,idim)
        dd(j,idim)=dd(j,idim)-id(j,idim)
        dg(j,idim)=1.0D0-dd(j,idim)
        ig(j,idim)=id(j,idim)-1
     end do
  end do

  ! Compute cloud volumes
#if NDIM==1
  do j=1,np
     vol(j,1)=dg(j,1)
     vol(j,2)=dd(j,1)
  end do
#endif
#if NDIM==2
  do j=1,np
     vol(j,1)=dg(j,1)*dg(j,2)
     vol(j,2)=dd(j,1)*dg(j,2)
     vol(j,3)=dg(j,1)*dd(j,2)
     vol(j,4)=dd(j,1)*dd(j,2)
  end do
#endif
#if NDIM==3
  do j=1,np
     vol(j,1)=dg(j,1)*dg(j,2)*dg(j,3)
     vol(j,2)=dd(j,1)*dg(j,2)*dg(j,3)
     vol(j,3)=dg(j,1)*dd(j,2)*dg(j,3)
     vol(j,4)=dd(j,1)*dd(j,2)*dg(j,3)
     vol(j,5)=dg(j,1)*dg(j,2)*dd(j,3)
     vol(j,6)=dd(j,1)*dg(j,2)*dd(j,3)
     vol(j,7)=dg(j,1)*dd(j,2)*dd(j,3)
     vol(j,8)=dd(j,1)*dd(j,2)*dd(j,3)
  end do
#endif

  ! Compute parent grids
  do idim=1,ndim
     do j=1,np
        igg(j,idim)=ig(j,idim)/2
        igd(j,idim)=id(j,idim)/2
     end do
  end do
#if NDIM==1
  do j=1,np
     kg(j,1)=1+igg(j,1)
     kg(j,2)=1+igd(j,1)
  end do
#endif
#if NDIM==2
  do j=1,np
     kg(j,1)=1+igg(j,1)+3*igg(j,2)
     kg(j,2)=1+igd(j,1)+3*igg(j,2)
     kg(j,3)=1+igg(j,1)+3*igd(j,2)
     kg(j,4)=1+igd(j,1)+3*igd(j,2)
  end do
#endif
#if NDIM==3
  do j=1,np
     kg(j,1)=1+igg(j,1)+3*igg(j,2)+9*igg(j,3)
     kg(j,2)=1+igd(j,1)+3*igg(j,2)+9*igg(j,3)
     kg(j,3)=1+igg(j,1)+3*igd(j,2)+9*igg(j,3)
     kg(j,4)=1+igd(j,1)+3*igd(j,2)+9*igg(j,3)
     kg(j,5)=1+igg(j,1)+3*igg(j,2)+9*igd(j,3)
     kg(j,6)=1+igd(j,1)+3*igg(j,2)+9*igd(j,3)
     kg(j,7)=1+igg(j,1)+3*igd(j,2)+9*igd(j,3)
     kg(j,8)=1+igd(j,1)+3*igd(j,2)+9*igd(j,3)
  end do
#endif
  do ind=1,twotondim
     do j=1,np
        if(nbors_father_cells(ind_grid_part(j),kg(j,ind)).gt.0) then
           igrid(j,ind)=son(nbors_father_cells(ind_grid_part(j),kg(j,ind)))
        else
           igrid(j,ind)=0
        endif
     end do
  end do

  ! Compute parent cell position
  do idim=1,ndim
     do j=1,np
        icg(j,idim)=ig(j,idim)-2*igg(j,idim)
        icd(j,idim)=id(j,idim)-2*igd(j,idim)
     end do
  end do
#if NDIM==1
  do j=1,np
     icell(j,1)=1+icg(j,1)
     icell(j,2)=1+icd(j,1)
  end do
#endif
#if NDIM==2
  do j=1,np
     icell(j,1)=1+icg(j,1)+2*icg(j,2)
     icell(j,2)=1+icd(j,1)+2*icg(j,2)
     icell(j,3)=1+icg(j,1)+2*icd(j,2)
     icell(j,4)=1+icd(j,1)+2*icd(j,2)
  end do
#endif
#if NDIM==3
  do j=1,np
     icell(j,1)=1+icg(j,1)+2*icg(j,2)+4*icg(j,3)
     icell(j,2)=1+icd(j,1)+2*icg(j,2)+4*icg(j,3)
     icell(j,3)=1+icg(j,1)+2*icd(j,2)+4*icg(j,3)
     icell(j,4)=1+icd(j,1)+2*icd(j,2)+4*icg(j,3)
     icell(j,5)=1+icg(j,1)+2*icg(j,2)+4*icd(j,3)
     icell(j,6)=1+icd(j,1)+2*icg(j,2)+4*icd(j,3)
     icell(j,7)=1+icg(j,1)+2*icd(j,2)+4*icd(j,3)
     icell(j,8)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
  end do
#endif

  ! Update mass density and number density fields
  do ind=1,twotondim

     ! Check if particles are entirely in level ilevel
     do j=1,np
        ok(j)=igrid(j,ind)>0
     end do

     ! Compute parent cell adress
     do j=1,np
        if(ok(j))then
           indp(j,ind)=ncoarse+(icell(j,ind)-1)*ngridmax+igrid(j,ind)
        end if
     end do

     do j=1,np
        ok(j)=ok(j).and.(idp(ind_part(j)).eq.1)
     end do

     ! Update hydro variables
     do j=1,np
        if(ok(j)) then
           ! Specific kinetic energy of the gas particle
           ethermal(j)=up(ind_part(j))
           ! Update hydro variable in CIC cells
           uold(indp(j,ind),1)=uold(indp(j,ind),1)+mp(ind_part(j))*vol(j,ind)/vol_loc(j)
           do idim=1,ndim
              uold(indp(j,ind),idim+1)=uold(indp(j,ind),idim+1)+mp(ind_part(j))*vol(j,ind)/vol_loc(j)*vp(ind_part(j),idim)
           end do
           uold(indp(j,ind),ndim+2)=uold(indp(j,ind),ndim+2)+mp(ind_part(j))*vol(j,ind)/vol_loc(j)*ethermal(j)
           if(metal) then
             uold(indp(j,ind),imetal)=uold(indp(j,ind),imetal)+mp(ind_part(j))*vol(j,ind)/vol_loc(j)*zp(ind_part(j))
           endif
        endif
        ! Update passive scalar mask
        if(ic_mask_ptype.gt.-1) then
            uold(indp(j,ind),ivar_refine)=uold(indp(j,ind),ivar_refine)+mp(ind_part(j))*vol(j,ind)/vol_loc(j)*maskp(ind_part(j))
        endif
     end do
  end do

end subroutine init_gas_cic

#if NDIM==3
subroutine init_gas_ngp(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use random
  use dice_commons
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !------------------------------------------------------------------
  ! This routine computes the initial density field at level ilevel using
  ! the NGP scheme. Only cells that are in level ilevel
  ! are updated by the input particle list.
  !------------------------------------------------------------------
  integer::i,j,idim,nx_loc
  real(dp)::dx,dx_loc,scale
  ! Grid based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  real(dp),dimension(1:nvector),save::ethermal
  ! Particle based arrays
  logical ,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector),save::vol_loc
  real(dp),dimension(1:nvector,1:ndim),save::x
  integer ,dimension(1:nvector,1:ndim),save::id,igd,icd
  integer ,dimension(1:nvector),save::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc

!$omp threadprivate(x0,ind_cell,nbors_father_cells,ok,id,igd,icd,igrid,icell,indp,kg,ethermal,vol_loc,x)

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc(1:nvector)=dx_loc**ndim

  ! Lower left corner of 3x3x3 grid-cube
  do idim=1,ndim
     do i=1,ng
        x0(i,idim)=xg(ind_grid(i),idim)-3.0D0*dx
     end do
  end do

  ! Gather 27 neighboring father cells (should be present anytime !)
  do i=1,ng
     ind_cell(i)=father(ind_grid(i))
  end do
  call get3cubefather(ind_cell,nbors_father_cells,ng,ilevel)

  ! Rescale position at level ilevel
  do idim=1,ndim
     do j=1,np
        x(j,idim)=xp(ind_part(j),idim)/scale+skip_loc(idim)
        x(j,idim)=x(j,idim)-x0(ind_grid_part(j),idim)
        x(j,idim)=x(j,idim)/dx
     end do
  end do

  ! NGP at level ilevel
  do idim=1,ndim
     do j=1,np
        id(j,idim)=x(j,idim)
     end do
  end do

  ! Compute parent grids
  do idim=1,ndim
     do j=1,np
        igd(j,idim)=id(j,idim)/2
     end do
  end do

  do j=1,np
     kg(j)=1+igd(j,1)+3*igd(j,2)+9*igd(j,3)
  end do

  do j=1,np
     if(nbors_father_cells(ind_grid_part(j),kg(j)).gt.0) then
        igrid(j)=son(nbors_father_cells(ind_grid_part(j),kg(j)))
     else
        igrid(j)=0
     endif
  end do

  ! Compute parent cell position
  do idim=1,ndim
     do j=1,np
        icd(j,idim)=id(j,idim)-2*igd(j,idim)
     end do
  end do

  do j=1,np
     icell(j)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
  end do

  ! Check if particles are entirely in level ilevel
  do j=1,np
     ok(j)=igrid(j)>0
  end do

  ! Compute parent cell adresses
  do j=1,np
     if(ok(j))then
        indp(j)=ncoarse+(icell(j)-1)*ngridmax+igrid(j)
     endif
  end do

  do j=1,np
     ok(j)=ok(j).and.(idp(ind_part(j)).eq.1)
  end do

  ! Update hydro variables
  do j=1,np
     if(ok(j))then
        ethermal(j)=up(ind_part(j))
        ! Update density in NGP cell
        uold(indp(j),1)=uold(indp(j),1)+mp(ind_part(j))/vol_loc(j)
        ! Update velocity in NGP cell
        do idim=1,ndim
           uold(indp(j),idim+1)=uold(indp(j),idim+1)+mp(ind_part(j))/vol_loc(j)*vp(ind_part(j),idim)
        end do
        ! Update temperature in NGP cell
        uold(indp(j),ndim+2)=uold(indp(j),ndim+2)+mp(ind_part(j))/vol_loc(j)*ethermal(j)
        ! Update passive hydro variables in NGP cell
        if(metal) then
           uold(indp(j),imetal)=uold(indp(j),imetal)+mp(ind_part(j))/vol_loc(j)*zp(ind_part(j))
        endif
     endif
     ! Update passive scalar mask
     if(ic_mask_ptype.gt.-1) then
       uold(indp(j),ivar_refine)=uold(indp(j),ivar_refine)+mp(ind_part(j))/vol_loc(j)*maskp(ind_part(j))
     endif
  end do

end subroutine init_gas_ngp
#endif

#ifdef SOLVERmhd
subroutine dice_mag_constant(ilevel)
  ! constant background magnetic field
  use amr_commons
  use hydro_commons
  use dice_commons
  implicit none
  integer::i,ind,iskip,ilevel

  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,active(ilevel)%ngrid
        uold(active(ilevel)%igrid(i)+iskip,6:8)           = ic_mag_const
        uold(active(ilevel)%igrid(i)+iskip,nvar+1:nvar+3) = ic_mag_const
     enddo
  enddo
end subroutine dice_mag_constant

subroutine dice_mag_compute(ilevel)
  use amr_commons
  !use pm_commons
  use hydro_commons
  use dice_commons
  !use random
  implicit none
  integer::i,j,ilevel,icell
  logical::nogal
  real(dp)::Axdl,Axdr,Axul,Axur
  real(dp)::Aydl,Aydr,Ayul,Ayur
  real(dp)::Azdl,Azdr,Azul,Azur
  real(dp)::Bxl,Bxr,Byl,Byr,Bzl,Bzr
  real(dp),dimension(1:3)::pos,cell_center

  real(dp)::dx,dxhalf,dxmin,dxminhalf,scale,dx_loc,vol_loc
  integer::nx_loc,ind,ix,iy,iz,iskip,nfine
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc

  nogal=.true.
  do i=1,MAXGAL
    ! check if galaxy has a magnetic field
    if (ic_mag_scale_B(i) .NE. 0.0) nogal=.false.
  enddo
  if (nogal) return

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  dxhalf = 0.5D0*dx
  dxmin = 0.5D0**nlevelmax
  dxminhalf = 0.5D0*dxmin
  nfine = 2**(nlevelmax-ilevel)

  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  ! compute field
  do ind=1,twotondim
    iskip=ncoarse+(ind-1)*ngridmax
    do j=1,active(ilevel)%ngrid
      icell=active(ilevel)%igrid(j)
      cell_center = xg(icell,:)+xc(ind,:)-skip_loc(:)

      ! edge coordinates
      ! Ax
      pos = cell_center
      pos(1) = pos(1) - dxhalf + dxminhalf
      pos(2) = pos(2) - dxhalf
      pos(3) = pos(3) - dxhalf
      Axdl=0.0;Axdr=0.0;Axul=0.0;Axur=0.0
      do i=1,nfine
        CALL dice_mag_toroidal(pos,1,Axdl)
        pos(1) = pos(1) + dxmin
      enddo
      pos(2) = pos(2) + dx
      do i=1,nfine
        pos(1) = pos(1) - dxmin
        CALL dice_mag_toroidal(pos,1,Axdr)
      enddo
      pos(3) = pos(3) + dx
      do i=1,nfine
        CALL dice_mag_toroidal(pos,1,Axur)
        pos(1) = pos(1) + dxmin
      enddo
      pos(2) = pos(2) - dx
      do i=1,nfine
        pos(1) = pos(1) - dxmin
        CALL dice_mag_toroidal(pos,1,Axul)
      enddo
      ! Ay
      pos = cell_center
      pos(1) = pos(1) - dxhalf
      pos(2) = pos(2) - dxhalf + dxminhalf
      pos(3) = pos(3) - dxhalf
      Aydl=0.0;Aydr=0.0;Ayul=0.0;Ayur=0.0
      do i=1,nfine
        CALL dice_mag_toroidal(pos,2,Aydl)
        pos(2) = pos(2) + dxmin
      enddo
      pos(1) = pos(1) + dx
      do i=1,nfine
        pos(2) = pos(2) - dxmin
        CALL dice_mag_toroidal(pos,2,Aydr)
      enddo
      pos(3) = pos(3) + dx
      do i=1,nfine
        CALL dice_mag_toroidal(pos,2,Ayur)
        pos(2) = pos(2) + dxmin
      enddo
      pos(1) = pos(1) - dx
      do i=1,nfine
        pos(2) = pos(2) - dxmin
        CALL dice_mag_toroidal(pos,2,Ayul)
      enddo
      ! Az
      pos = cell_center
      pos(1) = pos(1) - dxhalf
      pos(2) = pos(2) - dxhalf
      pos(3) = pos(3) - dxhalf + dxminhalf
      Azdl=0.0;Azdr=0.0;Azul=0.0;Azur=0.0
      do i=1,nfine
        CALL dice_mag_toroidal(pos,3,Azdl)
        pos(3) = pos(3) + dxmin
      enddo
      pos(1) = pos(1) + dx
      do i=1,nfine
        pos(3) = pos(3) - dxmin
        CALL dice_mag_toroidal(pos,3,Azdr)
      enddo
      pos(2) = pos(2) + dx
      do i=1,nfine
        CALL dice_mag_toroidal(pos,3,Azur)
        pos(3) = pos(3) + dxmin
      enddo
      pos(1) = pos(1) - dx
      do i=1,nfine
        pos(3) = pos(3) - dxmin
        CALL dice_mag_toroidal(pos,3,Azul)
      enddo

      ! Bx left
      Bxl = ((Azul - Azdl) - (Ayul - Aydl))/dx / nfine
      uold(icell+iskip,6)      = uold(icell+iskip,6) + Bxl
      ! By left
      Byl = ((Axul - Axdl) - (Azdr - Azdl))/dx / nfine
      uold(icell+iskip,7)      = uold(icell+iskip,7) + Byl
      ! Bz left
      Bzl = ((Aydr - Aydl) - (Axdr - Axdl))/dx / nfine
      uold(icell+iskip,8)      = uold(icell+iskip,8) + Bzl

      ! Bx right
      Bxr = ((Azur - Azdr) - (Ayur - Aydr))/dx / nfine
      uold(icell+iskip,nvar+1) = uold(icell+iskip,nvar+1) + Bxr
      ! By right
      Byr = ((Axur - Axdr) - (Azur - Azul))/dx / nfine
      uold(icell+iskip,nvar+2) = uold(icell+iskip,nvar+2) + Byr
      ! Bz right
      Bzr = ((Ayur - Ayul) - (Axur - Axul))/dx / nfine
      uold(icell+iskip,nvar+3) = uold(icell+iskip,nvar+3) + Bzr
    end do
  end do

end subroutine dice_mag_compute

subroutine dice_mag_toroidal(pos,dir,A)
  use dice_commons
  use amr_parameters, ONLY: boxlen
  implicit none
  real(dp)::r,h
  real(dp)::Ah,A
  integer::i,dir
  real(dp),dimension(1:3)::pos,gcenter,gaxis,xrel,grad
  real(dp)::sB,sR,sH

  do i=1,MAXGAL
    ! check if galaxy has a magnetic field
    if (ic_mag_scale_B(i) .EQ. 0.0) cycle

    gcenter(1) = 0.5d0 + ic_mag_center_x(i)/boxlen
    gcenter(2) = 0.5d0 + ic_mag_center_y(i)/boxlen
    gcenter(3) = 0.5d0 + ic_mag_center_z(i)/boxlen
    gaxis(1) = ic_mag_axis_x(i)
    gaxis(2) = ic_mag_axis_y(i)
    gaxis(3) = ic_mag_axis_z(i)
    sR = ic_mag_scale_R(i) / boxlen
    sH = ic_mag_scale_H(i) / boxlen
    sB = ic_mag_scale_B(i)

    ! coordinates in galaxy frame
    xrel = pos - gcenter
    h = DOT_PRODUCT(xrel,gaxis)
    grad = xrel - h*gaxis
    r = NORM2(grad)

    Ah = sB * sR * exp(-r/sR) * exp(-ABS(h)/sH)

    ! vector in cartesian frame
    A = A + Ah*gaxis(dir)
  end do
end subroutine
#endif



! Routines related to dumping the gas on the grid

#if NDIM==3
subroutine multipole_from_current_level(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ilevel
  !------------------------------------------------------------------
  ! This routine computes the density field at level ilevel using
  ! the CIC scheme from particles that are not entirely in
  ! level ilevel (boundary particles).
  ! Arrays flag1 and flag2 are used as temporary work space.
  ! Used for DICE initial conditions.
  !------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,idim,icpu,ind,iskip,ibound
  integer::j,ig,ip,npart1,npart2,next_part
  real(dp)::dx

  integer,dimension(1:nvector),save::ind_grid,ind_cell
  integer,dimension(1:nvector),save::ind_part,ind_grid_part
  real(dp),dimension(1:nvector,1:ndim),save::x0
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer ::nx_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  integer ::ix,iy,iz
  real(kind=8)::dx_loc,scale,vol_loc
  real(dp),dimension(1:3)::skip_loc
  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do
  !!!!!!!!!!!!!!!!!!!!!!!!!!!

  if(verbose)write(*,111)ilevel
  ! Mesh spacing in that level
  dx=0.5D0**ilevel


  ! Initialize unew field to zero
  do icpu=1,ncpu
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do idim=1,ndim+1
           do j=1,reception(icpu,ilevel)%ngrid
              unew(reception(icpu,ilevel)%igrid(j)+iskip,idim)=0.0D0
           end do
        end do
     end do
  end do
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do idim=1,ndim+1
        do j=1,active(ilevel)%ngrid
           unew(active(ilevel)%igrid(j)+iskip,idim)=0.0D0
        end do
     end do
  end do
  ! Reset unew in physical boundaries
  do ibound=1,nboundary
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do idim=1,ndim+1
           do j=1,boundary(ibound,ilevel)%ngrid
              unew(boundary(ibound,ilevel)%igrid(j)+iskip,idim)=0.0
           end do
        end do
     end do
  end do

  ! Loop over cpus
  do icpu=1,ncpu
     ! Loop over grids
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0

        ! Count gas particles
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if(idp(ipart).eq.1)then
                 npart2=npart2+1
              endif
              ipart=next_part  ! Go to next particle
           end do
        endif

        if(npart2>0)then
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              ! Select only gas particles
              if(idp(ipart).eq.1)then
                 if(ig==0)then
                    ig=1
                    ind_grid(ig)=igrid
                 end if
                 ip=ip+1
                 ind_part(ip)=ipart
                 ind_grid_part(ip)=ig
              endif
              if(ip==nvector)then
                 ! Lower left corner of 3x3x3 grid-cube
                 do idim=1,ndim
                    do j=1,ig
                       x0(j,idim)=xg(ind_grid(j),idim)-3.0D0*dx
                    end do
                 end do
                 do j=1,ig
                    ind_cell(j)=father(ind_grid(j))
                 end do
                 call ngp_amr_gas(ind_cell,ind_part,ind_grid_part,x0,ig,ip,ilevel)
                 ip=0
                 ig=0
              end if
              ipart=next_part  ! Go to next particle
           end do
           ! End loop over particles
        end if
        igrid=next(igrid)   ! Go to next grid
     end do
     ! End loop over grids

     if(ip>0)then
        ! Lower left corner of 3x3x3 grid-cube
        do idim=1,ndim
           do j=1,ig
              x0(j,idim)=xg(ind_grid(j),idim)-3.0D0*dx
           end do
        end do
        do j=1,ig
           ind_cell(j)=father(ind_grid(j))
        end do
        call ngp_amr_gas(ind_cell,ind_part,ind_grid_part,x0,ig,ip,ilevel)
     end if

  end do
  ! End loop over cpus

  ! Update boundaries
  do idim=1,ndim+1
     call make_virtual_reverse_dp(unew(1,idim),ilevel)
     call make_virtual_fine_dp(unew(1,idim),ilevel)
  end do

  ! Check for over-refinement
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do j=1,active(ilevel)%ngrid
        if(unew(active(ilevel)%igrid(j)+iskip,1)==0d0) then
           unew(active(ilevel)%igrid(j)+iskip,1)=smallr*vol_loc
           do idim=1,ndim
              unew(active(ilevel)%igrid(j)+iskip,idim+1)=(xg(active(ilevel)%igrid(j),idim)+xc(ind,idim)-skip_loc(idim))*scale &
                 & *unew(active(ilevel)%igrid(j)+iskip,1)
           end do
        endif
     end do
  end do

  do idim=1,ndim+1
     call make_virtual_fine_dp(unew(1,idim),ilevel)
  end do

  111 format('   Entering multipole_from_current_level for level',i2)

end subroutine multipole_from_current_level

subroutine ngp_amr_gas(ind_cell,ind_part,ind_grid_part,x0,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ng,np,ilevel
  integer ,dimension(1:nvector)::ind_cell,ind_grid_part,ind_part
  real(dp),dimension(1:nvector,1:ndim)::x0
  !------------------------------------------------------------------
  ! This routine computes the density field at level ilevel using
  ! the CIC scheme. Only cells that are in level ilevel
  ! are updated by the input particle list.
  !------------------------------------------------------------------
  integer::j,idim,nx_loc
  real(dp)::dx,dx_loc,scale
  ! Grid-based arrays
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  ! Particle-based arrays
  logical ,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector,1:ndim),save::x
  integer ,dimension(1:nvector,1:ndim),save::id,igd,icd
  integer ,dimension(1:nvector),save::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:nvector),save::vol_loc


  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc(1:nvector)=dx_loc**ndim

  call get3cubefather(ind_cell,nbors_father_cells,ng,ilevel)

  ! Rescale position at level ilevel
  do idim=1,ndim
     do j=1,np
        x(j,idim)=xp(ind_part(j),idim)/scale+skip_loc(idim)
        x(j,idim)=x(j,idim)-x0(ind_grid_part(j),idim)
        x(j,idim)=x(j,idim)/dx
     end do
  end do

  ! NGP at level ilevel
  do idim=1,ndim
     do j=1,np
        id(j,idim)=x(j,idim)
     end do
  end do

   ! Compute parent grids
  do idim=1,ndim
     do j=1,np
        igd(j,idim)=id(j,idim)/2
     end do
  end do
  do j=1,np
     kg(j)=1+igd(j,1)+3*igd(j,2)+9*igd(j,3)
     igrid(j)=son(nbors_father_cells(ind_grid_part(j),kg(j)))
  end do

  ! Check if particles are entirely in level ilevel
  ok(1:np)=.true.
  do j=1,np
     ok(j)=ok(j).and.igrid(j)>0
  end do

  ! Compute parent cell position
  do idim=1,ndim
     do j=1,np
        if(ok(j)) then
           icd(j,idim)=id(j,idim)-2*igd(j,idim)
        endif
     end do
  end do

  do j=1,np
     if(ok(j)) then
        icell(j)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
     endif
  end do

  ! Compute parent cell adresses
  do j=1,np
     if(ok(j))then
        indp(j)=ncoarse+(icell(j)-1)*ngridmax+igrid(j)
     else
        indp(j) = nbors_father_cells(ind_grid_part(j),kg(j))
     end if
  end do

  if(hydro)then
     do j=1,np
        unew(indp(j),1)=unew(indp(j),1)+mp(ind_part(j))
     end do
     do idim=1,ndim
        do j=1,np
           unew(indp(j),idim+1)=unew(indp(j),idim+1)+mp(ind_part(j))*xp(ind_part(j),idim)
        end do
     end do
  endif


end subroutine ngp_amr_gas
#endif
!###########################################################
!###########################################################
!###########################################################
!###########################################################

! Routines related to removing initial gas particles

subroutine kill_gas_part(ilevel)
  use pm_commons
  use amr_commons
  use mpi_mod
  implicit none
  integer::ilevel
  !--------------------------------------------------------
  ! This subroutine removes the gas particles
  ! initially present in the gadget1 DICE output
  !--------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,next_part
  integer::ig,ip,npart1,npart2,icpu,info
  integer,dimension(1:nvector)::ind_grid,ind_part,ind_grid_part
  logical,dimension(1:nvector)::ok=.true.
  integer::npart_all
  integer,dimension(1:ncpu)::npart_cpu,npart_cpu_all

  npart_cpu = 0
  npart_all = 0

  if(numbtot(1,ilevel)==0)return
  ! Gather gas particles.
  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0
        ! Count gas particles
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if(idp(ipart).eq.1)then
                 npart2=npart2+1
              endif
              ipart=next_part  ! Go to next particle
           end do
           npart_cpu(myid)=npart_cpu(myid)+npart2
        endif
        ! Gather gas particles
        if(npart2>0)then
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              ! Select only gas particles
              if(idp(ipart).eq.1)then
                 if(ig==0)then
                    ig=1
                    ind_grid(ig)=igrid
                 end if
                 ip=ip+1
                 ind_part(ip)=ipart
                 ind_grid_part(ip)=ig
              endif
              if(ip==nvector)then
                 call remove_list(ind_part,ind_grid_part,ok,ip)
                 call add_free_cond(ind_part,ok,ip)
                 ip=0
                 ig=0
              end if
              ipart=next_part  ! Go to next particle
           end do
           ! End loop over particles
        end if

        igrid=next(igrid)   ! Go to next grid
     end do

     ! End loop over grids
     if(ip>0)then
        call remove_list(ind_part,ind_grid_part,ok,ip)
        call add_free_cond(ind_part,ok,ip)
     end if
  end do

#ifndef WITHOUTMPI
  ! Give an array of number of gas on each cpu available to all cpus
  call MPI_ALLREDUCE(npart_cpu,npart_cpu_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#endif
  npart_all=sum(npart_cpu_all(1:ncpu))
  if(npart_all>0) then
     if(myid==1) then
        write(*,'(A50)')"__________________________________________________"
        write(*,'(A,I15)')' Gas particles deleted ->',npart_all
        write(*,'(A50)')"__________________________________________________"
     endif
  endif
npart_cpu(myid)=npart
#ifndef WITHOUTMPI
  ! Give an array of number of gas on each cpu available to all cpus
  call MPI_ALLREDUCE(npart_cpu,npart_cpu_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#endif
  npart_all=sum(npart_cpu_all(1:ncpu))
  if(npart_all>0) then
     if(myid==1) then
        write(*,'(A50)')"__________________________________________________"
        write(*,'(A,I15)')' Remaining particles ->',npart_all
        write(*,'(A50)')"__________________________________________________"
     endif
  endif

  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0
        ! Count gas particles
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if(idp(ipart).gt.0) idp(ipart)=idp(ipart)-1
              ipart=next_part  ! Go to next particle
           end do
           npart_cpu(myid)=npart_cpu(myid)+npart2
        endif
     end do
  end do


!---------------------------------------------
end subroutine
