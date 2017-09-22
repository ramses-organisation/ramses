!################################################################
!################################################################
!################################################################
!################################################################
subroutine flag_formation_sites
#if NDIM==3
  use amr_commons
  use pm_commons
  use clfind_commons
  use hydro_commons, only:uold
  use hydro_parameters, only:smallr
  use pm_parameters, only:mass_halo_AGN,mass_clump_AGN
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  !=============================================================================
  ! This routine flags (flag2 = 1)  the cells where a sink is going to be formed
  !=============================================================================

  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp)::factG
  real(dp),dimension(1:nvector,1:3)::pos
  real(dp),dimension(1:ndim)::rrel
  integer,dimension(1:nvector)::cell_index,cell_levl,cc
  integer::j,jj,i,nx_loc,idim
  integer::global_peak_id,local_peak_id
  integer::merge_to,local_halo_id
  logical::ok
  real(dp)::dx,dx_min,dist2,scale
  real(dp),dimension(1:npeaks)::peakd
  integer,dimension(1:npeaks)::ind_sort
  logical,dimension(1:ndim)::period

  period(1)=(nx==1)
  period(2)=(ny==1)
  period(3)=(nz==1)

  ! Grid spacing and physical scales
  dx=0.5D0**nlevelmax
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_min=dx*scale

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Gravitational constant
  factG=1d0
  if(cosmo)factG=3d0/8d0/3.1415926*omega_m*aexp

  ! Loop over sinks and mark all clumps which are already occupied by a sink
  allocate(occupied(1:npeaks_max))
  occupied=0
  pos=0.0
  if(myid==1 .and. clinfo .and. nsink>0)write(*,*)'looping over ',nsink,' sinks and marking their clumps'

  ! Block clumps (halo done later) that contain a sink for formation
  if (smbh)then
     do j=1,nsink
        pos(1,1:3)=xsink(j,1:3)
        call cmp_cpumap(pos,cc,1)
        if (cc(1) .eq. myid)then
           call get_cell_index(cell_index,cell_levl,pos,nlevelmax,1)
           global_peak_id=flag2(cell_index(1))
           if(global_peak_id.NE.0)then
              call get_local_peak_id(global_peak_id,local_peak_id)
              occupied(local_peak_id)=1
              if(verbose)write(*,"('CPU # ',I5,' blocked clump # ',I6,' for sink production because of sink # ',I6)")myid,global_peak_id,idsink(j)
           endif
        end if
     end do
  end if

  ! Block clumps whose centres are closer than twice R_accretion from existing sinks
  do j=1,nsink
     do i=1,npeaks
        rrel=xsink(j,1:ndim)-peak_pos(i,1:ndim)
        do idim=1,ndim
           if (period(idim) .and. rrel(idim)>boxlen*0.5)rrel(idim)=rrel(idim)-boxlen
           if (period(idim) .and. rrel(idim)<boxlen*(-0.5))rrel(idim)=rrel(idim)+boxlen
        end do
        dist2=sum(rrel**2)
        if (dist2<(2.*ir_cloud*dx_min/aexp)**2)then
           occupied(i)=1
           if(clinfo)write(*,*)'CPU # ',myid,'blocked clump # ',i+ipeak_start(myid),' for sink production because of sink # ',idsink(j)
        end if
     end do
  end do
  ! This is required because new sink positions might bring new peaks into local memory in the smbh case
  call build_peak_communicator

#ifndef WITHOUTMPI
  call virtual_peak_int(occupied,'max')
  call boundary_peak_int(occupied)
#endif

  ! Block halos that contain a blocked clump
  if(smbh)then
     do i=1,npeaks
        if(occupied(i)==1)then
           merge_to=ind_halo(i)
           call get_local_peak_id(merge_to,local_halo_id)
           occupied(local_halo_id)=1
           if(verbose)write(*,"('CPU # ',I5,' blocked halo # ',I6,' for sink production because of clump # ',I6)")myid,merge_to,i+ipeak_start(myid)
        endif
     enddo
#ifndef WITHOUTMPI
     call virtual_peak_int(occupied,'max')
     call boundary_peak_int(occupied)
#endif
  endif

  !------------------------------------------------------------------------------
  ! Determine whether a peak patch is allowed to form a new sink.
  ! if a new sink has to be created, flag2 is set to the clump number at the peak position
  ! -> criteria to be chosen depend on the physics
  ! -> this routine can be patched
  !------------------------------------------------------------------------------
  pos=0.0
  flag2=0

  ! Sort clumps by peak density in ascending order
  do i=1,npeaks
     peakd(i)=max_dens(i)
     ind_sort(i)=i
  end do
  call quick_sort_dp(peakd,ind_sort,npeaks)

  ! Compute and combine various sink formation criteria
  do j=npeaks,1,-1
     jj=ind_sort(j)
     ok=.true.
     if (smbh)then
        ! Clump has to be a halo
        ok=ok.and.(ind_halo(jj).EQ.jj+ipeak_start(myid))
        ! Halo must have no existing sink
        ok=ok.and.occupied(jj)==0
        ! Halo mass has to be larger than some threshold
        ok=ok.and.halo_mass(jj)>mass_halo_AGN*2d33/(scale_d*scale_l**3.0)
        ! 4-cell ball mass has to be larger than some threshold
        ok=ok.and.clump_mass4(jj)>mass_clump_AGN*2d33/(scale_d*scale_l**3.0)
        ! Peak density has to be larger than star formation thresold
        ok=ok.and.max_dens(jj)>n_star/scale_nH
        ! Then create a sink at the peak position
        if (ok)then
           pos(1,1:3)=peak_pos(jj,1:3)
           call cmp_cpumap(pos,cc,1)
           if (cc(1) .eq. myid)then
              call get_cell_index(cell_index,cell_levl,pos,nlevelmax,1)
              ! Allow sink formation only inside zoom region
              if(ivar_refine>0)then
                 if(uold(cell_index(1),ivar_refine)/max(uold(cell_index(1),1),smallr)>var_cut_refine)then
                    flag2(cell_index(1))=jj
                    write(*,"('CPU # ',I5,' produces a new sink for clump # ',I6)")myid,jj+ipeak_start(myid)
                 end if
              else
                 flag2(cell_index(1))=jj
                 write(*,"('CPU # ',I5,' produces a new sink for clump # ',I6)")myid,jj+ipeak_start(myid)
              end if
            end if
         end if
     else
        ! Clump has to be peaky enough
        ok=ok.and.relevance(jj)>0.
        ! Clump has to contain at least one cell
        ok=ok.and.n_cells(jj)>0.
        ! Clmup must have no existing sink
        ok=ok.and.occupied(jj)==0
        ! Peak has to be dense enough
        ok=ok.and.max_dens(jj)>d_sink
        ! Clump has to be massive enough
        ok=ok.and.clump_mass4(jj)>mass_sink_seed*2d33/(scale_d*scale_l**3.0)
!!$        ! Clump has to be contracting
!!$        ok=ok.and.contracting(jj)
!!$        ! Clump has to be virialized
!!$        ok=ok.and.Icl_dd(jj)<0.
        ! Avoid formation of sinks from gas which is only compressed by thermal pressure rather than gravity.
        ok=ok.and.(kinetic_support(jj)<abs(grav_term(jj)))
        ! Avoid formation of crazy spins
        ok=ok.and.(kinetic_support(jj)<factG*mass_sink_seed*2d33/(scale_d*scale_l**3.0)/(ir_cloud*dx_min/aexp))
        ! Then create a sink at the peak position
        if (ok)then
           pos(1,1:3)=peak_pos(jj,1:3)
           call cmp_cpumap(pos,cc,1)
           if (cc(1) .eq. myid)then
              call get_cell_index(cell_index,cell_levl,pos,nlevelmax,1)
              flag2(cell_index(1))=jj
              write(*,*)'cpu ',myid,' produces a new sink for clump number ',jj+ipeak_start(myid)
           end if
        end if
     end if
  end do

  deallocate(occupied)

#endif
end subroutine flag_formation_sites
!################################################################
!################################################################
!################################################################
!################################################################
#if NDIM==3
subroutine compute_clump_properties_round2
  use amr_commons
#if NENER>0
  use hydro_commons, ONLY:uold,gamma,nvar,nener,inener,smallr
#else
  use hydro_commons, ONLY:uold,gamma,nvar,smallr
#endif
  use poisson_commons, ONLY:f
  use clfind_commons
  use pm_commons, ONLY:cont_speed
  use pm_parameters
#ifdef RT
  use rt_parameters, only: nGroups,ev_to_erg,iGroups,group_egy,c_cgs
  use rt_hydro_commons, only:rtuold
  use rt_cooling_module, only:kappaAbs,kappaSc
#endif

  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  !----------------------------------------------------------------------------
  ! This subroutine performs another loop over all particles and collects
  ! more information like binding energies, etc, that can not be created by
  ! just summing up cell properties.
  !----------------------------------------------------------------------------
  integer::ipart,ilevel,i,peak_nr,global_peak_id,j,ii,jj
  integer::grid,nx_loc,ix,iy,iz,ind,idim
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp)::d,vol,ekk,err,etot,p,T2
  real(dp)::dx,dx_loc,scale,vol_loc,abs_err,A1=0.,A2=0.,A3=0.
  real(dp),dimension(1:nlevelmax)::volume
  real(dp),dimension(1:3)::vd,xcell,xpeak,rrel,vrel,fgrav,skip_loc,frad
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:3,1:3)::eigenv,a
  real(dp),dimension(1:npeaks,1:3)::contractions
  logical,dimension(1:ndim)::period
  real(dp)::emag,pmag
#ifdef SOLVERmhd
  real(dp),dimension(1:3)::B
#endif
#ifdef RT
  integer::iNp, ig
  real(dp)::scale_Np,scale_Fp,scale_kappa
  real(dp)::kappa, c_code, ev_to_uu
  real(dp),dimension(1:nGroups,1:ndim)::Fp
  real(dp),dimension(1:nGroups)::Np2Ep_flux
#endif

#if NENER>0
  integer :: irad,nener_offset
  nener_offset = inener-1
#endif

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

#ifdef RT
  call rt_units(scale_Np, scale_Fp)
  ev_to_uu=ev_to_erg/(scale_d * scale_l**3 * scale_v**2)
  do ig=1,nGroups
     Np2Ep_flux(ig)=(scale_Fp * group_egy(ig) * ev_to_erg)/(scale_d * scale_v**3)
  end do
  scale_kappa=(scale_d*scale_l)**(-1.d0)
  c_code=c_cgs/scale_v
#endif

  period(1)=(nx==1)
  period(2)=(ny==1)
  period(3)=(nz==1)

  call surface_int

  ! Initialize arrays
  clump_size=0.d0; clump_mass4=0.d0
  grav_term=0.d0; rad_term=0.d0
  kinetic_support=0.d0; thermal_support=0.d0; magnetic_support=0.d0
  Icl=0.d0; Icl_d=0.d0; Icl_dd=0.d0
  Icl_3by3=0.d0;  Icl_d_3by3=0.d0
  contracting=.false.

  !------------------------------------------
  ! Compute volume of a cell in a given level
  !------------------------------------------
  do ilevel=1,nlevelmax
     ! Mesh spacing in that level
     dx=0.5D0**ilevel
     nx_loc=(icoarse_max-icoarse_min+1)
     scale=boxlen/dble(nx_loc)
     dx_loc=dx*scale
     vol_loc=dx_loc**ndim
     volume(ilevel)=vol_loc
  end do

  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)

  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     xc(ind,1)=(dble(ix)-0.5D0)
     xc(ind,2)=(dble(iy)-0.5D0)
     xc(ind,3)=(dble(iz)-0.5D0)
  end do

#ifndef WITHOUTMPI
  call boundary_peak_dp(max_dens)
  call boundary_peak_dp(clump_mass)
  do i=1,ndim
     call boundary_peak_dp(peak_pos(1,i))
  end do
#endif

  !---------------------------------------------------------------------------
  ! Loop over all test particles to collect information from the cells
  !---------------------------------------------------------------------------
  do ipart=1,ntest
     global_peak_id=flag2(icellp(ipart))
     if (global_peak_id /=0 ) then
        call get_local_peak_id(global_peak_id,peak_nr)

        ! Peak position
        do i=1,ndim
           xpeak(i)=peak_pos(peak_nr,i)
        end do

        ! Cell coordinates
        ind=(icellp(ipart)-ncoarse-1)/ngridmax+1 ! cell position
        grid=icellp(ipart)-ncoarse-(ind-1)*ngridmax ! grid index
        dx=0.5D0**levp(ipart)
        xcell(1:ndim)=(xg(grid,1:ndim)+xc(ind,1:ndim)*dx-skip_loc(1:ndim))*scale

        ! Periodic boundary conditions relative to peak position
        do idim=1,ndim
           if (period(idim) .and. (xcell(idim)-xpeak(idim))>boxlen*0.5)xcell(idim)=xcell(idim)-boxlen
           if (period(idim) .and. (xcell(idim)-xpeak(idim))<boxlen*(-0.5))xcell(idim)=xcell(idim)+boxlen
        end do

        ! Cell mass density
        d=max(uold(icellp(ipart),1),smallr)

        ! Cell total energy density
        etot=uold(icellp(ipart),ndim+2)

        ! Cell velocity
        do i=1,ndim
           vd(i)=uold(icellp(ipart),i+1)
        end do
        vd(1:3)=vd(1:3)/d

        ! Cell kinetic energy
        ekk=0.d0
        do i=1,3
           ekk=ekk+0.5*d*vd(i)**2
        end do

        ! Cell radiation flux
#ifdef RT
        do ig=1,nGroups
           iNp=iGroups(ig)
           Fp(ig,1:ndim) = rtuold(icellp(ipart),iNp+1:iNp+ndim)* Np2Ep_flux(ig)
        end do
#endif

        ! Cell magnetic field
#ifdef SOLVERmhd
        do i=1,ndim
           B(i)=0.5d0*(uold(icellp(ipart),5+i)+uold(icellp(ipart),nvar+i))
        end do
#endif

        ! Cell magnetic energy and magnetic pressure
        emag=0.d0; pmag=0.0d0
#ifdef SOLVERmhd
        do i=1,3
           emag=emag+0.5d0*B(i)**2
        end do
        pmag=emag
#endif

        ! Cell non-themal energy
        err=0.d0
#if NENER>0
        do irad=1,nener
           err=err+uold(icellp(ipart),nener_offset+irad)
        end do
#endif

        ! Cell thermal pressure and temperature
        p=(etot-ekk-emag-err)*(gamma-1)
        T2=p/d*scale_T2

        ! Add radiation pressure by trapped photons
        p=p+err/3.d0

        ! Cell volume
        vol=volume(levp(ipart))

        ! Properties of the cell relative to center of mass
        rrel(1:3)=xcell(1:3)-center_of_mass(peak_nr,1:3)
        vrel(1:3)=vd(1:3)-clump_velocity(peak_nr,1:3)

        ! Size relative to center of mass
        do i=1,ndim
           clump_size(peak_nr,i)=clump_size(peak_nr,i)+rrel(i)**2 * vol
        end do

        ! Properties for regions close to peak (4 cells away)
        if (((xpeak(1)-xcell(1))**2.+(xpeak(2)-xcell(2))**2.+(xpeak(3)-xcell(3))**2.) .LE. 16.*volume(nlevelmax)**(2./3.)/aexp**2)then
           clump_mass4(peak_nr)=clump_mass4(peak_nr)+d*vol
        endif

        ! Cell gravitational acceleration
        fgrav(1:3)=f(icellp(ipart),1:3)

        ! Cell radiation acceleration
        frad=0.d0
#ifdef RT
        do ig=1,nGroups
           kappa = kappaSc(ig)/scale_kappa
           frad(1:3)  = frad(1:3) +  Fp(ig,1:3) * kappa / c_code
        end do
#endif

        ! Virial analysis volume terms
        magnetic_support(peak_nr)  = magnetic_support(peak_nr) + 3*vol*pmag
        thermal_support (peak_nr)  = thermal_support (peak_nr) + 3*vol*p
        do i=1,3
           kinetic_support(peak_nr)= kinetic_support(peak_nr)  + vrel(i)**2         * vol*d
           grav_term(peak_nr)      = grav_term(peak_nr)        + fgrav(i) * rrel(i) * vol*d
           rad_term(peak_nr)       = rad_term(peak_nr)         + frad(i)  * rrel(i) * vol*d
        end do

        ! Time derivatives of the moment of inertia
        do i=1,3
           Icl_d(peak_nr)        = Icl_d(peak_nr)     + vrel(i)  * rrel(i) * vol*d
           Icl(peak_nr)          = Icl(peak_nr)       + rrel(i)  * rrel(i) * vol*d
           do j=1,3
              Icl_d_3by3(peak_nr,i,j)=  Icl_d_3by3(peak_nr,i,j)   + ( vrel(j) * rrel(i)  +  vrel(i) * rrel(j) )   * vol*d
              Icl_3by3(peak_nr,i,j)  =  Icl_3by3(peak_nr,i,j)     +   rrel(j) * rrel(i)                           * vol*d
           end do
        end do

     end if
  end do

  !---------------------------------------------------------------------------
  ! MPI communication to collect the results from the different CPU
  !---------------------------------------------------------------------------
#ifndef WITHOUTMPI
  call virtual_peak_dp(thermal_support,'sum')
  call virtual_peak_dp(kinetic_support,'sum')
  call virtual_peak_dp(magnetic_support,'sum')
  call virtual_peak_dp(clump_mass4,'sum')
  call virtual_peak_dp(Icl,'sum')
  call virtual_peak_dp(Icl_d,'sum')
  call virtual_peak_dp(grav_term,'sum')
  call virtual_peak_dp(rad_term,'sum')
  do i=1,ndim
     call virtual_peak_dp(clump_size(1,i),'sum')
     do j=1,ndim
        call virtual_peak_dp(Icl_3by3(1,i,j),'sum')
        call virtual_peak_dp(Icl_d_3by3(1,i,j),'sum')
     end do
  end do
#endif

  ! Second time derivative of I
  Icl_dd(1:npeaks)=2.*(grav_term(1:npeaks)+rad_term(1:npeaks)&
       -Psurf(1:npeaks)-MagPsurf(1:npeaks)+MagTsurf(1:npeaks)&
       +kinetic_support(1:npeaks)+thermal_support(1:npeaks)+magnetic_support(1:npeaks))

  do j=npeaks,1,-1
     if (relevance(j)>0..and.n_cells(j)>0)then
        contracting(j)=.true.
        if (n_cells(j)>1)then ! Only if more than one cell...
           ! Compute eigenvalues and eigenvectors of Icl_d_3by3
           a=Icl_3by3(j,1:3,1:3)
           abs_err=1.d-8*Icl(j)**2+1.d-40
           call jacobi(a,eigenv,abs_err)
           A1=a(1,1); A2=a(2,2); A3=a(3,3)

           ! Compute the contraction rates along the eigenvectors of Icl
           contractions(j,1:3)=0._dp
           do ii=1,3
              do jj=1,3
                 contractions(j,1)=contractions(j,1)+Icl_d_3by3(j,ii,jj)*eigenv(1,ii)*eigenv(1,jj)
                 contractions(j,2)=contractions(j,2)+Icl_d_3by3(j,ii,jj)*eigenv(2,ii)*eigenv(2,jj)
                 contractions(j,3)=contractions(j,3)+Icl_d_3by3(j,ii,jj)*eigenv(3,ii)*eigenv(3,jj)
              end do
           end do

           ! Check wether clump is contracting fast enough along all axis
           contracting(j)=contracting(j) .and. contractions(j,1)/(A1+tiny(0.d0)) < cont_speed
           contracting(j)=contracting(j) .and. contractions(j,2)/(A2+tiny(0.d0)) < cont_speed
           contracting(j)=contracting(j) .and. contractions(j,3)/(A3+tiny(0.d0)) < cont_speed
        end if
     endif
  end do

  ! Write to the log file some information that could be of interest for debugging etc.
  if(clinfo .and. (.not. smbh) .and. sink)then
     if (myid==1.and.npeaks>0)then
        write(*,'(135A)')'======================================================================================================================================'
        write(*,'(135A)')'Clump ID   Ncell   Mass     Mass4    Ekin     Eth    Emag     Egrav    Erad   '
        write(*,'(135A)')'======================================================================================================================================'
     end if
     do j=npeaks,1,-1
        if (relevance(j)>0..and.n_cells(j)>0)then
           write(*,'(I4,2X,I8,2x,3(L2,2X),(L2,6X),8(E9.2E2,3X))')j+ipeak_start(myid)&
                & ,n_cells(j)&
                & ,clump_mass(j),clump_mass4(j)&
                & ,kinetic_support(j),thermal_support(j),magnetic_support(j)&
                & ,grav_term(j),rad_term(j)
        end if
     end do
  end if

end subroutine compute_clump_properties_round2
#endif
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine trim_clumps
#if NDIM==3
  use amr_commons
  use clfind_commons
  use pm_commons, only:ir_cloud
  implicit none

  !---------------------------------------------------------------------------
  ! this routine trims the clumps down to the intersection of the clump with
  ! the accretion zone of the sink. Cells that are too far away from the peak
  ! are removed from the clump by setting flag2 to 0.
  !---------------------------------------------------------------------------
#ifndef WITHOUTMPI
  integer::ilevel
#endif
  integer::ipart,nx_loc,ind,idim
  real(dp)::dx,scale,dx_loc,r2
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  integer ::ix,iy,iz,grid,peak_nr,glob_peak_nr

  real(dp),dimension(1:3)::skip_loc,xcell,rrel
  real(dp),dimension(1:twotondim,1:3)::xc
  logical,dimension(1:ndim)::period

  period(1)=(nx==1)
  period(2)=(ny==1)
  period(3)=(nz==1)

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Mesh spacing in max level
  dx=0.5D0**nlevelmax
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale/aexp

  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     xc(ind,1)=(dble(ix)-0.5D0)
     xc(ind,2)=(dble(iy)-0.5D0)
     xc(ind,3)=(dble(iz)-0.5D0)
  end do

#ifndef WITHOUTMPI
  do idim=1,ndim
     call boundary_peak_dp(peak_pos(1,idim))
  end do
#endif

  ! Update flag 2
  do ipart=1,ntest
     glob_peak_nr=flag2(icellp(ipart))
     if (glob_peak_nr /=0 ) then
        call get_local_peak_id(glob_peak_nr,peak_nr)
        ind=(icellp(ipart)-ncoarse-1)/ngridmax+1 ! cell position
        grid=icellp(ipart)-ncoarse-(ind-1)*ngridmax ! grid index
        dx=0.5D0**levp(ipart)
        xcell(1:ndim)=(xg(grid,1:ndim)+xc(ind,1:ndim)*dx-skip_loc(1:ndim))*scale
        rrel=xcell(1:ndim)-peak_pos(peak_nr,1:ndim)
        do idim=1,ndim
           if (period(idim) .and. rrel(idim)>boxlen*0.5)rrel(idim)=rrel(idim)-boxlen
           if (period(idim) .and. rrel(idim)<boxlen*(-0.5))rrel(idim)=rrel(idim)+boxlen
        end do
        r2=sum(rrel(1:ndim)**2)
        if (r2 > (ir_cloud*dx_loc)**2.)then
           flag2(icellp(ipart))=0 ! Remove cell from clump
        end if
     end if
  end do

#ifndef WITHOUTMPI
  do ilevel=levelmin,nlevelmax
     call make_virtual_fine_int(flag2(1),ilevel)
  end do
#endif

#endif
end subroutine trim_clumps
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine jacobi(A,x,err2)
  use amr_commons, only:dp
  implicit none
  real(dp)::err2
  real(dp),dimension(3,3)::A,x

  !---------------------------------------------------------------------------
  ! Compute eigenvalues and eigenvectors using the jacobi-Method
  ! as for example described in Numerical Recipes.
  ! Returns eigenvalues as diagonal elements of A
  !---------------------------------------------------------------------------

  integer::n
  integer::i,j,k
  real(dp)::b2, bar,dummy
  real(dp)::beta, coeff, c, s, cs, sc

  n=3
  ! x is identity matrix initially
  x = 0.0
  do i=1,n
     x(i,i) = 1.0
  end do

  ! sum all squared off-diagonal elements
  b2 = 0.0
  do i=1,n
     do j=1,n
        if (i.ne.j) b2 = b2 + A(i,j)**2
     end do
  end do

  !return if already diagonal "enough"
  if (b2 <= err2) then
     return
  endif

  ! average for off-diagonal elements /2
  bar = 0.5*b2/9.

  do while (b2 > err2)
     do i=1,n-1
        do j=i+1,n
           if (A(j,i)**2 <= bar) cycle  ! do not touch small elements
           dummy=b2
           b2 = b2 - 2.0*A(j,i)**2
           bar = max(0.5*b2/9.,0.) !deal with weird optimized arithmetics...
           ! calculate coefficient c and s for Givens matrix
           beta = (A(j,j)-A(i,i))/(2.0*A(j,i))
           coeff = 0.5*beta*(1.0+beta**2)**(-0.5)
           s = (max(0.5+coeff,0.0))**0.5
           c = (max(0.5-coeff,0.0))**0.5
           ! update rows i and j
           do k=1,n
              cs =  c*A(i,k)+s*A(j,k)
              sc = -s*A(i,k)+c*A(j,k)
              A(i,k) = cs
              A(j,k) = sc
           end do
           ! find new matrix A_{k+1}
           do k=1,n
              cs =  c*A(k,i)+s*A(k,j)
              sc = -s*A(k,i)+c*A(k,j)
              A(k,i) = cs
              A(k,j) = sc
              cs =  c*x(k,i)+s*x(k,j)
              sc = -s*x(k,i)+c*x(k,j)
              x(k,i) = cs
              x(k,j) = sc
           end do
        end do
     end do
  end do
end subroutine jacobi
!################################################################
!################################################################
!################################################################
!################################################################
#if NDIM==3
subroutine surface_int
  use amr_commons
  use hydro_commons
  use pm_commons
  use clfind_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  !---------------------------------------------------------------
  ! Compute all the surface terms for virial analysis.
  !---------------------------------------------------------------

  integer::ipart,ip,ilevel,next_level
  integer,dimension(1:nvector)::ind_cell

  Psurf=0.d0; MagPsurf=0.d0; MagTsurf=0.d0

  ! loop cells that belong to peak patches
  ip=0
  do ipart=1,ntest
     ip=ip+1
     ilevel=levp(testp_sort(ipart))
     if (verbose.and.ilevel/=nlevelmax)print*,'not all particles in max level',ilevel
     next_level=0
     if(ipart<ntest)next_level=levp(testp_sort(ipart+1))
     ind_cell(ip)=icellp(testp_sort(ipart))
     if(ip==nvector .or. next_level /= ilevel)then
        call surface_int_np(ind_cell,ip,ilevel)
        ip=0
     endif
  end do
  if (ip>0)then
     call surface_int_np(ind_cell,ip,ilevel)
  endif

#ifndef WITHOUTMPI
  call virtual_peak_dp(Psurf,'sum')
  call virtual_peak_dp(MagPsurf,'sum')
  call virtual_peak_dp(MagTsurf,'sum')
#endif

end subroutine surface_int
#endif
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
#if NDIM==3
subroutine surface_int_np(ind_cell,np,ilevel)
  use amr_commons
#ifdef SOLVERmhd
  use clfind_commons, ONLY: center_of_mass,Psurf,MagPsurf,MagTsurf
#else
  use clfind_commons, ONLY: center_of_mass,Psurf
#endif
#if NENER>0
  use hydro_commons, ONLY: uold,gamma,nvar,nener,inener,smallr
#else
  use hydro_commons, ONLY: uold,gamma,nvar,smallr
#endif
  implicit none
  integer::np,ilevel
  integer,dimension(1:nvector)::ind_grid,ind_cell

  !------------------------------------------------------------
  ! This routine constructs all neighboring leaf cells that
  ! have a common cell surface at levels
  ! ilevel-1, ilevel, ilevel+1. Then, it computes the pressure
  ! pressure onto these surfaces and integrates over the surface
  ! of the clumps.
  !------------------------------------------------------------
  integer::j,ind,nx_loc,i2,j2,k2,ix,iy,iz,idim,jdim,i3,j3,k3
  real(dp)::dx,dx_loc,scale,vol_loc
  integer ,dimension(1:nvector)::cell_index,cell_levl,clump_nr,indv,neigh_cl
  integer ,dimension(1:nvector)::loc_clump_nr
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim)::xtest,r
  real(dp),dimension(1:nvector)::ekk_cell,emag_cell,P_cell,r_dot_n,err_cell
  real(dp)::P_neigh,err_neigh,emag_neigh,ekk_neigh
  real(dp),dimension(1:3)::skip_loc,n
  logical ,dimension(1:nvector)::ok
  logical,dimension(1:ndim)::period

#ifdef SOLVERmhd
  real(dp),dimension(1:nvector)::B_dot_n,B_dot_r,B2
  real(dp),dimension(1:nvector,1:3)::B
#endif

#if NENER>0
  integer::irad
  integer::nener_offset
  nener_offset = inener-1
#endif

  period(1)=(nx==1)
  period(2)=(ny==1)
  period(3)=(nz==1)

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**3

  ! Cells center position relative to grid center position
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     xc(ind,1)=(dble(ix)-0.5D0)*dx
     xc(ind,2)=(dble(iy)-0.5D0)*dx
     xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  emag_cell=0.d0; ekk_cell=0.d0; err_cell=0.d0

  ! some preliminary action...
  do j=1,np
     indv(j)=(ind_cell(j)-ncoarse-1)/ngridmax+1           ! cell position in grid
     ind_grid(j)=ind_cell(j)-ncoarse-(indv(j)-1)*ngridmax ! grid index
     clump_nr(j)=flag2(ind_cell(j))                       ! save clump number

     ! compute pressure in cell
     do jdim=1,ndim
        ekk_cell(j)=ekk_cell(j)+0.5*uold(ind_cell(j),jdim+1)**2
     end do
     ekk_cell(j)=ekk_cell(j)/max(uold(ind_cell(j),1),smallr)

#ifdef SOLVERmhd
     do idim=1,3
        emag_cell(j)=emag_cell(j)+0.125d0*(uold(ind_cell(j),idim+5)+uold(ind_cell(j),idim+nvar))**2
     end do
#endif
#if NENER>0
     do irad=1,nener
        err_cell(j)=err_cell(j)+uold(ind_cell(j),nener_offset+irad)
     end do
#endif
     P_cell(j)=(gamma-1.0)*(uold(ind_cell(j),ndim+2)-ekk_cell(j)-err_cell(j)-emag_cell(j))
  end do

  do j=1,np
     if (clump_nr(j) .ne. 0)call get_local_peak_id(clump_nr(j),loc_clump_nr(j))
  end do

  !================================
  ! generate neighbors at level ilevel and ilevel-1
  !================================
  ! Generate 3x3x3 neighboring cells
  do k2=0,2
     do j2=0,2
        do i2=0,2
           if((k2-1.)**2+(j2-1.)**2+(i2-1.)**2==1)then !check whether common face exists

              !construct outward facing normal vector
              n=0.d0
              if (k2==0)n(3)=-1.d0
              if (k2==2)n(3)=1.d0
              if (j2==0)n(2)=-1.d0
              if (j2==2)n(2)=1.d0
              if (i2==0)n(1)=-1.d0
              if (i2==2)n(1)=1.d0
              if (n(1)**2+n(2)**2+n(3)**2/=1)print*,'n has wrong lenght'


              r=0.d0
              do j=1,np
                 xtest(j,1)=(xg(ind_grid(j),1)+xc(indv(j),1)-skip_loc(1))*scale+(i2-1)*dx_loc
                 xtest(j,2)=(xg(ind_grid(j),2)+xc(indv(j),2)-skip_loc(2))*scale+(j2-1)*dx_loc
                 xtest(j,3)=(xg(ind_grid(j),3)+xc(indv(j),3)-skip_loc(3))*scale+(k2-1)*dx_loc

                 if (clump_nr(j)>0)then
                    r(j,1)=(xg(ind_grid(j),1)+xc(indv(j),1)-skip_loc(1))*scale+(i2-1)*dx_loc*0.5&
                         -center_of_mass(loc_clump_nr(j),1)

                    if (period(1) .and. r(j,1)>boxlen*0.5)r(j,1)=r(j,1)-boxlen
                    if (period(1) .and. r(j,1)<boxlen*(-0.5))r(j,1)=r(j,1)+boxlen

                    r(j,2)=(xg(ind_grid(j),2)+xc(indv(j),2)-skip_loc(2))*scale+(j2-1)*dx_loc*0.5&
                         -center_of_mass(loc_clump_nr(j),2)

                    if (period(2) .and. r(j,2)>boxlen*0.5)r(j,2)=r(j,2)-boxlen
                    if (period(2) .and. r(j,2)<boxlen*(-0.5))r(j,2)=r(j,2)+boxlen

                    r(j,3)=(xg(ind_grid(j),3)+xc(indv(j),3)-skip_loc(3))*scale+(k2-1)*dx_loc*0.5&
                         -center_of_mass(loc_clump_nr(j),3)

                    if (period(3) .and. r(j,3)>boxlen*0.5)r(j,3)=r(j,3)-boxlen
                    if (period(3) .and. r(j,3)<boxlen*(-0.5))r(j,3)=r(j,3)+boxlen

                 endif
              end do

              call get_cell_index(cell_index,cell_levl,xtest,ilevel,np)

              do j=1,np
                 ok(j)=(son(cell_index(j))==0)
              end do

              do j=1,np
                 neigh_cl(j)=flag2(cell_index(j))          ! number of clump the neighboring cell is in
                 ok(j)=ok(j).and. neigh_cl(j)/=clump_nr(j) ! neighboring cell is in another clump
                 ok(j)=ok(j).and. 0/=clump_nr(j)           ! clump number is not zero
              end do

              r_dot_n=0.d0
              do idim=1,3
                 do j=1,np
                    r_dot_n(j)=r_dot_n(j)+n(idim)*r(j,idim)
                 end do
              end do

#ifdef SOLVERmhd
              ! compute B field at the face by averaging the two cell_center values (all components)
              do j=1,np
                 B(j,1)=0.25d0*(uold(ind_cell(j),6)+uold(ind_cell(j),nvar+1) + &
                      uold(cell_index(j),6)+uold(cell_index(j),nvar+1) )
                 B(j,2)=0.25d0*(uold(ind_cell(j),7)+uold(ind_cell(j),nvar+2) + &
                      uold(cell_index(j),7)+uold(cell_index(j),nvar+2) )
                 B(j,3)=0.25d0*(uold(ind_cell(j),8)+uold(ind_cell(j),nvar+3) + &
                      uold(cell_index(j),8)+uold(cell_index(j),nvar+3) )
              end do

              B_dot_n=0.d0
              do idim=1,3
                 do j=1,np
                    B_dot_n(j)=B_dot_n(j)+B(j,idim)*n(idim)
                 end do
              end do

              B_dot_r=0.d0
              do idim=1,3
                 do j=1,np
                    B_dot_r(j)=B_dot_r(j)+B(j,idim)*r(j,idim)
                 end do
              end do

              B2=0.d0
              do idim=1,3
                 do j=1,np
                    B2(j)=B2(j)+B(j,idim)**2
                 end do
              end do
#endif

              do j=1,np
                 if (ok(j))then

                    ekk_neigh=0.d0
                    do jdim=1,ndim
                       ekk_neigh=ekk_neigh+0.5*uold(cell_index(j),jdim+1)**2
                    end do
                    ekk_neigh=ekk_neigh/max(uold(cell_index(j),1),smallr)

                    emag_neigh=0.d0
#ifdef SOLVERmhd
                    do jdim=1,ndim
                       emag_neigh=emag_neigh+0.125d0*(uold(cell_index(j),jdim+5)+uold(cell_index(j),jdim+nvar))**2
                    end do
#endif
                    err_neigh=0.d0
#if NENER>0
                    do irad=1,nener
                       err_neigh=err_neigh+uold(cell_index(j),nener_offset+irad)
                    end do
#endif
                    P_neigh=(gamma-1.0)*(uold(cell_index(j),ndim+2)-ekk_neigh-emag_neigh-err_neigh)

                    ! add to the actual terms for the virial analysis
                    Psurf(loc_clump_nr(j))    = Psurf(loc_clump_nr(j))    + 0.5d0*(P_neigh + P_cell(j)) * r_dot_n(j) * dx_loc**2
#ifdef SOLVERmhd
                    MagPsurf(loc_clump_nr(j)) = MagPsurf(loc_clump_nr(j)) + 0.5d0*B2(j)                 * r_dot_n(j) * dx_loc**2
                    MagTsurf(loc_clump_nr(j)) = MagTsurf(loc_clump_nr(j)) + B_dot_r(j)                  * B_dot_n(j) * dx_loc**2
#endif
                 endif
              end do
           endif
        end do
     end do
  end do


  !===================================
  ! generate neighbors at level ilevel+1
  !====================================
  if(ilevel<nlevelmax)then
     ! Generate 4x4x4 neighboring cells at level ilevel+1
     do k3=0,3
        do j3=0,3
           do i3=0,3
              if((k3-1.5)**2+(j3-1.5)**2+(i3-1.5)**2==2.75)then !check whether common face exists

                 n=0.
                 if (k3==0)n(3)=-1.d0
                 if (k3==3)n(3)=1.d0
                 if (j3==0)n(2)=-1.d0
                 if (j3==3)n(2)=1.d0
                 if (i3==0)n(1)=-1.d0
                 if (i3==3)n(1)=1.d0
                 if (n(1)**2+n(2)**2+n(3)**2/=1)print*,'n has wrong lenght'

                 r=0.d0
                 do j=1,np

                    xtest(j,1)=(xg(ind_grid(j),1)+xc(indv(j),1)-skip_loc(1))*scale+(i3-1.5)*dx_loc/2.0
                    xtest(j,2)=(xg(ind_grid(j),2)+xc(indv(j),2)-skip_loc(2))*scale+(j3-1.5)*dx_loc/2.0
                    xtest(j,3)=(xg(ind_grid(j),3)+xc(indv(j),3)-skip_loc(3))*scale+(k3-1.5)*dx_loc/2.0

                    if (clump_nr(j)>0)then
                       r(j,1)=(xg(ind_grid(j),1)+xc(indv(j),1)-skip_loc(1))*scale+(i3-1.5)*dx_loc/2.0*0.5&
                            -center_of_mass(loc_clump_nr(j),1)

                       if (period(1) .and. r(j,1)>boxlen*0.5)r(j,1)=r(j,1)-boxlen
                       if (period(1) .and. r(j,1)<boxlen*(-0.5))r(j,1)=r(j,1)+boxlen

                       r(j,2)=(xg(ind_grid(j),2)+xc(indv(j),2)-skip_loc(2))*scale+(j3-1.5)*dx_loc/2.0*0.5&
                            -center_of_mass(loc_clump_nr(j),2)

                       if (period(2) .and. r(j,2)>boxlen*0.5)r(j,2)=r(j,2)-boxlen
                       if (period(2) .and. r(j,2)<boxlen*(-0.5))r(j,2)=r(j,2)+boxlen

                       r(j,3)=(xg(ind_grid(j),3)+xc(indv(j),3)-skip_loc(3))*scale+(k3-1.5)*dx_loc/2.0*0.5&
                            -center_of_mass(loc_clump_nr(j),3)

                       if (period(3) .and. r(j,3)>boxlen*0.5)r(j,3)=r(j,3)-boxlen
                       if (period(3) .and. r(j,3)<boxlen*(-0.5))r(j,3)=r(j,3)+boxlen

                    endif
                 end do
                 call get_cell_index(cell_index,cell_levl,xtest,ilevel+1,np)

                 ok=.false.
                 do j=1,np
                    ! check wether neighbor is in a leaf cell at the right level
                    if(son(cell_index(j))==0.and.cell_levl(j)==(ilevel+1))ok(j)=.true.
                 end do

                 do j=1,np
                    neigh_cl(j)=flag2(cell_index(j))          ! number of clump the neighboring cell is in
                    ok(j)=ok(j).and. neigh_cl(j)/=clump_nr(j) ! neighboring cell is in another clump
                    ok(j)=ok(j).and. 0/=clump_nr(j)           ! clump number is not zero
                 end do

                 r_dot_n=0.d0
                 do idim=1,3
                    do j=1,np
                       r_dot_n(j)=r_dot_n(j)+n(idim)*r(j,idim)
                    end do
                 end do

#ifdef SOLVERmhd
                 ! compute B field at the face by averaging the two cell_center values (all components)
                 do j=1,np
                    B(j,1)=0.25d0*(uold(ind_cell(j),6)+uold(ind_cell(j),nvar+1) + &
                         uold(cell_index(j),6)+uold(cell_index(j),nvar+1) )
                    B(j,2)=0.25d0*(uold(ind_cell(j),7)+uold(ind_cell(j),nvar+2) + &
                         uold(cell_index(j),7)+uold(cell_index(j),nvar+2) )
                    B(j,3)=0.25d0*(uold(ind_cell(j),8)+uold(ind_cell(j),nvar+3) + &
                         uold(cell_index(j),8)+uold(cell_index(j),nvar+3) )
                 end do

                 B_dot_n=0.d0
                 do idim=1,3
                    do j=1,np
                       B_dot_n(j)=B_dot_n(j)+B(j,idim)*n(idim)
                    end do
                 end do

                 B_dot_r=0.d0
                 do idim=1,3
                    do j=1,np
                       B_dot_r(j)=B_dot_r(j)+B(j,idim)*r(j,idim)
                    end do
                 end do

                 B2=0.d0
                 do idim=1,3
                    do j=1,np
                       B2(j)=B2(j)+B(j,idim)**2
                    end do
                 end do
#endif
                 do j=1,np
                    if (ok(j))then

                       ekk_neigh=0.d0
                       do jdim=1,ndim
                          ekk_neigh=ekk_neigh+0.5*uold(cell_index(j),jdim+1)**2
                       end do
                       ekk_neigh=ekk_neigh/max(uold(cell_index(j),1),smallr)

                       emag_neigh=0.d0
#ifdef SOLVERmhd
                       do jdim=1,ndim
                          emag_neigh=emag_neigh+0.125d0*(uold(cell_index(j),jdim+5)+uold(cell_index(j),jdim+nvar))**2
                       end do
#endif
                       err_neigh=0.d0
#if NENER>0
                       do irad=1,nener
                          err_neigh=err_neigh+uold(cell_index(j),nener_offset+irad)
                       end do
#endif
                       P_neigh=(gamma-1.0)*(uold(cell_index(j),ndim+2)-ekk_neigh-emag_neigh-err_neigh)

                       ! add to the actual terms for the virial analysis
                       Psurf(loc_clump_nr(j))    = Psurf(loc_clump_nr(j))    + 0.5d0*(P_neigh + P_cell(j)) * r_dot_n(j) * 0.25d0 * dx_loc**2
#ifdef SOLVERmhd
                       MagPsurf(loc_clump_nr(j)) = MagPsurf(loc_clump_nr(j)) + 0.5d0*B2(j)                 * r_dot_n(j) * 0.25d0 * dx_loc**2
                       MagTsurf(loc_clump_nr(j)) = MagTsurf(loc_clump_nr(j)) + B_dot_r(j)                  * B_dot_n(j) * 0.25d0 * dx_loc**2
#endif
                    endif
                 end do
              endif
           end do
        end do
     end do
  endif

end subroutine surface_int_np
#endif

