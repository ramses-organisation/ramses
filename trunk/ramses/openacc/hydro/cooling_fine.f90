!######################################################################
! Subroutine that computes cooling for fine levels
!######################################################################
subroutine cooling_fine(ilevel)
  use amr_commons
  use hydro_commons
  use cooling_module
#ifdef RT
  use rt_parameters, only: rt_UV_hom,rt_isDiffuseUVsrc, nGroups, iGroups
  use rt_hydro_commons
  use rt_cooling_module, only: update_UVrates, n_U,iNpU,iFpU,rt_solve_cooling
  use UV_module
#endif
#ifdef ATON
  use radiation_commons, ONLY: Erad
#endif
!---------------------------------------------------------------------
! Declaration of variables
!---------------------------------------------------------------------
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  integer::ncache,ioct,igrid,ngrid,info
  integer,dimension(active(ilevel)%ngrid)::ind_grid
#ifndef DERIVED_DATATYPES
  ! Auxiliary variable to avoid explicit reference to the derived data type
  integer,dimension(active(ilevel)%ngrid)::active_ilevel_igrid
#endif
  integer::nleaf
  integer::icell,iskip,idim,nx_loc,ix,iy,iz,ivar,tot_cells
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(kind=8)::dtcool,nISM,nCOM,damp_factor,cooling_switch,t_blast
  real(dp)::polytropic_constant,Fpnew,Npnew
  integer,dimension(active(ilevel)%ngrid)::ind_cell
#ifdef RT
  real(dp)::scale_Np,scale_Fp
  logical,dimension(1:nvector),save::cooling_on=.true.
  real(dp),dimension(1:nvector,n_U),save::U,U_old
  real(dp),dimension(1:nvector,nGroups),save::Fp, Fp_precool
  real(dp),dimension(1:nvector,nGroups),save::dNpdt=0., dFpdt=0.
#endif
  real(dp),dimension(1:3)::skip_loc
  real(kind=8)::dx,dx_loc,scale,vol_loc
  !----------------------------------------------------------------------
  ! Variables to be allocated:
  integer,dimension(:),allocatable::ind_leaf
  real(dp),dimension(:,:),allocatable::uaux
  real(kind=8),dimension(:),allocatable::nH,T2,delta_T2,ekk,T2min,Zsolar,boost
  !----------------------------------------------------------------------
  ! Variables used only for debugging purposes:
  real(kind=8)::max_energy_difference
  real(dp),dimension(:),allocatable::energy_before,energy_after,energy_difference
  !----------------------------------------------------------------------

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Save an output file for debugging purposes:
  open(unit=12, file="aoutput_GPU.txt", action="write") !, status="replace")


!######################################################################
! Calculations of quantities that are invariants for the current level
!######################################################################
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

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
#ifdef RT
  call rt_units(scale_Np, scale_Fp)
#endif

  ! Typical ISM density in H/cc
  nISM = n_star; nCOM=0d0
  if(cosmo)then
     nCOM = del_star*omega_b*rhoc*(h0/100.)**2/aexp**3*X/mH
  endif
  nISM = MAX(nCOM,nISM)

  ! Polytropic constant for Jeans length related polytropic EOS
  if(jeans_ncells>0)then
     polytropic_constant=2d0*(boxlen*jeans_ncells*0.5d0**dble(nlevelmax)*scale_l/aexp)**2/ &
          & (twopi)*6.67e-8*scale_d*(scale_t/scale_l)**2
  endif


#ifndef DERIVED_DATATYPES
  active_ilevel_igrid = active(ilevel)%igrid
#endif

  ncache=active(ilevel)%ngrid
  !print*, "ncache = ", ncache
  tot_cells = ncache * twotondim

  ! Allocate array ind_leaf
  allocate(ind_leaf(1:tot_cells)) 


!######################################################################
! Loop over ALL cells: Gathering of the leaf cells
!######################################################################
nleaf=0
do ioct=1,ncache

#ifdef DERIVED_DATATYPES
        ind_grid(ioct)=active(ilevel)%igrid(ioct)
#else
        ind_grid(ioct)=active_ilevel_igrid(ioct)
#endif

  do icell=1,twotondim
     !print*, "icell =", icell, " out of ", twotondim
     iskip=ncoarse+(icell-1)*ngridmax
     ind_cell(ioct)=iskip+ind_grid(ioct)
   
     ! Gather leaf cells (creates the array ind_leaf
     ! containing the indices of the leaf cells for that level)
     if(son(ind_cell(ioct))==0)then
        nleaf=nleaf+1
        ind_leaf(nleaf)=ind_cell(ioct)
     end if
  end do
  if(nleaf.eq.0)cycle

end do
!######################################################################
! End of loop over ALL cells
!######################################################################
 
! print*, "Number of leaf cells:", nleaf, "out of a total number of cells of:", tot_cells

!---------------------------------------
! Allocate arrays:
  allocate(energy_before(1:nleaf)) 
  allocate(energy_after(1:nleaf)) 
  allocate(energy_difference(1:nleaf)) 
  
  allocate(nH(1:nleaf))
  allocate(T2(1:nleaf))
  allocate(delta_T2(1:nleaf))
  allocate(ekk(1:nleaf))
  allocate(T2min(1:nleaf))
  allocate(Zsolar(1:nleaf))
  allocate(boost(1:nleaf))
 
  allocate(uaux(1:nleaf,1:7))
!---------------------------------------

!-------------------------------------------------------------------------------
! Define quantities in the auxiliary variable uaux
! Note (Marco): I had to do this in order to avoid a runtime error
! when update device(uaux(:,idelay) and update device(uaux(:,imetal) were used.
do ivar=1,5
   do icell=1,nleaf
      uaux(icell,ivar) = uold(ind_leaf(icell),ivar)
   end do
end do

   do icell=1,nleaf
      uaux(icell,6) = uold(ind_leaf(icell),idelay)
   end do

   do icell=1,nleaf
      uaux(icell,7) = uold(ind_leaf(icell),imetal)
   end do
!-------------------------------------------------------------------------------


!-------------- for debugging purposes --------------
do icell=1,nleaf
   energy_before(icell) = uaux(icell,5)
end do
!----------------------------------------------------


!**********************************************************************
!   Start the data region
!$acc data copyin(uaux(:,1:7)) &
!$acc pcreate(nH, Zsolar, T2, ekk, boost, T2min, delta_T2)
!**********************************************************************
!######################################################################
! Loop over LEAF cells:
! All the calculations.
! Note: uold(:,ndim+2) and uold(:,idelay) are modified by these
!       calculations and therefore must be brought back to the host 
!######################################################################
!$acc parallel loop gang vector
do icell=1,nleaf
 
     !Compute density nH
        nH(icell)=MAX(uaux(icell,1),smallr)
     
     
     ! Compute metallicity in solar units   
     if(metal)then
           Zsolar(icell)=uaux(icell,7)/nH(icell)/0.02
     else
           Zsolar(icell)=z_ave
     endif

     
     ! Compute pressure
        T2(icell)=uaux(icell,ndim+2)
        ekk(icell)=0.0d0
   
 
     do idim=1,ndim
           ekk(icell)=ekk(icell)+0.5*uaux(icell,idim+1)**2/nH(icell)
     end do


     ! Compute T2=T/mu in Kelvin
        T2(icell)=(gamma-1.0)*(T2(icell)-ekk(icell))/nH(icell)*scale_T2
 

     ! Compute nH in H/cc
        nH(icell)=nH(icell)*scale_nH


     ! Compute radiation boost factor
     if(self_shielding)then
           boost(icell)=exp(-nH(icell)/0.01)
#ifdef ATON
     else if (aton) then
        do ioct=1,nleaf(ivs,icell)
           boost(aa-1+**nvector+ioct)=MAX(Erad(ind_leaf(ivs,ioct,icell))/J0simple(aexp), &
                &                   J0min/J0simple(aexp) )
        end do
#endif
     else
           boost(icell)=1.0
     endif


     !==========================================
     ! Compute temperature from polytrope EOS
     !==========================================
     if(jeans_ncells>0)then
           T2min(icell) = nH(icell)*polytropic_constant*scale_T2
     else
           T2min(icell) = T2_star*(nH(icell)/nISM)**(g_star-1.0)
     endif
     !==========================================
     ! You can put your own polytrope EOS here
     !==========================================


     if(cooling)then
        ! Compute thermal temperature by substracting polytrope
           T2(icell) = max(T2(icell)-T2min(icell),T2_min_fix)
     endif

     
     ! Compute cooling time step in second
     dtcool = dtnew(ilevel)*scale_t


#ifdef RT
     if(neq_chem) then
        ! Get gas thermal temperature
        do ioct=1,nleaf(ivs,icell)
           U(ioct,1) = T2(aa-1+**nvector+ioct)
        end do

        ! Get the ionization fractions
        do ivar=0,nIons-1
           do ioct=1,nleaf(ivs,icell)
              U(ioct,2+ivar) = uaux(icell,iIons+ivar)/uaux(icell,1)
           end do
        end do

        ! Get photon densities and flux magnitudes
        do ivar=1,nGroups
           do ioct=1,nleaf(ivs,icell)
              U(ioct,iNpU(ivar)) = scale_Np * rtuaux(icell,iGroups(ivar))
              U(ioct,iFpU(ivar)) = scale_Fp &
                   * sqrt(sum((rtuaux(icell,iGroups(ivar)+1:iGroups(ivar)+ndim))**2))
           enddo
           if(rt_smooth) then                           ! Smooth RT update
              do ioct=1,nleaf(ivs,icell) !Calc addition per sec to Np, Fp for current dt
                 Npnew = scale_Np * rtunew(ind_leaf(ivs,ioct,icell),iGroups(ivar))
                 Fpnew = scale_Fp &
                      * sqrt(sum((rtunew(ind_leaf(ivs,ioct,icell),iGroups(ivar)+1:iGroups(ivar)+ndim))**2))
                 dNpdt(ioct,ivar) = (Npnew - U(ioct,iNpU(ivar))) / dtcool
                 dFpdt(ioct,ivar) = (Fpnew - U(ioct,iFpU(ivar))) / dtcool ! Change in magnitude
                 ! Update flux vector to get the right direction
                 rtuold(ind_leaf(ivs,ioct,icell),iGroups(ivar)+1:iGroups(ivar)+ndim) = &
                      rtunew(ind_leaf(ivs,ioct,icell),iGroups(ivar)+1:iGroups(ivar)+ndim)
                 Fp_precool(ioct,ivar)=Fpnew           ! For update after solve_cooling
              end do
           else
              do ioct=1,nleaf(ivs,icell)
                 Fp_precool(ioct,ivar)=U(ioct,iFpU(ivar)) ! For update after solve_cooling
              end do
           end if
        end do

        if(cooling .and. delayed_cooling) then
           cooling_on(1:nleaf(ivs,icell))=.true.
           do ioct=1,nleaf(ivs,icell)
              if(uold(ind_leaf(ivs,ioct,icell),idelay)/uold(ind_leaf(ivs,ioct,icell),1) .gt. 1d-3) &
                   cooling_on(ioct)=.false.
           end do
        end if
        if(isothermal)cooling_on(1:nleaf(ivs,icell))=.false.
     endif
#endif

end do
!$acc end parallel loop
!######################################################################
    
!-------------------------------------------------------------------------
! Note: this function call must to be outside of the loop over the cells 
!-------------------------------------------------------------------------
! Compute net cooling delta_T2 at constant nH
if(cooling.and..not.neq_chem)then
   call solve_cooling(nH,T2,Zsolar,boost,dtcool,delta_T2,nleaf)
endif
!-------------------------------------------------------------------------

#ifdef RT
     if(neq_chem) then
        U_old=U
        call rt_solve_cooling(U, dNpdt, dFpdt, nH, cooling_on, Zsolar, dtcool, aexp, nleaf(ivs,icell))
        do ioct=1,nleaf(ivs,icell)
           delta_T2(aa-1+**nvector+ioct) = U(ioct,1) - T2(aa-1+**nvector+ioct)
        end do
     endif
#endif


!######################################################################
!$acc parallel loop gang vector
do icell=1,nleaf
     ! Convert density nH
        nH(icell) = nH(icell)/scale_nH


     ! Compute net energy sink
     if(cooling.or.neq_chem)then
           delta_T2(icell) = delta_T2(icell)*nH(icell)/scale_T2/(gamma-1.0)
        ! Turn off cooling in blast wave regions
        if(delayed_cooling)then
              cooling_switch = uaux(icell,6)/uaux(icell,1)
              if(cooling_switch < 1d-3) cycle
                 delta_T2(icell) = MAX(delta_T2(icell),real(0,kind=dp))
        endif
     endif


        ! Compute minimal total energy from polytrope
        T2min(icell) = T2min(icell)*nH(icell)/scale_T2/(gamma-1.0) + ekk(icell)
        ! Update total fluid energy
        T2(icell) = uaux(icell,ndim+2)

     
     if(cooling.or.neq_chem)then
           T2(icell) = T2(icell)+delta_T2(icell)
     endif


     if(isothermal)then
           uaux(icell,ndim+2) = T2min(icell)
     else
           uaux(icell,ndim+2) = max(T2(icell),T2min(icell))
     endif


     ! Update delayed cooling switch
     if(delayed_cooling)then
        t_blast=20d0*1d6*(365.*24.*3600.)
        damp_factor=exp(-dtcool/t_blast)
           uaux(icell,6)=uaux(icell,6)*damp_factor
     endif


#ifdef RT
     if(neq_chem) then
        ! Update ionization fraction
        do ivar=0,nIons-1
           do ioct=1,nleaf(ivs,icell)
              uold(ind_leaf(ivs,ioct,icell),iIons+ivar) = U(ioct,2+ivar)*nH(aa-1+**nvector+ioct)
           end do
        end do
     endif
     if(rt) then
        ! Update photon densities and flux magnitudes
        do ivar=1,nGroups
           do ioct=1,nleaf(ivs,icell)
              rtuold(ind_leaf(ivs,ioct,icell),iGroups(ivar)) = U(ioct,iNpU(ivar)) /scale_Np
              if(Fp_precool(ioct,ivar) .gt. 0.d0)then
                 rtuold(ind_leaf(ivs,ioct,icell),iGroups(ivar)+1:iGroups(ivar)+ndim) = U(ioct,iFpU(ivar))/Fp_precool(ioct,ivar) &
                      & *rtuold(ind_leaf(ivs,ioct,icell),iGroups(ivar)+1:iGroups(ivar)+ndim)
              endif
           enddo
        end do
     endif
#endif
  

end do
!$acc end parallel loop
!######################################################################
! End loop over LEAF cells
!######################################################################
!-------------------------------------------------------------------------
! Bring back to the host uaux entries corresponding to energy and idelay
!$acc update host(uaux(:,5:6))
!-------------------------------------------------------------------------
!**********************************************************************
!    End of data region
!$acc end data
!**********************************************************************


! Update entries corresponding to energy and idelay back in uold
   do icell=1,nleaf
      uold(ind_leaf(icell),5) = uaux(icell,5)
   end do

   do icell=1,nleaf
      uold(ind_leaf(icell),idelay) = uaux(icell,6)
   end do

!######################################################################
! Loop for debugging:
! Write to a file for debugging purposes (this is done by the CPU)
! Note: this can be commented out when no more useful.
!######################################################################
!print*, "Writing results to file..."
max_energy_difference = 0.0
do icell=1,nleaf
   energy_after(icell) = uaux(icell,5)
   energy_difference(icell) = abs(energy_after(icell) - energy_before(icell))

   if(energy_difference(icell).ge.max_energy_difference)then
      max_energy_difference=energy_difference(icell)
   end if
end do
! Write to file the max energy difference for the current level
write(12,*) max_energy_difference
!######################################################################
! End loop for debugging
!######################################################################


  if((cooling.and..not.neq_chem).and.ilevel==levelmin.and.cosmo)then
     if(myid==1)write(*,*)'Computing new cooling table'
     call set_table(dble(aexp))
  endif
#ifdef RT
  if(neq_chem.and.ilevel==levelmin) then
     if(cosmo)call update_rt_c
     if(cosmo .and. rt_UV_hom)call update_UVrates
     if(cosmo .and. rt_isDiffuseUVsrc)call update_UVsrc
     if(ilevel==levelmin) call output_rt_stats
  endif
#endif

!-------------------------------------------------------------------
! Deallocate arrays:
  deallocate(ind_leaf,energy_before,energy_after,energy_difference) 

  deallocate(nH)
  deallocate(T2)
  deallocate(delta_T2)
  deallocate(ekk)
  deallocate(T2min)
  deallocate(Zsolar)
  deallocate(boost)

  deallocate(uaux)
!-------------------------------------------------------------------

111 format('   Entering cooling_fine for level',i2)

end subroutine cooling_fine
!######################################################################
!######################################################################
!######################################################################
!######################################################################
