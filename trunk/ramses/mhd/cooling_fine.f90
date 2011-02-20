subroutine cooling_fine(ilevel)
  use amr_commons
  use hydro_commons
  use cooling_module
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !-------------------------------------------------------------------
  ! Compute cooling for fine levels
  !-------------------------------------------------------------------
  integer::ncache,i,igrid,ngrid,info,isink
  integer,dimension(1:nvector),save::ind_grid

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Compute sink accretion rates
  if(sink)call compute_accretion_rate(0)

  ! Operator splitting step for cooling source term
  ! by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     call coolfine1(ind_grid,ngrid,ilevel)
  end do

  if(cooling.and.ilevel==levelmin.and.cosmo)then
     if(myid==1)write(*,*)'Computing new cooling table'
     call set_table(dble(aexp))
  end if

111 format('   Entering cooling_fine for level',i2)

end subroutine cooling_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine coolfine1(ind_grid,ngrid,ilevel)
  use amr_commons
  use hydro_commons
  use cooling_module
  implicit none
  integer::ilevel,ngrid
  integer,dimension(1:nvector)::ind_grid
  !-------------------------------------------------------------------
  !-------------------------------------------------------------------
  integer::i,ind,iskip,idim,nleaf,neul=5
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(kind=8)::dtcool,nISM,nCOM
  integer,dimension(1:nvector),save::ind_cell,ind_leaf
  real(kind=8),dimension(1:nvector),save::nH,T2,delta_T2,ekk,emag
  real(kind=8),dimension(1:nvector),save::T2min,Zsolar,boost

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Typical ISM density in H/cc
  nISM = n_star; nCOM=0d0
  if(cosmo)then
     nCOM = del_star*omega_b*rhoc*(h0/100.)**2/aexp**3*X/mH
  endif
  nISM = MAX(nCOM,nISM)

  ! Loop over cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,ngrid
        ind_cell(i)=iskip+ind_grid(i)
     end do

     ! Gather leaf cells
     nleaf=0
     do i=1,ngrid
        if(son(ind_cell(i))==0)then
           nleaf=nleaf+1
           ind_leaf(nleaf)=ind_cell(i)
        end if
     end do

     ! Compute rho
     do i=1,nleaf
        nH(i)=MAX(uold(ind_leaf(i),1),smallr)
     end do
     
     ! Compute metallicity in solar units
     if(metal)then
        do i=1,nleaf
           Zsolar(i)=uold(ind_leaf(i),neul+4)/nH(i)/0.02
        end do
     else
        do i=1,nleaf
           Zsolar(i)=z_ave
        end do
     endif

     ! Compute pressure
     do i=1,nleaf
        T2(i)=uold(ind_leaf(i),neul)
     end do
     do i=1,nleaf
        ekk(i)=0.0d0
     end do
     do idim=1,3
        do i=1,nleaf
           ekk(i)=ekk(i)+0.5*uold(ind_leaf(i),idim+1)**2/nH(i)
        end do
     end do
     do i=1,nleaf
        emag(i)=0.0d0
     end do
     do idim=1,3
        do i=1,nleaf
           emag(i)=emag(i)+0.125d0*(uold(ind_leaf(i),idim+neul)+uold(ind_leaf(i),idim+nvar))**2
        end do
     end do
     do i=1,nleaf
        T2(i)=(gamma-1.0)*(T2(i)-ekk(i)-emag(i))
     end do

     ! Compute T2=T/mu in Kelvin
     do i=1,nleaf
        T2(i)=T2(i)/nH(i)*scale_T2
     end do

     ! Compute nH in H/cc
     do i=1,nleaf
        nH(i)=nH(i)*scale_nH
     end do

     ! Compute radiation boost factor
     if(self_shielding)then
        do i=1,nleaf
           boost(i)=exp(-nH(i)/0.01)
        end do
     else
        do i=1,nleaf
           boost(i)=1.0
        end do
     endif

     !==========================================
     ! Compute temperature from polytrope EOS
     !==========================================
     do i=1,nleaf
        T2min(i) = T2_star*(nH(i)/nISM)**(g_star-1.0)
     end do
     !==========================================
     ! You can put your own polytrope EOS here
     !==========================================

     ! Compute cooling time step in second
     dtcool = dtnew(ilevel)*scale_t

     ! Compute net cooling at constant nH
     if(cooling)then
        ! Compute "thermal" temperature by substracting polytrope
        do i=1,nleaf
           T2(i)=max(T2(i)-T2min(i),T2_min_fix)
        end do
        call solve_cooling(nH,T2,Zsolar,boost,dtcool,delta_T2,nleaf)
     endif

     ! Compute rho
     do i=1,nleaf
        nH(i) = nH(i)/scale_nH
     end do

     ! Compute net energy sink
     if(cooling)then
        do i=1,nleaf
           delta_T2(i) = delta_T2(i)*nH(i)/scale_T2/(gamma-1.0)
        end do
     endif

     ! Compute minimal total energy from polytrope
     do i=1,nleaf
        T2min(i) = T2min(i)*nH(i)/scale_T2/(gamma-1.0) + ekk(i) + emag(i)
     end do

     ! Update total fluid energy
     do i=1,nleaf
        T2(i) = uold(ind_leaf(i),neul)
     end do
     if(cooling)then
        do i=1,nleaf
           T2(i) = T2(i)+delta_T2(i)
        end do
     endif
     if(isothermal)then
        do i=1,nleaf
           uold(ind_leaf(i),neul) = T2min(i)
        end do
     else
        do i=1,nleaf
           uold(ind_leaf(i),neul) = max(T2(i),T2min(i))
        end do
     endif

  end do
  ! End loop over cells

end subroutine coolfine1



