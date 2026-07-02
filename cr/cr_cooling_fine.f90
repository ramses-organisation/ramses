SUBROUTINE cr_cooling_fine(ilevel)
  ! Collisional cosmic-ray cooling (Coulomb + hadronic losses):
  ! Armillotta et al. 2021, 2024; Fitz Axen et al. 2024. Exponential decay
  ! of the CR energy and flux at the rate lambda_cr * n_H.
  !
  ! Called from crmom_step (after cr_set_uold) when cr_cooling=.true.
  use amr_commons
  use hydro_commons              ! gas uold, smallr
  use cr_parameters              ! Ecr_idx, ncr_groups, zeta_cr, cr_ne, cr_fneut, cr_c_fraction
  use cr_hydro_commons           ! cruold
  implicit none
  integer::ilevel
  !---------------------------------------------------------
  integer::i,ind,iskip,igrid,ngrid,ncache,nleaf,il,iGrp,icrE
  integer,dimension(1:nvector),save::ind_grid,ind_cell,ind_leaf
  real(dp)::dtcool,lambda_cr,decay
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(kind=8),dimension(1:nvector),save::nH

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,*)'Entering cr_cooling_fine: ilevel=',ilevel

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  dtcool=dtnew(ilevel)*scale_t
  lambda_cr=zeta_cr*(1d0+0.22d0*cr_ne+0.125d0*cr_fneut)   ! [cm^3 s^-1]

  ! Loop over myid grids by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do

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
        if(nleaf==0)cycle

        ! Gas density
        do i=1,nleaf
           nH(i)=MAX(uold(ind_leaf(i),1),smallr)
        end do

        ! Exponential decay of E_cr and F_cr per group (Fitz Axen et al. 2024)
        do iGrp=1,ncr_groups
           icrE=Ecr_idx(iGrp)
           do i=1,nleaf
              il=ind_leaf(i)
              if(cr_c_fraction==1.0d0)then
                 ! At cr_c_fraction==1, cr_c_fraction**2==1 and x*1d0==x
                 ! exactly, so the flux exponent equals the energy exponent:
                 ! compute the decay factor once and reuse for both.
                 decay=EXP(-lambda_cr*nH(i)*dtcool)
                 cruold(il,icrE)=cruold(il,icrE)*decay
                 cruold(il,icrE+1:icrE+ndim)=cruold(il,icrE+1:icrE+ndim)*decay
              else
                 cruold(il,icrE)=cruold(il,icrE) &
                      *EXP(-lambda_cr*nH(i)*dtcool)
                 cruold(il,icrE+1:icrE+ndim)=cruold(il,icrE+1:icrE+ndim) &
                      *EXP(-lambda_cr*nH(i)*dtcool*cr_c_fraction**2)
              endif
           end do
        end do
     end do
  end do

  if(verbose)write(*,*)'Exiting cr_cooling_fine: ilevel=',ilevel
END SUBROUTINE cr_cooling_fine
