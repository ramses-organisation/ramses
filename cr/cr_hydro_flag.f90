! CR refinement pass -- the cosmic-ray counterpart of rt_hydro_flag
! (rt/rt_hydro_flag.f90). Flags cells whose CR energy gradient exceeds
! err_grad_crmom, run as a SEPARATE pass after hydro_flag and OR-ing into the
! shared set-once flag1 (exactly as RT runs rt_hydro_flag separately). The grid/
! neighbour machinery mirrors rt_hydro_flag; the errcr gradient block and the
! floor_crmom constant are ported VERBATIM from the former #ifdef CRPHYS block in
! hydro/hydro_flag.f90, so the gas refinement routine stays CR-free and the
! refinement decision (hence the AMR grid) is byte-identical: flag1 is a 0/1
! set-once flag and A.or.B is order-independent across the two passes.
subroutine cr_hydro_flag(ilevel)
  use amr_commons
  use cr_parameters, only: ncr,ncrvars,iCRu,err_grad_crmom
  use cr_hydro_commons, only: cruold
  implicit none
  integer::ilevel
  ! -------------------------------------------------------------------
  ! Flag for refinement cells whose CR energy gradient at level ilevel
  ! exceeds err_grad_crmom. Separate pass; results OR into flag1.
  ! -------------------------------------------------------------------
  integer::i,j,ncache,nok,iskip
  integer::igrid,ind,idim,ngrid
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  integer,dimension(1:nvector,0:twondim),save::igridn
  integer,dimension(1:nvector,1:twondim),save::indn

  logical,dimension(1:nvector),save::ok

  real(dp),dimension(1:nvector,1:ncrvars),save::ucrg,ucrm,ucrd
  real(dp)::pcrg,pcrm,pcrd,errcr
  integer::icr,icrE
  ! Matches cral's floor_prad in the CR-gradient denominator (mhd/godunov_utils:381)
  real(dp),parameter::floor_crmom=1d-10

  if(ilevel==nlevelmax)return
  if(numbtot(1,ilevel)==0)return
  ! No CR refinement requested -> nothing to do (analogue of rt_err_grad_cn==-1).
  if(all(err_grad_crmom(1:ncr) < 0d0))return

  ! Loop over active grids
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector

     ! Gather nvector grids
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do

     ! Gather neighboring offsets
     call getnborgrids(ind_grid,igridn,ngrid)

     ! Loop over cells
     do ind=1,twotondim

        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do

        ! Initialize refinement to false
        do i=1,ngrid
           ok(i)=.false.
        end do

        ! Gather neighboring cells
        call getnborcells(igridn,ind,indn,ngrid)

        ! If a neighbor cell does not exist,
        ! replace it by its father cell
        do j=1,twondim
           do i=1,ngrid
              if(indn(i,j)==0)then
                 indn(i,j)=nbor(ind_grid(i),j)
              end if
           end do
        end do

        ! Loop over dimensions
        do idim=1,ndim
           ! CR-energy gradient refinement (err_grad_crmom). Ported verbatim
           ! from the former CR block in hydro/hydro_flag.f90: cral evaluates
           ! this inside hydro_refine on the embedded CR slots; the separated
           ! module reads the CR energy from cruold here.
           do icr=1,ncr
              if(err_grad_crmom(icr) >= 0.)then
                 icrE=iCRu+(ndim+1)*(icr-1)
                 do i=1,ngrid
                    ucrg(i,icrE)=cruold(indn(i,2*idim-1),icrE)
                    ucrm(i,icrE)=cruold(ind_cell(i     ),icrE)
                    ucrd(i,icrE)=cruold(indn(i,2*idim  ),icrE)
                 end do
                 do i=1,ngrid
                    pcrg=ucrg(i,icrE); pcrm=ucrm(i,icrE); pcrd=ucrd(i,icrE)
                    errcr=2.0d0*MAX( &
                         & ABS((pcrd-pcrm)/(pcrd+pcrm+floor_crmom)), &
                         & ABS((pcrm-pcrg)/(pcrm+pcrg+floor_crmom)) )
                    ok(i) = ok(i) .or. errcr > err_grad_crmom(icr)
                 end do
              end if
           end do
        end do

        ! Count newly flagged cells
        nok=0
        do i=1,ngrid
           if(flag1(ind_cell(i))==0.and.ok(i))then
              nok=nok+1
           end if
        end do

        do i=1,ngrid
           if(ok(i))flag1(ind_cell(i))=1
        end do

        nflag=nflag+nok
     end do
     ! End loop over cells

  end do
  ! End loop over grids

end subroutine cr_hydro_flag
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
