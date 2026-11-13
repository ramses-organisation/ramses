! CR refinement pass: flags cells whose CR energy gradient exceeds
! err_grad_cr, OR-ing the result into the shared flag1.
subroutine cr_hydro_flag(ilevel)
  use amr_commons
  use cr_parameters, only: ncr_groups,ncrvar,Ecr_idx,err_grad_cr
  use cr_hydro_commons, only: cruold
  implicit none
  integer::ilevel
  ! -------------------------------------------------------------------
  ! Flag for refinement cells whose CR energy gradient at level ilevel
  ! exceeds err_grad_cr. Separate pass; results OR into flag1.
  ! -------------------------------------------------------------------
  integer::i,j,ncache,nok,iskip
  integer::igrid,ind,idim,ngrid
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  integer,dimension(1:nvector,0:twondim),save::igridn
  integer,dimension(1:nvector,1:twondim),save::indn

  logical,dimension(1:nvector),save::ok

  real(dp),dimension(1:nvector,1:ncrvar),save::ucrg,ucrm,ucrd
  real(dp)::pcrg,pcrm,pcrd,errcr
  integer::icr,icrE
  ! Floor in the CR-gradient denominator (avoids divide-by-zero)
  real(dp),parameter::floor_crmom=1d-10

  if(ilevel==nlevelmax)return
  if(numbtot(1,ilevel)==0)return
  ! No CR refinement requested -> nothing to do.
  if(all(err_grad_cr(1:ncr_groups) < 0d0))return

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
           ! CR-energy gradient refinement (err_grad_cr).
           do icr=1,ncr_groups
              if(err_grad_cr(icr) >= 0.)then
                 icrE=Ecr_idx(icr)
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
                    ok(i) = ok(i) .or. errcr > err_grad_cr(icr)
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
