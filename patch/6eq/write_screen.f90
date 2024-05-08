subroutine write_screen
  use amr_commons
  use hydro_commons
  use pm_commons
  use poisson_commons
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::info
#endif
  integer::igrid,jgrid,ind,icpu
  integer::i,icell,ncell,ilevel,ncache
  integer::icellmin,imat,lll
  real(dp)::dx,scale,ddd

  integer     ,dimension(:),allocatable::ind_grid,ind_cell,ind_sort,ll,ll_all
  real(kind=8),dimension(:),allocatable::rr,rr_all
  real(kind=8),dimension(:,:),allocatable::qq,qq_all,ff,ff_all,gg,gg_all,pp
  real(dp),dimension(1:nvector),save::ppp_mat,ccc,dtot
  real(dp),dimension(1:nvector,1:nmat),save::fff,ggg
  real(dp),dimension(1:nvector,1:npri),save::qqq
  logical::inv

  integer,dimension(1:ncpu)::iskip,ncell_loc,ncell_all
  real(dp)::uuu,rrr

  if(ndim>1)return

#ifndef WITHOUTMPI
  call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif

  ncell=0
  do ilevel=1,nlevelmax
     ncache=numbl(myid,ilevel)
     if(ncache > 0)then
        allocate(ind_grid(1:ncache),ind_cell(1:ncache))
        ! Gather all grids
        igrid=headl(myid,ilevel)
        do jgrid=1,ncache
           ind_grid(jgrid)=igrid
           igrid=next(igrid)
        end do
        ! Count leaf cells
        do ind=1,twotondim
           do i=1,ncache
              ind_cell(i)=ncoarse+(ind-1)*ngridmax+ind_grid(i)
           end do
           do i=1,ncache
              if(son(ind_cell(i))== 0)then
                 ncell=ncell+1
              end if
           end do
        end do
        deallocate(ind_grid, ind_cell)
     end if
  end do

  ncell_loc=0
  ncell_all=0
  ncell_loc(myid)=ncell
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(ncell_loc,ncell_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  ncell_all=ncell_loc
#endif

  ncell=0
  iskip=0
  do icpu=1,ncpu
     iskip(icpu)=ncell
     ncell=ncell+ncell_all(icpu)
  end do

  if(myid==1)write(*,114)ncell

  if(ncell>0)then

  allocate(rr(1:ncell)       ,rr_all(1:ncell))
  allocate(ll(1:ncell)       ,ll_all(1:ncell))
  allocate(ff(1:ncell,1:nmat),ff_all(1:ncell,1:nmat))
  allocate(gg(1:ncell,1:nmat),gg_all(1:ncell,1:nmat))
  allocate(qq(1:ncell,1:npri),qq_all(1:ncell,1:npri))
  rr    =0.0D0; ll    =0; ff    =0.0; gg    =0.0; qq    =0.0
  rr_all=0.0D0; ll_all=0; ff_all=0.0; gg_all=0.0; qq_all=0.0

  icell=iskip(myid)
  do ilevel=1,nlevelmax
     icellmin=icell
     ncache=numbl(myid,ilevel)
     if(ncache > 0)then
        dx=0.5D0**ilevel
        allocate(ind_grid(1:ncache),ind_cell(1:ncache))
        ! Gather all grids
        igrid=headl(myid,ilevel)
        do jgrid=1,ncache
           ind_grid(jgrid)=igrid
           igrid=next(igrid)
        end do
        ! Gather variables
        icell=icellmin
        do ind=1,twotondim
           do i=1,ncache
              ind_cell(i)=ncoarse+(ind-1)*ngridmax+ind_grid(i)
           end do
           do i=1,ncache
              if(son(ind_cell(i))==0)then
                 icell=icell+1
                 rr(icell)=xg(ind_grid(i),1)+(dble(ind)-1.5D0)*dx
                 ll(icell)=ilevel
              end if
           end do
        end do
        if(hydro)then
           icell=icellmin
           do ind=1,twotondim
              do i=1,ncache
                 ind_cell(i)=ncoarse+(ind-1)*ngridmax+ind_grid(i)
              end do
              dtot(1:ncache) = 0.0
              do i=1,ncache
                do imat=1,nmat
                  dtot(i)    = dtot(i) + uold(ind_cell(i),nmat+imat)
                end do
                if(son(ind_cell(i))==0)then
                  icell=icell+1
                  qq(icell,1)                = uold(ind_cell(i),2*nmat+1)/dtot(i)
                  do imat=1,nmat
                    ff(icell,imat)           = uold(ind_cell(i),imat)
                    gg(icell,imat)           = uold(ind_cell(i),nmat+imat)/ff(icell,imat)
                    qq(icell,ndim+nmat+imat) = uold(ind_cell(i),2*nmat+ndim+imat)/ff(icell,imat) - 0.5*gg(icell,imat)*qq(icell,1)**2
                    inv=.false.
                    call eos(gg(icell,imat),qq(icell,ndim+nmat+imat),ppp_mat,ccc,imat,inv,1)
                    qq(icell,ndim+imat)      = ppp_mat(1)
                  end do
                end if
              end do
           end do
        end if
        deallocate(ind_grid, ind_cell)
     end if
  end do


#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(ff,ff_all,ncell*nmat,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(gg,gg_all,ncell*nmat,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(qq,qq_all,ncell*npri,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(ll,ll_all,ncell     ,MPI_INTEGER         ,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(rr,rr_all,ncell     ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  ff=ff_all; gg=gg_all; qq=qq_all; rr=rr_all; ll=ll_all
#endif

  if(myid==1)then
     write(*,*)'================================================'
     write(*,*)'lev      x           f1         f2         d1         d2          u         P1         P2         e1         e2'
     ! Sort radius
     allocate(ind_sort(1:ncell))
     call quick_sort(rr,ind_sort,ncell)


     ! Write results to screen
     scale=boxlen/dble(icoarse_max-icoarse_min+1)
     do i=1,ncell

        fff(1,1:nmat)=ff(ind_sort(i),1:nmat)
        ggg(1,1:nmat)=gg(ind_sort(i),1:nmat)
        qqq(1,1:npri)=qq(ind_sort(i),1:npri)
        lll=ll(ind_sort(i))
        rrr=(rr(i)-dble(icoarse_min))*scale
        uuu=qqq(1,1)

        if(ABS(uuu).lt.1d-15)uuu=0.0
        write(*,113) &
             & lll,  &
             & rrr,  &
             & (fff(1,imat),imat=1,nmat) , &
             & (ggg(1,imat),imat=1,nmat) , &
             & uuu, &
             & (qqq(1,ndim+imat),imat=1,nmat), &
             & (qqq(1,ndim+nmat+imat),imat=1,nmat)
     end do
     deallocate(ind_sort)
  end if

  ! Deallocate local arrays
  deallocate(ff,qq,rr,ll)
  deallocate(ff_all,qq_all,rr_all,ll_all)

  end if

#ifndef WITHOUTMPI
  call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif

111 format(2(1pe12.5,1x))
112 format(i3,1x,1pe10.3,1x,8(1pe10.3,1x))
113 format(i3,1x,1pe12.5,1x,9(1pe10.3,1x))
114 format(' Output ',i5,' cells')
115 format(' Output ',i5,' parts')

end subroutine write_screen
