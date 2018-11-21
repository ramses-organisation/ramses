subroutine write_screen
  use amr_commons
  use hydro_commons
  use pm_commons
  use mpi_mod
  implicit none
  !
  integer::igrid,jgrid,ind,icpu,info
  integer::i,icell,ncell,ilevel,ncache
  integer::icellmin,nx_loc
  real(dp)::dx,scale,smallp,ddd,ppp

  integer     ,dimension(:),allocatable::ind_grid,ind_cell,ind_sort,ll,ll_all
  real(dp),dimension(:),allocatable::rr,et,ei,dd,uu,vv,ww,mm,gg,dtot!!
  real(dp),dimension(:),allocatable::rr_all,et_all,ei_all!!
  real(dp),dimension(:),allocatable::dd_all,uu_all,mm_all,gg_all,dtot_all!!
  real(dp),allocatable,dimension(:,:)::qq ! primitive variables
  integer,dimension(1:ncpu)::iskip,ncell_loc,ncell_all
  real(dp) :: lor,entho ! Lorentz factor
  real(dp) :: D,M,E,Mx,My,Mz,u2,Xsi,R
  real(dp) ::rho,p,vpar,vx,vy,vz
  real(dp) :: small_bigD=1e-12

  if(ndim>1)return
  smallp = smallc**2/gamma

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
  call MPI_ALLREDUCE(ncell_loc,ncell_all,ncpu,MPI_INTEGER,MPI_SUM,&
       & MPI_COMM_WORLD,info)
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

  allocate(rr(1:ncell),mm(1:ncell),dd(1:ncell),dtot(1:ncell))
  allocate(et(1:ncell),ei(1:ncell))
  allocate(uu(1:ncell),ll(1:ncell),gg(1:ncell),vv(1:ncell),ww(1:ncell))
  allocate(rr_all(1:ncell),mm_all(1:ncell),dd_all(1:ncell),dtot_all(1:ncell))
  allocate(et_all(1:ncell),ei_all(1:ncell))
  allocate(uu_all(1:ncell),ll_all(1:ncell),gg_all(1:ncell))
  allocate(qq(1:ncell,1:nvar))
  rr=0.0D0; mm=0.0D0; dd=0.0D0; dtot=0.0D0; et=0.0D0
  ei=0.0D0; uu=0.0D0; vv=0.0d0 ; ww=0.0d0 ; gg=0.0D0; ll=0
  rr_all=0.0D0; mm_all=0.0D0; dd_all=0.0D0; dtot_all=0.0D0; et_all=0.0D0
  ei_all=0.0D0; uu_all=0.0D0; gg_all=0.0D0; ll_all=0
  qq=0.0d0

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
        ! get radii
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
              do i=1,ncache
                 if(son(ind_cell(i))==0)then
                    icell=icell+1

                    !convert to primitive variables

                    ! Compute density
                    D = uold(ind_cell(i),1)
                    ! Compute momentum
                    Mx=uold(ind_cell(i),2) ; My=uold(ind_cell(i),3) ; Mz=uold(ind_cell(i),4)
                    M = sqrt(Mx**2+My**2+Mz**2)
                    ! Compute total energy
                    E = uold(ind_cell(i),5)
                    !Method from Mignone,McKinney,2007. Same as BS2011 except one uses E'=U-D and u^2=Lor^2*v^2

                    if (D<0) then
                       write(*,*),'D<0 in ctoprimbis'
                       D=small_bigD
                       uold(ind_cell(i),1)=D
                    endif
                    if (E<0) then
                       E=sqrt(M**2+d**2+1d-8)
                       write(*,*),'E<0 in ctoprimbis'
                       uold(ind_cell(i),5)=E
                    endif

                    if (E**2<M**2+D**2) then

                       write (*,*) 'Switch...ctoprimbis'
                       E=sqrt(M**2+d**2+1d-8)
                       uold(ind_cell(i),5)=E
                    endif

                    if ( M .eq. 0) then
                       qq(icell,1) = D
                       qq(icell,2) = 0d0
                       qq(icell,3) = 0d0
                       qq(icell,4) = 0d0
                       if (eos .eq. 'TM') then
                          qq(icell,5) =(E**2-D**2)/3d0/E
                       else
                          qq(icell,5)=(E-D)*(gamma-1d0)
                       endif
                       lor=1d0
                    else

                       call Newton_Raphson_Mignone(D,M,E,gamma,R)

                       ! Compute the Lorentz factor
                       u2  = M**2.0d0/(R**2.0d0-M**2.0d0)
                       lor = (1.0d0+u2)**(1d0/2d0)

                       ! Compute the density
                       qq(icell,1) = D/lor

                       ! compute velocities
                       qq(icell,2) = Mx/R
                       qq(icell,3) = My/R
                       qq(icell,4) = Mz/R

                       ! Compute pressure
                       Xsi=((R-D)-u2/(lor+1d0)*D)/lor**2
                       if (eos .eq. 'TM') then
                          rho=qq(icell,1)
                          qq(icell,5)=(2d0*xsi*(xsi+2d0*rho))/(5d0*(xsi+rho)+sqrt(9d0*xsi**2+18d0*rho*xsi+25d0*rho**2))
                       else
                          qq(icell,5)=(gamma-1d0)/gamma*Xsi
                       endif



                    endif
                    if ((qq(icell,1)<0d0).or.(qq(icell,5)<0d0).or.E<0d0 ) then
                       write(*,*) 'negative pressure or density output'
                       stop
                    endif
                    dd(icell)=qq(icell,1)
                    uu(icell)=qq(icell,2)
                    vv(icell)=qq(icell,3)
                    ww(icell)=qq(icell,4)
                    ei(icell)=qq(icell,5)
                 endif
              end do
           end do
        end if
        deallocate(ind_grid, ind_cell)
     end if
  end do

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(rr,rr_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(mm,mm_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(dd,dd_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(dtot,dtot_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(et,et_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(ei,ei_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(uu,uu_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(gg,gg_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(ll,ll_all,ncell,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(rr,rr_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  rr=rr_all; mm=mm_all; dd=dd_all; dtot=dtot_all; et=et_all
  ei=ei_all; uu=uu_all; gg=gg_all; ll=ll_all
#endif

  if(myid==1)then
     write(*,*)'====================================================================='
        write(*,*)'lev      x           d          u         v         w          P'
     endif
     ! Sort radius
     allocate(ind_sort(1:ncell))
     call quick_sort(rr,ind_sort,ncell)
     ! Write results to screen
     smallp=smallc**2/gamma
     nx_loc=icoarse_max-icoarse_min+1
     scale=boxlen/dble(nx_loc)
     ! Prevent underflow for velocity
     do i=1,ncell
        if(ABS(uu(i))<smallc)uu(i)=0.0D0
        if(ABS(vv(i))<smallc)vv(i)=0.0D0
        if(ABS(ww(i))<smallc)ww(i)=0.0D0
     end do
     do i=1,ncell
        ddd=MAX(dd(ind_sort(i)),smallr)
        ppp=ei(ind_sort(i))
        write(*,116) &
             & ll(ind_sort(i)),  &
             & (rr(i)-dble(icoarse_min))*scale, &
             & ddd , &
             & uu(ind_sort(i)), vv(ind_sort(i)), ww(ind_sort(i)),&
             & ei(ind_sort(i))
     end do
     deallocate(ind_sort)
     write(*,*)'================================================'
  end if

  ! Deallocate local arrays
  deallocate(mm,rr,dd,dtot,et,ei,uu,vv,ww,ll,gg)
  deallocate(mm_all,rr_all,dd_all,dtot_all,et_all,ei_all,uu_all,ll_all,gg_all)

!end if

#ifndef WITHOUTMPI
  call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif

111 format(2(1pe12.5,1x))
112 format(i3,1x,1pe10.3,1x,8(1pe10.3,1x))
113 format(i3,1x,1pe12.5,1x,9(1pe10.3,1x))
114 format(' Output ',i5,' cells')
115 format(' Output ',i5,' parts')
116 format(i3,1x,1pe12.5,1x,9(1pe15.8,1x))!!!
end subroutine write_screen
