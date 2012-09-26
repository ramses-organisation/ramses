subroutine make_sink(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use cooling_module, ONLY: XH=>X, rhoc, mH 
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !----------------------------------------------------------------------
  ! Description: This subroutine create sink particle in cells where some
  ! density threshold is crossed. It also removes from the gas the
  ! corresponding particle mass. On exit, all fluid variables in the cell
  ! are modified. This is done only in leaf cells.
  ! Romain Teyssier, October 7th, 2007
  !----------------------------------------------------------------------
  ! local constants
  real(dp),dimension(1:twotondim,1:3)::xc
  integer ::ncache,nnew,ivar,ngrid,icpu,index_sink,index_sink_tot,icloud
  integer ::igrid,ix,iy,iz,ind,i,j,n,iskip,isink,inew,nx_loc
  integer ::ii,jj,kk,ind_cloud,ncloud
  integer ::ntot,ntot_all,info
  logical ::ok_free
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m
  real(dp)::d,x,y,z,u,v,w,e,temp,zg,factG
  real(dp)::dxx,dyy,dzz,drr
  real(dp)::msink_max2,rsink_max2
  real(dp)::velc,uc,vc,wc,l_jeans,d_jeans,d_thres,d_sink
  real(dp)::birth_epoch,xx,yy,zz,rr
  real(dp),dimension(1:3)::skip_loc
  real(dp)::dx,dx_loc,scale,vol_loc,dx_min,vol_min
  real(dp)::bx1,bx2,by1,by2,bz1,bz2

  integer ,dimension(1:nvector),save::ind_grid,ind_cell
  integer ,dimension(1:nvector),save::ind_grid_new,ind_cell_new,ind_part
  integer ,dimension(1:nvector),save::ind_part_cloud,ind_grid_cloud
  logical ,dimension(1:nvector),save::ok,ok_new=.true.,ok_true=.true.
  integer ,dimension(1:ncpu)::ntot_sink_cpu,ntot_sink_all
  


  if(numbtot(1,ilevel)==0) return
  if(.not. hydro)return
  if(ndim.ne.3)return

  if(verbose)write(*,*)' Entering make_sink for level ',ilevel

  ! Conversion factor from user units to cgs units                              
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_m=scale_d*scale_l**3

  ! Minimum radius to create a new sink from any other
  rsink_max2=(rsink_max*3.08d21/scale_l)**2

  ! Maximum value for the initial sink mass
  msink_max2=msink_max*2d33/scale_m
  
  ! Gravitational constant
  factG=1d0
  if(cosmo)factG=3d0/8d0/3.1415926*omega_m*aexp

  ! Density threshold for sink particle creation
  d_sink=n_sink/scale_nH

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
  dx_min=(0.5D0**nlevelmax)*scale
  vol_min=dx_min**ndim

  ! Birth epoch
  birth_epoch=t

  ! Cells center position relative to grid center position
  do ind=1,twotondim  
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     xc(ind,1)=(dble(ix)-0.5D0)*dx
     xc(ind,2)=(dble(iy)-0.5D0)*dx
     xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  xx=0.0; yy=0.0;zz=0.0
  ncloud=0
  do kk=-2*ir_cloud,2*ir_cloud
     zz=dble(kk)*0.5
     do jj=-2*ir_cloud,2*ir_cloud
        yy=dble(jj)*0.5
        do ii=-2*ir_cloud,2*ir_cloud
           xx=dble(ii)*0.5
           rr=sqrt(xx*xx+yy*yy+zz*zz)
           if(rr<=dble(ir_cloud))ncloud=ncloud+1
        end do
     end do
  end do
  ncloud_sink=ncloud

  ! Set new sink variables to zero
  msink_new=0d0; tsink_new=0d0; delta_mass_new=0d0; xsink_new=0d0; vsink_new=0d0; oksink_new=0d0; idsink_new=0

#if NDIM==3

  !------------------------------------------------
  ! Convert hydro variables to primitive variables
  !------------------------------------------------
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     do ind=1,twotondim  
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        do i=1,ngrid
           d=uold(ind_cell(i),1)
           u=uold(ind_cell(i),2)/d
           v=uold(ind_cell(i),3)/d
           w=uold(ind_cell(i),4)/d
           e=uold(ind_cell(i),5)/d
#ifdef SOLVERmhd
           bx1=uold(ind_cell(i),6)
           by1=uold(ind_cell(i),7)
           bz1=uold(ind_cell(i),8)
           bx2=uold(ind_cell(i),nvar+1)
           by2=uold(ind_cell(i),nvar+2)
           bz2=uold(ind_cell(i),nvar+3)
           e=e-0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)/d
#endif
           e=e-0.5d0*(u**2+v**2+w**2)
           uold(ind_cell(i),1)=d
           uold(ind_cell(i),2)=u
           uold(ind_cell(i),3)=v
           uold(ind_cell(i),4)=w
           uold(ind_cell(i),5)=e
        end do
        do ivar=imetal,nvar
           do i=1,ngrid
              d=uold(ind_cell(i),1)
              w=uold(ind_cell(i),ivar)/d
              uold(ind_cell(i),ivar)=w
           end do
        enddo
     end do
  end do

  !-------------------------------
  ! Determine SMBH formation sites
  !-------------------------------
  if(smbh)call quenching(ilevel)

  !----------------------------
  ! Compute number of new sinks
  !----------------------------
  ntot=0
  ! Loop over grids
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     ! Density threshold crossed ---> logical array ok(i)
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do

        ! Flag leaf cells
        do i=1,ngrid
           ok(i)=son(ind_cell(i))==0
        end do

        ! Create new sink if the gas density exceed some threshold
        do i=1,ngrid
           d=uold(ind_cell(i),1)

           ! User defined density threshold
           d_thres=d_sink

           ! Jeans length related density threshold
           if(d_thres<0.0)then
              temp=max(uold(ind_cell(i),5)*(gamma-1.0),smallc**2)
              d_jeans=temp*3.1415926/(4.0*dx_loc)**2/factG
              d_thres=d_jeans
           endif

           ! Density criterion
           if(d < d_thres)ok(i)=.false.

           ! Quenching criterion
           if(smbh.and.flag2(ind_cell(i))==1)ok(i)=.false.
           
           ! Geometrical criterion
           if(ivar_refine>0)then
              d=uold(ind_cell(i),ivar_refine)
              if(d<=var_cut_refine)ok(i)=.false.
           endif

           ! Proximity criterion
           if(ok(i).and.rsink_max>0d0)then
              x=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
              y=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
              z=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
              do isink=1,nsink
                 dxx=x-xsink(isink,1)
                 if(dxx> 0.5*scale)then
                    dxx=dxx-scale
                 endif
                 if(dxx<-0.5*scale)then
                    dxx=dxx+scale
                 endif
                 dyy=y-xsink(isink,2)
                 if(dyy> 0.5*scale)then
                    dyy=dyy-scale
                 endif
                 if(dyy<-0.5*scale)then
                    dyy=dyy+scale
                 endif
                 dzz=z-xsink(isink,3)
                 if(dzz> 0.5*scale)then
                    dzz=dzz-scale
                 endif
                 if(dzz<-0.5*scale)then
                    dzz=dzz+scale
                 endif
                 drr=dxx*dxx+dyy*dyy+dzz*dzz
                 if(drr.le.rsink_max2)ok(i)=.false.
              enddo
           endif

        end do

        ! Calculate number of new sinks in each cell
        do i=1,ngrid
           flag2(ind_cell(i))=0
           if(ok(i))then
              ntot=ntot+1
              flag2(ind_cell(i))=1
           endif
        enddo
     end do
  end do

  !---------------------------------
  ! Check for free particle memory
  !--------------------------------
  ok_free=(numbp_free-ntot)>=0
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(numbp_free,numbp_free_tot,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  numbp_free_tot=numbp_free
#endif
  if(.not. ok_free)then
     write(*,*)'No more free memory for particles'
     write(*,*)'New sink particles',ntot
     write(*,*)'Increase npartmax'
#ifndef WITHOUTMPI
    call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
    stop
#endif
  end if
  
  !---------------------------------
  ! Compute global sink statistics
  !---------------------------------
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(ntot,ntot_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  ntot_all=ntot
#endif
#ifndef WITHOUTMPI
  ntot_sink_cpu=0; ntot_sink_all=0
  ntot_sink_cpu(myid)=ntot
  call MPI_ALLREDUCE(ntot_sink_cpu,ntot_sink_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  ntot_sink_cpu(1)=ntot_sink_all(1)
  do icpu=2,ncpu
     ntot_sink_cpu(icpu)=ntot_sink_cpu(icpu-1)+ntot_sink_all(icpu)
  end do
#endif
  nsink=nsink+ntot_all  
  nindsink=nindsink+ntot_all
  if(myid==1)then
     if(ntot_all.gt.0)then
        write(*,'(" Level = ",I6," New sink before merging= ",I6," Tot =",I8)')&
             & ilevel,ntot_all,nsink
     endif
  end if

  !------------------------------
  ! Create new sink particles
  !------------------------------
  ! Starting identity number
  if(myid==1)then
     index_sink=nsink-ntot_all
     index_sink_tot=nindsink-ntot_all
  else
     index_sink=nsink-ntot_all+ntot_sink_cpu(myid-1)
     index_sink_tot=nindsink-ntot_all+ntot_sink_cpu(myid-1)
  end if

  ! Loop over grids
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
        
        ! Gather cells with a new sink
        nnew=0
        do i=1,ngrid
           if (flag2(ind_cell(i))>0)then
              nnew=nnew+1
              ind_grid_new(nnew)=ind_grid(i)
              ind_cell_new(nnew)=ind_cell(i)
           end if
        end do
        
        ! Create new sink particles
        do i=1,nnew
           index_sink=index_sink+1
           index_sink_tot=index_sink_tot+1

           ! Get gas variables
           d=uold(ind_cell_new(i),1)
           u=uold(ind_cell_new(i),2)
           v=uold(ind_cell_new(i),3)
           w=uold(ind_cell_new(i),4)
           e=uold(ind_cell_new(i),5)

           ! Get gas cell position
           x=(xg(ind_grid_new(i),1)+xc(ind,1)-skip_loc(1))*scale
           y=(xg(ind_grid_new(i),2)+xc(ind,2)-skip_loc(2))*scale
           z=(xg(ind_grid_new(i),3)+xc(ind,3)-skip_loc(3))*scale

           ! User defined density threshold
           d_thres=d_sink
           
           ! Jeans length related density threshold
           if(d_thres<0.0)then
              temp=max(e*(gamma-1.0),smallc**2)
              d_jeans=temp*3.1415926/(4.0*dx_loc)**2/factG
              d_thres=d_jeans
           endif

           ! Mass of the new sink
           msink_new(index_sink)=min((d-d_thres)*vol_loc,msink_max2)
           delta_mass_new(index_sink)=0d0

           ! Global index of the new sink
           oksink_new(index_sink)=1d0
           idsink_new(index_sink)=index_sink_tot

           ! Update linked list
           ind_grid_cloud(1)=ind_grid_new(i)
           call remove_free(ind_part_cloud,1)
           call add_list(ind_part_cloud,ind_grid_cloud,ok_true,1)
           ind_cloud=ind_part_cloud(1)

           ! Set new sink particle variables
           tp(ind_cloud)=birth_epoch     ! Birth epoch
           mp(ind_cloud)=msink_new(index_sink) ! Mass
           levelp(ind_cloud)=ilevel      ! Level
           idp(ind_cloud)=-index_sink    ! Identity
           xp(ind_cloud,1)=x
           xp(ind_cloud,2)=y
           xp(ind_cloud,3)=z
           vp(ind_cloud,1)=u
           vp(ind_cloud,2)=v
           vp(ind_cloud,3)=w
           
           ! Store properties of the new sink
           tsink_new(index_sink)=birth_epoch
           xsink_new(index_sink,1)=x
           xsink_new(index_sink,2)=y
           xsink_new(index_sink,3)=z
           vsink_new(index_sink,1)=u
           vsink_new(index_sink,2)=v
           vsink_new(index_sink,3)=w
           
           uold(ind_cell_new(i),1)=uold(ind_cell_new(i),1)-msink_new(index_sink)/vol_loc
           
        end do
        ! End loop over new sink particle cells

     end do
     ! End loop over cells
  end do
  ! End loop over grids

  !---------------------------------------------------------
  ! Convert hydro variables back to conservative variables
  !---------------------------------------------------------
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     do ind=1,twotondim  
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        do i=1,ngrid
           d=uold(ind_cell(i),1)
           u=uold(ind_cell(i),2)
           v=uold(ind_cell(i),3)
           w=uold(ind_cell(i),4)
           e=uold(ind_cell(i),5)
#ifdef SOLVERmhd
           bx1=uold(ind_cell(i),6)
           by1=uold(ind_cell(i),7)
           bz1=uold(ind_cell(i),8)
           bx2=uold(ind_cell(i),nvar+1)
           by2=uold(ind_cell(i),nvar+2)
           bz2=uold(ind_cell(i),nvar+3)
           e=e+0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)/d
#endif
           e=e+0.5d0*(u**2+v**2+w**2)
           uold(ind_cell(i),1)=d
           uold(ind_cell(i),2)=d*u
           uold(ind_cell(i),3)=d*v
           uold(ind_cell(i),4)=d*w
           uold(ind_cell(i),5)=d*e
        end do
        do ivar=imetal,nvar
           do i=1,ngrid
              d=uold(ind_cell(i),1)
              w=uold(ind_cell(i),ivar)
              uold(ind_cell(i),ivar)=d*w
           end do
        end do
     end do
  end do

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(oksink_new,oksink_all,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(idsink_new,idsink_all,nsinkmax,MPI_INTEGER         ,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(msink_new ,msink_all ,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(tsink_new ,tsink_all ,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(xsink_new ,xsink_all ,nsinkmax*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(vsink_new ,vsink_all ,nsinkmax*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(delta_mass_new,delta_mass_all,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#else
  oksink_all=oksink_new
  idsink_all=idsink_new
  msink_all=msink_new
  tsink_all=tsink_new
  xsink_all=xsink_new
  vsink_all=vsink_new
  delta_mass_all=delta_mass_new
#endif
  do isink=1,nsink
     if(oksink_all(isink)==1d0)then
        idsink(isink)=idsink_all(isink)
        msink(isink)=msink_all(isink)
        tsink(isink)=tsink_all(isink)
        xsink(isink,1:ndim)=xsink_all(isink,1:ndim)
        vsink(isink,1:ndim)=vsink_all(isink,1:ndim)
        delta_mass(isink)=delta_mass_all(isink)
     endif
  end do

#endif

end subroutine make_sink
!################################################################
!################################################################
!################################################################
!################################################################
subroutine merge_sink(ilevel)
  use pm_commons
  use amr_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !------------------------------------------------------------------------
  ! This routine merges sink usink the FOF algorithm.
  ! It keeps only the group centre of mass and remove other sinks.
  !------------------------------------------------------------------------
  integer::j,isink,ii,jj,kk,ind,idim,new_sink
  real(dp)::dx_loc,scale,dx_min,xx,yy,zz,rr,rmax2,rmax
  integer::igrid,jgrid,ipart,jpart,next_part,info
  integer::i,ig,ip,npart1,npart2,icpu,nx_loc
  integer::igrp,icomp,gndx,ifirst,ilast,indx
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part
  integer,dimension(:),allocatable::psink,gsink
  real(dp),dimension(1:3)::xbound,skip_loc


  if(numbtot(1,ilevel)==0)return
  if(nsink==0)return
  if(verbose)write(*,111)ilevel

  ! Mesh spacing in that level
  dx_loc=0.5D0**ilevel
  xbound(1:3)=(/dble(nx),dble(ny),dble(nz)/)
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax/aexp
  rmax=dble(ir_cloud)*dx_min ! Linking length in physical units
  rmax2=rmax*rmax

  allocate(psink(1:nsink),gsink(1:nsink))
  
  !-------------------------------
  ! Merge sinks using FOF
  !-------------------------------
  do isink=1,nsink
     psink(isink)=isink
     gsink(isink)=0
  end do
  
  igrp=0
  icomp=1
  ifirst=2
  do while(icomp.le.nsink)
     gndx=psink(icomp)
     if(gsink(gndx)==0)then
        igrp=igrp+1
        gsink(gndx)=igrp
     endif
     ilast=nsink
     do while((ilast-ifirst+1)>0)
        indx=psink(ifirst)
        xx=xsink(indx,1)-xsink(gndx,1)
        if(xx>scale*xbound(1)/2.0)then
           xx=xx-scale*xbound(1)
        endif
        if(xx<-scale*xbound(1)/2.0)then
           xx=xx+scale*xbound(1)
        endif
        rr=xx**2
        yy=xsink(indx,2)-xsink(gndx,2)
        if(yy>scale*xbound(2)/2.0)then
           yy=yy-scale*xbound(2)
        endif
        if(yy<-scale*xbound(2)/2.0)then
           yy=yy+scale*xbound(2)
        endif
        rr=yy**2+rr
        zz=xsink(indx,3)-xsink(gndx,3)
        if(zz>scale*xbound(3)/2.0)then
           zz=zz-scale*xbound(3)
        endif
        if(zz<-scale*xbound(3)/2.0)then
           zz=zz+scale*xbound(3)
        endif
        rr=zz**2+rr
        if(rr.le.rmax2)then
           ifirst=ifirst+1
           gsink(indx)=igrp
        else
           psink(ifirst)=psink(ilast)
           psink(ilast)=indx
           ilast=ilast-1
        endif
     end do
     icomp=icomp+1
  end do
  new_sink=igrp

  if(myid==1)then
     write(*,*)'Number of sinks after merging',new_sink
!     do isink=1,nsink
!        write(*,'(3(I3,1x),3(1PE10.3))')isink,psink(isink),gsink(isink),xsink(isink,1:ndim)
!     end do
  endif
  
  !----------------------------------------------------
  ! Compute group centre of mass and average velocty
  !----------------------------------------------------
  xsink_new=0d0; vsink_new=0d0; msink_new=0d0; tsink_new=0d0; delta_mass_new=0d0
  oksink_all=0d0; oksink_new=0d0; idsink_all=0d0; idsink_new=0d0
  do isink=1,nsink
     igrp=gsink(isink)
     if(oksink_new(igrp)==0d0)then
        oksink_all(isink)=igrp
        oksink_new(igrp)=isink
     endif
     msink_new(igrp)=msink_new(igrp)+msink(isink)
     delta_mass_new(igrp)=delta_mass_new(igrp)+delta_mass(isink)
     if(tsink_new(igrp)==0d0)then
        tsink_new(igrp)=tsink(isink)
     else
        tsink_new(igrp)=min(tsink_new(igrp),tsink(isink))
     endif
     if(idsink_new(igrp)==0)then
        idsink_new(igrp)=idsink(isink)
     else
        idsink_new(igrp)=min(idsink_new(igrp),idsink(isink))
     endif

     xx=xsink(isink,1)-xsink(int(oksink_new(igrp)),1)
     if(xx>scale*xbound(1)/2.0)then
        xx=xx-scale*xbound(1)
     endif
     if(xx<-scale*xbound(1)/2.0)then
        xx=xx+scale*xbound(1)
     endif
     xsink_new(igrp,1)=xsink_new(igrp,1)+msink(isink)*xx
     vsink_new(igrp,1)=vsink_new(igrp,1)+msink(isink)*vsink(isink,1)
     yy=xsink(isink,2)-xsink(int(oksink_new(igrp)),2)
     if(yy>scale*xbound(2)/2.0)then
        yy=yy-scale*xbound(2)
     endif
     if(yy<-scale*xbound(2)/2.0)then
        yy=yy+scale*xbound(2)
     endif
     xsink_new(igrp,2)=xsink_new(igrp,2)+msink(isink)*yy
     vsink_new(igrp,2)=vsink_new(igrp,2)+msink(isink)*vsink(isink,2)
     zz=xsink(isink,3)-xsink(int(oksink_new(igrp)),3)
     if(zz>scale*xbound(3)/2.0)then
        zz=zz-scale*xbound(3)
     endif
     if(zz<-scale*xbound(3)/2.0)then
        zz=zz+scale*xbound(3)
     endif
     xsink_new(igrp,3)=xsink_new(igrp,3)+msink(isink)*zz
     vsink_new(igrp,3)=vsink_new(igrp,3)+msink(isink)*vsink(isink,3)
  end do
  do isink=1,new_sink
     xsink_new(isink,1)=xsink_new(isink,1)/msink_new(isink)+xsink(int(oksink_new(isink)),1)
     vsink_new(isink,1)=vsink_new(isink,1)/msink_new(isink)
     xsink_new(isink,2)=xsink_new(isink,2)/msink_new(isink)+xsink(int(oksink_new(isink)),2)
     vsink_new(isink,2)=vsink_new(isink,2)/msink_new(isink)
     xsink_new(isink,3)=xsink_new(isink,3)/msink_new(isink)+xsink(int(oksink_new(isink)),3)
     vsink_new(isink,3)=vsink_new(isink,3)/msink_new(isink)
  end do
  nsink=new_sink
  msink(1:nsink)=msink_new(1:nsink)
  tsink(1:nsink)=tsink_new(1:nsink)
  idsink(1:nsink)=idsink_new(1:nsink)
  delta_mass(1:nsink)=delta_mass_new(1:nsink)
  xsink(1:nsink,1:ndim)=xsink_new(1:nsink,1:ndim)
  vsink(1:nsink,1:ndim)=vsink_new(1:nsink,1:ndim)

  ! Periodic boundary conditions
  do isink=1,nsink
     xx=xsink(isink,1)
     if(xx<-scale*skip_loc(1))then
        xx=xx+scale*(xbound(1)-skip_loc(1))
     endif
     if(xx>scale*(xbound(1)-skip_loc(1)))then
        xx=xx-scale*(xbound(1)-skip_loc(1))
     endif
     xsink(isink,1)=xx
     yy=xsink(isink,2)
     if(yy<-scale*skip_loc(2))then
        yy=yy+scale*(xbound(2)-skip_loc(2))
     endif
     if(yy>scale*(xbound(2)-skip_loc(2)))then
        yy=yy-scale*(xbound(2)-skip_loc(2))
     endif
     xsink(isink,2)=yy
     zz=xsink(isink,3)
     if(zz<-scale*skip_loc(3))then
        zz=zz+scale*(xbound(3)-skip_loc(3))
     endif
     if(zz>scale*(xbound(3)-skip_loc(3)))then
        zz=zz-scale*(xbound(3)-skip_loc(3))
     endif
     xsink(isink,3)=zz
  enddo

  deallocate(psink,gsink)
  
  !-----------------------------------------------------
  ! Remove sink particles that are part of a FOF group.
  !-----------------------------------------------------
  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0
        
        ! Count sink particles
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if(idp(ipart).lt.0)then
                 npart2=npart2+1
              endif
              ipart=next_part  ! Go to next particle
           end do
        endif
        
        ! Gather sink particles
        if(npart2>0)then        
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              ! Select only sink particles
              if(idp(ipart).lt.0)then
                 if(ig==0)then
                    ig=1
                    ind_grid(ig)=igrid
                 end if
                 ip=ip+1
                 ind_part(ip)=ipart
                 ind_grid_part(ip)=ig   
              endif
              if(ip==nvector)then
                 call kill_sink(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
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
     if(ip>0)call kill_sink(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
  end do 
  ! End loop over cpus

111 format('   Entering merge_sink for level ',I2)

end subroutine merge_sink
!################################################################
!################################################################
!################################################################
!################################################################
subroutine kill_sink(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !-----------------------------------------------------------------------
  ! This routine is called by subroutine merge_sink
  ! It removes sink particles that are part of a FOF group.
  !-----------------------------------------------------------------------
  integer::j,isink,ii,jj,kk,ind,idim,isink_new
  real(dp)::dx_loc,scale,dx_min,xx,yy,zz,rr,rmax
  ! Particle-based arrays
  logical ,dimension(1:nvector),save::ok

  do j=1,np
     isink=-idp(ind_part(j))
     ok(j)=(oksink_all(isink)==0)
     if(.not. ok(j))then
        isink_new=oksink_all(isink)
        idp(ind_part(j))=-isink_new
        mp(ind_part(j))=msink(isink_new)
        xp(ind_part(j),1)=xsink(isink_new,1)
        vp(ind_part(j),1)=vsink(isink_new,1)
        xp(ind_part(j),2)=xsink(isink_new,2)
        vp(ind_part(j),2)=vsink(isink_new,2)
        xp(ind_part(j),3)=xsink(isink_new,3)
        vp(ind_part(j),3)=vsink(isink_new,3)
     endif
  end do

  ! Remove particles from parent linked list
  call remove_list(ind_part,ind_grid_part,ok,np)
  call add_free_cond(ind_part,ok,np)

end subroutine kill_sink
!################################################################
!################################################################
!################################################################
!################################################################

subroutine kill_cloud(ilevel)
  use pm_commons
  use amr_commons
  implicit none
  integer::ilevel
  !------------------------------------------------------------------------
  ! This routine removes from the list cloud particles and keeps only
  ! sink particles. 
  !------------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,next_part
  integer::i,ig,ip,npart1,npart2,icpu,nx_loc
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Gather sink and cloud particles.

  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0
        
        ! Count sink and cloud particles
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if(idp(ipart).lt.0)then
                 npart2=npart2+1
              endif
              ipart=next_part  ! Go to next particle
           end do
        endif
        
        ! Gather sink and cloud particles
        if(npart2>0)then        
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              ! Select only sink particles
              if(idp(ipart).lt.0)then
                 if(ig==0)then
                    ig=1
                    ind_grid(ig)=igrid
                 end if
                 ip=ip+1
                 ind_part(ip)=ipart
                 ind_grid_part(ip)=ig   
              endif
              if(ip==nvector)then
                 call rm_cloud(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
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
     if(ip>0)call rm_cloud(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
  end do

111 format('   Entering kill_cloud for level ',I2)

end subroutine kill_cloud
!################################################################
!################################################################
!################################################################
!################################################################
subroutine rm_cloud(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !-----------------------------------------------------------------------
  ! This routine is called by subroutine kill_cloud.
  !-----------------------------------------------------------------------
  logical::error
  integer::j,isink,ii,jj,kk,ind,idim,nx_loc
  real(dp)::dx_loc,scale,dx_min,xx,yy,zz,rr,r2,r2_eps
  ! Particle-based arrays
  logical,dimension(1:nvector),save::ok

  ! Mesh spacing in that level
  dx_loc=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax/aexp
  r2_eps=(1d-15*dx_min)**2

  do j=1,np
     isink=-idp(ind_part(j))
     r2=0d0
     do idim=1,ndim
        r2=r2+(xp(ind_part(j),idim)-xsink(isink,idim))**2
     end do
     ok(j)=r2>r2_eps
  end do
  
  ! Remove particles from parent linked list
  call remove_list(ind_part,ind_grid_part,ok,np)
  call add_free_cond(ind_part,ok,np)

end subroutine rm_cloud
!################################################################
!################################################################
!################################################################
!################################################################
subroutine agn_feedback
  use amr_commons
  use pm_commons
  use hydro_commons
  use cooling_module, ONLY: XH=>X, rhoc, mH 
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::nSN_tot_all
  integer,dimension(1:ncpu)::nSN_icpu_all
  real(dp),dimension(:),allocatable::mSN_all,sSN_all,ZSN_all
  real(dp),dimension(:,:),allocatable::xSN_all,vSN_all
#endif
  !----------------------------------------------------------------------
  ! Description: This subroutine checks SN events in cells where a
  ! star particle has been spawned.
  ! Yohan Dubois
  !----------------------------------------------------------------------
  ! local constants
  integer::ip,icpu,igrid,jgrid,npart1,npart2,ipart,jpart,next_part
  integer::nSN,nSN_loc,nSN_tot,info,isink,ilevel,ivar
  integer,dimension(1:ncpu)::nSN_icpu
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,t0
  real(dp)::scale,dx_min,vol_min,nISM,nCOM,d0,mstar,temp_blast
  real(dp)::T2_AGN,T2_min,T2_max,delta_mass_max
  integer::nx_loc
  integer,dimension(:),allocatable::ind_part,ind_grid
  logical,dimension(:),allocatable::ok_free

  if(.not. hydro)return
  if(ndim.ne.3)return

  if(verbose)write(*,*)'Entering make_sn'
  
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Mesh spacing in that level
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_min=(0.5D0**nlevelmax)*scale
  vol_min=dx_min**ndim

  ! AGN specific energy
  T2_AGN=0.15*1d12 ! in Kelvin

  ! Minimum specific energy
  T2_min=1d7  ! in Kelvin

  ! Maximum specific energy
  T2_max=1d9 ! in Kelvin

  ! Compute the grid discretization effects
  call average_AGN

  ! Check if sink goes into blast wave mode
  ok_blast_agn(1:nsink)=.false.
  do isink=1,nsink
     ! Compute estimated average temperature in the blast
     temp_blast=0.0
     if(vol_gas_agn(isink)>0.0)then
        temp_blast=T2_AGN*delta_mass(isink)/mass_gas_agn(isink)
     else
        if(ind_blast_agn(isink)>0)then
           temp_blast=T2_AGN*delta_mass(isink)/mass_blast_agn(isink)
        endif
     endif
     if(temp_blast>T2_min)then
        ok_blast_agn(isink)=.true.
     endif
  end do
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(ok_blast_agn,ok_blast_agn_all,nsink,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,info)
  ok_blast_agn=ok_blast_agn_all
#endif

  ! Modify hydro quantities to account for the AGN blast
  call AGN_blast

  ! Reset accreted mass
  do isink=1,nsink
     if(ok_blast_agn(isink))then
        if(myid==1)then
           write(*,'("***BLAST***",I4,1X,2(1PE12.5,1X))')isink &
                & ,msink(isink)*scale_d*scale_l**3/2d33 &  
                & ,delta_mass(isink)*scale_d*scale_l**3/2d33
        endif
        ! Compute estimated average temperature in the blast
        temp_blast=0.0
        if(vol_gas_agn(isink)>0.0)then
           temp_blast=T2_AGN*delta_mass(isink)/mass_gas_agn(isink)
        else
           if(ind_blast_agn(isink)>0)then
              temp_blast=T2_AGN*delta_mass(isink)/mass_blast_agn(isink)
           endif
        endif
        if(temp_blast<T2_max)then
           delta_mass(isink)=0.0
        else
           if(vol_gas_agn(isink)>0.0)then
              delta_mass_max=T2_max/T2_AGN*mass_gas_agn(isink)
           else
              if(ind_blast_agn(isink)>0)then
                 delta_mass_max=T2_max/T2_AGN*mass_blast_agn(isink)
              endif
           endif
           delta_mass(isink)=max(delta_mass(isink)-delta_mass_max,0.0_dp)
        endif
     endif
  end do

  ! Update hydro quantities for split cells
  do ilevel=nlevelmax,levelmin,-1
     call upload_fine(ilevel)
     do ivar=1,nvar
        call make_virtual_fine_dp(uold(1,ivar),ilevel)
     enddo
  enddo

end subroutine agn_feedback
!################################################################
!################################################################
!################################################################
!################################################################
subroutine average_AGN
  use pm_commons
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !------------------------------------------------------------------------
  ! This routine average the hydro quantities inside the SN bubble
  !------------------------------------------------------------------------
  integer::ilevel,ncache,nSN,j,isink,ind,ix,iy,iz,ngrid,iskip
  integer::i,nx_loc,igrid,info
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp)::x,y,z,drr,d,u,v,w,ek,u2,v2,w2,dr_cell
  real(dp)::scale,dx,dxx,dyy,dzz,dx_min,dx_loc,vol_loc,rmax2,rmax
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  logical ,dimension(1:nvector),save::ok

  if(nsink==0)return
  if(verbose)write(*,*)'Entering average_AGN'

  ! Mesh spacing in that level
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax

  ! Maximum radius of the ejecta
  rmax=4.0d0*dx_min/aexp
  rmax2=rmax*rmax

  ! Initialize the averaged variables
  vol_gas_agn=0.0;vol_blast_agn=0.0;mass_gas_agn=0.0;ind_blast_agn=-1

  ! Loop over levels
  do ilevel=levelmin,nlevelmax
     ! Computing local volume (important for averaging hydro quantities) 
     dx=0.5D0**ilevel 
     dx_loc=dx*scale
     vol_loc=dx_loc**ndim
     ! Cells center position relative to grid center position
     do ind=1,twotondim  
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        xc(ind,1)=(dble(ix)-0.5D0)*dx
        xc(ind,2)=(dble(iy)-0.5D0)*dx
        xc(ind,3)=(dble(iz)-0.5D0)*dx
     end do

     ! Loop over grids
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

           ! Flag leaf cells
           do i=1,ngrid
              ok(i)=son(ind_cell(i))==0
           end do

           do i=1,ngrid
              if(ok(i))then
                 ! Get gas cell position
                 x=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
                 y=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
                 z=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale

                 do isink=1,nsink
                    ! Check if the cell lies within the sink radius
                    dxx=x-xsink(isink,1)
                    if(dxx> 0.5*scale)then
                       dxx=dxx-scale
                    endif
                    if(dxx<-0.5*scale)then
                       dxx=dxx+scale
                    endif
                    dyy=y-xsink(isink,2)
                    if(dyy> 0.5*scale)then
                       dyy=dyy-scale
                    endif
                    if(dyy<-0.5*scale)then
                       dyy=dyy+scale
                    endif
                    dzz=z-xsink(isink,3)
                    if(dzz> 0.5*scale)then
                       dzz=dzz-scale
                    endif
                    if(dzz<-0.5*scale)then
                       dzz=dzz+scale
                    endif
                    drr=dxx*dxx+dyy*dyy+dzz*dzz
                    dr_cell=MAX(ABS(dxx),ABS(dyy),ABS(dzz))
                    if(drr.lt.rmax2)then
                       vol_gas_agn(isink)=vol_gas_agn(isink)+vol_loc
                       mass_gas_agn(isink)=mass_gas_agn(isink)+vol_loc*uold(ind_cell(i),1)
                    endif
                    if(dr_cell.le.dx_loc/2.0)then
                       ind_blast_agn(isink)=ind_cell(i)
                       vol_blast_agn(isink)=vol_loc
                       mass_blast_agn(isink)=vol_loc*uold(ind_cell(i),1)
                    endif
                 end do
              endif
           end do
           
        end do
        ! End loop over cells
     end do
     ! End loop over grids
  end do
  ! End loop over levels

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(vol_gas_agn,vol_gas_agn_all,nsink,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(mass_gas_agn,mass_gas_agn_all,nsink,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  vol_gas_agn=vol_gas_agn_all
  mass_gas_agn=mass_gas_agn_all
#endif

  if(verbose)write(*,*)'Exiting average_AGN'

end subroutine average_AGN
!################################################################
!################################################################
!################################################################
!################################################################
subroutine AGN_blast
  use pm_commons
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !------------------------------------------------------------------------
  ! This routine merges SN using the FOF algorithm.
  !------------------------------------------------------------------------
  integer::ilevel,j,isink,nSN,ind,ix,iy,iz,ngrid,iskip
  integer::i,nx_loc,igrid,info,ncache
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp)::x,y,z,dx,dxx,dyy,dzz,drr,d,u,v,w,ek,u_r,ESN
  real(dp)::scale,dx_min,dx_loc,vol_loc,rmax2,rmax,T2_AGN,T2_max
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  logical ,dimension(1:nvector),save::ok


  if(nsink==0)return
  if(verbose)write(*,*)'Entering AGN_blast'

  ! Mesh spacing in that level
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Maximum radius of the ejecta
  rmax=4.0d0*dx_min/aexp
  rmax2=rmax*rmax
  
  ! AGN specific energy
  T2_AGN=0.15*1d12 ! in Kelvin
  T2_AGN=T2_AGN/scale_T2 ! in code units

  ! Maximum specific energy
  T2_max=1d9 ! in Kelvin
  T2_max=T2_max/scale_T2 ! in code units

  do isink=1,nsink
     if(ok_blast_agn(isink))then
        if(vol_gas_agn(isink)>0d0)then
           p_agn(isink)=MIN(delta_mass(isink)*T2_AGN/vol_gas_agn(isink), &
                &         mass_gas_agn(isink)*T2_max/vol_gas_agn(isink)  )
        else
           p_agn(isink)=MIN(delta_mass(isink)*T2_AGN, &
                &       mass_blast_agn(isink)*T2_max  )
        endif
     endif
  end do

  ! Loop over levels
  do ilevel=levelmin,nlevelmax
     ! Computing local volume (important for averaging hydro quantities) 
     dx=0.5D0**ilevel 
     dx_loc=dx*scale
     vol_loc=dx_loc**ndim
     ! Cells center position relative to grid center position
     do ind=1,twotondim  
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        xc(ind,1)=(dble(ix)-0.5D0)*dx
        xc(ind,2)=(dble(iy)-0.5D0)*dx
        xc(ind,3)=(dble(iz)-0.5D0)*dx
     end do

     ! Loop over grids
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

           ! Flag leaf cells
           do i=1,ngrid
              ok(i)=son(ind_cell(i))==0
           end do

           do i=1,ngrid
              if(ok(i))then
                 ! Get gas cell position
                 x=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
                 y=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
                 z=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
                 do isink=1,nsink
                    ! Check if sink is in blast mode
                    if(ok_blast_agn(isink))then
                       ! Check if the cell lies within the sink radius
                       dxx=x-xsink(isink,1)
                       if(dxx> 0.5*scale)then
                          dxx=dxx-scale
                       endif
                       if(dxx<-0.5*scale)then
                          dxx=dxx+scale
                       endif
                       dyy=y-xsink(isink,2)
                       if(dyy> 0.5*scale)then
                          dyy=dyy-scale
                       endif
                       if(dyy<-0.5*scale)then
                          dyy=dyy+scale
                       endif
                       dzz=z-xsink(isink,3)
                       if(dzz> 0.5*scale)then
                          dzz=dzz-scale
                       endif
                       if(dzz<-0.5*scale)then
                          dzz=dzz+scale
                       endif
                       drr=dxx*dxx+dyy*dyy+dzz*dzz
                       
                       if(drr.lt.rmax2)then
                          ! Update the total energy of the gas
                          uold(ind_cell(i),5)=uold(ind_cell(i),5)+p_agn(isink)
                       endif
                    endif
                 end do
              endif
           end do
           
        end do
        ! End loop over cells
     end do
     ! End loop over grids
  end do
  ! End loop over levels

  do isink=1,nsink
     if(ok_blast_agn(isink).and.vol_gas_agn(isink)==0d0)then
        if(ind_blast_agn(isink)>0)then
           uold(ind_blast_agn(isink),5)=uold(ind_blast_agn(isink),5)+p_agn(isink)/vol_blast_agn(isink)
        endif
     endif
  end do

  if(verbose)write(*,*)'Exiting AGN_blast'

end subroutine AGN_blast
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine quenching(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !------------------------------------------------------------------------
  ! This routine selects regions which are eligible for SMBH formation.
  ! It is based on a stellar density threshold and on a stellar velocity
  ! dispersion threshold.
  ! On exit, flag2 array is set to 0 for AGN sites and to 1 otherwise.
  !------------------------------------------------------------------------
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::dx,dx_loc,scale,vol_loc
  real(dp)::str_d,tot_m,ave_u,ave_v,ave_w,sig_u,sig_v,sig_w
  integer::igrid,jgrid,ipart,jpart,next_part,ind_cell,iskip,ind
  integer::i,ig,ip,npart1,npart2,icpu,nx_loc
  real(dp),dimension(1:3)::skip_loc
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Conversion factor from user units to cgs units
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
  vol_loc=dx_loc**ndim

#if NDIM==3
  ! Gather star particles only.

  ! Loop over grids
  do i=1,active(ilevel)%ngrid
     igrid=active(ilevel)%igrid(i)
     ! Number of particles in the grid
     npart1=numbp(igrid)
     npart2=0
     
     ! Reset velocity moments
     str_d=0.0
     tot_m=0.0
     ave_u=0.0
     ave_v=0.0
     ave_w=0.0
     sig_u=0.0
     sig_v=0.0
     sig_w=0.0
     
     ! Count star particles
     if(npart1>0)then
        ipart=headp(igrid)
        ! Loop over particles
        do jpart=1,npart1
           ! Save next particle   <--- Very important !!!
           next_part=nextp(ipart)
           if(idp(ipart).gt.0.and.tp(ipart).ne.0)then
              npart2=npart2+1
              tot_m=tot_m+mp(ipart)
              ave_u=ave_u+mp(ipart)*vp(ipart,1)
              ave_v=ave_v+mp(ipart)*vp(ipart,2)
              ave_w=ave_w+mp(ipart)*vp(ipart,3)
              sig_u=sig_u+mp(ipart)*vp(ipart,1)**2
              sig_v=sig_v+mp(ipart)*vp(ipart,2)**2
              sig_w=sig_w+mp(ipart)*vp(ipart,3)**2
           endif
           ipart=next_part  ! Go to next particle
        end do
     endif
     
     ! Normalize velocity moments
     if(npart2.gt.0)then
        ave_u=ave_u/tot_m
        ave_v=ave_v/tot_m
        ave_w=ave_w/tot_m
        sig_u=sqrt(sig_u/tot_m-ave_u**2)*scale_v/1d5
        sig_v=sqrt(sig_v/tot_m-ave_v**2)*scale_v/1d5
        sig_w=sqrt(sig_w/tot_m-ave_w**2)*scale_v/1d5
        str_d=tot_m/(2**ndim*vol_loc)*scale_nH
     endif
     
     ! Loop over cells
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        ind_cell=iskip+igrid
        ! AGN formation sites
        ! if n_star>0.1 H/cc and v_disp>100 km/s
        if(str_d>0.1.and.MAX(sig_u,sig_v,sig_w)>100.)then
           flag2(ind_cell)=0
        else
           flag2(ind_cell)=1
        end if
     end do
  end do
  ! End loop over grids

#endif

111 format('   Entering quenching for level ',I2)

end subroutine quenching
