!=======================================================================
!=======================================================================
!=======================================================================
!=======================================================================
subroutine output_frame()
  use amr_commons
  use pm_commons
  use hydro_commons
#ifdef RT
  use rt_parameters
  use rt_hydro_commons
#endif
  implicit none
#ifndef WITHOUTMPI
  include "mpif.h"
#endif
  
  integer::dummy_io,info
  integer,parameter::tag=100

  character(len=5) :: istep_str
  character(len=100) :: moviedir, moviecmd, infofile, sinkfile
#ifdef SOLVERmhd
  character(len=100),dimension(0:NVAR+6) :: moviefiles
#else
  character(len=100),dimension(0:NVAR+2) :: moviefiles
#endif
  integer::icell,ncache,iskip,ngrid,nlevelmax_frame
  integer::ilun,nx_loc,ipout,npout,npart_out,ind,ix,iy,iz
  integer::imin,imax,jmin,jmax,ii,jj,kk,ll
  character(LEN=80)::fileloc
  character(LEN=5)::nchar,dummy
  real(dp)::scale,scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::xcen,ycen,zcen,delx,dely,delz
  real(dp)::xleft_frame,xright_frame,yleft_frame,yright_frame,zleft_frame,zright_frame
  real(dp)::xleft,xright,yleft,yright,zleft,zright
  real(dp)::xxleft,xxright,yyleft,yyright,zzleft,zzright
  real(dp)::xpf,ypf,zpf
  real(dp)::dx_frame,dy_frame,dx,dx_loc,dx_min
  real(dp)::dx_cell,dy_cell,dz_cell,dvol
  real(kind=8)::cell_value
  integer ,dimension(1:nvector)::ind_grid,ind_cell
  logical,dimension(1:nvector)::ok
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim)::xx
  real(kind=8),dimension(:,:,:),allocatable::data_frame,data_frame_all
  real(kind=8),dimension(:,:),allocatable::dens,dens_all,vol,vol_all
  real(kind=4),dimension(:,:),allocatable::data_single
  real(kind=8) :: z1,z2,om0in,omLin,hubin,Lbox
  real(kind=8) :: observer(3),thetay,thetaz,theta,phi,temp,ekk
  integer::igrid,jgrid,ipart,jpart,idim,icpu,ilevel,next_part
  integer::i,j,ig,ip,npart1
  integer::nalloc1,nalloc2
  integer::proj_ind,l,nh_temp,nw_temp
  real(kind=4)::ratio

  integer,dimension(1:nvector),save::ind_part,ind_grid_part
  logical::opened

  character(len=1)::temp_string

#ifdef RT
  character(len=100),dimension(1:NGROUPS) :: rt_moviefiles
  real(kind=8),dimension(:,:,:),allocatable::rt_data_frame,rt_data_frame_all
#endif  

  
  nh_temp = nh_frame
  nw_temp = nw_frame

  
 do proj_ind=1,LEN(trim(proj_axis)) 
  opened=.false.

#if NDIM > 1
  if(imov<1)imov=1
  if(imov>imovout)return

  ! Determine the filename, dir, etc
  if(myid==1)write(*,*)'Computing and dumping movie frame'

  call title(imov, istep_str)
  write(temp_string,'(I1)') proj_ind
  moviedir = 'movie'//trim(temp_string)//'/'
  moviecmd = 'mkdir -p '//trim(moviedir)
  if(myid==1) write(*,*) "Writing frame ", istep_str
  if(.not.withoutmkdir) then 
#ifdef NOSYSTEM
     if(myid==1)call PXFMKDIR(TRIM(moviedir),LEN(TRIM(moviedir)),O'755',info)  
#else
     if(myid==1)call system(moviecmd)
#endif
  endif
  
  infofile = trim(moviedir)//'info_'//trim(istep_str)//'.txt'
  if(myid==1)call output_info(infofile)
  
  moviefiles(0) = trim(moviedir)//'temp_'//trim(istep_str)//'.map'
  moviefiles(1) = trim(moviedir)//'dens_'//trim(istep_str)//'.map'
  moviefiles(2) = trim(moviedir)//'vx_'//trim(istep_str)//'.map'
  moviefiles(3) = trim(moviedir)//'vy_'//trim(istep_str)//'.map'
#if NDIM>2
  moviefiles(4) = trim(moviedir)//'vz_'//trim(istep_str)//'.map'
#endif
#if NDIM==2
  moviefiles(4) = trim(moviedir)//'pres_'//trim(istep_str)//'.map'
#endif
#if NDIM>2
  moviefiles(5) = trim(moviedir)//'pres_'//trim(istep_str)//'.map'
#endif
#if NVAR>5
  do ll=6,NVAR
    write(dummy,'(I3.1)') ll
    moviefiles(ll) = trim(moviedir)//'var'//trim(adjustl(dummy))//'_'//trim(istep_str)//'.map'
 end do
#endif
#ifdef SOLVERmhd
  moviefiles(6) = trim(moviedir)//'bxl_'//trim(istep_str)//'.map'
  moviefiles(7) = trim(moviedir)//'byl_'//trim(istep_str)//'.map'
  moviefiles(8) = trim(moviedir)//'bzl_'//trim(istep_str)//'.map'
  moviefiles(NVAR+1) = trim(moviedir)//'bxr_'//trim(istep_str)//'.map'
  moviefiles(NVAR+2) = trim(moviedir)//'byr_'//trim(istep_str)//'.map'
  moviefiles(NVAR+3) = trim(moviedir)//'bzr_'//trim(istep_str)//'.map'
  moviefiles(NVAR+4) = trim(moviedir)//'pmag_'//trim(istep_str)//'.map'
  moviefiles(NVAR+5) = trim(moviedir)//'dm_'//trim(istep_str)//'.map'
  moviefiles(NVAR+6) = trim(moviedir)//'stars_'//trim(istep_str)//'.map'
#else
  moviefiles(NVAR+1) = trim(moviedir)//'dm_'//trim(istep_str)//'.map'
  moviefiles(NVAR+2) = trim(moviedir)//'stars_'//trim(istep_str)//'.map'
#endif

#ifdef RT
  ! Can generate mass weighted averages of cN_i for each group i
  if(rt) then
     do ll=1,NGROUPS
        write(dummy,'(I3.1)') ll
        rt_moviefiles(ll) = trim(moviedir)//'Fp'//trim(adjustl(dummy))//'_'//trim(istep_str)//'.map'
     end do
  endif
#endif

  ! sink filename
  if(sink)then
    sinkfile = trim(moviedir)//'sink_'//trim(istep_str)//'.txt'
    if(myid==1.and.proj_ind==1) call output_sink_csv(sinkfile)
  endif
  
  if(levelmax_frame==0)then
     nlevelmax_frame=nlevelmax
  else if (levelmax_frame.gt.nlevelmax)then
     nlevelmax_frame=nlevelmax
  else
     nlevelmax_frame=levelmax_frame
  endif

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Local constants
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)

  ! Compute frame boundaries
  if(proj_axis(proj_ind:proj_ind).eq.'x')then
    xcen=ycentre_frame(proj_ind*4-3)+ycentre_frame(proj_ind*4-2)*aexp+ycentre_frame(proj_ind*4-1)*aexp**2+ycentre_frame(proj_ind*4)*aexp**3
    ycen=zcentre_frame(proj_ind*4-3)+zcentre_frame(proj_ind*4-2)*aexp+zcentre_frame(proj_ind*4-1)*aexp**2+zcentre_frame(proj_ind*4)*aexp**3
    zcen=xcentre_frame(proj_ind*4-3)+xcentre_frame(proj_ind*4-2)*aexp+xcentre_frame(proj_ind*4-1)*aexp**2+xcentre_frame(proj_ind*4)*aexp**3
    delx=min(2*min(xcen,ycen,zcen,boxlen-xcen,boxlen-ycen,boxlen-zcen),deltay_frame(proj_ind*2-1)+deltay_frame(proj_ind*2)/aexp) !+deltax_frame(3)*aexp**2+deltax_frame(4)*aexp**3  !Essentially comoving or physical
    dely=min(2*min(xcen,ycen,zcen,boxlen-xcen,boxlen-ycen,boxlen-zcen),deltaz_frame(proj_ind*2-1)+deltaz_frame(proj_ind*2)/aexp) !+deltay_frame(3)*aexp**2+deltay_frame(4)*aexp**3
    delz=min(2*min(xcen,ycen,zcen,boxlen-xcen,boxlen-ycen,boxlen-zcen),deltax_frame(proj_ind*2-1)+deltax_frame(proj_ind*2)/aexp) !+deltaz_frame(3)*aexp**2+deltaz_frame(4)*aexp**3
  elseif(proj_axis(proj_ind:proj_ind).eq.'y')then
    xcen=xcentre_frame(proj_ind*4-3)+xcentre_frame(proj_ind*4-2)*aexp+xcentre_frame(proj_ind*4-1)*aexp**2+xcentre_frame(proj_ind*4)*aexp**3
    ycen=zcentre_frame(proj_ind*4-3)+zcentre_frame(proj_ind*4-2)*aexp+zcentre_frame(proj_ind*4-1)*aexp**2+zcentre_frame(proj_ind*4)*aexp**3
    zcen=ycentre_frame(proj_ind*4-3)+ycentre_frame(proj_ind*4-2)*aexp+ycentre_frame(proj_ind*4-1)*aexp**2+ycentre_frame(proj_ind*4)*aexp**3
    delx=min(2*min(xcen,ycen,zcen,boxlen-xcen,boxlen-ycen,boxlen-zcen),deltax_frame(proj_ind*2-1)+deltax_frame(proj_ind*2)/aexp) !+deltax_frame(3)*aexp**2+deltax_frame(4)*aexp**3  !Essentially comoving or physical
    dely=min(2*min(xcen,ycen,zcen,boxlen-xcen,boxlen-ycen,boxlen-zcen),deltaz_frame(proj_ind*2-1)+deltaz_frame(proj_ind*2)/aexp) !+deltay_frame(3)*aexp**2+deltay_frame(4)*aexp**3
    delz=min(2*min(xcen,ycen,zcen,boxlen-xcen,boxlen-ycen,boxlen-zcen),deltay_frame(proj_ind*2-1)+deltay_frame(proj_ind*2)/aexp) !+deltaz_frame(3)*aexp**2+deltaz_frame(4)*aexp**3
  else
    xcen=xcentre_frame(proj_ind*4-3)+xcentre_frame(proj_ind*4-2)*aexp+xcentre_frame(proj_ind*4-1)*aexp**2+xcentre_frame(proj_ind*4)*aexp**3
    ycen=ycentre_frame(proj_ind*4-3)+ycentre_frame(proj_ind*4-2)*aexp+ycentre_frame(proj_ind*4-1)*aexp**2+ycentre_frame(proj_ind*4)*aexp**3
    zcen=zcentre_frame(proj_ind*4-3)+zcentre_frame(proj_ind*4-2)*aexp+zcentre_frame(proj_ind*4-1)*aexp**2+zcentre_frame(proj_ind*4)*aexp**3
    delx=min(2*min(xcen,ycen,zcen,boxlen-xcen,boxlen-ycen,boxlen-zcen),deltax_frame(proj_ind*2-1)+deltax_frame(proj_ind*2)/aexp) !+deltax_frame(3)*aexp**2+deltax_frame(4)*aexp**3  !Essentially comoving or physical
    dely=min(2*min(xcen,ycen,zcen,boxlen-xcen,boxlen-ycen,boxlen-zcen),deltay_frame(proj_ind*2-1)+deltay_frame(proj_ind*2)/aexp) !+deltay_frame(3)*aexp**2+deltay_frame(4)*aexp**3
    delz=min(2*min(xcen,ycen,zcen,boxlen-xcen,boxlen-ycen,boxlen-zcen),deltaz_frame(proj_ind*2-1)+deltaz_frame(proj_ind*2)/aexp) !+deltaz_frame(3)*aexp**2+deltaz_frame(4)*aexp**3
  endif

  ratio = delx/dely
  if(ratio.gt.1)then
    nw_frame=nh_temp*ratio
  else
    nh_frame=nw_temp/ratio
  endif

  xleft_frame=xcen-delx/2.
  xright_frame=xcen+delx/2.
  yleft_frame=ycen-dely/2.
  yright_frame=ycen+dely/2.
  zleft_frame=zcen-delz/2.
  zright_frame=zcen+delz/2.
  
  ! Allocate image
#ifdef SOLVERmhd
  allocate(data_frame(1:nw_frame,1:nh_frame,0:NVAR+6))
#else
  allocate(data_frame(1:nw_frame,1:nh_frame,0:NVAR+2))
#endif
#ifdef RT
  if(rt) then
     allocate(rt_data_frame(1:nw_frame,1:nh_frame,1:NGROUPS))
     rt_data_frame(:,:,:) = 0d0
  endif
#endif
  allocate(dens(1:nw_frame,1:nh_frame))
  allocate(vol(1:nw_frame,1:nh_frame))
  data_frame=0d0
  dens=0d0
  vol=0d0
  dx_frame=delx/dble(nw_frame)
  dy_frame=dely/dble(nh_frame)

  if(hydro) then
     ! Loop over levels
     do ilevel=levelmin,nlevelmax_frame
   
        ! Mesh size at level ilevel in coarse cell units
        dx=0.5D0**ilevel
        
        ! Set position of cell centres relative to grid centre
        do ind=1,twotondim
           iz=(ind-1)/4
           iy=(ind-1-4*iz)/2
           ix=(ind-1-2*iy-4*iz)
           if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
           if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
           if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
        end do
     
        dx_loc=dx*scale
        dx_min=0.5D0**nlevelmax*scale
        ncache=active(ilevel)%ngrid
   
        ! Loop over grids by vector sweeps
        do igrid=1,ncache,nvector
           ngrid=MIN(nvector,ncache-igrid+1)
           do i=1,ngrid
              ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
           end do
           ! Loop over cells
           do ind=1,twotondim
              ! Gather cell indices
              iskip=ncoarse+(ind-1)*ngridmax
              do i=1,ngrid
                 ind_cell(i)=iskip+ind_grid(i)
              end do
              ! Gather cell centre positions
              do idim=1,ndim
                 do i=1,ngrid
                    xx(i,idim)=xg(ind_grid(i),idim)+xc(ind,idim)
                 end do
              end do
              ! Rescale position from code units to user units
              do idim=1,ndim
                 do i=1,ngrid
                    xx(i,idim)=(xx(i,idim)-skip_loc(idim))*scale
                 end do
              end do
              
              ! Check if cell is to be considered
              do i=1,ngrid
                 ok(i)=son(ind_cell(i))==0.or.ilevel==nlevelmax_frame
              end do
   
              do i=1,ngrid
                 if(ok(i))then
                    ! Check if the cell intersect the domain
#if NDIM>2                 
                    if(proj_axis(proj_ind:proj_ind).eq.'x')then
                      xleft=xx(i,2)-dx_loc/2.
                      xright=xx(i,2)+dx_loc/2.
                      yleft=xx(i,3)-dx_loc/2.
                      yright=xx(i,3)+dx_loc/2.
                    elseif(proj_axis(proj_ind:proj_ind).eq.'y')then
                      xleft=xx(i,1)-dx_loc/2.
                      xright=xx(i,1)+dx_loc/2.
                      yleft=xx(i,3)-dx_loc/2.
                      yright=xx(i,3)+dx_loc/2.
                    else
                      xleft=xx(i,1)-dx_loc/2.
                      xright=xx(i,1)+dx_loc/2.
                      yleft=xx(i,2)-dx_loc/2.
                      yright=xx(i,2)+dx_loc/2.
                    endif
                    
                    if(proj_axis(proj_ind:proj_ind).eq.'x')then
                      zleft=xx(i,1)-dx_loc/2.
                      zright=xx(i,1)+dx_loc/2.
                    elseif(proj_axis(proj_ind:proj_ind).eq.'y')then
                      zleft=xx(i,2)-dx_loc/2.
                      zright=xx(i,2)+dx_loc/2.
                    else
                      zleft=xx(i,3)-dx_loc/2.
                      zright=xx(i,3)+dx_loc/2.
                    endif
                    if(    xright.lt.xleft_frame.or.xleft.ge.xright_frame.or.&
                         & yright.lt.yleft_frame.or.yleft.ge.yright_frame.or.&
                         & zright.lt.zleft_frame.or.zleft.ge.zright_frame)cycle
#else
                    xleft=xx(i,1)-dx_loc/2.
                    xright=xx(i,1)+dx_loc/2.
                    yleft=xx(i,2)-dx_loc/2.
                    yright=xx(i,2)+dx_loc/2.
   
                    if(    xright.lt.xleft_frame.or.xleft.ge.xright_frame.or.&
                         & yright.lt.yleft_frame.or.yleft.ge.yright_frame)cycle
#endif
                    ! Compute map indices for the cell
                    if(xleft>xleft_frame)then
                       imin=min(int((xleft-xleft_frame)/dx_frame)+1,nw_frame)
                    else
                       imin=1
                    endif
                    imax=min(int((xright-xleft_frame)/dx_frame)+1,nw_frame)
                    if(yleft>yleft_frame)then
                       jmin=min(int((yleft-yleft_frame)/dy_frame)+1,nh_frame) ! change
                    else
                       jmin=1
                    endif
                    jmax=min(int((yright-yleft_frame)/dy_frame)+1,nh_frame) ! change
                    
                    ! Fill up map with projected mass
#if NDIM>2                 
                    dz_cell=min(zright_frame,zright)-max(zleft_frame,zleft) ! change
#endif
                    do ii=imin,imax
                       xxleft=xleft_frame+dble(ii-1)*dx_frame
                       xxright=xxleft+dx_frame
                       dx_cell=min(xxright,xright)-max(xxleft,xleft)
                       do jj=jmin,jmax
                          yyleft=yleft_frame+dble(jj-1)*dy_frame
                          yyright=yyleft+dy_frame
                          dy_cell=min(yyright,yright)-max(yyleft,yleft)
                          ! Intersection volume
                          dvol=dx_cell*dy_cell
#if NDIM>2                 
                          dvol=dvol*dz_cell
#endif
                          dens(ii,jj)=dens(ii,jj)+dvol*max(uold(ind_cell(i),1),smallr)
                          vol(ii,jj)=vol(ii,jj)+dvol
                          
                          data_frame(ii,jj,1)=data_frame(ii,jj,1)+dvol*max(uold(ind_cell(i),1),smallr)**2
#ifdef SOLVERmhd
                          do kk=2,NVAR+3
#else                       
                          do kk=2,NVAR
#endif
                             if(movie_vars(kk).eq.1) data_frame(ii,jj,kk)=data_frame(ii,jj,kk)+dvol*uold(ind_cell(i),kk)
                          end do
   
#ifdef RT
                          if(rt) then
                             do kk=1,NGROUPS
                                if(rt_movie_vars(kk).eq.1) then
                                   rt_data_frame(ii,jj,kk) = rt_data_frame(ii,jj,kk) &
                                        + dvol * rtuold(ind_cell(i), 1+(kk-1)*(ndim+1)) * rt_c_cgs &
                                               * max(uold(ind_cell(i),1),smallr) ! mass-weighted
                                endif
                             end do
                          endif
#endif
   
                          
                          if (movie_vars(0).eq.1)then
                            !Get temperature
                            ekk=0.0d0
                            do idim=1,3
                               ekk=ekk+0.5*uold(ind_cell(i),idim+1)**2/max(uold(ind_cell(i),1),smallr)
                            enddo
                            temp=(gamma-1.0)*(uold(ind_cell(i),5)-ekk) !pressure
                            temp=max(temp/max(uold(ind_cell(i),1),smallr),smallc**2)*scale_T2 !temperature in K
   
                            data_frame(ii,jj,0)=data_frame(ii,jj,0)+dvol*max(uold(ind_cell(i),1),smallr)*temp !mass weighted temperature
                          end if
   
#ifdef SOLVERmhd
                          if (movie_vars(NVAR+4).eq.1)then
                                  data_frame(ii,jj,NVAR+4)=data_frame(ii,jj,NVAR+4)+ dvol*0.125*(&
                                      uold(ind_cell(i),6)**2 + uold(ind_cell(i),7)**2 + uold(ind_cell(i),8)**2 &
                                      + uold(ind_cell(i),NVAR+1)**2 + uold(ind_cell(i),NVAR+2)**2 + uold(ind_cell(i),NVAR+3)**2)
                          end if
#endif
   
                       end do
                    end do
                 end if
              end do
   
           end do
           ! End loop over cells
   
        end do
        ! End loop over grids
     end do
  ! End loop over levels
  end if

  ! Loop over particles
  do j=1,npartmax
#if NDIM>2                 
     if(proj_axis(proj_ind:proj_ind).eq.'x')then
       xpf  = xp(j,2)
       ypf  = xp(j,3)
     elseif(proj_axis(proj_ind:proj_ind).eq.'y')then
       xpf  = xp(j,1)
       ypf  = xp(j,3)
     else
       xpf  = xp(j,1)
       ypf  = xp(j,2)
     endif
     
     if(proj_axis(proj_ind:proj_ind).eq.'x')then
       zpf  = xp(j,1)
     elseif(proj_axis(proj_ind:proj_ind).eq.'y')then
       zpf  = xp(j,2)
     else
       zpf  = xp(j,3)
     endif
     if(    xpf.lt.xleft_frame.or.xpf.ge.xright_frame.or.&
          & ypf.lt.yleft_frame.or.ypf.ge.yright_frame.or.&
          & zpf.lt.zleft_frame.or.zpf.ge.zright_frame)cycle
#else
     xpf  = xp(j,1)
     ypf  = xp(j,2)
     
     if(    xpf.lt.xleft_frame.or.xpf.ge.xright_frame.or.&
          & ypf.lt.yleft_frame.or.ypf.ge.yright_frame)cycle
#endif
     ! Compute map indices for the cell
     ii = min(int((xpf-xleft_frame)/dx_frame)+1,nw_frame)
     jj = min(int((ypf-yleft_frame)/dy_frame)+1,nh_frame)
     
     ! Fill up map with projected mass
#ifdef SOLVERmhd
     if(star) then
        if(tp(j).eq.0.) then
           if(mass_cut_refine>0.0) then
              if(mp(j)<mass_cut_refine) data_frame(ii,jj,NVAR+5)=data_frame(ii,jj,NVAR+5)+mp(j)
           else
              data_frame(ii,jj,NVAR+5)=data_frame(ii,jj,NVAR+5)+mp(j)
           endif
        else
           data_frame(ii,jj,NVAR+6)=data_frame(ii,jj,NVAR+6)+mp(j)
        endif
     else
        if(mass_cut_refine>0.0) then
           if(mp(j)<mass_cut_refine) data_frame(ii,jj,NVAR+5)=data_frame(ii,jj,NVAR+5)+mp(j)
        else
           data_frame(ii,jj,NVAR+5)=data_frame(ii,jj,NVAR+5)+mp(j)
        endif
     endif
#else
     if(star) then
        if(tp(j).eq.0.) then
           if(mass_cut_refine>0.0) then
              if(mp(j)<mass_cut_refine) data_frame(ii,jj,NVAR+1)=data_frame(ii,jj,NVAR+1)+mp(j)
           else
              data_frame(ii,jj,NVAR+1)=data_frame(ii,jj,NVAR+1)+mp(j)
           endif
        else
           data_frame(ii,jj,NVAR+2)=data_frame(ii,jj,NVAR+2)+mp(j)
        endif
     else
        if(mass_cut_refine>0.0) then
           if(mp(j)<mass_cut_refine) data_frame(ii,jj,NVAR+1)=data_frame(ii,jj,NVAR+1)+mp(j)
        else
           data_frame(ii,jj,NVAR+1)=data_frame(ii,jj,NVAR+1)+mp(j)
        endif
     endif
#endif
  end do
  ! End loop over particles

#ifndef WITHOUTMPI
#ifdef SOLVERmhd
  allocate(data_frame_all(1:nw_frame,1:nh_frame,0:NVAR+6))
  call MPI_ALLREDUCE(data_frame,data_frame_all,nw_frame*nh_frame*(NVAR+6+1),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#else
  allocate(data_frame_all(1:nw_frame,1:nh_frame,0:NVAR+2))
  call MPI_ALLREDUCE(data_frame,data_frame_all,nw_frame*nh_frame*(NVAR+2+1),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
  allocate(dens_all(1:nw_frame,1:nh_frame))
  call MPI_ALLREDUCE(dens,dens_all,nw_frame*nh_frame,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  allocate(vol_all(1:nw_frame,1:nh_frame))
  call MPI_ALLREDUCE(vol,vol_all,nw_frame*nh_frame,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  data_frame=data_frame_all
  dens=dens_all
  vol=vol_all
  deallocate(data_frame_all)
  deallocate(dens_all)
  deallocate(vol_all)
#ifdef RT
  if(rt) then
     allocate(rt_data_frame_all(1:nw_frame,1:nh_frame,1:NGROUPS))
     rt_data_frame_all(:,:,:)=0d0
     call MPI_ALLREDUCE(rt_data_frame,rt_data_frame_all        &
          ,nw_frame*nh_frame*NGROUPS,MPI_DOUBLE_PRECISION      &
          ,MPI_SUM,MPI_COMM_WORLD,info)
     rt_data_frame=rt_data_frame_all
     deallocate(rt_data_frame_all)
  endif
#endif
#endif
  ! Convert into mass weighted                                                                                                         
  do ii=1,nw_frame
    do jj=1,nh_frame
#ifdef SOLVERmhd
      do kk=0,5
        if(movie_vars(kk).eq.1) data_frame(ii,jj,kk)=data_frame(ii,jj,kk)/dens(ii,jj)
      end do
      do kk=6,8
        if(movie_vars(kk).eq.1) data_frame(ii,jj,kk)=data_frame(ii,jj,kk)/vol(ii,jj)
      end do
      do kk=9,NVAR
        if(movie_vars(kk).eq.1) data_frame(ii,jj,kk)=data_frame(ii,jj,kk)/dens(ii,jj)
      end do
      do kk=NVAR+1,NVAR+4
        if(movie_vars(kk).eq.1) data_frame(ii,jj,kk)=data_frame(ii,jj,kk)/vol(ii,jj)
      end do
#else
      do kk=0,NVAR
        if(movie_vars(kk).eq.1) data_frame(ii,jj,kk)=data_frame(ii,jj,kk)/dens(ii,jj)
      end do
#endif
#ifdef RT
      if(rt) then
         do kk=1,NGROUPS
            if(rt_movie_vars(kk).eq.1) &
                 rt_data_frame(ii,jj,kk)=rt_data_frame(ii,jj,kk)/dens(ii,jj)
         end do
      endif
#endif

    end do
  end do
  deallocate(dens)
  deallocate(vol)

  if(myid==1)then
     ilun=10
     allocate(data_single(1:nw_frame,1:nh_frame))
     ! Output mass weighted density
#ifdef SOLVERmhd
     do kk=0, NVAR+6
#else
     do kk=0, NVAR+2
#endif
       if (movie_vars(kk).eq.1)then
         open(ilun,file=TRIM(moviefiles(kk)),form='unformatted')
         data_single=data_frame(:,:,kk)
         rewind(ilun)  
         if(tendmov>0)then
            write(ilun)t,delx,dely,delz
         else
            write(ilun)aexp,delx,dely,delz
         endif
         write(ilun)nw_frame,nh_frame
         write(ilun)data_single
         close(ilun)
       end if
     end do

#ifdef RT
      if(rt) then
         do kk=1, NGROUPS
            if (rt_movie_vars(kk).eq.1) then
               open(ilun,file=TRIM(rt_moviefiles(kk)),form='unformatted')
               data_single(:,:)=0.
               data_single=rt_data_frame(:,:,kk)
               rewind(ilun)  
               if(tendmov>0)then
                  write(ilun)t,delx,dely,delz
               else
                  write(ilun)aexp,delx,dely,delz
               endif
               write(ilun)nw_frame,nh_frame
               write(ilun)data_single
               close(ilun)
            end if
         end do
      endif
#endif
     
     deallocate(data_single)
  endif

  deallocate(data_frame)
#ifdef RT
  if(rt) deallocate(rt_data_frame)
#endif
#endif
  ! Update counter
  if(proj_ind.eq.len(trim(proj_axis))) then 
     ! Increase counter and skip frames if timestep is large
     imov=imov+1
     do while((amovout(imov)<aexp.or.tmovout(imov)<t).and.(imov.lt.imovout))
        imov=imov+1
     end do
  endif

  nw_frame = nw_temp
  nh_frame = nh_temp
 enddo
end subroutine output_frame

subroutine set_movie_vars()
  use amr_commons
  ! This routine sets the movie vars from textual form
  integer::ll
  character(LEN=5)::dummy

  if(ANY(movie_vars_txt=='temp ')) movie_vars(0)=1
  if(ANY(movie_vars_txt=='dens ')) movie_vars(1)=1
  if(ANY(movie_vars_txt=='vx   ')) movie_vars(2)=1
  if(ANY(movie_vars_txt=='vy   ')) movie_vars(3)=1
#if NDIM>2
  if(ANY(movie_vars_txt=='vz   ')) movie_vars(4)=1
#endif
#if NDIM==2
  if(ANY(movie_vars_txt=='pres ')) movie_vars(4)=1
#endif
#if NDIM>2
  if(ANY(movie_vars_txt=='pres ')) movie_vars(5)=1
#endif
#if NVAR>5
  do ll=6,NVAR
    write(dummy,'(I3.1)') ll
    if(ANY(movie_vars_txt=='var'//trim(adjustl(dummy))//' ')) movie_vars(ll)=1
 end do
#endif
#ifdef SOLVERmhd
  if(ANY(movie_vars_txt=='bxl  ')) movie_vars(6)=1
  if(ANY(movie_vars_txt=='byl  ')) movie_vars(7)=1
  if(ANY(movie_vars_txt=='bzl  ')) movie_vars(8)=1
  if(ANY(movie_vars_txt=='bxr  ')) movie_vars(NVAR+1)=1
  if(ANY(movie_vars_txt=='byr  ')) movie_vars(NVAR+2)=1
  if(ANY(movie_vars_txt=='bzr  ')) movie_vars(NVAR+3)=1
  if(ANY(movie_vars_txt=='pmag ')) movie_vars(NVAR+4)=1
  if(ANY(movie_vars_txt=='dm   ')) movie_vars(NVAR+5)=1
  if(ANY(movie_vars_txt=='stars')) movie_vars(NVAR+6)=1
#else
  if(ANY(movie_vars_txt=='dm   ')) movie_vars(NVAR+1)=1
  if(ANY(movie_vars_txt=='stars')) movie_vars(NVAR+2)=1
#endif
end subroutine set_movie_vars
