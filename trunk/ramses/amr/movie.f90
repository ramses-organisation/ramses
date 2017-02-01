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
  
  integer::dummy_io,info,ierr,iframe
  integer,parameter::tag=100

  character(len=5) :: istep_str
  character(len=100) :: moviedir, moviecmd, infofile, sinkfile
#ifdef SOLVERmhd
  character(len=100),dimension(0:NVAR+6) :: moviefiles
#else
  character(len=100),dimension(0:NVAR+2) :: moviefiles
#endif
  integer::icell,ncache,iskip,irad,ngrid,nlevelmax_frame
  integer::nframes,rt_nframes,imap,ipart_start
  integer::ilun,nx_loc,ipout,npout,npart_out,ind,ix,iy,iz
  integer::imin,imax,jmin,jmax,ii,jj,kk,ll
  character(LEN=80)::fileloc
  character(LEN=5)::nchar,dummy
  real(dp)::scale,scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::xcen,ycen,zcen,delx,dely,delz
  real(dp)::xtmp,ytmp,ztmp,smooth,theta_cam,phi_cam,alpha,beta,smooth_theta,fov_camera,dist_cam
  real(dp)::xleft_frame,xright_frame,yleft_frame,yright_frame,zleft_frame,zright_frame,rr
  real(dp)::xleft,xright,yleft,yright,zleft,zright,xcentre,ycentre,zcentre
  real(dp)::xxleft,xxright,yyleft,yyright,zzleft,zzright,xxcentre,yycentre,zzcentre
  real(dp)::xpf,ypf,zpf
  real(dp)::dx_frame,dy_frame,dx,dx_loc,dx_min,pers_corr
  real(dp)::dx_cell,dy_cell,dz_cell,dvol,dx_proj,weight
  real(kind=8)::cell_value
  integer ,dimension(1:nvector)::ind_grid,ind_cell
  logical,dimension(1:nvector)::ok
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim)::xx
  real(dp),dimension(1:nvector,1:ndim)::xx2
  real(kind=8),dimension(:,:,:),allocatable::data_frame
  real(kind=8),dimension(:,:),allocatable::weights
  real(kind=8),dimension(:),allocatable::data_single,data_single_all
  real(kind=8) :: z1,z2,om0in,omLin,hubin,Lbox
  real(kind=8) :: observer(3),thetay,thetaz,theta,phi,temp,e,uvar
  real(kind=8) :: pi=3.14159265359
  real(dp),dimension(8)::xcube,ycube,zcube
  integer::igrid,jgrid,ipart,jpart,idim,icpu,ilevel,next_part,icube,iline
  integer::i,j,ig,ip,npart1
  integer::nalloc1,nalloc2
  integer::proj_ind,l,nh_temp,nw_temp
  real(dp)::minx,maxx,miny,maxy,minz,maxz,xpc,ypc,zpc,d1,d2,d3,d4,l1,l2,l3,l4
  integer,dimension(1:nvector),save::ind_part,ind_grid_part
  logical::opened,cube_face
  character(len=1)::temp_string
  integer,dimension(6,8)::lind = reshape((/1, 2, 3, 4, 1, 3, 2, 4,    &
                                           5, 6, 7, 8, 5, 7, 6, 8,    &
                                           1, 5, 2, 6, 1, 2, 5, 6,    &
                                           3, 7, 4, 8, 3, 4, 7, 8,    &
                                           1, 3, 5, 7, 1, 5, 3, 7,    &
                                           2, 4, 6, 8, 2, 6, 4, 8 /)  &
                                           ,shape(lind),order=(/2,1/))
#ifdef RT
  character(len=100),dimension(1:NGROUPS) :: rt_moviefiles
  real(kind=8),dimension(:,:,:),allocatable::rt_data_frame
#endif  

  nh_temp    = nh_frame
  nw_temp    = nw_frame

  ! Only one projection available in 2D
  if((ndim.eq.2).and.(trim(proj_axis).ne.'z')) proj_axis = 'z'
  
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
  if(.not.withoutmkdir) then 
#ifdef NOSYSTEM
     if(myid==1)call PXFMKDIR(TRIM(moviedir),LEN(TRIM(moviedir)),O'755',info)  
#else
     if(myid==1)then
        call EXECUTE_COMMAND_LINE(moviecmd,exitstat=ierr,wait=.true.)
     endif
#ifndef WITHOUTMPI
     call MPI_BCAST(ierr,1,MPI_INTEGER,0,MPI_COMM_WORLD,info)
     if(ierr.ne.0 .and. ierr.ne.127)then
        write(*,*) 'Error - Could not create ',trim(moviedir)
        call MPI_ABORT(MPI_COMM_WORLD,1,info)
        stop
     endif
#endif
#endif
  endif
  
  infofile = trim(moviedir)//'info_'//trim(istep_str)//'.txt'
  if(myid==1)call output_info(infofile)
#ifndef WITHOUTMPI
  call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
  
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

  nframes = 0
#ifdef SOLVERmhd
  do kk=0,NVAR+6
#else                       
  do kk=0,NVAR+2
#endif
     if(movie_vars(kk).eq.1) nframes = nframes+1
  enddo
  rt_nframes = 0
#ifdef RT
  if(rt)then
     do kk=1,NGROUPS
        if(rt_movie_vars(kk).eq.1) rt_nframes = rt_nframes+1
     enddo
  endif
#endif
  if(nframes+rt_nframes==0) continue

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Local constants
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)

  if(xcentre_frame(proj_ind*4-3).eq.0d0) xcentre_frame(proj_ind*4-3) = boxlen/2d0
  if(ycentre_frame(proj_ind*4-3).eq.0d0) ycentre_frame(proj_ind*4-3) = boxlen/2d0
  if(zcentre_frame(proj_ind*4-3).eq.0d0) zcentre_frame(proj_ind*4-3) = boxlen/2d0
  ! Compute frame boundaries
  if(proj_axis(proj_ind:proj_ind).eq.'x')then
    xcen=ycentre_frame(proj_ind*4-3)+ycentre_frame(proj_ind*4-2)*aexp+ycentre_frame(proj_ind*4-1)*aexp**2+ycentre_frame(proj_ind*4)*aexp**3
    ycen=zcentre_frame(proj_ind*4-3)+zcentre_frame(proj_ind*4-2)*aexp+zcentre_frame(proj_ind*4-1)*aexp**2+zcentre_frame(proj_ind*4)*aexp**3
    zcen=xcentre_frame(proj_ind*4-3)+xcentre_frame(proj_ind*4-2)*aexp+xcentre_frame(proj_ind*4-1)*aexp**2+xcentre_frame(proj_ind*4)*aexp**3
  elseif(proj_axis(proj_ind:proj_ind).eq.'y')then
    xcen=xcentre_frame(proj_ind*4-3)+xcentre_frame(proj_ind*4-2)*aexp+xcentre_frame(proj_ind*4-1)*aexp**2+xcentre_frame(proj_ind*4)*aexp**3
    ycen=zcentre_frame(proj_ind*4-3)+zcentre_frame(proj_ind*4-2)*aexp+zcentre_frame(proj_ind*4-1)*aexp**2+zcentre_frame(proj_ind*4)*aexp**3
    zcen=ycentre_frame(proj_ind*4-3)+ycentre_frame(proj_ind*4-2)*aexp+ycentre_frame(proj_ind*4-1)*aexp**2+ycentre_frame(proj_ind*4)*aexp**3
  else
    xcen=xcentre_frame(proj_ind*4-3)+xcentre_frame(proj_ind*4-2)*aexp+xcentre_frame(proj_ind*4-1)*aexp**2+xcentre_frame(proj_ind*4)*aexp**3
    ycen=ycentre_frame(proj_ind*4-3)+ycentre_frame(proj_ind*4-2)*aexp+ycentre_frame(proj_ind*4-1)*aexp**2+ycentre_frame(proj_ind*4)*aexp**3
    zcen=zcentre_frame(proj_ind*4-3)+zcentre_frame(proj_ind*4-2)*aexp+zcentre_frame(proj_ind*4-1)*aexp**2+zcentre_frame(proj_ind*4)*aexp**3
  endif
  if(deltax_frame(proj_ind*2-1).eq.0d0 .and. deltay_frame(proj_ind*2-1).gt.0d0)then
     deltax_frame(proj_ind*2-1)=deltay_frame(proj_ind*2-1)*float(nw_frame)/float(nh_frame)
  endif
  if(deltay_frame(proj_ind*2-1).eq.0d0 .and. deltax_frame(proj_ind*2-1).gt.0d0)then
     deltay_frame(proj_ind*2-1)=deltax_frame(proj_ind*2-1)*float(nh_frame)/float(nw_frame)
  endif
  if(deltaz_frame(proj_ind*2-1).eq.0d0)then
     deltaz_frame(proj_ind*2-1)=boxlen
  endif
  delx=deltax_frame(proj_ind*2-1)+deltax_frame(proj_ind*2)/aexp
  dely=deltay_frame(proj_ind*2-1)+deltay_frame(proj_ind*2)/aexp
  delz=deltaz_frame(proj_ind*2-1)+deltaz_frame(proj_ind*2)/aexp
  if(dist_camera(proj_ind).le.0D0) dist_camera(proj_ind) = boxlen
   
  ! Camera properties
  if(cosmo) then
     if(tend_theta_camera(proj_ind).le.0d0) tend_theta_camera(proj_ind) = aendmov
     if(tend_phi_camera(proj_ind).le.0d0) tend_phi_camera(proj_ind) = aendmov
     theta_cam  = theta_camera(proj_ind)*pi/180.                                                                                 &
                +min(max(aexp-tstart_theta_camera(proj_ind),0d0),tend_theta_camera(proj_ind))*dtheta_camera(proj_ind)*pi/180./(aendmov-astartmov)
     phi_cam    = phi_camera(proj_ind)*pi/180.                                                                                   &
                +min(max(aexp-tstart_theta_camera(proj_ind),0d0),tend_phi_camera(proj_ind))*dphi_camera(proj_ind)*pi/180./(aendmov-astartmov)
     dist_cam   = dist_camera(proj_ind)+min(max(aexp-tstart_theta_camera(proj_ind),0d0),tend_theta_camera(proj_ind))*ddist_camera(proj_ind)/(aendmov-astartmov)
  else
     if(tend_theta_camera(proj_ind).le.0d0) tend_theta_camera(proj_ind) = tendmov
     if(tend_phi_camera(proj_ind).le.0d0) tend_phi_camera(proj_ind) = tendmov
     theta_cam  = theta_camera(proj_ind)*pi/180.                                                                                 &
                +min(max(t-tstart_theta_camera(proj_ind),0d0),tend_theta_camera(proj_ind))*dtheta_camera(proj_ind)*pi/180./(tendmov-tstartmov)
     phi_cam    = phi_camera(proj_ind)*pi/180.                                                                                   &
                +min(max(t-tstart_phi_camera(proj_ind),0d0),tend_phi_camera(proj_ind))*dphi_camera(proj_ind)*pi/180./(tendmov-tstartmov)
     dist_cam   = dist_camera(proj_ind)+min(max(t-tstart_theta_camera(proj_ind),0d0),tend_theta_camera(proj_ind))*ddist_camera(proj_ind)/(tendmov-tstartmov)
  endif
  
  if((focal_camera(proj_ind).le.0D0).or.(focal_camera(proj_ind).gt.dist_camera(proj_ind))) focal_camera(proj_ind) = dist_cam
  fov_camera = atan((delx/2d0)/focal_camera(proj_ind))
#if NDIM>2                 
  if(myid==1) write(*,'(5A,F6.1,A,F6.1,A,F6.3,A,F4.1)') ' Writing frame ', istep_str,' los=',proj_axis(proj_ind:proj_ind),   &
  &                                              ' theta=',theta_cam*180./pi,' phi=',phi_cam*180./pi,' d=',dist_cam,' fov=',fov_camera*180./pi
#else
  if(myid==1) write(*,'(3A,F6.1)') " Writing frame ", istep_str,' theta=',theta_cam*180./pi
#endif
  ! Frame boundaries
  xleft_frame  = xcen-delx/2.
  xright_frame = xcen+delx/2.
  yleft_frame  = ycen-dely/2.
  yright_frame = ycen+dely/2.
  zleft_frame  = zcen-delz/2.
  zright_frame = zcen+delz/2.

  ! No cubic shader for 2D simulations
  if((ndim.eq.2).and.(shader_frame(proj_ind).eq.'cube')) shader_frame(proj_ind) = 'square'

  ! Allocate image
  allocate(data_frame(1:nw_frame,1:nh_frame,1:nframes),stat=ierr)
if(ierr .ne. 0)then
   write(*,*) 'Error - Movie frame allocation failed'
#ifndef WITHOUTMPI
   call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
   stop
#endif
endif
#ifdef RT
  if(rt) then
     allocate(rt_data_frame(1:nw_frame,1:nh_frame,1:rt_nframes),stat=ierr)
     if(ierr .ne. 0)then
        write(*,*) 'Error - Movie frame allocation failed'
#ifndef WITHOUTMPI
        call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
#else
        stop
#endif
     endif
     rt_data_frame(:,:,:) = 0d0
  endif
#endif
  allocate(weights(1:nw_frame,1:nh_frame),stat=ierr)
  if(ierr .ne. 0)then
     write(*,*) 'Error - Movie frame allocation failed'
#ifndef WITHOUTMPI
     call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
     stop
#endif
  endif
  if(ierr .ne. 0)then
     write(*,*) 'Error - Movie frame allocation failed'
#ifndef WITHOUTMPI
     call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
     stop
#endif
  endif
  if(method_frame(proj_ind).eq.'min')then
     data_frame(:,:,:) = 1e-3*huge(0.0)
  elseif(method_frame(proj_ind).eq.'max')then
     data_frame(:,:,:) = -1e-3*huge(0.0)
  else
     data_frame(:,:,:) = 0.0
  endif
  weights(:,:) = 0d0
  dx_frame = delx/dble(nw_frame)
  dy_frame = dely/dble(nh_frame)

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

        dx_proj = (dx_loc/2.0)*smooth_frame(proj_ind)
#if NDIM>2                 
        if((shader_frame(proj_ind).eq.'cube').and.(.not.perspective_camera(proj_ind)))then
           xcube = (/-dx_proj,-dx_proj,-dx_proj,-dx_proj, dx_proj, dx_proj, dx_proj, dx_proj/)
           ycube = (/-dx_proj,-dx_proj, dx_proj, dx_proj,-dx_proj,-dx_proj, dx_proj, dx_proj/)
           zcube = (/-dx_proj, dx_proj,-dx_proj, dx_proj,-dx_proj, dx_proj,-dx_proj, dx_proj/)
           do icube=1,8
              xtmp         = cos(theta_cam)*xcube(icube)+sin(theta_cam)*ycube(icube)
              ytmp         = cos(theta_cam)*ycube(icube)-sin(theta_cam)*xcube(icube)
              xcube(icube) = xtmp
              ycube(icube) = ytmp
              ytmp         = cos(phi_cam)*ycube(icube)+sin(phi_cam)*zcube(icube)
              ztmp         = cos(phi_cam)*zcube(icube)-sin(phi_cam)*ycube(icube)
              ycube(icube) = ytmp
              zcube(icube) = ztmp
           enddo
           minx = minval(xcube)
           maxx = maxval(xcube)
           miny = minval(ycube)
           maxy = maxval(ycube)
           minz = minval(zcube)
           maxz = maxval(zcube)
        endif
#endif
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
                 if(ivar_refine>0.and.zoom_only_frame(proj_ind)) then
                   ok(i)=ok(i).and. &
                      & (uold(ind_cell(i),ivar_refine)/uold(ind_cell(i),1) &
                      & > var_cut_refine)
                 endif
                 if((ok(i)).and.(ivar_frame(proj_ind)>0).and.(ivar_frame(proj_ind)<=nvar))then
                    uvar = uold(ind_cell(i),ivar_frame(proj_ind))
                    ! Scale temperature to K
                    if(ivar_frame(proj_ind)==ndim+2)then
                       e = 0.0d0
                       do idim=1,ndim
                          e = e+0.5*uold(ind_cell(i),idim+1)**2/uold(ind_cell(i),1)
                       enddo
#if NENER>0
                       do irad=0,nener-1
                          e = e+uold(ind_cell(i),inener+irad)
                       enddo
#endif
#ifdef SOLVERmhd
                       do idim=1,ndim 
                          e = e+0.125d0*(uold(ind_cell(i),idim+ndim+2)+uold(ind_cell(i),idim+nvar))**2
                       enddo
#endif
                       ! Pressure
                       uvar = (gamma-1.0)*(uold(ind_cell(i),ivar_frame(proj_ind))-e)*scale_T2
                    endif
                    ! Switch to primitive variables
                    if(ivar_frame(proj_ind)>1) uvar = uvar/uold(ind_cell(i),1)
                    ! Scale density to cm**-3
                    if(ivar_frame(proj_ind)==1) uvar = uvar*scale_nH
                    ! Scale velocities to km/s
                    if(ivar_frame(proj_ind)>1.and.ivar_frame(proj_ind)<ndim+2) uvar = uvar*scale_v/1e5
                    ok(i) = ok(i).and.(uvar.ge.varmin_frame(proj_ind))
                    ok(i) = ok(i).and.(uvar.le.varmax_frame(proj_ind))
                 endif
              end do
   
              do i=1,ngrid
                 if(ok(i))then
#if NDIM>2                 
                    ! Centering
                    xx(i,1) = xx(i,1)-xcen
                    xx(i,2) = xx(i,2)-ycen
                    xx(i,3) = xx(i,3)-zcen
                    ! Rotating
                    xtmp    = cos(theta_cam)*xx(i,1)+sin(theta_cam)*xx(i,2)
                    ytmp    = cos(theta_cam)*xx(i,2)-sin(theta_cam)*xx(i,1)
                    xx(i,1) = xtmp
                    xx(i,2) = ytmp
                    ytmp    = cos(phi_cam)*xx(i,2)+sin(phi_cam)*xx(i,3)
                    ztmp    = cos(phi_cam)*xx(i,3)-sin(phi_cam)*xx(i,2)
                    xx(i,2) = ytmp
                    xx(i,3) = ztmp
                    ! Perspective correction factor
                    pers_corr = 1.0 
                    if(proj_axis(proj_ind:proj_ind).eq.'x')then
                      if(dist_cam-xx(i,1).lt.0d0) continue
                      if(perspective_camera(proj_ind))then
                         alpha  = atan(xx(i,2)/(dist_cam-xx(i,1)))
                         beta   = atan(xx(i,3)/(dist_cam-xx(i,1)))
                         if(abs(alpha)/2d0.gt.fov_camera) continue
                         if(abs(beta)/2d0.gt.fov_camera) continue
                         pers_corr = focal_camera(proj_ind)/(dist_cam-xx(i,1))
                         xx(i,2)   = xx(i,2)*pers_corr
                         xx(i,3)   = xx(i,3)*pers_corr
                         dx_proj   = (dx_loc/2.0)*pers_corr*smooth_frame(proj_ind)
                      endif
                      xcentre = xx(i,2)+ycen
                      ycentre = xx(i,3)+zcen
                      zcentre = xx(i,1)+xcen
                    elseif(proj_axis(proj_ind:proj_ind).eq.'y')then
                      if(dist_cam-xx(i,2).lt.0d0) continue
                      if(perspective_camera(proj_ind))then
                         alpha  = atan(xx(i,1)/(dist_cam-xx(i,2)))
                         beta   = atan(xx(i,3)/(dist_cam-xx(i,2)))
                         if(abs(alpha)/2d0.gt.fov_camera) continue
                         if(abs(beta)/2d0.gt.fov_camera) continue
                         pers_corr = focal_camera(proj_ind)/(dist_cam-xx(i,2))
                         xx(i,1)   = xx(i,1)*pers_corr
                         xx(i,3)   = xx(i,3)*pers_corr
                         dx_proj   = (dx_loc/2.0)*pers_corr*smooth_frame(proj_ind)
                      endif
                      xcentre = xx(i,1)+xcen
                      ycentre = xx(i,3)+zcen
                      zcentre = xx(i,2)+ycen
                    else
                      if(dist_cam-xx(i,3).lt.0d0) continue
                      if(perspective_camera(proj_ind))then
                         alpha  = atan(xx(i,1)/(dist_cam-xx(i,3)))
                         beta   = atan(xx(i,2)/(dist_cam-xx(i,3)))
                         if(abs(alpha)/2d0.gt.fov_camera) continue
                         if(abs(beta)/2d0.gt.fov_camera) continue
                         pers_corr = focal_camera(proj_ind)/(dist_cam-xx(i,3))
                         xx(i,1)   = xx(i,1)*pers_corr
                         xx(i,2)   = xx(i,2)*pers_corr
                         dx_proj   = (dx_loc/2.0)*pers_corr*smooth_frame(proj_ind)
                      endif
                      xcentre = xx(i,1)+xcen
                      ycentre = xx(i,2)+ycen
                      zcentre = xx(i,3)+zcen
                    endif
                    ! Rotating the cube shader
                    if(shader_frame(proj_ind).eq.'cube'.and.perspective_camera(proj_ind))then
                       xcube = (/-dx_proj,-dx_proj,-dx_proj,-dx_proj, dx_proj, dx_proj, dx_proj, dx_proj/)
                       ycube = (/-dx_proj,-dx_proj, dx_proj, dx_proj,-dx_proj,-dx_proj, dx_proj, dx_proj/)
                       zcube = (/-dx_proj, dx_proj,-dx_proj, dx_proj,-dx_proj, dx_proj,-dx_proj, dx_proj/)
                       do icube=1,8
                          xtmp         = cos(theta_cam)*xcube(icube)+sin(theta_cam)*ycube(icube)
                          ytmp         = cos(theta_cam)*ycube(icube)-sin(theta_cam)*xcube(icube)
                          xcube(icube) = xtmp
                          ycube(icube) = ytmp
                          ytmp         = cos(phi_cam)*ycube(icube)+sin(phi_cam)*zcube(icube)
                          ztmp         = cos(phi_cam)*zcube(icube)-sin(phi_cam)*ycube(icube)
                          ycube(icube) = ytmp
                          zcube(icube) = ztmp
                          ! Additional coordinate dependent rotation for perspective effect
                          xtmp         = cos(alpha)*xcube(icube)+sin(alpha)*zcube(icube)
                          ztmp         = cos(alpha)*zcube(icube)-sin(alpha)*xcube(icube)
                          xcube(icube) = xtmp
                          zcube(icube) = ztmp
                          ytmp         = cos(beta)*ycube(icube)+sin(beta)*zcube(icube)
                          ztmp         = cos(beta)*zcube(icube)-sin(beta)*ycube(icube)
                          ycube(icube) = ytmp
                          zcube(icube) = ztmp
                          pers_corr    = zcentre/(zcentre-zcube(icube))
                          xcube(icube) = xcube(icube)*pers_corr
                          ycube(icube) = ycube(icube)*pers_corr
                       enddo
                       minx = minval(xcube)
                       maxx = maxval(xcube)
                       miny = minval(ycube)
                       maxy = maxval(ycube)
                       minz = minval(zcube)
                       maxz = maxval(zcube)
                    endif
                    if(shader_frame(proj_ind).ne.'cube')then
                       minx = -dx_proj
                       maxx =  dx_proj
                       miny = -dx_proj
                       maxy =  dx_proj
                       minz = -dx_proj
                       maxz =  dx_proj
                    endif
                    xleft   = xcentre+minx
                    xright  = xcentre+maxx
                    yleft   = ycentre+miny
                    yright  = ycentre+maxy
                    zleft   = zcentre+minz
                    zright  = zcentre+maxz
                    if(    xright.lt.xleft_frame.or.xleft.ge.xright_frame.or.&
                         & yright.lt.yleft_frame.or.yleft.ge.yright_frame.or.&
                         & zright.lt.zleft_frame.or.zleft.ge.zright_frame)cycle
#else
                    xx(i,1) = xx(i,1)-xcen
                    xx(i,2) = xx(i,2)-ycen
                    ! Rotating
                    xtmp    = cos(theta_cam)*xx(i,1)+sin(theta_cam)*xx(i,2)
                    ytmp    = cos(theta_cam)*xx(i,2)-sin(theta_cam)*xx(i,1)
                    xx(i,1) = xtmp
                    xx(i,2) = ytmp
                    xcentre = xx(i,1)+xcen
                    ycentre = xx(i,2)+ycen
                    xleft   = xcentre-dx_proj
                    xright  = xcentre+dx_proj
                    yleft   = ycentre-dx_proj
                    yright  = ycentre+dx_proj
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
                    do ii=imin,imax
                       ! Pixel x-axis position
                       xxleft      = xleft_frame+dble(ii-1)*dx_frame
                       xxright     = xxleft+dx_frame
                       xxcentre    = xxleft+0.5*dx_frame
                       dx_cell     = dx_frame/dx_proj
                       do jj=jmin,jmax
                          ! Pixel y-axis position
                          yyleft   = yleft_frame+dble(jj-1)*dy_frame
                          yyright  = yyleft+dy_frame
                          yycentre = yyleft+0.5*dy_frame
                          dy_cell  = dx_frame/dx_proj
                          xpc      = xxcentre-xcentre
                          ypc      = yycentre-ycentre
                          if(abs(xxcentre-xleft).lt.1d-2*dx_frame) xpc = xpc-1e-2*dx_frame
                          if(abs(yycentre-yleft).lt.1d-2*dx_frame) ypc = ypc-1e-2*dx_frame
                          if(abs(xxcentre-xright).lt.1d-2*dx_frame) xpc = xpc-1e-2*dx_frame
                          if(abs(yycentre-yright).lt.1d-2*dx_frame) ypc = ypc-1e-2*dx_frame
#if NDIM>2                 
                          if(shader_frame(proj_ind).eq.'cube')then
                             cube_face = .false.
                             if(sqrt(xpc**2+ypc**2).gt.dx_proj*sqrt(3.0)) goto 666
                             ! Filling the 6 cube shader faces
                             do iline=1,6
                                l1 = (ycube(lind(iline,2))-ycube(lind(iline,1)))**2+(xcube(lind(iline,2))-xcube(lind(iline,1)))**2
                                l2 = (ycube(lind(iline,4))-ycube(lind(iline,3)))**2+(xcube(lind(iline,4))-xcube(lind(iline,3)))**2
                                if(l1.eq.0d0) cycle
                                if(l2.eq.0d0) cycle
                                d1 = ((ycube(lind(iline,2))-ycube(lind(iline,1)))*xpc-(xcube(lind(iline,2))-xcube(lind(iline,1)))*ypc+xcube(lind(iline,2))*ycube(lind(iline,1))-ycube(lind(iline,2))*xcube(lind(iline,1)))/l1
                                d2 = ((ycube(lind(iline,4))-ycube(lind(iline,3)))*xpc-(xcube(lind(iline,4))-xcube(lind(iline,3)))*ypc+xcube(lind(iline,4))*ycube(lind(iline,3))-ycube(lind(iline,4))*xcube(lind(iline,3)))/l2
                                if(d1.eq.-sign(d1,d2)) cube_face=.true.
                                if(.not.cube_face) cycle
                                l3 = (ycube(lind(iline,6))-ycube(lind(iline,5)))**2+(xcube(lind(iline,6))-xcube(lind(iline,5)))**2
                                l4 = (ycube(lind(iline,8))-ycube(lind(iline,7)))**2+(xcube(lind(iline,8))-xcube(lind(iline,7)))**2
                                if(l3.eq.0d0) cycle
                                if(l4.eq.0d0) cycle
                                d3 = ((ycube(lind(iline,6))-ycube(lind(iline,5)))*xpc-(xcube(lind(iline,6))-xcube(lind(iline,5)))*ypc+xcube(lind(iline,6))*ycube(lind(iline,5))-ycube(lind(iline,6))*xcube(lind(iline,5)))/l3
                                d4 = ((ycube(lind(iline,8))-ycube(lind(iline,7)))*xpc-(xcube(lind(iline,8))-xcube(lind(iline,7)))*ypc+xcube(lind(iline,8))*ycube(lind(iline,7))-ycube(lind(iline,8))*xcube(lind(iline,7)))/l4
                                ! Within the projected face?
                                if(d3.eq.sign(d3,d4)) cube_face=.false.
                                if(cube_face) exit
                             enddo
666                          continue
                          endif
#endif
                          if((shader_frame(proj_ind).eq.'cube'          &
                             .and.(cube_face))                          &
                             .or.(shader_frame(proj_ind).eq.'sphere'    &
                             .and.sqrt(xpc**2+ypc**2).le.dx_proj)       &
                             .or.(shader_frame(proj_ind).eq.'square'    &
                             .and.(abs(xpc).lt.dx_proj)                 &
                             .and.(abs(ypc).lt.dx_proj)))then
                             ! Intersection volume
                             dvol      = dx_cell*dy_cell
                             if(method_frame(proj_ind).eq.'mean_mass')then
                                weight = dvol*uold(ind_cell(i),1)*dx_loc**3
                             elseif(method_frame(proj_ind).eq.'mean_dens')then
                                weight = dvol*uold(ind_cell(i),1)
                             elseif(method_frame(proj_ind).eq.'mean_vol')then
                                weight = dvol*dx_loc**3
                             elseif(method_frame(proj_ind).eq.'sum')then
                                weight = 1.0
                             endif
                             ! Update weights map
                             if(method_frame(proj_ind)(1:4).eq.'mean')then
                                weights(ii,jj) = weights(ii,jj)+weight
                             endif

			     imap = 1
#ifdef SOLVERmhd
                             do kk=0,NVAR+4
#else                       
                             do kk=0,NVAR
#endif
                                if(movie_vars(kk).eq.1)then
                                   ! Temperature map case
                                   if(kk==0)then
                                      e = 0.0d0
                                      do idim=1,ndim
                                         e = e+0.5*uold(ind_cell(i),idim+1)**2/max(uold(ind_cell(i),1),smallr)
                                      enddo
#if NENER>0
                                      do irad=0,nener-1
                                         e = e+uold(ind_cell(i),inener+irad)
                                      enddo
#endif
#ifdef SOLVERmhd
                                      do idim=1,ndim 
                                         e = e+0.125d0*(uold(ind_cell(i),idim+ndim+2)+uold(ind_cell(i),idim+nvar))**2
                                      enddo
#endif
                                      uvar = (gamma-1.0)*(uold(ind_cell(i),ndim+2)-e)
                                      uvar = uvar/uold(ind_cell(i),1)*scale_T2
                                   endif
                                   ! Density map case
                                   if(kk==1) then
                                      uvar = uold(ind_cell(i),kk)
                                   endif 
                                   ! Other scalars map
                                   if(kk>1)then
                                      uvar = uold(ind_cell(i),kk)/max(uold(ind_cell(i),1),smallr)
                                   endif
#ifdef SOLVERmhd
                                   ! Magnetic energy map case
                                   if(kk.eq.NVAR+4)then
                                      uvar = 0.125*(uold(ind_cell(i),6)**2+uold(ind_cell(i),7)**2+uold(ind_cell(i),8)**2 &
                                      &    + uold(ind_cell(i),NVAR+1)**2+uold(ind_cell(i),NVAR+2)**2+uold(ind_cell(i),NVAR+3)**2)
                                   endif
#endif
                                   ! Frame update
                                   if(method_frame(proj_ind).eq.'min')then
                                      data_frame(ii,jj,imap) = min(data_frame(ii,jj,imap),uvar)
                                   elseif(method_frame(proj_ind).eq.'max')then
                                      data_frame(ii,jj,imap) = max(data_frame(ii,jj,imap),uvar)
                                   else
                                      data_frame(ii,jj,imap) = data_frame(ii,jj,imap)+weight*uvar
                                   endif
                                   imap = imap+1
                                endif
                             end do
#ifdef RT
                             if(rt) then
			        imap = 1
                                do kk=1,NGROUPS
                                   if(rt_movie_vars(kk).eq.1) then
                                      if(method_frame(proj_ind).eq.'min')then
                                         rt_data_frame(ii,jj,imap) = &
                                         &   min(rt_data_frame(ii,jj,imap),rtuold(ind_cell(i),1+(kk-1)*(ndim+1))*rt_c_cgs*uold(ind_cell(i),1)
                                      elseif(method_frame(proj_ind).eq.'max')then
                                         rt_data_frame(ii,jj,imap) = &
                                         &   max(rt_data_frame(ii,jj,imap),rtuold(ind_cell(i),1+(kk-1)*(ndim+1))*rt_c_cgs*uold(ind_cell(i),1)
                                      else
                                         rt_data_frame(ii,jj,imap) = rt_data_frame(ii,jj,imap)+ &
                                         &   weight*rtuold(ind_cell(i),1+(kk-1)*(ndim+1))*rt_c_cgs*uold(ind_cell(i),1)
                                      endif
                                      imap = imap+1
                                   endif
                                end do
                             endif
#endif
                          endif
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
     xpf  = xp(j,1)-xcen
     ypf  = xp(j,2)-ycen
     zpf  = xp(j,3)-zcen
     ! Projection
     xtmp = cos(theta_cam)*xpf+sin(theta_cam)*ypf
     ytmp = cos(theta_cam)*ypf-sin(theta_cam)*xpf
     xpf  = xtmp
     ypf  = ytmp
     ytmp = cos(phi_cam)*ypf+sin(phi_cam)*zpf
     ztmp = cos(phi_cam)*zpf-sin(phi_cam)*ypf

     if(proj_axis(proj_ind:proj_ind).eq.'x')then
       xpf = ytmp
       ypf = ztmp
       zpf = xtmp
     elseif(proj_axis(proj_ind:proj_ind).eq.'y')then
       xpf = xtmp
       ypf = ztmp
       zpf = ytmp
     else
       xpf = xtmp
       ypf = ytmp
       zpf = ztmp
     endif
     if(perspective_camera(proj_ind))then
        xpf  = xpf*focal_camera(proj_ind)/(dist_cam-zpf)
        ypf  = ypf*focal_camera(proj_ind)/(dist_cam-zpf)
     endif
     xpf  = xpf+xcen
     ypf  = ypf+ycen
     zpf  = zpf+zcen
     
     if(    xpf.lt.xleft_frame.or.xpf.ge.xright_frame.or.&
          & ypf.lt.yleft_frame.or.ypf.ge.yright_frame.or.&
          & zpf.lt.zleft_frame.or.zpf.ge.zright_frame)cycle
#else
     xpf  = xp(j,1)-xcen
     ypf  = xp(j,2)-ycen
     xtmp = cos(theta_cam)*xpf+sin(theta_cam)*xpf
     ytmp = cos(theta_cam)*ypf-sin(theta_cam)*ypf
     xpf  = xtmp
     ypf  = ytmp
     xpf  = xpf+xcen
     ypf  = ypf+xcen
     if(    xpf.lt.xleft_frame.or.xpf.ge.xright_frame.or.&
          & ypf.lt.yleft_frame.or.ypf.ge.yright_frame)cycle
#endif
     ! Compute map indices for the cell
     ii = min(int((xpf-xleft_frame)/dx_frame)+1,nw_frame)
     jj = min(int((ypf-yleft_frame)/dy_frame)+1,nh_frame)
     
     ! Fill up map with projected mass
#ifdef SOLVERmhd
     ipart_start = NVAR+5
#else
     ipart_start = NVAR+1
#endif
     imap = 1
     do kk=0,ipart_start+1
        if(movie_vars(kk).eq.1)then
           if(star) then
              ! DM particles
              if((tp(j).eq.0d0).and.(kk.eq.ipart_start)) then
                 if(mass_cut_refine>0.0.and.zoom_only_frame(proj_ind)) then
                    if(mp(j)<mass_cut_refine) data_frame(ii,jj,imap)=data_frame(ii,jj,imap)+mp(j)
                 else
                    data_frame(ii,jj,imap)=data_frame(ii,jj,imap)+mp(j)
                 endif
              endif
              ! Star particles
              if((tp(j).ne.0d0).and.(kk.eq.ipart_start+1)) then
                 data_frame(ii,jj,imap)=data_frame(ii,jj,imap)+mp(j)
              endif
           else
              ! DM particles only
              if(kk.eq.ipart_start) then
                 if(mass_cut_refine>0d0.and.zoom_only_frame(proj_ind)) then
                    if(mp(j)<mass_cut_refine) data_frame(ii,jj,imap)=data_frame(ii,jj,imap)+mp(j)
                 else
                    data_frame(ii,jj,imap)=data_frame(ii,jj,imap)+mp(j)
                 endif
              endif
           endif
           imap = imap+1
        endif
     enddo
  end do
  ! End loop over particles

#ifndef WITHOUTMPI
  ! Maps communication
  allocate(data_single(1:nw_frame*nh_frame),stat=ierr)
  if(ierr .ne. 0)then
     write(*,*) 'Error - Movie frame allocation failed'
     call MPI_ABORT(MPI_COMM_WORLD,1,info)
  endif
  allocate(data_single_all(1:nw_frame*nh_frame),stat=ierr)
  if(ierr .ne. 0)then
     write(*,*) 'Error - Movie frame allocation failed'
     call MPI_ABORT(MPI_COMM_WORLD,1,info)
  endif
  ! Loop over maps
  do iframe=1,nframes
     i=1
     ! Load data in comm arrays
     do ii=1,nw_frame
        do jj=1,nh_frame
           data_single(i) = data_frame(ii,jj,iframe)
           i = i+1
        enddo
     enddo
     if(method_frame(proj_ind).eq.'min')then
        call MPI_REDUCE(data_single,data_single_all,nw_frame*nh_frame,MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,info)
     elseif(method_frame(proj_ind).eq.'max')then
        call MPI_REDUCE(data_single,data_single_all,nw_frame*nh_frame,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,info)
     else
        call MPI_REDUCE(data_single,data_single_all,nw_frame*nh_frame,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,info)
     endif
     ! Fill master pid frame
     if(myid==1)then
        i=1
        do ii=1,nw_frame
           do jj=1,nh_frame
              data_frame(ii,jj,iframe) = data_single_all(i)
              i = i+1
           enddo
        enddo
     endif
  enddo
  if(info.ne.MPI_SUCCESS)then
     if(myid==1) write(*,*) 'MPI error - map reduce failed'
     call MPI_ABORT(MPI_COMM_WORLD,1,info)
  endif
  ! Weights communication
  if(method_frame(proj_ind)(1:4).eq.'mean')then
     i=1
     do ii=1,nw_frame
        do jj=1,nh_frame
           data_single(i) = weights(ii,jj)
           i = i+1
        enddo
     enddo
     call MPI_REDUCE(data_single,data_single_all,nw_frame*nh_frame,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,info)
     if(myid==1)then
        i=1
        do ii=1,nw_frame
           do jj=1,nh_frame
              weights(ii,jj) = data_single_all(i)
              i = i+1
           enddo
        enddo
     endif
     if(info.ne.MPI_SUCCESS)then
        if(myid==1) write(*,*) 'MPI error - weigths reduce failed'
        call MPI_ABORT(MPI_COMM_WORLD,1,info)
     endif
  endif
#ifdef RT
  if(rt) then
     if(ierr .ne. 0)then
        call MPI_ABORT(MPI_COMM_WORLD,1,info)
     endif
     do iframe=1,rt_nframes
        i=1
        do ii=1,nw_frame
           do jj=1,nh_frame
              data_single(i) = rt_data_frame(ii,jj,iframe)
              i = i+1
           enddo
        enddo
        if(method_frame(proj_ind).eq.'min')then
           call MPI_REDUCE(data_single,data_single_all,nw_frame*nh_frame,MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,info)
        elseif(method_frame(proj_ind).eq.'max')then
           call MPI_REDUCE(data_single,data_single_all,nw_frame*nh_frame,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,info)
        else
           call MPI_REDUCE(data_single,data_single_all,nw_frame*nh_frame,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,info)
        endif
        if(myid==1)then
           i=1
           do ii=1,nw_frame
              do jj=1,nh_frame
                 rt_data_frame(ii,jj,iframe) = data_single_all(i)
                 i = i+1
              enddo
           enddo
        endif
     enddo
  endif
#endif
  deallocate(data_single)
  deallocate(data_single_all)
#endif
  if(myid==1)then
     if(method_frame(proj_ind)(1:4).ne.'sum')then
        ! Convert into mass weighted                                                                                                         
        do ii=1,nw_frame
          do jj=1,nh_frame
            imap = 1
            do kk=0,NVAR
              if(movie_vars(kk).eq.1)then
                 if((method_frame(proj_ind)(1:4).eq.'mean').and.(weights(ii,jj).gt.0d0))then
                    data_frame(ii,jj,imap) = data_frame(ii,jj,imap)/weights(ii,jj)
                 endif
                 if(method_frame(proj_ind)(1:4).eq.'min'.and.data_frame(ii,jj,imap).ge.1e-3*huge(0.0))then
                    data_frame(ii,jj,imap) = 0.0
                 endif
                 if(method_frame(proj_ind)(1:4).eq.'max'.and.data_frame(ii,jj,imap).le.-1e-3*huge(0.0))then
                    data_frame(ii,jj,imap) = 0.0
                 endif
                 imap = imap+1
              endif
            end do
#ifdef RT
            if(rt) then
               imap = 1
               do kk=1,NGROUPS
                  if(rt_movie_vars(kk).eq.1)then
                     if((method_frame(proj_ind)(1:4).eq.'mean').and.(weights(ii,jj).gt.0d0))then
                        rt_data_frame(ii,jj,imap) = data_frame(ii,jj,imap)/weights(ii,jj)
                     endif
                     if(method_frame(proj_ind)(1:4).eq.'min'.and.rt_data_frame(ii,jj,imap).ge.1e-3*huge(0.0))then
                        rt_data_frame(ii,jj,imap) = 0.0
                     endif
                     if(method_frame(proj_ind)(1:4).eq.'max'.and.rt_data_frame(ii,jj,imap).le.-1e-3*huge(0.0))then
                        rt_data_frame(ii,jj,imap) = 0.0
                     endif
                     imap = imap+1
                  endif
               end do
            endif
#endif
          end do
        end do
     endif
     if(method_frame(proj_ind)(1:4).eq.'min')then
     endif
  endif
  if(method_frame(proj_ind)(1:4).eq.'mean') deallocate(weights)

  if(myid==1)then
     ilun = 10
     if(ierr .ne. 0)then
        write(*,*) 'Error - Cannot alllocate movie frame'
#ifndef WITHOUTMPI
        call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
        stop
#endif
     endif
     ! Output mass weighted density
     imap = 1
#ifdef SOLVERmhd
     do kk=0, NVAR+6
#else
     do kk=0, NVAR+2
#endif
       if (movie_vars(kk).eq.1)then
         open(ilun,file=TRIM(moviefiles(kk)),form='unformatted',iostat=ierr)
         if(ierr .ne. 0)then
            write(*,*) 'Error - Could not open ',TRIM(moviefiles(kk))
#ifndef WITHOUTMPI
            call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
            stop
#endif
	 endif
         rewind(ilun)  
         if(tendmov>0)then
            write(ilun)t,delx,dely,delz
         else
            write(ilun)aexp,delx,dely,delz
         endif
         write(ilun)nw_frame,nh_frame
         write(ilun) real(data_frame(:,:,imap),4)
         close(ilun)
         ilun = ilun+1
         imap = imap+1
       end if
     end do

#ifdef RT
     if(rt) then
        imap = 1
        do kk=1, NGROUPS
           if (rt_movie_vars(kk).eq.1) then
              open(ilun,file=TRIM(rt_moviefiles(kk)),form='unformatted',iostat=ierr)
              if(ierr .ne. 0)then
                 write(*,*) 'Error - Could not open ',TRIM(moviefiles(kk))
#ifndef WITHOUTMPI
                 call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
                 stop
#endif
              endif
              rewind(ilun)  
              if(tendmov>0)then
                 write(ilun)t,delx,dely,delz
              else
                 write(ilun)aexp,delx,dely,delz
              endif
              write(ilun)nw_frame,nh_frame
              write(ilun)real(rt_data_frame(:,:,imap),4)
              close(ilun)
              ilun = ilun+1
              imap = imap+1
           end if
        end do
     endif
#endif
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
#ifdef RT
  do ll=1,NGROUPS
    write(dummy,'(I3.1)') ll
    if(ANY(movie_vars_txt=='Fp'//trim(adjustl(dummy))//' ')) rt_movie_vars(ll)=1
  enddo
#endif

end subroutine set_movie_vars
