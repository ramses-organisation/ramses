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
  use constants, only: pi, c_cgs, L_sun, M_sun, yr2sec
  use mpi_mod
  use file_module, ONLY: mkdir
  implicit none
#if NDIM > 1
#ifndef WITHOUTMPI
  integer::info
  real(kind=8),dimension(:),allocatable::data_single,data_single_all
#endif

#ifdef NOSYSTEM
  integer::info2
#endif

  integer::ierr
  integer,parameter::tag=100

  character(len=5)::istep_str
  character(len=100)::moviedir,moviecmd,infofile,sinkfile,filename
#if NENER>0
  integer::irad
#endif
#ifdef RT
  character(len=100)::rt_infofile
#endif
  integer::ncache,iskip,ngrid,nlevelmax_frame
  integer::idim,ilun,nx_loc,ind,ix,iy,iz
  integer::imin,imax,jmin,jmax,ii,jj,kk
  real(dp)::scale,scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::xcen,ycen,zcen,delx,dely,delz,timer
  real(dp)::xtmp,ytmp,theta_cam,phi_cam,fov_camera,dist_cam
  real(dp)::xleft_frame,xright_frame,yleft_frame,yright_frame,zleft_frame,zright_frame
  real(dp)::xleft,xright,yleft,yright,xcentre,ycentre
  real(dp)::xxleft,xxright,yyleft,yyright,xxcentre,yycentre
  real(dp)::xpf,ypf
  real(dp)::dx_frame,dy_frame,dx,dx_loc,dx_min
  real(dp)::dx_cell,dy_cell,dvol,dx_proj,weight=0
  integer,dimension(1:nvector)::ind_grid,ind_cell
  logical,dimension(1:nvector)::ok
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim)::xx
  real(kind=8),dimension(:,:,:),allocatable::data_frame
  real(kind=8),dimension(:,:),allocatable::weights
  real(kind=8)::e,uvar
  integer, parameter :: mode = int(O'755')
  integer::igrid,ilevel
  integer::i,j,ivar
  integer::proj_ind,nh_temp,nw_temp,proj_ax
  real(dp)::xpc,ypc
  logical::opened,cube_face=.true.,ok_frame
  logical::is_cube, is_sphere, is_square, is_mean_mass, is_mean_dens
  logical::is_mean_vol, is_sum, is_mean, is_min, is_max
  character(len=1)::temp_string
  integer,dimension(6,8)::lind = reshape((/1, 2, 3, 4, 1, 3, 2, 4,    &
                                           5, 6, 7, 8, 5, 7, 6, 8,    &
                                           1, 5, 2, 6, 1, 2, 5, 6,    &
                                           3, 7, 4, 8, 3, 4, 7, 8,    &
                                           1, 3, 5, 7, 1, 5, 3, 7,    &
                                           2, 4, 6, 8, 2, 6, 4, 8 /)  &
                                           ,shape(lind),order=(/2,1/))
  integer::npoly
  real(dp)::msol,lumsol,year,log_lum
  real(dp),dimension(46)::lum_poly = (/-5.6801548098d+05,5.9946628460d+05,-2.3731255261d+05,3.8560177896d+04,   &
                                       -5.6808026741d+02,-4.1194896342d+02,-5.4456684492d+00,4.4307467151d+00,  &
                                       4.2381242791d-01,-1.0831618977d-02,-6.2879957495d-03,-6.4348257709d-04,  &
                                       -8.7868432521d-06,6.6149289040d-06,1.0739840644d-06,8.1516389787d-08,    &
                                       -4.9767430769d-10,-1.0869960505d-09,-1.7184492747d-10,-1.5168404854d-11, &
                                       -4.3154482167d-13,1.1515161236d-13,2.5427591874d-14,3.0240690024d-15,    &
                                       2.1644677501d-16,5.4516441993d-19,-2.6382406892d-18,-4.8388838089d-19,   &
                                       -5.4310814625d-20,-3.7196590586d-21,1.6452219840d-23,5.0868990738d-23,   &
                                       9.0947242158d-24,1.0014439287d-24,6.3725561865d-26,-1.6122919369d-27,    &
                                       -1.1379163662d-27,-1.8221651341d-28,-1.6791574329d-29,-4.0655271535d-31, &
                                       1.6771769041d-31,3.3366386465d-32,2.6122101904d-33,-1.7749944262d-34,    &
                                       -6.6026967339d-35,3.9721180887d-36 /)
#if NDIM>2
  real(dp)::ztmp,alpha=0,beta=0,pers_corr
  real(dp)::zleft,zright,zcentre,zpf
  real(dp),dimension(8)::xcube,ycube,zcube
  integer::icube,iline
  real(dp)::minx=0,maxx=0,miny=0,maxy=0,minz=0,maxz=0,d1,d2,d3,d4,l1,l2,l3,l4
#endif

 nh_temp = nh_frame
 nw_temp = nw_frame

 if(n_movie_vars==0) return

 ! Only one projection available in 2D
 if((ndim.eq.2).and.(trim(proj_axis).ne.'z')) proj_axis = 'z'

 ! Loop over projections
 do proj_ind=1,LEN(trim(proj_axis))

    opened=.false.

    if(imov<1)imov=1
    if(imov>imovout)return

    ! No cubic shader for 2D simulations
    if((ndim.eq.2).and.(shader_frame(proj_ind).eq.'cube')) shader_frame(proj_ind) = 'square'

    ! Some booleans for speedup
    is_cube=.false.
    if(shader_frame(proj_ind).eq.'cube')      is_cube=.true.
    is_sphere=.false.
    if(shader_frame(proj_ind).eq.'sphere')    is_sphere=.true.
    is_square=.false.
    if(shader_frame(proj_ind).eq.'square')    is_square=.true.
    is_mean_mass=.false.
    if(method_frame(proj_ind).eq.'mean_mass') is_mean_mass=.true.
    is_mean_dens=.false.
    if(method_frame(proj_ind).eq.'mean_dens') is_mean_dens=.true.
    is_mean_vol=.false.
    if(method_frame(proj_ind).eq.'mean_vol')  is_mean_vol=.true.
    is_mean_vol=.false.
    if(method_frame(proj_ind).eq.'sum')       is_sum=.true.
    is_mean=.false.
    if(method_frame(proj_ind)(1:4).eq.'mean') is_mean=.true.
    is_min=.false.
    if(method_frame(proj_ind).eq.'min')       is_min=.true.
    is_max=.false.
    if(method_frame(proj_ind).eq.'max')       is_max=.true.

    if(proj_axis(proj_ind:proj_ind).eq.'x') then
       proj_ax=1
    else if(proj_axis(proj_ind:proj_ind).eq.'y') then
       proj_ax=2
    else
       proj_ax=3
    endif

    ! Determine the filename, dir, etc
    if(myid==1)write(*,*)'Computing and dumping movie frame'

    call title(imov, istep_str)
    write(temp_string,'(I1)') proj_ind
    moviedir = 'movie'//trim(temp_string)//'/'
    moviecmd = 'mkdir -p '//trim(moviedir)
    if(.not.withoutmkdir) then
#ifdef NOSYSTEM
       if(myid==1)call PXFMKDIR(TRIM(moviedir),LEN(TRIM(moviedir)),O'755',info2)
#else
       if(myid==1)then
          ierr=1
          ! call system(moviecmd,ierr)
          ! call EXECUTE_COMMAND_LINE(moviecmd,exitstat=ierr,wait=.true.)
          call mkdir(trim(moviedir),mode,ierr)
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

#ifdef RT
    if(rt)then
       rt_infofile = trim(moviedir)//'rt_info_'//trim(istep_str)//'.txt'
       if(myid==1.and.proj_ind==1) call output_rtinfo(rt_infofile)
    endif
#endif

    ! Conversion factor from user units to cgs units
    call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
    msol   = M_sun/(scale_d*scale_l**3)
    lumsol = L_sun*(scale_t/(c_cgs*scale_v)**2)
    year   = yr2sec/scale_t
    if(cosmo) year = year*aexp**2

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
    if(cosmo) then
       timer = aexp
    else
       timer = t
    endif

    ! Compute frame boundaries
    xcen=xcentre_frame(proj_ind*4-3)+xcentre_frame(proj_ind*4-2)*timer+xcentre_frame(proj_ind*4-1)*timer**2+xcentre_frame(proj_ind*4)*timer**3
    ycen=ycentre_frame(proj_ind*4-3)+ycentre_frame(proj_ind*4-2)*timer+ycentre_frame(proj_ind*4-1)*timer**2+ycentre_frame(proj_ind*4)*timer**3
    zcen=zcentre_frame(proj_ind*4-3)+zcentre_frame(proj_ind*4-2)*timer+zcentre_frame(proj_ind*4-1)*timer**2+zcentre_frame(proj_ind*4)*timer**3
    if(deltax_frame(proj_ind*2-1).eq.0d0 .and. deltay_frame(proj_ind*2-1).gt.0d0)then
       deltax_frame(proj_ind*2-1)=deltay_frame(proj_ind*2-1)*float(nw_frame)/float(nh_frame)
    endif
    if(deltay_frame(proj_ind*2-1).eq.0d0 .and. deltax_frame(proj_ind*2-1).gt.0d0)then
       deltay_frame(proj_ind*2-1)=deltax_frame(proj_ind*2-1)*float(nh_frame)/float(nw_frame)
    endif
    delx=deltax_frame(proj_ind*2-1)+deltax_frame(proj_ind*2)/aexp
    dely=deltay_frame(proj_ind*2-1)+deltay_frame(proj_ind*2)/aexp
    delz=deltaz_frame(proj_ind*2-1)+deltaz_frame(proj_ind*2)/aexp
    if(dist_camera(proj_ind).le.0D0) dist_camera(proj_ind) = boxlen

    ! Camera properties
    if(cosmo) then
       if(tend_theta_camera(proj_ind).le.0d0) tend_theta_camera(proj_ind) = aendmov
       if(tend_phi_camera(proj_ind).le.0d0) tend_phi_camera(proj_ind) = aendmov
       theta_cam  = theta_camera(proj_ind)*pi/180d0                                                                                 &
            +min(max(aexp-tstart_theta_camera(proj_ind),0d0),tend_theta_camera(proj_ind))*dtheta_camera(proj_ind)*pi/180./(aendmov-astartmov)
       phi_cam    = phi_camera(proj_ind)*pi/180d0                                                                                   &
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
    if(myid==1) then
       write(*,'(5A,F6.1,A,F6.1)',advance='no') ' Writing frame ', istep_str,' los=',proj_axis(proj_ind:proj_ind),   &
            &                                              ' theta=',theta_cam*180./pi,' phi=',phi_cam*180./pi
       if(perspective_camera(proj_ind))then
          write(*,'(A,F8.2,A,F6.1)') ' d=',dist_cam,' fov=',fov_camera*180./pi
       else
          write(*,'(A)') ''
       endif
    endif
#else
    if(myid==1) write(*,'(3A,F6.1)') " Writing frame ", istep_str,' theta=',theta_cam*180./pi
#endif

    ! Frame boundaries
    if(proj_ax.eq.1) then ! x-projection
       xleft_frame  = ycen-dely/2.
       xright_frame = ycen+dely/2.
       yleft_frame  = zcen-delz/2.
       yright_frame = zcen+delz/2.
       zleft_frame  = xcen-delx/2.
       zright_frame = xcen+delx/2.

       dx_frame = dely/dble(nw_frame)
       dy_frame = delz/dble(nh_frame)
    elseif(proj_ax.eq.2) then ! y-projection
       xleft_frame  = xcen-delx/2.
       xright_frame = xcen+delx/2.
       yleft_frame  = zcen-delz/2.
       yright_frame = zcen+delz/2.
       zleft_frame  = ycen-dely/2.
       zright_frame = ycen+dely/2.

       dx_frame = delx/dble(nw_frame)
       dy_frame = delz/dble(nh_frame)
    else                      ! z-projection
       xleft_frame  = xcen-delx/2.
       xright_frame = xcen+delx/2.
       yleft_frame  = ycen-dely/2.
       yright_frame = ycen+dely/2.
       zleft_frame  = zcen-delz/2.
       zright_frame = zcen+delz/2.

       dx_frame = delx/dble(nw_frame)
       dy_frame = dely/dble(nh_frame)
    endif

    ! Allocate image
    allocate(data_frame(1:nw_frame,1:nh_frame,1:n_movie_vars),stat=ierr)
    if(ierr .ne. 0)then
       write(*,*) 'Error - Movie frame allocation failed'
#ifndef WITHOUTMPI
       call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
       stop
#endif
    endif

    allocate(weights(1:nw_frame,1:nh_frame),stat=ierr)
    if(ierr .ne. 0)then
       write(*,*) 'Error - Movie frame allocation failed'
#ifndef WITHOUTMPI
       call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
       stop
#endif
    endif

    if(is_min)then
       data_frame(:,:,:) = 1e-3*huge(0.0)
    elseif(is_max)then
       data_frame(:,:,:) = -1e-3*huge(0.0)
    else
       data_frame(:,:,:) = 0.0
    endif
    weights(:,:) = 0d0

    ! Deal with hydro variables
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
          if(is_cube .and. (.not.perspective_camera(proj_ind)))then
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
#ifdef SOLVERmhd
                      ! Scale temperature to K
                      if(ivar_frame(proj_ind)==5)then
                         e = 0.0d0
                         do idim=1,3
#else
                      ! Scale temperature to K
                      if(ivar_frame(proj_ind)==ndim+2)then
                         e = 0.0d0
                         do idim=1,ndim
#endif
                            e = e+0.5*uold(ind_cell(i),idim+1)**2/uold(ind_cell(i),1)
                         enddo
#if NENER>0
                         do irad=0,nener-1
                            e = e+uold(ind_cell(i),inener+irad)
                         enddo
#endif
#ifdef SOLVERmhd
                         do idim=1,3
                            e = e+0.125d0*(uold(ind_cell(i),idim+5)+uold(ind_cell(i),idim+nvar))**2
                         enddo
#endif
                         ! Pressure
                         uvar = (gamma-1.0)*(uold(ind_cell(i),ivar_frame(proj_ind))-e)*scale_T2
                      endif
                      ! Switch to primitive variables
                      if(ivar_frame(proj_ind)>1) uvar = uvar/uold(ind_cell(i),1)
                      ! Scale density to cm**-3
                      if(ivar_frame(proj_ind)==1) uvar = uvar*scale_nH
#ifdef SOLVERmhd
                      ! Scale velocities to km/s
                      if(ivar_frame(proj_ind)>1.and.ivar_frame(proj_ind)<5) uvar = uvar*scale_v/1e5
#else
                      ! Scale velocities to km/s
                      if(ivar_frame(proj_ind)>1.and.ivar_frame(proj_ind)<ndim+2) uvar = uvar*scale_v/1e5
#endif
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
                      if(proj_ax.eq.1) then ! x-projection
                         if(dist_cam-xx(i,1).lt.0d0) cycle
                         if(perspective_camera(proj_ind))then
                            alpha  = atan(xx(i,2)/(dist_cam-xx(i,1)))
                            beta   = atan(xx(i,3)/(dist_cam-xx(i,1)))
                            if(abs(alpha)/2d0.gt.fov_camera) cycle
                            if(abs(beta)/2d0.gt.fov_camera) cycle
                            pers_corr = focal_camera(proj_ind)/(dist_cam-xx(i,1))
                            xx(i,2)   = xx(i,2)*pers_corr
                            xx(i,3)   = xx(i,3)*pers_corr
                            dx_proj   = (dx_loc/2.0)*pers_corr*smooth_frame(proj_ind)
                         endif
                         xcentre = xx(i,2)+ycen
                         ycentre = xx(i,3)+zcen
                         zcentre = xx(i,1)+xcen
                      elseif(proj_ax.eq.2) then ! y-projection
                         if(dist_cam-xx(i,2).lt.0d0) cycle
                         if(perspective_camera(proj_ind))then
                            alpha  = atan(xx(i,1)/(dist_cam-xx(i,2)))
                            beta   = atan(xx(i,3)/(dist_cam-xx(i,2)))
                            if(abs(alpha)/2d0.gt.fov_camera) cycle
                            if(abs(beta)/2d0.gt.fov_camera) cycle
                            pers_corr = focal_camera(proj_ind)/(dist_cam-xx(i,2))
                            xx(i,1)   = xx(i,1)*pers_corr
                            xx(i,3)   = xx(i,3)*pers_corr
                            dx_proj   = (dx_loc/2.0)*pers_corr*smooth_frame(proj_ind)
                         endif
                         xcentre = xx(i,1)+xcen
                         ycentre = xx(i,3)+zcen
                         zcentre = xx(i,2)+ycen
                      else
                         if(dist_cam-xx(i,3).lt.0d0) cycle
                         if(perspective_camera(proj_ind))then
                            alpha  = atan(xx(i,1)/(dist_cam-xx(i,3)))
                            beta   = atan(xx(i,2)/(dist_cam-xx(i,3)))
                            if(abs(alpha)/2d0.gt.fov_camera) cycle
                            if(abs(beta)/2d0.gt.fov_camera) cycle
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
                      if(is_cube .and. perspective_camera(proj_ind))then
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
                      if(.not. is_cube)then
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
                            if(is_cube)then
                               cube_face = .false.
                               if(sqrt(xpc**2+ypc**2).gt.dx_proj*sqrt(3.0)) goto 666
                               ! Filling the 6 cube shader faces
                               do iline=1,6
                                  l1 = (ycube(lind(iline,2))-ycube(lind(iline,1)))**2+(xcube(lind(iline,2))-xcube(lind(iline,1)))**2
                                  l2 = (ycube(lind(iline,4))-ycube(lind(iline,3)))**2+(xcube(lind(iline,4))-xcube(lind(iline,3)))**2
                                  if(l1.eq.0d0) cycle
                                  if(l2.eq.0d0) cycle
                                  d1 = ((ycube(lind(iline,2))-ycube(lind(iline,1)))*xpc-(xcube(lind(iline,2))-xcube(lind(iline,1)))*ypc &
                                       & +xcube(lind(iline,2))*ycube(lind(iline,1))-ycube(lind(iline,2))*xcube(lind(iline,1)))/l1
                                  d2 = ((ycube(lind(iline,4))-ycube(lind(iline,3)))*xpc-(xcube(lind(iline,4))-xcube(lind(iline,3)))*ypc &
                                       & +xcube(lind(iline,4))*ycube(lind(iline,3))-ycube(lind(iline,4))*xcube(lind(iline,3)))/l2
                                  if(d1.eq.-sign(d1,d2)) cube_face=.true.
                                  if(.not.cube_face) cycle
                                  l3 = (ycube(lind(iline,6))-ycube(lind(iline,5)))**2+(xcube(lind(iline,6))-xcube(lind(iline,5)))**2
                                  l4 = (ycube(lind(iline,8))-ycube(lind(iline,7)))**2+(xcube(lind(iline,8))-xcube(lind(iline,7)))**2
                                  if(l3.eq.0d0) cycle
                                  if(l4.eq.0d0) cycle
                                  d3 = ((ycube(lind(iline,6))-ycube(lind(iline,5)))*xpc-(xcube(lind(iline,6))-xcube(lind(iline,5)))*ypc &
                                       & +xcube(lind(iline,6))*ycube(lind(iline,5))-ycube(lind(iline,6))*xcube(lind(iline,5)))/l3
                                  d4 = ((ycube(lind(iline,8))-ycube(lind(iline,7)))*xpc-(xcube(lind(iline,8))-xcube(lind(iline,7)))*ypc &
                                       & +xcube(lind(iline,8))*ycube(lind(iline,7))-ycube(lind(iline,8))*xcube(lind(iline,7)))/l4
                                  ! Within the projected face?
                                  if(d3.eq.sign(d3,d4)) cube_face=.false.
                                  if(cube_face) exit
                               enddo
666                            continue
                            endif
#endif
                            if((is_cube .and. (cube_face))                           &
                                 .or.(is_sphere .and. sqrt(xpc**2+ypc**2).le.dx_proj)  &
                                 .or.(is_square .and.(abs(xpc).lt.dx_proj)             &
                                 .and.(abs(ypc).lt.dx_proj)))                       then
                               ! Intersection volume
                               dvol      = dx_cell*dy_cell
                               if(is_mean_mass) then
                                  weight = dvol*uold(ind_cell(i),1)*dx_loc**3
                               elseif(is_mean_dens)then
                                  weight = dvol*uold(ind_cell(i),1)
                               elseif(is_mean_vol)then
                                  weight = dvol*dx_loc**3
                               elseif(is_sum)then
                                  weight = 1.0
                               endif
                               ! Update weights map
                               if(is_mean) weights(ii,jj) = weights(ii,jj)+weight

                               do kk=1,n_movie_vars
                                  ok_frame=.false.
                                  ! Temperature map case
                                  if(movie_vars(kk).eq.i_mv_temp)then
                                     ok_frame=.true.
                                     e = 0.0d0
#ifdef SOLVERmhd
                                     do idim=1,3
#else
                                     do idim=1,ndim
#endif
                                        e = e+0.5*uold(ind_cell(i),idim+1)**2/max(uold(ind_cell(i),1),smallr)
                                     enddo
#if NENER>0
                                     do irad=0,nener-1
                                        e = e+uold(ind_cell(i),inener+irad)
                                     enddo
#endif
#ifdef SOLVERmhd
                                     do idim=1,3
                                        e = e+0.125d0*(uold(ind_cell(i),idim+5)+uold(ind_cell(i),idim+nvar))**2
                                     enddo
                                     uvar = (gamma-1.0)*(uold(ind_cell(i),5)-e)
#else
                                     uvar = (gamma-1.0)*(uold(ind_cell(i),ndim+2)-e)
#endif
                                     uvar = uvar/uold(ind_cell(i),1)*scale_T2

                                     ! Density map
                                  else if(movie_vars(kk).eq.i_mv_dens)then
                                     ok_frame=.true.
                                     uvar = uold(ind_cell(i),1)

                                     ! Pressure map
                                  else if(movie_vars(kk).eq.i_mv_p)then
                                     ok_frame=.true.
                                     uvar = uold(ind_cell(i),ndim+2)

                                     ! Speed map
                                  else if(movie_vars(kk).eq.i_mv_speed)then
                                     ok_frame=.true.
                                     uvar=0.0d0
                                     do idim=1,ndim
                                        uvar = uvar + uold(ind_cell(i),idim+1)**2 / max(uold(ind_cell(i),1),smallr)
                                     end do
                                     uvar = sqrt(uvar)

                                     ! Velocity map
                                  else if(movie_vars(kk).eq.i_mv_vx)then
                                     ok_frame=.true.
                                     uvar = uold(ind_cell(i),2) / max(uold(ind_cell(i),1),smallr)
                                  else if(movie_vars(kk).eq.i_mv_vy)then
                                     ok_frame=.true.
                                     uvar = uold(ind_cell(i),3) / max(uold(ind_cell(i),1),smallr)
#if NDIM>2
                                  else if(movie_vars(kk).eq.i_mv_vz)then
                                     ok_frame=.true.
                                     uvar = uold(ind_cell(i),4) / max(uold(ind_cell(i),1),smallr)
#endif

                                     ! Metallicity map
                                  else if(movie_vars(kk).eq.i_mv_metallicity)then
                                     ok_frame=.true.
                                     uvar = uold(ind_cell(i),imetal)/max(uold(ind_cell(i),1),smallr)

                                     ! Any scalars map
                                  else if(movie_vars(kk).eq.i_mv_var)then
                                     ok_frame=.true.
                                     ivar = movie_var_number(kk)
                                     uvar = uold(ind_cell(i),ivar)/max(uold(ind_cell(i),1),smallr)
#ifdef SOLVERmhd
                                     ! Magnetic energy map case
                                  else if(movie_vars(kk).eq.i_mv_pmag)then
                                     ok_frame=.true.
                                     uvar=0.125*( (uold(ind_cell(i),6)+uold(ind_cell(i),NVAR+1))**2 &
                                          & + (uold(ind_cell(i),7)+uold(ind_cell(i),NVAR+2))**2 &
                                          & + (uold(ind_cell(i),8)+uold(ind_cell(i),NVAR+3))**2 )
#endif
                                  endif
#ifdef RT
                                  ! Ionization fraction map
                                  if(movie_vars(kk).eq.i_mv_xhi)then
                                     ok_frame=.true.
                                     if(ixhi .ne. 0) then
                                        uvar = uold(ind_cell(i),ichem-1+ixhi)/max(uold(ind_cell(i),1),smallr)
                                     else
                                        uvar = 1d0-uold(ind_cell(i),ichem-1+ixhii)/max(uold(ind_cell(i),1),smallr)
                                     endif
                                  else if(movie_vars(kk).eq.i_mv_xhii)then
                                     ok_frame=.true.
                                     uvar = uold(ind_cell(i),ichem-1+ixhii)/max(uold(ind_cell(i),1),smallr)
                                  else if(movie_vars(kk).eq.i_mv_xheii .and. ixheii .ne. 0)then
                                     ok_frame=.true.
                                     uvar = uold(ind_cell(i),ichem-1+ixheii)/max(uold(ind_cell(i),1),smallr)
                                  else if(movie_vars(kk).eq.i_mv_xheiii .and. ixheiii .ne. 0)then
                                     ok_frame=.true.
                                     uvar = uold(ind_cell(i),ichem-1+ixheiii)/max(uold(ind_cell(i),1),smallr)
                                  else if(movie_vars(kk).eq.i_mv_xh2 .and. ixhi .ne. 0)then
                                     ok_frame=.true.
                                     uvar = 1.-(uold(ind_cell(i),ichem-1+ixhi)+uold(ind_cell(i),ichem-1+ixhii))/max(uold(ind_cell(i),1),smallr)
                                  endif

                                  ! Photon map
                                  if(rt) then
                                     if(movie_vars(kk).eq.i_mv_fp)then
                                        ok_frame=.true.
                                        ivar = movie_var_number(kk)
                                        uvar=rtuold(ind_cell(i),1+(ivar-1)*(ndim+1))*rt_c
                                     endif
                                  endif ! if(rt)
#endif
                                  ! Frame update
                                  if(ok_frame) then
                                     if(is_min)then
                                        data_frame(ii,jj,kk) = min(data_frame(ii,jj,kk),uvar)
                                     elseif(is_max)then
                                        data_frame(ii,jj,kk) = max(data_frame(ii,jj,kk),uvar)
                                     else
                                        data_frame(ii,jj,kk) = data_frame(ii,jj,kk)+weight*uvar
                                     endif
                                  endif

                               end do ! loop over kk
                            endif ! if(shader is cube)
                         end do ! jj=jmin,jmax
                      end do ! ii=imin,imax
                   end if ! if(ok)
                end do ! i=1,ngrid

             end do
             ! End loop over cells

          end do
          ! End loop over grids
       end do
       ! End loop over levels
    end if
    ! End block if hydro

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

       if(proj_ax.eq.1) then ! x-projection
          xpf = ytmp
          ypf = ztmp
          zpf = xtmp
       elseif(proj_ax.eq.2) then ! y-projection
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
       if(proj_ax.eq.1) then ! x-projection
          xpf  = xpf+ycen
          ypf  = ypf+zcen
          zpf  = zpf+xcen
       elseif(proj_ax.eq.2) then ! y-projection
          xpf  = xpf+xcen
          ypf  = ypf+zcen
          zpf  = zpf+ycen
       else
          xpf  = xpf+xcen
          ypf  = ypf+ycen
          zpf  = zpf+zcen
       endif

       ! Check if particle is in front of camera
       if(dist_cam-zpf.lt.0) cycle

       ! Check if particle is in the movie box
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
       ypf  = ypf+ycen
       if(    xpf.lt.xleft_frame.or.xpf.ge.xright_frame.or.&
            & ypf.lt.yleft_frame.or.ypf.ge.yright_frame)cycle
#endif
       ! Compute map indices for the particle
       ii = min(int((xpf-xleft_frame)/dx_frame)+1,nw_frame)
       jj = min(int((ypf-yleft_frame)/dy_frame)+1,nh_frame)

       do kk=1,n_movie_vars
          if(movie_vars(kk) .eq. i_mv_dm .and. is_DM(typep(j)))then
             ! DM particles
             if(mass_cut_refine>0.0.and.zoom_only_frame(proj_ind)) then
                if(mp(j)<mass_cut_refine) data_frame(ii,jj,kk)=data_frame(ii,jj,kk)+mp(j)
             else
                data_frame(ii,jj,kk)=data_frame(ii,jj,kk)+mp(j)
             endif
          else if(star .and. movie_vars(kk) .eq. i_mv_stars .and. is_star(typep(j)))then
             ! Star particles
             data_frame(ii,jj,kk)=data_frame(ii,jj,kk)+mp(j)
          else if(star .and. movie_vars(kk) .eq. i_mv_lum .and. is_star(typep(j)))then
             ! Star particles luminosity in code units (luminosity over speed of light squared)
             ! The polynome is fitted on Starburst99 instantaneous bolometric magnitude
             ! for  Z = 0.04, alpha = 2.35, M_up = 100 Msol
             ! http://www.stsci.edu/science/starburst99/data/bol_inst_a.dat
             ! Polynome is poorly constrained on high and low ends
             if(log10((texp-tp(j))/year)<6)then
                log_lum = 3.2d0
             else if(log10((texp-tp(j))/year)>9)then
                log_lum = log10((texp-tp(j))/year)*(-9.79362D-01)+9.08855D+00
             else
                log_lum = 0d0
                do npoly=1,size(lum_poly)
                   log_lum = log_lum+lum_poly(npoly)*(log10((texp-tp(j))/year))**(npoly-1)
                enddo
             endif
             data_frame(ii,jj,kk)=data_frame(ii,jj,kk)+(10d0**(log_lum))*(mp(j)/msol)*lumsol
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
    do kk=1,n_movie_vars
       i=1
       ! Load data in comm arrays
       do ii=1,nw_frame
          do jj=1,nh_frame
             data_single(i) = data_frame(ii,jj,kk)
             i = i+1
          enddo
       enddo
       if(is_min)then
          call MPI_REDUCE(data_single,data_single_all,nw_frame*nh_frame,MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,info)
       elseif(is_max)then
          call MPI_REDUCE(data_single,data_single_all,nw_frame*nh_frame,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,info)
       else
          call MPI_REDUCE(data_single,data_single_all,nw_frame*nh_frame,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,info)
       endif
       ! Fill master pid frame
       if(myid==1)then
          i=1
          do ii=1,nw_frame
             do jj=1,nh_frame
                data_frame(ii,jj,kk) = data_single_all(i)
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
    if(is_mean)then
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
    deallocate(data_single)
    deallocate(data_single_all)
#endif
    if(myid==1)then
       if(method_frame(proj_ind)(1:3).ne.'sum')then
          ! Convert into mass weighted
          do kk=1,n_movie_vars
             if(movie_vars(kk) .eq. i_mv_dm) cycle
             if(movie_vars(kk) .eq. i_mv_stars) cycle
             if(movie_vars(kk) .eq. i_mv_lum) cycle
             do ii=1,nw_frame
                do jj=1,nh_frame
                   if((is_mean).and.(weights(ii,jj).gt.0d0))then
                      data_frame(ii,jj,kk) = data_frame(ii,jj,kk)/weights(ii,jj)
                   endif
                   if(method_frame(proj_ind)(1:3).eq.'min'.and.data_frame(ii,jj,kk).ge.1e-3*huge(0.0))then
                      data_frame(ii,jj,kk) = 0.0
                   endif
                   if(method_frame(proj_ind)(1:3).eq.'max'.and.data_frame(ii,jj,kk).le.-1e-3*huge(0.0))then
                      data_frame(ii,jj,kk) = 0.0
                   endif
                end do ! loop over jj
             end do ! loop over ii
          end do ! loop over kk
       endif
    endif
    if(is_mean) deallocate(weights)

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
       ! Write the frames to files
       do kk=1, n_movie_vars
          filename = trim(moviedir)//trim(movie_vars_txt(kk))//'_'//trim(istep_str)//'.map'
          open(ilun,file=TRIM(filename),form='unformatted',iostat=ierr)
          if(ierr .ne. 0)then
             write(*,*) 'Error - Could not open ',TRIM(filename)
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
          !print*,'Frame is ',movie_vars_txt(kk),minval(data_frame(:,:,kk)),maxval(data_frame(:,:,kk))
          write(ilun) real(data_frame(:,:,kk),4)
          close(ilun)
          ilun = ilun+1
       end do

    endif

    deallocate(data_frame)

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
 ! End loop over projections

#endif
end subroutine output_frame

subroutine set_movie_vars()
  use amr_commons
  implicit none
  integer::kk, ivar
  ! This routine sets up movie_vars to draw the correct
  ! variables

  ! Determine the number of frames to output
  n_movie_vars=0
  do
     if(TRIM(movie_vars_txt(n_movie_vars+1)).eq.'') exit
     n_movie_vars = n_movie_vars + 1
  end do

  do kk=1,n_movie_vars
     if(movie_vars_txt(kk) .eq. 'dm') then
        if(i_mv_dm .eq. -1) i_mv_dm = kk
        movie_vars(kk) = i_mv_dm

     else if (movie_vars_txt(kk) .eq. 'stars') then
        if(i_mv_stars .eq. -1) i_mv_stars = kk
        movie_vars(kk) = i_mv_stars

     else if (movie_vars_txt(kk) .eq. 'lum') then
        if(i_mv_lum .eq. -1) i_mv_lum = kk
        movie_vars(kk) = i_mv_lum

     else if (movie_vars_txt(kk) .eq. 'temp') then
        if(i_mv_temp .eq. -1) i_mv_temp = kk
        movie_vars(kk) = i_mv_temp

     else if (movie_vars_txt(kk) .eq. 'dens') then
        if(i_mv_dens .eq. -1) i_mv_dens = kk
        movie_vars(kk) = i_mv_dens

     else if (movie_vars_txt(kk) .eq. 'P') then
        if(i_mv_dens .eq. -1) i_mv_p = kk
        movie_vars(kk) = i_mv_p

     else if (movie_vars_txt(kk) .eq. '|v|') then
        if(i_mv_dens .eq. -1) i_mv_speed = kk
        movie_vars(kk) = i_mv_speed

    else if (movie_vars_txt(kk) .eq. 'vx') then
        if(i_mv_vx .eq. -1) i_mv_vx = kk
        movie_vars(kk) = i_mv_vx

    else if (movie_vars_txt(kk) .eq. 'vy') then
        if(i_mv_vy .eq. -1) i_mv_vy = kk
        movie_vars(kk) = i_mv_vy

    else if (movie_vars_txt(kk) .eq. 'vz') then
        if(i_mv_vz .eq. -1) i_mv_vz = kk
        movie_vars(kk) = i_mv_vz

     else if (movie_vars_txt(kk) .eq. 'Z') then
        if(i_mv_metallicity .eq. -1) i_mv_metallicity = kk
        movie_vars(kk) = i_mv_metallicity

     else if (movie_vars_txt(kk) .eq. 'pmag') then
        if(i_mv_pmag .eq. -1) i_mv_pmag = kk
        movie_vars(kk) = i_mv_pmag

     else if (movie_vars_txt(kk) .eq. 'xH2') then
        if(i_mv_xh2 .eq. -1) i_mv_xh2 = kk
        movie_vars(kk) = i_mv_xh2

     else if (movie_vars_txt(kk) .eq. 'xHI') then
        if(i_mv_xhi .eq. -1) i_mv_xhi = kk
        movie_vars(kk) = i_mv_xhi

     else if (movie_vars_txt(kk) .eq. 'xHII') then
        if(i_mv_xhii .eq. -1) i_mv_xhii = kk
        movie_vars(kk) = i_mv_xhii

     else if (movie_vars_txt(kk) .eq. 'xHeII') then
        if(i_mv_xheii .eq. -1) i_mv_xheii = kk
        movie_vars(kk) = i_mv_xheii

     else if (movie_vars_txt(kk) .eq. 'xHeIII') then
        if(i_mv_xheiii .eq. -1) i_mv_xheiii = kk
        movie_vars(kk) = i_mv_xheiii

     else if (index(movie_vars_txt(kk), 'var') .eq. 1) then
        if(i_mv_var .eq. -1) i_mv_var = kk
        movie_vars(kk) = i_mv_var
        ! Find which numbered variable to show
        read( movie_vars_txt(kk)(4:5), '(i1)' ) ivar
        movie_var_number(kk) = ivar

     else if (index(movie_vars_txt(kk), 'Fp') .eq. 1) then
        if(i_mv_fp .eq. -1) i_mv_fp = kk
        movie_vars(kk) = i_mv_fp
        ! Find which photon group to show
        read( movie_vars_txt(kk)(3:4), '(i1)' ) ivar
        movie_var_number(kk) = ivar

     endif

  end do

end subroutine set_movie_vars
