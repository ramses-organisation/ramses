subroutine synchro_fine(ilevel)
  use pm_commons
  use amr_commons
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::info
#endif
  integer::ilevel,xtondim
  !--------------------------------------------------------------------
  ! This routine synchronizes particle velocity with particle
  ! position for ilevel particle only. If particle sits entirely
  ! in level ilevel, then use inverse CIC at fine level to compute
  ! the force. Otherwise, use coarse level force and coarse level CIC.
  !--------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart
  integer::ig,ip,npart1,isink
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  if(sink)then
     fsink_new=0
  endif

#ifdef TSC
    xtondim=threetondim
#else
    xtondim=twotondim
#endif

  ! Synchronize velocity using CIC
  ig=0
  ip=0
  ! Loop over grids
  igrid=headl(myid,ilevel)
  do jgrid=1,numbl(myid,ilevel)
     npart1=numbp(igrid)  ! Number of particles in the grid
     if(npart1>0)then
        ig=ig+1
        ind_grid(ig)=igrid
        ipart=headp(igrid)
        ! Loop over particles
        do jpart=1,npart1
           if(ig==0)then
              ig=1
              ind_grid(ig)=igrid
           end if
           ip=ip+1
           ind_part(ip)=ipart
           ind_grid_part(ip)=ig
           if(ip==nvector)then
              call sync(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel,xtondim)
              ip=0
              ig=0
           end if
           ipart=nextp(ipart)  ! Go to next particle
        end do
        ! End loop over particles
     end if
     igrid=next(igrid)   ! Go to next grid
  end do
  ! End loop over grids
  if(ip>0)call sync(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel,xtondim)

  !sink cloud particles are used to average the grav. acceleration
  if(sink)then
     if(nsink>0)then
#ifndef WITHOUTMPI
        call MPI_ALLREDUCE(fsink_new,fsink_all,nsinkmax*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#else
        fsink_all=fsink_new
#endif
     endif
     do isink=1,nsink
        if (.not. direct_force_sink(isink))then
           fsink_partial(isink,1:ndim,ilevel)=fsink_all(isink,1:ndim)
        end if
     end do
  endif

111 format('   Entering synchro_fine for level ',I2)

end subroutine synchro_fine
!####################################################################
!####################################################################
!####################################################################
!####################################################################
subroutine synchro_fine_static(ilevel)
  use pm_commons
  use amr_commons
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::info
#endif
  integer::ilevel,xtondim
  !--------------------------------------------------------------------
  ! This routine synchronizes particle velocity with particle
  ! position for ilevel particle only. If particle sits entirely
  ! in level ilevel, then use inverse CIC at fine level to compute
  ! the force. Otherwise, use coarse level force and coarse level CIC.
  !--------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart
  integer::ig,ip,next_part,npart1,npart2,isink
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  if(sink)then
     fsink_new=0
     fsink_all=0
  endif

#ifdef TSC
    xtondim=threetondim
#else
    xtondim=twotondim
#endif

  ! Synchronize velocity using CIC
  ig=0
  ip=0
  ! Loop over grids
  igrid=headl(myid,ilevel)
  do jgrid=1,numbl(myid,ilevel)
     npart1=numbp(igrid)  ! Number of particles in the grid
     npart2=0

     ! Count particles
     if(npart1>0)then
        ipart=headp(igrid)
        ! Loop over particles
        do jpart=1,npart1
           ! Save next particle   <--- Very important !!!
           next_part=nextp(ipart)
           if(star) then
              if ( (.not. static_DM .and. is_DM(typep(ipart))) .or. &
                   & (.not. static_stars .and. (is_star(typep(ipart)) .or. is_debris(typep(ipart))) )  ) then
                 ! FIXME: there should be a static_sink as well
                 npart2=npart2+1
              endif
           else
              if(.not.static_dm) then
                 npart2=npart2+1
              endif
           endif
           ipart=next_part  ! Go to next particle
        end do
     endif

     ! Gather star particles
     if(npart2>0)then
        ig=ig+1
        ind_grid(ig)=igrid
        ipart=headp(igrid)
        ! Loop over particles
        do jpart=1,npart1
           ! Save next particle   <--- Very important !!!
           next_part=nextp(ipart)
           ! Select particles
           if(star) then
              if ( (.not. static_DM .and. is_DM(typep(ipart))) .or. &
                   & (.not. static_stars .and. (is_star(typep(ipart)) .or. is_debris(typep(ipart))) )  ) then
                 ! FIXME: what about sinks?
                 if(ig==0)then
                    ig=1
                    ind_grid(ig)=igrid
                 end if
                 ip=ip+1
                 ind_part(ip)=ipart
                 ind_grid_part(ip)=ig
              endif
           else
              if(.not.static_dm) then
                 if(ig==0)then
                    ig=1
                    ind_grid(ig)=igrid
                 end if
                 ip=ip+1
                 ind_part(ip)=ipart
                 ind_grid_part(ip)=ig
              endif
           endif
           if(ip==nvector)then
              call sync(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel,xtondim)
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
  if(ip>0)call sync(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel,xtondim)

  !sink cloud particles are used to average the grav. acceleration
  if(sink)then
     if(nsink>0)then
#ifndef WITHOUTMPI
        call MPI_ALLREDUCE(fsink_new,fsink_all,nsinkmax*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#else
        fsink_all=fsink_new
#endif
     endif
     do isink=1,nsink
        if (.not. direct_force_sink(isink))then
           fsink_partial(isink,1:ndim,ilevel)=fsink_all(isink,1:ndim)
        end if
     end do
  endif

111 format('   Entering synchro_fine for level ',I2)

end subroutine synchro_fine_static
!####################################################################
!####################################################################
!####################################################################
!####################################################################
subroutine sync(ind_grid,ind_part,ind_grid_part,ng,np,ilevel,xtondim)
  use amr_commons
  use pm_commons
  use poisson_commons
  implicit none
  integer::ng,np,ilevel,xtondim
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !
  !
  !
  logical::error
  integer::i,j,ind,idim,nx_loc,isink
  real(dp)::dx,scale
  ! Grid-based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  logical ,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector),save::dteff
  real(dp),dimension(1:3)::skip_loc
  ! Particle-based arrays
#ifndef TSC
  real(dp),dimension(1:nvector,1:ndim),save::x,ff,new_vp,dd,dg
  integer ,dimension(1:nvector,1:ndim),save::ig,id,igg,igd,icg,icd
  real(dp),dimension(1:nvector,1:twotondim),save::vol
  integer ,dimension(1:nvector,1:twotondim),save::igrid,icell,indp,kg
#else
  real(dp),dimension(1:nvector,1:ndim),save::x,ff,new_vp
  real(dp),dimension(1:nvector,1:ndim),save::cl,cr,cc,wl,wr,wc
  integer ,dimension(1:nvector,1:ndim),save::igl,igr,igc,icl,icr,icc
  real(dp),dimension(1:nvector,1:threetondim),save::vol
  integer ,dimension(1:nvector,1:threetondim),save::igrid,icell,indp,kg
#endif


  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)

  ! Lower left corner of 3x3x3 grid-cube
  do idim=1,ndim
     do i=1,ng
        x0(i,idim)=xg(ind_grid(i),idim)-3.0D0*dx
     end do
  end do

  ! Gather 27 neighboring father cells (should be present anytime !)
  do i=1,ng
     ind_cell(i)=father(ind_grid(i))
  end do
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ng,ilevel)

  ! Rescale position at level ilevel
  do idim=1,ndim
     do j=1,np
        x(j,idim)=xp(ind_part(j),idim)/scale+skip_loc(idim)
     end do
  end do
  do idim=1,ndim
     do j=1,np
        x(j,idim)=x(j,idim)-x0(ind_grid_part(j),idim)
     end do
  end do
  do idim=1,ndim
     do j=1,np
        x(j,idim)=x(j,idim)/dx
     end do
  end do
#ifdef TSC
    if (ndim .ne. 3)then
       write(*,*)'TSC not supported for ndim neq 3'
       call clean_stop
    end if

    xtondim=threetondim

    ! Check for illegal moves
    error=.false.
    do idim=1,ndim
       do j=1,np
          if(x(j,idim)<1.0D0.or.x(j,idim)>5.0D0)error=.true.
       end do
    end do
    if(error)then
       write(*,*)'problem in tsc_fine'
       do idim=1,ndim
          do j=1,np
             if(x(j,idim)<1.0D0.or.x(j,idim)>5.0D0)then
                write(*,*)x(j,1:ndim)
             endif
          end do
       end do
       stop
    end if

    ! TSC at level ilevel; a particle contributes
    !     to three cells in each dimension
    ! cl: position of leftmost cell centre
    ! cc: position of central cell centre
    ! cr: position of rightmost cell centre
    ! wl: weighting function for leftmost cell
    ! wc: weighting function for central cell
    ! wr: weighting function for rightmost cell
    do idim=1,ndim
       do j=1,np
          cl(j,idim)=dble(int(x(j,idim)))-0.5D0
          cc(j,idim)=dble(int(x(j,idim)))+0.5D0
          cr(j,idim)=dble(int(x(j,idim)))+1.5D0
          wl(j,idim)=0.50D0*(1.5D0-abs(x(j,idim)-cl(j,idim)))**2
          wc(j,idim)=0.75D0-          (x(j,idim)-cc(j,idim)) **2
          wr(j,idim)=0.50D0*(1.5D0-abs(x(j,idim)-cr(j,idim)))**2
       end do
    end do

    ! Compute parent grids
    do idim=1,ndim
       do j=1,np
          igl(j,idim)=(int(cl(j,idim)))/2
          igc(j,idim)=(int(cc(j,idim)))/2
          igr(j,idim)=(int(cr(j,idim)))/2
       end do
    end do
    ! #if NDIM==3
    do j=1,np
       kg(j,1 )=1+igl(j,1)+3*igl(j,2)+9*igl(j,3)
       kg(j,2 )=1+igc(j,1)+3*igl(j,2)+9*igl(j,3)
       kg(j,3 )=1+igr(j,1)+3*igl(j,2)+9*igl(j,3)
       kg(j,4 )=1+igl(j,1)+3*igc(j,2)+9*igl(j,3)
       kg(j,5 )=1+igc(j,1)+3*igc(j,2)+9*igl(j,3)
       kg(j,6 )=1+igr(j,1)+3*igc(j,2)+9*igl(j,3)
       kg(j,7 )=1+igl(j,1)+3*igr(j,2)+9*igl(j,3)
       kg(j,8 )=1+igc(j,1)+3*igr(j,2)+9*igl(j,3)
       kg(j,9 )=1+igr(j,1)+3*igr(j,2)+9*igl(j,3)
       kg(j,10)=1+igl(j,1)+3*igl(j,2)+9*igc(j,3)
       kg(j,11)=1+igc(j,1)+3*igl(j,2)+9*igc(j,3)
       kg(j,12)=1+igr(j,1)+3*igl(j,2)+9*igc(j,3)
       kg(j,13)=1+igl(j,1)+3*igc(j,2)+9*igc(j,3)
       kg(j,14)=1+igc(j,1)+3*igc(j,2)+9*igc(j,3)
       kg(j,15)=1+igr(j,1)+3*igc(j,2)+9*igc(j,3)
       kg(j,16)=1+igl(j,1)+3*igr(j,2)+9*igc(j,3)
       kg(j,17)=1+igc(j,1)+3*igr(j,2)+9*igc(j,3)
       kg(j,18)=1+igr(j,1)+3*igr(j,2)+9*igc(j,3)
       kg(j,19)=1+igl(j,1)+3*igl(j,2)+9*igr(j,3)
       kg(j,20)=1+igc(j,1)+3*igl(j,2)+9*igr(j,3)
       kg(j,21)=1+igr(j,1)+3*igl(j,2)+9*igr(j,3)
       kg(j,22)=1+igl(j,1)+3*igc(j,2)+9*igr(j,3)
       kg(j,23)=1+igc(j,1)+3*igc(j,2)+9*igr(j,3)
       kg(j,24)=1+igr(j,1)+3*igc(j,2)+9*igr(j,3)
       kg(j,25)=1+igl(j,1)+3*igr(j,2)+9*igr(j,3)
       kg(j,26)=1+igc(j,1)+3*igr(j,2)+9*igr(j,3)
       kg(j,27)=1+igr(j,1)+3*igr(j,2)+9*igr(j,3)
    end do

    do ind=1,threetondim
       do j=1,np
          igrid(j,ind)=son(nbors_father_cells(ind_grid_part(j),kg(j,ind)))
       end do
    end do

    ! Check if particles are entirely in level ilevel
    ok(1:np)=.true.
    do ind=1,threetondim
       do j=1,np
          ok(j)=ok(j).and.igrid(j,ind)>0
       end do
    end do

    ! If not, rescale position at level ilevel-1
    do idim=1,ndim
       do j=1,np
          if(.not.ok(j))then
             x(j,idim)=x(j,idim)/2.0D0
          end if
       end do
    end do
    ! If not, redo TSC at level ilevel-1
    do idim=1,ndim
       do j=1,np
          if(.not.ok(j))then
            cl(j,idim)=dble(int(x(j,idim)))-0.5D0
            cc(j,idim)=dble(int(x(j,idim)))+0.5D0
            cr(j,idim)=dble(int(x(j,idim)))+1.5D0
            wl(j,idim)=0.50D0*(1.5D0-abs(x(j,idim)-cl(j,idim)))**2
            wc(j,idim)=0.75D0-          (x(j,idim)-cc(j,idim)) **2
            wr(j,idim)=0.50D0*(1.5D0-abs(x(j,idim)-cr(j,idim)))**2
          end if
       end do
    end do

    ! Compute parent cell position
    do idim=1,ndim
       do j=1,np
         if(ok(j))then
          icl(j,idim)=int(cl(j,idim))-2*igl(j,idim)
          icc(j,idim)=int(cc(j,idim))-2*igc(j,idim)
          icr(j,idim)=int(cr(j,idim))-2*igr(j,idim)
         else ! ERM: this else may or may not be correct? But I believe it is.
          icl(j,idim)=int(cl(j,idim))
          icc(j,idim)=int(cc(j,idim))
          icr(j,idim)=int(cr(j,idim))
         endif
       end do
    end do

    ! #if NDIM==3
    do j=1,np
      if(ok(j))then
       icell(j,1 )=1+icl(j,1)+2*icl(j,2)+4*icl(j,3)
       icell(j,2 )=1+icc(j,1)+2*icl(j,2)+4*icl(j,3)
       icell(j,3 )=1+icr(j,1)+2*icl(j,2)+4*icl(j,3)
       icell(j,4 )=1+icl(j,1)+2*icc(j,2)+4*icl(j,3)
       icell(j,5 )=1+icc(j,1)+2*icc(j,2)+4*icl(j,3)
       icell(j,6 )=1+icr(j,1)+2*icc(j,2)+4*icl(j,3)
       icell(j,7 )=1+icl(j,1)+2*icr(j,2)+4*icl(j,3)
       icell(j,8 )=1+icc(j,1)+2*icr(j,2)+4*icl(j,3)
       icell(j,9 )=1+icr(j,1)+2*icr(j,2)+4*icl(j,3)
       icell(j,10)=1+icl(j,1)+2*icl(j,2)+4*icc(j,3)
       icell(j,11)=1+icc(j,1)+2*icl(j,2)+4*icc(j,3)
       icell(j,12)=1+icr(j,1)+2*icl(j,2)+4*icc(j,3)
       icell(j,13)=1+icl(j,1)+2*icc(j,2)+4*icc(j,3)
       icell(j,14)=1+icc(j,1)+2*icc(j,2)+4*icc(j,3)
       icell(j,15)=1+icr(j,1)+2*icc(j,2)+4*icc(j,3)
       icell(j,16)=1+icl(j,1)+2*icr(j,2)+4*icc(j,3)
       icell(j,17)=1+icc(j,1)+2*icr(j,2)+4*icc(j,3)
       icell(j,18)=1+icr(j,1)+2*icr(j,2)+4*icc(j,3)
       icell(j,19)=1+icl(j,1)+2*icl(j,2)+4*icr(j,3)
       icell(j,20)=1+icc(j,1)+2*icl(j,2)+4*icr(j,3)
       icell(j,21)=1+icr(j,1)+2*icl(j,2)+4*icr(j,3)
       icell(j,22)=1+icl(j,1)+2*icc(j,2)+4*icr(j,3)
       icell(j,23)=1+icc(j,1)+2*icc(j,2)+4*icr(j,3)
       icell(j,24)=1+icr(j,1)+2*icc(j,2)+4*icr(j,3)
       icell(j,25)=1+icl(j,1)+2*icr(j,2)+4*icr(j,3)
       icell(j,26)=1+icc(j,1)+2*icr(j,2)+4*icr(j,3)
       icell(j,27)=1+icr(j,1)+2*icr(j,2)+4*icr(j,3)
     else
       icell(j,1 )=1+icl(j,1)+3*icl(j,2)+9*icl(j,3)
       icell(j,2 )=1+icc(j,1)+3*icl(j,2)+9*icl(j,3)
       icell(j,3 )=1+icr(j,1)+3*icl(j,2)+9*icl(j,3)
       icell(j,4 )=1+icl(j,1)+3*icc(j,2)+9*icl(j,3)
       icell(j,5 )=1+icc(j,1)+3*icc(j,2)+9*icl(j,3)
       icell(j,6 )=1+icr(j,1)+3*icc(j,2)+9*icl(j,3)
       icell(j,7 )=1+icl(j,1)+3*icr(j,2)+9*icl(j,3)
       icell(j,8 )=1+icc(j,1)+3*icr(j,2)+9*icl(j,3)
       icell(j,9 )=1+icr(j,1)+3*icr(j,2)+9*icl(j,3)
       icell(j,10)=1+icl(j,1)+3*icl(j,2)+9*icc(j,3)
       icell(j,11)=1+icc(j,1)+3*icl(j,2)+9*icc(j,3)
       icell(j,12)=1+icr(j,1)+3*icl(j,2)+9*icc(j,3)
       icell(j,13)=1+icl(j,1)+3*icc(j,2)+9*icc(j,3)
       icell(j,14)=1+icc(j,1)+3*icc(j,2)+9*icc(j,3)
       icell(j,15)=1+icr(j,1)+3*icc(j,2)+9*icc(j,3)
       icell(j,16)=1+icl(j,1)+3*icr(j,2)+9*icc(j,3)
       icell(j,17)=1+icc(j,1)+3*icr(j,2)+9*icc(j,3)
       icell(j,18)=1+icr(j,1)+3*icr(j,2)+9*icc(j,3)
       icell(j,19)=1+icl(j,1)+3*icl(j,2)+9*icr(j,3)
       icell(j,20)=1+icc(j,1)+3*icl(j,2)+9*icr(j,3)
       icell(j,21)=1+icr(j,1)+3*icl(j,2)+9*icr(j,3)
       icell(j,22)=1+icl(j,1)+3*icc(j,2)+9*icr(j,3)
       icell(j,23)=1+icc(j,1)+3*icc(j,2)+9*icr(j,3)
       icell(j,24)=1+icr(j,1)+3*icc(j,2)+9*icr(j,3)
       icell(j,25)=1+icl(j,1)+3*icr(j,2)+9*icr(j,3)
       icell(j,26)=1+icc(j,1)+3*icr(j,2)+9*icr(j,3)
       icell(j,27)=1+icr(j,1)+3*icr(j,2)+9*icr(j,3)
     endif
    end do

    ! Compute parent cell adress
    do ind=1,threetondim
       do j=1,np
         if(ok(j))then
          indp(j,ind)=ncoarse+(icell(j,ind)-1)*ngridmax+igrid(j,ind)
         else ! ERM: for AMR(?) there may be an issue with ind_grid_part(j) being used here.
           indp(j,ind)=nbors_father_cells(ind_grid_part(j),icell(j,ind))
         endif
       end do
    end do

    ! Compute cloud volumes (NDIM==3)
    do j=1,np
       vol(j,1 )=wl(j,1)*wl(j,2)*wl(j,3)
       vol(j,2 )=wc(j,1)*wl(j,2)*wl(j,3)
       vol(j,3 )=wr(j,1)*wl(j,2)*wl(j,3)
       vol(j,4 )=wl(j,1)*wc(j,2)*wl(j,3)
       vol(j,5 )=wc(j,1)*wc(j,2)*wl(j,3)
       vol(j,6 )=wr(j,1)*wc(j,2)*wl(j,3)
       vol(j,7 )=wl(j,1)*wr(j,2)*wl(j,3)
       vol(j,8 )=wc(j,1)*wr(j,2)*wl(j,3)
       vol(j,9 )=wr(j,1)*wr(j,2)*wl(j,3)
       vol(j,10)=wl(j,1)*wl(j,2)*wc(j,3)
       vol(j,11)=wc(j,1)*wl(j,2)*wc(j,3)
       vol(j,12)=wr(j,1)*wl(j,2)*wc(j,3)
       vol(j,13)=wl(j,1)*wc(j,2)*wc(j,3)
       vol(j,14)=wc(j,1)*wc(j,2)*wc(j,3)
       vol(j,15)=wr(j,1)*wc(j,2)*wc(j,3)
       vol(j,16)=wl(j,1)*wr(j,2)*wc(j,3)
       vol(j,17)=wc(j,1)*wr(j,2)*wc(j,3)
       vol(j,18)=wr(j,1)*wr(j,2)*wc(j,3)
       vol(j,19)=wl(j,1)*wl(j,2)*wr(j,3)
       vol(j,20)=wc(j,1)*wl(j,2)*wr(j,3)
       vol(j,21)=wr(j,1)*wl(j,2)*wr(j,3)
       vol(j,22)=wl(j,1)*wc(j,2)*wr(j,3)
       vol(j,23)=wc(j,1)*wc(j,2)*wr(j,3)
       vol(j,24)=wr(j,1)*wc(j,2)*wr(j,3)
       vol(j,25)=wl(j,1)*wr(j,2)*wr(j,3)
       vol(j,26)=wc(j,1)*wr(j,2)*wr(j,3)
       vol(j,27)=wr(j,1)*wr(j,2)*wr(j,3)
    end do

#else

  ! Check for illegal moves
  error=.false.
  do idim=1,ndim
     do j=1,np
        if(x(j,idim)<0.5D0.or.x(j,idim)>5.5D0)error=.true.
     end do
  end do
  if(error)then
     write(*,*)'problem in sync'
     do idim=1,ndim
        do j=1,np
           if(x(j,idim)<0.5D0.or.x(j,idim)>5.5D0)then
              write(*,*)x(j,1:ndim)
           endif
        end do
     end do
     stop
  end if

  ! CIC at level ilevel (dd: right cloud boundary; dg: left cloud boundary)
  do idim=1,ndim
     do j=1,np
        dd(j,idim)=x(j,idim)+0.5D0
        id(j,idim)=int(dd(j,idim))
        dd(j,idim)=dd(j,idim)-id(j,idim)
        dg(j,idim)=1.0D0-dd(j,idim)
        ig(j,idim)=id(j,idim)-1
     end do
  end do

   ! Compute parent grids
  do idim=1,ndim
     do j=1,np
        igg(j,idim)=ig(j,idim)/2
        igd(j,idim)=id(j,idim)/2
     end do
  end do
#if NDIM==1
  do j=1,np
     kg(j,1)=1+igg(j,1)
     kg(j,2)=1+igd(j,1)
  end do
#endif
#if NDIM==2
  do j=1,np
     kg(j,1)=1+igg(j,1)+3*igg(j,2)
     kg(j,2)=1+igd(j,1)+3*igg(j,2)
     kg(j,3)=1+igg(j,1)+3*igd(j,2)
     kg(j,4)=1+igd(j,1)+3*igd(j,2)
  end do
#endif
#if NDIM==3
  do j=1,np
     kg(j,1)=1+igg(j,1)+3*igg(j,2)+9*igg(j,3)
     kg(j,2)=1+igd(j,1)+3*igg(j,2)+9*igg(j,3)
     kg(j,3)=1+igg(j,1)+3*igd(j,2)+9*igg(j,3)
     kg(j,4)=1+igd(j,1)+3*igd(j,2)+9*igg(j,3)
     kg(j,5)=1+igg(j,1)+3*igg(j,2)+9*igd(j,3)
     kg(j,6)=1+igd(j,1)+3*igg(j,2)+9*igd(j,3)
     kg(j,7)=1+igg(j,1)+3*igd(j,2)+9*igd(j,3)
     kg(j,8)=1+igd(j,1)+3*igd(j,2)+9*igd(j,3)
  end do
#endif
  do ind=1,twotondim
     do j=1,np
        igrid(j,ind)=son(nbors_father_cells(ind_grid_part(j),kg(j,ind)))
     end do
  end do

  ! Check if particles are entirely in level ilevel
  ok(1:np)=.true.
  do ind=1,twotondim
     do j=1,np
        ok(j)=ok(j).and.igrid(j,ind)>0
     end do
  end do

  ! If not, rescale position at level ilevel-1
  do idim=1,ndim
     do j=1,np
        if(.not.ok(j))then
           x(j,idim)=x(j,idim)/2.0D0
        end if
     end do
  end do
  ! If not, redo CIC at level ilevel-1
  do idim=1,ndim
     do j=1,np
        if(.not.ok(j))then
           dd(j,idim)=x(j,idim)+0.5D0
           id(j,idim)=int(dd(j,idim))
           dd(j,idim)=dd(j,idim)-id(j,idim)
           dg(j,idim)=1.0D0-dd(j,idim)
           ig(j,idim)=id(j,idim)-1
        end if
     end do
  end do

 ! Compute parent cell position
  do idim=1,ndim
     do j=1,np
        if(ok(j))then
           icg(j,idim)=ig(j,idim)-2*igg(j,idim)
           icd(j,idim)=id(j,idim)-2*igd(j,idim)
        else
           icg(j,idim)=ig(j,idim)
           icd(j,idim)=id(j,idim)
        end if
     end do
  end do
#if NDIM==1
  do j=1,np
     icell(j,1)=1+icg(j,1)
     icell(j,2)=1+icd(j,1)
  end do
#endif
#if NDIM==2
  do j=1,np
     if(ok(j))then
        icell(j,1)=1+icg(j,1)+2*icg(j,2)
        icell(j,2)=1+icd(j,1)+2*icg(j,2)
        icell(j,3)=1+icg(j,1)+2*icd(j,2)
        icell(j,4)=1+icd(j,1)+2*icd(j,2)
     else
        icell(j,1)=1+icg(j,1)+3*icg(j,2)
        icell(j,2)=1+icd(j,1)+3*icg(j,2)
        icell(j,3)=1+icg(j,1)+3*icd(j,2)
        icell(j,4)=1+icd(j,1)+3*icd(j,2)
     end if
  end do
#endif
#if NDIM==3
  do j=1,np
     if(ok(j))then
        icell(j,1)=1+icg(j,1)+2*icg(j,2)+4*icg(j,3)
        icell(j,2)=1+icd(j,1)+2*icg(j,2)+4*icg(j,3)
        icell(j,3)=1+icg(j,1)+2*icd(j,2)+4*icg(j,3)
        icell(j,4)=1+icd(j,1)+2*icd(j,2)+4*icg(j,3)
        icell(j,5)=1+icg(j,1)+2*icg(j,2)+4*icd(j,3)
        icell(j,6)=1+icd(j,1)+2*icg(j,2)+4*icd(j,3)
        icell(j,7)=1+icg(j,1)+2*icd(j,2)+4*icd(j,3)
        icell(j,8)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
     else
        icell(j,1)=1+icg(j,1)+3*icg(j,2)+9*icg(j,3)
        icell(j,2)=1+icd(j,1)+3*icg(j,2)+9*icg(j,3)
        icell(j,3)=1+icg(j,1)+3*icd(j,2)+9*icg(j,3)
        icell(j,4)=1+icd(j,1)+3*icd(j,2)+9*icg(j,3)
        icell(j,5)=1+icg(j,1)+3*icg(j,2)+9*icd(j,3)
        icell(j,6)=1+icd(j,1)+3*icg(j,2)+9*icd(j,3)
        icell(j,7)=1+icg(j,1)+3*icd(j,2)+9*icd(j,3)
        icell(j,8)=1+icd(j,1)+3*icd(j,2)+9*icd(j,3)
     end if
  end do
#endif

  ! Compute parent cell adresses
  do ind=1,twotondim
     do j=1,np
        if(ok(j))then
           indp(j,ind)=ncoarse+(icell(j,ind)-1)*ngridmax+igrid(j,ind)
        else
           indp(j,ind)=nbors_father_cells(ind_grid_part(j),icell(j,ind))
        end if
     end do
  end do

  ! Compute cloud volumes
#if NDIM==1
  do j=1,np
     vol(j,1)=dg(j,1)
     vol(j,2)=dd(j,1)
  end do
#endif
#if NDIM==2
  do j=1,np
     vol(j,1)=dg(j,1)*dg(j,2)
     vol(j,2)=dd(j,1)*dg(j,2)
     vol(j,3)=dg(j,1)*dd(j,2)
     vol(j,4)=dd(j,1)*dd(j,2)
  end do
#endif
#if NDIM==3
  do j=1,np
     vol(j,1)=dg(j,1)*dg(j,2)*dg(j,3)
     vol(j,2)=dd(j,1)*dg(j,2)*dg(j,3)
     vol(j,3)=dg(j,1)*dd(j,2)*dg(j,3)
     vol(j,4)=dd(j,1)*dd(j,2)*dg(j,3)
     vol(j,5)=dg(j,1)*dg(j,2)*dd(j,3)
     vol(j,6)=dd(j,1)*dg(j,2)*dd(j,3)
     vol(j,7)=dg(j,1)*dd(j,2)*dd(j,3)
     vol(j,8)=dd(j,1)*dd(j,2)*dd(j,3)
  end do
#endif

!#else
!#include "tsc_fine.F90"
#endif
  ! Gather 3-force
  ff(1:np,1:ndim)=0.0D0
  if(poisson)then
     do ind=1,xtondim
        do idim=1,ndim
           do j=1,np
              ff(j,idim)=ff(j,idim)+f(indp(j,ind),idim)*vol(j,ind)
           end do
        end do
     end do
  endif

  ! For sink particle only, store contribution to the sink force
  if(sink)then
     do idim=1,ndim
        do j=1,np
           if ( is_cloud(typep(ind_part(j))) ) then
              isink=-idp(ind_part(j))
              if(.not. direct_force_sink(isink))then
                 fsink_new(isink,idim)=fsink_new(isink,idim)+ff(j,idim)
              endif
           endif
        end do
     end do
  end if

  ! Compute individual time steps
  do j=1,np
     if(levelp(ind_part(j))>=ilevel)then
        dteff(j)=dtnew(levelp(ind_part(j)))
     else
        dteff(j)=dtold(levelp(ind_part(j)))
     endif
  end do

  ! Update particles level
  do j=1,np
     levelp(ind_part(j))=ilevel
  end do

  ! Update 3-velocity
  do idim=1,ndim
     if(static)then
        do j=1,np
           new_vp(j,idim)=ff(j,idim)
        end do
     else
        do j=1,np
           new_vp(j,idim)=vp(ind_part(j),idim)+ff(j,idim)*0.5D0*dteff(j)
        end do
     endif
  end do
  do idim=1,ndim
     do j=1,np
        vp(ind_part(j),idim)=new_vp(j,idim)
     end do
  end do

  ! Extra acceleration forces are added here.
  if((accel_gr(1).ne.0).or.(accel_gr(2).ne.0).or.(accel_gr(3).ne.0))then
     do idim=1,ndim
        do j=1,np
           vp(ind_part(j),idim)=vp(ind_part(j),idim)+dteff(j)*accel_gr(idim)
        end do
     end do
  endif

  ! For sink particle only, overwrite cloud particle velocity with sink velocity
  if(sink)then
     do idim=1,ndim
        do j=1,np
           if ( is_cloud(typep(ind_part(j))) ) then
              isink=-idp(ind_part(j))
              ! Remember that vsink is half time step older than other particles
              vp(ind_part(j),idim)=vsink(isink,idim)
           endif
        end do
     end do
  end if

end subroutine sync
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine init_dust_fine(ilevel)
  use pm_commons
  use amr_commons
  use hydro_commons
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::info
#endif
  integer::ilevel,xtondim
  ! First, reset uold to zero.
  ! Can remove gravity and sink particle related things.
  ! Can remove synchro_fine_static as well.
  ! In "sync", want to remove the gravity...
  ! Syncing up the velocity, get rid of that too.
  !--------------------------------------------------------------------
  ! This routine synchronizes particle velocity with particle
  ! position for ilevel particle only. If particle sits entirely
  ! in level ilevel, then use inverse CIC at fine level to compute
  ! the force. Otherwise, use coarse level force and coarse level CIC.
  !--------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart
  integer::ig,ip,npart1,isink
  integer::i,iskip,icpu,ind,ibound,ivar,ivar_dust
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel
  ivar_dust=9

#ifdef TSC
    xtondim=threetondim
#else
    xtondim=twotondim
#endif

  ! Reset unew to zero for dust ``stopping time''
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
        do i=1,active(ilevel)%ngrid
           unew(active(ilevel)%igrid(i)+iskip,ivar_dust)=0.0D0
        end do
  end do
  do icpu=1,ncpu
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
           do i=1,reception(icpu,ilevel)%ngrid
              unew(reception(icpu,ilevel)%igrid(i)+iskip,ivar_dust)=0.0D0
           end do
     end do
  end do

  ! Reset uold to zero for dust mass and momentum densities
  do icpu=1,ncpu
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do ivar=ivar_dust,ivar_dust+3
           do i=1,reception(icpu,ilevel)%ngrid
              uold(reception(icpu,ilevel)%igrid(i)+iskip,ivar)=0.0D0
           end do
        end do
     end do
  end do
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do ivar=ivar_dust,ivar_dust+3
        do i=1,active(ilevel)%ngrid
           uold(active(ilevel)%igrid(i)+iskip,ivar)=0.0D0
        end do
     end do
  end do

!!!!!!!!!!!!!!!!!!!!!!! NEW !!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Reset dust variables in physical boundaries
   do ibound=1,nboundary
      do ind=1,twotondim
         iskip=ncoarse+(ind-1)*ngridmax
         do ivar=ivar_dust,ivar_dust+ndim
            do i=1,boundary(ibound,ilevel)%ngrid
               uold(boundary(ibound,ilevel)%igrid(i)+iskip,ivar)=0.0D0
               ! unew(boundary(ibound,ilevel)%igrid(i)+iskip,ivar)=&
               ! &uold(boundary(ibound,ilevel)%igrid(i)+iskip,ivar)
            end do
         end do
      end do
   end do

   do ibound=1,nboundary
      do ind=1,twotondim
         iskip=ncoarse+(ind-1)*ngridmax
         do ivar=ivar_dust,ivar_dust+ndim
            do i=1,boundary(ibound,ilevel)%ngrid
               unew(boundary(ibound,ilevel)%igrid(i)+iskip,ivar)=0.0D0
               ! unew(boundary(ibound,ilevel)%igrid(i)+iskip,ivar)=&
               ! &uold(boundary(ibound,ilevel)%igrid(i)+iskip,ivar)
            end do
         end do
      end do
   end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Synchronize velocity using CIC (No longer need velocity to be synced.)
  ig=0
  ip=0
  ! Loop over grids
  igrid=headl(myid,ilevel)
  do jgrid=1,numbl(myid,ilevel)
     npart1=numbp(igrid)  ! Number of particles in the grid
     if(npart1>0)then
        ig=ig+1
        ind_grid(ig)=igrid
        ipart=headp(igrid)
        ! Loop over particles
        do jpart=1,npart1
           if(ig==0)then
              ig=1
              ind_grid(ig)=igrid
           end if
           ip=ip+1
           ind_part(ip)=ipart
           ind_grid_part(ip)=ig
           if(ip==nvector)then
              call init_dust(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel,xtondim)
              ip=0
              ig=0
           end if
           ipart=nextp(ipart)  ! Go to next particle
        end do
        ! End loop over particles
     end if
     igrid=next(igrid)   ! Go to next grid
  end do
  ! End loop over grids
  if(ip>0)call init_dust(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel,xtondim)

  ! Update MPI boundary conditions for uold for dust mass and momentum densities
  do ivar=ivar_dust,ivar_dust+ndim
     call make_virtual_reverse_dp(uold(1,ivar),ilevel)
     call make_virtual_fine_dp   (uold(1,ivar),ilevel)
  end do

  ! Update MPI boundary conditions for unew for dust "mass" and "momentum" densities
  call make_virtual_reverse_dp(unew(1,ivar_dust),ilevel)
  call make_virtual_fine_dp   (unew(1,ivar_dust),ilevel)



111 format('   Entering init_dust_fine for level ',I2)

end subroutine init_dust_fine
!####################################################################
!####################################################################
!####################################################################
!####################################################################
subroutine init_dust(ind_grid,ind_part,ind_grid_part,ng,np,ilevel,xtondim)
  use amr_commons
  !use amr_parameters ERM
  use pm_commons
  use poisson_commons
  use hydro_commons, ONLY: uold,unew,smallr,nvar,gamma
  implicit none
  integer::ng,np,ilevel,xtondim
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !
  !
  !
  logical::error
  integer::i,j,ind,idim,nx_loc,isink,ivar_dust
  real(dp)::dx,scale,dx_loc,vol_loc
  real(dp)::ctm ! ERM: recommend 1.15D3
  real(dp)::ts !ERM: recommend 2.2D-1
  real(dp)::rd ! ERM: Grain size parameter

  ! Grid-based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle-based arrays
  logical ,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector),save::mmm,dteff,nu_stop
  real(dp),dimension(1:nvector),save::dgr,tss,mm ! ERM: density, (non-constant) stopping times
  real(dp),dimension(1:nvector,1:ndim),save::uu,bb,vv ! ERM: Added these arrays
#ifndef TSC
  real(dp),dimension(1:nvector,1:ndim),save::x,ff,new_vp,dd,dg
  integer ,dimension(1:nvector,1:ndim),save::ig,id,igg,igd,icg,icd
  real(dp),dimension(1:nvector,1:twotondim),save::vol
  integer ,dimension(1:nvector,1:twotondim),save::igrid,icell,indp,kg
#else
  real(dp),dimension(1:nvector,1:ndim),save::x,ff,new_vp,cl,cr,cc,wl,wr,wc
  integer ,dimension(1:nvector,1:ndim),save::igl,igr,igc,icl,icr,icc
  real(dp),dimension(1:nvector,1:threetondim),save::vol
  integer ,dimension(1:nvector,1:threetondim),save::igrid,icell,indp,kg
#endif
  real(dp),dimension(1:3)::skip_loc

  ctm = charge_to_mass
  ts = t_stop
  rd = 0.62665706865775*grain_size ! constant for epstein drag law. used to havesqrt(gamma)*

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

  ! Lower left corner of 3x3x3 grid-cube
  do idim=1,ndim
     do i=1,ng
        x0(i,idim)=xg(ind_grid(i),idim)-3.0D0*dx
     end do
  end do

  ! Gather 27 neighboring father cells (should be present anytime !)
  do i=1,ng
     ind_cell(i)=father(ind_grid(i))
  end do
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ng,ilevel)

  ! Rescale position at level ilevel
  do idim=1,ndim
     do j=1,np
        x(j,idim)=xp(ind_part(j),idim)/scale+skip_loc(idim)
     end do
  end do
  do idim=1,ndim
     do j=1,np
        x(j,idim)=x(j,idim)-x0(ind_grid_part(j),idim)
     end do
  end do
  do idim=1,ndim
     do j=1,np
        x(j,idim)=x(j,idim)/dx
     end do
  end do

#ifndef TSC
  ! Check for illegal moves
  error=.false.
  do idim=1,ndim
     do j=1,np
        if(x(j,idim)<0.5D0.or.x(j,idim)>5.5D0)error=.true.
     end do
  end do
  if(error)then
     write(*,*)'problem in init_dust'
     do idim=1,ndim
        do j=1,np
           if(x(j,idim)<0.5D0.or.x(j,idim)>5.5D0)then
              write(*,*)x(j,1:ndim)
           endif
        end do
     end do
     stop
  end if

  ! CIC at level ilevel (dd: right cloud boundary; dg: left cloud boundary)
  do idim=1,ndim
     do j=1,np
        dd(j,idim)=x(j,idim)+0.5D0
        id(j,idim)=int(dd(j,idim))
        dd(j,idim)=dd(j,idim)-id(j,idim)
        dg(j,idim)=1.0D0-dd(j,idim)
        ig(j,idim)=id(j,idim)-1
     end do
  end do

   ! Compute parent grids
  do idim=1,ndim
     do j=1,np
        igg(j,idim)=ig(j,idim)/2
        igd(j,idim)=id(j,idim)/2
     end do
  end do
#if NDIM==1
  do j=1,np
     kg(j,1)=1+igg(j,1)
     kg(j,2)=1+igd(j,1)
  end do
#endif
#if NDIM==2
  do j=1,np
     kg(j,1)=1+igg(j,1)+3*igg(j,2)
     kg(j,2)=1+igd(j,1)+3*igg(j,2)
     kg(j,3)=1+igg(j,1)+3*igd(j,2)
     kg(j,4)=1+igd(j,1)+3*igd(j,2)
  end do
#endif
#if NDIM==3
  do j=1,np
     kg(j,1)=1+igg(j,1)+3*igg(j,2)+9*igg(j,3)
     kg(j,2)=1+igd(j,1)+3*igg(j,2)+9*igg(j,3)
     kg(j,3)=1+igg(j,1)+3*igd(j,2)+9*igg(j,3)
     kg(j,4)=1+igd(j,1)+3*igd(j,2)+9*igg(j,3)
     kg(j,5)=1+igg(j,1)+3*igg(j,2)+9*igd(j,3)
     kg(j,6)=1+igd(j,1)+3*igg(j,2)+9*igd(j,3)
     kg(j,7)=1+igg(j,1)+3*igd(j,2)+9*igd(j,3)
     kg(j,8)=1+igd(j,1)+3*igd(j,2)+9*igd(j,3)
  end do
#endif
  do ind=1,twotondim
     do j=1,np
        igrid(j,ind)=son(nbors_father_cells(ind_grid_part(j),kg(j,ind)))
     end do
  end do

  ! Check if particles are entirely in level ilevel
  ok(1:np)=.true.
  do ind=1,twotondim
     do j=1,np
        ok(j)=ok(j).and.igrid(j,ind)>0
     end do
  end do

  ! If not, rescale position at level ilevel-1
  do idim=1,ndim
     do j=1,np
        if(.not.ok(j))then
           x(j,idim)=x(j,idim)/2.0D0
        end if
     end do
  end do
  ! If not, redo CIC at level ilevel-1
  do idim=1,ndim
     do j=1,np
        if(.not.ok(j))then
           dd(j,idim)=x(j,idim)+0.5D0
           id(j,idim)=int(dd(j,idim))
           dd(j,idim)=dd(j,idim)-id(j,idim)
           dg(j,idim)=1.0D0-dd(j,idim)
           ig(j,idim)=id(j,idim)-1
        end if
     end do
  end do

 ! Compute parent cell position
  do idim=1,ndim
     do j=1,np
        if(ok(j))then
           icg(j,idim)=ig(j,idim)-2*igg(j,idim)
           icd(j,idim)=id(j,idim)-2*igd(j,idim)
        else
           icg(j,idim)=ig(j,idim)
           icd(j,idim)=id(j,idim)
        end if
     end do
  end do
#if NDIM==1
  do j=1,np
     icell(j,1)=1+icg(j,1)
     icell(j,2)=1+icd(j,1)
  end do
#endif
#if NDIM==2
  do j=1,np
     if(ok(j))then
        icell(j,1)=1+icg(j,1)+2*icg(j,2)
        icell(j,2)=1+icd(j,1)+2*icg(j,2)
        icell(j,3)=1+icg(j,1)+2*icd(j,2)
        icell(j,4)=1+icd(j,1)+2*icd(j,2)
     else
        icell(j,1)=1+icg(j,1)+3*icg(j,2)
        icell(j,2)=1+icd(j,1)+3*icg(j,2)
        icell(j,3)=1+icg(j,1)+3*icd(j,2)
        icell(j,4)=1+icd(j,1)+3*icd(j,2)
     end if
  end do
#endif
#if NDIM==3
  do j=1,np
     if(ok(j))then
        icell(j,1)=1+icg(j,1)+2*icg(j,2)+4*icg(j,3)
        icell(j,2)=1+icd(j,1)+2*icg(j,2)+4*icg(j,3)
        icell(j,3)=1+icg(j,1)+2*icd(j,2)+4*icg(j,3)
        icell(j,4)=1+icd(j,1)+2*icd(j,2)+4*icg(j,3)
        icell(j,5)=1+icg(j,1)+2*icg(j,2)+4*icd(j,3)
        icell(j,6)=1+icd(j,1)+2*icg(j,2)+4*icd(j,3)
        icell(j,7)=1+icg(j,1)+2*icd(j,2)+4*icd(j,3)
        icell(j,8)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
     else
        icell(j,1)=1+icg(j,1)+3*icg(j,2)+9*icg(j,3)
        icell(j,2)=1+icd(j,1)+3*icg(j,2)+9*icg(j,3)
        icell(j,3)=1+icg(j,1)+3*icd(j,2)+9*icg(j,3)
        icell(j,4)=1+icd(j,1)+3*icd(j,2)+9*icg(j,3)
        icell(j,5)=1+icg(j,1)+3*icg(j,2)+9*icd(j,3)
        icell(j,6)=1+icd(j,1)+3*icg(j,2)+9*icd(j,3)
        icell(j,7)=1+icg(j,1)+3*icd(j,2)+9*icd(j,3)
        icell(j,8)=1+icd(j,1)+3*icd(j,2)+9*icd(j,3)
     end if
  end do
#endif

  ! Compute parent cell adresses
  do ind=1,twotondim
     do j=1,np
        if(ok(j))then
           indp(j,ind)=ncoarse+(icell(j,ind)-1)*ngridmax+igrid(j,ind)
        else
           indp(j,ind)=nbors_father_cells(ind_grid_part(j),icell(j,ind))
        end if
     end do
  end do

  ! Compute cloud volumes
#if NDIM==1
  do j=1,np
     vol(j,1)=dg(j,1)
     vol(j,2)=dd(j,1)
  end do
#endif
#if NDIM==2
  do j=1,np
     vol(j,1)=dg(j,1)*dg(j,2)
     vol(j,2)=dd(j,1)*dg(j,2)
     vol(j,3)=dg(j,1)*dd(j,2)
     vol(j,4)=dd(j,1)*dd(j,2)
  end do
#endif
#if NDIM==3
  do j=1,np
     vol(j,1)=dg(j,1)*dg(j,2)*dg(j,3)
     vol(j,2)=dd(j,1)*dg(j,2)*dg(j,3)
     vol(j,3)=dg(j,1)*dd(j,2)*dg(j,3)
     vol(j,4)=dd(j,1)*dd(j,2)*dg(j,3)
     vol(j,5)=dg(j,1)*dg(j,2)*dd(j,3)
     vol(j,6)=dd(j,1)*dg(j,2)*dd(j,3)
     vol(j,7)=dg(j,1)*dd(j,2)*dd(j,3)
     vol(j,8)=dd(j,1)*dd(j,2)*dd(j,3)
  end do
#endif

#else
!#include "tsc_fine.F90"
    if (ndim .ne. 3)then
       write(*,*)'TSC not supported for ndim neq 3'
       call clean_stop
    end if

    xtondim=threetondim

    ! Check for illegal moves
    error=.false.
    do idim=1,ndim
       do j=1,np
          if(x(j,idim)<1.0D0.or.x(j,idim)>5.0D0)error=.true.
       end do
    end do
    if(error)then
       write(*,*)'problem in tsc_fine'
       do idim=1,ndim
          do j=1,np
             if(x(j,idim)<1.0D0.or.x(j,idim)>5.0D0)then
                write(*,*)x(j,1:ndim)
             endif
          end do
       end do
       stop
    end if

    ! TSC at level ilevel; a particle contributes
    !     to three cells in each dimension
    ! cl: position of leftmost cell centre
    ! cc: position of central cell centre
    ! cr: position of rightmost cell centre
    ! wl: weighting function for leftmost cell
    ! wc: weighting function for central cell
    ! wr: weighting function for rightmost cell
    do idim=1,ndim
       do j=1,np
          cl(j,idim)=dble(int(x(j,idim)))-0.5D0
          cc(j,idim)=dble(int(x(j,idim)))+0.5D0
          cr(j,idim)=dble(int(x(j,idim)))+1.5D0
          wl(j,idim)=0.50D0*(1.5D0-abs(x(j,idim)-cl(j,idim)))**2
          wc(j,idim)=0.75D0-          (x(j,idim)-cc(j,idim)) **2
          wr(j,idim)=0.50D0*(1.5D0-abs(x(j,idim)-cr(j,idim)))**2
       end do
    end do

    ! Compute parent grids
    do idim=1,ndim
       do j=1,np
          igl(j,idim)=(int(cl(j,idim)))/2
          igc(j,idim)=(int(cc(j,idim)))/2
          igr(j,idim)=(int(cr(j,idim)))/2
       end do
    end do
    ! #if NDIM==3
    do j=1,np
       kg(j,1 )=1+igl(j,1)+3*igl(j,2)+9*igl(j,3)
       kg(j,2 )=1+igc(j,1)+3*igl(j,2)+9*igl(j,3)
       kg(j,3 )=1+igr(j,1)+3*igl(j,2)+9*igl(j,3)
       kg(j,4 )=1+igl(j,1)+3*igc(j,2)+9*igl(j,3)
       kg(j,5 )=1+igc(j,1)+3*igc(j,2)+9*igl(j,3)
       kg(j,6 )=1+igr(j,1)+3*igc(j,2)+9*igl(j,3)
       kg(j,7 )=1+igl(j,1)+3*igr(j,2)+9*igl(j,3)
       kg(j,8 )=1+igc(j,1)+3*igr(j,2)+9*igl(j,3)
       kg(j,9 )=1+igr(j,1)+3*igr(j,2)+9*igl(j,3)
       kg(j,10)=1+igl(j,1)+3*igl(j,2)+9*igc(j,3)
       kg(j,11)=1+igc(j,1)+3*igl(j,2)+9*igc(j,3)
       kg(j,12)=1+igr(j,1)+3*igl(j,2)+9*igc(j,3)
       kg(j,13)=1+igl(j,1)+3*igc(j,2)+9*igc(j,3)
       kg(j,14)=1+igc(j,1)+3*igc(j,2)+9*igc(j,3)
       kg(j,15)=1+igr(j,1)+3*igc(j,2)+9*igc(j,3)
       kg(j,16)=1+igl(j,1)+3*igr(j,2)+9*igc(j,3)
       kg(j,17)=1+igc(j,1)+3*igr(j,2)+9*igc(j,3)
       kg(j,18)=1+igr(j,1)+3*igr(j,2)+9*igc(j,3)
       kg(j,19)=1+igl(j,1)+3*igl(j,2)+9*igr(j,3)
       kg(j,20)=1+igc(j,1)+3*igl(j,2)+9*igr(j,3)
       kg(j,21)=1+igr(j,1)+3*igl(j,2)+9*igr(j,3)
       kg(j,22)=1+igl(j,1)+3*igc(j,2)+9*igr(j,3)
       kg(j,23)=1+igc(j,1)+3*igc(j,2)+9*igr(j,3)
       kg(j,24)=1+igr(j,1)+3*igc(j,2)+9*igr(j,3)
       kg(j,25)=1+igl(j,1)+3*igr(j,2)+9*igr(j,3)
       kg(j,26)=1+igc(j,1)+3*igr(j,2)+9*igr(j,3)
       kg(j,27)=1+igr(j,1)+3*igr(j,2)+9*igr(j,3)
    end do

    do ind=1,threetondim
       do j=1,np
          igrid(j,ind)=son(nbors_father_cells(ind_grid_part(j),kg(j,ind)))
       end do
    end do

    ! Check if particles are entirely in level ilevel
    ok(1:np)=.true.
    do ind=1,threetondim
       do j=1,np
          ok(j)=ok(j).and.igrid(j,ind)>0
       end do
    end do

    ! If not, rescale position at level ilevel-1
    do idim=1,ndim
       do j=1,np
          if(.not.ok(j))then
             x(j,idim)=x(j,idim)/2.0D0
          end if
       end do
    end do
    ! If not, redo TSC at level ilevel-1
    do idim=1,ndim
       do j=1,np
          if(.not.ok(j))then
            cl(j,idim)=dble(int(x(j,idim)))-0.5D0
            cc(j,idim)=dble(int(x(j,idim)))+0.5D0
            cr(j,idim)=dble(int(x(j,idim)))+1.5D0
            wl(j,idim)=0.50D0*(1.5D0-abs(x(j,idim)-cl(j,idim)))**2
            wc(j,idim)=0.75D0-          (x(j,idim)-cc(j,idim)) **2
            wr(j,idim)=0.50D0*(1.5D0-abs(x(j,idim)-cr(j,idim)))**2
          end if
       end do
    end do

    ! Compute parent cell position
    do idim=1,ndim
       do j=1,np
         if(ok(j))then
          icl(j,idim)=int(cl(j,idim))-2*igl(j,idim)
          icc(j,idim)=int(cc(j,idim))-2*igc(j,idim)
          icr(j,idim)=int(cr(j,idim))-2*igr(j,idim)
         else ! ERM: this else may or may not be correct? But I believe it is.
          icl(j,idim)=int(cl(j,idim))
          icc(j,idim)=int(cc(j,idim))
          icr(j,idim)=int(cr(j,idim))
         endif
       end do
    end do

    ! #if NDIM==3
    do j=1,np
      if(ok(j))then
       icell(j,1 )=1+icl(j,1)+2*icl(j,2)+4*icl(j,3)
       icell(j,2 )=1+icc(j,1)+2*icl(j,2)+4*icl(j,3)
       icell(j,3 )=1+icr(j,1)+2*icl(j,2)+4*icl(j,3)
       icell(j,4 )=1+icl(j,1)+2*icc(j,2)+4*icl(j,3)
       icell(j,5 )=1+icc(j,1)+2*icc(j,2)+4*icl(j,3)
       icell(j,6 )=1+icr(j,1)+2*icc(j,2)+4*icl(j,3)
       icell(j,7 )=1+icl(j,1)+2*icr(j,2)+4*icl(j,3)
       icell(j,8 )=1+icc(j,1)+2*icr(j,2)+4*icl(j,3)
       icell(j,9 )=1+icr(j,1)+2*icr(j,2)+4*icl(j,3)
       icell(j,10)=1+icl(j,1)+2*icl(j,2)+4*icc(j,3)
       icell(j,11)=1+icc(j,1)+2*icl(j,2)+4*icc(j,3)
       icell(j,12)=1+icr(j,1)+2*icl(j,2)+4*icc(j,3)
       icell(j,13)=1+icl(j,1)+2*icc(j,2)+4*icc(j,3)
       icell(j,14)=1+icc(j,1)+2*icc(j,2)+4*icc(j,3)
       icell(j,15)=1+icr(j,1)+2*icc(j,2)+4*icc(j,3)
       icell(j,16)=1+icl(j,1)+2*icr(j,2)+4*icc(j,3)
       icell(j,17)=1+icc(j,1)+2*icr(j,2)+4*icc(j,3)
       icell(j,18)=1+icr(j,1)+2*icr(j,2)+4*icc(j,3)
       icell(j,19)=1+icl(j,1)+2*icl(j,2)+4*icr(j,3)
       icell(j,20)=1+icc(j,1)+2*icl(j,2)+4*icr(j,3)
       icell(j,21)=1+icr(j,1)+2*icl(j,2)+4*icr(j,3)
       icell(j,22)=1+icl(j,1)+2*icc(j,2)+4*icr(j,3)
       icell(j,23)=1+icc(j,1)+2*icc(j,2)+4*icr(j,3)
       icell(j,24)=1+icr(j,1)+2*icc(j,2)+4*icr(j,3)
       icell(j,25)=1+icl(j,1)+2*icr(j,2)+4*icr(j,3)
       icell(j,26)=1+icc(j,1)+2*icr(j,2)+4*icr(j,3)
       icell(j,27)=1+icr(j,1)+2*icr(j,2)+4*icr(j,3)
     else
       icell(j,1 )=1+icl(j,1)+3*icl(j,2)+9*icl(j,3)
       icell(j,2 )=1+icc(j,1)+3*icl(j,2)+9*icl(j,3)
       icell(j,3 )=1+icr(j,1)+3*icl(j,2)+9*icl(j,3)
       icell(j,4 )=1+icl(j,1)+3*icc(j,2)+9*icl(j,3)
       icell(j,5 )=1+icc(j,1)+3*icc(j,2)+9*icl(j,3)
       icell(j,6 )=1+icr(j,1)+3*icc(j,2)+9*icl(j,3)
       icell(j,7 )=1+icl(j,1)+3*icr(j,2)+9*icl(j,3)
       icell(j,8 )=1+icc(j,1)+3*icr(j,2)+9*icl(j,3)
       icell(j,9 )=1+icr(j,1)+3*icr(j,2)+9*icl(j,3)
       icell(j,10)=1+icl(j,1)+3*icl(j,2)+9*icc(j,3)
       icell(j,11)=1+icc(j,1)+3*icl(j,2)+9*icc(j,3)
       icell(j,12)=1+icr(j,1)+3*icl(j,2)+9*icc(j,3)
       icell(j,13)=1+icl(j,1)+3*icc(j,2)+9*icc(j,3)
       icell(j,14)=1+icc(j,1)+3*icc(j,2)+9*icc(j,3)
       icell(j,15)=1+icr(j,1)+3*icc(j,2)+9*icc(j,3)
       icell(j,16)=1+icl(j,1)+3*icr(j,2)+9*icc(j,3)
       icell(j,17)=1+icc(j,1)+3*icr(j,2)+9*icc(j,3)
       icell(j,18)=1+icr(j,1)+3*icr(j,2)+9*icc(j,3)
       icell(j,19)=1+icl(j,1)+3*icl(j,2)+9*icr(j,3)
       icell(j,20)=1+icc(j,1)+3*icl(j,2)+9*icr(j,3)
       icell(j,21)=1+icr(j,1)+3*icl(j,2)+9*icr(j,3)
       icell(j,22)=1+icl(j,1)+3*icc(j,2)+9*icr(j,3)
       icell(j,23)=1+icc(j,1)+3*icc(j,2)+9*icr(j,3)
       icell(j,24)=1+icr(j,1)+3*icc(j,2)+9*icr(j,3)
       icell(j,25)=1+icl(j,1)+3*icr(j,2)+9*icr(j,3)
       icell(j,26)=1+icc(j,1)+3*icr(j,2)+9*icr(j,3)
       icell(j,27)=1+icr(j,1)+3*icr(j,2)+9*icr(j,3)
     endif
    end do

    ! Compute parent cell adress
    do ind=1,threetondim
       do j=1,np
         if(ok(j))then
          indp(j,ind)=ncoarse+(icell(j,ind)-1)*ngridmax+igrid(j,ind)
         else ! ERM: for AMR(?) there may be an issue with ind_grid_part(j) being used here.
           indp(j,ind)=nbors_father_cells(ind_grid_part(j),icell(j,ind))
         endif
       end do
    end do

    ! Compute cloud volumes (NDIM==3)
    do j=1,np
       vol(j,1 )=wl(j,1)*wl(j,2)*wl(j,3)
       vol(j,2 )=wc(j,1)*wl(j,2)*wl(j,3)
       vol(j,3 )=wr(j,1)*wl(j,2)*wl(j,3)
       vol(j,4 )=wl(j,1)*wc(j,2)*wl(j,3)
       vol(j,5 )=wc(j,1)*wc(j,2)*wl(j,3)
       vol(j,6 )=wr(j,1)*wc(j,2)*wl(j,3)
       vol(j,7 )=wl(j,1)*wr(j,2)*wl(j,3)
       vol(j,8 )=wc(j,1)*wr(j,2)*wl(j,3)
       vol(j,9 )=wr(j,1)*wr(j,2)*wl(j,3)
       vol(j,10)=wl(j,1)*wl(j,2)*wc(j,3)
       vol(j,11)=wc(j,1)*wl(j,2)*wc(j,3)
       vol(j,12)=wr(j,1)*wl(j,2)*wc(j,3)
       vol(j,13)=wl(j,1)*wc(j,2)*wc(j,3)
       vol(j,14)=wc(j,1)*wc(j,2)*wc(j,3)
       vol(j,15)=wr(j,1)*wc(j,2)*wc(j,3)
       vol(j,16)=wl(j,1)*wr(j,2)*wc(j,3)
       vol(j,17)=wc(j,1)*wr(j,2)*wc(j,3)
       vol(j,18)=wr(j,1)*wr(j,2)*wc(j,3)
       vol(j,19)=wl(j,1)*wl(j,2)*wr(j,3)
       vol(j,20)=wc(j,1)*wl(j,2)*wr(j,3)
       vol(j,21)=wr(j,1)*wl(j,2)*wr(j,3)
       vol(j,22)=wl(j,1)*wc(j,2)*wr(j,3)
       vol(j,23)=wc(j,1)*wc(j,2)*wr(j,3)
       vol(j,24)=wr(j,1)*wc(j,2)*wr(j,3)
       vol(j,25)=wl(j,1)*wr(j,2)*wr(j,3)
       vol(j,26)=wc(j,1)*wr(j,2)*wr(j,3)
       vol(j,27)=wr(j,1)*wr(j,2)*wr(j,3)
    end do

#endif

  ! Gather 3-force
  ! ERM: deleted.
  ! Update 3-velocity
  ! ERM: Block 2. Modifying vp instead of new_vp.
  ! ERM: deleted.

  ! Acceleration forces will be added here. ERM: deleted

  ! Update old dust mass and momentum density variables
  ivar_dust=9
  if(nvar<ivar_dust+ndim)then
     write(*,*)'You need to compile ramses with nvar=',ivar_dust+ndim
     stop
  endif

  do idim=1,ndim ! set vv equal to the velocity.
     do j=1,np
        vv(j,idim)=vp(ind_part(j),idim)
     end do
  end do

  do ind=1,xtondim
     do j=1,np ! deposit the dust mass density.
        if(ok(j))then
           uold(indp(j,ind),ivar_dust)=uold(indp(j,ind),ivar_dust)&
           &+mp(ind_part(j))*vol(j,ind)/vol_loc
        end if
     end do
     do idim=1,ndim
        do j=1,np ! deposit the dust momentum density
           if(ok(j))then
              uold(indp(j,ind),ivar_dust+idim)=uold(indp(j,ind),ivar_dust+idim)&
              &+mp(ind_part(j))*vp(ind_part(j),idim)*vol(j,ind)/vol_loc
           end if
        end do
     end do
  end do

  ! I don't think we actually want to do this until after the Lorentz kick
  call InitStoppingRate(np,dtnew(ilevel),indp,vol,vv,nu_stop,xtondim)
  ! Will have to deposit the grain charges eventually too
  do ind=1,xtondim
     do j=1,np ! deposit the dust mass weighted stopping time.
        if(ok(j))then
           unew(indp(j,ind),ivar_dust)=unew(indp(j,ind),ivar_dust)+&
           &(mp(ind_part(j))*vol(j,ind)/vol_loc)*&!rho^d_ij
           &nu_stop(j)
        end if
     end do
  end do
  !
  ! do ind=1,twotondim
  !    do j=1,np !deposit the dust "mass" density.
  !       if(ok(j))then
  !          unew(indp(j,ind),ivar_dust)=unew(indp(j,ind),ivar_dust)+&
  !          &(mp(ind_part(j))*vol(j,ind)/vol_loc)*&!rho^d_ij
  !          &nu_stop(j)*dtnew(ilevel)/(1.+nu_stop(j)*dtnew(ilevel))
  !       end if
  !    end do
  !    do idim=1,ndim
  !       do j=1,np ! deposit the dust "momentum" density
  !          if(ok(j))then
  !             unew(indp(j,ind),ivar_dust+idim)=unew(indp(j,ind),ivar_dust+idim)+&
  !             &(mp(ind_part(j))*vp(ind_part(j),idim)*vol(j,ind)/vol_loc)*&
  !             &vv(j,idim)*nu_stop(j)*dtnew(ilevel)/(1.+nu_stop(j)*dtnew(ilevel))
  !          end if
  !       end do
  !    end do
  ! end do

  ! Deposit initial dust momentum to new gas momentum.
  ! This variable will collect changes in dust momentum so we can subtract
  ! from the gas at the end of move_fine.
  ! This is unew's gas momentum slot.
  ! Try just adding in the momentum CHANGES.
  ! do ind=1,twotondim
  !    do idim=1,ndim
  !       do j=1,np
  !          if(ok(j))then
  !             unew(indp(j,ind),1+idim)=unew(indp(j,ind),1+idim)+mp(ind_part(j))*vp(ind_part(j),idim)*vol(j,ind)/vol_loc
  !          end if
  !       end do
  !    end do
  ! end do

  ! Gather 3-force
  ! ERM: block 1
  !ff(1:np,1:ndim)=0.0D0
  !if(poisson)then
  !   do ind=1,twotondim
  !      do idim=1,ndim
  !         do j=1,np
  !            ff(j,idim)=ff(j,idim)+f(indp(j,ind),idim)*vol(j,ind)
  !         end do
  !      end do
  !   end do
  !endif


  ! ERM: I don't think we need to set new_vp = vp... In fact, I don't think we
  ! need to do anything with vp at all in this routine.
  !do idim=1,ndim
  !   do j=1,np
  !      new_vp(j,idim)=vp(ind_part(j),idim)
  !   end do
  !end do

  ! No longer Perform the second electric kick

!!$  if(boris.and.constant_t_stop)then
!!$     vv(1:np,1:ndim)=new_vp(1:np,1:ndim)
!!$     call ThirdBorisKick(np,dteff,ctm,ts,bb,uu,vv)
!!$     new_vp(1:np,1:ndim)=vv(1:np,1:ndim)
!!$  endif
!!$
!!$  if(.not.constant_t_stop.and.boris)then
!!$     vv(1:np,1:ndim)=new_vp(1:np,1:ndim)
!!$     ! ERM: Compute the stopping time assuming a constant (and unit) sound speed
!!$     if(.not.constant_t_stop.and.boris)then
!!$       do j=1,np
!!$         tss(j)=rd/(dgr(j)*&
!!$         sqrt(1.0+0.2209*gamma*&
!!$         ((vv(j,1)-uu(j,1))**2+(vv(j,2)-uu(j,2))**2+(vv(j,3)-uu(j,3))**2)))
!!$       end do
!!$     endif
!!$     call ThirdBorisKickWithVarTs(np,dteff,ctm,tss,bb,uu,vv)
!!$     new_vp(1:np,1:ndim)=vv(1:np,1:ndim)
!!$  endif

  !do idim=1,ndim
  !   do j=1,np
  !      vp(ind_part(j),idim)=new_vp(j,idim)
  !   end do
  !end do

end subroutine init_dust
!#########################################################################
!#########################################################################
subroutine InitStoppingRate(nn,dt,indp,vol,v,nu,xtondim)
  ! The following subroutine will alter its last argument, nu
  ! to be a half-step advanced. Because we are operator splitting,
  ! one must use the updated dust and gas velocities.
  ! "Large dust fractions can prevent the propagation of soundwaves"
  ! Above is a paper that we should use to test our code at high mu
  use amr_parameters
  use hydro_parameters
  use hydro_commons, ONLY: uold,unew,smallr,nvar,gamma
  implicit none
  integer ::nn,xtondim ! number of cells
  integer ::ivar_dust ! cell-centered dust variables start.
  real(dp) ::dt ! timestep.
  real(dp)::rd,cs! ERM: Grain size parameter
  real(dp),dimension(1:nvector) ::nu,c
  real(dp),dimension(1:nvector,1:xtondim)::vol
  integer ,dimension(1:nvector,1:xtondim)::indp
  real(dp),dimension(1:nvector),save ::dgr! gas density at grain.
  real(dp),dimension(1:nvector,1:ndim) ::v! grain velocity
  real(dp),dimension(1:nvector,1:xtondim,1:ndim)::big_v
  real(dp),dimension(1:nvector,1:ndim),save ::w! drift at half step.
  integer ::i,j,idim,ind
  ivar_dust=9
  rd = 0.62665706865775*grain_size !constant for epstein drag law. #used to have *sqrt(gamma)
   ! isothermal sound speed... Need to get this right. This works for now,
         ! but only if you have scaled things so that the sound speed is 1.

     if ((constant_t_stop).and.(stopping_rate .lt. 0.0))then ! add a "constant_nu_stop" option so you can turn drag totally off.
       nu(1:nvector)=1./t_stop ! Or better yet, add pre-processor directives to turn drag off.
     else if ((constant_t_stop) .and. (stopping_rate .ge. 0.0))then
       nu(1:nvector)=stopping_rate
     else
     dgr(1:nn) = 0.0D0
     if(boris)then
        do ind=1,xtondim
            do j=1,nn
               dgr(j)=dgr(j)+uold(indp(j,ind),1)*vol(j,ind)
               c(j)=c(j)+vol(j,ind)*&
               & sqrt((uold(indp(j,ind),5) - &
               &0.125D0*(uold(indp(j,ind),1+5)+uold(indp(j,ind),1+nvar))*(uold(indp(j,ind),1+5)+uold(indp(j,ind),1+nvar)) &
               &-0.125D0*(uold(indp(j,ind),2+5)+uold(indp(j,ind),2+nvar))*(uold(indp(j,ind),2+5)+uold(indp(j,ind),2+nvar)) &
               &-0.125D0*(uold(indp(j,ind),3+5)+uold(indp(j,ind),3+nvar))*(uold(indp(j,ind),3+5)+uold(indp(j,ind),3+nvar)) &
               &- 0.5D0*(&
               & uold(indp(j,ind),1+1)*uold(indp(j,ind),1+1)+ &
               & uold(indp(j,ind),2+1)*uold(indp(j,ind),2+1)+ &
               & uold(indp(j,ind),3+1)*uold(indp(j,ind),3+1) &
               &)/max(uold(indp(j,ind),1),smallr))&! Subtract from total energy agnetic and kinetic energies
               &*(gamma-1.0D0)/max(uold(indp(j,ind),1),smallr))
           end do
        end do
     endif

     w(1:nn,1:ndim) = 0.0D0 ! Set to the drift velocity post-Lorentz force
     if(boris .and. supersonic_drag)then
        do ind=1,xtondim
          do idim=1,ndim
            do j=1,nn
               w(j,idim)=w(j,idim)+vol(j,ind)*&
               &(v(j,idim)-uold(indp(j,ind),1+idim)/&
               &max(uold(indp(j,ind),1),smallr))
           end do
         end do
        end do
     endif
     do i=1,nn
       nu(i)=(dgr(i)*c(i)/rd)*sqrt(1.+&
       &0.22089323345553233*&
       &(w(i,1)**2+w(i,2)**2+w(i,3)**2)&
       &/(c(i)*c(i)))
     end do
  endif
end subroutine InitStoppingRate
