subroutine move_fine(ilevel)
  use amr_commons
  use hydro_commons, ONLY: uold,unew,smallr,nvar,gamma
  use pm_commons
  use mpi_mod
  implicit none
  integer::ilevel,xtondim
  !----------------------------------------------------------------------
  ! Update particle position and time-centred velocity at level ilevel.
  ! If particle sits entirely in level ilevel, then use fine grid force
  ! for CIC interpolation. Otherwise, use coarse grid (ilevel-1) force.
  !----------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,next_part,ig,ip,npart1,icpu,ind,iskip,ivar,i
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part
  character(LEN=80)::filename,fileloc
  character(LEN=5)::nchar

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  if(trajectories(1)>0)then
   filename='trajectory.dat'
   call title(myid,nchar)
   fileloc=TRIM(filename)//TRIM(nchar)
   open(25+myid, file = fileloc, status = 'unknown', access = 'append')
  endif

#ifdef TSC
    xtondim=threetondim
#else
    xtondim=twotondim
#endif

  ! Set unew = uold in the active region
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do ivar=2,4
        do i=1,active(ilevel)%ngrid
           unew(active(ilevel)%igrid(i)+iskip,ivar)=&
           &uold(active(ilevel)%igrid(i)+iskip,ivar)
        end do
     end do
  end do
  ! Set unew reception cells to zero
  do icpu=1,ncpu
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do ivar=2,4
           do i=1,reception(icpu,ilevel)%ngrid
              unew(reception(icpu,ilevel)%igrid(i)+iskip,ivar)=0.0D0
           end do
        end do
     end do
  end do

  ! Update particles position and velocity
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
           ! Save next particle  <---- Very important !!!
           next_part=nextp(ipart)
           if(ig==0)then
              ig=1
              ind_grid(ig)=igrid
           end if
           ip=ip+1
           ind_part(ip)=ipart
           ind_grid_part(ip)=ig
           if(ip==nvector)then
              call move1(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel,xtondim)
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
  if(ip>0)call move1(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel,xtondim)

  ! Update MPI boundary conditions for unew for dust mass and momentum densities
  do ivar=2,4 ! Gas momentum indices
     call make_virtual_reverse_dp(unew(1,ivar),ilevel)
  end do

  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do ivar=2,4
        do i=1,active(ilevel)%ngrid
           uold(active(ilevel)%igrid(i)+iskip,ivar)=&
           &unew(active(ilevel)%igrid(i)+iskip,ivar)
        end do
     end do
  end do

  do ivar=2,4 ! Gas momentum indices
     call make_virtual_fine_dp   (uold(1,ivar),ilevel)
  end do

  close(25+myid)
!!!!!!!!!!!!!!!!!!!!!!! NEW !!!!!!!!!!!!!!!!!!!!!!!!!!!
if(simple_boundary)call make_boundary_hydro(ilevel)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

111 format('   Entering move_fine for level ',I2)

end subroutine move_fine
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine move_fine_static(ilevel)
  use amr_commons
  use pm_commons
  use mpi_mod
  implicit none
  integer::ilevel,xtondim
  !----------------------------------------------------------------------
  ! Update particle position and time-centred velocity at level ilevel.
  ! If particle sits entirely in level ilevel, then use fine grid force
  ! for CIC interpolation. Otherwise, use coarse grid (ilevel-1) force.
  !----------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,next_part,ig,ip,npart1,npart2
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel
#ifdef TSC
    xtondim=threetondim
#else
    xtondim=twotondim
#endif

  ! Update particles position and velocity
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
                   & (.not. static_stars .and. is_not_DM(typep(ipart)) )  ) then
                 ! FIXME: there should be a static_sink as well
                 ! FIXME: what about debris?
                 npart2=npart2+1
              endif
           else
              if(.not.static_DM) then
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
                   & (.not. static_stars .and. is_not_DM(typep(ipart)) )  ) then
                 ! FIXME: there should be a static_sink as well
                 ! FIXME: what about debris?
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
              call move1(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel,xtondim)
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
  if(ip>0)call move1(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel,xtondim)

111 format('   Entering move_fine for level ',I2)

end subroutine move_fine_static
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine move1(ind_grid,ind_part,ind_grid_part,ng,np,ilevel,xtondim)
  use amr_commons
  use pm_commons
  use poisson_commons
  use hydro_commons, ONLY: uold,unew,smallr,nvar,gamma
  implicit none
  integer::ng,np,ilevel,xtondim
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !------------------------------------------------------------
  ! This routine computes the force on each particle B2
  ! inverse CIC and computes new positions for all particles.
  ! If particle sits entirely in fine level, then CIC is performed
  ! at level ilevel. Otherwise, it is performed at level ilevel-1.
  ! This routine is called B2 move_fine.
  !------------------------------------------------------------
  logical::error
  integer::i,j,ind,idim,nx_loc,isink,index_part,ivar_dust,iskip,icpu
  real(dp)::dx,dx_loc,scale,vol_loc
  real(dp)::ctm! ERM: recommend 1.15D3

  ! Grid-based arrays
  integer ,dimension(1:nvector),save::father_cell
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  real(dp),dimension(1:nvector),save:: grain_sizes,grain_charges! ERM: fluid density interpolated to grain pos. and stopping times
  ! Particle-based arrays
#ifndef TSC
  logical ,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector,1:ndim),save::x,ff,new_xp,new_vp,dd,dg
  real(dp),dimension(1:nvector,1:ndim),save::vv
  real(dp),dimension(1:nvector,1:ndim),save::bb,uu
  real(dp),dimension(1:nvector,1:twotondim,1:ndim),save::big_vv,big_ww
  real(dp),dimension(1:nvector),save:: nu_stop,mov,dgr,ddgr,ciso
  integer ,dimension(1:nvector,1:ndim),save::ig,id,igg,igd,icg,icd
  real(dp),dimension(1:nvector,1:twotondim),save::vol
  integer ,dimension(1:nvector,1:twotondim),save::igrid,icell,indp,kg
#else
  logical ,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector,1:ndim),save::x,ff,new_xp,new_vp
  real(dp),dimension(1:nvector,1:ndim),save::vv,cl,cr,cc,wl,wr,wc
  real(dp),dimension(1:nvector,1:ndim),save::bb,uu
  real(dp),dimension(1:nvector,1:threetondim,1:ndim),save::big_vv,big_ww
  real(dp),dimension(1:nvector),save:: nu_stop,mov,dgr,ddgr,ciso ! ERM: fluid density interpolated to grain pos. and stopping times
  integer ,dimension(1:nvector,1:ndim),save::igl,igr,igc,icl,icr,icc
  real(dp),dimension(1:nvector,1:threetondim),save::vol
  integer ,dimension(1:nvector,1:threetondim),save::igrid,icell,indp,kg
#endif
  real(dp),dimension(1:3)::skip_loc
  real(dp)::den_dust,den_gas,mom_dust,mom_gas,velocity_com
  ! ERM: w is the cell dust-gas drift, B the mag field.
  ctm = charge_to_mass
  !ts = t_stop!  ERM: Not used if constant_t_stop==.false.
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

  ! Gather neighboring father cells (should be present anytime !)
  do i=1,ng
     father_cell(i)=father(ind_grid(i))
  end do
  call get3cubefather(father_cell,nbors_father_cells,nbors_father_grids,&
       & ng,ilevel)

  ! Rescale particle position at level ilevel
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
  ! Check for illegal moves. Is this different for CIC vs TSC?
  error=.false.
  do idim=1,ndim
     do j=1,np
        if(x(j,idim)<0.5D0.or.x(j,idim)>5.5D0)error=.true.
     end do
  end do
  if(error)then
     write(*,*)'problem in move'
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
  ! Gather center of mass 3-velocity
  ivar_dust=9
  if(nvar<ivar_dust+ndim)then
     write(*,*)'You need to compile ramses with nvar=',ivar_dust+ndim
     stop
  endif
  ! ERM: probably going to delete this stuff.
  !vcom(1:np,1:ndim)=0.0D0 ! Will probably want to break things off and do them separately....
  !if(boris.and.hydro)then
  !   do ind=1,twotondim
  !      do idim=1,ndim
  !         do j=1,np
  !            den_gas=uold(indp(j,ind),1)
  !            mom_gas=uold(indp(j,ind),1+idim)
  !            den_dust=uold(indp(j,ind),ivar_dust)
  !            mom_dust=uold(indp(j,ind),ivar_dust+idim)
  !            velocity_com=(mom_gas*(1.0d0+dtnew(ilevel)/ts)+mom_dust*dtnew(ilevel)/ts)/(den_gas*(1.0d0+dtnew(ilevel)/ts)+den_dust*dtnew(ilevel)/ts)
  !            vcom(j,idim)=vcom(j,idim)+velocity_com*vol(j,ind)
!              write(*,*)idim,vcom(j,idim),den_gas,mom_gas,den_dust,mom_dust
  !         end do
  !      end do
  !   end do
  !endif

  ! Various fields interpolated to particle positions
  ! Gather 3-velocity and 3-magnetic field
  uu(1:np,1:ndim)=0.0D0
  bb(1:np,1:ndim)=0.0D0
  if(boris.and.hydro)then
     do ind=1,xtondim
        do idim=1,ndim
           do j=1,np
              uu(j,idim)=uu(j,idim)+uold(indp(j,ind),idim+1)/max(uold(indp(j,ind),1),smallr)*vol(j,ind)
              bb(j,idim)=bb(j,idim)+0.5D0*(uold(indp(j,ind),idim+5)+uold(indp(j,ind),idim+nvar))*vol(j,ind)
           end do
        end do
     end do
  endif

  dgr(1:np) = 0.0D0 ! Gas density. Only used for trajectory file.
  ddgr(1:np) = 0.0D0 ! Dust density.
  ciso(1:np) = 0.0D0 ! Isothermal sound speed (squared).
  if(boris)then
     do ind=1,xtondim
         do j=1,np
            dgr(j)=dgr(j)+uold(indp(j,ind),1)*vol(j,ind)
            ddgr(j)=ddgr(j)+uold(indp(j,ind),ivar_dust)*vol(j,ind)
            ciso(j)=ciso(j)+vol(j,ind)*&
            & sqrt((uold(indp(j,ind),5) - &
            &0.125D0*(uold(indp(j,ind),1+5)+uold(indp(j,ind),1+nvar))*(uold(indp(j,ind),1+5)+uold(indp(j,ind),1+nvar)) &
            &-0.125D0*(uold(indp(j,ind),2+5)+uold(indp(j,ind),2+nvar))*(uold(indp(j,ind),2+5)+uold(indp(j,ind),2+nvar)) &
            &-0.125D0*(uold(indp(j,ind),3+5)+uold(indp(j,ind),3+nvar))*(uold(indp(j,ind),3+5)+uold(indp(j,ind),3+nvar)) &
            &- 0.5D0*(&
            & uold(indp(j,ind),1+1)*uold(indp(j,ind),1+1)+ &
            & uold(indp(j,ind),2+1)*uold(indp(j,ind),2+1)+ &
            & uold(indp(j,ind),3+1)*uold(indp(j,ind),3+1) &
            &)/max(uold(indp(j,ind),1),smallr))&! Subtract from total energy agnetic and kinetic energies
            &*(gamma-1.0D0)/max(uold(indp(j,ind),1),smallr))! P/rho = ciso**2, to be interpolated.
        end do
     end do
  endif

  ! if(trajectories(1)>0)then!Various fields interpolated to particle positions
  !    do index_part=trajectories(1),trajectories(2)
  !       do j=1,np
  !          if(idp(ind_part(j)).EQ.index_part)then
  !             write(25+myid,*)t-dtnew(ilevel),idp(ind_part(j)),dgr(j),ddgr(j), & ! Old time
  !                  & xp(ind_part(j),1),xp(ind_part(j),2),xp(ind_part(j),3),& ! Old particle position
  !                  & vp(ind_part(j),1),vp(ind_part(j),2),vp(ind_part(j),3),& ! Old particle velocity
  !                  &  uu(j,1),uu(j,2),uu(j,3),& ! Old fluid velocity
  !                  &  bb(j,1),bb(j,2),bb(j,3)! Old magnetic field.
  !                  ! & new_vp(j,1),new_vp(j,2),new_vp(j,3) ! NEW particle velocity (for comparison)
  !          endif
  !       end do
  !    end do
  ! endif

 if (trajectories(1)>0)then
        !Various fields interpolated to particle positions
         do j=1,np
             if( ANY(trajectories .eq. idp(ind_part(j))) )then
                write(25+myid,*)t-dtnew(ilevel),idp(ind_part(j)),dgr(j),ddgr(j), & ! Old time
                     & xp(ind_part(j),1),xp(ind_part(j),2),xp(ind_part(j),3),& ! Old particle position
                     & vp(ind_part(j),1),vp(ind_part(j),2),vp(ind_part(j),3),& ! Old particle velocity
                     &  uu(j,1),uu(j,2),uu(j,3),& ! Old fluid velocity
                     &  bb(j,1),bb(j,2),bb(j,3)! Old magnetic field.
                     ! & new_vp(j,1),new_vp(j,2),new_vp(j,3) ! NEW particle velocity (for comparison)
             endif
       end do
 endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! ERM: Block here is only used for computing variable stopping times.

  if(boris)then ! ERM: may not be needed.
    do j=1,np
        mov(j) = mp(ind_part(j))/vol_loc
    end do
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Gather 3-velocity
  ff(1:np,1:ndim)=0.0D0
  if(tracer.and.hydro)then
     do ind=1,xtondim
        do idim=1,ndim
           do j=1,np
              ff(j,idim)=ff(j,idim)+uold(indp(j,ind),idim+1)/max(uold(indp(j,ind),1),smallr)*vol(j,ind)
           end do
        end do
     end do
  endif
  ! Gather 3-force
  if(poisson)then
     do ind=1,xtondim
        do idim=1,ndim
           do j=1,np
              ff(j,idim)=ff(j,idim)+f(indp(j,ind),idim)*vol(j,ind)
           end do
        end do
#ifdef OUTPUT_PARTICLE_POTENTIAL
        do j=1,np
           ptcl_phi(ind_part(j)) = phi(indp(j,ind))
        end do
#endif
     end do
  endif

  ! Update velocity
  do idim=1,ndim
     if(static.or.tracer)then
        do j=1,np
           new_vp(j,idim)=ff(j,idim)
        end do
     else
        do j=1,np
           new_vp(j,idim)=vp(ind_part(j),idim)+ff(j,idim)*0.5D0*dtnew(ilevel)
        end do
     endif
  end do

  if(boris.and.hydro)then
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! GRAIN SIZES AND CHARGES
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (ddex.ne.0.0)then
      do j=1,np ! construct charges and grain sizes
        grain_sizes(j)=grain_size*1.0d1**(ddex*idp(ind_part(j))*0.5d0**(3.*levelmin)/ndust)
        grain_charges(j)=charge_to_mass*(grain_sizes(j)/grain_size)**charge_slope
      end do
    else
      grain_sizes(1:np)=grain_size
      grain_charges(1:np)=charge_to_mass
    endif
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! STOPPING RATE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! We compute this before the EM kick for the sole reason that we had to do
  ! that in init_dust_fine. For second order accuracy, things will be more
  ! complicated.
  call StoppingRate(np,dtnew(ilevel),indp,vol,vv,nu_stop,ciso,dgr,grain_sizes,xtondim)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! LORENTZ KICK
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! ERM: Determining the evolved versions of the drift, then interpolating those
  ! For now, do the 8N calculations, but in the future, would be good to loop
  ! over CELLS rather than particles (esp. when doing TSC).
  ! Set unew's dust momentum slots to be the gas velocity.

    !call ResetUnewToFluidVel(ilevel)
    !call reset_unew(ilevel)
    big_vv(1:np,1:xtondim,1:ndim)=0.0D0 ! contains actual sub-cloud velocities.
    big_ww(1:np,1:xtondim,1:ndim)=0.0D0 ! Contains net mean drift velocity
    ! might want a "big_ww"? I think that's how I'll approach it.
    ! We want to evolve each of the subclouds. Knowing the new w will
    ! allow us to compute the evolution of the sub-clouds with the drag too.
    vv(1:np,1:ndim)=new_vp(1:np,1:ndim)
    call EMKick(np,dtnew(ilevel),indp,grain_charges,ok,vol,mov,vv,big_vv,big_ww,xtondim)
    ! big_vv now contains changes to sub-cloud velocities. vv is still the old
    ! velocity. As well, unew's dust slot contains u**n+du**EM
    !write(*,*)'big_vv=',big_vv(1,1,1),big_vv(1,1,2),big_vv(1,1,3)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! DRAG KICK
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call DragKick(np,dtnew(ilevel),indp,ok,vol,nu_stop,big_vv,big_ww,vv,xtondim)
    ! Was just "DragKick"
    !write(*,*)'big_vv+=',big_vv(1,1,1),big_vv(1,1,2),big_vv(1,1,3)
    ! DragKick will modify big_ww as well as big_vv, but not vv.
    ! Now kick the dust given these quantities.
    vv(1:np,1:ndim)=0.0D0
    do ind=1,xtondim
      do idim=1,ndim
        do j=1,np
          vv(j,idim)=vv(j,idim)+vol(j,ind)*big_vv(j,ind,idim)
        end do
      end do
    end do

    !call DragKick(np,dtnew(ilevel),indp,ok,vol,mov,nu_stop,big_vv,vv)
    new_vp(1:np,1:ndim)=vv(1:np,1:ndim)
    ! big_vv is not actually modified in this process:
    ! Rather, we go straight to interpolating onto vv.
  endif

  ! For sink cloud particle only
  if(sink)then
     ! Overwrite cloud particle velocity with sink velocity
     do idim=1,ndim
        do j=1,np
           if( is_cloud(typep(ind_part(j))) ) then
              isink=-idp(ind_part(j))
              new_vp(j,idim)=vsnew(isink,idim,ilevel)
           end if
        end do
     end do
  end if

  !
  ! ! Output data to trajectory file
  ! ! May have to think more carefully about when and where this is placed
  ! if((boris.or.tracer).and.constant_t_stop)then
  !   !Various fields interpolated to particle positions
  !   !Gather 3-velocity and 3-magnetic field
  !   ! uu(1:10,1:ndim)=0.0D0
  !   ! bb(1:10,1:ndim)=0.0D0
  !   ! if(boris.and.hydro)then
  !   !    do ind=1,twotondim
  !   !       do idim=1,ndim
  !   !         do index_part=1,10
  !   !          do j=1,np
  !   !            if(idp(ind_part(j)).EQ.index_part)then
  !   !             uu(index_part,idim)=uu(index_part,idim)+uold(indp(j,ind),idim+1)/max(uold(indp(j,ind),1),smallr)*vol(j,ind)
  !   !             bb(index_part,idim)=bb(index_part,idim)+0.5D0*(uold(indp(j,ind),idim+5)+uold(indp(j,ind),idim+nvar))*vol(j,ind)
  !   !            endif
  !   !          end do
  !   !         end do
  !   !       end do
  !   !    end do
  !   ! endif
  !    do index_part=1,10
  !       do j=1,np
  !          if(idp(ind_part(j)).EQ.index_part)then
  !             write(25+myid,*)t-dtnew(ilevel),idp(ind_part(j)),& ! Old time
  !                  & xp(ind_part(j),1),xp(ind_part(j),2),xp(ind_part(j),3),& ! Old particle position
  !                  & vp(ind_part(j),1),vp(ind_part(j),2),vp(ind_part(j),3) ! Old particle velocity
  !                   ! &  uu(index_part,1),uu(index_part,2),uu(index_part,3),& ! Old fluid velocity
  !                   ! &  bb(index_part,1),bb(index_part,2),bb(index_part,3) ! Old magnetic field.
  !          endif
  !       end do
  !    end do
  ! endif

  ! Update position BEFORE setting new velocity using trapezoidal rule.
  do idim=1,ndim
     if(static)then
        do j=1,np
           new_xp(j,idim)=xp(ind_part(j),idim)
        end do
     else
        do j=1,np
           new_xp(j,idim)=xp(ind_part(j),idim)+0.5*(new_vp(j,idim)+vp(ind_part(j),idim))*dtnew(ilevel)
        end do
     endif
  end do
  do idim=1,ndim
     do j=1,np
        xp(ind_part(j),idim)=new_xp(j,idim) ! This is where you'd do the monte-carlo bit?
     end do
  end do

  ! Deposit minus final dust momentum to new gas momentum
  do ind=1,xtondim
     do idim=1,ndim
        do j=1,np
           if(ok(j))then
              unew(indp(j,ind),1+idim)=unew(indp(j,ind),1+idim)&
              &-mp(ind_part(j))*(new_vp(j,idim)-vp(ind_part(j),idim))*vol(j,ind)/vol_loc
           end if
        end do
     end do
  end do

  ! Deposit minus final dust energy onto gas energy.
  do ind=1,xtondim
        do j=1,np
           if(ok(j))then
              unew(indp(j,ind),5)=unew(indp(j,ind),5)&
              &-0.5*mp(ind_part(j))*&
              &(new_vp(j,idim)*new_vp(j,idim)-vp(ind_part(j),idim)*vp(ind_part(j),idim))&
              &*vol(j,ind)/vol_loc
           end if
        end do
  end do

  ! Store new velocity
  do idim=1,ndim
     do j=1,np
        vp(ind_part(j),idim)=new_vp(j,idim)
     end do
  end do

end subroutine move1
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine EMKick(nn,dt,indp,ctms,ok,vol,mov,v,big_v,big_w,xtondim)
  ! This subroutine will compute changes to sub-cloud velocity in big_v,
  ! as well as set unew's dust momentum slot to being u+du**EM.
  use amr_parameters
  use hydro_parameters
  use hydro_commons, ONLY: uold,unew,smallr,nvar,gamma
  implicit none
  integer ::ivar_dust,xtondim ! cell-centered dust variables start.
  integer ::nn ! number of cells
  real(dp) ::dt,ctmav! timestep, average charge to mass ratio
  logical ,dimension(1:nvector)::ok
  real(dp),dimension(1:nvector)::mov,ctms !mass over volume, charge to mass ratios
  real(dp),dimension(1:nvector,1:xtondim)::vol
  integer ,dimension(1:nvector,1:xtondim)::indp
  real(dp),dimension(1:nvector,1:xtondim,1:ndim)::big_v,big_w
  real(dp),dimension(1:nvector,1:ndim) ::v! grain velocity
  real(dp) ::den_dust,den_gas,mu
  real(dp),dimension(1:3) ::vtemp,w,B
  integer ::i,j,ind,idim! Just an -index

  ivar_dust=9

  do ind=1,xtondim
     do i=1,nn
        den_gas=uold(indp(i,ind),1)
        den_dust=uold(indp(i,ind),ivar_dust)
        mu=den_dust/max(den_gas,smallr)
        ctmav=unew(indp(i,ind),ivar_dust+1)
        do idim=1,3
          B(idim)=0.5D0*(uold(indp(i,ind),idim+5)+uold(indp(i,ind),idim+nvar))
          w(idim)=uold(indp(i,ind),idim+ivar_dust)/max(uold(indp(i,ind),ivar_dust),smallr)&
          &-uold(indp(i,ind),idim+1)/max(uold(indp(i,ind),1),smallr)
        end do

        big_w(i,ind,1)=-1.*& !velocity changes to drift
        &(ctmav*dt*(2.*B(2)**2*ctmav*dt*(1.+mu)**2*w(1)+&
        &B(2)*(-2.*B(1)*ctmav*dt*(1.+mu)**2*w(2) + 4.*(1.+mu)*w(3)) +&
        & B(3)*(2.*B(3)*ctmav*dt*(1.+mu)**2*w(1)-4.*(1.+mu)*w(2) - 2*B(1)*ctmav*dt*(1.+mu)**2*w(3))))/&
        &((4.+(B(1)**2+B(2)**2+B(3)**2)*ctmav**2*dt**2*(1.+mu)**2))

        big_w(i,ind,2)=-1.*& !velocity changes to drift
        &(ctmav*dt*(2.*B(3)**2*ctmav*dt*(1.+mu)**2*w(2)+&
        &B(3)*(-2.*B(2)*ctmav*dt*(1.+mu)**2*w(3) + 4.*(1.+mu)*w(1)) +&
        & B(1)*(2.*B(1)*ctmav*dt*(1.+mu)**2*w(2)-4.*(1.+mu)*w(3) - 2*B(2)*ctmav*dt*(1.+mu)**2*w(1))))/&
        &((4.+(B(1)**2+B(2)**2+B(3)**2)*ctmav**2*dt**2*(1.+mu)**2))

        big_w(i,ind,3)=-1.*& !velocity changes to drift
        &(ctmav*dt*(2.*B(1)**2*ctmav*dt*(1.+mu)**2*w(3)+&
        &B(1)*(-2.*B(3)*ctmav*dt*(1.+mu)**2*w(1) + 4.*(1.+mu)*w(2)) +&
        & B(2)*(2.*B(2)*ctmav*dt*(1.+mu)**2*w(3)-4.*(1.+mu)*w(1) - 2*B(3)*ctmav*dt*(1.+mu)**2*w(2))))/&
        &((4.+(B(1)**2+B(2)**2+B(3)**2)*ctmav**2*dt**2*(1.+mu)**2))

        do idim=1,ndim
          vtemp(idim) = v(i,idim)-uold(indp(i,ind),1+idim)/max(uold(indp(i,ind),1),smallr)&
          &+0.5*mu*big_w(i,ind,idim)/(1.+mu) ! v^n - u^(n+1)
          big_w(i,ind,idim)=w(idim)+big_w(i,ind,idim) !w^(n+1)
        end do

        big_v(i,ind,1)=v(i,1)+& ! subcloud velocity update
        &(ctms(i)*dt*(-2.*B(2)**2*ctms(i)*dt*vtemp(1) + B(2)*(2.*B(1)*ctms(i)*dt*vtemp(2) - 4.*vtemp(3)) +&
        & B(3)*(-2.*B(3)*ctms(i)*dt*vtemp(1) + 4.*vtemp(2) + 2.*B(1)*ctms(i)*dt*vtemp(3))))/&
        &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctms(i)**2*dt**2)

        big_v(i,ind,2)=v(i,2)+&
        &(ctms(i)*dt*(-2.*B(3)**2*ctms(i)*dt*vtemp(2) + B(3)*(2.*B(2)*ctms(i)*dt*vtemp(3) - 4.*vtemp(1)) +&
        & B(1)*(-2.*B(1)*ctms(i)*dt*vtemp(2) + 4.*vtemp(3) + 2.*B(2)*ctms(i)*dt*vtemp(1))))/&
        &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctms(i)**2*dt**2)

        big_v(i,ind,3)=v(i,3)+&
        &(ctms(i)*dt*(-2.*B(1)**2*ctms(i)*dt*vtemp(3) + B(1)*(2.*B(3)*ctms(i)*dt*vtemp(1) - 4.*vtemp(2)) +&
        & B(2)*(-2.*B(2)*ctms(i)*dt*vtemp(3) + 4.*vtemp(1) + 2.*B(3)*ctms(i)*dt*vtemp(2))))/&
        &(4. + (B(1)**2 + B(2)**2 + B(3)**2)*ctms(i)**2*dt**2)
     end do
  end do
end subroutine EMKick
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine StoppingRate(nn,dt,indp,vol,v,nu,c,dgr,gs,xtondim)
  ! The following subroutine will alter its last argument, nu
  ! to be a half-step advanced. Because we are operator splitting,
  ! one must use the updated dust and gas velocities.
  ! "Large dust fractions can prevent the propagation of soundwaves"
  ! Above is a paper that we should use to test our code at high mu
  use amr_parameters
  use hydro_parameters
  use hydro_commons, ONLY: uold,unew,smallr,nvar,gamma
  implicit none
  integer ::nn ! number of cells
  integer ::ivar_dust,xtondim  ! cell-centered dust variables start.
  real(dp) ::dt ! timestep.
  real(dp)::rd,cs! ERM: Grain size parameter
  real(dp),dimension(1:nvector) ::nu,c
  real(dp),dimension(1:nvector,1:xtondim)::vol
  integer ,dimension(1:nvector,1:xtondim)::indp
  real(dp),dimension(1:nvector) ::dgr,gs! gas density at grain, grain size array
  real(dp),dimension(1:nvector,1:ndim) ::v! grain velocity
  real(dp),dimension(1:nvector,1:xtondim,1:ndim)::big_v
  real(dp),dimension(1:nvector,1:ndim),save ::wh! drift at half step.
  integer ::i,j,idim,ind
  ivar_dust=9
  rd = 0.62665706865775 !*sqrt(gamma) constant for epstein drag law.

  if ((constant_t_stop).and.(stopping_rate .lt. 0.0))then ! add a "constant_nu_stop" option so you can turn drag totally off.
    nu(1:nvector)=(1./t_stop)*grain_size/gs(1:nvector) ! Or better yet, add pre-processor directives to turn drag off.
  else if ((constant_t_stop) .and. (stopping_rate .ge. 0.0))then
    nu(1:nvector)=stopping_rate*grain_size/gs(1:nvector)
  else
     !dgr(1:nn) = 0.0D0 ! I don't have to do this twice... It's done previously...
     !if(boris)then
    !    do ind=1,xtondim
  !          do j=1,nn
!               dgr(j)=dgr(j)+uold(indp(j,ind),1)*vol(j,ind)
    !       end do
    !    end do
    ! endif

     wh(1:nn,1:ndim) = 0.0D0 ! Set to the drift velocity post-Lorentz force
     if(boris .and. supersonic_drag)then
        do ind=1,xtondim
          do idim=1,ndim
            do j=1,nn
               wh(j,idim)=wh(j,idim)+vol(j,ind)*&
               &(v(j,idim)-uold(indp(j,ind),1+idim)/&
               &max(uold(indp(j,ind),1),smallr))
           end do
         end do
        end do
     endif
     do i=1,nn
       nu(i)=(dgr(i)*c(i)/(rd*gs(i)))*sqrt(1.+&
       &0.22089323345553233*&
       &(wh(i,1)**2+wh(i,2)**2+wh(i,3)**2)&
       &/(c(i)*c(i)))
     end do
  endif
end subroutine StoppingRate
!#########################################################################
!#########################################################################
! !#########################################################################
! !#########################################################################
! subroutine OldDragKick(nn,dt,indp,ok,vol,nu,big_v,big_w,v) ! mp is actually mov
!   use amr_parameters
!   use hydro_parameters
!   use hydro_commons, ONLY: uold,unew,smallr,nvar,gamma
!   implicit none
!   integer::nn
!   integer ::ivar_dust ! cell-centered dust variables start.  integer ::nn ! number of cells
!   real(dp) ::dt ! timestep
!   real(dp) ::vol_loc ! cloud volume
!   real(dp),dimension(1:nvector) ::nu,mp
!   logical ,dimension(1:nvector)::ok
!   real(dp),dimension(1:nvector,1:twotondim)::vol
!   integer ,dimension(1:nvector,1:twotondim)::indp
!   real(dp),dimension(1:nvector,1:ndim) ::v ! grain velocity
!   real(dp),dimension(1:nvector,1:twotondim,1:ndim) ::big_v,big_w
!   real(dp) ::den_dust,den_gas,mu,nuj,vo
!   integer ::i,j,ind,idim! Just an index
!   ivar_dust=9
!
!   do ind=1,twotondim
!      do i=1,nn
!         den_gas=uold(indp(i,ind),1)
!         den_dust=uold(indp(i,ind),ivar_dust)
!         mu=den_dust/max(den_gas,smallr)
!         nuj=(1.+mu)*unew(indp(i,ind),ivar_dust)/max(uold(indp(i,ind),ivar_dust),smallr)
!         do idim=1,ndim
!           ! w = &
!           ! &uold(indp(i,ind),ivar_dust+idim)/max(uold(indp(i,ind),ivar_dust),smallr)&
!           ! &-uold(indp(i,ind),1+idim)/max(uold(indp(i,ind),1),smallr)
!
!           big_w(i,ind,idim)=big_w(i,ind,idim)&
!           &/(1.+nuj*dt+0.5*nuj*nuj*dt*dt)
!
!           vo = -big_v(i,ind,idim) +(uold(indp(i,ind),1+idim)&
!           &+uold(indp(i,ind),ivar_dust+idim))/&
!           &max(uold(indp(i,ind),1)+uold(indp(i,ind),ivar_dust),smallr)
!
!           big_v(i,ind,idim)=big_v(i,ind,idim)+&
!           &((dt*nu(i)+0.5*dt*dt*nu(i)*nu(i))*vo&
!           &-dt*nu(i)*mu*(1.+0.5*dt*(nu(i)-nuj))*big_w(i,ind,idim)/(1.+mu))&
!           &/(1.+dt*nu(i)+0.5*dt*dt*nu(i)*nu(i))
!         end do
!         ! big_w corresponds directly to a change in the gas velocity.
!      end do
!   end do
! end subroutine OldDragKick
!
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine DragKick(nn,dt,indp,ok,vol,nu,big_v,big_w,v,xtondim) ! mp is actually mov
  use amr_parameters
  use hydro_parameters
  use hydro_commons, ONLY: uold,unew,smallr,nvar,gamma
  implicit none
  integer::nn
  integer ::ivar_dust, xtondim ! cell-centered dust variables start.  integer ::nn ! number of cells
  real(dp) ::dt ! timestep
  real(dp) ::vol_loc ! cloud volume
  real(dp),dimension(1:nvector) ::nu,mp
  logical ,dimension(1:nvector)::ok
  real(dp),dimension(1:nvector,1:xtondim)::vol
  integer ,dimension(1:nvector,1:xtondim)::indp
  real(dp),dimension(1:nvector,1:ndim) ::v ! grain velocity
  real(dp),dimension(1:nvector,1:xtondim,1:ndim) ::big_v,big_w
  real(dp) ::den_dust,den_gas,mu,nuj,up
  integer ::i,j,ind,idim! Just an index
  ivar_dust=9

  if (second_order)then
    do ind=1,xtondim
       do i=1,nn
          den_gas=uold(indp(i,ind),1)
          den_dust=uold(indp(i,ind),ivar_dust)
          mu=den_dust/max(den_gas,smallr)
          nuj=(1.+mu)*unew(indp(i,ind),ivar_dust)/max(uold(indp(i,ind),ivar_dust),smallr)
          do idim=1,ndim
            ! w = &
            ! &uold(indp(i,ind),ivar_dust+idim)/max(uold(indp(i,ind),ivar_dust),smallr)&
            ! &-uold(indp(i,ind),1+idim)/max(uold(indp(i,ind),1),smallr)

            big_w(i,ind,idim)=big_w(i,ind,idim)&
            &/(1.+nuj*dt+0.5*nuj*nuj*dt*dt)

            up = -mu*big_w(i,ind,idim)*(1.+0.5*dt*nuj)/(1.+mu)+(uold(indp(i,ind),1+idim)&
            &+uold(indp(i,ind),ivar_dust+idim))/&
            &max(uold(indp(i,ind),1)+uold(indp(i,ind),ivar_dust),smallr)

            big_v(i,ind,idim)=(big_v(i,ind,idim)+dt*nu(i)*(1.+0.5*dt*nu(i))*up)&
            &/(1.+dt*nu(i)*(1.+0.5*dt*nu(i)))
            ! up = -mu*big_w(i,ind,idim)*(1.+0.5*dt*nuj)/(1.+mu)+(uold(indp(i,ind),1+idim)&
            ! &+uold(indp(i,ind),ivar_dust+idim))/&
            ! &max(uold(indp(i,ind),1)+uold(indp(i,ind),ivar_dust),smallr)
            !
            ! big_v(i,ind,idim)=(big_v(i,ind,idim)+dt*nu(i)*(1.+0.5*dt*nu(i))*up)&
            ! &/(1.+dt*nu(i)*(1.+0.5*dt*nu(i)))
          end do
          ! big_w corresponds directly to a change in the gas velocity.
       end do
    end do
   else
     do ind=1,xtondim
        do i=1,nn
           den_gas=uold(indp(i,ind),1)
           den_dust=uold(indp(i,ind),ivar_dust)
           mu=den_dust/max(den_gas,smallr)
           nuj=(1.+mu)*unew(indp(i,ind),ivar_dust)/max(uold(indp(i,ind),ivar_dust),smallr)
           do idim=1,ndim
             ! w = &
             ! &uold(indp(i,ind),ivar_dust+idim)/max(uold(indp(i,ind),ivar_dust),smallr)&
             ! &-uold(indp(i,ind),1+idim)/max(uold(indp(i,ind),1),smallr)

             big_w(i,ind,idim)=big_w(i,ind,idim)&
             &/(1.+nuj*dt)

             up = -mu*big_w(i,ind,idim)/(1.+mu)+(uold(indp(i,ind),1+idim)&
             &+uold(indp(i,ind),ivar_dust+idim))/&
             &max(uold(indp(i,ind),1)+uold(indp(i,ind),ivar_dust),smallr)

             big_v(i,ind,idim)=(big_v(i,ind,idim)+dt*nu(i)*up)/(1.+dt*nu(i))
           end do
           ! big_w corresponds directly to a change in the gas velocity.
        end do
     end do
   endif
end subroutine DragKick


subroutine DragKickAlt(nn,dt,indp,ok,vol,nu,big_v,big_w,v,xtondim) ! mp is actually mov
  use amr_parameters
  use hydro_parameters
  use hydro_commons, ONLY: uold,unew,smallr,nvar,gamma
  implicit none
  integer::nn
  integer ::ivar_dust, xtondim ! cell-centered dust variables start.  integer ::nn ! number of cells
  real(dp) ::dt ! timestep
  real(dp) ::vol_loc ! cloud volume
  real(dp),dimension(1:nvector) ::nu,mp
  logical ,dimension(1:nvector)::ok
  real(dp),dimension(1:nvector,1:xtondim)::vol
  integer ,dimension(1:nvector,1:xtondim)::indp
  real(dp),dimension(1:nvector,1:ndim) ::v ! grain velocity
  real(dp),dimension(1:nvector,1:xtondim,1:ndim) ::big_v,big_w
  real(dp) ::den_dust,den_gas,mu,nuj,up
  integer ::i,j,ind,idim! Just an index
  ivar_dust=9

  if (second_order)then
    do ind=1,xtondim
       do i=1,nn
          den_gas=uold(indp(i,ind),1)
          den_dust=uold(indp(i,ind),ivar_dust)
          mu=den_dust/max(den_gas,smallr)
          nuj=(1.+mu)*unew(indp(i,ind),ivar_dust)/max(uold(indp(i,ind),ivar_dust),smallr)
          do idim=1,ndim
            ! w = &
            ! &uold(indp(i,ind),ivar_dust+idim)/max(uold(indp(i,ind),ivar_dust),smallr)&
            ! &-uold(indp(i,ind),1+idim)/max(uold(indp(i,ind),1),smallr)

            big_w(i,ind,idim)=big_w(i,ind,idim)&
            &/(1.+nuj*dt+0.5*nuj*nuj*dt*dt)

            up = -mu*big_w(i,ind,idim)/(1.+mu)+(uold(indp(i,ind),1+idim)&
            &+uold(indp(i,ind),ivar_dust+idim))/&
            &max(uold(indp(i,ind),1)+uold(indp(i,ind),ivar_dust),smallr)

            big_v(i,ind,idim)= up + (big_v(i,ind,idim)-up)&
            &/(1.+dt*nu(i)*(1.+0.5*dt*nu(i)))
            !big_v(i,ind,idim)=(big_v(i,ind,idim)+dt*nu(i)*(1.+0.5*dt*nu(i))*up)&
            !&/(1.+dt*nu(i)*(1.+0.5*dt*nu(i))) ! suspected error begins here?
          end do
          ! big_w corresponds directly to a change in the gas velocity.
       end do
    end do
   else ! Change this to be the original algorithm
     do ind=1,xtondim
        do i=1,nn
           den_gas=uold(indp(i,ind),1)
           den_dust=uold(indp(i,ind),ivar_dust)
           mu=den_dust/max(den_gas,smallr)
           nuj=(1.+mu)*unew(indp(i,ind),ivar_dust)/max(uold(indp(i,ind),ivar_dust),smallr)
           do idim=1,ndim
             ! w = &
             ! &uold(indp(i,ind),ivar_dust+idim)/max(uold(indp(i,ind),ivar_dust),smallr)&
             ! &-uold(indp(i,ind),1+idim)/max(uold(indp(i,ind),1),smallr)

             big_w(i,ind,idim)=big_w(i,ind,idim)&
             &/(1.+nuj*dt)

             up = -mu*big_w(i,ind,idim)/(1.+mu)+(uold(indp(i,ind),1+idim)&
             &+uold(indp(i,ind),ivar_dust+idim))/&
             &max(uold(indp(i,ind),1)+uold(indp(i,ind),ivar_dust),smallr)

             big_v(i,ind,idim)=(big_v(i,ind,idim)+dt*nu(i)*up)/(1.+dt*nu(i))
           end do
           ! big_w corresponds directly to a change in the gas velocity.
        end do
     end do
   endif
end subroutine DragKickAlt
!# Later on, can implement the better algorithm that does everything
! self-consistently.
!#########################################################################
!#########################################################################

! subroutine StoppingRateMidpt(nn,twodt,indp,ok,vol,mov,v,big_v,nu)
!   ! The following subroutine will alter its last argument, nu
!   ! to be a half-step advanced. Because we are operator splitting,
!   ! one must use the updated dust and gas velocities.
!   ! "Large dust fractions can prevent the propagation of soundwaves"
!   ! Above is a paper that we should use to test our code at high mu
!   use amr_parameters
!   use hydro_parameters
!   use hydro_commons, ONLY: uold,unew,smallr,nvar,gamma
!   implicit none
!   integer ::nn ! number of cells
!   integer ::ivar_dust ! cell-centered dust variables start.
!   real(dp) ::dt,twodt ! half-timestep, full timestep.
!   real(dp)::rd,cs! ERM: Grain size parameter
!   real(dp),dimension(1:nvector) ::nu
!   real(dp),dimension(1:nvector) ::mov
!   real(dp),dimension(1:nvector,1:twotondim)::vol
!   integer ,dimension(1:nvector,1:twotondim)::indp
!   logical ,dimension(1:nvector)::ok
!   real(dp),dimension(1:nvector),save ::dgr! gas density at grain.
!   real(dp),dimension(1:nvector,1:ndim) ::v! grain velocity
!   real(dp),dimension(1:nvector,1:twotondim,1:ndim)::big_v
!   real(dp),dimension(1:nvector,1:ndim),save ::wh! drift at half step.
!   integer ::i,j,idim,ind
!   ivar_dust=9
!   dt=0.5*twodt
!   rd = sqrt(gamma)*0.62665706865775*grain_size !constant for epstein drag law.
!   cs=1.0 ! isothermal sound speed... Need to get this right. This works for now,
!          ! but only if you have scaled things so that the sound speed is 1.
!
!   if (constant_t_stop)then
!     nu(1:nvector)=1./t_stop
!   else
!      dgr(1:nn) = 0.0D0
!      if(boris)then
!         do ind=1,twotondim
!             do j=1,nn
!                dgr(j)=dgr(j)+uold(indp(j,ind),1)*vol(j,ind)
!            end do
!         end do
!      endif
!
!      wh(1:nn,1:ndim) = 0.0D0 ! Set to the drift velocity post-Lorentz force
!      if(boris)then
!         do ind=1,twotondim
!           do idim=1,ndim
!             do j=1,nn
!                wh(j,idim)=wh(j,idim)+&
!                &(v(j,idim)+big_v(j,ind,idim)-unew(indp(j,ind),ivar_dust+idim))&
!                &*vol(j,ind)
!            end do
!          end do
!         end do
!      endif
!      ! Initial stopping time used to compute half-step stopping time.
!      do i=1,nn
!        nu(i)=(dgr(i)*cs/rd)*sqrt(1.+&
!        &0.22089323345553233*&
!        &(wh(i,1)**2+wh(i,2)**2+wh(i,3)**2)&
!        &/(cs*cs))
!      end do
!      ! Deposit relevant quantities.
!      do ind=1,twotondim
!         do j=1,nn !deposit first order effective dust mass
!            if(ok(j))then
!               uold(indp(j,ind),ivar_dust)=uold(indp(j,ind),ivar_dust)+&
!               &mov(j)*dt*nu(j)*vol(j,ind)/(1.+dt*nu(j))
!            end if
!         end do
!         do idim=1,ndim
!            do j=1,nn ! deposit first order effective dust momentum.
!               if(ok(j))then
!                  uold(indp(j,ind),ivar_dust+idim)=uold(indp(j,ind),ivar_dust+idim)&
!                  &+mov(j)*(v(j,idim)+big_v(j,ind,idim))*dt*nu(j)*vol(j,ind)/(1.+dt*nu(j))
!               end if
!            end do
!         end do
!      end do
!      ! Do a half-step update here in order to compute the w**(n+1/2) drift mag.
!      wh(1:nn,1:ndim)=0.0D0
!      do ind=1,twotondim
!        do idim=1,ndim
!          do j=1,nn
!            wh(j,idim)=wh(j,idim)+vol(j,ind)*&
!            &((v(j,idim)+big_v(j,ind,idim)-unew(indp(j,ind),ivar_dust+idim))/(1.+dt*nu(j))-&
!            &(uold(indp(j,ind),ivar_dust+idim)-uold(indp(j,ind),ivar_dust)*unew(indp(j,ind),ivar_dust+idim))&
!            &/((1.+dt*nu(j))*max(uold(indp(j,ind),1)+uold(indp(j,ind),ivar_dust),smallr)))
!          end do
!        end do
!      end do
!     do i=1,nn
!        nu(i)=(dgr(i)*cs/rd)*sqrt(1.+&
!        &0.22089323345553233*(wh(i,1)**2+wh(i,2)**2+wh(i,3)**2)/(cs*cs)) ! In principle, will also depend on the sound speed.
!     end do
!   endif
! end subroutine StoppingRateMidpt
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
