subroutine move_fine(ilevel)
  use amr_commons
  use pm_commons
  use mpi_mod
  implicit none
  integer::ilevel
  !----------------------------------------------------------------------
  ! Update particle position and time-centred velocity at level ilevel.
  ! If particle sits entirely in level ilevel, then use fine grid force
  ! for CIC interpolation. Otherwise, use coarse grid (ilevel-1) force.
  !----------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,next_part,ig,ip,npart1
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part
  character(LEN=80)::filename,fileloc
  character(LEN=5)::nchar

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  filename='trajectory.dat'
  call title(myid,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)
  open(25+myid, file = fileloc, status = 'unknown', access = 'append')

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
              call move1(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
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
  if(ip>0)call move1(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)

  close(25+myid)

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
  integer::ilevel
  !----------------------------------------------------------------------
  ! Update particle position and time-centred velocity at level ilevel.
  ! If particle sits entirely in level ilevel, then use fine grid force
  ! for CIC interpolation. Otherwise, use coarse grid (ilevel-1) force.
  !----------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,next_part,ig,ip,npart1,npart2
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

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
              call move1(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
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
  if(ip>0)call move1(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)

111 format('   Entering move_fine for level ',I2)

end subroutine move_fine_static
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine move1(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use poisson_commons
  use hydro_commons, ONLY: uold,unew,smallr,nvar,gamma
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !------------------------------------------------------------
  ! This routine computes the force on each particle by
  ! inverse CIC and computes new positions for all particles.
  ! If particle sits entirely in fine level, then CIC is performed
  ! at level ilevel. Otherwise, it is performed at level ilevel-1.
  ! This routine is called by move_fine.
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
  ! Particle-based arrays
  logical ,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector,1:ndim),save::x,ff,new_xp,new_vp,dd,dg
  real(dp),dimension(1:nvector,1:ndim),save::vv
  real(dp),dimension(1:10,1:ndim),save ::bb,uu
  real(dp),dimension(1:nvector,1:twotondim,1:ndim),save::big_vv
  real(dp),dimension(1:nvector),save:: nu_stop,mov ! ERM: fluid density interpolated to grain pos. and stopping times
  integer ,dimension(1:nvector,1:ndim),save::ig,id,igg,igd,icg,icd
  real(dp),dimension(1:nvector,1:twotondim),save::vol
  integer ,dimension(1:nvector,1:twotondim),save::igrid,icell,indp,kg
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
  vol_loc=dx_loc**3

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

  ! Check for illegal moves
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
     do ind=1,twotondim
        do idim=1,ndim
           do j=1,np
              ff(j,idim)=ff(j,idim)+uold(indp(j,ind),idim+1)/max(uold(indp(j,ind),1),smallr)*vol(j,ind)
           end do
        end do
     end do
  endif
  ! Gather 3-force
  if(poisson)then
     do ind=1,twotondim
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

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! LORENTZ KICK
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! ERM: Determining the evolved versions of the drift, then interpolating those
  ! For now, do the 8N calculations, but in the future, would be good to loop
  ! over CELLS rather than particles (esp. when doing TSC).

  if(boris.and.hydro)then
    ! Set unew's dust momentum slots to be the gas velocity.

    call ResetUnewToFluidVel()
    big_vv(1:np,1:twotondim,1:ndim)=0.0D0 ! collects velocity changes to sub-clouds
    vv(1:np,1:ndim)=new_vp(1:np,1:ndim)
    call EMKick(np,dtnew(ilevel),indp,ctm,ok,vol,mov,vv,big_vv)
    ! big_vv now contains changes to sub-cloud velocities. vv is still the old
    ! velocity. As well, unew's dust slot contains u^n+du^EM
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! DRAG KICK
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Compute the stopping rates at a half-timestep advanced.
  ! If the stopping time is a constant, the array will be constant.
  if (boris.and.hydro)then
    ! Reset dust variables to zero
    call ResetUoldToZero()
    call StoppingRate(np,dtnew(ilevel),indp,ok,vol,mov,vv,big_vv,nu_stop)
    ! Compute stopping rates at time n+1/2 using backward Euler.
    ! If we have a constant stopping rate, it just sets nu_stop=1/t_stop.

    call ResetUoldToZero()
    vv(1:np,1:ndim)=new_vp(1:np,1:ndim)
    call DragKick(np,dtnew(ilevel),indp,ok,vol,mov,nu_stop,big_vv,vv)
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

  if(boris.or.tracer)then
     do index_part=1,10
        do j=1,np
           if(idp(ind_part(j)).EQ.index_part)then
              write(25+myid,*)t-dtnew(ilevel),idp(ind_part(j)),& ! Old time
                   & xp(ind_part(j),1),xp(ind_part(j),2),xp(ind_part(j),3),& ! Old particle position
                   & vp(ind_part(j),1),vp(ind_part(j),2),vp(ind_part(j),3) ! Old particle velocity
                   ! &  uu(j,1),uu(j,2),uu(j,3),& ! Old fluid velocity
                   ! &  bb(j,1),bb(j,2),bb(j,3),& ! Old magnetic field.
                   ! & new_vp(j,1),new_vp(j,2),new_vp(j,3) ! NEW particle velocity (for comparison)
           end if
        end do
     end do
  endif
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
        xp(ind_part(j),idim)=new_xp(j,idim)
     end do
  end do

  ! Store new velocity
  do idim=1,ndim
     do j=1,np
        vp(ind_part(j),idim)=new_vp(j,idim)
     end do
  end do

  ! Deposit minus final dust momentum to new gas momentum
  do ind=1,twotondim
     do idim=1,ndim
        do j=1,np
           if(ok(j))then
              unew(indp(j,ind),1+idim)=unew(indp(j,ind),1+idim)-mov(j)*vp(ind_part(j),idim)*vol(j,ind)
           end if
        end do
     end do
  end do


end subroutine move1
!#########################################################################
!#########################################################################

!#########################################################################
!#########################################################################

subroutine EMKick(nn,dt,indp,ctm,ok,vol,mov,v,big_v)
  ! This subroutine will compute changes to sub-cloud velocity in big_v,
  ! as well as set unew's dust momentum slot to being u+du^EM.
  use amr_parameters
  use hydro_parameters
  use hydro_commons, ONLY: uold,unew,smallr,nvar,gamma
  implicit none
  integer ::ivar_dust ! cell-centered dust variables start.
  integer ::nn ! number of cells
  real(dp) ::dt! timestep
  real(dp) ::ctm ! charge-to-mass ratio
  logical ,dimension(1:nvector)::ok
  real(dp),dimension(1:nvector)::mov
  real(dp),dimension(1:nvector,1:twotondim)::vol
  integer ,dimension(1:nvector,1:twotondim)::indp
  real(dp),dimension(1:nvector,1:twotondim,1:ndim)::big_v
  real(dp),dimension(1:nvector,1:ndim) ::v,u! grain velocity,magnetic field,velocity
  real(dp) ::w1,w2,w3,B1,B2,B3,den_dust,den_gas,mu
  real(dp),dimension(1:3) ::vtemp
  integer ::i,j,ind,idim! Just an -index

  ivar_dust=9

  do ind=1,twotondim
     do i=1,nn
        den_gas=uold(indp(i,ind),1)
        den_dust=uold(indp(i,ind),ivar_dust)
        mu=den_dust/max(den_gas,smallr)
        B1=0.5D0*(uold(indp(i,ind),1+5)+uold(indp(i,ind),1+nvar))
        B2=0.5D0*(uold(indp(i,ind),2+5)+uold(indp(i,ind),2+nvar))
        B3=0.5D0*(uold(indp(i,ind),3+5)+uold(indp(i,ind),3+nvar))
        w1=uold(indp(i,ind),1+ivar_dust)/max(uold(indp(i,ind),ivar_dust),smallr)&
        &-uold(indp(i,ind),1+1)/max(uold(indp(i,ind),1),smallr)
        w2=uold(indp(i,ind),2+ivar_dust)/max(uold(indp(i,ind),ivar_dust),smallr)&
        &-uold(indp(i,ind),2+1)/max(uold(indp(i,ind),1),smallr)
        w3=uold(indp(i,ind),3+ivar_dust)/max(uold(indp(i,ind),ivar_dust),smallr)&
        &-uold(indp(i,ind),3+1)/max(uold(indp(i,ind),1),smallr)

        big_v(i,ind,1)=& ! First contains gas momentum changes
        &((uold(indp(i,ind),1+ivar_dust)-mu*uold(indp(i,ind),1+1))-&
        &den_dust*& !dust mass
        &(((4. + (B1*B1 - B2*B2 - B3*B3)*ctm*(1.+mu)*ctm*(1.+mu)*dt*dt)*w1 &
        &+ 2.*ctm*(1.+mu)*dt*(2.*B3*w2 + &
        &B1*B2*ctm*(1.+mu)*dt*w2 - 2.*B2*w3 + B1*B3*ctm*(1.+mu)*dt*w3))&
        &/(4. + (B1*B1 + B2*B2+ B3*B3)*ctm*(1.+mu)*ctm*(1.+mu)*dt*dt)))&
        &/max(den_dust+den_gas,smallr)

        big_v(i,ind,2)=&
        &((uold(indp(i,ind),2+ivar_dust)-mu*uold(indp(i,ind),2+1))-&
        &den_dust*& !dust mass
        &(((4. + (B2*B2 - B3*B3 - B1*B1)*ctm*(1.+mu)*ctm*(1.+mu)*dt*dt)*w2 &
        &+ 2.*ctm*(1.+mu)*dt*(2.*B1*w3 + &
        &B2*B3*ctm*(1.+mu)*dt*w3 - 2.*B3*w1 + B2*B1*ctm*(1.+mu)*dt*w1))&
        &/(4. + (B1*B1 + B2*B2+ B3*B3)*ctm*(1.+mu)*ctm*(1.+mu)*dt*dt)))&
        &/max(den_dust+den_gas,smallr)

        big_v(i,ind,3)=&
        &((uold(indp(i,ind),3+ivar_dust)-mu*uold(indp(i,ind),3+1))-&
        &den_dust*& !dust mass
        &(((4. + (B3*B3 - B1*B1 - B2*B2)*ctm*(1.+mu)*ctm*(1.+mu)*dt*dt)*w3 &
        &+ 2.*ctm*(1.+mu)*dt*(2.*B2*w1 + &
        &B3*B1*ctm*(1.+mu)*dt*w1 - 2.*B1*w2 + B3*B2*ctm*(1.+mu)*dt*w2))&
        &/(4. + (B1*B1 + B2*B2+ B3*B3)*ctm*(1.+mu)*ctm*(1.+mu)*dt*dt)))&
        &/max(den_dust+den_gas,smallr)

        do idim=1,ndim
          vtemp(idim) = v(i,idim)-uold(indp(i,ind),1+idim)/uold(indp(i,ind),1)&
          &-0.5*big_v(i,ind,idim)
        end do

        big_v(i,ind,1)=& ! Now let it be the change in sub-cloud velocity.
        &(2.*ctm*ctm*dt*dt*(B1*B2*vtemp(2) + B1*B3*vtemp(3) &
        &- (B3*B3 + B2*B2)*vtemp(1)) &
        &+4.*ctm*dt*(B3*vtemp(2) - B2*vtemp(3)))/&
        &(4. + (B1*B1 + B2*B2 + B3*B3)*ctm*ctm*dt*dt)
        big_v(i,ind,2)=&
        &(2.*ctm*ctm*dt*dt*(B2*B3*vtemp(3) + B2*B1*vtemp(1) &
        &- (B1*B1 + B3*B3)*vtemp(2)) &
        &+4.*ctm*dt*(B1*vtemp(3) - B3*vtemp(1)))/&
        &(4. + (B1*B1 + B2*B2 + B3*B3)*ctm*ctm*dt*dt)
        big_v(i,ind,3)=&
        &(2.*ctm*ctm*dt*dt*(B3*B1*vtemp(1) + B3*B2*vtemp(2) &
        &- (B2*B2 + B1*B1)*vtemp(3)) &
        &+4.*ctm*dt*(B2*vtemp(1) - B1*vtemp(2)))/&
        &(4. + (B1*B1 + B2*B2 + B3*B3)*ctm*ctm*dt*dt)
     end do
  end do

  ! do ind=1,twotondim
  !    do idim=1,ndim
  !       do j=1,nn ! Add up velocity changes to for particles.
  !          v(j,idim)=v(j,idim)+big_vv(j,ind,idim)*vol(j,ind)
  !       end do
  !    end do
  ! end do

  ! Deposit dust sub-cloud momentum onto the unew ``dust'' slot.
  ! This is still velocity here.
  do ind=1,twotondim
     do idim=1,ndim
        do j=1,nn
           if(ok(j))then
              unew(indp(j,ind),ivar_dust+idim)=unew(indp(j,ind),ivar_dust+idim)-&
              &mov(j)*big_v(j,ind,idim)*vol(j,ind)/uold(indp(j,ind),1)
           end if
        end do
     end do
  end do
end subroutine EMKick

!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine StoppingRate(nn,twodt,indp,ok,vol,mov,v,big_v,nu)
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
  integer ::ivar_dust ! cell-centered dust variables start.
  real(dp) ::dt,twodt ! half-timestep, full timestep.
  real(dp)::rd,cs! ERM: Grain size parameter
  real(dp),dimension(1:nvector) ::nu
  real(dp),dimension(1:nvector) ::mov
  real(dp),dimension(1:nvector,1:twotondim)::vol
  integer ,dimension(1:nvector,1:twotondim)::indp
  logical ,dimension(1:nvector)::ok
  real(dp),dimension(1:nvector),save ::dgr! gas density at grain.
  real(dp),dimension(1:nvector,1:ndim) ::v! grain velocity
  real(dp),dimension(1:nvector,1:twotondim,1:ndim)::big_v
  real(dp),dimension(1:nvector,1:ndim),save ::wh! drift at half step.
  integer ::i,j,idim,ind
  ivar_dust=9
  dt=0.5*twodt
  rd = sqrt(gamma)*0.62665706865775*grain_size !constant for epstein drag law.
  cs=1.0 ! isothermal sound speed... Need to get this right. This works for now,
         ! but only if you have scaled things so that the sound speed is 1.

  if (constant_t_stop)then
    nu(1:nvector)=1./t_stop
  else
     dgr(1:nn) = 0.0D0
     if(boris)then
        do ind=1,twotondim
            do j=1,nn
               dgr(j)=dgr(j)+uold(indp(j,ind),1)*vol(j,ind)
           end do
        end do
     endif

     wh(1:nn,1:ndim) = 0.0D0 ! Set to the drift velocity post-Lorentz force
     if(boris)then
        do ind=1,twotondim
          do idim=1,ndim
            do j=1,nn
               wh(j,idim)=wh(j,idim)+&
               &(v(j,idim)+big_v(j,ind,idim)-unew(indp(j,ind),ivar_dust+idim))&
               &*vol(j,ind)
           end do
         end do
        end do
     endif
     ! Initial stopping time used to compute half-step stopping time.
     do i=1,nn
       nu(i)=(dgr(i)*cs/rd)*sqrt(1.+&
       &0.22089323345553233*&
       &(wh(i,1)**2+wh(i,2)**2+wh(i,3)**2)&
       &/(cs*cs))
     end do
     ! Deposit relevant quantities.
     do ind=1,twotondim
        do j=1,nn !deposit first order effective dust mass
           if(ok(j))then
              uold(indp(j,ind),ivar_dust)=uold(indp(j,ind),ivar_dust)+&
              &mov(j)*dt*nu(j)*vol(j,ind)/(1.+dt*nu(j))
           end if
        end do
        do idim=1,ndim
           do j=1,nn ! deposit first order effective dust momentum.
              if(ok(j))then
                 uold(indp(j,ind),ivar_dust+idim)=uold(indp(j,ind),ivar_dust+idim)&
                 &+mov(j)*(v(j,idim)+big_v(j,ind,idim))*dt*nu(j)*vol(j,ind)/(1.+dt*nu(j))
              end if
           end do
        end do
     end do
     ! Do a half-step update here in order to compute the w^(n+1/2) drift mag.
     wh(1:nn,1:ndim)=0.0D0
     do ind=1,twotondim
       do idim=1,ndim
         do j=1,nn
           wh(j,idim)=wh(j,idim)+vol(j,ind)*&
           &((v(j,idim)+big_v(j,ind,idim)-unew(indp(j,ind),ivar_dust+idim))/(1.+dt*nu(j))-&
           &(uold(indp(j,ind),ivar_dust+idim)-uold(indp(j,ind),ivar_dust)*unew(indp(j,ind),ivar_dust+idim))&
           &/((1.+dt*nu(j))*max(uold(indp(j,ind),1)+uold(indp(j,ind),ivar_dust),smallr)))
         end do
       end do
     end do
    do i=1,nn
       nu(i)=(dgr(i)*cs/rd)*sqrt(1.+&
       &0.22089323345553233*(wh(i,1)**2+wh(i,2)**2+wh(i,3)**2)/(cs*cs)) ! In principle, will also depend on the sound speed.
    end do
  endif
end subroutine StoppingRate
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################

subroutine DragKick(nn,dt,indp,ok,vol,mp,nu,big_v,v) ! mp is actually mov
  use amr_parameters
  use hydro_parameters
  use hydro_commons, ONLY: uold,unew,smallr,nvar,gamma
  implicit none
  integer::nn
  integer ::ivar_dust ! cell-centered dust variables start.  integer ::nn ! number of cells
  real(dp) ::dt ! timestep
  real(dp) ::vol_loc ! cloud volume
  real(dp),dimension(1:nvector) ::nu,mp
  logical ,dimension(1:nvector)::ok
  real(dp),dimension(1:nvector,1:twotondim)::vol
  integer ,dimension(1:nvector,1:twotondim)::indp
  real(dp),dimension(1:nvector,1:ndim) ::v ! grain velocity
  real(dp),dimension(1:nvector,1:twotondim,1:ndim) ::big_v
  real(dp),dimension(1:nvector,1:ndim),save ::vo ! grain velocity "new" (?)
  real(dp) ::den_dust,den_gas,mu
  integer ::i,j,ind,idim! Just an index
  ivar_dust=9

  do ind=1,twotondim
    !deposit effective dust mass
     do j=1,nn
        if(ok(j))then
           uold(indp(j,ind),ivar_dust)=uold(indp(j,ind),ivar_dust)+&
           &mp(j)*0.5*(dt*nu(j)+dt*dt*nu(j)*nu(j))*vol(j,ind)/&
           &(1.+dt*nu(j)+0.5*dt*dt*nu(j)*nu(j))
        end if
     end do
     do idim=1,ndim
        do j=1,nn ! deposit vector for change in gas velocity. (V0)
           if(ok(j))then
              uold(indp(j,ind),ivar_dust+idim)=uold(indp(j,ind),ivar_dust+idim)&
              &+mp(j)*(dt*nu(j)+0.5*dt*dt*nu(j)*nu(j))*&
              &(v(j,idim)+big_v(j,ind,idim)-unew(indp(j,ind),ivar_dust+idim))*&
              &vol(j,ind)/((1.+dt*nu(j)+0.5*dt*dt*nu(j)*nu(j))*&
              &max(uold(indp(j,ind),1)+uold(indp(j,ind),ivar_dust),smallr))
           end if
        end do
     end do
  end do

  ! Now kick the dust given these quantities.
  do ind=1,twotondim
    do idim=1,ndim
      do j=1,nn
        v(j,idim)=v(j,idim)+vol(j,ind)*&
        &(big_v(j,ind,idim)+((dt*nu(j)+0.5*dt*dt*nu(j)*nu(j))*&
        &(unew(indp(j,ind),ivar_dust+idim)-v(j,idim)-big_v(j,ind,idim))+&
        &0.5*(dt*nu(j)+dt*dt*nu(j)*nu(j))*uold(indp(j,ind),ivar_dust+idim)&
        &)/(1.+dt*nu(j)+0.5*dt*dt*nu(j)*nu(j)))
      end do
    end do
  end do
end subroutine DragKick
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine ResetUoldToZero()
  use pm_commons
  use amr_commons
  use hydro_commons
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::info
#endif
  integer::ilevel
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

  do icpu=1,ncpu
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do ivar=ivar_dust,ivar_dust+ndim
           do i=1,reception(icpu,ilevel)%ngrid
              uold(reception(icpu,ilevel)%igrid(i)+iskip,ivar)=0.0D0
           end do
        end do
     end do
  end do

  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do ivar=ivar_dust,ivar_dust+ndim
        do i=1,active(ilevel)%ngrid
           uold(active(ilevel)%igrid(i)+iskip,ivar)=0.0D0
        end do
     end do
  end do

  ! Reset rho in physical boundaries
  do ibound=1,nboundary
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do ivar=ivar_dust,ivar_dust+ndim
           do i=1,boundary(ibound,ilevel)%ngrid
              uold(boundary(ibound,ilevel)%igrid(i)+iskip,ivar)=0.0D0
           end do
        end do
     end do
  end do
end subroutine ResetUoldToZero
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine ResetUnewToFluidVel()
  use pm_commons
  use amr_commons
  use hydro_commons
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::info
#endif
  integer::ilevel
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
  ivar_dust=9
  do icpu=1,ncpu
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do ivar=1,ndim
           do i=1,reception(icpu,ilevel)%ngrid
              unew(reception(icpu,ilevel)%igrid(i)+iskip,ivar+ivar_dust)=&
              &uold(reception(icpu,ilevel)%igrid(i)+iskip,ivar+1)/&
              &max(uold(active(ilevel)%igrid(i)+iskip,1),smallr)
           end do
        end do
     end do
  end do

  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do ivar=1,ndim
        do i=1,active(ilevel)%ngrid
           unew(active(ilevel)%igrid(i)+iskip,ivar+ivar_dust)=&
           &uold(active(ilevel)%igrid(i)+iskip,1+ivar)/&
           &max(uold(active(ilevel)%igrid(i)+iskip,1),smallr)
        end do
     end do
  end do

  ! Reset rho in physical boundaries
  do ibound=1,nboundary
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do ivar=1,ndim
           do i=1,boundary(ibound,ilevel)%ngrid
              unew(boundary(ibound,ilevel)%igrid(i)+iskip,ivar+ivar_dust)=&
              &uold(boundary(ibound,ilevel)%igrid(i)+iskip,ivar+1)/&
              &max(uold(active(ilevel)%igrid(i)+iskip,1),smallr)
           end do
        end do
     end do
  end do
end subroutine ResetUnewToFluidVel

!#########################################################################
!#########################################################################
!#########################################################################
