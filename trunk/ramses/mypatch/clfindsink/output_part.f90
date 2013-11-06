subroutine backup_part(filename)
  use amr_commons
  use pm_commons
  implicit none
  character(LEN=80)::filename

  integer::i,idim,ilun,ipart
  character(LEN=80)::fileloc
  character(LEN=5)::nchar
  real(dp),allocatable,dimension(:)::xdp
  integer ,allocatable,dimension(:)::ii
  integer ,allocatable,dimension(:)::ll

  if(verbose)write(*,*)'Entering backup_part'
  
  ilun=2*ncpu+myid+10

  call title(myid,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)
  open(unit=ilun,file=TRIM(fileloc),form='unformatted')
  rewind(ilun)
  ! Write header
  write(ilun)ncpu
  write(ilun)ndim
  write(ilun)npart
  write(ilun)localseed
  write(ilun)nstar_tot   
  write(ilun)mstar_tot   
  write(ilun)mstar_lost
  write(ilun)nsink
  ! Write position
  allocate(xdp(1:npart))
  do idim=1,ndim
     ipart=0
     do i=1,npartmax
        if(levelp(i)>0)then
           ipart=ipart+1
           xdp(ipart)=xp(i,idim)
        end if
     end do
     write(ilun)xdp
  end do
  ! Write velocity
  do idim=1,ndim
     ipart=0
     do i=1,npartmax
        if(levelp(i)>0)then
           ipart=ipart+1
           xdp(ipart)=vp(i,idim)
        end if
     end do
     write(ilun)xdp
  end do
  ! Write mass
  ipart=0
  do i=1,npartmax
     if(levelp(i)>0)then
        ipart=ipart+1
        xdp(ipart)=mp(i)
     end if
  end do
  write(ilun)xdp
  deallocate(xdp)
  ! Write identity
  allocate(ii(1:npart))
  ipart=0
  do i=1,npartmax
     if(levelp(i)>0)then
        ipart=ipart+1
        ii(ipart)=idp(i)
     end if
  end do
  write(ilun)ii
  deallocate(ii)
  ! Write level
  allocate(ll(1:npart))
  ipart=0
  do i=1,npartmax
     if(levelp(i)>0)then
        ipart=ipart+1
        ll(ipart)=levelp(i)
     end if
  end do
  write(ilun)ll
  deallocate(ll)
  ! Write birth epoch
  if(star.or.sink)then
     allocate(xdp(1:npart))
     ipart=0
     do i=1,npartmax
        if(levelp(i)>0)then
           ipart=ipart+1
           xdp(ipart)=tp(i)
        end if
     end do
     write(ilun)xdp
     ! Write metallicity
     if(metal)then
     ipart=0
     do i=1,npartmax
        if(levelp(i)>0)then
           ipart=ipart+1
           xdp(ipart)=zp(i)
        end if
     end do
     write(ilun)xdp
     end if
     deallocate(xdp)
  end if
  !======================
  ! Write sink properties
  !======================
  if(sink)then
     if(nsink>0)then
     allocate(xdp(1:nsink))
     do i=1,nsink
        xdp(i)=msink(i)
     end do
     write(ilun)xdp ! Write sink mass
     do i=1,nsink
        xdp(i)=tsink(i)
     end do
     write(ilun)xdp ! Write sink birth epoch
     do idim=1,ndim
        do i=1,nsink
           xdp(i)=xsink(i,idim)
        end do
        write(ilun)xdp ! Write sink position
     enddo
     do idim=1,ndim
        do i=1,nsink
           xdp(i)=vsink(i,idim)
        end do
        write(ilun)xdp ! Write sink velocity
     enddo
     do i=1,nsink
        xdp(i)=delta_mass(i)
     end do
     write(ilun)xdp ! Write sink accumulated rest mass energy
     deallocate(xdp)
     allocate(ii(1:nsink))
     do i=1,nsink
        ii(i)=idsink(i)
     end do
     write(ilun)ii ! Write sink index
     deallocate(ii)
     endif
  endif

  close(ilun)

end subroutine backup_part

subroutine output_sink(filename)
  use amr_commons
  use pm_commons
  implicit none
  character(LEN=80)::filename

  integer::i,idim,ipart,isink
  integer::nx_loc,ny_loc,nz_loc,ilun,icpu,idom
  real(dp)::scale
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m
  character(LEN=80)::fileloc
  character(LEN=5)::nchar

  if(verbose)write(*,*)'Entering output_amr'

  ilun=myid+10

  ! Conversion factor from user units to cgs units                                                                   
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_m=scale_d*scale_l**3d0

  if(verbose)write(*,*)'Entering output_sink'
  
  ilun=2*ncpu+myid+10

  fileloc=TRIM(filename)
  open(unit=ilun,file=TRIM(fileloc),form='formatted')
  !======================
  ! Write sink properties
  !======================
  write(ilun,*)'Number of sink = ',nsink
  write(ilun,'(" ========================================================================================================================================= ")')
  write(ilun,'(" Id     Mass(Msol)              x              y              z             vx             vy             vz      new_born ")')
  write(ilun,'(" ========================================================================================================================================= ")')
  do isink=1,nsink
     write(ilun,'(I6,7(1X,1PE14.7),X,X,X,X,X,I1)')idsink(isink),msink(isink)*scale_m/2d33,xsink(isink,1:ndim),vsink(isink,1:ndim),new_born_all(isink)
  end do
  write(ilun,'(" ========================================================================================================================================= ")')

  close(ilun)

end subroutine output_sink



