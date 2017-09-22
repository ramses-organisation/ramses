subroutine backup_part(filename, file_desc)
  use amr_commons
  use pm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  character(LEN=80), intent(in)::filename, file_desc

  integer::i,idim,ilun,ipart
  character(LEN=80)::fileloc
  character(LEN=5)::nchar
  real(dp),allocatable,dimension(:)::xdp
  integer(i8b),allocatable,dimension(:)::ii8
  integer(2),allocatable,dimension(:)::ishort
  integer,allocatable,dimension(:)::ll
  integer,parameter::tag=1122
  integer::dummy_io,info2
  integer(1),allocatable,dimension(:)::ii1

  integer :: unit_info, ivar

  character(len=1),dimension(1:3) :: dim_keys = (/"x", "y", "z"/)

  if(verbose)write(*,*)'Entering backup_part'

  ! Set ivar to 1 for first variable
  ivar = 1

  ! Wait for the token
#ifndef WITHOUTMPI
  if(IOGROUPSIZE>0) then
     if (mod(myid-1,IOGROUPSIZE)/=0) then
        call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
             & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
     end if
  endif
#endif


  ilun=2*ncpu+myid+10

  call title(myid,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)
  open(unit=ilun,file=TRIM(fileloc),form='unformatted')
  if (myid == 1) then
     open(newunit=unit_info,file=trim(file_desc),form='formatted')
     write(unit_info, *) 'version: ', 1
  end if

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
     if (myid == 1) call dump_var("position_"//dim_keys(idim), ivar, unit_info)
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
     if (myid == 1) call dump_var("velocity_"//dim_keys(idim), ivar, unit_info)
  end do
  ! Write mass
  ipart=0
  do i=1,npartmax
     if(levelp(i)>0)then
        ipart=ipart+1
        xdp(ipart)=mp(i)
     end if
  end do
  if (myid == 1) call dump_var("mass", ivar, unit_info)
  write(ilun)xdp
  deallocate(xdp)
  ! Write identity
  allocate(ii8(1:npart))
  ipart=0
  do i=1,npartmax
     if(levelp(i)>0)then
        ipart=ipart+1
        ii8(ipart)=idp(i)
     end if
  end do
  write(ilun)ii8
  if (myid == 1) call dump_var("identity", ivar, unit_info)
  deallocate(ii8)

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
  if (myid == 1) call dump_var("levelp", ivar, unit_info)

  deallocate(ll)

  ! Write family
  allocate(ii1(1:npart))
  ipart=0
  do i=1,npartmax
     if(levelp(i)>0)then
        ipart=ipart+1
        ii1(ipart)=typep(i)%family
     end if
  end do
  write(ilun)ii1
  if (myid == 1) call dump_var("family", ivar, unit_info)

  ! Write tag
  ipart=0
  do i=1,npartmax
     if(levelp(i)>0)then
        ipart=ipart+1
        ii1(ipart)=typep(i)%tag
     end if
  end do
  write(ilun)ii1
  if (myid == 1) call dump_var("tag", ivar, unit_info)

  deallocate(ii1)

#ifdef OUTPUT_PARTICLE_POTENTIAL
  ! Write potential (added by AP)
  allocate(xdp(1:npart))
  ipart=0
  do i=1, npartmax
     if(levelp(i)>0) then
        ipart=ipart+1
        xdp(ipart)=ptcl_phi(i)
     end if
  end do
  write(ilun)xdp
  if (myid == 1) call dump_var("potential", ivar, unit_info)

  deallocate(xdp)
#endif

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
     if (myid == 1) call dump_var("birth_time", ivar, unit_info)
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
        if (myid == 1) call dump_var("metallicity", ivar, unit_info)
     end if
     deallocate(xdp)
  end if

  close(ilun)
  if (myid == 1) close(unit_info)

  ! Send the token
#ifndef WITHOUTMPI
  if(IOGROUPSIZE>0) then
     if(mod(myid,IOGROUPSIZE)/=0 .and.(myid.lt.ncpu))then
        dummy_io=1
        call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag, &
             & MPI_COMM_WORLD,info2)
     end if
  endif
#endif

contains
  subroutine dump_var(varname, ivar, unit)
    character(len=*), intent(in) :: varname
    integer, intent(in) :: unit
    integer, intent(inout) :: ivar

    if (myid == 1) &
         write(unit, '("variable #",I2,": ",a)') ivar, trim(varname)

    ivar = ivar + 1
  end subroutine dump_var

end subroutine backup_part
