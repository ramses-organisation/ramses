
subroutine writemem(usedmem)
  real(kind=4)::usedmem

  usedmem=real(usedmem)*4096

  if(usedmem>1024.**4.)then
     write(*,999)usedmem/1024.**4.
  else if (usedmem>1024.**3.) then
     write(*,998)usedmem/1024.**3.
  else if (usedmem>1024.**2.) then
     write(*,997)usedmem/1024.**2.
  else if (usedmem>1024.) then
     write(*,996)usedmem/1024.
  endif

996 format(' Used memory:',F9.1,' kB')
997 format(' Used memory:',F9.1,' MB')
998 format(' Used memory:',F9.3,' GB')
999 format(' Used memory:',F9.3,' TB')

end subroutine writemem

subroutine getmem(outmem)
  use amr_commons,only:myid
#ifndef WITHOUTMPI
  use amr_commons,only:IOGROUPSIZE
  use amr_commons,only:ncpu
#endif
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::dummy_io,info2
#endif
  real(kind=4)::outmem
  character(len=300) :: dir, dir2, file
  integer::read_status
  integer,parameter::tag=1134
  integer::nmem,ind,j
  logical::file_exists

  file='/proc/self/stat'
#ifndef WITHOUTMPI
  if(IOGROUPSIZE>0) then
     if (mod(myid-1,IOGROUPSIZE)/=0) then
        call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
             & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
     end if
  endif
#endif

  inquire(file=file, exist=file_exists)
  if (file_exists) then
     open(unit=12,file=file,form='formatted')
     read(12,'(A300)',IOSTAT=read_status)dir
     close(12)
  else
     read_status=-1000
  endif

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


  if (read_status < 0)then
     outmem=0
     if (myid==1 .and. read_status .ne. -1000)write(*,*)'Problem in checking free memory'
  else
     ind=300
     j=0
     do while (j<23)
        ind=index(dir,' ')
        dir2=dir(ind+1:300)
        j=j+1
        dir=dir2
     end do
     ind=index(dir,' ')
     dir2=dir(1:ind)
     read(dir2,'(I12)')nmem
     outmem=real(nmem,kind=4)
  end if

end subroutine getmem
