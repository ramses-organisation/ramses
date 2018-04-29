subroutine backup_radiation(filename)
  use amr_commons
  use hydro_commons
  use radiation_commons, ONLY: Erad
  use mpi_mod
  implicit none
  character(LEN=80)::filename

  integer::i,ivar,ncache,ind,ilevel,igrid,iskip,ilun,istart,ibound,nvar_rad
  integer,allocatable,dimension(:)::ind_grid
  real(dp),allocatable,dimension(:)::xdp
  character(LEN=5)::nchar
  character(LEN=80)::fileloc
  integer,parameter::tag=1124
  integer::dummy_io,info2

  if(verbose)write(*,*)'Entering backup_radiation'

  ilun=ncpu+myid+10
  nvar_rad=1

  ! Wait for the token
#ifndef WITHOUTMPI
  if(IOGROUPSIZE>0) then
     if (mod(myid-1,IOGROUPSIZE)/=0) then
        call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
             & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
     end if
  endif
#endif

  call title(myid,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)
  open(unit=ilun,file=fileloc,form='unformatted')
  write(ilun)ncpu
  write(ilun)nvar_rad
  write(ilun)ndim
  write(ilun)nlevelmax
  write(ilun)nboundary

  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        write(ilun)ilevel
        write(ilun)ncache
        if(ncache>0)then
           allocate(ind_grid(1:ncache),xdp(1:ncache))
           ! Loop over level grids
           igrid=istart
           do i=1,ncache
              ind_grid(i)=igrid
              igrid=next(igrid)
           end do
           ! Loop over cells
           do ind=1,twotondim
              iskip=ncoarse+(ind-1)*ngridmax
              ! Output the radiation photon density.
              do i=1,ncache
                 xdp(i)=Erad(ind_grid(i)+iskip)
              end do
              write(ilun)xdp
           end do
           deallocate(ind_grid, xdp)
        end if
     end do
  end do
  close(ilun)

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


end subroutine backup_radiation

subroutine store_radiation(filename)
  use data_common
  use radiation_commons
  use mpi_mod
  implicit none

  character(LEN=80)::filename
  integer::ilun
  character(LEN=5)::nchar
  character(LEN=80)::fileloc
  integer,parameter::tag=1125
  integer::dummy_io,info2

  ilun=ncpu+myid+10
  call title(myid,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)

  ! Wait for the token
#ifndef WITHOUTMPI
  if(IOGROUPSIZE>0) then
     if (mod(myid-1,IOGROUPSIZE)/=0) then
        call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
             & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
     end if
  endif
#endif

  open(unit=ilun,file=fileloc,form='unformatted')
  write(ilun)grid_size_x,grid_size_y,grid_size_z
  write(ilun)cpu_e
  write(ilun)cpu_f
  close(ilun)

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

end subroutine store_radiation

subroutine restore_radiation
  use data_common
  use radiation_commons
  use mpi_mod
  implicit none

  integer::ilun
  character(LEN=5)::nchar,ncharcpu
  character(LEN=80)::fileloc
  integer::x,y,z
  integer,parameter::tag=1126
  integer::dummy_io,info2

  ilun=ncpu+myid+10
  call title(nrestart,nchar)

  if(IOGROUPSIZEREP>0)then
     call title(((myid-1)/IOGROUPSIZEREP)+1,ncharcpu)
     fileloc='output_'//TRIM(nchar)//'/group_'//TRIM(ncharcpu)//'/radgpu_'//TRIM(nchar)//'.out'
  else
     fileloc='output_'//TRIM(nchar)//'/radgpu_'//TRIM(nchar)//'.out'
  endif
 call title(myid,nchar)
  fileloc=TRIM(fileloc)//TRIM(nchar)

  if (myid.eq.1) then
     write(*,*)"Restoring radiation data: ", fileloc
  end if

   ! Wait for the token
#ifndef WITHOUTMPI
  if(IOGROUPSIZE>0) then
     if (mod(myid-1,IOGROUPSIZE)/=0) then
        call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
             & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
     end if
  endif
#endif


  open(unit=ilun,file=fileloc,form='unformatted')
  read(ilun)x,y,z
  if ((x.ne.grid_size_x).or.(y.ne.grid_size_y).or.(z.ne.grid_size_z)) then
     write(*,*)"Error: Radiation grid size mismatch from ",fileloc
     write(*,*)"  expected:",grid_size_x,grid_size_y,grid_size_z
     write(*,*)"  actual:",x,y,z
     close(ilun)
     call clean_stop
  end if
  read(ilun)cpu_e
  read(ilun)cpu_f
  close(ilun)

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

end subroutine restore_radiation
