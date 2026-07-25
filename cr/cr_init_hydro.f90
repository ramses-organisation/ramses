!*************************************************************************
subroutine cr_init_hydro
  ! Allocate the cosmic-ray state arrays (separate from the hydro
  ! uold/unew), initialise them to zero, and read them back on restart.
  use amr_commons
  use cr_parameters
  use cr_hydro_commons
  use hydro_parameters, only: gamma
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::dummy_io,info,info2
  integer,parameter::tag=1133
#endif
  integer::ncell,ncache,iskip,igrid,i,ilevel,ind,igroup,idim
  integer::ilevel2,numbl2,ilun,ibound,istart
  integer::ncpu2,ndim2,nlevelmax2,nboundary2,ncrvar2
  integer ,dimension(:),allocatable::ind_grid
  real(dp),dimension(:),allocatable::xx
  real(dp)::gamma2
  character(LEN=80)::fileloc
  character(LEN=5)::nchar,ncharcpu
  logical::ok

  if(verbose)write(*,*)'Entering cr_init_hydro'

  ncell=ncoarse+twotondim*ngridmax
  allocate(cruold(1:ncell,1:ncrvar))
  allocate(crunew(1:ncell,1:ncrvar))
  cruold=0d0 ; crunew=0d0

  if(verbose)write(*,*)'Allocate done for',ncrvar,'CR variables'

  if(nrestart==0)return
  if(.not.cr_advect)return

  !--------------------------------
  ! For a restart, read the CR file
  !--------------------------------
  ilun=ncpu+myid+113
  call title(nrestart,nchar)

  if(IOGROUPSIZEREP>0)then
     call title(((myid-1)/IOGROUPSIZEREP)+1,ncharcpu)
     fileloc='output_'//TRIM(nchar)//'/group_'//TRIM(ncharcpu)//'/cr_'//TRIM(nchar)//'.out'
  else
     fileloc='output_'//TRIM(nchar)//'/cr_'//TRIM(nchar)//'.out'
  endif

  call title(myid,nchar)
  fileloc=TRIM(fileloc)//TRIM(nchar)
  inquire(file=fileloc, exist=ok)
  if(.not.ok)then
     if(myid==1)write(*,*) &
          'Could not read CR output, but continuing in case of postprocessing.'
     return
  endif

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
  read(ilun)ncpu2
  read(ilun)ncrvar2
  read(ilun)ndim2
  read(ilun)nlevelmax2
  read(ilun)nboundary2
  read(ilun)gamma2
  if(ncrvar2.gt.ncrvar .and. myid==1)then ! Not ok to drop CR variables
     write(*,*)'File cr.tmp is not compatible'
     write(*,*)'Found ncrvar  =',ncrvar2
     write(*,*)'Expected=',ncrvar
     call clean_stop
  end if
  do ilevel=1,nlevelmax2
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        read(ilun)ilevel2
        read(ilun)numbl2
        if(numbl2.ne.ncache)then
           write(*,*)'File cr.tmp is not compatible'
           write(*,*)'Found   =',numbl2,' for level ',ilevel2
           write(*,*)'Expected=',ncache,' for level ',ilevel
        end if
        if(ncache>0)then
           allocate(ind_grid(1:ncache))
           allocate(xx(1:ncache))
           ! Loop over level grids
           igrid=istart
           do i=1,ncache
              ind_grid(i)=igrid
              igrid=next(igrid)
           end do
           ! Loop over cells
           do ind=1,twotondim
              iskip=ncoarse+(ind-1)*ngridmax
              do igroup=1,ncr_groups
                 ! CR energy density (raw, no unit conversion)
                 read(ilun)xx
                 do i=1,ncache
                    cruold(ind_grid(i)+iskip,Ecr_idx(igroup))=xx(i)
                 end do
                 ! CR flux components
                 do idim=1,ndim
                    read(ilun)xx
                    do i=1,ncache
                       cruold(ind_grid(i)+iskip,Ecr_idx(igroup)+idim)=xx(i)
                    end do
                 end do
              end do
           end do
           deallocate(ind_grid,xx)
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

#ifndef WITHOUTMPI
  if(debug)write(*,*)'cr.tmp read for processor ',myid
  call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
  if(verbose)write(*,*)'CR backup files read completed'

end subroutine cr_init_hydro
