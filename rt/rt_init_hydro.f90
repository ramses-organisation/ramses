!*************************************************************************
subroutine rt_init_hydro

!  Allocate rt variables, initialise for RT or non-equilibrium chemistry,
!  initialise ionisation states if needed, and read rt vars if restarting.
!-------------------------------------------------------------------------
  use amr_commons
  use hydro_commons
  use rt_hydro_commons
  use rt_parameters
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::dummy_io,info,info2
#endif
  integer::ncell,ncache,iskip,igrid,i,ilevel,ind,ivar
  integer::ilevel2,numbl2,ilun,ibound,istart
  integer::ncpu2,ndim2,nlevelmax2,nboundary2,idim
  integer ,dimension(:),allocatable::ind_grid
  real(dp),dimension(:),allocatable::xx
  real(dp)::gamma2
  character(LEN=80)::fileloc
  character(LEN=5)::nchar,ncharcpu
  integer::nRTvar2=0
  logical::ok
  integer,parameter::tag=1130

  if(verbose)write(*,*)'Entering init_rt'
  !------------------------------------------------------
  ! Allocate conservative, cell-centered variables arrays
  !------------------------------------------------------
  ncell=ncoarse+twotondim*ngridmax
  allocate(rtuold(1:ncell,1:nrtvar))
  allocate(rtunew(1:ncell,1:nrtvar))
  rtuold=smallNp ; rtunew=smallNp

  if(verbose)write(*,*)'Allocate done for nrtvar'

  call rt_init

  if(nrestart .eq. 0) return

  if(rt_is_init_xion) then
     if(myid==1) &
          write(*,*) 'Initializing ionization states from T profile'
     do ilevel=nlevelmax,1,-1
        call rt_init_xion(ilevel)
        call upload_fine(ilevel)
     end do
  end if

  if(.not.rt) return

  !--------------------------------
  ! For a restart, read rt file
  !--------------------------------
  ilun=ncpu+myid+103
  call title(nrestart,nchar)

  if(IOGROUPSIZEREP>0)then
     call title(((myid-1)/IOGROUPSIZEREP)+1,ncharcpu)
     fileloc='output_'//TRIM(nchar)//'/group_'//TRIM(ncharcpu)//'/rt_'//TRIM(nchar)//'.out'
  else
     fileloc='output_'//TRIM(nchar)//'/rt_'//TRIM(nchar)//'.out'
  endif

  call title(myid,nchar)
  fileloc=TRIM(fileloc)//TRIM(nchar)
  inquire(file=fileloc, exist=ok)
  if(.not.ok) then
     if(myid==1) write(*,*) &
          'Could not read RT output, but continuing in case of postprocessing.'
     return ! May be post-processing of normal RAMSES
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
  read(ilun)nrtvar2
  read(ilun)ndim2
  read(ilun)nlevelmax2
  read(ilun)nboundary2
  read(ilun)gamma2
  if(nrtvar2.gt.nrtvar .and. myid==1)then ! Not ok to drop RT variables
     write(*,*)'File rt.tmp is not compatible (1)'
     write(*,*)'Found nrtvar  =',nrtvar2
     write(*,*)'Expected=',nrtvar
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
           write(*,*)'File rt.tmp is not compatible'
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
              ! Loop over RT variables
              do ivar=1,nGroups
                 ! Read photon density in flux units
                 read(ilun)xx
                 do i=1,ncache
                    rtuold(ind_grid(i)+iskip,iGroups(ivar))=xx(i)/rt_c
                 end do
                 ! Read photon flux
                 do idim=1,ndim
                    read(ilun)xx
                    do i=1,ncache
                       rtuold(ind_grid(i)+iskip,iGroups(ivar)+idim)=xx(i)
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
  if(debug)write(*,*)'rt.tmp read for processor ',myid
  call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
  if(verbose)write(*,*)'RT backup files read completed'

end subroutine rt_init_hydro
