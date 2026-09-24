subroutine backup_poisson(filename)
  use amr_commons
  use poisson_commons
  use mpi_mod
  implicit none
  character(LEN=80)::filename

  integer::i,ivar,ncache,ind,ilevel,igrid,iskip,ilun,istart,ibound
  integer,allocatable,dimension(:)::ind_grid
  real(dp),allocatable,dimension(:)::xdp
  character(LEN=5)::nchar
  character(LEN=80)::fileloc

#ifndef WITHOUTMPI
  integer,parameter::tag=1123
  integer::dummy_io,info2
#endif

  if(verbose)write(*,*)'Entering backup_poisson'

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
  write(ilun)ncpu
#ifdef OUTPUT_PARTICLE_DENSITY
  write(ilun)ndim+2
#else
  write(ilun)ndim+1
#endif
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
#ifdef OUTPUT_PARTICLE_DENSITY
              ! Write density
              do i=1,ncache
                 xdp(i)=rho(ind_grid(i)+iskip)
              end do
              write(ilun)xdp
#endif
              ! Write potential
              do i=1,ncache
                 xdp(i)=phi(ind_grid(i)+iskip)
              end do
              write(ilun)xdp
              ! Write force
              do ivar=1,ndim
                 do i=1,ncache
                    xdp(i)=f(ind_grid(i)+iskip,ivar)
                 end do
                 write(ilun)xdp
              end do
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

end subroutine backup_poisson

!################################################################
!################################################################
!################################################################
!################################################################
subroutine backup_safe_mode(filename)
  use amr_commons        ,only:verbose,nlevelmax,ncpu,myid
  use poisson_parameters ,only:mg_safe_mode_reset
  use poisson_commons    ,only:safe_mode,safe_mode_countdown
  implicit none
  character(LEN=80)::filename
  !--------------------------------------------------------------
  ! Output the multigrid safe_mode variable so that it can be
  ! restored on restart.
  !--------------------------------------------------------------
  integer::ilevel,ilun

  if(verbose)write(*,*)'Entering backup_safe_mode'

  ilun=ncpu+myid+10

  open(unit=ilun,file=TRIM(filename),form='formatted')
  write(ilun,'(A)')'# Multigrid safe mode state'
  write(ilun,'(A8,i5)')'levelmax',nlevelmax
  write(ilun,'(A18,1X,i8)')'mg_safe_mode_reset',mg_safe_mode_reset
  write(ilun,'(A)')'# level  safe_mode  solves_left'
  do ilevel=1,nlevelmax
     write(ilun,'(i5,L9,i12)')ilevel,safe_mode(ilevel),safe_mode_countdown(ilevel)
  end do
  close(ilun)

end subroutine backup_safe_mode
!################################################################
!################################################################
!################################################################
!################################################################
subroutine restore_safe_mode
  use amr_commons
  use poisson_parameters ,only:mg_safe_mode_reset
  use poisson_commons    ,only:safe_mode,safe_mode_countdown
  use mpi_mod
  implicit none
  !--------------------------------------------------------------
  ! Read backed-up safe_mode for restart. If the file is not present
  ! safe_mode is set to .false.
  !--------------------------------------------------------------
  integer::ilevel,ilun,idummy,nlevelmax2,ierr,mg_safe_mode_reset_old,reset_diff
  logical::ok
  character(LEN=80)::fileloc,line
  character(LEN=5)::nchar
#ifndef WITHOUTMPI
  integer::info
#endif

  call title(nrestart,nchar)
  fileloc='output_'//TRIM(nchar)//'/safe_mode_'//TRIM(nchar)//'.txt'

  if(myid==1)then
     inquire(file=TRIM(fileloc),exist=ok)
     if(ok)then
        ilun=ncpu+myid+10
        open(unit=ilun,file=TRIM(fileloc),form='formatted')
        read(ilun,'(A)')line
        read(ilun,'(A8,i5)')line,nlevelmax2
        read(ilun,'(A18,1X,i8)')line,mg_safe_mode_reset_old
        read(ilun,'(A)')line
        do ilevel=1,nlevelmax2
           if(ilevel>nlevelmax)exit
           read(ilun,'(i5,L9,i12)',iostat=ierr)idummy,safe_mode(ilevel),safe_mode_countdown(ilevel)
           if(ierr/=0)then
              write(*,*)'Truncated safe mode file, ignoring levels from ',ilevel
              safe_mode(ilevel:nlevelmax)=.false.
              safe_mode_countdown(ilevel:nlevelmax)=0
              exit
           end if
        end do
        close(ilun)

        ! Verify the restart state is not in contradition with the current namelist
        if(mg_safe_mode_reset>0 .and. mg_safe_mode_reset.ne.mg_safe_mode_reset_old)then
           reset_diff = mg_safe_mode_reset - mg_safe_mode_reset_old !positive or negative
           do ilevel=1,nlevelmax
              if(safe_mode(ilevel))then
                 safe_mode_countdown(ilevel) = max(safe_mode_countdown(ilevel)+reset_diff, 0)
                 if(safe_mode_countdown(ilevel)==0)safe_mode(ilevel)=.false.
              end if
           end do
        endif

     else
        write(*,*)'No safe_mode file found in the restart output.'
        write(*,*)'All levels restart with safe_mode=.false.'
     end if
  end if


#ifndef WITHOUTMPI
  call MPI_BCAST(safe_mode          ,nlevelmax,MPI_LOGICAL,0,MPI_COMM_WORLD,info)
  call MPI_BCAST(safe_mode_countdown,nlevelmax,MPI_INTEGER,0,MPI_COMM_WORLD,info)
#endif

  if(verbose)write(*,*)'SAFE MODE backup file read completed'

end subroutine restore_safe_mode
