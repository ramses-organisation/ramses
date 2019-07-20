subroutine init_hydro
  use amr_commons
  use hydro_commons
  implicit none
  integer::ncell,ncache,iskip,igrid,i,ilevel,ind,ivar
  integer::nvar2,ilevel2,numbl2,ilun
  integer ,dimension(:),allocatable::ind_grid
  real(dp),dimension(:),allocatable::xx
  character(LEN=80)::fileloc
  character(LEN=5)::nchar

  if(verbose)write(*,*)'Entering init_hydro'

  ncell=ncoarse+twotondim*ngridmax
  
  ! Allocate conservative, cell-centered variables arrays
  allocate(uold(1:ncell,1:nvar))
  allocate(unew(1:ncell,1:nvar))
  uold=0.0d0; unew=0.0d0 
  allocate(divu(1:ncell))
  divu=0.0d0  

  if(nrestart>0)then

     ilun=ncpu+myid+10
     call title(nrestart,nchar)
     fileloc='backup_'//TRIM(nchar)//'/hydro_'//TRIM(nchar)//'.bak'
     call title(myid,nchar)
     fileloc=TRIM(fileloc)//TRIM(nchar)
     open(unit=ilun,file=fileloc,form='unformatted')
     read(ilun)nvar2
     if(nvar2.ne.nvar)then
        write(*,*)'File hydro.tmp is not compatible'
        write(*,*)'Found   =',nvar2
        write(*,*)'Expected=',nvar
        call clean_stop
     end if
     do ilevel=1,nlevelmax
        read(ilun)ilevel2
        read(ilun)numbl2
        if(numbl2.ne.numbl(myid,ilevel))then
           write(*,*)'File hydro.tmp is not compatible for level',ilevel
           write(*,*)'Found   =',numbl2,ilevel2
           write(*,*)'Expected=',numbl(myid,ilevel)
        end if
        if(numbl(myid,ilevel)>0)then
           ncache=numbl(myid,ilevel)
           allocate(ind_grid(1:ncache))
           allocate(xx(1:ncache))
           ! Loop over level grids
           igrid=headl(myid,ilevel)
           do i=1,ncache
              ind_grid(i)=igrid
              igrid=next(igrid)
           end do
           ! Loop over cells
           do ind=1,twotondim
              iskip=ncoarse+(ind-1)*ngridmax
              ! Loop over conservative variables
              do ivar=1,nvar
                 read(ilun)xx
                 do i=1,ncache
                    uold(ind_grid(i)+iskip,ivar)=xx(i)
                 end do
              end do
           end do
           deallocate(ind_grid,xx)
        end if
     end do
     close(ilun)

     ! Update boundaries
     do ilevel=1,nlevelmax
        do ivar=1,nvar
           call make_virtual_fine_dp(uold(1,ivar),ilevel)
        end do
        if(simple_boundary)call make_boundary_hydro(ilevel)
     end do

  end if

end subroutine init_hydro




