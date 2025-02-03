#if NDIM==3
!################################################################
!################################################################
!################################################################
subroutine write_clump_field
  use amr_commons
  use clfind_commons
  use mpi_mod
  implicit none

  !---------------------------------------------------------------------------
  ! This routine writes the clump field to disk as integer AMR field. This AMR
  ! field can be read by the usual postprocessing tools.
  !---------------------------------------------------------------------------
  character(LEN=80)::filename
  integer::i,ncache,ind,ilevel,igrid,iskip,ilun,istart,ibound
  integer, allocatable, dimension(:) :: ind_grid, id_clump
  character(LEN=5)::myidstring,nchar,ncharcpu
  integer,parameter::tag=1102
  integer::dummy_io,info2

  if(verbose)write(*,*)'Entering write_clump_field'

  ilun=ncpu+myid+10

  ! Wait for the token
#ifndef WITHOUTMPI
  if(IOGROUPSIZE>0) then
     if (mod(myid-1,IOGROUPSIZE)/=0) then
        call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
             & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
     end if
  endif
#endif

  call title(ifout,nchar)
  call title(myid,myidstring)

  ! write file with description of the fields
  if (myid==1)then
     filename = 'output_'//TRIM(nchar)//'/clump_field_file_descriptor.txt'
     call file_descriptor_clump(filename)
  end if

  if(IOGROUPSIZEREP>0)then
     call title(((myid-1)/IOGROUPSIZEREP)+1,ncharcpu)
     open(unit=ilun, file=TRIM('output_'//TRIM(nchar)//'/group_'//TRIM(ncharcpu)//'/clump_field_'//TRIM(nchar)//'.out'//TRIM(myidstring)),form='unformatted')
  else
     open(unit=ilun, file=TRIM('output_'//TRIM(nchar)//'/clump_field_'//TRIM(nchar)//'.out'//TRIM(myidstring)),form='unformatted')
  endif

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
           allocate(ind_grid(1:ncache),id_clump(1:ncache))
           ! Loop over level grids
           igrid=istart
           do i=1,ncache
              ind_grid(i)=igrid
              igrid=next(igrid)
           end do
           ! Loop over cells
           do ind=1,twotondim
              iskip=ncoarse+(ind-1)*ngridmax
              ! Write clump index
              do i=1,ncache
                 id_clump(i)=flag2(ind_grid(i)+iskip)
              end do
              write(ilun)id_clump
           end do
           deallocate(ind_grid, id_clump)
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

end subroutine write_clump_field
!################################################################
!################################################################
!################################################################
subroutine file_descriptor_clump(filename)
  use amr_commons
  use hydro_commons
  implicit none

  ! Pretty much a dummy file at the moment since the clump file contains
  ! only one AMR field. However, this might change at some point...

  character(LEN=80)::filename
  character(LEN=80)::fileloc
  integer::ivar,ilun

  if(verbose)write(*,*)'Entering file_descriptor_clump'
  ilun=11

  ! Open file
  fileloc=TRIM(filename)
  open(unit=ilun,file=fileloc,form='formatted')

  ! Write run parameters
  write(ilun,'("nvar        = 1")')
  ivar=1
  write(ilun,'("variable #",I2,": Level 0 clump id")')ivar
  close(ilun)

end subroutine file_descriptor_clump
!################################################################
!################################################################
!################################################################
#endif
