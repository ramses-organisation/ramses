subroutine backup_hydro(filename)
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'  
#endif

  character(LEN=80)::filename

  integer::i,ivar,ncache,ind,ilevel,igrid,iskip,ilun,istart,ibound,irad
  integer,allocatable,dimension(:)::ind_grid
  real(dp),allocatable,dimension(:)::xdp
  character(LEN=5)::nchar
  character(LEN=80)::fileloc
  integer,parameter::tag=1121
  integer::dummy_io,info2

  if(verbose)write(*,*)'Entering backup_hydro'

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
  write(ilun)nvar
  write(ilun)ndim
  write(ilun)nlevelmax
  write(ilun)nboundary
  write(ilun)gamma
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
              do ivar=1,ndim+1
                 if(ivar==1)then
                    ! Write density
                    do i=1,ncache
                       xdp(i)=uold(ind_grid(i)+iskip,1)
                    end do
                 else if(ivar>=2.and.ivar<=ndim+1)then
                    ! Write velocity field
                    do i=1,ncache
                       xdp(i)=uold(ind_grid(i)+iskip,ivar)/max(uold(ind_grid(i)+iskip,1),smallr)
                    end do
                 endif
                 write(ilun)xdp
              end do
#if NENER>0
              ! Write non-thermal pressures
              do ivar=ndim+3,ndim+2+nener
                 do i=1,ncache
                    xdp(i)=(gamma_rad(ivar-ndim-2)-1d0)*uold(ind_grid(i)+iskip,ivar)
                 end do
                 write(ilun)xdp
              end do
#endif
              ! Write thermal pressure
              do i=1,ncache
                 xdp(i)=uold(ind_grid(i)+iskip,ndim+2)
                 xdp(i)=xdp(i)-0.5d0*uold(ind_grid(i)+iskip,2)**2/max(uold(ind_grid(i)+iskip,1),smallr)
#if NDIM>1
                 xdp(i)=xdp(i)-0.5d0*uold(ind_grid(i)+iskip,3)**2/max(uold(ind_grid(i)+iskip,1),smallr)
#endif
#if NDIM>2
                 xdp(i)=xdp(i)-0.5d0*uold(ind_grid(i)+iskip,4)**2/max(uold(ind_grid(i)+iskip,1),smallr)
#endif
#if NENER>0
                 do irad=1,nener
                    xdp(i)=xdp(i)-uold(ind_grid(i)+iskip,ndim+2+irad)
                 end do
#endif
                 xdp(i)=(gamma-1d0)*xdp(i)
              end do
              write(ilun)xdp
#if NVAR>NDIM+2+NENER
              ! Write passive scalars
              do ivar=ndim+3+nener,nvar
                 do i=1,ncache
                    xdp(i)=uold(ind_grid(i)+iskip,ivar)/max(uold(ind_grid(i)+iskip,1),smallr)
                 end do
                 write(ilun)xdp
              end do
#endif
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
  
  
end subroutine backup_hydro





