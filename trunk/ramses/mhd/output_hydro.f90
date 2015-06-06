subroutine backup_hydro(filename)
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'  
#endif
  character(LEN=80)::filename

  integer::i,ivar,ncache,ind,ilevel,igrid,iskip,ilun,istart,ibound,irad
  real(dp)::d,u,v,w,A,B,C,e
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
  write(ilun)nvar+3
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
              do ivar=1,4
                 if(ivar==1)then ! Write density
                    do i=1,ncache
                       xdp(i)=uold(ind_grid(i)+iskip,1)
                    end do
                 else ! Write velocity field
                    do i=1,ncache
                       xdp(i)=uold(ind_grid(i)+iskip,ivar)/max(uold(ind_grid(i)+iskip,1),smallr)
                    end do
                 endif
                 write(ilun)xdp
              end do
              do ivar=6,8 ! Write left B field
                 do i=1,ncache
                    xdp(i)=uold(ind_grid(i)+iskip,ivar)
                 end do
                 write(ilun)xdp
              end do
              do ivar=nvar+1,nvar+3 ! Write right B field
                 do i=1,ncache
                    xdp(i)=uold(ind_grid(i)+iskip,ivar)
                 end do
                 write(ilun)xdp
              end do
#if NENER>0
              ! Write non-thermal pressures
              do ivar=9,8+nener
                 do i=1,ncache
                    xdp(i)=(gamma_rad(ivar-8)-1d0)*uold(ind_grid(i)+iskip,ivar)
                 end do
                 write(ilun)xdp
              end do
#endif
              do i=1,ncache ! Write thermal pressure
                 d=max(uold(ind_grid(i)+iskip,1),smallr)
                 u=uold(ind_grid(i)+iskip,2)/d
                 v=uold(ind_grid(i)+iskip,3)/d
                 w=uold(ind_grid(i)+iskip,4)/d
                 A=0.5*(uold(ind_grid(i)+iskip,6)+uold(ind_grid(i)+iskip,nvar+1))
                 B=0.5*(uold(ind_grid(i)+iskip,7)+uold(ind_grid(i)+iskip,nvar+2))
                 C=0.5*(uold(ind_grid(i)+iskip,8)+uold(ind_grid(i)+iskip,nvar+3))
                 e=uold(ind_grid(i)+iskip,5)-0.5*d*(u**2+v**2+w**2)-0.5*(A**2+B**2+C**2)
#if NENER>0
                 do irad=1,nener
                    e=e-uold(ind_grid(i)+iskip,8+irad)
                 end do
#endif
                 xdp(i)=(gamma-1d0)*e
              end do
              write(ilun)xdp
#if NVAR > 8+NENER
              do ivar=9+nener,nvar ! Write passive scalars if any
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





