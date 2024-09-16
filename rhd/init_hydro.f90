subroutine init_hydro
  use amr_commons
  use hydro_commons
  use mpi_mod
  implicit none
  integer::ncell,ncache,iskip,igrid,i,ilevel,ind,ivar
  integer::nvar2,ilevel2,numbl2,ilun,ibound,istart,info
  integer::ncpu2,ndim2,nlevelmax2,nboundary2
  integer ,dimension(:),allocatable::ind_grid
  real(dp),dimension(:),allocatable::rr,vx,vy,vz,pp,xx
  real(dp)::gamma2
  character(LEN=80)::fileloc
  character(LEN=5)::nchar
  real(dp)::lor,tau,h

  if(verbose)write(*,*)'Entering init_hydro'

  !------------------------------------------------------
  ! Allocate conservative, cell-centered variables arrays
  !------------------------------------------------------
  ncell=ncoarse+twotondim*ngridmax
  allocate(uold(1:ncell,1:nvar))
  allocate(unew(1:ncell,1:nvar))
  uold=0.0d0; unew=0.0d0
  if(pressure_fix)then
     allocate(divu(1:ncell))
     allocate(enew(1:ncell))
     divu=0.0d0; enew=0.0d0
  end if

  !--------------------------------
  ! For a restart, read hydro file
  !--------------------------------
  if(nrestart>0)then
     ilun=ncpu+myid+10
     call title(nrestart,nchar)
     fileloc='output_'//TRIM(nchar)//'/hydro_'//TRIM(nchar)//'.out'
     call title(myid,nchar)
     fileloc=TRIM(fileloc)//TRIM(nchar)
     open(unit=ilun,file=fileloc,form='unformatted')
     read(ilun)ncpu2
     read(ilun)nvar2
     read(ilun)ndim2
     read(ilun)nlevelmax2
     read(ilun)nboundary2
     read(ilun)gamma2
     if(nvar2.ne.nvar)then
        write(*,*)'File hydro.tmp is not compatible'
        write(*,*)'Found   =',nvar2
        write(*,*)'Expected=',nvar
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
              write(*,*)'File hydro.tmp is not compatible'
              write(*,*)'Found   =',numbl2,' for level ',ilevel2
              write(*,*)'Expected=',ncache,' for level ',ilevel
           end if
           if(ncache>0)then
              allocate(ind_grid(1:ncache))
              allocate(xx(1:ncache))
              allocate(rr(1:ncache))
              allocate(vx(1:ncache))
              allocate(vy(1:ncache))
              allocate(vz(1:ncache))
              allocate(pp(1:ncache))
              ! Loop over level grids
              igrid=istart
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
                    if(ivar==1)then
                       rr=xx
                    else if(ivar==2)then
                       vx=xx
                    else if(ivar==3)then
                       vy=xx
                    else if(ivar==4)then
                       vz=xx
                    else if(ivar==5)then
                       pp=xx
                    else
                       !passive scalars
                       do i=1,ncache
                          uold(ind_grid(i)+iskip,ivar)=xx(i)!*uold(ind_grid(i)+iskip,1)
                       enddo
                    endif
                 enddo
                 do i=1,ncache
                    ! convert to conservative variables
                    lor=(1-vx(i)**2-vy(i)**2-vz(i)**2)**(-1./2.)
                    tau=pp(i)/rr(i)
                    if (eos .eq. 'TM') then
                       h=5d0/2d0*tau+3d0/2d0*sqrt(tau**2+4d0/9d0)
                    else
                       h=1d0+tau*gamma/(gamma-1d0)
                    endif
                    uold(ind_grid(i)+iskip,1)=rr(i)*lor
                    uold(ind_grid(i)+iskip,2)=rr(i)*h*lor**2*vx(i)
                    uold(ind_grid(i)+iskip,3)=rr(i)*h*lor**2*vy(i)
                    uold(ind_grid(i)+iskip,4)=rr(i)*h*lor**2*vz(i)
                    uold(ind_grid(i)+iskip,5)=rr(i)*h*lor**2-pp(i)

                    do ivar=6,nvar
                       uold(ind_grid(i)+iskip,ivar)= uold(ind_grid(i)+iskip,ivar)*uold(ind_grid(i)+iskip,1)*lor
                    enddo
                 end do
              end do
              deallocate(ind_grid,xx,rr,vx,vy,vz,pp)
           end if
        end do
     end do
     close(ilun)
#ifndef WITHOUTMPI
     if(debug)write(*,*)'hydro.tmp read for processor ',myid
     call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
     if(verbose)write(*,*)'HYDRO backup files read completed'
  end if

end subroutine init_hydro
