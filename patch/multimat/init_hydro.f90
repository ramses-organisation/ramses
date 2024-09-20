subroutine init_hydro
  use amr_commons
  use hydro_commons
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::info,info2,dummy_io
#endif
  integer::ncell,ncache,iskip,igrid,i,ilevel,ind,ivar,idim,imat
  integer::nvar2,ilevel2,numbl2,ilun,ibound,istart
  integer::ncpu2,ndim2,nlevelmax2,nboundary2
  integer ,dimension(:),allocatable::ind_grid
  real(dp),dimension(:),allocatable::xx
  real(dp)::gamma2
  character(LEN=80)::fileloc
  character(LEN=5)::nchar,ncharcpu
  real(dp)::dtot,ekin,erad
  real(dp),dimension(1:nvector,1:nmat),save::ff,gg,kk_mat
  real(dp),dimension(1:nvector,1:npri),save::qq
  real(dp),dimension(1:nvector),save::ee,cc,kk_hat

  if(verbose)write(*,*)'Entering init_hydro'

  ncell=ncoarse+twotondim*ngridmax

  ! Allocate conservative, cell-centered variables arrays
  allocate(uold(1:ncell,1:nvar))
  allocate(unew(1:ncell,1:nvar))
  uold=0.0d0; unew=0.0d0
  allocate(divu(1:ncell))
  divu=0.0d0

  if(nrestart>0) then

     ilun=ncpu+myid+103
     call title(nrestart,nchar)

     fileloc='output_'//TRIM(nchar)//'/hydro_'//TRIM(nchar)//'.out'

     call title(myid,nchar)
     fileloc=TRIM(fileloc)//TRIM(nchar)

     open(unit=ilun,file=fileloc,form='unformatted')
     read(ilun)ncpu2
     read(ilun)nvar2
     if(strict_equilibrium>0)nvar2=nvar2-2
     read(ilun)ndim2
     read(ilun)nlevelmax2
     read(ilun)nboundary2
     read(ilun)gamma2
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
              ! Loop over level grids
              igrid=istart
              do i=1,ncache
                 ind_grid(i)=igrid
                 igrid=next(igrid)
              end do
              ! Loop over cells
              do ind=1,twotondim
                 iskip=ncoarse+(ind-1)*ngridmax
                 ! Read density and velocities --> density and momenta
                 do ivar=1,ndim+1
                    read(ilun)xx
                    if(ivar==1)then
                       do i=1,ncache
                          uold(ind_grid(i)+iskip,1)=xx(i)
                       end do
                    else if(ivar>=2.and.ivar<=ndim+1)then
                       do i=1,ncache
                          uold(ind_grid(i)+iskip,ivar)=xx(i)*max(uold(ind_grid(i)+iskip,1),smallr)
                       end do
                    endif
                 end do
                 ! Read thermal pressure
                 read(ilun)xx
                 do i=1,ncache
                    uold(ind_grid(i)+iskip,ndim+2)=xx(i)
                 end do
                 ! Read volume fraction
                 do imat=1,nmat
                    ivar=ndim+2+nener+imat
                    read(ilun)xx
                    do i=1,ncache
                       uold(ind_grid(i)+iskip,ivar)=xx(i)
                    end do
                 end do
                 ! Read physical densities
                 do imat=1,nmat
                    ivar=ndim+2+nener+nmat+imat
                    read(ilun)xx
                    do i=1,ncache
                       uold(ind_grid(i)+iskip,ivar)=xx(i)
                    end do
                 end do
                 ! Convert pressure to total energy
                 do i = 1, ncache
                    dtot=max(uold(ind_grid(i)+iskip,1),smallr)
                    do imat=1,nmat
                       ff(1,imat)=uold(ind_grid(i)+iskip,imat+nener+npri)
                       gg(1,imat)=uold(ind_grid(i)+iskip,imat+nener+npri+nmat)
                    end do
                    qq(1,1)=dtot
                    ekin=0.0
                    do idim=1,ndim
                       qq(1,idim+1)=uold(ind_grid(i)+iskip,idim+1)/dtot
                       ekin=ekin+0.5d0*qq(1,idim+1)**2
                    end do
                    erad=00
#if NENER > 0
                    do irad = 1,nener
                       erad=erad+uold(ind_grid(i)+iskip,ndim+2+irad)
                    end do
#endif
                    qq(1,npri)=uold(ind_grid(i)+iskip,npri) ! Pressure
                    call eosinv(ff,gg,qq,ee,cc,kk_mat,kk_hat,1)
                    uold(ind_grid(i)+iskip,npri)=ee(1)+dtot*ekin+erad ! Total energy
                 end do

                 ! Read equilibrium density and pressure profiles
                 if(strict_equilibrium>0)then
                    read(ilun)xx
                    do i=1,ncache
                       rho_eq(ind_grid(i)+iskip)=xx(i)
                    end do
                    read(ilun)xx
                    do i=1,ncache
                       p_eq(ind_grid(i)+iskip)=xx(i)
                    end do
                 endif

              end do
              deallocate(ind_grid,xx)
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
