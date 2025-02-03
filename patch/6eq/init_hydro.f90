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
  real(dp)::ekin,erad
  logical ::inv
  real(dp),dimension(1:nvector,1:nmat),save::ff,gg
  real(dp),dimension(1:nvector,1:npri),save::qq
  real(dp),dimension(1:nvector),save::dtot,gg_mat,ee,pp_mat,cc

  if(verbose)write(*,*)'Entering init_hydro'

  ncell=ncoarse+twotondim*ngridmax

  ! Allocate conservative, cell-centered variables arrays
  allocate(uold(1:ncell,1:nvar))
  allocate(unew(1:ncell,1:nvar))
  uold=0.0d0; unew=0.0d0
  allocate(divu(1:ncell))
  allocate(dive(1:ncell,1:nmat))
  divu=0.0d0
  dive=0.0d0

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
                 ! Read volume fractions
                 do imat=1,nmat
                    ivar=imat
                    read(ilun)xx
                    do i=1,ncache
                       uold(ind_grid(i)+iskip,ivar)=xx(i)
                    end do
                 end do
                 ! Read physical densities
                 do imat=1,nmat
                    ivar=nmat+imat
                    read(ilun)xx
                    do i=1,ncache
                       uold(ind_grid(i)+iskip,ivar)=xx(i)
                    end do
                 end do
                 ! Calculate total density
                 dtot(1:ncache) = 0.0
                 do imat=1,nmat
                  do i=1,ncache
                    dtot(i) = dtot(i) + uold(ind_grid(i)+iskip,nmat+imat)
                  end do
                 end do
                 ! Read velocities -->  momenta
                 do ivar=2*nmat+1,2*nmat+ndim
                    read(ilun)xx
                    do i=1,ncache
                      uold(ind_grid(i)+iskip,ivar)=xx(i)*dtot(i)
                    end do
                 end do
                 ! Read thermal pressures
                 read(ilun)xx
                 do imat=1,nmat
                   do i=1,ncache
                      uold(ind_grid(i)+iskip,2*nmat+ndim+imat)=xx(i)              ! Saving the pressure into the energy slots
                      qq(i,ndim+imat) = uold(ind_grid(i)+iskip,2*nmat+ndim+imat)
                   end do
                 end do
                 ! Convert pressure to total energy
                 inv = .true.
                 do imat=1,nmat
                   do i=1,ncache
                     gg_mat(i) = gg(i,imat)
                     pp_mat(i) = qq(i,ndim+imat)
                   end do
                   call eos(gg_mat,ee,pp_mat,cc,imat,inv,ncache)
                   do i = 1, ncache
                      ff(i,imat)    = uold(ind_grid(i)+iskip,imat)
                      gg(i,imat)    = uold(ind_grid(i)+iskip,imat+nmat)/ff(i,imat)
                      ekin=0.0
                      do idim=1,ndim
                         qq(i,idim) = uold(ind_grid(i)+iskip,2*nmat+idim)/dtot(i)
                         ekin       = ekin + 0.5d0*qq(i,idim)**2
                      end do
                      erad=0.0
#if NENER > 0
                      do irad = 1,nener
                         erad       = erad + uold(ind_grid(i)+iskip,3*nmat+ndim+irad)
                      end do
#endif
                      uold(ind_grid(i)+iskip,2*nmat+ndim+imat) = (ee(i) + gg(i,imat)*ekin + erad)*ff(i,imat) ! f_k.E_k
                   end do
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
