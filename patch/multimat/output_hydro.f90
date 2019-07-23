subroutine output_hydro(filename)
  use amr_commons
  use hydro_commons
  implicit none
  character(LEN=80)::filename

  integer::i,ivar,idim,imat,ncache,ind,ilevel,igrid,iskip,ilun
  integer,allocatable,dimension(:)::ind_grid
  real(sp),allocatable,dimension(:)::xsp
  real(sp)::gamma_sp
  real(dp)::dtot,ekin
  real(dp),dimension(1:nvector,1:nmat),save::ff,gg
  real(dp),dimension(1:nvector,1:npri),save::qq
  real(dp),dimension(1:nvector),save::pp,cc
  character(LEN=5)::nchar
  character(LEN=80)::fileloc

  if(verbose)write(*,*)'Entering output_hydro'

  ilun=ncpu+myid+10

  gamma_sp=gamma
  call title(myid,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)
  open(unit=ilun,file=fileloc,form='unformatted')
  write(ilun)ncpu
  write(ilun)nvar
  write(ilun)ndim
  write(ilun)nlevelmax
  write(ilun)gamma_sp
  ! Output primitive variables
  do ilevel=1,nlevelmax
     write(ilun)ilevel
     write(ilun)numbl(myid,ilevel)
     if(numbl(myid,ilevel)>0)then
        ncache=numbl(myid,ilevel)
        allocate(ind_grid(1:ncache),xsp(1:ncache))
        ! Loop over level grids
        igrid=headl(myid,ilevel)
        do i=1,ncache
           ind_grid(i)=igrid
           igrid=next(igrid)
        end do
        ! Loop over cells
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ncache
              ind_grid(i)=ind_grid(i)+iskip
           end do
           ! Output total density
           do i=1,ncache
              xsp(i)=uold(ind_grid(i),1)
           end do
           write(ilun)xsp
           ! Output velocity (in comoving km s-1)
           do idim=1,ndim
              do i=1,ncache
                 dtot=uold(ind_grid(i),1)
                 xsp(i)=uold(ind_grid(i),idim+1)/dtot
              end do
              write(ilun)xsp
           end do
           ! Output pressure (in consistant units)
           do i=1,ncache
              dtot=uold(ind_grid(i),1)
              do imat=1,nmat
                 ff(1,imat)=uold(ind_grid(i),imat+npri)
                 gg(1,imat)=uold(ind_grid(i),imat+npri+nmat)
              end do
              qq(1,1)=dtot
              ekin=0.0
              do idim=1,ndim
                 qq(1,idim+1)=uold(ind_grid(i),idim+1)/dtot
                 ekin=ekin+0.5d0*qq(1,idim+1)**2
              end do
              qq(1,npri)=uold(ind_grid(i),npri)-dtot*ekin
              call eos(ff,gg,qq,pp,cc,1)
              xsp(i)=pp(1)       ! Pressure
           end do
           write(ilun)xsp
           ! Output volume fraction
           do imat=1,nmat
              do i=1,ncache
                 xsp(i)=uold(ind_grid(i),imat+npri)
              end do
              write(ilun)xsp
           end do
           ! Output fluid density
           do imat=1,nmat
              do i=1,ncache
                 xsp(i)=uold(ind_grid(i),imat+npri+nmat)
              end do
              write(ilun)xsp
           end do
           do i=1,ncache
              ind_grid(i)=ind_grid(i)-iskip
           end do
        end do
        deallocate(ind_grid,xsp)
     end if
  end do
  close(ilun)

end subroutine output_hydro

subroutine backup_hydro(filename)
  use amr_commons
  use hydro_commons
  implicit none
  character(LEN=80)::filename

  integer::i,ivar,ncache,ind,ilevel,igrid,iskip,ilun
  integer,allocatable,dimension(:)::ind_grid
  real(dp),allocatable,dimension(:)::xdp
  real(sp),allocatable,dimension(:)::xsp
  character(LEN=5)::nchar
  character(LEN=80)::fileloc

  if(verbose)write(*,*)'Entering backup_hydro'

  ilun=ncpu+myid+10
     
  call title(myid,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)
  open(unit=ilun,file=fileloc,form='unformatted')
  write(ilun)nvar
  do ilevel=1,nlevelmax
     write(ilun)ilevel
     write(ilun)numbl(myid,ilevel)
     if(numbl(myid,ilevel)>0)then
        ncache=numbl(myid,ilevel)
        allocate(ind_grid(1:ncache),xdp(1:ncache))
        ! Loop over level grids
        igrid=headl(myid,ilevel)
        do i=1,ncache
           ind_grid(i)=igrid
           igrid=next(igrid)
        end do
        ! Loop over cells
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ncache
              ind_grid(i)=ind_grid(i)+iskip
           end do
           do ivar=1,nvar
              do i=1,ncache
                 xdp(i)=uold(ind_grid(i),ivar)
              end do
              write(ilun)xdp
           end do
           do i=1,ncache
              ind_grid(i)=ind_grid(i)-iskip
           end do
        end do
        deallocate(ind_grid, xdp)
     end if
  end do
  close(ilun)
     
end subroutine backup_hydro





