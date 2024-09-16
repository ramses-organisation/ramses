subroutine backup_hydro(filename)
  use amr_commons
  use hydro_commons
  implicit none
  character(LEN=80)::filename

  integer::i,ivar,ncache,ind,ilevel,igrid,iskip,ilun,istart,ibound
  integer,allocatable,dimension(:)::ind_grid
  real(dp),allocatable,dimension(:)::xdp
  real(dp),allocatable,dimension(:,:)::qq ! primitive variables
  character(LEN=5)::nchar
  character(LEN=80)::fileloc
  integer::k,ncell
  real(dp) :: lor,entho ! Lorentz factor
  real(dp) :: D,M,E,Mx,My,Mz,u2,Xsi,R
  real(dp) ::rho,p,vpar,vx,vy,vz,smallp,tau
  real(dp) :: small_bigD=1e-12

  if(verbose)write(*,*)'Entering backup_hydro'
  ilun=ncpu+myid+10

  call title(myid,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)
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
           allocate(qq(1:ncache,1:nvar))
           ! Loop over level grids
           igrid=istart
           do i=1,ncache
              ind_grid(i)=igrid
              igrid=next(igrid)
           end do
           ! Loop over cells
           do ind=1,twotondim
              iskip=ncoarse+(ind-1)*ngridmax

              do i=1,ncache

              !convert to primitive variables
              ! Compute density
              D = uold(ind_grid(i)+iskip,1)
              ! Compute momentum
              Mx=uold(ind_grid(i)+iskip,2)
              My=uold(ind_grid(i)+iskip,3)
              Mz=uold(ind_grid(i)+iskip,4)
              M = sqrt(Mx**2+My**2+Mz**2)
              ! Compute total energy
              E = uold(ind_grid(i)+iskip,5)
              !Method from Mignone,McKinney,2007. Same as BS2011 except one uses E'=U-D and u^2=Lor^2*v^2

              if (D<0) then
                 write(*,*),'D<0 in ctoprimbis'
                 D=small_bigD
                 uold(ind_grid(i)+iskip,1)=D
              endif
              if (E<0) then
                 E=sqrt(M**2+d**2+1d-8)
                 write(*,*),'E<0 in ctoprimbis'
                 uold(ind_grid(i)+iskip,5)=E
              endif

              if (E**2<M**2+D**2) then

                 write (*,*) 'Switch...ctoprimbis'
                 E=sqrt(M**2+d**2+1d-8)
                 uold(ind_grid(i)+iskip,5)=E
              endif


              if ( M .eq. 0) then
                 qq(i,1) = D
                 qq(i,2) = 0d0
                 qq(i,3) = 0d0
                 qq(i,4) = 0d0
                 if (eos .eq. 'TM') then
                    qq(i,5) =(E**2-D**2)/3d0/E
                 else
                    qq(i,5)=(E-D)*(gamma-1d0)
                 endif
                 lor=1d0
              else

                 call Newton_Raphson_Mignone(D,M,E,gamma,R)

                 ! Compute the Lorentz factor
                 u2  = M**2.0d0/(R**2.0d0-M**2.0d0)
                 lor = (1.0d0+u2)**(1d0/2d0)

                 ! Compute the density
                 qq(i,1) = D/lor

                 ! compute velocities
                 qq(i,2) = Mx/R
                 qq(i,3) = My/R
                 qq(i,4) = Mz/R

                 ! Compute pressure
                 Xsi=((R-D)-u2/(lor+1d0)*D)/lor**2

                 if (eos .eq. 'TM') then
                    rho=qq(i,1)
                    qq(i,5)=(2d0*xsi*(xsi+2d0*rho))/(5d0*(xsi+rho)+sqrt(9d0*xsi**2+18d0*rho*xsi+25d0*rho**2))
                 else
                    qq(i,5)=(gamma-1d0)/gamma*Xsi
                 endif
              endif
              if ((qq(i,1)<0d0) .or.((qq(i,5)<0d0))) then

                 write(*,*) 'negative pressure or density output'
              endif

              do ivar=6,nvar
                 qq(i,ivar)=uold(ind_grid(i)+iskip,ivar)/uold(ind_grid(i)+iskip,1)
              enddo
           enddo
           do ivar=1,nvar
              if(ivar==1)then ! Write density
                 do i=1,ncache
                    xdp(i)=qq(i,1)
                 end do
              else if(ivar>=2.and.ivar<=4)then ! Write velocity field
                 do i=1,ncache
                    xdp(i)=qq(i,ivar)
                 end do
              else if(ivar==5)then ! Write pressure
                 do i=1,ncache
                    xdp(i)=qq(i,ivar)
                 end do
              else ! Write passive scalars if any
                 do i=1,ncache
                    xdp(i)=qq(i,ivar)
                 end do
              endif
              write(ilun)xdp
           end do
        end do
        deallocate(ind_grid, xdp,qq)
     end if
  end do
end do
close(ilun)

end subroutine backup_hydro
