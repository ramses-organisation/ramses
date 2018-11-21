!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine upload_fine(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !----------------------------------------------------------------------
  ! This routine performs a restriction operation (averaging down)
  ! for the hydro variables.
  !----------------------------------------------------------------------
  integer::i,ncache,igrid,ngrid,ind,iskip,nsplit,icell
  integer,dimension(1:nvector),save::ind_grid,ind_cell,ind_split
  logical,dimension(1:nvector),save::ok

  if(ilevel==nlevelmax)return
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Loop over active grids by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do

     ! Loop over cells
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do

        ! Gather split cells
        do i=1,ngrid
           ok(i)=son(ind_cell(i))>0
        end do

        ! Count split cells
        nsplit=0
        do i=1,ngrid
           if(ok(i))nsplit=nsplit+1
        end do

        ! Upload for selected cells
        if(nsplit>0)then
           icell=0
           do i=1,ngrid
              if(ok(i))then
                 icell=icell+1
                 ind_split(icell)=ind_cell(i)
              end if
           end do
           call upl(ind_split,nsplit)
        end if

     end do
     ! End loop over cells

  end do
  ! End loop over grids

111 format('   Entering upload_fine for level',i2)

end subroutine upload_fine
!##########################################################################
!##########################################################################
!##########################################################################
!##########################################################################
subroutine upl(ind_cell,ncell)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ncell
  integer,dimension(1:nvector)::ind_cell
  !---------------------------------------------------------------------
  ! This routine performs a restriction operation (averaging down)
  ! for the following variables:
  ! interpol_var=0: use D, Mx,My,Mz and E
  ! interpol_var=1: use D, Mx, My,Mz  and  epsilon
  !---------------------------------------------------------------------
  integer ::ivar,i,idim,ind_son,iskip_son
  integer ,dimension(1:nvector),save::igrid_son,ind_cell_son
  real(dp),dimension(1:nvector),save::getx,eint
  real(dp)::D,Mx,My,Mz,M,E,M2, lor, h,xsi,qd,qp,qx,qy,qz,r,v2
  real(dp) :: small_bigD=1e-12,tau

  ! Get child oct index
  do i=1,ncell
     igrid_son(i)=son(ind_cell(i))
  end do

  if(interpol_var==1)then

     do ivar=1,nvar
        if (ivar .ne. 5) then
           ! Average conservative variable
           getx(1:ncell)=0.0d0
           do ind_son=1,twotondim
              iskip_son=ncoarse+(ind_son-1)*ngridmax
              do i=1,ncell
                 ind_cell_son(i)=iskip_son+igrid_son(i)
              end do
              do i=1,ncell
                 getx(i)=getx(i)+uold(ind_cell_son(i),ivar)
              end do
           end do
           ! Scatter result to cells
           do i=1,ncell
              uold(ind_cell(i),ivar)=getx(i)/dble(twotondim)
           end do

        else
           ! Average internal energy
           getx(1:ncell)=0.0d0
           do ind_son=1,twotondim
              iskip_son=ncoarse+(ind_son-1)*ngridmax
              do i=1,ncell
                 ind_cell_son(i)=iskip_son+igrid_son(i)
              end do
              ! to compute the child internal energy epsilon=P/rho/(gamma-1), we need the P and rho
              do i=1,ncell
                 D = uold(ind_cell_son(i),1)
                 ! Compute momentum
                 Mx=uold(ind_cell_son(i),2) ; My=uold(ind_cell_son(i),3) ; Mz=uold(ind_cell_son(i),4)
                 M = sqrt(Mx**2+My**2+Mz**2)
                 ! Compute total energy
                 E = uold(ind_cell_son(i),5)
                 !Method from Mignone,McKinney,2007. Same as BS2011 except one uses E'=U-D and u^2=Lor^2*v^2

                 if (D<0) then
                    write(*,*),'D<0 switch'
                    D=small_bigD
                    uold(ind_cell_son(i),1)=D
                 endif
                 if (E<0) then
                    write(*,*),'E<0 switch'
                    E=sqrt((M**2+d**2)*(1+1d-8))
                    uold(ind_cell_son(i),5)=E
                 endif

                 if (E**2<M**2+D**2) then
                    write(*,*),'E<M switch'
                    E=sqrt((M**2+d**2)*(1+1d-8))
                    uold(ind_cell_son(i),5)=E
                 endif


                 call Newton_Raphson_Mignone(D,M,E,gamma,R)

                 ! Compute the Lorentz factor
                 v2  = M**2.0d0/(R**2.0d0-M**2.0d0)
                 lor = (1.0d0+v2)**(1d0/2d0)

                 ! Compute the density
                 qd = D/lor

                 ! Compute pressure
                 Xsi=((R-D)-v2/(lor+1d0)*D)/lor**2

                 if (eos .eq. 'TM') then
                    qp=(2d0*xsi*(xsi+2d0*qd))/(5d0*(xsi+qd)+sqrt(9d0*xsi**2+18d0*qd*xsi+25d0*qd**2))
                    ! Compute child internal energy (eint)
                    tau=qp/qd
                    eint(i)=3d0/2d0*tau+3d0/2d0*(tau**2+4d0/9d0)**(1d0/2d0)-1d0
                 else
                    qp=(gamma-1d0)/gamma*Xsi
                 ! Compute child internal energy (eint)
                    eint(i)= qp/qd/(gamma-1d0)
                 endif
                 if ((qd<0d0).or.(qp<0d0).or.E<0d0) then
                    write(*,*) 'negative pressure or density interpol hydro'
!                    write(*,*),qp,qd,D,M,E
                    stop
                 endif



              end do
              ! update the average for  eint
              do i=1,ncell
                 getx(i)=getx(i)+eint(i)
              end do
           end do

           ! compute the updated total energy for father cells
           do i=1,ncell
              getx(i)=getx(i)/dble(twotondim)
              if (eos .eq. 'TM') then
                 tau=getx(i)*(getx(i)+2d0)/(3d0*(getx(i)+1d0))
              else
                 tau=getx(i)*(gamma-1d0)
              endif
              ! for any EOS
              h=1d0+getx(i)+tau
              m2=uold(ind_cell(i),2)**2+uold(ind_cell(i),3)**2+uold(ind_cell(i),4)**2
              lor=sqrt(m2/uold(ind_cell(i),1)**2/h**2+1.0d0)
              uold(ind_cell(i),5)=uold(ind_cell(i),1)*h*lor-(h-1d0-getx(i))*uold(ind_cell(i),1)/lor
           enddo
        endif
     enddo
  else

     do ivar=1,nvar
        ! Average conservative variable
        getx(1:ncell)=0.0d0
        do ind_son=1,twotondim
           iskip_son=ncoarse+(ind_son-1)*ngridmax
           do i=1,ncell
              ind_cell_son(i)=iskip_son+igrid_son(i)
           end do
           do i=1,ncell
              getx(i)=getx(i)+uold(ind_cell_son(i),ivar)
           end do
        end do

     ! Scatter result to cells
        do i=1,ncell
           uold(ind_cell(i),ivar)=getx(i)/dble(twotondim)
        end do
     enddo

  endif



end subroutine upl
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine interpol_hydro(u1,g1,u2,g2,nn)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::nn
  real(dp),dimension(1:nvector,0:twondim  ,1:nvar)::u1
  real(dp),dimension(1:nvector,0:twondim  ,1:ndim)::g1
  real(dp),dimension(1:nvector,1:twotondim,1:nvar)::u2
  real(dp),dimension(1:nvector,1:twotondim,1:ndim)::g2
  !----------------------------------------------------------
  ! This routine performs a prolongation (interpolation)
  ! operation for newly refined cells or buffer cells.
  ! The interpolated variables are:
  ! interpol_var=0: D, M and E
  ! interpol_var=1: D,  M and epsilon
  ! The interpolation method is:
  ! interpol_type=0 straight injection
  ! interpol_type=1 linear interpolation with MinMod slope
  ! interpol_type=2 linear interpolation with Monotonized Central slope
  ! interpol_type=3 linear interpolation with unlimited Central slope
  ! The gravitational acceleration is also prolongated
  ! using straight injection only.
  !----------------------------------------------------------
  integer::i,j,ivar,idim,ind,ix,iy,iz,compteur

  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,0:twondim),save::a
  real(dp),dimension(1:nvector,1:ndim),save::w
  real(dp),dimension(1:nvector),save::ekin

  real(dp) :: lor ! Lorentz factor
  real(dp) :: D,M,Mx,My,Mz,v2,Xsi,R,vx,vy,vz
  real(dp) :: Dr,Mr,E,Mxr,Myr,Mzr
  real(dp),dimension(1:nvector,0:twotondim,1:nvar)::u_record
  real(dp),dimension(1:nvector,0:twotondim) ::Erecord ! total energy
  real(dp) ::qp,qd,qx,qy,qz,h,m2,tau
  real(dp) :: small_bigD=1e-12



  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)
  end do

  ! If necessary, convert father total energy into internal energy
  if(interpol_var==1)then
     do j=0,twondim
        do i=1,nn
           ! eint=p/rho/(gamma-1)
           D = u1(i,j,1)
           ! Compute momentum
           Mx=u1(i,j,2) ; My=u1(i,j,3) ; Mz=u1(i,j,4)
           M = sqrt(Mx**2+My**2+Mz**2)
           ! Compute total energy
           E = u1(i,j,5)
           !Method from Mignone,McKinney,2007. Same as BS2011 except one uses E'=U-D and u^2=Lor^2*v^2

           if (D<0) then
              write(*,*),'D<0 switch'
              D=small_bigD
              u1(i,j,1)=D
           endif
           if (E<0) then
              write(*,*),'E<0 switch'
              E=sqrt(M**2+d**2+1d-8)
              u1(i,j,5)=E
           endif
           if (E**2<M**2+D**2) then
              write(*,*),'E<M switch'
              E=sqrt(M**2+d**2+1d-8)
              u1(i,j,5)=E
           endif


           call Newton_Raphson_Mignone(D,M,E,gamma,R)

           ! Compute the Lorentz factor
           v2  = M**2.0d0/(R**2.0d0-M**2.0d0)
           lor = (1.0d0+v2)**(1d0/2d0)

           ! Compute the density
           qd = D/lor

           ! Compute pressure
           Xsi=((R-D)-v2/(lor+1d0)*D)/lor**2

           if (eos .eq. 'TM') then
              qp=(2d0*xsi*(xsi+2d0*qd))/(5d0*(xsi+qd)+sqrt(9d0*xsi**2+18d0*qd*xsi+25d0*qd**2))
              ! Compute child internal energy (eint)
              tau=qp/qd
              u1(i,j,5)=3d0/2d0*tau+3d0/2d0*(tau**2+4d0/9d0)**(1d0/2d0)-1d0
           else
              qp=(gamma-1d0)/gamma*Xsi
              ! Compute child internal energy (eint)
              u1(i,j,5)= qp/qd/(gamma-1d0)
           endif

           if ((qd<0d0).or.(qp<0d0).or. (E<0d0)) then
              write(*,*) 'negative pressure or density interpol hydro'
              stop
           endif

           !store total energie
           Erecord(i,j)=E
        end do
     end do
  end if

  ! Loop over interpolation variables
  do ivar=1,nvar

     ! Load father variable
     do j=0,twondim
        do i=1,nn
           a(i,j)=u1(i,j,ivar)
        end do
     end do


     ! Reset gradient
     w(1:nn,1:ndim)=0.0D0

     ! Compute gradient with chosen limiter
     if(interpol_type==1)call compute_limiter_minmod(a,w,nn)
     if(interpol_type==2)call compute_limiter_central(a,w,nn)
     if(interpol_type==3)call compute_central(a,w,nn)

     ! Interpolate over children cells
     do ind=1,twotondim
        u2(1:nn,ind,ivar)=a(1:nn,0)
        u_record(1:nn,ind,ivar)=a(1:nn,0)
        if (interpol_var == 1) then
           u_record(1:nn,ind,5)=Erecord(1:nn,0)
        endif
        do idim=1,ndim
           do i=1,nn
              u2(i,ind,ivar)=u2(i,ind,ivar)+w(i,idim)*xc(ind,idim)
           end do
        end do
     end do

  end do
  ! End loop over variables

  do j=0,twondim
     do i=1,nn
        if (u1(i,j,5) <0 )then
           write(*,*),'oops pere'
        endif
     end do
  end do


  ! If necessary, convert children internal energy into total energy
  if(interpol_var==1)then
     do ind=1,twotondim
        do i=1,nn
           if (eos .eq. 'TM') then
              tau=u2(i,ind,5)*(u2(i,ind,5)+2d0)/(3d0*u2(i,ind,5)+1d0)
           else
              tau=u2(i,ind,5)*(gamma-1d0)
              !              h=1.0d0+gamma*getx(i)/dble(twotondim)
           endif

           h=1+tau+u2(i,ind,5)
           m2=u2(i,ind,2)**2+u2(i,ind,3)**2+u2(i,ind,4)**2
           lor=sqrt(m2/u2(i,ind,1)**2/h**2+1.0d0)
           u2(i,ind,5)=u2(i,ind,1)*h*lor-(h-1d0-u2(i,ind,5))*u2(i,ind,1)/lor
        end do
     end do
  end if


! check if we get a physical state, otherwise, switch to first order
  compteur=0
  do ind=1,twotondim
     do idim=1,ndim
        do i=1,nn
           D = u2(i,ind,1)
           Mx=u2(i,ind,2) ; My=u2(i,ind,3) ; Mz=u2(i,ind,4)
           M = sqrt(Mx**2+My**2+Mz**2)
           E = u2(i,ind,5)
           if (M>E .or. (E**2<M**2+D**2).or.(E<0)) then
              compteur=compteur+1
           endif

       enddo
    enddo
 enddo
 if (compteur .ne. 0) then
    do ind=1,twotondim
       do idim=1,ndim
          do i=1,nn
             write(*,*),'back to first order'
             do ivar=1,nvar
                u2(i,ind,ivar)=u_record(i,ind,ivar)
             enddo
          end do
       end do
    end do
 endif




end subroutine interpol_hydro
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine compute_limiter_minmod(a,w,nn)
  use amr_commons
  use hydro_commons
  implicit none
  integer::nn
  real(dp),dimension(1:nvector,0:twondim)::a
  real(dp),dimension(1:nvector,1:ndim)::w
  !---------------
  ! MinMod slope
  !---------------
  integer::i,idim
  real(dp)::diff_left,diff_right,minmod

  do idim=1,ndim
     do i=1,nn
        diff_left=0.5*(a(i,2*idim)-a(i,0))
        diff_right=0.5*(a(i,0)-a(i,2*idim-1))
        if(diff_left*diff_right<=0.0)then
           minmod=0.0
        else
           minmod=MIN(ABS(diff_left),ABS(diff_right)) &
                &   *diff_left/ABS(diff_left)
        end if
        w(i,idim)=minmod
     end do
  end do

end subroutine compute_limiter_minmod
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine compute_limiter_central(a,w,nn)
  use amr_commons
  use hydro_commons
  implicit none
  integer::nn
  real(dp),dimension(1:nvector,0:twondim)::a
  real(dp),dimension(1:nvector,1:ndim)::w
  !---------------------------
  ! Monotonized Central slope
  !---------------------------
  integer::i,j,idim,ind,ix,iy,iz
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp)::xxc
  real(dp),dimension(1:nvector,1:twotondim),save::ac
  real(dp),dimension(1:nvector),save::corner,kernel,diff_corner,diff_kernel
  real(dp),dimension(1:nvector),save::max_limiter,min_limiter,limiter

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)
  end do

  ! Second order central slope
  do idim=1,ndim
     do i=1,nn
        w(i,idim)=0.25D0*(a(i,2*idim)-a(i,2*idim-1))
     end do
  end do

  ! Compute corner interpolated values
  do ind=1,twotondim
     do i=1,nn
        ac(i,ind)=a(i,0)
     end do
  end do
  do idim=1,ndim
     do ind=1,twotondim
        xxc = xc(ind,idim)
        do i=1,nn
           corner(i)=ac(i,ind)+2d0*w(i,idim)*xxc
        end do
        do i=1,nn
           ac(i,ind)=corner(i)
        end do
     end do
  end do

  ! Compute max of corners
  do i=1,nn
     corner(i)=ac(i,1)
  end do
  do j=2,twotondim
     do i=1,nn
        corner(i)=MAX(corner(i),ac(i,j))
     end do
  end do

  ! Compute max of gradient kernel
  do i=1,nn
     kernel(i)=a(i,1)
  end do
  do j=2,twondim
     do i=1,nn
        kernel(i)=MAX(kernel(i),a(i,j))
     end do
  end do

  ! Compute differences
  do i=1,nn
     diff_kernel(i)=a(i,0)-kernel(i)
     diff_corner(i)=a(i,0)-corner(i)
  end do

  ! Compute max_limiter
  max_limiter=0.0D0
  do i=1,nn
     if(diff_kernel(i)*diff_corner(i) > 0.0D0)then
        max_limiter(i)=MIN(1.0_dp,diff_kernel(i)/diff_corner(i))
     end if
  end do

  ! Compute min of corners
  do i=1,nn
     corner(i)=ac(i,1)
  end do
  do j=2,twotondim
     do i=1,nn
        corner(i)=MIN(corner(i),ac(i,j))
     end do
  end do

  ! Compute min of gradient kernel
  do i=1,nn
     kernel(i)=a(i,1)
  end do
  do j=2,twondim
     do i=1,nn
        kernel(i)=MIN(kernel(i),a(i,j))
     end do
  end do

  ! Compute differences
  do i=1,nn
     diff_kernel(i)=a(i,0)-kernel(i)
     diff_corner(i)=a(i,0)-corner(i)
  end do

  ! Compute max_limiter
  min_limiter=0.0D0
  do i=1,nn
     if(diff_kernel(i)*diff_corner(i) > 0.0D0)then
        min_limiter(i)=MIN(1.0_dp,diff_kernel(i)/diff_corner(i))
     end if
  end do

  ! Compute limiter
  do i=1,nn
     limiter(i)=MIN(min_limiter(i),max_limiter(i))
  end do

  ! Correct gradient with limiter
  do idim=1,ndim
     do i=1,nn
        w(i,idim)=w(i,idim)*limiter(i)
     end do
  end do

end subroutine compute_limiter_central
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine compute_central(a,w,nn)
  use amr_commons
  use hydro_commons
  implicit none
  integer::nn
  real(dp),dimension(1:nvector,0:twondim)::a
  real(dp),dimension(1:nvector,1:ndim)::w
  !---------------------------
  ! Unlimited Central slope
  !---------------------------
  integer::i,idim

  ! Second order central slope
  do idim=1,ndim
     do i=1,nn
        w(i,idim)=0.25D0*(a(i,2*idim)-a(i,2*idim-1))
     end do
  end do

end subroutine compute_central
