program icmask
  use utils
  !--------------------------------------------------------------------------
  ! Ce programme calcule la carte de densite surfacique projetee
  ! des particules de matiere noire d'une simulation RAMSES.
  ! Version F90 par R. Teyssier le 01/04/01.
  !--------------------------------------------------------------------------
  implicit none
  integer::ncpu,ndim,npart,i,j,k,icpu,ipos,nstar,nparttot
  integer::ncpu2,npart2,ndim2,levelmin,levelmax,ismooth,ico
  integer::ix,iy,iz,ix1,iy1,iz1,ncpu_read,smt=0,rsm=1
  real(KIND=8)::mtot,t,btime,unit_l,aexp,unit_t
  real(KIND=8)::metal=-1
  real(KIND=8)::xmin=0,xmax=1,ymin=0,ymax=1,zmin=0,zmax=1,r,xc=0.5,yc=0.5,zc=0.5,rad=-1
  integer::ipart
  real(KIND=8)::maxdisp
  real(KIND=8),dimension(:)  ,allocatable::x
  real(KIND=8),dimension(:)  ,allocatable::y
  real(KIND=8),dimension(:)  ,allocatable::z
  real(KIND=8),dimension(:)  ,allocatable::vx
  real(KIND=8),dimension(:)  ,allocatable::vy
  real(KIND=8),dimension(:)  ,allocatable::vz
  real(KIND=8),dimension(:)  ,allocatable::m
  real(KIND=8),dimension(:)  ,allocatable::tempx,tempy,tempz,tempvx,tempvy,tempvz,tempm,tempbt
  real(KIND=8),dimension(:)  ,allocatable::bt
  integer ,allocatable,dimension(:)::tempid,temp2,indtempid
  integer ,allocatable,dimension(:)::id
  integer ,allocatable,dimension(:)::idpart,indidpart
  real ,allocatable,dimension(:,:,:)::imark
  character(LEN=1)::proj='z'
  character(LEN=5)::nchar,ncharcpu
  character(LEN=80)::ordering
  character(LEN=128)::nomfich,repository,filetype='bin',grafic
  logical::ok,periodic=.false.,okerode=.false.
  logical::gid=.false.,fid=.false.
  integer::impi,maxid,idd
  real(KIND=8)::vfact
  real(kind=8),dimension(:),allocatable::bound_key
  logical,dimension(:),allocatable::cpu_read
  integer,dimension(:),allocatable::cpu_list
  integer(kind=4)::np1,np2,np3
  real::dx,x1o,x2o,x3o,astart,omegam,omegav,h0,x1or,x2or,x3or,dxor,omegak

  call read_params

  !-----------------------------------------------
  ! Lecture du fichier particules au format RAMSES
  !-----------------------------------------------
  ipos=INDEX(repository,'output_')
  nchar=repository(ipos+7:ipos+13)
  nomfich=TRIM(repository)//'/part_'//TRIM(nchar)//'.out00001'
  inquire(file=nomfich, exist=ok) ! verify input file
  if ( .not. ok ) then
     print *,TRIM(nomfich)//' not found.'
     stop
  endif

  nomfich=TRIM(repository)//'/info_'//TRIM(nchar)//'.txt'
  inquire(file=nomfich, exist=ok) ! verify input file
  if ( .not. ok ) then
     print *,TRIM(nomfich)//' not found.'
     stop
  endif
  open(unit=10,file=nomfich,form='formatted',status='old')
  read(10,'("ncpu        =",I11)')ncpu
  read(10,'("ndim        =",I11)')ndim
  read(10,'("levelmin    =",I11)')levelmin
  read(10,'("levelmax    =",I11)')levelmax
  read(10,*)
  read(10,*)
  read(10,*)
  write(*,*)ncpu,ndim,levelmin,levelmax

  read(10,*)
  read(10,'("time        =",E23.15)')t
  read(10,'("aexp        =",E23.15)')aexp
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,'("unit_l      =",E23.15)')unit_l
  read(10,*)
  read(10,'("unit_t      =",E23.15)')unit_t

  read(10,*)
  read(10,'("ordering type=",A80)') ordering
  read(10,*)
  write(*,'(" ordering type=",A20)') TRIM(ordering)
  allocate(cpu_list(1:ncpu))
  if(TRIM(ordering).eq.'hilbert')then
     allocate(bound_key(0:ncpu))
     allocate(cpu_read(1:ncpu))
     cpu_read=.false.
     do impi=1,ncpu
        read(10,'(I8,1X,E23.15,1X,E23.15)')i,bound_key(impi-1),bound_key(impi)
     end do
  endif
  close(10)

  ncpu_read=ncpu
  do j=1,ncpu
     cpu_list(j)=j
  end do

  write(*,*) 'Reading ID file...'
  open(19,file='partID.dat',form='formatted',status='old')
  read(19,*) npart,nparttot,maxid
  write(*,*) 'Number of particles in selected region =',npart
  allocate(idpart(1:npart))
  allocate(indidpart(1:npart))
  do i=1,npart
     read(19,*) idd
     !idpart(idd)=1
     idpart(i)=idd
  enddo
  close(19)

  CALL quick_sort(idpart, indidpart, npart)

  allocate(m(1:npart))
  allocate(x(1:npart))
  allocate(y(1:npart))
  allocate(z(1:npart))
  allocate(vx(1:npart))
  allocate(vy(1:npart))
  allocate(vz(1:npart))
  allocate(id(1:npart))
  allocate(bt(1:npart))

  !-----------------------------------------------
  ! Compute projected mass using CIC smoothing
  !----------------------------------------------
  mtot=0.0d0
  !nstart=1
  ipart=0
  do k=1,ncpu_read
     icpu=cpu_list(k)
     call title(icpu,ncharcpu)
     nomfich=TRIM(repository)//'/part_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
     open(unit=1,file=nomfich,status='old',form='unformatted')
     write(*,*)'Processing file '//TRIM(nomfich)
     read(1)ncpu2
     read(1)ndim2
     read(1)npart2
     read(1)
     read(1)
     read(1)
     read(1)
     read(1)
     allocate(tempx(1:npart2))
     allocate(tempy(1:npart2))
     allocate(tempz(1:npart2))
     allocate(tempvx(1:npart2))
     allocate(tempvy(1:npart2))
     allocate(tempvz(1:npart2))
     allocate(tempm(1:npart2))
     allocate(tempid(1:npart2))
     allocate(indtempid(1:npart2))
     allocate(temp2(1:npart2))
     allocate(tempbt(1:npart2))
     ! Read positions
     read(1)tempx
     read(1)tempy
     read(1)tempz
     ! Read velocity
     read(1)tempvx
     read(1)tempvy
     read(1)tempvz
     ! Read mass
     read(1)tempm
     ! Read identity
     read(1)tempid
     ! Read level
     read(1)temp2
     if(nstar.gt.0) then
        ! Read BT
        read(1)tempbt
     else
        tempbt=0
     endif

     close(1)

     call quick_sort(tempid, indtempid, npart2)

     ico=1
     i=1
     do while (i.le.npart2.and.ico.le.npart)
        if(tempid(i).lt.idpart(ico))then
           i=i+1
        else
           if(tempid(i).gt.idpart(ico))then
              ico=ico+1
           else
              if(tempid(i).eq.idpart(ico).and.tempbt(i).ne.0)then
                 i=i+1
              else
                 if(tempid(i).eq.idpart(ico).and.tempbt(i).eq.0)then
                    ipart=ipart+1
                    m(ipart)=tempm(indtempid(i))
                    x(ipart)=tempx(indtempid(i))
                    y(ipart)=tempy(indtempid(i))
                    z(ipart)=tempz(indtempid(i))
                    vx(ipart)=tempvx(indtempid(i))
                    vy(ipart)=tempvy(indtempid(i))
                    vz(ipart)=tempvz(indtempid(i))
                    id(ipart)=tempid(i)
                    bt(ipart)=tempbt(indtempid(i))
                    i=i+1
                    ico=ico+1
                 end if
              end if
           end if
        end if
     end do

!     maxxid=1
!     minnid=1e30
!     do i=1,npart2
!     maxxid=max(maxxid,tempid(i))
!     minnid=min(minnid,tempid(i))
!     end do
!     write(*,*)minnid,maxxid
!     icomin=1
!     icomax=npart
!     do i=1,npart
!        if(idpart(i)<=minnid.and.i>=icomin)then
!           icomin=i
!        end if
!        if(idpart(i)>=maxxid.and.i<=icomax)then
!           icomax=i
!        end if
!     end do
!     write(*,*)icomin,icomax,1,npart
!     istart=icomin
!     do i=1,npart2
!        if(nstar==0) then
!           ico=istart
!           do while(ico.le.icomax.and.idpart(ico).ne.tempid(i))
!              ico=ico+1
!           end do
!           if(idpart(ico).eq.tempid(i))then
!              istart=ico+1
!              ipart=ipart+1
!              m(ipart)=tempm(indtempid(i))
!              x(ipart)=tempx(indtempid(i))
!              y(ipart)=tempy(indtempid(i))
!              z(ipart)=tempz(indtempid(i))
!              vx(ipart)=tempvx(indtempid(i))
!              vy(ipart)=tempvy(indtempid(i))
!              vz(ipart)=tempvz(indtempid(i))
!              id(ipart)=tempid(i)
!              bt(ipart)=0
!           end if
!        end if
!        if(nstar>0.and.tempbt(indtempid(i)).eq.0.and.tempid(i)>0) then
!           ico=istart
!           do while(ico.le.icomax.and.idpart(ico).ne.tempid(i))
!              ico=ico+1
!           end do
!           if(idpart(ico).eq.tempid(i))then
!              istart=ico+1
!              ipart=ipart+1
!              m(ipart)=tempm(indtempid(i))
!              x(ipart)=tempx(indtempid(i))
!              y(ipart)=tempy(indtempid(i))
!              z(ipart)=tempz(indtempid(i))
!              vx(ipart)=tempvx(indtempid(i))
!              vy(ipart)=tempvy(indtempid(i))
!              vz(ipart)=tempvz(indtempid(i))
!              id(ipart)=tempid(i)
!              bt(ipart)=tempbt(indtempid(i))
!           end if
!        end if
!     end do
! ----------------------------
     !nstart=nstart+npart2  !Fill up the next set
     deallocate(tempx)
     deallocate(tempy)
     deallocate(tempz)
     deallocate(tempvx)
     deallocate(tempvy)
     deallocate(tempvz)
     deallocate(tempm)
     deallocate(tempid)
     deallocate(indtempid)
     deallocate(temp2)
     deallocate(tempbt)
  enddo

  !do i=1,npart
  !   write(56,*)i,npart,m(i),x(i),y(i),z(i)
  !end do

  grafic =TRIM(grafic)//'/ic_deltab'
  open(10,file=grafic,form='unformatted')
  read (10)np1,np2,np3,dx,x1o,x2o,x3o,astart,omegam,omegav,h0
  close(10)
  write(*,*)'Input array size is:',np1,np2,np3
  write(*,*)'Current box size is',np1*dx,' Mpc'
  allocate(imark(1:np1,1:np2,1:np3))

  x1or=x1o/(unit_l/aexp/3.08e24)
  x2or=x2o/(unit_l/aexp/3.08e24)
  x3or=x3o/(unit_l/aexp/3.08e24)
  dxor=dx/(unit_l/aexp/3.08e24)

  write(*,*)'Total box size ',(unit_l/aexp/3.08e24)

  omegak=1.0-omegam-omegav

  vfact=aexp*fpeebl(aexp) & ! Same scale factor as in grafic1
       & *sqrt(omegam/aexp+omegav*aexp*aexp+omegak)

  write(*,*) 'vfact =', vfact,aexp

  maxdisp=0.
  do i=1,npart
     if(nstar.eq.0) then  !Only DM particles
        btime=0
     else
        btime=bt(i)
     endif
     if(btime.eq.0) then
        ix=int((x(i)-x1or)/dxor) !Lay down perturbed mark
        iy=int((y(i)-x2or)/dxor)
        iz=int((z(i)-x3or)/dxor)
        if(ix.GT.0.and.ix.LE.np1.and.iy.GT.0.and.iy.LE.np2.and.iz.GT.0.and.iz.LE.np3) then
           imark(ix,iy,iz)=1.0
        endif
        !              write(*,*) z(i), vz(i)/vfact/dxor
        x(i)=x(i)-vx(i)/vfact   !Trace back the Zeldovich approx.
        y(i)=y(i)-vy(i)/vfact
        z(i)=z(i)-vz(i)/vfact
        maxdisp=max(maxdisp,abs(vx(i)/vfact),abs(vy(i)/vfact),abs(vz(i)/vfact))
        if(x(i)<0.)x(i)=1.0+x(i)
        if(y(i)<0.)y(i)=1.0+y(i)
        if(z(i)<0.)z(i)=1.0+z(i)
        if(x(i)>1.)x(i)=x(i)-1.0
        if(y(i)>1.)y(i)=y(i)-1.0
        if(z(i)>1.)z(i)=z(i)-1.0
        ix=int((x(i)-x1or)/dxor)
        iy=int((y(i)-x2or)/dxor)
        iz=int((z(i)-x3or)/dxor)
        if(ix.GT.0.and.ix.LE.np1.and.iy.GT.0.and.iy.LE.np2.and.iz.GT.0.and.iz.LE.np3) then
           imark(ix,iy,iz)=1.0
        endif
     endif
  enddo
  write(*,*) 'Max. displacement ',maxdisp/dxor,' cells'
  write(*,*) 'Outputting Grafics file ic_ref_ini'
  open(33,file='ic_ref_ini',form='unformatted')
  write(33) np1,np2,np3,dx,x1o,x2o,x3o,astart,omegam,omegav,h0
  do k=1,np3
     write(33) ((imark(i,j,k),i=1,np1),j=1,np2)
  enddo
  close(33)



  if(smt>0) then !Dilation
     write(*,*) rsm,smt
     do ismooth=1,rsm !Number of smoothing loops
        write(*,*) 'Smoothing ',ismooth
        do i=1,np1
           do j=1,np2
              do k=1,np3
                 if(imark(i,j,k)==1.0) then
                    do ix=i-smt,i+smt
                       ix1=ix
                       if(ix1.LT.1)ix1=ix1+np1
                       if(ix1.GT.np1)ix1=ix1-np1
                       do iy=j-smt,j+smt
                          iy1=iy
                          if(iy1.LT.1)iy1=iy1+np2
                          if(iy1.GT.np2)iy1=iy1-np2
                          do iz=k-smt,k+smt
                             iz1=iz
                             if(iz1.LT.1)iz1=iz1+np3
                             if(iz1.GT.np3)iz1=iz1-np3
                             if(periodic)then
                                r=(i-ix)**2+(j-iy)**2+(k-iz)**2
                                r=r**0.5
                                if(r.LE.smt.and.imark(ix1,iy1,iz1)==0.0) then
                                   imark(ix1,iy1,iz1)=2.0
                                endif
                             else
                                if(ix.GT.0.and.ix.LE.np1.and.iy.GT.0.and.iy.LE.np2.and.iz.GT.0.and.iz.LE.np3) then
                                   r=(i-ix)**2+(j-iy)**2+(k-iz)**2
                                   r=r**0.5
                                   if(r.LE.smt.and.imark(ix,iy,iz)==0.0) then
                                      imark(ix,iy,iz)=2.0
                                   endif
                                endif
                             endif
                          enddo
                       enddo
                    enddo
                 endif
              enddo
           enddo
        enddo
        do i=1,np1
           do j=1,np2
              do k=1,np3
                 if(imark(i,j,k)>0) then
                    imark(i,j,k)=1.0
                 endif
              enddo
           enddo
        end do
     enddo
  endif

  write(*,*) 'Outputting Grafics file ic_ref_dil'
  open(33,file='ic_ref_dil',form='unformatted')
  write(33) np1,np2,np3,dx,x1o,x2o,x3o,astart,omegam,omegav,h0
  do k=1,np3
     write(33) ((imark(i,j,k),i=1,np1),j=1,np2)
  enddo
  close(33)
  if(smt>0) then !Erosion
     write(*,*) rsm,smt
     do ismooth=1,rsm !Number of smoothing loops
        write(*,*) 'Erosion ',ismooth
        do i=1,np1
           do j=1,np2
              do k=1,np3
                 if(imark(i,j,k)==1.0) then
                    okerode=.false.
                    do ix=i-smt,i+smt
                       ix1=ix
                       if(ix1.LE.0)ix1=ix1+np1
                       if(ix1.GT.np1)ix1=ix1-np1
                       do iy=j-smt,j+smt
                          iy1=iy
                          if(iy1.LE.0)iy1=iy1+np2
                          if(iy1.GT.np2)iy1=iy1-np2
                          do iz=k-smt,k+smt
                             iz1=iz
                             if(iz1.LE.0)iz1=iz1+np3
                             if(iz1.GT.np3)iz1=iz1-np3
                             if(periodic)then
                                r=(i-ix)**2+(j-iy)**2+(k-iz)**2
                                r=r**0.5
                                if(r.LE.smt.and.imark(ix1,iy1,iz1)==0.0) then
                                   okerode=.true.
                                endif
                             else
                                if(ix.GT.0.and.ix.LE.np1.and.iy.GT.0.and.iy.LE.np2.and.iz.GT.0.and.iz.LE.np3) then
                                   r=(i-ix)**2+(j-iy)**2+(k-iz)**2
                                   r=r**0.5
                                   if(r.LE.smt.and.imark(ix,iy,iz)==0.0) then
                                      okerode=.true.
                                   endif
                                endif
                             endif
                          enddo
                       enddo
                    enddo
                    if(okerode)imark(i,j,k)=-1
                 endif
              enddo
           enddo
        enddo
        do i=1,np1
           do j=1,np2
              do k=1,np3
                 if(imark(i,j,k)<=0) then
                    imark(i,j,k)=0.0
                 endif
              enddo
           enddo
        end do
     enddo
  endif
  if(smt>0) then !Erosion
     write(*,*) rsm,smt
     write(*,*) 'Final erosion to remove isolated particles'
     smt=1
     do i=1,np1
        do j=1,np2
           do k=1,np3
              if(imark(i,j,k)==1.0) then
                 okerode=.false.
                 do ix=i-smt,i+smt
                    ix1=ix
                    if(ix1.LE.0)ix1=ix1+np1
                    if(ix1.GT.np1)ix1=ix1-np1
                    do iy=j-smt,j+smt
                       iy1=iy
                       if(iy1.LE.0)iy1=iy1+np2
                       if(iy1.GT.np2)iy1=iy1-np2
                       do iz=k-smt,k+smt
                          iz1=iz
                          if(iz1.LE.0)iz1=iz1+np3
                          if(iz1.GT.np3)iz1=iz1-np3
                          if(periodic)then
                             r=(i-ix)**2+(j-iy)**2+(k-iz)**2
                             r=r**0.5
                             if(r.LE.smt.and.imark(ix1,iy1,iz1)==0.0) then
                                okerode=.true.
                             endif
                          else
                             if(ix.GT.0.and.ix.LE.np1.and.iy.GT.0.and.iy.LE.np2.and.iz.GT.0.and.iz.LE.np3) then
                                r=(i-ix)**2+(j-iy)**2+(k-iz)**2
                                r=r**0.5
                                if(r.LE.smt.and.imark(ix,iy,iz)==0.0) then
                                   okerode=.true.
                                endif
                             endif
                          endif
                       enddo
                    enddo
                 enddo
                 if(okerode)imark(i,j,k)=-1
              endif
           enddo
        enddo
     enddo
     do i=1,np1
        do j=1,np2
           do k=1,np3
              if(imark(i,j,k)<=0) then
                 imark(i,j,k)=0.0
              endif
           enddo
        enddo
     end do
     write(*,*) 'Final dilatation to restore initial region'
     do i=1,np1
        do j=1,np2
           do k=1,np3
              if(imark(i,j,k)==1.0) then
                 do ix=i-smt,i+smt
                    ix1=ix
                    if(ix1.LE.0)ix1=ix1+np1
                    if(ix1.GT.np1)ix1=ix1-np1
                    do iy=j-smt,j+smt
                       iy1=iy
                       if(iy1.LE.0)iy1=iy1+np2
                       if(iy1.GT.np2)iy1=iy1-np2
                       do iz=k-smt,k+smt
                          iz1=iz
                          if(iz1.LE.0)iz1=iz1+np3
                          if(iz1.GT.np3)iz1=iz1-np3
                          if(periodic)then
                             r=(i-ix)**2+(j-iy)**2+(k-iz)**2
                             r=r**0.5
                             if(r.LE.smt.and.imark(ix1,iy1,iz1)==0.0) then
                                imark(ix1,iy1,iz1)=2.0
                             endif
                          else
                             if(ix.GT.0.and.ix.LE.np1.and.iy.GT.0.and.iy.LE.np2.and.iz.GT.0.and.iz.LE.np3) then
                                r=(i-ix)**2+(j-iy)**2+(k-iz)**2
                                r=r**0.5
                                if(r.LE.smt.and.imark(ix,iy,iz)==0.0) then
                                   imark(ix,iy,iz)=2.0
                                endif
                             endif
                          endif
                       enddo
                    enddo
                 enddo
              endif
           enddo
        enddo
     enddo
     do i=1,np1
        do j=1,np2
           do k=1,np3
              if(imark(i,j,k)>0) then
                 imark(i,j,k)=1.0
              endif
           enddo
        enddo
     end do
  endif


  write(*,*) 'Outputting Grafics file ic_refmap'
  if(metal>0.)then
     open(33,file='ic_pvar_00001',form='unformatted')
     write(33) np1,np2,np3,dx,x1o,x2o,x3o,astart,omegam,omegav,h0
     do k=1,np3
        write(33) ((real(metal,kind=4)*imark(i,j,k),i=1,np1),j=1,np2)
     enddo
     close(33)
  endif
  open(33,file='ic_refmap',form='unformatted')
  write(33) np1,np2,np3,dx,x1o,x2o,x3o,astart,omegam,omegav,h0
  do k=1,np3
     write(33) ((imark(i,j,k),i=1,np1),j=1,np2)
  enddo
  close(33)

contains

  subroutine read_params

      implicit none

      integer       :: i,n

      character(len=4)   :: opt
      character(len=128) :: arg

      n = command_argument_count()
      if (n < 4) then
         print *, 'usage: geticmask  -inp  input_dir'
         print *, '                 [-smt smt] '
         print *, '                 [-rsm rsm] '
         print *, '                 [-gfc grafic] '
         print *, '                 [-xc xc] '
         print *, '                 [-yc yc] '
         print *, '                 [-zc zc] '
         print *, '                 [-rad rad] '
         print *, '                 [-per per] '
         print *, 'ex: geticmask -inp output_00001 -xc 0.5 -yc 0.5 -zc 0.5 -rad 0.1 -smt 4'// &
              &   ' -rsm 12 -gfc ic_files/AqC_boxlen12p5_n256 -per true'
         stop
      end if

      do i = 1,n,2
         call get_command_argument(i,opt)
         if (i == n) then
            print '("option ",a2," has no argument")', opt
            stop 2
         end if
         call get_command_argument(i+1,arg)
         select case (opt)
         case ('-inp')
            repository = trim(arg)
         case ('-dir')
            proj = trim(arg)
         case ('-xmi')
            read (arg,*) xmin
         case ('-xma')
            read (arg,*) xmax
         case ('-ymi')
            read (arg,*) ymin
         case ('-yma')
            read (arg,*) ymax
         case ('-zmi')
            read (arg,*) zmin
         case ('-zma')
            read (arg,*) zmax
         case ('-xc')
            read (arg,*) xc
         case ('-yc')
            read (arg,*) yc
         case ('-zc')
            read (arg,*) zc
         case ('-rad')
            read (arg,*) rad
         case ('-per')
            read (arg,*) periodic
         case ('-gid')
            read (arg,*) gid
         case ('-fid')
            read (arg,*) fid
         case ('-smt')
            read (arg,*) smt
         case ('-met')
            read (arg,*) metal
         case ('-rsm')
            read (arg,*) rsm
         case ('-gfc')
            grafic = trim(arg)
         case ('-fil')
            filetype = trim(arg)
         case default
            print '("unknown option ",a2," ignored")', opt
         end select
      end do

      return

    end subroutine read_params
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function fy(a)
    implicit none
    !      Computes the integrand
    real(kind=8)::fy
    real(kind=8)::y,a

    y=omegam*(1.d0/a-1.d0) + omegav*(a*a-1.d0) + 1.d0
    fy=1.d0/y**1.5d0

    return
  end function fy

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function fpeebl(a)
    implicit none
    real(kind=8) :: fpeebl,a
    !     Computes the growth factor f=d\log D1/d\log a.
    real(kind=8) :: fact,y,eps

    eps=1.0d-6
    y=omegam*(1.d0/a-1.d0) + omegav*(a*a-1.d0) + 1.d0
    fact=rombint(eps,a,eps)
    fpeebl=(omegav*a*a-0.5d0*omegam/a)/y - 1.d0 + a*fy(a)/fact
    return
  end function fpeebl
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function rombint(a,b,tol)
    implicit none
    real(kind=8)::rombint
    !
    !     Rombint returns the integral from a to b of f(x)dx using Romberg
    !     integration. The method converges provided that f(x) is continuous
    !     in (a,b). The function f must be double precision and must be
    !     declared external in the calling routine.
    !     tol indicates the desired relative accuracy in the integral.
    !
    integer::maxiter=16,maxj=5
    real(kind=8),dimension(100):: g
    real(kind=8)::a,b,tol,fourj
    real(kind=8)::h,error,gmax,g0,g1
    integer::nint,i,j,k,jmax

    h=0.5d0*(b-a)
    gmax=h*(fy(a)+fy(b))
    g(1)=gmax
    nint=1
    error=1.0d20
    i=0
10  i=i+1
    if(.not.  (i>maxiter.or.(i>5.and.abs(error)<tol)))then
       !  Calculate next trapezoidal rule approximation to integral.

       g0=0.0d0
       do k=1,nint
          g0=g0+fy(a+(k+k-1)*h)
       end do
       g0=0.5d0*g(1)+h*g0
       h=0.5d0*h
       nint=nint+nint
       jmax=min(i,maxj)
       fourj=1.0d0

       do j=1,jmax
          ! Use Richardson extrapolation.
          fourj=4.0d0*fourj
          g1=g0+(g0-g(j))/(fourj-1.0d0)
          g(j)=g0
          g0=g1
       enddo
       if (abs(g0).gt.tol) then
          error=1.0d0-gmax/g0
       else
          error=gmax
       end if
       gmax=g0
       g(jmax+1)=g0
       go to 10
    end if
    rombint=g0
    if (i>maxiter.and.abs(error)>tol) &
         &    write(*,*) 'Rombint failed to converge; integral, error=', &
         &    rombint,error
    return
  end function rombint
end program icmask
