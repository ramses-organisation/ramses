program get_music_refmask
  !--------------------------------------------------------------------------
  ! Extract a mask to use in combination with MUSIC.
  ! Version F90 by R. Teyssier on 01/04/01.
  ! Modifier by C. Cadiou on 27/11/17.
  !--------------------------------------------------------------------------
  use utils
  implicit none
  integer::ncpu,ndim,npart,i,j,k,icpu,ipos,nstar,nstart,ico
  integer::ncpu2,npart2,ndim2,levelmin,levelmax,ilevel
  integer::ncpu_read
  real(KIND=8)::mtot,t,mass,btime,unit_l,aexp,unit_t
  real(KIND=8)::xmin=0,xmax=1,ymin=0,ymax=1,zmin=0,zmax=1,r,xc=0.5,yc=0.5,zc=0.5,rad=-1
  integer::imin,imax,jmin,jmax,kmin,kmax,lmin,ipart
  real(KIND=8)::deltax
  real(KIND=8),dimension(:)  ,allocatable::x
  real(KIND=8),dimension(:)  ,allocatable::y
  real(KIND=8),dimension(:)  ,allocatable::z
  real(KIND=8),dimension(:)  ,allocatable::vx
  real(KIND=8),dimension(:)  ,allocatable::vy
  real(KIND=8),dimension(:)  ,allocatable::vz
  real(KIND=8),dimension(:)  ,allocatable::m
  real(KIND=8),dimension(:)  ,allocatable::temp
  real(KIND=8),dimension(:)  ,allocatable::bt
  real(KIND=8),dimension(:)  ,allocatable::tempx,tempy,tempz,tempvx,tempvy,tempvz,tempm,tempbt
  integer(i8b),allocatable,dimension(:)::tempid
  integer,allocatable,dimension(:)::temp2,indtempid
  integer(i8b),allocatable,dimension(:)::id
  integer(i8b)::maxid
  integer(i8b),allocatable,dimension(:)::idpart
  integer,allocatable,dimension(:)::indidpart
  integer::outputmode
  character(LEN=1)::proj='z'
  character(LEN=5)::nchar,ncharcpu
  character(LEN=80)::ordering
  character(LEN=128)::nomfich,repository,filetype='bin',grafic
  character(LEN=128)::repository2, outputname
  logical::ok,ok_part,periodic=.false.
  integer::impi,ndom,bit_length,maxdom
  integer,dimension(1:8)::idom,jdom,kdom,cpu_min,cpu_max
  real(KIND=8),dimension(1:8)::bounding_min,bounding_max
  real(KIND=8)::dkey,order_min,dmax
  real(kind=8),dimension(:),allocatable::bound_key
  logical,dimension(:),allocatable::cpu_read
  integer,dimension(:),allocatable::cpu_list

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

  if(rad>0) then
     xmin=xc-rad
     xmax=xc+rad
     ymin=yc-rad
     ymax=yc+rad
     zmin=zc-rad
     zmax=zc+rad
  endif

  if(TRIM(ordering).eq.'hilbert')then

     dmax=max(xmax-xmin,ymax-ymin,zmax-zmin)
     do ilevel=1,levelmax
        deltax=0.5d0**ilevel
        if(deltax.lt.dmax)exit
     end do
     lmin=ilevel
     bit_length=lmin-1
     maxdom=2**bit_length
     imin=0; imax=0; jmin=0; jmax=0; kmin=0; kmax=0
     if(bit_length>0)then
        imin=int(xmin*dble(maxdom))
        imax=imin+1
        jmin=int(ymin*dble(maxdom))
        jmax=jmin+1
        kmin=int(zmin*dble(maxdom))
        kmax=kmin+1
     endif

     dkey=(dble(2**(levelmax+1)/dble(maxdom)))**ndim
     ndom=1
     if(bit_length>0)ndom=8
     idom(1)=imin; idom(2)=imax
     idom(3)=imin; idom(4)=imax
     idom(5)=imin; idom(6)=imax
     idom(7)=imin; idom(8)=imax
     jdom(1)=jmin; jdom(2)=jmin
     jdom(3)=jmax; jdom(4)=jmax
     jdom(5)=jmin; jdom(6)=jmin
     jdom(7)=jmax; jdom(8)=jmax
     kdom(1)=kmin; kdom(2)=kmin
     kdom(3)=kmin; kdom(4)=kmin
     kdom(5)=kmax; kdom(6)=kmax
     kdom(7)=kmax; kdom(8)=kmax

     do i=1,ndom
        if(bit_length>0)then
           call hilbert3d(idom(i),jdom(i),kdom(i),order_min,bit_length,1)
        else
           order_min=0.0d0
        endif
        bounding_min(i)=(order_min)*dkey
        bounding_max(i)=(order_min+1.0D0)*dkey
     end do
     cpu_min=0; cpu_max=0
     do impi=1,ncpu
        do i=1,ndom
           if (   bound_key(impi-1).le.bounding_min(i).and.&
                & bound_key(impi  ).gt.bounding_min(i))then
              cpu_min(i)=impi
           endif
           if (   bound_key(impi-1).lt.bounding_max(i).and.&
                & bound_key(impi  ).ge.bounding_max(i))then
              cpu_max(i)=impi
           endif
        end do
     end do

     ncpu_read=0
     do i=1,ndom
        do j=cpu_min(i),cpu_max(i)
           if(.not. cpu_read(j))then
              ncpu_read=ncpu_read+1
              cpu_list(ncpu_read)=j
              cpu_read(j)=.true.
           endif
        enddo
     enddo
  else
     ncpu_read=ncpu
     do j=1,ncpu
        cpu_list(j)=j
     end do
  end  if

  npart=0
  do k=1,ncpu_read
     write(*,*) 'CPU=',k
     icpu=cpu_list(k)
     call title(icpu,ncharcpu)
     nomfich=TRIM(repository)//'/part_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
     open(unit=1,file=nomfich,status='old',form='unformatted')
     read(1)ncpu2
     read(1)ndim2
     read(1)npart2
     read(1)
     read(1)nstar
     close(1)
     npart=npart+npart2
  end do
  write(*,*) npart,' particles in the region'
  allocate(m(1:npart))
  allocate(x(1:npart))
  allocate(y(1:npart))
  allocate(z(1:npart))
  allocate(vx(1:npart))
  allocate(vy(1:npart))
  allocate(vz(1:npart))
  allocate(id(1:npart))
  if(nstar>0) then
     allocate(bt(1:npart))
  endif

  !-----------------------------------------------
  ! Compute projected mass using CIC smoothing
  !----------------------------------------------
  mtot=0.0d0
  nstart=1
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
     allocate(temp(1:npart2))
     allocate(tempid(1:npart2))
     allocate(temp2(1:npart2))
     ! Read positions
     read(1)temp
     x(nstart:nstart+npart2-1)=temp
     read(1)temp
     y(nstart:nstart+npart2-1)=temp
     read(1)temp
     z(nstart:nstart+npart2-1)=temp
     ! Read velocity
     read(1)temp
     vx(nstart:nstart+npart2-1)=temp
     read(1)temp
     vy(nstart:nstart+npart2-1)=temp
     read(1)temp
     vz(nstart:nstart+npart2-1)=temp
!     Read mass
     read(1)temp
     m(nstart:nstart+npart2-1)=temp
     !Read identity
     read(1)tempid
     id(nstart:nstart+npart2-1)=tempid
     !Read level
     read(1)temp2
     if(nstar>0) then
        ! Read BT
        read(1)temp
        bt(nstart:nstart+npart2-1)=temp
     endif
! ----------------------------
     nstart=nstart+npart2  !Fill up the next set
     deallocate(temp)
     deallocate(tempid)
     deallocate(temp2)
  enddo

50 format(2I16)
  !Outputs IDs of selected particles
  ipart=0
  mass=0.0
  write(*,*) 'Getting IDs...'
  open(18,file='partID.dat',form='formatted')
  maxid=0
  do i=1,npart !To get maximum identity of the particle
     if(nstar.eq.0) then  !Only DM particles
        btime=0
     else
        btime=bt(i)
     endif
     if(btime.eq.0) then
        ok_part=(x(i)>=xmin.and.x(i)<=xmax.and. &
             &   y(i)>=ymin.and.y(i)<=ymax.and. &
             &   z(i)>=zmin.and.z(i)<=zmax)
        if(rad>0) then
           r=(x(i)-xc)**2+(y(i)-yc)**2+(z(i)-zc)**2
           ok_part=(sqrt(r)<=rad)
        endif
        if(ok_part) then
           maxid=max(maxid,id(i))
           ipart=ipart+1
           mass=mass+m(i)
        endif
     endif
  enddo
  write(*,*) 'We have',ipart,' particles in selected region'
  write(*,*) 'Total mass =', mass

30 format(i16)
  write(18,50) ipart,npart,maxid

  allocate(idpart(1:ipart))
  j=1
  do i=1,npart  !Start finding the IDs
     if(nstar.eq.0) then  !Only DM particles
        btime=0
     else
        btime=bt(i)
     endif
     if(btime.eq.0) then
        ok_part=(x(i)>=xmin.and.x(i)<=xmax.and. &
             &   y(i)>=ymin.and.y(i)<=ymax.and. &
             &   z(i)>=zmin.and.z(i)<=zmax)
        if(rad>0) then
           r=(x(i)-xc)**2+(y(i)-yc)**2+(z(i)-zc)**2
           ok_part=sqrt(r)<=rad
        endif
        if(ok_part) then
           write(18,30) id(i)   !Write IDs
           idpart(j) = id(i)
           j = j+1
        endif
     endif
  enddo

  npart = ipart
  allocate(indidpart(1:npart))
  call quick_sort(idpart, indidpart, npart)

  !------- read IC data and match -----------

  deallocate(m)
  deallocate(x)
  deallocate(y)
  deallocate(z)
  deallocate(vx)
  deallocate(vy)
  deallocate(vz)
  deallocate(id)
  if(nstar>0) then
     deallocate(bt)
  endif


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

  ipos=INDEX(repository2,'output_')
  nchar=repository2(ipos+7:ipos+13)

  do icpu=1,ncpu
     call title(icpu,ncharcpu)
     nomfich=TRIM(repository2)//'/part_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
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

  end do

  open(20,file=TRIM(outputname),form='formatted')

  if( outputmode .eq. 1 ) then
     do i=1,ipart
        write(20,1001) x(i), y(i), z(i), vx(i), vy(i), vz(i)
     end do
  else
      do i=1,ipart
        write(20,1002) x(i), y(i), z(i)
     end do
  end if

  close(20)

  write(*,*)'Wrote data to file ',trim(outputname)

1001 format (f16.8,f16.8,f16.8,e16.7,e16.7,e16.7)
1002 format (f16.8,f16.8,f16.8)

contains

  subroutine read_params

      implicit none

      integer       :: i,n
      character(len=4)   :: opt
      character(len=128) :: arg

      n = command_argument_count()
      if (n < 4) then
         print *, 'usage: get_music_refmask  -inf  input_dir_final_snapshot'
         print *, '                          -ini  input_dir_first_snapshot'
         print *, '                          [-out output_name] '
         print *, '                          [-xc xc] '
         print *, '                          [-yc yc] '
         print *, '                          [-zc zc] '
         print *, '                          [-rad rad] '
         print *, '                          [-vel 0|1] output also velocity data'
         print *, 'ex: get_music_refmask -inf output_00010 -ini output_00001 -xc 0.5 -yc 0.5 -zc 0.5 -rad 0.1'
         stop
      end if

      outputname = 'music_region_file.txt'
      outputmode = 0

      do i = 1,n,2
         call get_command_argument(i, opt)
         if (i == n) then
            print '("option ",a2," has no argument")', opt
            stop 2
         end if
         call get_command_argument(i+1, arg)
         select case (opt)
         case ('-inf')
            repository = trim(arg)
         case ('-ini')
            repository2 = trim(arg)
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
         case ('-out')
            outputname = trim(arg)
         case ('-vel')
            read (arg,*) outputmode
         case ('-per')
            read (arg,*) periodic
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

end program get_music_refmask
