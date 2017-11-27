program halo_evolution
  implicit none

  integer::ncpu,ndim,npart,i,j,n,n_2,icpu,ipos,nstar,ngroup,halo_nr,k,cou,ida,test_1,test_2,ngroup_2,group_ida_2
  integer::ncpu2,npart2,ndim2,id_group,h,npart3,npart4,npart5,npart6,n_list,npart7,outf
  integer,dimension(:),allocatable::group_id,n_group,group_id_2,group_id_3,list,group_id_4,npart_halo
  real(kind=8),dimension(:),allocatable::xref_group,yref_group,zref_group
  logical,dimension(:),allocatable::flag_group
  real(KIND=8)::mtot,mcut=1000000,sthres1,sthres2,aexp
  real(KIND=8)::period=1,xmin=0,xmax=1,ymin=0,ymax=1,zmin=0,zmax=1,rvir
  real(kind=8),dimension(:),allocatable::m_group,x_group,y_group,z_group,u_group,v_group,w_group,mpure_group,contamine
  real(KIND=8),dimension(:,:),allocatable::x,v
  real(KIND=8),dimension(:),allocatable::m,age
  integer,dimension(:),allocatable::id,id_new,id_new_2,ind_new_1,ind_new_2,ind_new_3
  character(LEN=5)::nstring,ncharcpu,noutput
  character(LEN=80)::directory,file_groupe_in,hop_name
  character(LEN=128)::nomfich,repository
  logical::ok,first
  integer::verbose=0

  !1. Specify Halo number in final output and mass threshold.

  call read_params

  !2. Extract particle id of particles in this halo, from final output and also the halo proprties (ramses output, .tag and .pos needed)

  !2.a) Reading group numbers

  call title(outf,noutput)
  nomfich=TRIM(hop_name)//TRIM(noutput)//'.tag'
  inquire(file=nomfich, exist=ok) ! verify input file
  if ( .not. ok ) then
     print *,TRIM(nomfich)//' not found'
     stop
  endif
  !write(*,*)'Reading file '//TRIM(nomfich)
#ifdef DIR
  open(unit=15,file=trim(nomfich),form='unformatted',status='old',access='direct',recl=1)
  read(15,rec=1)npart
  read(15,rec=2)ngroup
#else
  Open(15,file=nomfich,FORM='UNFORMATTED')
  Read(15)npart, ngroup
#endif

  if(verbose==1)write(*,*)npart
  allocate(group_id(1:npart))

#ifdef DIR
  do j=1,npart
     read(15,rec=j+2)group_id(j)
     if(mod(j,1000000).eq.0.and.verbose==1)write(*,*)j,group_id(j)
  end do
#else
  Read(15)(group_id(j),j=1,npart)
#endif
  Close(15)
  group_id=group_id+1  ! WARNING: convert from (0 to ng-1) to (1 to ng)

  id_group=0
  npart3=0

  do j=1,npart
     id_group=group_id(j)

     if(id_group.EQ.halo_nr) then
        npart3=npart3+1
     end if
  end do

  if(verbose==1)write(*,*) 'Number of particles in the selected halo =',npart3

  if(verbose==1)write(*,*)'Reading initial particle ids'
  nomfich='output_'//TRIM(noutput)//'/info_'//TRIM(noutput)//'.txt'
  inquire(file=nomfich, exist=ok) ! verify input file
  if ( .not. ok ) then
     print *,TRIM(nomfich)//' not found'
     stop
  endif
  open(unit=10,file=nomfich,form='formatted',status='old')
  read(10,'(13X,I11)')ncpu
  read(10,'(13X,I11)')ndim
  read(10,'(13X,I11)')
  read(10,'(13X,I11)')
  read(10,'(13X,I11)')
  read(10,'(13X,I11)')
  read(10,'(13X,I11)')
  read(10,'(13X,I11)')
  read(10,'(13X,I11)')
  read(10,'(13X,E23.15)')aexp
  close(10)

  !2.b)Extract halo properties for this list of halos from the .pos file

  allocate(group_id_4(1:ngroup))
  allocate(n_group(1:ngroup))
  allocate(m_group(1:ngroup))
  allocate(contamine(1:ngroup))
  allocate(x_group(1:ngroup))
  allocate(y_group(1:ngroup))
  allocate(z_group(1:ngroup))
  allocate(u_group(1:ngroup))
  allocate(v_group(1:ngroup))
  allocate(w_group(1:ngroup))

  nomfich=TRIM(hop_name)//TRIM(noutput)//'.pos'
  inquire(file=nomfich, exist=ok) ! verify input file
  if ( .not. ok ) then
     print *,TRIM(nomfich)//' not found'
     stop
  endif
  Open(16,file=nomfich,FORM='FORMATTED',status='old')
  read(16,'(A100)')
  do i=1,ngroup
     read(16,*) &
          & group_id_4(i), n_group(i), m_group(i), contamine(i), &
          & x_group(i),y_group(i),z_group(i), &
          & u_group(i),v_group(i),w_group(i)
  end do
  Close(16)

  !2.c Output
  !write(*,*)
  !write(*,*)"Working with output number ", noutput
  !write(*,'(A119)')'    #     aexp       frac     npart     mass    cont.frac      xc         yc         zc        uc         vc         wc'
  do i=1,ngroup
     if(group_id_4(i).EQ.halo_nr) then
        sthres2=1.0
        if(sthres2>sthres1) then
           write(*,'(I5,A,1E10.2,A,1E10.3,A,I7,A,1PE10.3,A,1E10.3,A,1E10.3,A,1E10.3,A,1E10.3,A,1E10.3,A,1E10.3,A,1E10.3)') &
                & group_id_4(i),' ', aexp,' ',sthres2,' ', n_group(i),' ', m_group(i),' ', contamine(i),' ', x_group(i),' ', &
                      & y_group(i),' ', z_group(i),' ', u_group(i),' ' , v_group(i),' ', w_group(i)
        end if
        cycle
     endif
  end do
  deallocate(group_id_4)
  deallocate(n_group)
  deallocate(m_group)
  deallocate(contamine)
  deallocate(x_group)
  deallocate(y_group)
  deallocate(z_group)
  deallocate(u_group)
  deallocate(v_group)
  deallocate(w_group)

  !2.d) Extracting particle ids

  npart=0
  do icpu=1,ncpu
     call title(icpu,ncharcpu)
     nomfich='output_'//TRIM(noutput)//'/part_'//TRIM(noutput)//'.out'//TRIM(ncharcpu)
     inquire(file=nomfich, exist=ok) ! verify input file
     if ( .not. ok ) then
        print *,TRIM(nomfich)//' not found'
        stop
     endif
     open(unit=1,file=nomfich,status='old',form='unformatted')
     read(1)ncpu2
     read(1)ndim2
     read(1)npart2
     read(1)
     read(1)

     close(1)
     npart=npart+npart2

  end do

  allocate(id_new(1:npart3))
  allocate(ind_new_1(1:npart3))
  id_group=0
  npart=0
  npart4=0
  k=1
  do icpu=1,ncpu
     call title(icpu,ncharcpu)
     nomfich='output_'//TRIM(noutput)//'/part_'//TRIM(noutput)//'.out'//TRIM(ncharcpu)
     inquire(file=nomfich, exist=ok) ! verify input file
     if ( .not. ok ) then
        print *,TRIM(nomfich)//' not found'
        stop
     endif
     open(unit=1,file=nomfich,status='old',form='unformatted')
     read(1)ncpu2
     read(1)ndim2
     read(1)npart2
     read(1)
     read(1)
     read(1)
     read(1)
     read(1)
     allocate(m(1:npart2))
     allocate(id(1:npart2))
     !Read position
     do i=1,ndim
        read(1)m
     end do
     !Read velocity
     do i=1,ndim
        read(1)m
     end do
     !Read mass
     read(1)m
     read(1)id
     read(1)
     close(1)

     do i=1,npart2
        npart=npart+1
        id_group=group_id(npart)
        if(id_group.EQ.halo_nr) then
           id_new(k)=id(i)
           k=k+1
           npart4=npart4+1
        end if
     end do

     deallocate(m)
     deallocate(id)
  end do

  deallocate(group_id)

  !3. Sort particles in the desired halo by the number of their id

  if(verbose==1)write(*,*)'Sorting initial particle ids'
  !call quicksort(id_new,1,npart3)
  call quick_sort_int(id_new,ind_new_1,npart3)
  if(verbose==1)write(*,*)'Sort done.'

  !4. Search previous outputs for particles with these ids

  !4.a) Read in group numbers

  do n=1,(outf-1)
     n_2=(outf-n)
     call title(n_2,noutput)
     !write(*,*)
     !write(*,*)"Working with output number ", noutput
     nomfich=TRIM(hop_name)//TRIM(noutput)//'.tag'
     inquire(file=nomfich, exist=ok) ! verify input file
     if ( .not. ok ) then
        print *,TRIM(nomfich)//' not found'
        cycle
     endif

#ifdef DIR
     open(unit=15,file=trim(nomfich),form='unformatted',status='old',access='direct',recl=1)
     read(15,rec=1)npart5
     read(15,rec=2)ngroup_2
#else
     Open(15,file=trim(nomfich),FORM='UNFORMATTED')
     Read(15)npart5,ngroup_2
#endif

     allocate(group_id_2(1:npart5))

#ifdef DIR
     do j=1,npart5
        read(15,rec=j+2)group_id_2(j)
        if(mod(j,1000000).eq.0.and.verbose==1)write(*,*)j,group_id_2(j)
     end do
#else
     Read(15)(group_id_2(j),j=1,npart5)
#endif

     Close(15)
     group_id_2=group_id_2+1  ! WARNING: convert from (0 to ng-1) to (1 to ng)

     !4.b) Extract particle ids

     if(verbose==1)write(*,*)'Reading current snapshot particle ids'
     nomfich='output_'//TRIM(noutput)//'/info_'//TRIM(noutput)//'.txt'
     inquire(file=nomfich, exist=ok) ! verify input file
     if ( .not. ok ) then
        print *,TRIM(nomfich)//' not found'
        goto 1001
     endif
     open(unit=10,file=nomfich,form='formatted',status='old')
     read(10,'(13X,I11)')ncpu
     read(10,'(13X,I11)')ndim
     read(10,'(13X,I11)')
     read(10,'(13X,I11)')
     read(10,'(13X,I11)')
     read(10,'(13X,I11)')
     read(10,'(13X,I11)')
     read(10,'(13X,I11)')
     read(10,'(13X,I11)')
     read(10,'(13X,E23.15)')aexp
     close(10)

     npart=0

     do icpu=1,ncpu
        call title(icpu,ncharcpu)
        nomfich='output_'//TRIM(noutput)//'/part_'//TRIM(noutput)//'.out'//TRIM(ncharcpu)
        inquire(file=nomfich, exist=ok) ! verify input file
        if ( .not. ok ) then
           print *,TRIM(nomfich)//' not found'
           goto 1001
        endif
        open(unit=1,file=nomfich,status='old',form='unformatted')
        read(1)ncpu2
        read(1)ndim2
        read(1)npart2
        read(1)
        read(1)
        close(1)
        npart=npart+npart2
     end do
     if(verbose==1)write(*,*)'Found ',npart,' particles.'

     allocate(id_new_2(1:npart))
     allocate(ind_new_2(1:npart))
     k=0

     do icpu=1,ncpu
        call title(icpu,ncharcpu)
        nomfich='output_'//TRIM(noutput)//'/part_'//TRIM(noutput)//'.out'//TRIM(ncharcpu)
        inquire(file=nomfich, exist=ok) ! verify input file
        if ( .not. ok ) then
           print *,TRIM(nomfich)//' not found'
           goto 1002
        endif
        open(unit=1,file=nomfich,status='old',form='unformatted')
        read(1)ncpu2
        read(1)ndim2
        read(1)npart2
        read(1)
        read(1)
        read(1)
        read(1)
        read(1)
        allocate(m(1:npart2))
        allocate(id(1:npart2))
        !Read position
        do i=1,ndim
           read(1)m
        end do
        !Read velocity
        do i=1,ndim
           read(1)m
        end do
        !Read mass
        read(1)m
        read(1)id
        close(1)

        do i=1,npart2
           k=k+1
           id_new_2(k)=id(i)
        end do

        deallocate(m)
        deallocate(id)

     end do

     !4.c) Sort array with particle ids and halo number according to particle ids

     if(verbose==1)write(*,*)'Sorting current snapshot particle ids'
     call quick_sort_int(id_new_2,ind_new_2,npart)
     if(verbose==1)write(*,*)'Sort done.'

     !5. Extract the halos in which they are and the properties of those

     !5.a)Pick the halo numbers of the particles in the selected final halo

     i=1
     allocate(group_id_3(1:npart3))
     allocate(ind_new_3(1:npart3))
     npart6=0
     do j=1,npart
        if(id_new_2(j).EQ.id_new(i)) then
           group_id_3(i)=group_id_2(ind_new_2(j))
           npart6=npart6+1
           i=i+1
        end if
     end do

     !5.b)Sort the halo numbers

     if(verbose==1)write(*,*)'Sorting current halo number'
     call quick_sort_int(group_id_3,ind_new_3,npart3)
     if(verbose==1)write(*,*)'Sort done.'

     !5.c)Check over how many halos the particles in the final halo are distributed in this output and list their numbers

     n_list=1

     do i=2,npart3
        if(group_id_3(i).NE.group_id_3(i-1)) then
           n_list=n_list+1
        end if
     end do

     allocate(list(1:n_list))
     j=1

     do i=1,(npart3-1)
        if(group_id_3(i).NE.group_id_3(i+1)) then
           list(j)=group_id_3(i)
           j=j+1
        end if
     end do

     list(j)=group_id_3(npart3)
     allocate(npart_halo(1:n_list))

     npart_halo=0
     j=1
     npart_halo(1)=1

     do i=2,npart3
        if(group_id_3(i).EQ.group_id_3(i-1)) then
           npart_halo(j)=npart_halo(j)+1
        else
           j=j+1
           npart_halo(j)=npart_halo(j)+1
        end if
     end do

     !5.d)Extract halo properties for this list of halos from the .pos file

     allocate(group_id_4(1:ngroup_2))
     allocate(n_group(1:ngroup_2))
     allocate(m_group(1:ngroup_2))
     allocate(contamine(1:ngroup_2))
     allocate(x_group(1:ngroup_2))
     allocate(y_group(1:ngroup_2))
     allocate(z_group(1:ngroup_2))
     allocate(u_group(1:ngroup_2))
     allocate(v_group(1:ngroup_2))
     allocate(w_group(1:ngroup_2))

     nomfich=TRIM(hop_name)//TRIM(noutput)//'.pos'
     inquire(file=nomfich, exist=ok) ! verify input file
     if ( .not. ok ) then
        print *,TRIM(nomfich)//' not found'
        goto 1003
     endif
     Open(16,file=nomfich,FORM='FORMATTED',status='old')
     read(16,'(A100)')
     do i=1,ngroup_2
        read(16,*) &
             & group_id_4(i), n_group(i), m_group(i), contamine(i),x_group(i),y_group(i),z_group(i), &
             & u_group(i),v_group(i),w_group(i)
     end do
     Close(16)

     !6. Output

     j=1
first=.true.
     do i=1,ngroup_2
        if(j<(n_list+1)) then
           if(list(j).EQ.0) then
              !write(*,'(A119)')'    #     aexp       frac     npart     mass    cont.frac      xc         yc         zc        uc         vc         wc'
              j=j+1
           end if
           if(group_id_4(i).EQ.list(j)) then
              sthres2=REAL(npart_halo(j))/REAL(n_group(i))
              if(sthres2>sthres1.and.first.eq..true.) then
                 write(*,'(I5,A,1E10.2,A,1E10.3,A,I7,A,1PE10.3,A,1E10.3,A,1E10.3,A,1E10.3,A,1E10.3,A,1E10.3,A,1E10.3,A,1E10.3)') &
                      & group_id_4(i),' ', aexp,' ',sthres2,' ', n_group(i),' ', m_group(i),' ', contamine(i),' ', x_group(i),' ', &
                      & y_group(i),' ', z_group(i),' ', u_group(i),' ' , v_group(i),' ', w_group(i)
first=.false.
exit
              end if
              j=j+1
           end if
        end if
     end do

1003 deallocate(group_id_4)
     deallocate(n_group)
     deallocate(m_group)
     deallocate(contamine)
     deallocate(x_group)
     deallocate(y_group)
     deallocate(z_group)
     deallocate(u_group)
     deallocate(v_group)
     deallocate(w_group)

     deallocate(npart_halo)

     deallocate(list)

     deallocate(group_id_3,ind_new_3)

1002 deallocate(id_new_2,ind_new_2)

1001 deallocate(group_id_2)

  end do

  deallocate(id_new)
  deallocate(ind_new_1)

contains

  subroutine read_params

    implicit none

    integer       :: i,n
    
    character(len=4)   :: opt
    character(len=128) :: arg

    n = command_argument_count()
    if (n.le.4) then
       print *, 'Too few arguments'
       stop
    end if

    do i = 1,n,2
       call get_command_argument(i,opt)
       if (i == n) then
          print '("option ",a4," has no argument")', opt
          stop 2
       end if

       call get_command_argument(i+1,arg)

       select case (opt)

       case ('-hal')
          read (arg,*) halo_nr
       case ('-thr')
          read (arg,*) sthres1
       case ('-ouf')
          read (arg,*) outf
       case ('-ver')
          read (arg,*) verbose
       case ('-hop')
          hop_name = trim(arg)
       case default
          print '("unknown option ",a2," ignored")', opt

       end select
    end do

    return

  end subroutine read_params

end program halo_evolution


subroutine title(n,nstring)

  implicit none
  integer::n
  character*5::nstring

  character*1::nchar1
  character*2::nchar2
  character*3::nchar3
  character*4::nchar4
  character*5::nchar5

  if(n.ge.10000)then
    write(nchar5,'(i5)') n
    nstring = nchar5
  elseif(n.ge.1000)then
    write(nchar4,'(i4)') n
    nstring = '0'//nchar4
  elseif(n.ge.100)then
    write(nchar3,'(i3)') n
    nstring = '00'//nchar3
  elseif(n.ge.10)then
    write(nchar2,'(i2)') n
    nstring = '000'//nchar2
  else
    write(nchar1,'(i1)') n
    nstring = '0000'//nchar1
  endif

end subroutine title


subroutine swap(A,C)

  implicit none

  real :: A,C,Anew
  Anew=A
  A=C
  C=Anew

end subroutine


recursive subroutine quicksort(D,start,ent)

  implicit none

  integer :: start,ent
  real,dimension(1:ent)::D
  integer :: i,p
  real :: pivot

  if(ent.LE.start) then
     return
  end if

  pivot=D(start)

  p=start

  do i=(start+1),ent
     if(D(i)<pivot) then
        p=p+1
        call swap(D(i),D(p))
     end if
  end do

  call swap(D(start),D(p))

  call quicksort(D,start,p-1)
  call quicksort(D,p+1,ent)

end subroutine quicksort


recursive subroutine quicksort2(D,E,start,ent)

  implicit none

  integer :: start,ent
  real,dimension(1:ent)::D,E
  integer :: i,p
  real :: pivot

  if(ent.LE.start) then
     return
  end if

  pivot=D(start)

  p=start

  do i=(start+1),ent
     if(D(i)<pivot) then
        p=p+1
        call swap(D(i),D(p))
        call swap(E(i),E(p))
     end if
  end do

  call swap(D(start),D(p))
  call swap(E(start),E(p))

  call quicksort2(D,E,start,p-1)
  call quicksort2(D,E,p+1,ent)

end subroutine
!================================================================
!================================================================
!================================================================
!================================================================
SUBROUTINE quick_sort_int(list, order, n)
  ! Quick sort routine from:
  ! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
  ! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
  ! Modified by Alan Miller to include an associated integer array which gives
  ! the positions of the elements in the original order.
  IMPLICIT NONE
  INTEGER :: n
  INTEGER, DIMENSION (1:n), INTENT(INOUT)  :: list
  INTEGER, DIMENSION (1:n), INTENT(OUT)  :: order

  ! Local variable
  INTEGER :: i

  DO i = 1, n
     order(i) = i
  END DO

  CALL quick_sort_int_1(1, n)

CONTAINS

  RECURSIVE SUBROUTINE quick_sort_int_1(left_end, right_end)

    INTEGER, INTENT(IN) :: left_end, right_end

    !     Local variables
    INTEGER             :: i, j, itemp
    INTEGER             :: reference, temp
    INTEGER, PARAMETER  :: max_simple_sort_size = 6

    IF (right_end < left_end + max_simple_sort_size) THEN
       ! Use interchange sort for small lists
       CALL interchange_sort_int(left_end, right_end)

    ELSE
       ! Use partition ("quick") sort
       reference = list((left_end + right_end)/2)
       i = left_end - 1; j = right_end + 1

       DO
          ! Scan list from left end until element >= reference is found
          DO
             i = i + 1
             IF (list(i) >= reference) EXIT
          END DO
          ! Scan list from right end until element <= reference is found
          DO
             j = j - 1
             IF (list(j) <= reference) EXIT
          END DO


          IF (i < j) THEN
             ! Swap two out-of-order elements
             temp = list(i); list(i) = list(j); list(j) = temp
             itemp = order(i); order(i) = order(j); order(j) = itemp
          ELSE IF (i == j) THEN
             i = i + 1
             EXIT
          ELSE
             EXIT
          END IF
       END DO

       IF (left_end < j) CALL quick_sort_int_1(left_end, j)
       IF (i < right_end) CALL quick_sort_int_1(i, right_end)
    END IF

  END SUBROUTINE quick_sort_int_1


  SUBROUTINE interchange_sort_int(left_end, right_end)

    INTEGER, INTENT(IN) :: left_end, right_end

    !     Local variables
    INTEGER             :: i, j, itemp
    INTEGER             :: temp

    DO i = left_end, right_end - 1
       DO j = i+1, right_end
          IF (list(i) > list(j)) THEN
             temp = list(i); list(i) = list(j); list(j) = temp
             itemp = order(i); order(i) = order(j); order(j) = itemp
          END IF
       END DO
    END DO

  END SUBROUTINE interchange_sort_int

END SUBROUTINE quick_sort_int
