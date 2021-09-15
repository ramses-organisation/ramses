!----------------------------------------------------------------------------------------------
! basic sparse matrix package for the use in the RAMSES clumpfinder
! every line of the matrix is saved as a linked list
! issues/improvements:
!    -write a routine for quicker "maxmerging" two lines (outer loop calling "get_value" at each position
!     is very iniefficient)
!    -disconnect a value which is set to zero rather than just writing zero into memory (not too bad since
!    -have short lifetime)
!    -reuse disconnected space
!----------------------------------------------------------------------------------------------

module sparse_matrix
  use amr_commons
  implicit none
  type sparse_mat
     real(dp),allocatable,dimension(:)::val,maxval
     integer,allocatable,dimension(:)::next,col,first,maxloc
     integer::used,n,m
  end type sparse_mat
  integer,parameter::NSPARSEMAX=10000000
contains
  !----------------------------------------------------------------------------------------------
  subroutine sparse_initialize(m,mat)
    type(sparse_mat)::mat
    integer::m !size of the array
    allocate(mat%val(1:NSPARSEMAX))
    mat%val=0
    allocate(mat%next(1:NSPARSEMAX))
    mat%next=0
    allocate(mat%col(0:NSPARSEMAX))
    mat%col=0
    mat%col(0)=huge(1)
    allocate(mat%first(1:m))
    mat%first=0
    allocate(mat%maxval(1:m))
    mat%maxval=0
    allocate(mat%maxloc(1:m))
    mat%maxloc=0
    mat%used=0
  end subroutine sparse_initialize
  !----------------------------------------------------------------------------------------------

  !----------------------------------------------------------------------------------------------
  subroutine sparse_kill(mat)
    type(sparse_mat)::mat
    deallocate(mat%val,mat%maxval)
    deallocate(mat%next,mat%col,mat%first,mat%maxloc)
  end subroutine sparse_kill
  !----------------------------------------------------------------------------------------------

  !----------------------------------------------------------------------------------------------
  subroutine set_value(i,j,new_value,mat)
    type(sparse_mat)::mat
    integer::i,j !new entry
    real(dp)::new_value
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! set new_value at position i,j in mat
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer::current,save_next

    if(mat%used.eq.NSPARSEMAX)then
       write(*,*)'Maximum size reached',mat%used
       stop
    endif

    ! if corresponding line is empty
    if (mat%first(i)==0) then
       if(new_value.NE.0.)then
          mat%used=mat%used+1
          mat%first(i)=mat%used
          mat%val(mat%used)=new_value
          mat%col(mat%used)=j
          mat%next(mat%used)=0
       end if
       return
    end if

    ! if element needs to be added to start of the list
    if (mat%col(mat%first(i))>j) then
       if(new_value.NE.0.)then
          mat%used=mat%used+1
          save_next=mat%first(i)
          mat%first(i)=mat%used
          mat%val(mat%used)=new_value
          mat%col(mat%used)=j
          mat%next(mat%used)=save_next
       endif
       return
    end if

    ! if first element needs to be replaced -> update max if needed
    if (mat%col(mat%first(i)) == j) then
       mat%val(mat%first(i))=new_value
       return
    end if

    ! walk the line...
    current=mat%first(i)
    do  while( mat%col(mat%next(current)) < j )
       current=mat%next(current)
    end do

    ! there is already an existing value in place: overwrite and update max if needed
    if ( mat%col(mat%next(current)) == j)then
       mat%val(mat%next(current))=new_value
       return
    end if

    ! next points to zero -> add to the end
    if ( mat%next(current) == 0)then
       if(new_value.NE.0.)then
          mat%used=mat%used+1
          mat%next(current)=mat%used
          mat%next(mat%used)=0
          mat%col(mat%used)=j
          mat%val(mat%used)=new_value
       endif
       return
    end if

    ! next points not to zero -> link in between
    if ( mat%next(current) > 0)then
       if(new_value.NE.0.)then
          mat%used=mat%used+1
          save_next=mat%next(current)
          mat%next(current)=mat%used
          mat%next(mat%used)=save_next
          mat%col(mat%used)=j
          mat%val(mat%used)=new_value
       endif
       return
    end if

    write(*,*)'ooops, I should not be here!'
  end subroutine set_value
  !----------------------------------------------------------------------------------------------

  !----------------------------------------------------------------------------------------------
  function get_value(i,j,mat)
    type(sparse_mat)::mat
    integer::i,j
    real(dp)::get_value
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! gets the value of mat at position i,j
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer::current

    ! empty line
    get_value=0
    if (mat%first(i)==0)then
       return
    end if

    ! walk the line...
    current=mat%first(i)
    do  while( mat%col(current) < j )
       current=mat%next(current)
    end do

    ! we are sitting at the right spot
    if ( mat%col(current) == j)then
       get_value= mat%val(current)
       return
    end if

    ! we are too far (means there was no entry for column j)
    if ( mat%col(current) > j)then
       get_value=0
       return
    end if

    write(*,*)'ooops, I should not be here!'
  end function get_value
  !----------------------------------------------------------------------------------------------

end module sparse_matrix
