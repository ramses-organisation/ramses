module sparse_matrix
  use amr_commons

  type sparse_mat
     real(dp),allocatable,dimension(:)::val,maxval
     integer,allocatable,dimension(:)::next,col,first,maxloc
     integer::used,n,m
  end type sparse_mat


contains

  subroutine sparse_initialize(m,n,mat)
    type(sparse_mat)::mat 
    integer::m,n !size of the array
    allocate(mat%val(1:100000))
    mat%val=0.
    allocate(mat%next(1:100000))
    mat%next=0
    allocate(mat%col(0:100000))
    mat%col=0
    mat%col(0)=huge(1)
    allocate(mat%first(1:m))
    mat%first=0
    allocate(mat%maxval(1:m))
    mat%maxval=0.
    allocate(mat%maxloc(1:m))
    mat%maxloc=0
    mat%used=0
    

  end subroutine sparse_initialize

  subroutine set_value(i,j,new_value,mat)
    type(sparse_mat)::mat 
    integer::i,j !new entry
    real(dp)::new_value

    integer::current,save_next

    ! update maximum and its location
    if (new_value>mat%maxval(i))then
       mat%maxval(i)=new_value
       mat%maxloc(i)=j
    end if
    ! in case of equal values, take smaller index (for comparability)
    if (new_value==mat%maxval(i))then
       mat%maxloc(i)=min(mat%maxloc(i),j)
    end if

    ! start new line
    if (mat%first(i)==0) then
       mat%used=mat%used+1
       mat%first(i)=mat%used
       mat%val(mat%used)=new_value
       mat%col(mat%used)=j
       mat%next(mat%used)=0
       return
    end if

    ! if element needs to be added to start
    if (mat%col(mat%first(i))>j) then
       mat%used=mat%used+1
       save_next=mat%first(i)
       mat%first(i)=mat%used
       mat%val(mat%used)=new_value
       mat%col(mat%used)=j
       mat%next(mat%used)=save_next
       return
    end if

    ! if first element needs to be replaced -> update max if needed
    if (mat%col(mat%first(i)) == j) then
       mat%val(mat%first(i))=new_value
       if (j==mat%maxloc(i))call get_max(i,mat)
       return
    end if

    !walk the line...
    current=mat%first(i)
    do  while( mat%col(mat%next(current)) < j )
       current=mat%next(current)
    end do
    
    !there is already an existing value in place: overwrite and update max if needed
    if ( mat%col(mat%next(current)) == j)then
       mat%val(mat%next(current))=new_value
       if (j==mat%maxloc(i))call get_max(i,mat)
       return
    end if
    
    !next points to zero -> add to the end
    if ( mat%next(current) == 0)then
       mat%used=mat%used+1
       mat%next(current)=mat%used
       mat%next(mat%used)=0
       mat%col(mat%used)=j
       mat%val(mat%used)=new_value       
       return
    end if

    !next point not to zero -> link in between
    if ( mat%next(current) > 0)then
       mat%used=mat%used+1
       save_next=mat%next(current)
       mat%next(current)=mat%used
       mat%next(mat%used)=save_next
       mat%col(mat%used)=j
       mat%val(mat%used)=new_value
       return
    end if

    write(*,*)'ooops, I should not be here!'

  end subroutine set_value


  function get_value(i,j,mat)
    type(sparse_mat)::mat 
    integer::i,j 
    real(dp)::get_value

    integer::current


    !empty line
    if (mat%first(i)==0)then
       get_value=0.
       return
    end if

    !walk the line...
    current=mat%first(i)
    do  while( mat%col(current) < j )
       current=mat%next(current)
    end do
    
    !we are sitting at the right spot
    if ( mat%col(current) == j)then
       get_value= mat%val(current)
       return
    end if
    
    !we are too far (means there was no entry for column j)
    if ( mat%col(current) > j)then
       get_value=0.
       return
    end if

    write(*,*)'ooops, I should not be here!'

  end function get_value

  subroutine get_max(i,mat)
    type(sparse_mat)::mat 
    integer::i

    integer::current

    mat%maxval(i)=0.
    mat%maxloc(i)=0


    !walk the line...
    current=mat%first(i)
    do  while( current /= 0 )
       if(mat%maxval(i)<mat%val(current))then
          mat%maxval(i)=mat%val(current)
          mat%maxloc(i)=mat%col(current)
       end if
       current=mat%next(current)
    end do
    
  end subroutine get_max




  subroutine sparse_kill(mat)
    type(sparse_mat)::mat 
    deallocate(mat%val,mat%maxval)
    deallocate(mat%next,mat%col,mat%first,mat%maxloc)
  end subroutine sparse_kill


end module sparse_matrix
