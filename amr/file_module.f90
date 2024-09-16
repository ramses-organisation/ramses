module file_module

  interface
    integer(c_int) function c_mkdir(name,mode) BIND(C,NAME="mkdir")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_CHAR, C_INT
      character(kind=c_char), dimension(*) :: name
      integer(C_INT), VALUE :: mode
    end function c_mkdir
  end interface

  contains

    ! Make a directory
    subroutine mkdir(name,mode,errno,parent)
      use iso_c_binding, only: c_char, c_null_char
      character(len=*), intent(in) :: name
      integer, value :: mode
      integer, optional :: errno
      logical, optional :: parent
      integer :: rc,i
      character(kind=c_char), dimension(len(name)+1) :: c_name

      if (present(errno)) errno = 0
      do i = 1, len(name)
        if (present(parent)) then
          if (i.gt.1 .and. name(i:i).eq.'/') then
            c_name(i) = c_null_char
            rc = c_mkdir(c_name,mode)
            if (rc .ne. 0 .and. ierrno() .ne. 17) then
              if (present(errno)) errno = ierrno()
            end if
          endif
        end if
        c_name(i) = name(i:i)
      end do
      c_name(len(name)+1) = c_null_char
      rc = c_mkdir(c_name,mode)
      if (rc .ne. 0 .and. ierrno() .ne. 17) then
        if (present(errno)) errno = ierrno()
      end if
    end subroutine mkdir

end module file_module
