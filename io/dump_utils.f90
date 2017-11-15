module dump_utils
  use iso_fortran_env
  implicit none

  private

  interface generic_dump
     module procedure logicaldump
     module procedure int8dump
     module procedure int16dump
     module procedure int32dump
     module procedure int64dump
     module procedure real32dump
     module procedure real64dump
  end interface generic_dump

  character(len=1), dimension(1:3), parameter :: dim_keys = ["x", "y", "z"]

  public :: generic_dump, dump_header_info, dim_keys
contains

  subroutine logicaldump(varname, ivar, data, unit_out, dump_info, unit_info)
    character(len=*), intent(in) :: varname
    logical, intent(in), dimension(:) :: data
    integer, intent(inout) :: ivar
    integer, intent(in) :: unit_out, unit_info
    logical, intent(in) :: dump_info

    character(len=1) :: kind

    write(unit_out) data
    kind = '?'
    if (dump_info) call dump_var_info(varname, ivar, kind, unit_info)
    ivar = ivar + 1
  end subroutine logicaldump

  subroutine int8dump(varname, ivar, data, unit_out, dump_info, unit_info)
    character(len=*), intent(in) :: varname
    integer(int8), intent(in), dimension(:) :: data
    integer, intent(inout) :: ivar
    integer, intent(in) :: unit_out, unit_info
    logical, intent(in) :: dump_info

    character(len=1) :: kind

    write(unit_out) data
    kind = 'b'
    if (dump_info) call dump_var_info(varname, ivar, kind, unit_info)
    ivar = ivar + 1
  end subroutine int8dump

  subroutine int16dump(varname, ivar, data, unit_out, dump_info, unit_info)
    character(len=*), intent(in) :: varname
    integer(int16), intent(in), dimension(:) :: data
    integer, intent(inout) :: ivar
    integer, intent(in) :: unit_out, unit_info
    logical, intent(in) :: dump_info

    character(len=1) :: kind

    write(unit_out) data
    kind = 'h'
    if (dump_info) call dump_var_info(varname, ivar, kind, unit_info)
    ivar = ivar + 1
  end subroutine int16dump

  subroutine int32dump(varname, ivar, data, unit_out, dump_info, unit_info)
    character(len=*), intent(in) :: varname
    integer(int32), intent(in), dimension(:) :: data
    integer, intent(inout) :: ivar
    integer, intent(in) :: unit_out, unit_info
    logical, intent(in) :: dump_info

    character(len=1) :: kind

    write(unit_out) data
    kind = 'i'
    if (dump_info) call dump_var_info(varname, ivar, kind, unit_info)
    ivar = ivar + 1
  end subroutine int32dump

  subroutine int64dump(varname, ivar, data, unit_out, dump_info, unit_info)
    character(len=*), intent(in) :: varname
    integer(int64), intent(in), dimension(:) :: data
    integer, intent(inout) :: ivar
    integer, intent(in) :: unit_out, unit_info
    logical, intent(in) :: dump_info

    character(len=1) :: kind

    write(unit_out) data
    kind = 'q'
    if (dump_info) call dump_var_info(varname, ivar, kind, unit_info)
    ivar = ivar + 1
  end subroutine int64dump

  subroutine real32dump(varname, ivar, data, unit_out, dump_info, unit_info)
    character(len=*), intent(in) :: varname
    real(real32), intent(in), dimension(:) :: data
    integer, intent(inout) :: ivar
    integer, intent(in) :: unit_out, unit_info
    logical, intent(in) :: dump_info

    character(len=1) :: kind

    write(unit_out) data
    kind = 'f'
    if (dump_info) call dump_var_info(varname, ivar, kind, unit_info)
    ivar = ivar + 1
  end subroutine real32dump

  subroutine real64dump(varname, ivar, data, unit_out, dump_info, unit_info)
    character(len=*), intent(in) :: varname
    real(real64), intent(in), dimension(:) :: data
    integer, intent(inout) :: ivar
    integer, intent(in) :: unit_out, unit_info
    logical, intent(in) :: dump_info

    character(len=1) :: kind

    write(unit_out) data
    kind = 'd'
    if (dump_info) call dump_var_info(varname, ivar, kind, unit_info)
    ivar = ivar + 1
  end subroutine real64dump

  subroutine dump_var_info(varname, ivar, kind, unit_info)
    character(len=*), intent(in) :: varname, kind
    integer, intent(in) :: unit_info
    integer, intent(in) :: ivar

    write(unit_info, '(I2,", ", a, ", ", a)') ivar, trim(varname), trim(kind)
  end subroutine dump_var_info

  subroutine dump_header_info(unit_info)
    integer, intent(in) :: unit_info
    write(unit_info, '("# version: ", i2)') 1
    write(unit_info, '("# ", a, ", ", a, ", ", a)') 'ivar', 'variable_name', 'variable_type'
  end subroutine dump_header_info

end module dump_utils
