subroutine output_patch(filename)
  character(LEN=80)::filename
  character(LEN=80)::fileloc
  character(LEN=30)::format
  integer::ilun

  ilun=11

  fileloc=TRIM(filename)
  format="(A)"
  open(unit=ilun,file=fileloc,form='formatted')
  write(ilun,format)"no patches"
  write(ilun,format)" "
  close(ilun)
end subroutine output_patch
