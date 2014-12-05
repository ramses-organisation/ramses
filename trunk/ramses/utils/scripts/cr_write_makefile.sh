#! /bin/bash
#######################################
# cr_write_makefile.sh [makefile name]
# creates a .f90 file with code to
# write the Makefile content to disk
#######################################


if [ $# == 0 ]
  then
  exit
fi

MAKEFILE=$1

sed 's/\$/ /g' ${MAKEFILE} | sed "s/\"/'/g" | cat -e | sed 's/\$/\"/' | sed 's/^/  write(ilun,*)\"/' > .test_middle.f90
  
cat << EOF > .test_after.f90

  close(ilun)
end subroutine output_makefile
EOF

cat << EOF > .test_before.f90
subroutine output_makefile(filename)
  
  character(LEN=80)::filename
  character(LEN=80)::fileloc
  integer::ilun
  
  fileloc=TRIM(filename)
  open(newunit=ilun,file=fileloc,form='formatted')
EOF

cat .test_before.f90 .test_middle.f90 .test_after.f90 > write_makefile.f90

rm .test_before.f90 .test_middle.f90 .test_after.f90

