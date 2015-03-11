#! /bin/bash

#######################################
# cr_write_patch.sh [patch dir]
# creates a .f90 file with code to
# write the patch dir content to disk
#######################################

if [ $# == 0 ]
    then
    echo "no patches" > .tmp_output.txt
    echo " " >> .tmp_output.txt
    else
    PATCHDIRS=$1
    PATCHDIRS=(${PATCHDIRS//:/ })
    for PATCHDIR in "${PATCHDIRS[@]}"; do
	for filename in ${PATCHDIR}/*.f90; do
            echo "$filename" 
            cat "$filename"
	done 
    done > .tmp_output.txt
fi

sed 's/\$/ /g' .tmp_output.txt | sed "s/\"/'/g" | cat -e | sed 's/\$/\"/' | sed 's/^/  write(ilun,format)"/' > .test_middle.f90
  
cat << EOF > .test_after.f90
  close(ilun)
end subroutine output_patch
EOF

cat << EOF > .test_before.f90
subroutine output_patch(filename)
  character(LEN=80)::filename
  character(LEN=80)::fileloc
  character(LEN=30)::format
  integer::ilun

  ilun=11

  fileloc=TRIM(filename)
  format="(A)"
  open(unit=ilun,file=fileloc,form='formatted')
EOF

cat .test_before.f90 .test_middle.f90 .test_after.f90 > write_patch.f90

rm .tmp_output.txt .test_before.f90 .test_middle.f90 .test_after.f90

