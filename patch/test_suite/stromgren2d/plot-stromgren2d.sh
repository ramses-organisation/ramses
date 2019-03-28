#!/bin/bash

../visu/amr2map -inp output_00046 -out data1.dat -typ 1 -fil ascii;

# create the md5 hash
echo $(md5sum data1.dat | cut -d ' ' -f 1) > data.dat

exit;
