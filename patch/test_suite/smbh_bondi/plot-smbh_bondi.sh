#!/bin/bash

../visu/amr2map -inp output_00002 -out data1.dat -typ 1 -fil ascii;

echo $(md5sum output_00002/sink_00002.csv | cut -d ' ' -f 1) > data.dat

exit;
