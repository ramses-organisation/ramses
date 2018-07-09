#!/bin/bash

../visu/amr2map -inp output_00002 -out data1.dat -typ 1 -fil ascii;

echo $(python smbh-bondi-comp.py) > data.dat

exit;
