#!/bin/bash

# Extract time from log
s=$(grep Fine log | tail -n 1);
a=( $s );
echo ${a[4]} > time.dat;

# Extract data from log
l1=$(grep -n "===" log | tail -n 2 | head -n 1 | cut -d ':' -f1);
l2=$(grep -n "===" log | tail -n 1 | head -n 1 | cut -d ':' -f1);

line1=`expr $l1 + 2`;
line2=`expr $l2 - 1`;

sed -n "${line1},${line2}p" log > data.dat;
