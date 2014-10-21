#!/bin/bash
# This script takes all snapshots made with movie scheme of RAMSES, 
# produces png files with python/pyplot and converts them to
# a movie with ffmpeg
#
# Author: P. Biernacki
# biernack@physik.uzh.ch

# defining command-line arguments
dir=$1
kind=$2
step=$3
proj=$4
# path to where your directories with simulations are
path='/zbox/user/biernack/'
 
if [ "${dir:(-1)}" != "/" ]; then # check if dir was given with /
  dir=$dir/
fi

echo "Type: $kind, step: $step"
if [ "$#" -ne 4 ]; then
  printf "Usage: ./map2mp4.sh <dir> <kind> <step> <proj_index>, where:\n<kind>: dens, temp, metal\n <dir>"
else
  echo "Preparing pngs..."
  eval "../py/map2img.py -s $step -c RdBu_r -l True -p $proj -d ${path}${dir} -k $kind -m 1e4 -M 1e8"
  echo "Creating the movie: ${path}${dir}$kind$proj.mp4"
  eval "ffmpeg -loglevel quiet -i ${path}${dir}movie${proj}/pngs/${kind}_%05d.png -y -vcodec libx264 -vpre slow -crf 10 -r 25 ${path}${dir}$kind$proj.mp4"
  echo "Cleaning up!"
  eval "rm ${path}${dir}movie${proj}/pngs/$kind*.png"
  eval "rmdir ${path}${dir}movie${proj}/pngs"
  eval "chmod a+r ${path}${dir}${kind}${proj}.mp4"
fi
