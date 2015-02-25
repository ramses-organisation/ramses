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
max=$5
# path to where your directories with simulations are
path=$ZUSER/
 
if [ "${dir:(-1)}" != "/" ]; then # check if dir was given with /
  dir=$dir/
fi

echo "Type: $kind, step: $step"
if [ "$#" -lt 4 ]; then
  printf "Usage: ./map2mp4.sh <dir> <kind> <step> <proj_index> [<max_iter>], where:\n<kind>: dens, temp, metal\n <dir>"
else
	if [ "$#" -eq 4 ]; then
		max=-1
	fi
  echo "Preparing pngs..."
  eval "~/scripts/map2img.py -s $step -c bone -l True -p $proj -d ${path}${dir} -k $kind -F $max" # -m 1e5 -M 1e8"
  #eval "~/scripts/map2img.py -s $step -c Spectral_r -l True -p $proj -d ${path}${dir} -k $kind -F $max  -m 1e5 -M 1e8"
  echo "Creating the movie: ${path}${dir}$kind$proj.mp4"
  eval "~/scripts/ffmpeg -loglevel quiet -i ${path}${dir}movie${proj}/pngs/${kind}_%05d.png -y -vcodec h264 -pix_fmt yuv420p  -r 25 -qp 5 ${path}${dir}$kind$proj.mp4"
  echo "Cleaning up!"
  eval "rm ${path}${dir}movie${proj}/pngs -r"
  eval "chmod a+r ${path}${dir}${kind}${proj}.mp4"
fi
