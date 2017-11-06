#!/bin/bash
# This script takes all snapshots made with movie scheme of RAMSES, 
# produces png files with python/pyplot and converts them to
# a movie with ffmpeg
#
# Author: P. Biernacki
# biernack@physik.uzh.ch

# defining command-line arguments
dir=$1
kind1=$2
kind2=$3
step=$4
proj1=$5
proj2=$6
max=$7
# path to where your directories with simulations are
path=/zbox/data/wombat/
 
if [ "${dir:(-1)}" != "/" ]; then # check if dir was given with /
  dir=$dir/
fi

echo "Type: $kind, step: $step"
if [ "$#" -lt 7 ]; then
  printf "Usage: ./map2qmp4.sh <dir> <kind1> <kind2> <step> <proj_index1> <proj_index2> <max_iter> [<min1> <max1> <min2> <max2>], where:\n<kind>: dens, temp, metal <dir>\n"
else
  echo "Preparing pngs..."
  if [ "$#" -eq 11 ]; then
    eval "map2img.py -s $step -c bone -l True -p $proj1 -d ${path}${dir} -k $kind1 -F $max -m ${8} -M ${9}"
    eval "map2img.py -s $step -c Spectral_r -l True -p $proj1 -d ${path}${dir} -k $kind2 -F $max -m ${10} -M ${11}"
    eval "map2img.py -s $step -c bone -l True -p $proj2 -d ${path}${dir} -k $kind1 -F $max -m ${8} -M ${9}"
    eval "map2img.py -s $step -c Spectral_r -l True -p $proj2 -d ${path}${dir} -k $kind2 -F $max -m ${10} -M ${11}"
  else
    eval "map2img.py -s $step -c bone -l True -p $proj1 -d ${path}${dir} -k $kind1 -F $max"
    eval "map2img.py -s $step -c Spectral_r -l True -p $proj1 -d ${path}${dir} -k $kind2 -F $max"
    eval "map2img.py -s $step -c bone -l True -p $proj2 -d ${path}${dir} -k $kind1 -F $max"
    eval "map2img.py -s $step -c Spectral_r -l True -p $proj2 -d ${path}${dir} -k $kind2 -F $max"
  fi
  for j in `seq 1 $max`; do
		i=$(printf %05d $j)
		montage movie$proj1/pngs/${kind1}_$i.png movie$proj2/pngs/${kind1}_$i.png movie$proj1/pngs/${kind2}_$i.png movie$proj2/pngs/${kind2}_$i.png -geometry 620x640 out_$i.png
	done
  echo "Creating the movie: ${path}${dir}multi.mp4"
  eval "ffmpeg -i out_%05d.png -y -pix_fmt yuv420p -vcodec h264 -r 25 -qp 5 multi.mp4"
  echo "Cleaning up!"
  eval "rm ${path}${dir}movie${proj1}/pngs -r"
  eval "rm ${path}${dir}movie${proj2}/pngs -r"
  eval "rm ${path}${dir}out_*.png"
  eval "chmod a+r multi.mp4"
fi
