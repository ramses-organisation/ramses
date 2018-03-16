#! /bin/bash

bin=$HOME'/bin'

a=$1
count='0000'$a
if [ $a -ge 10 ]; then
    count='000'$a
fi
if [ $a -ge 100 ]; then
    count='00'$a
fi
if [ $a -ge 1000 ]; then
    count='0'$a
fi
echo output_$count
aexp=`grep aexp output_$count/info_* | cut -d = -f 2`
echo 'aexp='$aexp
redshift=`echo $aexp|awk '{OFMT="%.2f"; print 1/$1-1}'`
echo 'z='$redshift
unit_l=`grep unit_l output_$count/info_* | cut -d = -f 2`
echo 'unit_l='$unit_l
boxlen1=`grep boxlen output_$count/info_* | cut -d = -f 2`
echo 'boxlen='$boxlen1
boxlen=`echo $boxlen1 $unit_l|awk '{OFMT="%.7e"; print $1*$2/3.08e21}'`
echo 'boxsize='$boxlen
lmax=`grep levelmax output_$count/info_* | cut -d = -f 2`
echo 'levelmax='$lmax

source params_$count.sh

r200=`echo $rvir $redshift $boxlen|awk '{print $1*(1+$2)/($3)}'`
echo 'r200c='$r200
rgal=`echo $r200|awk '{print $1/10}'`
hgal=`echo $rgal|awk '{print $1}'`

dir=$2

imsize=$r200
xmi=`echo $xc $imsize|awk '{print $1-$2}'`
xma=`echo $xc $imsize|awk '{print $1+$2}'`
ymi=`echo $yc $imsize|awk '{print $1-$2}'`
yma=`echo $yc $imsize|awk '{print $1+$2}'`
zmi=`echo $zc $imsize|awk '{print $1-$2}'`
zma=`echo $zc $imsize|awk '{print $1+$2}'`

part2map -inp output_$count -out dark_${count}_dir${dir}.map  -xmi $xmi -xma $xma -ymi $ymi -yma $yma -zmi $zmi -zma $zma -nx 1024 -ny 1024 -dir ${dir}
amr2map -inp output_$count -typ 1 -out dens_${count}_dir${dir}.map -xmi $xmi -xma $xma -ymi $ymi -yma $yma -zmi $zmi -zma $zma -nx 1024 -ny 1024 -lma ${lmax} -dir ${dir}
amr2map -inp output_$count -typ 0 -out level_${count}_dir${dir}.map -xmi $xmi -xma $xma -ymi $ymi -yma $yma -zmi $zmi -zma $zma -nx 1024 -ny 1024 -lma ${lmax} -dir ${dir} -act 2

#amr2map -inp output_$count -typ 5 -out temp_${count}_vir_dir${dir}.map -xmi $xmi -xma $xma -ymi $ymi -yma $yma -zmi $zmi -zma $zma -nx 1024 -ny 1024 -lma ${lmax} -dir ${dir}
#amr2map -inp output_$count -typ 6 -out met_${count}_vir_dir${dir}.map  -xmi $xmi -xma $xma -ymi $ymi -yma $yma -zmi $zmi -zma $zma -nx 1024 -ny 1024 -lma ${lmax} -dir ${dir}
#amr2map -inp output_$count -typ 7 -out blast_${count}_vir_dir${dir}.map  -xmi $xmi -xma $xma -ymi $ymi -yma $yma -zmi $zmi -zma $zma -nx 1024 -ny 1024 -lma ${lmax} -dir ${dir}
#part2map -inp output_$count -out star_${count}_vir_dir${dir}.map  -xmi $xmi -xma $xma -ymi $ymi -yma $yma -zmi $zmi -zma $zma -nx 1024 -ny 1024 -dir ${dir} -str true -den hop_star/hop$count.den
#part2map -inp output_$count -out star_${count}_vir_dir${dir}.map  -xmi $xmi -xma $xma -ymi $ymi -yma $yma -zmi $zmi -zma $zma -nx 1024 -ny 1024 -dir ${dir} -str true -age true




