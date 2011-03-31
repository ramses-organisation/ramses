count='00420'
xc=0.499055
yc=0.500313
zc=0.500152
uc=-85.
vc=144.
wc=157.
jx=0.402
jy=-0.644
jz=0.657

#part2prof -inp output_$count -out snap_zoom0.prof -xce $xc -yce $yc -zce $zc -uce $uc -vce $vc -wce $wc  -rma 0.02 -nra 50
amr2prof -inp output_$count -out snap_zoom0.prof -xce $xc -yce $yc -zce $zc -uce $uc -vce $vc -wce $wc  -rma 0.02 -nra 50

#part2prof -inp output_$count -out snap_zoom1.prof -xce $xc -yce $yc -zce $zc -uce $uc -vce $vc -wce $wc  -rma 0.004 -nra 50
amr2prof -inp output_$count -out snap_zoom1.prof -xce $xc -yce $yc -zce $zc -uce $uc -vce $vc -wce $wc  -rma 0.004 -nra 50

#part2prof -inp output_$count -out snap_zoom2.prof -xce $xc -yce $yc -zce $zc -uce $uc -vce $vc -wce $wc  -rma 0.0008 -nra 50
amr2prof -inp output_$count -out snap_zoom2.prof -xce $xc -yce $yc -zce $zc -uce $uc -vce $vc -wce $wc  -rma 0.0008 -nra 50

#part2prof -inp output_$count -out snap_zoom3.prof -xce $xc -yce $yc -zce $zc -uce $uc -vce $vc -wce $wc  -rma 0.00016 -nra 50
amr2prof -inp output_$count -out snap_zoom3.prof -xce $xc -yce $yc -zce $zc -uce $uc -vce $vc -wce $wc  -rma 0.00016 -nra 50

#cat snap_zoom3.prof.dark snap_zoom2.prof.dark snap_zoom1.prof.dark snap_zoom0.prof.dark > snap_$count.prof.dark
#cat snap_zoom3.prof.star snap_zoom2.prof.star snap_zoom1.prof.star snap_zoom0.prof.star > snap_$count.prof.star
cat snap_zoom3.prof.gas  snap_zoom2.prof.gas  snap_zoom1.prof.gas  snap_zoom0.prof.gas  > snap_$count.prof.gas

#part2cylprof -inp output_$count -out snap_$count.cyl -xce $xc -yce $yc -zce $zc -uce $uc -vce $vc -wce $wc -jx $jx -jy $jy -jz $jz -rma 0.003 -nra 200 -hma 0.0002





