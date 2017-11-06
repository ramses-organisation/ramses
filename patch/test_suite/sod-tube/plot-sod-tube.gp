reset

scale_d = 1.0
scale_l = 1.0
scale_v = 1.0

lwidth = 2

file1 = 'data.dat'
file2 = 'sod-tube-ana.dat'

lmin = int(system(sprintf("grep levelmin sod-tube.nml | cut -d '=' -f2")))
lmax = int(system(sprintf("grep levelmax sod-tube.nml | cut -d '=' -f2")))

set term post enh color portrait
set output 'sod-tube.ps'

set xlabel 'Distance x (cm)'

nxpanel=1
nypanel=3
panlgap=0.08

marginx1=0.10
marginx2=0.90
marginy1=0.07
marginy2=0.99

dxpanel=(marginx2-marginx1-(nxpanel-1)*panlgap)/nxpanel
dypanel=(marginy2-marginy1-(nypanel-1)*panlgap)/nypanel

set multiplot
unset key

set lmargin at screen marginx1
set rmargin at screen marginx2
set bmargin at screen (marginy2-dypanel)
set tmargin at screen marginy2

set y2label 'AMR level' offset -2.0
set ytics nomirror
set y2tics
set autoscale  y
set y2tics lmin,1,lmax

set ylabel 'Density (g/cm3)' offset 2.0
plot file1 u 2:3 lc 0 lt 6, file2 u 2:4 w l lw lwidth lt 1, file1 u 2:1 w l lt 0 lc 0 axes x1y2

unset y2label
unset y2tics

set bmargin at screen (marginy1+dypanel+panlgap)
set tmargin at screen (marginy2-dypanel-panlgap)

set ylabel 'Velocity (cm/s)' offset 2.0
plot file1 u 2:4 lc 0 lt 6, file2 u 2:3 w l lw lwidth lt 1

set bmargin at screen marginy1
set tmargin at screen (marginy1+dypanel)

set ylabel 'Pressure (g/cm/s2)' offset 2.0
plot file1 u 2:5 lc 0 lt 6, file2 u 2:5 w l lw lwidth lt 1

unset multiplot

unset output
