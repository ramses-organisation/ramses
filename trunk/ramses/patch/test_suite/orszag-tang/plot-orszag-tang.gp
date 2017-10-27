reset

file1 = 'data1.dat'

set palette defined (0 "#000090",1 "#000fff",2 "#0090ff",3 "#0fffee",4 "#90ff70",5 "#ffee00",6 "#ff7000",7 "#ee0000",8 "#7f0000")

set term post enh color landscape size 9in,9in
set output 'orszag-tang.ps'

set size square
set view map
set xlabel 'Distance x'
set ylabel 'Distance y'
set cblabel 'Density'
unset key
splot file1 u 1:2:3 w image

unset output
