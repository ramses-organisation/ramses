reset

file1 = 'data4plot.dat'

set term post enh color landscape size 8in,7in
set output 'smbh-bondi.ps'

# set size square
set view map
set xrange [0:60]
set yrange [1e-3:3e-3]
set format y "%.0s x 10^{%S}"
set ytics (0.001,0.002,0.003)
set logscale y
set xlabel 'Time [kyr]' font ",20"
set ylabel 'dot(M_{Bondi}) [Msol/yr]' font ",20"
unset key
plot file1 u 1:2 w lp lw 3 lt rgb "#1B62A5"

unset output
