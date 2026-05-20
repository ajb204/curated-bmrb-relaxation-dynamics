set size square
set border 3
set tics nomirror
set term pdf
set key outside
set title 'Model free'
set xlabel 'Rate(exp)'
set ylabel 'Rate(calc)'
set output 'figGlobFreeAx1.pdf'
plot 'outy.globModelFreeAx.out' i 0 u 1:2 ti 'R1'w points pt 7 ps 0.2,x noti lc -1 lw 0.05,'outy.globModelFreeAx.out' i 1 u 1:2 ti 'R2'w points pt 7 ps 0.2,x noti lc -1 lw 0.05,'outy.globModelFreeAx.out' i 2 u 1:2 ti 'NOE'w points pt 7 ps 0.2,x noti lc -1 lw 0.05
set border 11
set tics nomirror
set xlabel 'residue'
set ylabel 'S^2'
set y2label ' taue'
set y2tics
set logscale y2
set format y2 "10^{%L}"
set output 'figGlobFreeAx2.pdf'
plot 'outy.globModelFreeParsAx.out' u 1:2 ti 'S2'w points pt 7 ps 0.2 lc 1,'outy.globModelFreeParsAx.out' u 1:3 ti 'taue'w points pt 7 ps 0.2 lc 2 axes x1y2,'outy.globModelFreeParsAx.out' u 1:4 ti 'tauc' w li lc 2 axes x1y2 
set border 3
unset y2tics
set output 'figGlobFreeAx4.pdf'
set xlabel 'residue'
set ylabel 'error (s-1)'
unset y2label
plot 'outy.globModelFreeParsAx.out' u 1:(abs($6-$8)+abs($9-$11)+abs($12-$14)) noti w points pt 7 ps 0.2
set border 3
unset y2tics
set output 'figGlobFreeAx3.pdf'
set xlabel 'S^2'
set ylabel 'taue ns'
set logscale y
set format y "10^{%L}"
set cblabel 'residue'
plot 'outy.globModelFreeParsAx.out' u 2:($3*1E9):1 noti w points pt 7 ps 0.2 lc palette
