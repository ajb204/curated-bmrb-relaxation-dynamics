set term pdf enh color solid
set size square
set key outside
set logscale 
set yrange[*:*]
set xlabel '{/Symbol t}_c(s)'
set ylabel 'Rate(s^{-1})'
set border 3
set xtics nomirror
set ytics nomirror
set title '500'
set output 'tmp/figNH_1.pdf'
plot 'tmp/outy.out' u 1:2 ti 'N1p'w li lw 5 lt 1,'tmp/outy.out' u 1:3 ti 'N1_z'w li lw 5 lt 2,'tmp/outy.out' u 1:4 ti 'H2_z'w li lw 5 lt 3
