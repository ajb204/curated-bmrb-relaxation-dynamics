set term pdf
set border 3
set size square
set tics nomirror
set xlabel 'tc(fit) (ns)'
set ylabel 'tc(estimated) (ns)'
unset key
set output 'tc1.pdf'
set title 'R2'
plot 'test.out' u ($3*1E9):((1E-3/(1.0*1E6))*($2)*(1/($1*1.381E-23*6.022E23))*1E9) lc 1,x
set output 'tc2.pdf'
set title 'NOE'
plot 'test.out' u ($4*1E9):((1E-3/(1.0*1E6))*($2)*(1/($1*1.381E-23*6.022E23))*1E9) lc 2,x
set output 'tc3.pdf'
set title 'R1'
plot 'test.out' u ($5*1E9):((1E-3/(1.0*1E6))*($2)*(1/($1*1.381E-23*6.022E23))*1E9) lc 3,x

