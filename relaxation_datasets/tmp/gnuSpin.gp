set term pdf enh color solid
set ticslevel 0
set output 'tmp/figNHgeo.1.pdf'
set xlabel 'x'
set ylabel 'y'
set zlabel 'z(A)'
unset colorbox
unset xtics
unset ytics
unset key
set view 80,73
set view equal xyz
set size 1.2,1.2
set bmargin at screen 0.0
set xrange[-1.530000:1.530000]
set yrange[-1.530000:1.530000]
set zrange[-1.530000:1.530000]
set palette defined(0"blue",1"red")
set label sprintf("NH") at graph 0.1,0.1,0.9
splot 'tmp/dipA.out' u 1:2:3:4  noti  w li lw 2 lc palette,'tmp/dipV.out' u 1:2:3:4:5:6:7 noti w vectors lw 2 lc palette,'tmp/csaA.out' i 0 u 1:2:3 noti w li lw 0.5 lc 2,'tmp/csaA.out' i 1 u 1:2:3 noti w li lw 0.5 lc 3,'tmp/vpts.out' i 0 u 1:2:3 noti w points lc 2 ps 1 pt 7,'' u 1:2:3:4 noti w labels,'tmp/vpts.out' i 1 u 1:2:3 noti w points lc 3 ps 1 pt 7,'' u 1:2:3:4 noti w labels
unset label
set output 'tmp/figNHgeo.2.pdf'
set label sprintf("NH map xy") at graph 0.1,1.05,0.9
set view 0,0
unset zlabel
splot 'tmp/dipA.out' u 1:2:3:4  noti  w li lw 2 lc palette,'tmp/dipV.out' u 1:2:3:4:5:6:7 noti w vectors lw 2 lc palette,'tmp/csaA.out' i 0 u 1:2:3 noti w li lw 0.5 lc 2,'tmp/csaA.out' i 1 u 1:2:3 noti w li lw 0.5 lc 3,'tmp/vpts.out' i 0 u 1:2:3 noti w points lc 2 ps 1 pt 7,'' u 1:2:3:4 noti w labels,'tmp/vpts.out' i 1 u 1:2:3 noti w points lc 3 ps 1 pt 7,'' u 1:2:3:4 noti w labels
unset label
set output 'tmp/figNHgeo.2a.pdf'
set label sprintf("NH map xz") at graph 0.1,0.9,1.05
set view 90,0
set zlabel 'Z(A)'
splot 'tmp/dipA.out' u 1:2:3:4  noti  w li lw 2 lc palette,'tmp/dipV.out' u 1:2:3:4:5:6:7 noti w vectors lw 2 lc palette,'tmp/csaA.out' i 0 u 1:2:3 noti w li lw 0.5 lc 2,'tmp/csaA.out' i 1 u 1:2:3 noti w li lw 0.5 lc 3,'tmp/vpts.out' i 0 u 1:2:3 noti w points lc 2 ps 1 pt 7,'' u 1:2:3:4 noti w labels,'tmp/vpts.out' i 1 u 1:2:3 noti w points lc 3 ps 1 pt 7,'' u 1:2:3:4 noti w labels
unset label
set output 'tmp/figNHgeo.3.pdf'
set label sprintf("NH dipolar") at graph 0.1,0.1,0.9
set view 80,73
splot 'tmp/dipA.out' u 1:2:3:4  noti  w li lw 2 lc palette,'tmp/dipV.out' u 1:2:3:4:5:6:7 noti w vectors lw 2 lc palette,'tmp/vpts.out' i 0 u 1:2:3 noti w points lc 2 ps 1 pt 7,'' u 1:2:3:4 noti w labels,'tmp/vpts.out' i 1 u 1:2:3 noti w points lc 3 ps 1 pt 7,'' u 1:2:3:4 noti w labels
unset label
set output 'tmp/figNHgeo.4.pdf'
set label sprintf("NH CSA") at graph 0.1,0.1,0.9
set view 80,73
splot 'tmp/csaA.out' i 0 u 1:2:3 noti w li lw 0.5 lc 2,'tmp/csaA.out' i 1 u 1:2:3 noti w li lw 0.5 lc 3,'tmp/vpts.out' i 0 u 1:2:3 noti w points lc 2 ps 1 pt 7,'' u 1:2:3:4 noti w labels,'tmp/vpts.out' i 1 u 1:2:3 noti w points lc 3 ps 1 pt 7,'' u 1:2:3:4 noti w labels
unset label
