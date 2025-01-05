
set terminal cairo size 11cm,15cm standalone color colortext header \
   "\\renewcommand\\familydefault{\\sfdefault}\\usepackage{cmbright}"


set output 'mom_9040.tex'

set tmargin 0.5
set bmargin 0.2
set lmargin 7
set rmargin 2
unset xtics


set style line 2 lc rgb '#0060ad' lt 1 lw 5 pt 5 ps 0.5   #  blue
set style line 3 lc rgb '#dd181f' lt 1 lw 5 pt 7 ps .5   # red
set style line 4 lc rgb '#66A61E' lt 1 lw 5 pt 9 ps .5   #  green
set style line 6 lc rgb '#ff8800' lt 1 lw 5 pt 11 ps .5   #  orange
set style line 7 lc rgb '#9FAFDF' lt 1 lw 5 pt 13 ps 0.5   #royal blue
set style line 8 lc rgb '#1E90FF' lt 1 lw 5 pt 6 ps 0.5   #  slate blue
set style line 1 lc rgb '#F897C5' lt 1 lw 5 pt 8 ps 0.5   # taffy
set style line 5 lc rgb 'black' lt 1 lw 5 pt 10 ps 0.5   # black
set style line 9 lc rgb '#f1c40f' lt 1 lw 5 pt 12 ps 0.5   #yellow


set multiplot 

set origin 0.02,0.775
set size 1,0.225

set xtics
set format x ''

set xrange[0.0:255]
set yrange[9:2200]

#set ytics 0,100,401

set logscale y

unset xlabel

set ylabel "MFPT $\\mu / \\tau$" offset 1.2,0

set key at graph 1.,0.94 samplen 3  horizontal

set label '(a)' at screen 0.01,0.98
set label '(b)' at screen 0.01,0.77
set label '(c)' at screen 0.01,0.545
set label '(d)' at screen 0.01,0.320

plot \
"../../data/Urew/mom1_cts_9040.dat" u 1:3 w l ls 3 lw 5 notitle, \
"../../data/Urew/mom1_cts_9040.dat" u 1:6 w l ls 6 lw 5 notitle, \
"../../data/Urew/mom1_cts_9040.dat" u 1:13 w l ls 2 lw 5 notitle, \
"../../data/Urew/mom1_cts_9040.dat" u 1:16 w l ls 7 lw 5 notitle, \
"../../data/Urew/mom1_9040.dat" u 1:3 ls 3 lw 3 title "$A \\rightarrow B$",\
"../../data/Urew/mom1_9040.dat" u 1:6 ls 6 lw 3 title "$B \\rightarrow A$" ,\
"../../data/Urew/mom1_9040.dat" u 1:13 ls 2 lw 3 title "$C \\rightarrow D$" ,\
"../../data/Urew/mom1_9040.dat" u 1:16 ls 7 lw 3 title "$D \\rightarrow C$"

unset key
set origin 0.02,0.55
set size 1,0.225
set tmargin 0.1

set yrange [9:2200]
#set ytics 0,100,390
set ylabel "STD $\\sigma /\\tau$" offset 1.2,0
plot \
"../../data/Urew/mom2_cts_9040.dat" u 1:3 w l ls 3 lw 5 notitle, \
"../../data/Urew/mom2_cts_9040.dat" u 1:6 w l ls 6 lw 5 notitle, \
"../../data/Urew/mom2_cts_9040.dat" u 1:13 w l ls 2 lw 5 notitle, \
"../../data/Urew/mom2_cts_9040.dat" u 1:16 w l ls 7 lw 5 notitle, \
"../../data/Urew/mom2_9040.dat" u 1:3 ls 3 lw 3 title "A-B",\
"../../data/Urew/mom2_9040.dat" u 1:6 ls 6 lw 3 title "B-A" ,\
"../../data/Urew/mom2_9040.dat" u 1:13 ls 2 lw 3 title "B-C" ,\
"../../data/Urew/mom2_9040.dat" u 1:16 ls 7 lw 3 title "C-B",\


unset logscale 
set origin 0.02,0.325
set size 1,0.225

unset yrange
set yrange[1.:5]
set ytics 2,1,4
set ylabel "skewness $\\kappa$" offset -0.3,0

plot \
"../../data/Urew/mom3_cts_9040.dat" u 1:3 w l ls 3 lw 5 notitle, \
"../../data/Urew/mom3_cts_9040.dat" u 1:6 w l ls 6 lw 5 notitle, \
"../../data/Urew/mom3_cts_9040.dat" u 1:13 w l ls 2 lw 5 notitle, \
"../../data/Urew/mom3_cts_9040.dat" u 1:16 w l ls 7 lw 5 notitle, \
"../../data/Urew/mom3_9040.dat" u 1:3 ls 3 lw 3 title "$A \\rightarrow B$",\
"../../data/Urew/mom3_9040.dat" u 1:6 ls 6 lw 3 title "$B \\rightarrow A$" ,\
"../../data/Urew/mom3_9040.dat" u 1:13 ls 2 lw 3 title "$A \\rightarrow C$" ,\
"../../data/Urew/mom3_9040.dat" u 1:16 ls 7 lw 3 title "$C \\rightarrow B$"


set origin 0.02,0.1
set size 1,0.225

set xtics nomirror


set yrange[0:0.31]
set ytics 0.05,0.05,0.25
set xlabel "$f~/ ( \\epsilon / \\mathcal{L} )$"
set ylabel "occupation $\\Pi$" offset 0.0,0,0
set format x

set key at graph 0.9,0.9

plot \
"../../data/Urew/ps_cts_9040.dat" u 1:2 w l ls 1 lw 5 notitle, \
"../../data/Urew/ps_cts_9040.dat" u 1:3 w l ls 5 lw 5 notitle, \
"../../data/Urew/ps_cts_9040.dat" u 1:4 w l ls 9 lw 5 notitle, \
"../../data/Urew/ps_cts_9040.dat" u 1:5 w l ls 4 lw 5 notitle, \
"../../data/Urew/P_9040.dat" u 1:2 ls 1 lw 3 title "$A$", \
"../../data/Urew/P_9040.dat" u 1:3 ls 5 lw 3 title "$B$", \
"../../data/Urew/P_9040.dat" u 1:4 ls 9 lw 3 title "$C$",\
"../../data/Urew/P_9040.dat" u 1:5 ls 4 lw 3 title "$D$"


