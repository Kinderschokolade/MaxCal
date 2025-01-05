
set terminal cairo size 11cm,15cm standalone color colortext header \
   "\\renewcommand\\familydefault{\\sfdefault}\\usepackage{cmbright}"


set output 'mom_1030.tex'

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

set origin 0,0.775
set size 1,0.225

set xtics
set format x ''

set xrange[0.0:9.1]
set yrange[10:750]

set ytics 100,100,601


unset xlabel

set ylabel "MFPT $\\mu /\\tau$" offset -0,0

set key at graph 1.,0.91 samplen 3  horizontal

set label '(a)' at screen 0.01,0.98
set label '(b)' at screen 0.01,0.77
set label '(c)' at screen 0.01,0.545
set label '(d)' at screen 0.01,0.320

plot \
"../../data/Frew/mom1_cts_1030.dat" u 1:3 w l ls 2 lw 5 notitle, \
"../../data/Frew/mom1_cts_1030.dat" u 1:4 w l ls 3 lw 5 notitle, \
"../../data/Frew/mom1_cts_1030.dat" u 1:5 w l ls 4 lw 5 notitle, \
"../../data/Frew/mom1_cts_1030.dat" u 1:7 w l ls 6 lw 5 notitle, \
"../../data/Frew/mom1_cts_1030.dat" u 1:8 w l ls 7 lw 5 notitle, \
"../../data/Frew/mom1_cts_1030.dat" u 1:9 w l ls 8 lw 5 notitle, \
"../../data/Frew/mom1_1030.dat" u 1:3 ls 2 lw 3 title "$A \\rightarrow B$",\
"../../data/Frew/mom1_1030.dat" u 1:5 ls 4 lw 3 title "$B \\rightarrow A$" ,\
"../../data/Frew/mom1_1030.dat" u 1:4 ls 3 lw 3 title "$A \\rightarrow C$" ,\
"../../data/Frew/mom1_1030.dat" u 1:8 ls 7 lw 3 title "$C \\rightarrow A$",\
"../../data/Frew/mom1_1030.dat" u 1:7 ls 6 lw 3 title "$B \\rightarrow C$",\
"../../data/Frew/mom1_1030.dat" u 1:9 ls 8 lw 3 title "$C \\rightarrow B$",\
"../../data/Frew/mom1_cts_5010.dat" u 1:3 w l ls 2 lw 3 dt 2 notitle, \
"../../data/Frew/mom1_cts_5010.dat" u 1:4 w l ls 3 lw 3 dt 2 notitle, \
"../../data/Frew/mom1_cts_5010.dat" u 1:5 w l ls 4 lw 3 dt 2 notitle, \
"../../data/Frew/mom1_cts_5010.dat" u 1:7 w l ls 6 lw 3 dt 2 notitle, \
"../../data/Frew/mom1_cts_5010.dat" u 1:8 w l ls 7 lw 3 dt 2 notitle, \
"../../data/Frew/mom1_cts_5010.dat" u 1:9 w l ls 8 lw 3 dt 2 notitle, 


unset key
set origin 0,0.55
set size 1,0.225
set tmargin 0.1

set yrange [10:650]
set ytics 100,100,601
set ylabel "STD $\\sigma /\\tau$" offset 0.2,0
plot \
"../../data/Frew/mom2_cts_1030.dat" u 1:3 w l ls 2 lw 5 notitle, \
"../../data/Frew/mom2_cts_1030.dat" u 1:4 w l ls 3 lw 5 notitle, \
"../../data/Frew/mom2_cts_1030.dat" u 1:5 w l ls 4 lw 5 notitle, \
"../../data/Frew/mom2_cts_1030.dat" u 1:7 w l ls 6 lw 5 notitle, \
"../../data/Frew/mom2_cts_1030.dat" u 1:8 w l ls 7 lw 5 notitle, \
"../../data/Frew/mom2_cts_1030.dat" u 1:9 w l ls 8 lw 5 notitle, \
"../../data/Frew/mom2_1030.dat" u 1:3 ls 2 lw 3 title "A \rightarrow B",\
"../../data/Frew/mom2_1030.dat" u 1:5 ls 4 lw 3 title "B \rightarrow A" ,\
"../../data/Frew/mom2_1030.dat" u 1:4 ls 3 lw 3 title "A \rightarrow C" ,\
"../../data/Frew/mom2_1030.dat" u 1:8 ls 7 lw 3 title "C \rightarrow A",\
"../../data/Frew/mom2_1030.dat" u 1:7 ls 6 lw 3 title "B \rightarrow C",\
"../../data/Frew/mom2_1030.dat" u 1:9 ls 8 lw 3 title "C \rightarrow B",\
"../../data/Frew/mom2_cts_5010.dat" u 1:3 w l ls 2 lw 3 dt 2 notitle, \
"../../data/Frew/mom2_cts_5010.dat" u 1:4 w l ls 3 lw 3 dt 2 notitle, \
"../../data/Frew/mom2_cts_5010.dat" u 1:5 w l ls 4 lw 3 dt 2 notitle, \
"../../data/Frew/mom2_cts_5010.dat" u 1:7 w l ls 6 lw 3 dt 2 notitle, \
"../../data/Frew/mom2_cts_5010.dat" u 1:8 w l ls 7 lw 3 dt 2 notitle, \
"../../data/Frew/mom2_cts_5010.dat" u 1:9 w l ls 8 lw 3 dt 2 notitle, 



set origin 0,0.325
set size 1,0.225

unset yrange
set yrange[1.5:3]
set ytics 1.7,0.2,2.8
set ylabel "skewness $\\kappa$" offset 1.1,0,0

plot \
"../../data/Frew/mom3_cts_1030.dat" u 1:3 w l ls 2 lw 5 notitle, \
"../../data/Frew/mom3_cts_1030.dat" u 1:4 w l ls 3 lw 5 notitle, \
"../../data/Frew/mom3_cts_1030.dat" u 1:5 w l ls 4 lw 5 notitle, \
"../../data/Frew/mom3_cts_1030.dat" u 1:7 w l ls 6 lw 5 notitle, \
"../../data/Frew/mom3_cts_1030.dat" u 1:8 w l ls 7 lw 5 notitle, \
"../../data/Frew/mom3_cts_1030.dat" u 1:9 w l ls 8 lw 5 notitle, \
"../../data/Frew/mom3_1030.dat" u 1:3 ls 2 lw 3 title "$A \\rightarrow B$",\
"../../data/Frew/mom3_1030.dat" u 1:5 ls 4 lw 3 title "$B \\rightarrow A$" ,\
"../../data/Frew/mom3_1030.dat" u 1:4 ls 3 lw 3 title "$A \\rightarrow C$" ,\
"../../data/Frew/mom3_1030.dat" u 1:8 ls 7 lw 3 title "$C \\rightarrow A$",\
"../../data/Frew/mom3_1030.dat" u 1:7 ls 6 lw 3 title "$B \\rightarrow C$",\
"../../data/Frew/mom3_1030.dat" u 1:9 ls 8 lw 3 title "$C \\rightarrow B$",\
"../../data/Frew/mom3_cts_5010.dat" u 1:3 w l ls 2 lw 3 dt 2 notitle, \
"../../data/Frew/mom3_cts_5010.dat" u 1:4 w l ls 3 lw 3 dt 2 notitle, \
"../../data/Frew/mom3_cts_5010.dat" u 1:5 w l ls 4 lw 3 dt 2 notitle, \
"../../data/Frew/mom3_cts_5010.dat" u 1:7 w l ls 6 lw 3 dt 2 notitle, \
"../../data/Frew/mom3_cts_5010.dat" u 1:8 w l ls 7 lw 3 dt 2 notitle, \
"../../data/Frew/mom3_cts_5010.dat" u 1:9 w l ls 8 lw 3 dt 2 notitle, 


set origin 0,0.1
set size 1,0.225

set xtics nomirror

set yrange[0:0.78]
set ytics 0.1,0.1,0.7
set xlabel "$f~/ ( \\epsilon / \\mathcal{L} )$"
set ylabel "occupation $\\Pi$" offset 0.0,0,0
set format x

set key at graph 1,0.92

plot \
"../../data/Frew/ps_cts_1030.dat" u 1:2 w l ls 1 lw 5 notitle, \
"../../data/Frew/ps_cts_1030.dat" u 1:3 w l ls 5 lw 5 notitle, \
"../../data/Frew/ps_cts_1030.dat" u 1:4 w l ls 9 lw 5 notitle, \
"../../data/Frew/P_1030.dat" u 1:2 ls 1 lw 3 title "$A$", \
"../../data/Frew/P_1030.dat" u 1:3 ls 5 lw 3 title "$B$", \
"../../data/Frew/P_1030.dat" u 1:4 ls 9 lw 3 title "$C$",\
"../../data/Frew/P_5010.dat" u 1:2 w l ls 1 lw 3 dt 2 notitle, \
"../../data/Frew/P_5010.dat" u 1:3 w l ls 5 lw 3 dt 2 notitle, \
"../../data/Frew/P_5010.dat" u 1:4 w l ls 9 lw 3 dt 2 notitle, \



#inset of skewness
set origin 0.33,0.45
set size 0.54,0.11
set bmargin 0.1
set tmargin 1
unset key
unset grid
set ytics 2,3,10
set yrange [1.9:*]
set format x ''
unset xlabel
unset ylabel
unset xtics

plot \
"../../data/Frew/mom3_1030.dat" u 1:3 ls 2 lw 5 notitle,\
"../../data/Frew/mom3_cts_1030.dat" u 1:3 w l ls 2 lw 5 notitle, \
"../../data/Frew/mom3_cts_5010.dat" u 1:3 w l ls 2 lw 3 dt 2 notitle, \


