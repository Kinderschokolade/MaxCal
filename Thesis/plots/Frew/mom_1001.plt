
set terminal cairo size 11cm,17cm standalone color colortext header \
   "\\renewcommand\\familydefault{\\sfdefault}\\usepackage{cmbright}"


set output 'mom_1001.tex'

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

set origin 0.02,0.8
set size 0.98,0.19

set xtics
set format x ''

set xrange[-10:10]
set yrange[10:90]

set ytics 20,20,90

unset xlabel


set ylabel "MFPT $\\mu/\\tau$" offset -0,0

set key at graph 0.8,0.9 samplen 3  horizontal

set label '(a)' at screen 0.01,0.98
set label '(b)' at screen 0.01,0.79
set label '(c)' at screen 0.01,0.60
set label '(d)' at screen 0.01,0.41
set label '(e)' at screen 0.01,0.22

plot \
"../../data/Frew/mom1_cts_1001.dat" u 1:4 w l ls 2 lw 5 notitle, \
"../../data/Frew/mom1_cts_1001.dat" u 1:5 w l ls 3 lw 5 notitle, \
"../../data/Frew/mom1_cts_1001.dat" u 1:6 w l ls 4 lw 5 notitle, \
"../../data/Frew/mom1_cts_1001.dat" u 1:8 w l ls 6 lw 5 notitle, \
"../../data/Frew/mom1_cts_1001.dat" u 1:9 w l ls 7 lw 5 notitle, \
"../../data/Frew/mom1_cts_1001.dat" u 1:10 w l ls 8 lw 5 notitle, \
"../../data/Frew/mom1_1001.dat" u 1:4 ls 2 lw 3 title "H $\\rightarrow$ I$",\
"../../data/Frew/mom1_1001.dat" u 1:6 ls 4 lw 3 title "I $\\rightarrow$ H$" ,\
"../../data/Frew/mom1_1001.dat" u 1:5 ls 3 lw 3 title "H $\\rightarrow$ E$" ,\
"../../data/Frew/mom1_1001.dat" u 1:9 ls 7 lw 3 title "E $\\rightarrow$ H",\
"../../data/Frew/mom1_1001.dat" u 1:8 ls 6 lw 3 title "I $\\rightarrow$ E",\
"../../data/Frew/mom1_1001.dat" u 1:10 ls 8 lw 3 title "E $\\rightarrow$ I",\
"../../data/Frew/mom1_1002.dat" u 1:4 ls 2 lw 3 notitle ,\
"../../data/Frew/mom1_1002.dat" u 1:6 ls 4 lw 3 notitle ,\
"../../data/Frew/mom1_1002.dat" u 1:5 ls 3 lw 3 notitle ,\
"../../data/Frew/mom1_1002.dat" u 1:9 ls 7 lw 3 notitle ,\
"../../data/Frew/mom1_1002.dat" u 1:8 ls 6 lw 3 notitle ,\
"../../data/Frew/mom1_1002.dat" u 1:10 ls 8 lw 3 notitle


unset key
set origin 0.02,0.61
set size 0.98,0.19
set tmargin 0.1

set yrange [10:70]
set ytics 20,20,90
set ylabel "STD $\\sigma /\\tau$" offset 0.2,0
plot \
"../../data/Frew/mom2_cts_1001.dat" u 1:4 w l ls 2 lw 5 notitle, \
"../../data/Frew/mom2_cts_1001.dat" u 1:5 w l ls 3 lw 5 notitle, \
"../../data/Frew/mom2_cts_1001.dat" u 1:6 w l ls 4 lw 5 notitle, \
"../../data/Frew/mom2_cts_1001.dat" u 1:8 w l ls 6 lw 5 notitle, \
"../../data/Frew/mom2_cts_1001.dat" u 1:9 w l ls 7 lw 5 notitle, \
"../../data/Frew/mom2_cts_1001.dat" u 1:10 w l ls 8 lw 5 notitle, \
"../../data/Frew/mom2_1001.dat" u 1:4 ls 2 lw 3 title "A \rightarrow B",\
"../../data/Frew/mom2_1001.dat" u 1:6 ls 4 lw 3 title "B \rightarrow A" ,\
"../../data/Frew/mom2_1001.dat" u 1:5 ls 3 lw 3 title "A \rightarrow C" ,\
"../../data/Frew/mom2_1001.dat" u 1:9 ls 7 lw 3 title "C \rightarrow A",\
"../../data/Frew/mom2_1001.dat" u 1:8 ls 6 lw 3 title "B \rightarrow C",\
"../../data/Frew/mom2_1001.dat" u 1:10 ls 8 lw 3 title "C \rightarrow B",\
"../../data/Frew/mom2_1002.dat" u 1:4 ls 2 lw 3 title "A \rightarrow B",\
"../../data/Frew/mom2_1002.dat" u 1:6 ls 4 lw 3 title "B \rightarrow A" ,\
"../../data/Frew/mom2_1002.dat" u 1:5 ls 3 lw 3 title "A \rightarrow C" ,\
"../../data/Frew/mom2_1002.dat" u 1:9 ls 7 lw 3 title "C \rightarrow A",\
"../../data/Frew/mom2_1002.dat" u 1:8 ls 6 lw 3 title "B \rightarrow C",\
"../../data/Frew/mom2_1002.dat" u 1:10 ls 8 lw 3 title "C \rightarrow B",\



set origin 0.02,0.42
set size 0.98,0.19

unset yrange
set yrange[1.8:2.7]
set ytics 1.8,0.2,2.8
set ylabel "skewness $\\kappa$" offset 1.1,0,0

plot \
"../../data/Frew/mom3_cts_1001.dat" u 1:4 w l ls 2 lw 5 notitle, \
"../../data/Frew/mom3_cts_1001.dat" u 1:5 w l ls 3 lw 5 notitle, \
"../../data/Frew/mom3_cts_1001.dat" u 1:6 w l ls 4 lw 5 notitle, \
"../../data/Frew/mom3_cts_1001.dat" u 1:8 w l ls 6 lw 5 notitle, \
"../../data/Frew/mom3_cts_1001.dat" u 1:9 w l ls 7 lw 5 notitle, \
"../../data/Frew/mom3_cts_1001.dat" u 1:10 w l ls 8 lw 5 notitle, \
"../../data/Frew/mom3_1001.dat" u 1:4 ls 2 lw 3 title "$A \\rightarrow B$",\
"../../data/Frew/mom3_1001.dat" u 1:6 ls 4 lw 3 title "$B \\rightarrow A$" ,\
"../../data/Frew/mom3_1001.dat" u 1:5 ls 3 lw 3 title "$A \\rightarrow C$" ,\
"../../data/Frew/mom3_1001.dat" u 1:9 ls 7 lw 3 title "$C \\rightarrow A$",\
"../../data/Frew/mom3_1001.dat" u 1:8 ls 6 lw 3 title "$B \\rightarrow C$",\
"../../data/Frew/mom3_1001.dat" u 1:10 ls 8 lw 3 title "$C \\rightarrow B$",\
"../../data/Frew/mom3_1002.dat" u 1:4 ls 2 lw 3 title "$A \\rightarrow B$",\
"../../data/Frew/mom3_1002.dat" u 1:6 ls 4 lw 3 title "$B \\rightarrow A$" ,\
"../../data/Frew/mom3_1002.dat" u 1:5 ls 3 lw 3 title "$A \\rightarrow C$" ,\
"../../data/Frew/mom3_1002.dat" u 1:9 ls 7 lw 3 title "$C \\rightarrow A$",\
"../../data/Frew/mom3_1002.dat" u 1:8 ls 6 lw 3 title "$B \\rightarrow C$",\
"../../data/Frew/mom3_1002.dat" u 1:10 ls 8 lw 3 title "$C \\rightarrow B$",\


set origin 0.02,0.23
set size 0.98,0.19


set yrange[0:0.25]
set ytics 0.05,0.05,0.22
set ylabel "occupation $\\Pi$" offset 0.0,0,0

set key at graph 1,0.92

plot \
"../../data/Frew/ps_cts_1001.dat" u 1:3 w l ls 1 lw 5 notitle, \
"../../data/Frew/ps_cts_1001.dat" u 1:4 w l ls 5 lw 5 notitle, \
"../../data/Frew/ps_cts_1001.dat" u 1:5 w l ls 9 lw 5 notitle, \
"../../data/Frew/P_1001.dat" u 1:3 ls 1 lw 3 title "H", \
"../../data/Frew/P_1001.dat" u 1:5 ls 9 lw 3 title "E",\
"../../data/Frew/P_1001.dat" u 1:4 ls 5 lw 3 title "I", \
"../../data/Frew/P_1002.dat" u 1:3 ls 1 lw 3 notitle, \
"../../data/Frew/P_1002.dat" u 1:5 ls 9 lw 3 notitle ,\
"../../data/Frew/P_1002.dat" u 1:4 ls 5 lw 3 notitle



set origin 0.02,0.06
set size 0.98,0.17

set xtics nomirror

set yrange[100:1400]
set ytics 200,300,1300
set xlabel "$f~/ ( \\epsilon / nm )$"
set ylabel "$t_i /$ fs" offset 0.0,0.0
set format x

set key at graph 1,0.92
plot \
"../../data/Frew/timescale_cts_1001.dat" u 1:3 w l ls 1 lw 5 title "$t_1$", \
"../../data/Frew/timescale_cts_1001.dat" u 1:4 w l ls 5 lw 5 title "$t_2$", \
"../../data/Frew/lagv_1001.dat" u 1:3 ls 1 notitle, \
"../../data/Frew/lagv_1001.dat" u 1:4 ls 5 notitle, \
"../../data/Frew/lagv_1002.dat" u 1:3 ls 1 notitle, \
"../../data/Frew/lagv_1002.dat" u 1:4 ls 5 notitle, \




##inset of skewness
#set origin 0.33,0.45
#set size 0.54,0.11
#set bmargin 0.1
#set tmargin 1
#unset key
#unset grid
#set ytics 2,3,10
#set yrange [1.9:*]
#set format x ''
#unset xlabel
#unset ylabel
#unset xtics
#
#plot \
#"../../data/Frew/mom3_1001.dat" u 1:3 ls 2 lw 5 notitle,\
#"../../data/Frew/mom3_cts_1001.dat" u 1:3 w l ls 2 lw 5 notitle, \
#

