
set terminal cairo size 11cm,11cm standalone color colortext header \
   "\\renewcommand\\familydefault{\\sfdefault}\\usepackage{cmbright}"

set output 'mom_2010.tex'

set style line 2 lc rgb '#0060ad' lt 1 lw 5 pt 10 ps 0.5   #  blue
set style line 3 lc rgb '#dd181f' lt 1 lw 5 pt 5 ps .5   # red
set style line 4 lc rgb '#66A61E' lt 1 lw 5 pt 8 ps .5   #  green
set style line 6 lc rgb '#ff8800' lt 1 lw 5 pt 7 ps .5   #  orange
set style line 7 lc rgb '#9FAFDF' lt 1 lw 5 pt 8 ps 0.5   #royal blue
set style line 8 lc rgb '#1E90FF' lt 1 lw 5 pt 8 ps 0.5   #  slate blue
set style line 1 lc rgb '#F897C5' lt 1 lw 5 pt 8 ps 0.5   # taffy
set style line 5 lc rgb 'black' lt 1 lw 5 pt 5 ps 0.5   # black


set tmargin 0.5
set bmargin 0.2
set lmargin 7
set rmargin 2
unset xtics

set multiplot 

set origin 0,0.7
set size 1,0.3

set xtics
set format x ''

set xrange[0.0:9.1]
set yrange[10:465]

set ytics 100,100,401


unset xlabel
#set grid 

set ylabel "Mean $ / \\tau$" offset -0.0,0

set key at graph 1.,0.94 samplen 3  horizontal


set label '(a)' at screen 0.01,0.985
set label '(b)' at screen 0.01,0.69
set label '(c)' at screen 0.01,0.39

plot \
"../../data/Urew/mom1_2010-x.dat" u 0:4 w l ls 2 lw 5 notitle, \
"../../data/Urew/mom1_2010-x.dat" u 0:5 w l ls 3 lw 5 notitle, \
"../../data/Urew/mom1_2010-x.dat" u 0:6 w l ls 4 lw 5 notitle, \
"../../data/Urew/mom1_2010-x.dat" u 0:8 w l ls 6 lw 5 notitle, \
"../../data/Urew/mom1_2010-x.dat" u 0:9 w l ls 7 lw 5 notitle, \
"../../data/Urew/mom1_2010-x.dat" u 0:10 w l ls 8 lw 5 notitle, \
"../../data/Urew/mom1_2010.dat" u 0:3 ls 2 lw 5 title "A-B",\
"../../data/Urew/mom1_2010.dat" u 0:5 ls 4 lw 5 title "B-A" ,\
"../../data/Urew/mom1_2010.dat" u 0:4 ls 3 lw 5 title "A-C" ,\
"../../data/Urew/mom1_2010.dat" u 0:8 ls 7 lw 5 title "C-A",\
"../../data/Urew/mom1_2010.dat" u 0:7 ls 6 lw 5 title "B-C",\
"../../data/Urew/mom1_2010.dat" u 0:9 ls 8 lw 5 title "C-B"


unset key
set origin 0,0.4
set size 1,0.3
set tmargin 0.1

set yrange [10:400]
set ytics 100,100,301
set ylabel "Variance $/ \\tau^2$" offset -0.0,0
plot \
"../../data/Urew/mom2_2010-x.dat" u 0:4 w l ls 2 lw 5 notitle, \
"../../data/Urew/mom2_2010-x.dat" u 0:5 w l ls 3 lw 5 notitle, \
"../../data/Urew/mom2_2010-x.dat" u 0:6 w l ls 4 lw 5 notitle, \
"../../data/Urew/mom2_2010-x.dat" u 0:8 w l ls 6 lw 5 notitle, \
"../../data/Urew/mom2_2010-x.dat" u 0:9 w l ls 7 lw 5 notitle, \
"../../data/Urew/mom2_2010-x.dat" u 0:10 w l ls 8 lw 5 notitle, \
"../../data/Urew/mom2_2010.dat" u 0:3 ls 2 lw 5 title "A-B",\
"../../data/Urew/mom2_2010.dat" u 0:5 ls 4 lw 5 title "B-A" ,\
"../../data/Urew/mom2_2010.dat" u 0:4 ls 3 lw 5 title "A-C" ,\
"../../data/Urew/mom2_2010.dat" u 0:8 ls 7 lw 5 title "C-A",\
"../../data/Urew/mom2_2010.dat" u 0:7 ls 6 lw 5 title "B-C",\
"../../data/Urew/mom2_2010.dat" u 0:9 ls 8 lw 5 title "C-B"


set origin 0,0.1
set size 1,0.3

set xtics nomirror


set yrange[1.7:2.5]
set ytics 1.7,0.2,2.41
set xlabel "$f / ( \\epsilon / \\mathcal{L} )$"
set ylabel "Skewness $/ \\tau^3$" offset 0.4,0,0
set format x

plot \
"../../data/Urew/mom3_2010-x.dat" u 0:4 w l ls 2 lw 5 notitle, \
"../../data/Urew/mom3_2010-x.dat" u 0:5 w l ls 3 lw 5 notitle, \
"../../data/Urew/mom3_2010-x.dat" u 0:6 w l ls 4 lw 5 notitle, \
"../../data/Urew/mom3_2010-x.dat" u 0:8 w l ls 6 lw 5 notitle, \
"../../data/Urew/mom3_2010-x.dat" u 0:9 w l ls 7 lw 5 notitle, \
"../../data/Urew/mom3_2010-x.dat" u 0:10 w l ls 8 lw 5 notitle, \
"../../data/Urew/mom3_2010.dat" u 0:3 ls 2 lw 5 title "$A \\rightarrow B$",\
"../../data/Urew/mom3_2010.dat" u 0:5 ls 4 lw 5 title "$B \\rightarrow A$" ,\
"../../data/Urew/mom3_2010.dat" u 0:4 ls 3 lw 5 title "$A \\rightarrow C$" ,\
"../../data/Urew/mom3_2010.dat" u 0:8 ls 7 lw 5 title "$C \\rightarrow A$",\
"../../data/Urew/mom3_2010.dat" u 0:7 ls 6 lw 5 title "$B \\rightarrow C$",\
"../../data/Urew/mom3_2010.dat" u 0:9 ls 8 lw 5 title "$C \\rightarrow B$"

set origin 0.35,0.27
set size 0.6,0.15
set bmargin 0.1
set tmargin 1
unset key
#unset grid
set ytics 2,1,4.9
set yrange [1.9:4.9]
set format x ''
unset xlabel
unset ylabel

plot \
"../../data/Urew/mom3_2010.dat" u 0:3 ls 2 lw 5 notitle,\
"../../data/Urew/mom3_2010-x.dat" u 0:4 w l ls 2 lw 5 notitle, \

