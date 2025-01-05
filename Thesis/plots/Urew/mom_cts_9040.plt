reset
#

set terminal cairo size 3.8,5.0 standalone color colortext header \
   "\\renewcommand\\familydefault{\\sfdefault}\\usepackage{cmbright}"



set output 'mom_9040.tex'

set style line 1 lc rgb '#0060ad' lt 1 lw 5 pt 8 ps 0.9   #  blue
set style line 2 lc rgb '#dd181f' lt 1 lw 5 pt 9 ps .9   # --- red
set style line 3 lc rgb '#66A61E' lt 1 lw 5 pt 6 ps .9   # ——— green
set style line 4 lc rgb '#ff8800' lt 1 lw 5 pt 4 ps .9   # ——— orange
set style line 5 lc rgb '#9FAFDF' lt 1 lw 5 pt 5 ps 0.9   # --- royal blue
set style line 6 lc rgb '#1E90FF' lt 1 lw 5 pt 7 ps 0.9   # --- slate blue
set style line 7 lc rgb '#39FF14' lt 1 lw 5 pt 12 ps 0.9   # --- neon green 
set style line 8 lc rgb '#CA1F7B' lt 1 lw 5 pt 13 ps 0.9   # --- magenta 


set multiplot layout 2,1

set xtics 25
set xrange[0.0:100]

set grid 

set ylabel "Mean $[\\tau]$" 

set yrange [10:2300]

plot \
"../data/mom1_cts_9040.dat" u 1:3 w l ls 1 lw 5 notitle, \
"../data/mom1_cts_9040.dat" u 1:5 w l ls 2 lw 5 notitle, \
"../data/mom1_cts_9040.dat" u 1:6 w l ls 3 lw 5 notitle, \
"../data/mom1_cts_9040.dat" u 1:8 w l ls 4 lw 5 notitle, \
"../data/mom1_cts_9040.dat" u 1:11 w l ls 5 lw 5 notitle, \
"../data/mom1_cts_9040.dat" u 1:13 w l ls 6 lw 5 notitle, \
"../data/mom1_cts_9040.dat" u 1:14 w l ls 7 lw 5 notitle, \
"../data/mom1_cts_9040.dat" u 1:16 w l ls 8 lw 5 notitle, \
"../data/mom1_9040.dat" u 1:3 ls 1 lw 5 title "A-B",\
"../data/mom1_9040.dat" u 1:5 ls 2 lw 5 title "A-D" ,\
"../data/mom1_9040.dat" u 1:6 ls 3 lw 5 title "B-A" ,\
"../data/mom1_9040.dat" u 1:8 ls 4 lw 5 title "B-C",\
"../data/mom1_9040.dat" u 1:11 ls 5 lw 5 title "C-B",\
"../data/mom1_9040.dat" u 1:13 ls 6 lw 5 title "C-D",\
"../data/mom1_9040.dat" u 1:14 ls 7 lw 5 title "D-A",\
"../data/mom1_9040.dat" u 1:16 ls 8 lw 5 title "D-C"


set xtics 25
set xrange[100:250]

set grid 

set ylabel "Mean $[\\tau]$" 
set xlabel  "$f$" 

set yrange [0:140]

plot \
"../data/mom1_cts_9040.dat" u 1:3 w l ls 1 lw 5 notitle, \
"../data/mom1_cts_9040.dat" u 1:5 w l ls 2 lw 5 notitle, \
"../data/mom1_cts_9040.dat" u 1:6 w l ls 3 lw 5 notitle, \
"../data/mom1_cts_9040.dat" u 1:8 w l ls 4 lw 5 notitle, \
"../data/mom1_cts_9040.dat" u 1:11 w l ls 5 lw 5 notitle, \
"../data/mom1_cts_9040.dat" u 1:13 w l ls 6 lw 5 notitle, \
"../data/mom1_cts_9040.dat" u 1:14 w l ls 7 lw 5 notitle, \
"../data/mom1_cts_9040.dat" u 1:16 w l ls 8 lw 5 notitle, \
"../data/mom1_9040.dat" u 1:3 ls 1 lw 5 notitle ,\
"../data/mom1_9040.dat" u 1:5 ls 2 lw 5 notitle ,\
"../data/mom1_9040.dat" u 1:6 ls 3 lw 5 notitle ,\
"../data/mom1_9040.dat" u 1:8 ls 4 lw 5 notitle ,\
"../data/mom1_9040.dat" u 1:11 ls 5 lw 5 notitle,\
"../data/mom1_9040.dat" u 1:13 ls 6 lw 5 notitle,\
"../data/mom1_9040.dat" u 1:14 ls 7 lw 5 notitle,\
"../data/mom1_9040.dat" u 1:16 ls 8 lw 5 notitle


pause -1
