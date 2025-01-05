set terminal cairo size 7cm,5cm standalone color colortext header \
  "\\renewcommand\\familydefault{\\sfdefault}\\usepackage{cmbright}"

set output 'potential_traj.tex'

set style line 2 lc rgb '#0060ad' lt 1 lw 5 pt 10 ps 0.5   #  blue
set style line 3 lc rgb '#dd181f' lt 1 lw 5 pt 5 ps .5   # red
set style line 4 lc rgb '#66A61E' lt 1 lw 5 pt 8 ps .5   #  green
set style line 6 lc rgb '#ff8800' lt 1 lw 5 pt 7 ps .5   #  orange
set style line 7 lc rgb '#9FAFDF' lt 1 lw 5 pt 8 ps 0.5   #royal blue
set style line 8 lc rgb '#1E90FF' lt 1 lw 5 pt 8 ps 0.5   #  slate blue
set style line 1 lc rgb '#F897C5' lt 1 lw 5 pt 8 ps 0.5   # taffy
set style line 5 lc rgb 'black' lt 1 lw 5 pt 5 ps 0.5   # black


set xlabel "$x / \\mathcal{L}$"

set xtics 0.2
set ytics 0.01

set xrange [0:1]
set ylabel "$U / \\epsilon$"
set ytics 2
set yrange [-1.2:5.7]

set arrow from 0.52,1.3 to 0.55,0.5 head lc rgb '#ff8800'  lw 3  size 2,0.5
set arrow from 0.645,1.2 to 0.62,0.5 head lc rgb '#66A61E' lw 3  size 2,0.5 

set arrow from 0.25,0 to 0.25,3 nohead ls 5 dt 2 lw 2
set arrow from 0.58,-1 to 0.58,3 nohead ls 5 dt 2 lw 2

set label "$x_0$" at 0.24,3.3
set label "$x_T$" at 0.57,3.3

plot \
"../../data/Jaynes/traj_p1.dat" u 1:($2-2) smooth csplines ls 4 lw 4 notitle,\
"../../data/Jaynes/traj_p2.dat" u 1:($2-2) smooth csplines ls 4 lw 4 notitle,\
"../../data/Jaynes/traj_p3.dat" u 1:($2-2) smooth csplines ls 6 lw 4 notitle,\
"../../data/MSM/potential_2010.dat" u ($1):($2) ls 5  lw 6 w l title "potential"


