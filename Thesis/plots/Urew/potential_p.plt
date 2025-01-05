set terminal cairo size 11cm,9cm standalone color colortext header \
   "\\renewcommand\\familydefault{\\sfdefault}\\usepackage{cmbright}"

set style line 2 lc rgb '#0060ad' lt 1 lw 5 pt 10 ps 0.5   #  blue
set style line 3 lc rgb '#dd181f' lt 1 lw 5 pt 5 ps .5   # red
set style line 4 lc rgb '#66A61E' lt 1 lw 5 pt 8 ps .5   #  green
set style line 6 lc rgb '#ff8800' lt 1 lw 5 pt 7 ps .5   #  orange
set style line 7 lc rgb '#9FAFDF' lt 1 lw 5 pt 8 ps 0.5   #royal blue
set style line 8 lc rgb '#1E90FF' lt 1 lw 5 pt 8 ps 0.5   #  slate blue
set style line 1 lc rgb '#F897C5' lt 1 lw 5 pt 8 ps 0.5   # taffy
set style line 5 lc rgb 'black' lt 1 lw 5 pt 5 ps 0.5   # black

set output "potential_p.tex"


set xlabel "$x / \\mathcal{L}$"
set format x2 


set xtics 0.2
set ytics 0.01

set xrange [0:1]
set y2range [-2.5:5.5]

set multiplot

set origin 0.0,0.0
set size 1,0.6

set lmargin 8
set rmargin 1
set tmargin 0.1
set bmargin 4

#set grid x

set key top left 
set ylabel "$p_s$"
set yrange [0:0.08]
set ytics 0.0,0.02,0.07

plot \
"../../data/Urew/ps_2010.dat" u ($0/30+1./120):1 every 2 w l ls 3 lw 5 dt 1 lc 7 title "$f_{\\textup{eq}}$",\
"../../data/Urew/ps_2910.dat" u ($0/30+1./120.):1 every 2 w l ls 7 lw 5 dt 1 lc 6 smooth csplines title "$f_{\\textup{neq}}$",\
"../../data/Urew/p_2910_2010.dat" u ($0/30+1./120.):1 every 2 pt 9 lw 2 ps .9 lc 7 title "$f_{\\textup{neq}} \\rightarrow f_{\\textup{eq}}$",\
"../../data/Urew/p_2010_2910.dat" u ($0/30+1./120.):1 every 2 pt 7 lw 2 ps .9 lc 6 title "$f_{\\textup{eq}} \\rightarrow f_{\\textup{neq}}$" ,\
#"../../data/Urew/ps_2910.dat" u ($0/60+1./120.):1  title "$f_{\\textup{neq}}$",\

set origin 0.0,0.6
set size 1,0.4
set bmargin 0.1
set tmargin 1
unset key
set ylabel "$U / \\epsilon$"
set ytics 2,2,6
set yrange [-0.:7.5]
set format x ''
unset xlabel

set label '$A$' at 0.24,2.6
set label '$B$' at 0.565,1.4
set label '$C$' at 0.9,2.6
set label '(b)' at screen 0.05,0.58
set label '(a)' at screen 0.05,0.95

set arrow from 0.52,3.3 to 0.55,2.5 head lc rgb '#66A61E' lw 3  size 2,0.5
set arrow from 0.645,3.2 to 0.62,2.5 head lc rgb '#ff8800' lw 3  size 2,0.5 

#set object 1 circle at first 0.23,4.0 size 1 fc rgb '#ff8800'
#set object 2 circle at first 0.27,4.0 size 1 fc rgb '#66A61E'
# does not work for some reasonset ytics 2,2,6

plot \
"../../data/Urew/traj_p1.dat" u 1:2 smooth csplines dt 1 lc rgb '#ff8800' lw 4,\
"../../data/Urew/traj_p2.dat" u 1:2 smooth csplines dt 1 lc rgb '#ff8800' lw 4,\
"../../data/Urew/traj_p3.dat" u 1:2 smooth csplines dt 1 lc rgb '#66A61E' lw 4,\
"../../data/Urew/potential_2010.dat" u ($1):($2+2) ls 5  lw 6 w l 



