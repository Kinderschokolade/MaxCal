reset

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

set output "potential_equ.tex"

set xlabel "$x~/\\mathcal{L}$"
set format x2 


set xtics 0.2
set ytics 0.01

set xrange [0.15:1]

set multiplot

set origin 0.0,0.0
set size 1,0.6

set lmargin 8
set rmargin 1
set tmargin 0.1
set bmargin 4

set grid x

set key top left 
set ylabel "$p_s$"
set yrange [0:0.15]
set ytics 0.0,0.03,0.14

plot \
"../../data/Urew/ps_3003.dat" u ($0/30+1./120):1 every 2 w l ls 3 lw 5 dt 1 lc 7 title "$U_1$",\
"../../data/Urew/ps_3903.dat" u ($0/30+1./120.):1 every 2 w l ls 7 lw 5 dt 1 lc 6 smooth csplines title "$U_2$",\
"../../data/Urew/p_3903_3003.dat" u ($0/30+1./120.):1 every 2 pt 9 lw 2 ps .9 lc 7 title "   $U_2 \\rightarrow U_1$",\
"../../data/Urew/p_3003_3903.dat" u ($0/30+1./120.):1 every 2 pt 7 lw 2 ps .9 lc 6 title "   $U_1 \\rightarrow U_2$" ,\

set origin 0.0,0.6
set size 1,0.4
set bmargin 0.1
set tmargin 1
set ylabel "$U~/\\epsilon$"
set ytics 2
set yrange [-5.5:5]
set format x ''
unset xlabel

set label '$A$' at 0.24,2.6
set label '$B$' at 0.565,1.4
set label '$C$' at 0.9,2.6
set label '(b)' at screen 0.05,0.58
set label '(a)' at screen 0.05,0.95


plot \
"../../data/Urew/potential_3003.dat" u ($1):($2+1) lc 7  lw 5 w l notitle  ,\
"../../data/Urew/potential_3003.dat" u ($1):($2-9*$1+3) lc 6  lw 5 w l notitle



