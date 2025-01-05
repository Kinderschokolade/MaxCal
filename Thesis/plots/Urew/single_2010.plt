reset

set terminal cairo size 15cm,9cm standalone color colortext header \
   "\\renewcommand\\familydefault{\\sfdefault}\\usepackage{cmbright}"

set output "single_2010.tex"

set style line 2 lc rgb '#0060ad' lt 1 lw 5 pt 5 ps 0.5   #  blue
set style line 3 lc rgb '#dd181f' lt 1 lw 5 pt 7 ps .5   # red
set style line 4 lc rgb '#66A61E' lt 1 lw 5 pt 9 ps .5   #  green
set style line 6 lc rgb '#ff8800' lt 1 lw 5 pt 11 ps .5   #  orange
set style line 7 lc rgb '#9FAFDF' lt 1 lw 5 pt 13 ps 0.5   #royal blue
set style line 8 lc rgb '#1E90FF' lt 1 lw 5 pt 6 ps 0.5   #  slate blue
set style line 1 lc rgb '#F897C5' lt 1 lw 5 pt 8 ps 0.5   # taffy
set style line 5 lc rgb 'black' lt 1 lw 5 pt 10 ps 0.5   # black
set style line 9 lc rgb '#f1c40f' lt 1 lw 5 pt 12 ps 0.5   #yellow


set xlabel "$x~/\\mathcal{L}$"
set format x2 


set xtics 0.2
set ytics 0.01

set xrange [0:1]

set multiplot

set origin 0.0,0.0
set size 0.5,0.6

set lmargin 8
set rmargin 1
set tmargin 0.1
set bmargin 4


set key top left 
set ylabel "$p_s$"
set yrange [0:0.08]
set ytics 0.0,0.02,0.079


set arrow from 0.6,5 to 0.7,5 lw 5 lc rgb '#dd181f'
set label 'f'  at 0.635,4.45 textcolor rgb '#dd181f'


plot \
"../../data/Urew/p_2910_2010.dat" u ($0/30+1./120.):1 every 2   ls 2 notitle,\
"../../data/Urew/p_2010_2910.dat" u ($0/30+1./120.):1 every 2   ls 3 notitle  ,\
"../../data/Urew/ps_2010.dat" u ($0/30+1./120):1 every 2 w l ls 2 lw 5  notitle ,\
"../../data/Urew/ps_2910.dat" u ($0/30+1./120.):1 every 2 w l ls 3 lw 5  smooth csplines notitle ,\




set origin 0.0,0.6
set size 0.5,0.4
set bmargin 0.1
set tmargin 1
set ylabel "$U~/\\epsilon$"
set ytics 2
set yrange [-1:6]
set format x ''
unset xlabel

set label '$A$' at 0.24,2.6
set label '$B$' at 0.565,1.4
set label '$C$' at 0.9,2.6
set label '(b)' at screen 0.01,0.58
set label '(a)' at screen 0.01,0.95
set label '(c)' at screen 0.51,0.79

unset key

plot \
"../../data/Urew/potential_2010.dat" u ($1):($2+1) ls 2  lw 5 w l notitle  ,\

set origin 0.5,0.14
set size 0.48,0.7

unset format

set logscale x
set xrange [1:1000]

set ytics 0,0.005,00.051
set xtics 1,10,1000

set yrange[0:0.017]


set ylabel "$p_{\\; \\textrm{FPT}}$" 
set xlabel "$t / \\tau$"


set key default
set key at screen 0.97,0.96

plot "../../data/Urew/hist_2910_2010.dat" u 0:6 ls 2 title "Reweight$: f=9~\\epsilon / \\mathcal{L} \\rightarrow f=0~\\epsilon / \\mathcal{L}$" ,\
"../../data/Urew/hist_2010_2910.dat" u 0:6 ls 3 title "Reweight$: f=0~\\epsilon / \\mathcal{L} \\rightarrow f=9~\\epsilon / \\mathcal{L}$"   ,\
"../../data/Urew/hist_2010_2910.dat" u 0:15 w l  ls 3 lw 5 title "Simulation$: f=0~\\epsilon / \\mathcal{L}$" ,\
"../../data/Urew/hist_2910_2010.dat" u 0:15 w l  ls 2 lw 5 title "Simulation$: f=9~\\epsilon / \\mathcal{L}$" ,\


#plot "../../data/Urew/hist_2910_2010.dat" u 0:6 ls 2 ps 1 notitle ,\
#"../../data/Urew/hist_2010_2910.dat" u 0:6 ls 3 ps 1 notitle  ,\
#"../../data/Urew/hist_2910_2010.dat" u 0:15 w l  ls 2 lw 5 notitle ,\
#"../../data/Urew/hist_2010_2910.dat" u 0:15 w l  ls 3 lw 5 notitle   



