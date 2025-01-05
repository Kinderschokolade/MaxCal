reset

set terminal cairo size 14cm,6cm standalone color colortext header \
   "\\renewcommand\\familydefault{\\sfdefault}\\usepackage{cmbright}"

set output "fpt_2010.tex"

set style line 2 lc rgb '#0060ad' lt 1 lw 5 pt 5 ps 0.5   #  blue
set style line 3 lc rgb '#dd181f' lt 1 lw 5 pt 7 ps .5   # red
set style line 4 lc rgb '#66A61E' lt 1 lw 5 pt 9 ps .5   #  green
set style line 6 lc rgb '#ff8800' lt 1 lw 5 pt 11 ps .5   #  orange
set style line 7 lc rgb '#9FAFDF' lt 1 lw 5 pt 13 ps 0.5   #royal blue
set style line 8 lc rgb '#1E90FF' lt 1 lw 5 pt 6 ps 0.5   #  slate blue
set style line 1 lc rgb '#F897C5' lt 1 lw 5 pt 8 ps 0.5   # taffy
set style line 5 lc rgb 'black' lt 1 lw 5 pt 10 ps 0.5   # black
set style line 9 lc rgb '#f1c40f' lt 1 lw 5 pt 12 ps 0.5   #yellow


set format x2 

set multiplot
set origin 0.0,0.0
set size 0.6,1

set ylabel "$p_{\\; \\textrm{FPT}}$" 
set xlabel "$t / \\tau$"

set ytics 0.002

set ylabel "$p_s$"
set yrange [0:0.008]

set logscale x

set key at screen 0.97,0.35

set xrange[1:1000]
set xtics 1,10,1000

plot "../../data/Urew/markovana_2010.dat" u 0:7 w lp ls 2 title "Simulation$: f=0~\\epsilon / \\mathcal{L}$" ,\
"../../data/Urew/markovana_2210.dat" u 0:7 w lp ls 3 title "Simulation$: f=2~\\epsilon / \\mathcal{L}$" ,\
"../../data/Urew/markovana_2910.dat" u 0:7 w lp ls 4 title "Simulation$: f=9~\\epsilon / \\mathcal{L}$" 

