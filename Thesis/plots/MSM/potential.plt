

set terminal cairo size 8cm,5cm standalone color colortext header \
   "\\renewcommand\\familydefault{\\sfdefault}\\usepackage{cmbright}"

set output 'potential_2010.tex'

set style line 2 lc rgb '#0060ad' lt 1 lw 5 pt 10 ps 0.5   #  blue
set style line 3 lc rgb '#dd181f' lt 1 lw 5 pt 5 ps .5   # red
set style line 4 lc rgb '#66A61E' lt 1 lw 5 pt 8 ps .5   #  green
set style line 6 lc rgb '#ff8800' lt 1 lw 5 pt 7 ps .5   #  orange
set style line 7 lc rgb '#9FAFDF' lt 1 lw 5 pt 8 ps 0.5   #royal blue
set style line 8 lc rgb '#1E90FF' lt 1 lw 5 pt 8 ps 0.5   #  slate blue
set style line 1 lc rgb '#F897C5' lt 1 lw 5 pt 8 ps 0.5   # taffy
set style line 5 lc rgb 'black' lt 1 lw 5 pt 5 ps 0.5   # black

set style fill transparent solid 0.25 # partial transparency
set style fill noborder # no separate top/bottom lines

set xlabel "$x / \\mathcal{L}$"

set xtics 0.2
set ytics 0.01

set xrange [0:1]
set ylabel "$U / \\epsilon$"
set ytics 1
set yrange [-1.2:4.2]



plot \
"../../data/MSM/potential_2010.dat" u 1:2  w l ls 3 lw 6 dt 1 lc 8 notitle,


pause -1
