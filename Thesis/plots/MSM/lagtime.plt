

set terminal cairo size 8cm,7cm standalone color colortext header \
   "\\renewcommand\\familydefault{\\sfdefault}\\usepackage{cmbright}"

set output 'lagtime_2010.tex'

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

set yrange[0.085:0.13]
set ytics 0.01
set xtics 0.002
set xrange [0:0.01]

set ylabel "$t_i / \\mathcal{T}$"
set xlabel "$\\tau / \\mathcal{T}$"



set arrow from 0.002,0.08500 to 0.002,0.1205 lw 3 nohead

plot "../../data/MSM/lagtime_2010.dat" u ($1*0.00001):($2*0.00001) ls 2  w lp title "$t_1$" ,\
"" u ($1*0.00001):($3*0.00001) ls 3  w lp title "$t_2$" ,\


pause -1
