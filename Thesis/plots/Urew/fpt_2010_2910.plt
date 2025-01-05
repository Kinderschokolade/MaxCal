
set terminal cairo size 11cm,7cm standalone color colortext header \
   "\\renewcommand\\familydefault{\\sfdefault}\\usepackage{cmbright}"


set style line 2 lc rgb '#0060ad' lt 1 lw 5 pt 10 ps 0.5   #  blue
set style line 3 lc rgb '#dd181f' lt 1 lw 5 pt 5 ps .5   # red
set style line 4 lc rgb '#66A61E' lt 1 lw 5 pt 8 ps .5   #  green
set style line 6 lc rgb '#ff8800' lt 1 lw 5 pt 7 ps .5   #  orange
set style line 7 lc rgb '#9FAFDF' lt 1 lw 5 pt 8 ps 0.5   #royal blue
set style line 8 lc rgb '#1E90FF' lt 1 lw 5 pt 8 ps 0.5   #  slate blue
set style line 1 lc rgb '#F897C5' lt 1 lw 5 pt 8 ps 0.5   # taffy
set style line 5 lc rgb 'black' lt 1 lw 5 pt 5 ps 0.5   # black



set output 'fpt_2010_2910.tex'

#set multiplot layout 1,2
#set origin 0,0
#set size 1,1

set logscale x
set xrange [1:1000]

set ytics 0,0.005,00.021
set xtics 1,10,1000

set yrange[0:0.021]


set ylabel "$p_{\\; \\textrm{FPT}}$" 
set xlabel "$t / \\tau$"


#set title "$f=0$"
#set title "equilibrium"

plot "../../data/Urew/hist_2910_2010.dat" u 0:6 pt 7 lc 6 ps .9 title "Reweight$: f=9~\\epsilon / \\mathcal{L} \\rightarrow f=0~\\epsilon / \\mathcal{L}$" ,\
"../../data/Urew/hist_2010_2910.dat" u 0:6 pt 6 lc 7 ps .9 title "Reweight$: f=0~\\epsilon / \\mathcal{L} \\rightarrow f=9~\\epsilon / \\mathcal{L}$"  ,\
"../../data/Urew/hist_2910_2010.dat" u 0:15 w l dt 1 lc 6 lw 5 title "Simulation$: f=0~\\epsilon / \\mathcal{L}$" ,\
"../../data/Urew/hist_2010_2910.dat" u 0:15 w l dt 1 lc 7 lw 5 title "Simulation$: f=9~\\epsilon / \\mathcal{L}$"  ,\


