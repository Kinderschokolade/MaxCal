
set terminal cairo size 15cm,7cm standalone color colortext header \
   "\\renewcommand\\familydefault{\\sfdefault}\\usepackage{cmbright}"

set output 'fpt_Sloc_2010_2910.tex'

set style line 2 lc rgb '#0060ad' lt 1 lw 5 pt 10 ps 0.5   #  blue
set style line 3 lc rgb '#dd181f' lt 1 lw 5 pt 5 ps .5   # red
set style line 4 lc rgb '#66A61E' lt 1 lw 5 pt 8 ps .5   #  green
set style line 6 lc rgb '#ff8800' lt 1 lw 5 pt 7 ps .5   #  orange
set style line 7 lc rgb '#9FAFDF' lt 1 lw 5 pt 8 ps 0.5   #royal blue
set style line 8 lc rgb '#1E90FF' lt 1 lw 5 pt 8 ps 0.5   #  slate blue
set style line 1 lc rgb '#F897C5' lt 1 lw 5 pt 8 ps 0.5   # taffy
set style line 5 lc rgb 'black' lt 1 lw 5 pt 5 ps 0.5   # black

set multiplot layout 1,2

set logscale x


set key outside bottom horizontal

set ylabel "$p_{\\; \\textrm{FPT}}$" 
set xlabel "t / $\\tau$" offset 0,0
set yrange[0:0.02]
set xrange[1:3000]
set ytics 0.005
set label "a)" at graph -0.3,1

plot "../../data/Jaynes/hist_Sloc_2910_2010.dat" u 0:6 pt 7 lc 6 ps .5  title "Reweight$: f=9~\\epsilon / \\mathcal{L} \\rightarrow f=0~\\epsilon / \\mathcal{L}$" ,\
"../../data/Jaynes/hist_Sloc_2010_2910.dat" u 0:6 pt 7 lc 7 ps .5  title "Reweight$: f=0~\\epsilon / \\mathcal{L} \\rightarrow f=9~\\epsilon / \\mathcal{L}$"  ,\
"../../data/Jaynes/hist_Sloc_2910_2010.dat" u 0:15 w l dt 1 lw 5 lc 6  title "Simulation$: f=0~\\epsilon / \\mathcal{L}$" ,\
"../../data/Jaynes/hist_Sloc_2010_2910.dat" u 0:15 w l dt 1 lw 5 lc 7  title "Simulation$: f=9~\\epsilon / \\mathcal{L}$"  ,\

unset label
unset logscale
unset yrange
unset xrange
set ytics 0.02
set xtics 0.2
set ylabel "$p_s$"
set xlabel "$x / \\mathcal{L}$" 
set label "b)" at graph -0.3,1

plot "../../data/Jaynes/p_Sloc_2910_2010.dat" u ($0/60.+1./120.):1 pt 7 lc 6 ps .5 notitle,\
"../../data/Jaynes/p_Sloc_2010_2910.dat" u  ($0/60.+1./120.):1 pt 7 lc 7 ps .5 notitle,\
"../../data/Jaynes/ps_2010.dat" u  ($0/60.+1./120.):1 w l dt 1 lw 5 lc 6  notitle ,\
"../../data/Jaynes/ps_2910.dat" u  ($0/60.+1./120.):1 w l dt 1 lw 5 lc 7  notitle


pause -1


