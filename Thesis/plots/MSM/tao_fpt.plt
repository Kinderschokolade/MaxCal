

set terminal cairo size 13cm,7cm standalone color colortext header \
   "\\renewcommand\\familydefault{\\sfdefault}\\usepackage{cmbright}"

set output 'tao_fpt.tex'

set style line 2 lc rgb '#0060ad' lt 1 lw 5 pt 5 ps 0.5   #  blue
set style line 3 lc rgb '#dd181f' lt 1 lw 5 pt 7 ps .5   # red
set style line 4 lc rgb '#66A61E' lt 1 lw 5 pt 9 ps .5   #  green
set style line 6 lc rgb '#ff8800' lt 1 lw 5 pt 11 ps .5   #  orange
set style line 7 lc rgb '#9FAFDF' lt 1 lw 5 pt 13 ps 0.5   #royal blue
set style line 8 lc rgb '#1E90FF' lt 1 lw 5 pt 6 ps 0.5   #  slate blue
set style line 1 lc rgb '#F897C5' lt 1 lw 5 pt 8 ps 0.5   # taffy
set style line 5 lc rgb 'black' lt 1 lw 5 pt 10 ps 0.5   # black
set style line 9 lc rgb '#f1c40f' lt 1 lw 5 pt 12 ps 0.5   #yellow

set logscale x

set xlabel "$t  / \\mathcal{T} $"
set ylabel "$p_{\\; \\textrm{FPT}} $"

set ytics 0.001
set xrange[0.001:5]
set yrange [0:0.0055]

set key outside right

plot "../../data/MSM/hist_2010_50.dat" u ($0*50*0.00001):($4*4) w lp ls 2 title "$\\tau = 0.0005\\; \\mathcal{T}$",\
"../../data/MSM/hist_2010_100.dat" u ($0*100*0.00001):($4*2) w lp ls 3  title "$\\tau = 0.001\\; \\mathcal{T}$",\
"../../data/MSM/hist_2010_200.dat" u ($0*200*0.00001):4 w lp ls 4  title "$\\tau = 0.002\\; \\mathcal{T}$",\
"../../data/MSM/hist_2010_500.dat" u ($0*500*0.00001):($4/2.5) w lp ls 6  title "$\\tau = 0.005\\; \\mathcal{T}$",\
"../../data/MSM/hist_2010_1000.dat" u ($0*1000*0.00001):($4/5) w lp ls 7 title "$\\tau = 0.01\\; \\mathcal{T}$",\

pause -1
