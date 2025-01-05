

set terminal cairo size 11cm,7cm standalone color colortext header \
   "\\renewcommand\\familydefault{\\sfdefault}\\usepackage{cmbright}"

set output 'tao_MFPT.tex'

set style line 2 lc rgb '#0060ad' lt 1 lw 5 pt 10 ps 0.5   #  blue
set style line 3 lc rgb '#dd181f' lt 1 lw 5 pt 5 ps .5   # red
set style line 4 lc rgb '#66A61E' lt 1 lw 5 pt 8 ps .5   #  green
set style line 6 lc rgb '#ff8800' lt 1 lw 5 pt 7 ps .5   #  orange
set style line 7 lc rgb '#9FAFDF' lt 1 lw 5 pt 8 ps 0.5   #royal blue
set style line 8 lc rgb '#1E90FF' lt 1 lw 5 pt 8 ps 0.5   #  slate blue
set style line 1 lc rgb '#F897C5' lt 1 lw 5 pt 8 ps 0.5   # taffy
set style line 5 lc rgb 'black' lt 1 lw 5 pt 5 ps 0.5   # black


set xlabel "$\\tau / \\mathcal{T}$"
set ylabel "MFPT  $ / \\mathcal{T}$"

set key outside right


set xtics 0.002 
set ytics 0.1

set xrange[0:0.01]

plot "../../data/MSM/mom1_2010.dat" u ($1*0.00001):($3*$1*0.00001) w lp ls 2 title "$A \\rightarrow B$" ,\
"" u ($1*0.00001):($4*$1*0.00001) w lp ls 3  title "$A \\rightarrow C$",\
"" u ($1*0.00001):($5*$1*0.00001) w lp ls 4  title "$B \\rightarrow A$",\
"" u ($1*0.00001):($7*$1*0.00001) w lp ls 6  title "$B \\rightarrow C$",\
"" u ($1*0.00001):($8*$1*0.00001) w lp ls 7  title "$C \\rightarrow A$",\
"" u ($1*0.00001):($9*$1*0.00001) w lp ls 8  title "$C \\rightarrow B$",\

pause -1
