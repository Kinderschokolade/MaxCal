reset

set terminal cairo size 8cm,6cm standalone color colortext header \
   "\\renewcommand\\familydefault{\\sfdefault}\\usepackage{cmbright}"

set output "traj.tex"

set style line 2 lc rgb '#0060ad' lt 1 lw 5 pt 5 ps 0.5   #  blue
set style line 3 lc rgb '#dd181f' lt 1 lw 5 pt 7 ps .5   # red
set style line 4 lc rgb '#66A61E' lt 1 lw 5 pt 9 ps .5   #  green
set style line 6 lc rgb '#ff8800' lt 1 lw 5 pt 11 ps .5   #  orange
set style line 7 lc rgb '#9FAFDF' lt 1 lw 5 pt 13 ps 0.5   #royal blue
set style line 8 lc rgb '#1E90FF' lt 1 lw 5 pt 6 ps 0.5   #  slate blue
set style line 1 lc rgb '#F897C5' lt 1 lw 5 pt 8 ps 0.5   # taffy
set style line 5 lc rgb 'black' lt 1 lw 5 pt 10 ps 0.5   # black
set style line 9 lc rgb '#f1c40f' lt 1 lw 5 pt 12 ps 0.5   #yellow


set xlabel "$x / \\mathcal{L}$"
set ylabel "$y / \\mathcal{L}$"

set xtics 0.85,0.05,1.01
set ytics 0.35,0.05,0.52

set xrange[0.84:1.01]
set yrange[0.33:0.52]

set object 1 rect from  0.8,0.4 to 0.9,0.5 lw 3 dt 2 back
set object 2 rect from 0.9,0.3 to 1,0.4 lw 3 dt 2 back

set object 3 rect from  0.8333,0.43333 to 0.8666,0.4666 lw 3 back
set object 4 rect from 0.9333,0.33333 to 0.9666,0.3666 lw 3 back



set label "$x_0$" at  0.86,0.48 textcolor rgb '#dd181f' front 
set label "$x_T$" at 0.955,0.355 textcolor rgb '#dd181f' front

set label "$x_0$" at  0.863,0.425 textcolor rgb '#0060ad' front 
set label "$x_T$" at  0.939,0.345 textcolor rgb '#0060ad' front 

plot "../../data/Complementary/path9_128_99.dat" w lp ls 2 ps 0.1 notitle ,\
"../../data/Complementary/path10_128_99.dat" w lp ls 3 ps 0.1 notitle


