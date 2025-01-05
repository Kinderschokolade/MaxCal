reset

set terminal cairo size 15cm,8cm standalone color colortext header \
   "\\renewcommand\\familydefault{\\sfdefault}\\usepackage{cmbright}"

set output "Shist.tex"

set style line 2 lc rgb '#0060ad' lt 1 lw 5 pt 5 ps 0.5   #  blue
set style line 3 lc rgb '#dd181f' lt 1 lw 5 pt 7 ps .5   # red
set style line 4 lc rgb '#66A61E' lt 1 lw 5 pt 9 ps .5   #  green
set style line 6 lc rgb '#ff8800' lt 1 lw 5 pt 11 ps .5   #  orange
set style line 7 lc rgb '#9FAFDF' lt 1 lw 5 pt 13 ps 0.5   #royal blue
set style line 8 lc rgb '#1E90FF' lt 1 lw 5 pt 6 ps 0.5   #  slate blue
set style line 1 lc rgb '#F897C5' lt 1 lw 5 pt 8 ps 0.5   # taffy
set style line 5 lc rgb 'black' lt 1 lw 5 pt 10 ps 0.5   # black
set style line 9 lc rgb '#f1c40f' lt 1 lw 5 pt 12 ps 0.5   #yellow


set multiplot

set origin 0.0,0.0
set size 0.5,1


set label "(a)" at screen 0.01,0.97
set label "(b)" at screen 0.51,0.97

set xlabel "$\\Delta S~/  (\\epsilon / K)$"
set ylabel "Probability"

set xrange[-0.2:1.5]
set yrange[0:0.1]

set arrow nohead ls 5 dt 1  lc  rgb 'yellow'from 0.0103,0 to 0.0103,0.1
set arrow nohead ls 5 dt 1  lc  rgb 'yellow'from 0.05366,0 to 0.05366,0.1
set arrow nohead ls 5 dt 1 lc  rgb 'yellow' from 0.2175,0 to 0.2175,0.1
set arrow nohead ls 5 dt 1 lc  rgb 'yellow' from 0.66234,0 to 0.66234,0.1
set arrow nohead ls 5 dt 1 lc  rgb 'yellow' from 1.325,0 to 1.325,0.1

set arrow nohead ls 5 dt 2 from 0.011,0 to 0.011,0.1
set arrow nohead ls 5 dt 2 from 0.0571,0 to 0.0571,0.1
set arrow nohead ls 5 dt 2 from 0.2189,0 to 0.2189,0.1
set arrow nohead ls 5 dt 2 from 0.6625,0 to 0.6625,0.1
set arrow nohead ls 5 dt 2 from 1.3189,0 to 1.3189,0.1



plot "../../data/Complementary/Shist_2015_15_16.dat" u 1:($2/14606100) w lp ls 2 notitle ,\
"../../data/Complementary/Shist_2015_15_17.dat" u 1:($2/11602900) w lp ls 3 notitle ,\
"../../data/Complementary/Shist_2015_15_18.dat" u 1:($2/8261230) w lp ls 4 notitle ,\
"../../data/Complementary/Shist_2015_15_19.dat" u 1:($2/4562600) w lp ls 6 notitle ,\
"../../data/Complementary/Shist_2015_15_20.dat" u 1:($2/1984370) w lp ls 7 notitle ,\


unset arrow

set origin 0.5,0.0
set size 0.5,1


set xrange[0:3.25]
set yrange[0:0.1]
set ytics 0.02


set arrow nohead ls 5 dt 1 lc  rgb 'yellow' from 0.621,0 to 0.621,0.1
set arrow nohead ls 5 dt 2 from 0.6258,0 to 0.6258,0.1

set arrow nohead ls 5 dt 1  lc  rgb 'yellow'from 1.479,0 to 1.479,0.1
set arrow nohead ls 5 dt 2 from 1.484,0 to 1.484,0.1

set arrow nohead ls 5 dt 1  lc  rgb 'yellow'from 2.164,0 to 2.164,0.1
set arrow nohead ls 5 dt 2 from 2.176,0 to 2.176,0.1

set arrow nohead ls 5 dt 1  lc  rgb 'yellow'from 2.447,0 to 2.447,0.1
set arrow nohead ls 5 dt 2 from 2.510,0 to 2.510,0.1




plot "../../data/Complementary/Shist_5099_125_126.dat" u (-$1):($2/5710490) w lp ls 2 notitle ,\
"../../data/Complementary/Shist_5099_125_127.dat" u (-$1):($2/1209580) w lp ls 3 notitle ,\
"../../data/Complementary/Shist_5099_125_128.dat" u (-$1):($2/235032) w lp ls 4 notitle ,\
"../../data/Complementary/Shist_5099_125_129.dat" u (-$1):($2/58452) w lp ls 6 notitle ,\
















