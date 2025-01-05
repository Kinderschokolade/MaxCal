
set terminal cairo size 15cm,7cm standalone color colortext header \
   "\\renewcommand\\familydefault{\\sfdefault}\\usepackage{cmbright}"

set output 'cktest_2010.tex'

set style line 2 lc rgb '#0060ad' lt 1 lw 5 pt 10 ps 0.5   #  blue
set style line 3 lc rgb '#dd181f' lt 1 lw 5 pt 5 ps .5   # red
set style line 4 lc rgb '#66A61E' lt 1 lw 5 pt 8 ps .5   #  green
set style line 6 lc rgb '#ff8800' lt 1 lw 5 pt 7 ps .5   #  orange
set style line 7 lc rgb '#9FAFDF' lt 1 lw 5 pt 8 ps 0.5   #royal blue
set style line 8 lc rgb '#1E90FF' lt 1 lw 5 pt 8 ps 0.5   #  slate blue
set style line 1 lc rgb '#F897C5' lt 1 lw 5 pt 8 ps 0.5   # taffy
set style line 5 lc rgb 'black' lt 1 lw 5 pt 5 ps 0.5   # black

set logscale x

set style fill transparent solid 0.25 # partial transparency
set style fill noborder # no separate top/bottom lines

set multiplot layout 1,3 margins 0.05,0.98,0.23,0.88 spacing 0.08,0.01


set ylabel "$p_{\\; \\textrm{CK}}$" 
set xlabel "$t / \\mathcal{T}$" 

set xrange[0.001:3.5]
set format x "10^{%L}"

set xtics(0.001, 0.01,0.1, 1) offset 0,0
set ytics 0,0.1,0.9

set title "State $A$"
set key at screen 0.3,0.05
plot  "../../data/MSM/pi_2010.dat" u ($0*50*0.00001):($1/0.000795) w l ls 5 lw 0 notitle , \
"" u ($0*50*0.00001):(($1-$4)/0.000795):(($1+$4)/0.000795) w filledcurves lc rgb 'black' notitle , \
"../../data/MSM/cktest_50_2010.dat" u ($1*50*0.00001):2 w l ls 2 lw 2 title "$\\tau = 0.0005 \\mathcal{T}$"  ,\
"../../data/MSM/cktest_200_2010.dat" u ($1*4*40*0.00001):2 w l ls 3 lw 2 notitle  ,\
 "../../data/MSM/cktest_1000_2010.dat" u ($1*20*40*0.00001):2 w l ls 4 lw 2 notitle ,\


unset yrange
unset ylabel
set title "State $B$"
set key at screen 0.63,0.05
plot  "../../data/MSM/pi_2010.dat" u ($0*50*0.00001):($2/0.0003405) w l ls 5 lw 0 notitle,\
"" u ($0*50*0.00001):(($2-$5)/0.0003405):(($2+$5)/0.0003405) w filledcurves lc rgb 'black' notitle,\
"../../data/MSM/cktest_50_2010.dat" u ($1*50*0.00001):3 w l ls 2 lw 2 notitle  ,\
"../../data/MSM/cktest_200_2010.dat" u ($1*4*50*0.00001):3 w l ls 3 lw 2  title "$\\tau = 0.002 \\mathcal{T}$" ,\
 "../../data/MSM/cktest_1000_2010.dat" u ($1*20*50*0.00001):3 w l ls 4 lw 2 notitle ,\

set title "State $C$"
set key at screen 0.96,0.05
plot "../../data/MSM/pi_2010.dat" u ($0*50*0.00001):($3/0.000775) w l ls 5 lw 0 notitle,\
"" u ($0*50*0.00001):(($3-$6)/0.000775):(($3+$6)/0.000775) w filledcurves lc rgb 'black' notitle,\
 "../../data/MSM/cktest_50_2010.dat" u ($1*50*0.00001):4 w l ls 2 lw 2 notitle  ,\
"../../data/MSM/cktest_200_2010.dat" u ($1*4*50*0.00001):4 w l ls 3 lw 2 notitle  ,\
 "../../data/MSM/cktest_1000_2010.dat" u ($1*20*50*0.00001):4 w l ls 4 lw 2  title "$\\tau = 0.01 \\mathcal{T}$" ,\

