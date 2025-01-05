reset

set terminal cairo size 15cm,8cm standalone color colortext header \
   "\\renewcommand\\familydefault{\\sfdefault}\\usepackage{cmbright}"

set output "Invariant_2010.tex"

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
set size 0.5,1.0
set xlabel "$U~/\\epsilon$"
set xtics 0,1,5
set xrange [-1:6]

set ytics 0.2



set yrange[0:1]

set ylabel "$y/\\mathcal{L}$"


set arrow nohead lt 8 lw 5 dt 2 front from screen 0.22,0.3 to screen 0.53,0.3
set arrow nohead lt 8 lw 5 dt 2 front from screen 0.22,0.38 to screen 0.57,0.38

set arrow nohead lt 8 lw 5 dt 2 front from screen 0.16,0.575 to screen 0.67,0.575
set arrow nohead lt 8 lw 5 dt 2 front from screen 0.16,0.655 to screen 0.71,0.655

set arrow nohead lt 8 lw 5 dt 2 front from screen 0.22,0.85 to screen 0.79,0.85
set arrow nohead lt 8 lw 5 dt 2 front from screen 0.22,0.93 to screen 0.83,0.93

#set label '(b)' at screen 0.01,0.58
#set label '(a)' at screen 0.01,0.95
#set label '(c)' at screen 0.51,0.79

unset key
set parametric
set style fill  transparent solid 0.50 noborder
plot \
"../../data/Complementary/potential_2010.dat" u ($2+1):($1+1./200) w filledcurves x1=-1  ls 2  lw 5  notitle  ,\


unset parametric
unset style


set origin 0.45,0.0
set size 0.52,1.0

set xrange[0:1]
set yrange[0:1]

set xtics 0.2

set xlabel "$x/\\mathcal{L}$"
set cblabel "$I_{ij} / 10^{-3}$" offset 1.3,0

unset ylabel
unset ytics

#set format cb "%s*10^{%S}"
set cbtics ("1" 0.001,"2" 0.002,"3" 0.003)

set palette defined ( 0 "white", 0.003  "#0060ad") 

plot "../../data/Complementary/Invariant_2010.dat" matrix u ($1/60.+1./120.):($2/60.+1./120.):3 w image


