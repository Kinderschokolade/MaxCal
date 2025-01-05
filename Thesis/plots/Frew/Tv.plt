reset

set terminal cairo size 15cm,8cm standalone color colortext header \
   "\\renewcommand\\familydefault{\\sfdefault}\\usepackage{cmbright}"

set output "Tv.tex"

set style line 2 lc rgb '#0060ad' lt 1 lw 5 pt 5 ps 0.5   #  blue
set style line 3 lc rgb '#dd181f' lt 1 lw 5 pt 7 ps .5   # red
set style line 4 lc rgb '#66A61E' lt 1 lw 5 pt 9 ps .5   #  green
set style line 6 lc rgb '#ff8800' lt 1 lw 5 pt 11 ps .5   #  orange
set style line 7 lc rgb '#9FAFDF' lt 1 lw 5 pt 13 ps 0.5   #royal blue
set style line 8 lc rgb '#1E90FF' lt 1 lw 5 pt 6 ps 0.5   #  slate blue
set style line 1 lc rgb '#F897C5' lt 1 lw 5 pt 8 ps 0.5   # taffy
set style line 5 lc rgb 'black' lt 1 lw 5 pt 10 ps 0.5   # black
set style line 9 lc rgb '#f1c40f' lt 1 lw 5 pt 12 ps 0.5   #yellow

msR= 17
msp= 15
phi(x) =  x/15.* 2.*pi - pi + 2*pi/30
R(x) = x /15 * 0.7 + 0.45 - 0.7/30.0

set xrange[-pi:pi]
set yrange[0.4:1.2]
set xlabel "$\\varphi /$ $(\\pi$ rad)"
set ylabel "$R_{14}$ / nm" 

set format x2 


set label '(c)' at screen 0.005,0.49
set label '(a)' at screen 0.005,0.97
set label '(b)' at screen 0.50,0.97
set label '(d)' at screen 0.50,0.49


set multiplot

set palette defined ( 0 "black", 0.05 "orange",  0.1 "red" ) 

set origin 0.0,0.5
set size 0.5,0.5

set cbtics 0,0.01,0.06
set cbrange[0:0.063]
set cblabel "$p_s$" offset 0.5

set xtics ("$-\\frac{2}{3}$" -2./3.*pi, "$-\\frac{1}{3}$" -1./3.*pi,"$0$" 0,"$\\frac{1}{3}$" 1./3.*pi ,"$\\frac{2}{3}$" 2./3.*pi)

set ytics 0.2

set object 1 circle front at first -2.1,0.8 radius char 0.3 fillcolor rgb 'green' fillstyle solid noborder

set arrow 3 front from 0.653,1.2 to 0.653,0.4 nohead lc rgb 'red' dt "." lw 5



plot "../../data/Frew/Tv_2001.dat" every :::42::42  u (phi($2)):(R($1)):($3>0.000?$3:NaN) w image  notitle,\
(x < -2.30?1.2:NaN) with filledcurves x1 lc "green" fs transparent solid 0.2 notitle,\
((x>-1.885 && x < 0.653)?1.2:NaN) with filledcurves x1 lc "blue" fs transparent solid 0.2 notitle,\
(x > 0.653?1.2:NaN) with filledcurves x1 lc "green" fs transparent solid 0.2 notitle



set origin 0.5,0.5
set size 0.5,0.5

plot "../../data/Frew/Tv_2901.dat" every :::42::42 u (phi($2)):(R($1)):($3>0.001?$3:NaN) w image notitle  ,\
(x < -2.30?1.2:NaN) with filledcurves x1 lc "green" fs transparent solid 0.2 notitle,\
((x>-1.885 && x < 0.653)?1.2:NaN) with filledcurves x1 lc "blue" fs transparent solid 0.2 notitle,\
(x > 0.653?1.2:NaN) with filledcurves x1 lc "green" fs transparent solid 0.2 notitle



set origin 0.0,0.0
set size 0.5,0.5

unset object 1
unset arrow 3
set object 2 circle front at first -0.85,0.565 radius char 0.3 fillcolor rgb 'green' fillstyle solid noborder

set arrow 4 front from 1.883,1.2 to 1.883,0.4 nohead lc rgb 'red' dt "." lw 5
set cbrange[0:0.11]
set cbtics 0,0.02,0.1


plot "../../data/Frew/Tv_2001.dat" every :::87::87 u (phi($2)):(R($1)):($3>0.000?$3:NaN) w image notitle  ,\
(x < -1.047?1.2:NaN) with filledcurves x1 lc "green" fs transparent solid 0.2 notitle,\
((x>-0.628 && x < 1.883)?1.2:NaN) with filledcurves x1 lc "blue" fs transparent solid 0.2 notitle,\
(x > 1.883 ?1.2:NaN) with filledcurves x1 lc "green" fs transparent solid 0.2 notitle



set origin 0.5,0.0
set size 0.5,0.5
plot "../../data/Frew/Tv_2901.dat" every :::87::87 u (phi($2)):(R($1)):($3>0.000?$3:NaN) w image notitle  ,\
(x < -1.047?1.2:NaN) with filledcurves x1 lc "green" fs transparent solid 0.2 notitle,\
((x>-0.628 && x < 1.883)?1.2:NaN) with filledcurves x1 lc "blue" fs transparent solid 0.2 notitle,\
(x > 1.883 ?1.2:NaN) with filledcurves x1 lc "green" fs transparent solid 0.2 notitle



