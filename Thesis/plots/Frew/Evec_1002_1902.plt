reset

set terminal cairo size 15cm,19.5cm standalone color colortext header \
   "\\renewcommand\\familydefault{\\sfdefault}\\usepackage{cmbright}"

set output "Evec_1002_1902.tex"

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

set multiplot



set xrange[-pi:pi]
set yrange[0.4:1.2]



set label "$f_R = 0$" at screen 0.19,0.96
set label "$f_R = 5$" at screen 0.45,0.96
set label "$f_R = 0 \\rightarrow f_R = 5$" at screen 0.7,0.96

#set label "$\\Phi_{0-2}$" at screen 0.2,0.95

set ylabel "$R_{14}/ nm$" 
set xlabel "$\\varphi / $( $\\pi$  rad )"

set origin 0.015,0.62
set size 0.37,0.31
unset colorbox

#set title "$\\Phi_0$"
set title "  "
set xtics ("$-\\frac{2}{3}$" -2./3.*pi, "$-\\frac{1}{3}$" -1./3.*pi,"$0$" 0,"$\\frac{1}{3}$" 1./3.*pi ,"$\\frac{2}{3}$" 2./3.*pi)


plot "../../data/Frew/Evec0_1002.dat" matrix u (phi($1)):(R($2)):3 w image notitle

set origin 0.33,0.62
set size 0.33,0.31
unset ytics
unset ylabel
set title "$\\Phi_0$"
plot "../../data/Frew/Evec0_1902.dat" matrix  u (phi($1)):(R($2)):3 w image notitle

set origin 0.627,0.62
set size 0.39,0.31
set colorbox
set cbtics 0.01

#set title "$\\Phi_0$"
set title "  "
plot "../../data/Frew/Evec0_1002_1902.dat" matrix  u (phi($1)):(R($2)):3 w image notitle


set palette defined (-1 "blue", 0 "white", 1 "red")

set cbrange[-1:1]
set cbtics 0.3 

set origin 0.015,0.31
set size 0.37,0.31
unset colorbox
set ytics

set ylabel "$R_{14}/ nm$" 
set title "  "
plot "../../data/Frew/Evec1_1002.dat" matrix u (phi($1)):(R($2)):3 w image notitle

set origin 0.33,0.31
set size 0.33,0.31
unset ytics
unset ylabel
set title "$\\Phi_1$"
plot "../../data/Frew/Evec1_1902.dat" matrix  u (phi($1)):(R($2)):3 w image notitle

set origin 0.627,0.31
set size 0.39,0.31
set colorbox


set title "  "
plot "../../data/Frew/Evec1_1002_1902.dat" matrix  u (phi($1)):(R($2)):(-$3) w image notitle



set origin 0.015,0.0
set size 0.37,0.31
unset colorbox

set ytics

set ylabel "$R_{14}/ nm$" 
set title "  "
plot "../../data/Frew/Evec2_1002.dat" matrix u (phi($1)):(R($2)):3 w image notitle

set origin 0.33,0.0
set size 0.33,0.31
unset ytics
unset ylabel
set title "$\\Phi_2$"
plot "../../data/Frew/Evec2_1902.dat" matrix  u (phi($1)):(R($2)):3 w image notitle

set origin 0.627,0.0
set size 0.39,0.31
set colorbox


set title "  "
plot "../../data/Frew/Evec2_1002_1902.dat" matrix  u (phi($1)):(R($2)):(-$3) w image notitle





#set label '(b)' at screen 0.01,0.58
#set label '(a)' at screen 0.01,0.95

