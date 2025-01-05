reset

set terminal cairo size 15cm,7cm standalone color colortext header \
   "\\renewcommand\\familydefault{\\sfdefault}\\usepackage{cmbright}"

set output "lag_1001.tex"

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

set label '(a)' at screen 0.01,0.95
set label '(b)' at screen 0.51,0.95

set origin 0.0,0.0
set size 0.5,0.9

set xrange[-pi:pi]
set yrange[0.4:1.2]
set xlabel "$\\varphi / ( \\pi$ rad $)$"
set ylabel "$R_{14}$ / nm" 
set cblabel "$F/ \\epsilon$" offset 2.5,0

set rmargin 4

set key horizontal at screen 0.42,0.97

set palette defined (-35 "black", -25 "red", -15 "orange", 0  "yellow")

#set arrow front from 1.87,0.94 to 2.13,0.94 lw 5 lc rgb '#dd181f'
#set label front 'f'  at 1.95,0.85 textcolor rgb '#dd181f'

set xtics ("$-\\frac{2}{3}$" -2./3.*pi, "$-\\frac{1}{3}$" -1./3.*pi,"$0$" 0,"$\\frac{1}{3}$" 1./3.*pi ,"$\\frac{2}{3}$" 2./3.*pi)

plot \
"../../data/Frew/F_1001.dat" matrix u (phi($1)):(R($2)):3 w image notitle,\
"../../data/Frew/mem_1001.dat" u (phi($1)):(R($2)) pt 5 lc 1 ps 0.3 title "H",\
"../../data/Frew/mem_1001.dat" u (phi($5)):(R($6)) pt 5 lc 3 ps 0.3 title "E",\
"../../data/Frew/mem_1001.dat" u (phi($3)):(R($4)) pt 5 lc 2 ps 0.3 title "I",\
#"../../data/Frew/out_1001.dat" u (phi($1)):(R($2)) pt 5 lc 8 ps 0.3 title "del"


set origin 0.5,0.0
set size 0.5,1

#set lmargin 8
#set rmargin 1
#set tmargin 0.1
#set bmargin 4

set xrange[100:300]
set yrange[100:1300]
set xlabel "$\\tau /$ fs"
set ylabel "$t_i /$ fs"

set cbrange[-35:-1]

set xtics auto 

set key horizontal at screen 0.94,0.92

plot \
"../../data/Frew/1001.dat" u 1:2 w lp ls 2 title "$t_1\\;\\;\\; f=0~ \\epsilon / nm$",\
"../../data/Frew/1001.dat" u 1:3 w lp ls 7 title "$t_2 \\;\\;\\;f=0~ \\epsilon / nm$",\
"../../data/Frew/1901.dat" u 1:2 w lp ls 3 title "$t_1 \\;\\;\\;f=-9~ \\epsilon / nm$",\
"../../data/Frew/1901.dat" u 1:3 w lp ls 6 title "$t_2 \\;\\;\\;f=-9~ \\epsilon / nm$",\
x with filledcurves x1 lc "black" fs transparent solid 0.2 notitle

set style fill transparent solid 0.35 noborder
replot x w filledcurves ls 5 notitle

