reset

set terminal cairo size 15cm,8cm standalone color colortext header \
   "\\renewcommand\\familydefault{\\sfdefault}\\usepackage{cmbright}"

set output "single_2001.tex"

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
set xlabel "$\\varphi /$ $ ( \\pi$ rad$)$"
set ylabel "$R_{14}$ / nm" 

set format x2 


set multiplot

#set palette rgb 3,11,6

set origin 0.0,0.5
set size 0.5,0.5

#set lmargin 8
#set rmargin 1
#set tmargin 8
#set bmargin 4
set cbtics 0,0.01,0.05
set cblabel "$p_s$" offset 0.5

set xtics ("$-\\frac{2}{3}$" -2./3.*pi, "$-\\frac{1}{3}$" -1./3.*pi,"$0$" 0,"$\\frac{1}{3}$" 1./3.*pi ,"$\\frac{2}{3}$" 2./3.*pi)

set ytics 0.2

#set title "Simulation$: f=5~\\epsilon / 2\\pi \\mathcal{L} "

plot "../../data/Frew/Evec0_2501.dat" matrix  u (phi($1)):(R($2)):3 w image  notitle,\


set origin 0.0,0.0
set size 0.5,0.5
#set tmargin 1
#set bmargin 0.1


#set title "Reweight$: f=0  \\rightarrow f=5~\\epsilon /$  rad"

plot "../../data/Frew/Evec0_2001_2501.dat" matrix u (phi($1)):(R($2)):3 w image notitle  ,\


set label '(b)' at screen 0.01,0.49
set label '(a)' at screen 0.01,0.97
set label '(c)' at screen 0.50,0.97


set origin 0.5,0.2
set size 0.5,0.8

unset format

unset title
set logscale x
set xrange [1:120]


set xtics 1,10,1000

set yrange[0:0.069]
set ytics 0.01

set ylabel "$p_{\\; \\textrm{FPT}}$" 
set xlabel "$t / \\tau$"


set key default
set key at screen 0.96,0.2

plot "../../data/Frew/hist_2501_2001.dat" u 0:3 ls 2 title "Reweight$: f=5~\\epsilon /$ rad $\\rightarrow f=0$" ,\
"../../data/Frew/hist_2001_2501.dat" u 0:3 ls 3 title "Reweight$: f=0 \\rightarrow f=5~\\epsilon /$ rad"   ,\
"../../data/Frew/hist_2501_2001.dat" u 0:12 w l  ls 2 lw 5 title "Simulation$: f=5~\\epsilon / $ rad" ,\
"../../data/Frew/hist_2001_2501.dat" u 0:12 w l  ls 3 lw 5 title "Simulation$: f=0$" ,\


