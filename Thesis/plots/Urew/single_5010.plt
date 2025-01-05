reset

set terminal cairo size 15cm,17cm standalone color colortext header \
   "\\renewcommand\\familydefault{\\sfdefault}\\usepackage{cmbright}"

set output "single_5010.tex"

set style line 2 lc rgb '#0060ad' lt 1 lw 5 pt 5 ps 0.5   #  blue
set style line 3 lc rgb '#dd181f' lt 1 lw 5 pt 7 ps .5   # red
set style line 4 lc rgb '#66A61E' lt 1 lw 5 pt 9 ps .5   #  green
set style line 6 lc rgb '#ff8800' lt 1 lw 5 pt 11 ps .5   #  orange
set style line 7 lc rgb '#9FAFDF' lt 1 lw 5 pt 13 ps 0.5   #royal blue
set style line 8 lc rgb '#1E90FF' lt 1 lw 5 pt 6 ps 0.5   #  slate blue
set style line 1 lc rgb '#F897C5' lt 1 lw 5 pt 8 ps 0.5   # taffy
set style line 5 lc rgb 'black' lt 1 lw 5 pt 10 ps 0.5   # black
set style line 9 lc rgb '#f1c40f' lt 1 lw 5 pt 12 ps 0.5   #yellow


set xlabel "$x~/\\mathcal{L}$"
set ylabel "$y~/\\mathcal{L}$"
set format x2 

#set xtics 0.2
#set ytics 0.01

#set xrange [0:1]

set multiplot

#set palette rgb 3,11,6

set origin 0.0,0.33
set size 0.5,0.33

#set lmargin 8
#set rmargin 1
#set tmargin 8
#set bmargin 4
set cbtics 0,0.02,0.08
set cblabel "$p_s$" offset 0.5

set title "Simulation$: f=9~\\epsilon / \\mathcal{L} "

plot "../../data/Urew/ps_5910.dat" matrix  u ($2/10+1./20.):($1/10+1./20.):3 w image  notitle,\


set origin 0.0,0.0
set size 0.5,0.33
#set tmargin 1
#set bmargin 0.1

set xrange[0:3]
set yrange[0:1]


set title "Reweight$: f=0~\\epsilon / \\mathcal{L} \\rightarrow f=9~\\epsilon / \\mathcal{L}$"

plot "../../data/Urew/p_5010_5910.dat" matrix u ($2/10+1./20.):($1/10+1./20.):3 w image notitle  ,\




set origin 0.0,0.66
set size 0.5,0.33
#set bmargin 0.1
#set tmargin 1
#set ytics 2
#set yrange [-1:6]
#set format x ''
#unset xlabel
unset title


set label 1 front '$A$' at 0.4,0.3
set label 2 front '$B$' at 1.4,0.3
set label 3 front '$C$' at 2.4,0.3
set label 111 '(a)' at screen 0.01,0.97
set label 112 '(c)' at screen 0.01,0.62
set label 113 '(b)' at screen 0.51,0.97
set label 114 '(d)' at screen 0.51,0.4


set cbtics -7,2,0
set palette rgb 21,22,23
set cblabel "$U / \\epsilon$" offset 1.4,0


set arrow 1 front from 1.87,0.94 to 2.13,0.94 lw 5 lc rgb '#dd181f'
set label 4 front 'f'  at 1.95,0.85 textcolor rgb '#dd181f'


plot \
"../../data/Urew/potential_ms_5010.dat" u 2:1:3 w image notitle  ,\
"../../data/Urew/mem_5010.dat" u 1:2 pt 5 lc 1 ps 0.5 notitle,\
"../../data/Urew/mem_5010.dat" u 3:4 pt 5 lc 2 ps 0.5 notitle,\
"../../data/Urew/mem_5010.dat" u 5:6 pt 5 lc 3 ps 0.5 notitle,\

unset label 4
unset arrow 1


set origin 0.5,0.1
set size 0.5,0.3

unset format

unset title
set logscale x
set xrange [1:1000]

set ytics 0,0.01,00.081
set xtics 1,10,1000

set yrange[0:*]


set ylabel "$p_{\\; \\textrm{FPT}}$" 
set xlabel "$t / \\tau$"


set key default
set key at screen 0.96,0.1

plot "../../data/Urew/hist_5910_5010.dat" u 0:6 ls 2 title "Reweight$: f=9~\\epsilon / \\mathcal{L} \\rightarrow f=0~\\epsilon / \\mathcal{L}$" ,\
"../../data/Urew/hist_5010_5910.dat" u 0:6 ls 3 title "Reweight$: f=0~\\epsilon / \\mathcal{L} \\rightarrow f=9~\\epsilon / \\mathcal{L}$"   ,\
"../../data/Urew/hist_5010_5910.dat" u 0:15 w l  ls 3 lw 5 title "Simulation$: f=0~\\epsilon / \\mathcal{L}$" ,\
"../../data/Urew/hist_5910_5010.dat" u 0:15 w l  ls 2 lw 5 title "Simulation$: f=9~\\epsilon / \\mathcal{L}$" ,\

set origin 0.5,0.4
set size 0.5,0.3

unset format

unset title
unset logscale
unset label 1
unset label 2
unset label 3

set label 10 "$f=0~\\epsilon / \\mathcal{L}$" at 0.2,0.05

set xrange [0:3]

set xtics 0.5

set logscale y
set yrange[0.0001:0.14]
set ytics auto

set ylabel "$p_s$" 
set xlabel "$x / \\mathcal{L}$"

set key default
set key outside top horizontal

plot "../../data/Urew/ps_5010.dat" u ($0/10.+ 1./20):5 w l ls 4 title "Sim $y=0.45\\mathcal{L}$",\
"../../data/Urew/ps_5010.dat" u ($0/10.+ 1./20):3 w l ls 6 title "$y=0.25 \\mathcal{L}$",\
"../../data/Urew/ps_5010.dat" u ($0/10.+ 1./20):1 w l ls 7 title "$y=0.05\\mathcal{L}$",\
"../../data/Urew/p_5910_5010.dat" u ($0/10.+ 1./20):1 w p ls 7 ps 0.3 notitle ,\
"../../data/Urew/p_5910_5010.dat" u ($0/10.+ 1./20):3 w p ls 6 ps 0.3 notitle ,\
"../../data/Urew/p_5910_5010.dat" u ($0/10.+ 1./20):5 w p ls 4 ps 0.3 notitle ,\

set origin 0.5,0.7
set size 0.5,0.3

unset format

unset label 10
set label 11 "$f=9~\\epsilon / \\mathcal{L}$" at 0.2,0.05

unset title
unset logscale
set xrange [0:3]

set xtics 0.5

set logscale y
set yrange[0.0001:0.14]
set ytics auto

set ylabel "$p_s$" 
set xlabel "$x / \\mathcal{L}$"

set key outside bottom horizontal

plot "../../data/Urew/ps_5910.dat" u ($0/10.+ 1./20):5 w l ls 4 notitle,\
"../../data/Urew/ps_5910.dat" u ($0/10.+ 1./20):3 w l ls 6  notitle,\
"../../data/Urew/ps_5910.dat" u ($0/10.+ 1./20):1 w l ls 7  notitle,\
"../../data/Urew/p_5010_5910.dat" u ($0/10.+ 1./20):5 w p ls 4 ps 0.3 title "Rew $y=0.45\\mathcal{L}$"  ,\
"../../data/Urew/p_5010_5910.dat" u ($0/10.+ 1./20):3 w p ls 6 ps 0.3 title "$y=0.25\\mathcal{L}$"  ,\
"../../data/Urew/p_5010_5910.dat" u ($0/10.+ 1./20):1 w p ls 7 ps 0.3 title "$y=0.05\\mathcal{L}$" ,\

