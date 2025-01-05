reset

set terminal cairo size 15cm,11cm standalone color colortext header \
   "\\renewcommand\\familydefault{\\sfdefault}\\usepackage{cmbright}"

set output "single_1030.tex"

set style line 2 lc rgb '#0060ad' lt 1 lw 5 pt 5 ps 0.5   #  blue
set style line 3 lc rgb '#dd181f' lt 1 lw 5 pt 7 ps .5   # red
set style line 4 lc rgb '#66A61E' lt 1 lw 5 pt 9 ps .5   #  green
set style line 6 lc rgb '#ff8800' lt 1 lw 5 pt 11 ps .5   #  orange
set style line 7 lc rgb '#9FAFDF' lt 1 lw 5 pt 13 ps 0.5   #royal blue
set style line 8 lc rgb '#1E90FF' lt 1 lw 5 pt 6 ps 0.5   #  slate blue
set style line 1 lc rgb '#F897C5' lt 1 lw 5 pt 8 ps 0.5   # taffy
set style line 5 lc rgb 'black' lt 1 lw 5 pt 10 ps 0.5   # black
set style line 9 lc rgb '#f1c40f' lt 1 lw 5 pt 12 ps 0.5   #yellow


set format x2 


set ytics 0.2

set xrange [0:3]

set multiplot

set origin 0.0,0.66666
set size 0.57,0.3333

set lmargin 8.1
set rmargin 0
set tmargin 1
set bmargin 0

set cbtics -7,2,0
set palette rgb 21,22,23
set cblabel "$U / \\epsilon$" offset 1.5,0

unset xtics

set arrow 1 front from 1.87,0.94 to 2.13,0.94 lw 5 lc rgb '#dd181f'
set label 4 front 'f'  at 1.95,0.85 textcolor rgb '#dd181f'


set ylabel "$y / \\mathcal{L}$"

plot \
"../../data/Urew/potential_ms_5010.dat" u 2:1:3 w image notitle  ,\
"../../data/Urew/mem_5010.dat" u 1:2 pt 5 lc 1 ps 0.5 notitle,\
"../../data/Urew/mem_5010.dat" u 3:4 pt 5 lc 2 ps 0.5 notitle,\
"../../data/Urew/mem_5010.dat" u 5:6 pt 5 lc 3 ps 0.5 notitle,\
((x>0.35 && x < 0.65)?1:NaN) with filledcurves x1 lc "black" fs transparent solid 0.2 notitle,\
((x>1.35 && x < 1.65)?1:NaN) with filledcurves x1 lc "black" fs transparent solid 0.2 notitle,\
((x>2.35 && x < 2.65)?1:NaN) with filledcurves x1 lc "black" fs transparent solid 0.2 notitle


unset label 4
unset arrow 1


set origin 0.0,0.0
set size 0.5,0.4

set lmargin 8
set rmargin 1
set tmargin 0.1
set bmargin 3


set key top left 
set ylabel "$p_s$"
set yrange [0:0.33]
set ytics 0.0,0.05,0.31

set xtics 0.5
set xlabel "$x / \\mathcal{L}$"
plot \
((x>0.35 && x < 0.65)?1:NaN) with filledcurves x1 lc "black" fs transparent solid 0.2 notitle,\
((x>1.35 && x < 1.65)?1:NaN) with filledcurves x1 lc "black" fs transparent solid 0.2 notitle,\
((x>2.35 && x < 2.65)?1:NaN) with filledcurves x1 lc "black" fs transparent solid 0.2 notitle,\
"../../data/Frew/p_1930_1030.dat" u ($0/30*3+3./60.):1  ls 2 notitle,\
"../../data/Frew/p_1030_1930.dat" u ($0/30*3+3./60.):1  ls 3 notitle  ,\
"../../data/Frew/ps_1030.dat" u ($0/30*3+3./60):1 w l ls 2 lw 5   smooth csplines notitle ,\
"../../data/Frew/ps_1930.dat" u ($0/30*3+3./60.):1 w l ls 3 lw 5  smooth csplines notitle 


unset xlabel


set origin 0.0,0.4
set size 0.5,0.266666
set bmargin 0.1
set tmargin 1
set ylabel "$U~/\\epsilon$"
set ytics 2
set yrange [-5:1]
set format x ''
unset xlabel

unset key

set label '$A$' at 0.45,-1 front
set label '$B$' at 1.45,-2.5 front 
set label '$C$' at 2.45,-4 front
set label 102 '(b)' at screen 0.01,0.63
set label 101 '(a)' front at screen 0.01,0.96
set label 103 '(c)' at screen 0.51,0.63
set label 104 '(d)' at screen 0.01,0.39


set arrow 2 from 0.5,-2 to 0.7,-2 lw 5 lc rgb '#dd181f'
set label 5 'f' front  at 0.58,-2.9 textcolor rgb '#dd181f'

plot \
((x>0.35 && x < 0.65)?1:NaN) with filledcurves x1 lc "black" fs transparent solid 0.2 notitle,\
((x>1.35 && x < 1.65)?1:NaN) with filledcurves x1 lc "black" fs transparent solid 0.2 notitle,\
((x>2.35 && x < 2.65)?1:NaN) with filledcurves x1 lc "black" fs transparent solid 0.2 notitle,\
"../../data/Frew/potential_ms_1030.dat" u ($1):($2-6.0) ls 2  lw 5 w l smooth csplines notitle  

unset arrow 2
unset label 5

set origin 0.5,0.09
set size 0.48,0.58

unset format

set logscale x
set xrange [1:1000]

set ytics 0,0.02,0.11
set xtics 1,10,1000

set yrange[0:0.11]


set ylabel "$p_{\\; \\textrm{FPT}}$" 
set xlabel "$t / \\tau$"


set key default
set key at screen 1.0,0.92

plot "../../data/Frew/hist_1930_1030.dat" u 0:6 ls 2 title "Rew 1D$: f=9~\\epsilon / \\mathcal{L} \\rightarrow f=0$" ,\
"../../data/Frew/hist_1030_1930.dat" u 0:6 ls 3 title "Rew 1D$: f=0 \\rightarrow f=9~\\epsilon / \\mathcal{L}$"   ,\
"../../data/Frew/hist_1030_1930.dat" u 0:15 w l   ls 3 lw 5 title "Simulation 1D$: f=0~\\epsilon / \\mathcal{L}$" ,\
"../../data/Frew/hist_1930_1030.dat" u 0:15 w l  ls 2 lw 5 title "Simulation 1D$: f=9~\\epsilon / \\mathcal{L}$" ,\
"../../data/Urew/hist_5910_5010.dat" u 0:6 ls 7 title "Rew 2D$: f=9~\\epsilon / \\mathcal{L} \\rightarrow f=0$" ,\
"../../data/Urew/hist_5010_5910.dat" u 0:6 ls 6 title "Rew 2D$: f=0 \\rightarrow f=9~\\epsilon / \\mathcal{L}$"   ,\
"../../data/Urew/hist_5910_5010.dat" u 0:15 w l  ls 7 dt 2 lw 5 title "Simulation 2D$: f=0~\\epsilon / \\mathcal{L}$" ,\
"../../data/Urew/hist_5010_5910.dat" u 0:15 w l ls 6 dt 2 lw 5 title "Simulation 2D$: f=9~\\epsilon / \\mathcal{L}$" ,\


#plot "../../data/Frew/hist_1930_1030.dat" u 0:6 ls 2 ps 1 notitle ,\
#"../../data/Frew/hist_1030_1930.dat" u 0:6 ls 3 ps 1 notitle  ,\
#"../../data/Frew/hist_1930_1030.dat" u 0:15 w l  ls 2 lw 5 notitle ,\
#"../../data/Frew/hist_1030_1930.dat" u 0:15 w l  ls 3 lw 5 notitle   



