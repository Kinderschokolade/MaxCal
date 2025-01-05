reset

set terminal cairo size 15cm,9cm standalone color colortext header \
   "\\renewcommand\\familydefault{\\sfdefault}\\usepackage{cmbright}"

set output "single_9040.tex"

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
set format x2 


set xtics 0.2
set ytics 0.01

set xrange [0:1]

set multiplot

set origin 0.,0
set size 0.5,0.85

set lmargin 8
set rmargin 1
set tmargin 0.1
set bmargin 4


set key top left 
set ylabel "$p_s$"
set yrange [0:0.065]
set ytics 0.0,0.01,0.061

set key at screen 0.1,0.96
plot \
"../../data/Urew/p_9840_9040.dat" u ($0/30+1./120.):1 every 2   ls 2 notitle,\
"../../data/Urew/p_9040_9840.dat" u ($0/30+1./120.):1 every 2   ls 3 notitle  ,\
"../../data/Urew/ps_9040.dat" u ($0/30+1./120):1 every 2 w l ls 2 lw 5  title "Simulation$: f=0~\\epsilon / \\mathcal{L}$",\
"../../data/Urew/ps_9840.dat" u ($0/30+1./120.):1 every 2 w l ls 3 lw 5  smooth csplines title "Simulation$: f=200~\\epsilon / \\mathcal{L}$",\

set label '(b)' at screen 0.51,0.85
set label '(a)' at screen 0.01,0.85

set origin 0.5,0
set size 0.48,0.85

unset format

set logscale x
set logscale y
set xrange [1:*]
set yrange[*:*]

set format y "10^{%S}"
set xtics 1,10,10000
set ytics 1e-17,100,1e-02



set ylabel "$p_{\\; \\textrm{FPT}}$" 
set xlabel "$t / \\tau$"


set key default
set key at screen 1,0.96

plot "../../data/Urew/hist_9840_9040.dat" u 0:3 ls 2 title "Reweight$: f=200~\\epsilon / \\mathcal{L} \\rightarrow f=0~\\epsilon / \\mathcal{L}$" ,\
"../../data/Urew/hist_9040_9840.dat" u 0:3 ls 3 title "Reweight$: f=0~\\epsilon / \\mathcal{L} \\rightarrow f=200~\\epsilon / \\mathcal{L}$"   ,\
"../../data/Urew/hist_9040_9840.dat" u 0:19 w l  ls 3 lw 5 notitle  ,\
"../../data/Urew/hist_9840_9040.dat" u 0:19 w l  ls 2 lw 5 notitle ,\



