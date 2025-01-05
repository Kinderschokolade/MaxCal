reset

set terminal cairo size 3.8,2.5 standalone color colortext header \
    "\\renewcommand\\familydefault{\\sfdefault}\\usepackage{cmbright}"
set output "potential_only.tex"


set xlabel "$x~[\\mathcal{L}]$"

set xtics 0.2
set ytics 0.01

set xrange [0:1]

set style line 1 lc rgb '#0060ad' lt 1 lw 5 pt 10 ps 0.9   # --- blue
set style line 2 lc rgb '#dd181f' lt 1 lw 5 pt 5 ps .9   # --- red
set style line 3 lc rgb '#66A61E' lt 1 lw 5 pt 8 ps .9   # ——— green
set style line 4 lc rgb '#ff8800' lt 1 lw 5 pt 7 ps .9   # ——— orange
set style line 5 lc rgb '#9FAFDF' lt 1 lw 5 pt 8 ps 0.9   # --- royal blue
set style line 6 lc rgb '#1E90FF' lt 1 lw 5 pt 8 ps 0.9   # --- slate blue

set style line 7 lc rgb '#0060ad' lt 1 lw 5 pt 7 ps 1.5   # --- blue
set style line 8 lc rgb 'black' lt 1 lw 5 pt 5 ps 1.5   # --- red

set grid x

set key top left 
set ylabel "$p_s$"
set yrange [0.006:0.0299]
set ytics 0.01,0.01,0.021

plot \
"../data/ps_6000.dat" u ($0/60+1./120):1 w l ls 3 lw 6 dt 1 lc 8 title "$p_{\\textup{equilibrium}}$",\

