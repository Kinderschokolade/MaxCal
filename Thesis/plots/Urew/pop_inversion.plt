reset

set terminal cairo size 3.5,2.5 standalone color colortext header \
   "\\renewcommand\\familydefault{\\sfdefault}\\usepackage{cmbright}"
set output "population_inversion.tex"


set style line 1 lc rgb '#0060ad' lt 1 lw 5 pt 10 ps 0.9   # --- blue
set style line 2 lc rgb '#dd181f' lt 1 lw 5 pt 5 ps .9   # --- red
set style line 3 lc rgb '#66A61E' lt 1 lw 5 pt 8 ps .9   # ——— green
set style line 4 lc rgb '#ff8800' lt 1 lw 5 pt 7 ps .9   # ——— orange
set style line 5 lc rgb '#9FAFDF' lt 1 lw 5 pt 8 ps 0.9   # --- royal blue
set style line 6 lc rgb '#1E90FF' lt 1 lw 5 pt 8 ps 0.9   # --- slate blue

set style line 7 lc rgb '#0060ad' lt 1 lw 5 pt 7 ps 1.5   # --- blue
set style line 8 lc rgb 'black' lt 1 lw 5 pt 5 ps 1.5   # --- red

set xlabel "$f~[\\epsilon / \\mathcal{L}]$
set ylabel "$\\pi_i$"

set xrange[0:250]


plot "../data/P_9040.dat" u 1:2 title "A" ls 1 ,\
"../data/P_9040.dat" u 1:3 title "B" ls 2 ,\
"../data/P_9040.dat" u 1:4 title "C" ls 3 ,\
"../data/P_9040.dat" u 1:5 title "D" ls 4 ,\
"../data/ps_cts_9040.dat" u 1:2 w l ls 1 notitle,\
"../data/ps_cts_9040.dat" u 1:3 w l ls 2 notitle,\
"../data/ps_cts_9040.dat" u 1:4 w l ls 3 notitle,\
"../data/ps_cts_9040.dat" u 1:5 w l ls 4 notitle

pause -1

