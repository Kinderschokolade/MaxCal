reset

set output "potential_laser.tex"

set terminal cairo size 8cm,6cm standalone color colortext header \
   "\\renewcommand\\familydefault{\\sfdefault}\\usepackage{cmbright}"

set xlabel "$x~/\\mathcal{L}$"
set format x2 

set xrange [0:1]

set style line 1 lc rgb '#0060ad' lt 1 lw 5 pt 10 ps 0.9   # --- blue
set style line 2 lc rgb '#dd181f' lt 1 lw 5 pt 5 ps .9   # --- red
set style line 3 lc rgb '#66A61E' lt 1 lw 5 pt 8 ps .9   # ——— green
set style line 4 lc rgb '#ff8800' lt 1 lw 5 pt 7 ps .9   # ——— orange
set style line 5 lc rgb '#9FAFDF' lt 1 lw 5 pt 8 ps 0.9   # --- royal blue
set style line 6 lc rgb '#1E90FF' lt 1 lw 5 pt 8 ps 0.9   # --- slate blue

set style line 7 lc rgb '#0060ad' lt 1 lw 5 pt 7 ps 1.5   # --- blue
set style line 8 lc rgb 'black' lt 1 lw 5 pt 5 ps 1.5   # --- red

set ylabel "$U~/\\epsilon$"
set y2label "$f~/(\\epsilon/\\mathcal{L})$"
set ytics 1 nomirror
set yrange [0.:6.]

set y2range[0:1.05]
set y2tics 0.25
set y2tics format ""

set label '$f_{\text{max}}$' at screen 0.85,0.903

f(x) = exp(-(x-0.125)*(x-0.125)/(2*0.02*0.02))

set label '$A$' at 0.02,0.25
set label '$B$' at 0.28,4.8
set label '$C$' at 0.54,3.8
set label '$D$' at 0.8,0.8

#set key at 0.37,4.7
set key top right

plot "../../data/Urew/potential_9040.dat" u ($1):($2+0.6) ls 8  lw 6 w l title "potential", f(x) ls 2 lw 6  title "force" axis x1y2 
pause -1

