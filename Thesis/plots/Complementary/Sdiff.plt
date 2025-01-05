reset

set terminal cairo size 15cm,8cm standalone color colortext header \
   "\\renewcommand\\familydefault{\\sfdefault}\\usepackage{cmbright}"

set output "Sdiff.tex"

set style line 2 lc rgb '#0060ad' lt 1 lw 5 pt 5 ps 0.5   #  blue
set style line 3 lc rgb '#dd181f' lt 1 lw 5 pt 7 ps .5   # red
set style line 4 lc rgb '#66A61E' lt 1 lw 5 pt 9 ps .5   #  green
set style line 6 lc rgb '#ff8800' lt 1 lw 5 pt 11 ps .5   #  orange
set style line 7 lc rgb '#9FAFDF' lt 1 lw 5 pt 13 ps 0.5   #royal blue
set style line 8 lc rgb '#1E90FF' lt 1 lw 5 pt 6 ps 0.5   #  slate blue
set style line 1 lc rgb '#F897C5' lt 1 lw 5 pt 8 ps 0.5   # taffy
set style line 5 lc rgb 'black' lt 1 lw 5 pt 10 ps 0.5   # black
set style line 9 lc rgb '#f1c40f' lt 1 lw 5 pt 12 ps 0.5   #yellow


set multiplot

set origin 0.0,0.0
set size 0.45,1


set label "(a)" at screen 0.01,0.97
set label "(b)" at screen 0.46,0.97


set xlabel "from $x / \\mathcal{L}$"
set ylabel "to $x / \\mathcal{L}$"

set xrange[0:1]
set yrange[0:1]

set palette defined (-0.1 "blue", 0 "black" , 0.1 "red")
set cbrange[-0.1:0.1]

unset colorbox

plot "../../data/Complementary/Sprod_diff_2010.dat" matrix u  ($1/60.+1./120.):($2/60.+1./120.):3 w image notitle

set colorbox
set origin 0.45,0.0
set size 0.54,1

set ytics 0.2

set xlabel "$x / \\mathcal{L}$"
set ylabel "$y / \\mathcal{L}$"



set xrange[0:1.2]
set yrange[0:1]

set cblabel "$\\Delta S~/  (\\epsilon / K)$"
set object 2 circle front at first 0.55,0.35 radius char 0.3 fillcolor rgb 'green' fillstyle solid noborder


plot "../../data/Complementary/Sprod_vdiff_5010.dat" u ($1/10.+1./20.):($2/10.+1./20.):3 w image notitle





