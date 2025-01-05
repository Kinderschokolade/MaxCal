

set terminal cairo size  9cm,10cm standalone color colortext header \
   "\\renewcommand\\familydefault{\\sfdefault}\\usepackage{cmbright}"

set output 'Evec1.tex'

set style line 2 lc rgb '#0060ad' lt 1 lw 5 pt 5 ps 0.5   #  blue
set style line 3 lc rgb '#dd181f' lt 1 lw 5 pt 5 ps .5   # red
set style line 4 lc rgb '#66A61E' lt 1 lw 5 pt 8 ps .5   #  green
set style line 6 lc rgb '#ff8800' lt 1 lw 5 pt 8 ps .5   #  orange
set style line 7 lc rgb '#9FAFDF' lt 1 lw 5 pt 8 ps 0.5   #royal blue
set style line 8 lc rgb '#1E90FF' lt 1 lw 5 pt 8 ps 0.5   #  slate blue
set style line 1 lc rgb '#F897C5' lt 1 lw 5 pt 8 ps 0.5   # taffy
set style line 5 lc rgb 'black' lt 1 lw 5 pt 5 ps 0.5   # black


set multiplot layout 2,1

set key at 16,0.078

set yrange[0:0.08]

set label "a)" at graph -0.185,1



set ytics 0.02
set ylabel "$p_s$"
plot ((x>(9.5) && x < 15.5)?1:NaN) with filledcurves x1 lc "black" fs transparent solid 0.2 notitle,\
((x>(29.5) && x < 39.5)?1:NaN) with filledcurves x1 lc "black" fs transparent solid 0.2 notitle,\
((x>(53.5) && x < 59.5)?1:NaN) with filledcurves x1 lc "black" fs transparent solid 0.2 notitle,\
"../../data/MSM/ps_2010.dat" u 0:($1) w lp ls 5 notitle,\
"" u ($0+10):1 every ::10::15  ls 4  title "$A$",\
"" u ($0+30):1 every ::30::39  ls 6 title "$B$",\
"" u ($0+54):1 every ::54::59  ls 7 title "$C$"

set ytics 0.2
set ylabel "$ \\Psi_i$"
set xlabel "\\# microstate"

set key at 16.5,0.29

unset label
set label "b)" at graph -0.185,1
set yrange[-0.4:0.31]
plot ((x>(9.5) && x < 15.5)?1:NaN) with filledcurves x1 lc "black" fs transparent solid 0.2 notitle,\
((x>(29.5) && x < 39.5)?1:NaN) with filledcurves x1 lc "black" fs transparent solid 0.2 notitle,\
((x>(53.5) && x < 59.5)?1:NaN) with filledcurves x1 lc "black" fs transparent solid 0.2 notitle,\
((x>(9.5) && x < 15.5)?-1:NaN) with filledcurves x1 lc "black" fs transparent solid 0.2 notitle,\
((x>(29.5) && x < 39.5)?-1:NaN) with filledcurves x1 lc "black" fs transparent solid 0.2 notitle,\
((x>(53.5) && x < 59.5)?-1:NaN) with filledcurves x1 lc "black" fs transparent solid 0.2 notitle,\
"../../data/MSM/Evec2_com_2010.dat" u 0:1 w lp ls 2  title "$\\Psi_1$",\
"../../data/MSM/Evec3_com_2010.dat" u 0:1 w lp ls 3  title "$\\Psi_2$"





pause -1
